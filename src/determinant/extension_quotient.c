/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_quotient.h"
#include "dixon_flint.h"
#include <setjmp.h>
#include <stdarg.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "fq_poly_mat_det.h"
#define DEGREE_QUOTIENT_EXTENSION
#include "degree_quotient_extension_field.h"
#define DEGREE_QUOTIENT_LIBRARY
#include "degree_quotient_impl.h"

/* Only dimension admission uses this escape; no matrices exist at that point.
 * All later rank/degree failures take the ordinary cleanup path. */
static jmp_buf admission;
static const char *admission_reason;
static void fail(const char *message)
{
    admission_reason=message;
    longjmp(admission,1);
}
static int admit_dimensions(const char *spec)
{
    if(setjmp(admission))return 0;
    init_dimensions(spec,0);
    return 1;
}
static void info(const char *fmt,...)
{
    if(g_dixon_verbose_level<1)return;
    va_list args;va_start(args,fmt);vprintf(fmt,args);va_end(args);fflush(stdout);
}
/* Share monomial keys, but copy the established quotient basis and allocate
 * private normal forms. No worker can change another worker's reduction. */
static degree_worker clone_worker(const degree_worker *source)
{
    degree_worker out={0};out.grades=flint_calloc(limit+1,sizeof(grade));
    for(slong d=0;d<=limit;d++) {
        const grade *src=source->grades+d;grade *dst=out.grades+d;
        dst->size=src->size;dst->h=src->h;dst->keys=src->keys;
        dst->free=flint_calloc(dst->h,sizeof(slong));
        if(d<=right_degree)memcpy(dst->free,src->free,dst->h*sizeof(slong));
        dq_mat_init(dst->reduce,dst->size,dst->h,prime);
    }
    return out;
}
static void clear_worker_copy(degree_worker *worker)
{
    for(slong d=0;d<=limit;d++) {
        flint_free(worker->grades[d].free);dq_mat_clear(worker->grades[d].reduce);
    }
    flint_free(worker->grades);
}
static dq_coeff **new_coefficients(void)
{
    dq_coeff **q=flint_calloc(n,sizeof(*q));
    for(slong f=0;f<n;f++)q[f]=dq_vec_new(spaces[degrees[f]].size);
    return q;
}
static void clear_coefficients(dq_coeff **q)
{
    for(slong f=0;f<n;f++)dq_vec_free(q[f],spaces[degrees[f]].size);flint_free(q);
}
/* Limit concurrency only by the requested threads and available samples. */
static slong sample_workers(slong requested,slong points)
{
    slong count=FLINT_MIN(FLINT_MAX(1,requested),points);
#ifdef _OPENMP
    if(omp_in_parallel())count=1;
#else
    count=1;
#endif
    return count;
}
static int sample_cores(dq_coeff *values,const dq_coeff *nodes,slong points,slong h,slong workers,slong *used)
{
    degree_worker *contexts=flint_calloc(workers,sizeof(*contexts));
    contexts[0]=active_worker;
    for(slong w=1;w<workers;w++)contexts[w]=clone_worker(contexts);
    int valid=1;
    #pragma omp parallel num_threads(workers) if(workers>1) reduction(&:valid)
    {
        slong id=0;
#ifdef _OPENMP
        id=omp_get_thread_num();
        if(id==0)*used=omp_get_num_threads();
#else
        *used=1;
#endif
        degree_worker saved=active_worker;active_worker=contexts[id];
        dq_coeff **q=new_coefficients();
        dq_mat_t core;dq_mat_init(core,h,h,prime);
        #pragma omp for schedule(dynamic,1)
        for(slong sample=0;sample<points;sample++) {
            specialize(q,nodes[sample]);
            if(!normal_forms(q,0)){valid=0;continue;}
            entry *entries=difference_matrix(q);bezout_core(core,entries,1);
            for(slong k=0;k<n*n;k++)flint_free(entries[k].terms);flint_free(entries);
            for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)
                dq_set(values[(sample*h+i)*h+j],dq_mat_entry(core,i,j));
        }
        dq_mat_clear(core);clear_coefficients(q);active_worker=saved;
    }
    for(slong w=1;w<workers;w++)clear_worker_copy(contexts+w);
    flint_free(contexts);return valid;
}
#include "quotient_extension_interpolation.h"

static int compute(unified_mpoly_struct *result,const unified_mpoly_struct *polys,
                   slong nv,slong np,slong requested_workers,mq_quotient_stats *stats)
{
    double begin=seconds();memset(stats,0,sizeof(*stats));
    stats->reason="extension quotient requires at least two equations over an extension field and one parameter";
    if(np!=1 || nv<1 || fq_nmod_ctx_degree(polys[0].ctx)<=1 ||
       g_field_equation_reduction || g_field_equation_final_only)return 0;
    stats->reason="subset workspace exceeds machine address capacity";
    if(nv>=FLINT_BITS-2)return 0;
    n=nv+1;socle=max_degree=0;input_count=0;
    spaces=NULL;input=NULL;degrees=prefix_degree=hilbert_dims=NULL;
    prime=fq_nmod_ctx_modulus(polys[0].ctx)->mod.n;
    nmod_init(&mod,prime);grouped_kernel=1;
    /* Enforce representation limits only; there is no memory/work budget. */
    slong ds[n],sum=0,terms=0;
    char spec[n*(3*sizeof(slong)+2)];size_t used=0;
    stats->reason="general-degree quotient requires positive degrees fitting machine address capacity";
    for(slong f=0;f<n;f++) {
        if(polys[f].nvars!=nv || polys[f].npars!=1)return 0;
        ds[f]=0;
        for(slong t=0;t<dr_mpoly_length(polys+f);t++) {
            DR_MPOLY_TERM(term,polys+f,t);
            slong degree=term.par_exp[0];
            if(degree<0)return 0;
            for(slong v=0;v<nv;v++) {
                if(term.var_exp[v]<0 || term.var_exp[v]>LONG_MAX-degree)return 0;
                degree+=term.var_exp[v];
            }
            ds[f]=FLINT_MAX(ds[f],degree);
        }
        if(!ds[f] || ds[f]>LONG_MAX-n-sum)return 0;
        sum+=ds[f]-1;
        if(dr_mpoly_length(polys+f)>LONG_MAX-terms)return 0;
        terms+=dr_mpoly_length(polys+f);
        used+=snprintf(spec+used,sizeof(spec)-used,"%s%ld",f?",":"",ds[f]);
    }
    stats->reason="field too small for interpolation and validation points";
    if(calculate_field_size(prime,fq_nmod_ctx_degree(polys[0].ctx))<((ulong)sum+5))return 0;
    int ok=0,saved_threads=flint_get_num_threads();
    const fq_nmod_ctx_struct *ctx=polys[0].ctx;
    dq_field_setup(ctx);
    dq_coeff *nodes=NULL,*values=NULL;
    slong sample_count=0,value_count=0;
    fq_nmod_poly_mat_t matrix;int matrix_live=0;
    fq_nmod_poly_t det;int det_live=0;
    fq_nmod_poly_struct *weights=NULL;
    fq_nmod_t scalar,check,node;
    fq_nmod_init(scalar,ctx);fq_nmod_init(check,ctx);fq_nmod_init(node,ctx);
    if(!admit_dimensions(spec)){stats->reason=admission_reason;goto done;}
    stats->reason="input storage exceeds machine address capacity";
    if((ulong)terms>SIZE_MAX/sizeof(*input) ||
       (ulong)terms>SIZE_MAX/(n*DQ_COEFF_BYTES))goto done;
    flint_set_num_threads(1);
    double t=seconds();
    init_spaces();
    /* Sparse affine input; shared specialization accumulates equal monomials. */
    input=flint_calloc(terms,sizeof(*input));
    for(slong f=0;f<n;f++)for(slong j=0;j<dr_mpoly_length(polys+f);j++) {
        DR_MPOLY_TERM(term,polys+f,j);
        input_term *out=input+input_count++;
        out->key=key_set(0,nv,term.par_exp[0]);
        for(slong v=0;v<nv;v++)out->key=key_set(out->key,v,term.var_exp[v]);
        out->coeff=dq_vec_new(n);
        fq_nmod_to_field_elem(out->coeff+f,term.coeff,&extension_field);
    }
    slong h=hilbert_dims[left_degree],samples=socle+3;
    stats->core_size=h;stats->samples=samples;
    sample_count=samples+2;value_count=h*h*sample_count;
    nodes=dq_vec_new(sample_count);values=dq_vec_new(value_count);
    /* Base-p polynomial digits enumerate distinct field elements, even in
     * characteristic two. Integer embedding would only enumerate F_p. */
    for(slong i=0;i<sample_count;i++) {
        ulong digits=i;fq_nmod_zero(node,ctx);
        for(slong j=0;digits;j++,digits/=prime)nmod_poly_set_coeff_ui(node,j,digits%prime);
        fq_nmod_to_field_elem(nodes+i,node,&extension_field);
    }
    info("\nStep 1: Prepare graded quotient bases for the Dixon Bezoutian\n");
    info("  Extension field; %s; grouped quotient recurrence; Core: %ld x %ld\n",
         field_ctx_get_description(&extension_field),h,h);
    /* Establish the basis once. Every sample, including those handled by a
     * freshly cloned worker, must match it rather than select its own basis. */
    dq_coeff **initial=new_coefficients();specialize(initial,nodes[0]);
    int basis_ok=normal_forms(initial,1);clear_coefficients(initial);
    stats->reason="non-generic Macaulay rank or changing quotient basis";
    if(!basis_ok)goto done;
    stats->leading=seconds()-t;dixon_maybe_print_step_time("Step 1",stats->leading);
    slong workers=sample_workers(requested_workers,samples+2);
    info("\nStep 2: Construct Dixon Bezoutian cores by grouped quotient recurrence\n");
    info("  Evaluation points: %ld; validation points: 2; requested threads: %ld\n",
         samples,FLINT_MAX(1,requested_workers));
    t=seconds();
    slong used_workers=1;
    int samples_ok=sample_cores(values,nodes,samples+2,h,workers,&used_workers);
    info("  sample workers: %ld\n",used_workers);
    if(!samples_ok)goto done;
    stats->evaluation=seconds()-t;dixon_maybe_print_step_time("Step 2",stats->evaluation);
    info("\nStep 3: Interpolate and validate the Dixon Bezoutian core\n");
    t=seconds();fq_nmod_poly_mat_init(matrix,h,h,ctx);matrix_live=1;
    weights=interpolation_weights(nodes,samples,ctx);
    stats->reason="core interpolation or weighted degree check failed";
    for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
        fq_nmod_poly_struct *p=fq_nmod_poly_mat_entry(matrix,i,j);
        for(slong sample=0;sample<samples;sample++) {
            field_elem_to_fq_nmod(scalar,values+(sample*h+i)*h+j,&extension_field);
            fq_nmod_poly_scalar_addmul_fq_nmod(p,weights+sample,scalar,ctx);
        }
        slong wi=left_degree-exponent(spaces[left_degree].keys[spaces[left_degree].free[i]],n-1);
        slong wj=right_degree-exponent(spaces[right_degree].keys[spaces[right_degree].free[j]],n-1);
        if(fq_nmod_poly_degree(p,ctx)>socle+1-wi-wj)goto done;
        for(slong sample=samples;sample<samples+2;sample++) {
            field_elem_to_fq_nmod(node,nodes+sample,&extension_field);
            field_elem_to_fq_nmod(scalar,values+(sample*h+i)*h+j,&extension_field);
            fq_nmod_poly_evaluate_fq_nmod(check,p,node,ctx);
            if(!fq_nmod_equal(check,scalar,ctx))goto done;
        }
    }
    stats->interpolation=seconds()-t;dixon_maybe_print_step_time("Step 3",stats->interpolation);
    info("\nStep 4: Compute Dixon resultant from the quotient core determinant\n");
    t=seconds();fq_nmod_poly_init(det,ctx);det_live=1;
    stats->reason="core determinant or expected resultant degree check failed";
    fq_nmod_poly_mat_det_iter(det,matrix,ctx);
    if(fq_nmod_poly_degree(det,ctx)!=expected_degree)goto done;
    fq_nmod_poly_make_monic(det,det,ctx);
    stats->determinant=seconds()-t;dixon_maybe_print_step_time("Step 4",stats->determinant);
    info("  Resultant degree: %ld; normalize to monic.\n",fq_nmod_poly_degree(det,ctx));
    dr_mpoly_init(result,0,1,ctx);
    for(slong i=0;i<fq_nmod_poly_length(det,ctx);i++) {
        fq_nmod_poly_get_coeff(scalar,det,i,ctx);
        if(!fq_nmod_is_zero(scalar,ctx))dr_mpoly_add_term_fast(result,NULL,&i,scalar);
    }
    dr_mpoly_normalize(result);stats->reason=NULL;ok=1;
done:
    if(det_live)fq_nmod_poly_clear(det,ctx);
    if(matrix_live)fq_nmod_poly_mat_clear(matrix,ctx);
    if(weights) {
        for(slong i=0;i<sample_count-2;i++)fq_nmod_poly_clear(weights+i,ctx);
        flint_free(weights);
    }
    dq_vec_free(nodes,sample_count);dq_vec_free(values,value_count);
    fq_nmod_clear(scalar,ctx);fq_nmod_clear(check,ctx);fq_nmod_clear(node,ctx);
    if(spaces)for(slong d=0;d<=limit;d++) {
        flint_free(spaces[d].keys);flint_free(spaces[d].free);dq_mat_clear(spaces[d].reduce);
    }
    flint_free(spaces);key_tables_clear();
    for(slong i=0;i<input_count;i++)dq_vec_free(input[i].coeff,n);flint_free(input);
    flint_free(degrees);flint_free(prefix_degree);flint_free(hilbert_dims);
    field_ctx_clear(&extension_field);
    flint_set_num_threads(saved_threads);stats->total=seconds()-begin;
    return ok;
}

int dr_extension_quotient_resultant(unified_mpoly_struct *result,
    const unified_mpoly_struct *polys,slong nv,slong np,slong requested_workers,mq_quotient_stats *stats)
{
    int ok;mq_quotient_stats local;if(!stats)stats=&local;
    #pragma omp critical(drsolve_quotient_context)
    {ok=compute(result,polys,nv,np,requested_workers,stats);}
    return ok;
}

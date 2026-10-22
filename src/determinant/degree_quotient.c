/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_quotient.h"
#include "dixon_flint.h"
#include <setjmp.h>
#include <stdarg.h>
#ifdef _OPENMP
#include <omp.h>
#endif
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
        nmod_mat_init(dst->reduce,dst->size,dst->h,prime);
    }
    return out;
}
static void clear_worker_copy(degree_worker *worker)
{
    for(slong d=0;d<=limit;d++) {
        flint_free(worker->grades[d].free);nmod_mat_clear(worker->grades[d].reduce);
    }
    flint_free(worker->grades);
}
static ulong **new_coefficients(void)
{
    ulong **q=flint_calloc(n,sizeof(*q));
    for(slong f=0;f<n;f++)q[f]=flint_malloc(spaces[degrees[f]].size*sizeof(ulong));
    return q;
}
static void clear_coefficients(ulong **q)
{
    for(slong f=0;f<n;f++)flint_free(q[f]);flint_free(q);
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
static int sample_cores(ulong *values,slong points,slong h,slong workers,slong *used)
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
        ulong **q=new_coefficients();
        nmod_mat_t core;nmod_mat_init(core,h,h,prime);
        #pragma omp for schedule(dynamic,1)
        for(slong sample=0;sample<points;sample++) {
            specialize(q,sample);
            if(!normal_forms(q,0)){valid=0;continue;}
            entry *entries=difference_matrix(q);bezout_core(core,entries,1);
            for(slong k=0;k<n*n;k++)flint_free(entries[k].terms);flint_free(entries);
            for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)
                values[(sample*h+i)*h+j]=nmod_mat_entry(core,i,j);
        }
        nmod_mat_clear(core);clear_coefficients(q);active_worker=saved;
    }
    for(slong w=1;w<workers;w++)clear_worker_copy(contexts+w);
    flint_free(contexts);return valid;
}
static int compute(unified_mpoly_struct *result,const unified_mpoly_struct *polys,
                   slong nv,slong np,slong requested_workers,mq_quotient_stats *stats)
{
    double begin=seconds();memset(stats,0,sizeof(*stats));
    stats->reason="general-degree quotient requires at least two equations over a prime field and one parameter";
    if(np!=1 || nv<1 || polys[0].field_id!=FIELD_ID_NMOD ||
       fq_nmod_ctx_degree(polys[0].ctx)!=1 ||
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
    if(prime<=((ulong)sum+4))return 0;
    int ok=0,saved_threads=flint_get_num_threads();
    ulong *xs=NULL,*values=NULL,*column=NULL;
    nmod_poly_mat_t matrix;int matrix_live=0;
    nmod_poly_t det;int det_live=0;
    if(!admit_dimensions(spec)){stats->reason=admission_reason;goto done;}
    stats->reason="input storage exceeds machine address capacity";
    if((ulong)terms>SIZE_MAX/sizeof(*input) ||
       (ulong)terms>SIZE_MAX/(n*sizeof(ulong)))goto done;
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
        out->coeff=flint_calloc(n,sizeof(ulong));
        out->coeff[f]=nmod_poly_get_coeff_ui(term.coeff,0);
    }
    slong h=hilbert_dims[left_degree],samples=socle+3;
    stats->core_size=h;stats->samples=samples;
    xs=flint_malloc(samples*sizeof(ulong));
    column=flint_malloc(samples*sizeof(ulong));
    values=flint_malloc(h*h*(samples+2)*sizeof(ulong));
    info("\nStep 1: Prepare graded quotient bases for the Dixon Bezoutian\n");
    info("  General degrees; prime field; grouped quotient recurrence; Core: %ld x %ld\n",h,h);
    /* Establish the basis once. Every sample, including those handled by a
     * freshly cloned worker, must match it rather than select its own basis. */
    ulong **initial=new_coefficients();specialize(initial,0);
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
    int samples_ok=sample_cores(values,samples+2,h,workers,&used_workers);
    info("  sample workers: %ld\n",used_workers);
    if(!samples_ok)goto done;
    for(slong sample=0;sample<samples;sample++)xs[sample]=sample;
    stats->evaluation=seconds()-t;dixon_maybe_print_step_time("Step 2",stats->evaluation);
    info("\nStep 3: Interpolate and validate the Dixon Bezoutian core\n");
    t=seconds();nmod_poly_mat_init(matrix,h,h,prime);matrix_live=1;
    stats->reason="core interpolation or weighted degree check failed";
    for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
        nmod_poly_struct *p=nmod_poly_mat_entry(matrix,i,j);
        for(slong sample=0;sample<samples;sample++)column[sample]=values[(sample*h+i)*h+j];
        nmod_poly_interpolate_nmod_vec(p,xs,column,samples);
        slong wi=left_degree-exponent(spaces[left_degree].keys[spaces[left_degree].free[i]],n-1);
        slong wj=right_degree-exponent(spaces[right_degree].keys[spaces[right_degree].free[j]],n-1);
        if(nmod_poly_degree(p)>socle+1-wi-wj)goto done;
        for(slong sample=samples;sample<samples+2;sample++)
            if(nmod_poly_evaluate_nmod(p,sample)!=values[(sample*h+i)*h+j])goto done;
    }
    stats->interpolation=seconds()-t;dixon_maybe_print_step_time("Step 3",stats->interpolation);
    info("\nStep 4: Compute Dixon resultant from the quotient core determinant\n");
    t=seconds();nmod_poly_init(det,prime);det_live=1;
    stats->reason="core determinant or expected resultant degree check failed";
#ifdef HAVE_PML
    if(!nmod_poly_mat_det_hnf(det,matrix))goto done;
#else
    nmod_poly_mat_det(det,matrix);
#endif
    if(nmod_poly_degree(det)!=expected_degree)goto done;
    nmod_poly_make_monic(det,det);
    stats->determinant=seconds()-t;dixon_maybe_print_step_time("Step 4",stats->determinant);
    info("  Resultant degree: %ld; normalize to monic.\n",nmod_poly_degree(det));
    dr_mpoly_init(result,0,1,polys[0].ctx);
    for(slong i=0;i<det->length;i++)if(det->coeffs[i]) {
        ulong e=i;
        nmod_mpoly_push_term_ui_ui(GET_NMOD_POLY(result),det->coeffs[i],&e,GET_NMOD_CTX(result->ctx_ptr));
    }
    nmod_mpoly_sort_terms(GET_NMOD_POLY(result),GET_NMOD_CTX(result->ctx_ptr));
    dr_mpoly_normalize(result);stats->reason=NULL;ok=1;
done:
    if(det_live)nmod_poly_clear(det);
    if(matrix_live)nmod_poly_mat_clear(matrix);
    flint_free(xs);flint_free(values);flint_free(column);
    if(spaces)for(slong d=0;d<=limit;d++) {
        flint_free(spaces[d].keys);flint_free(spaces[d].free);nmod_mat_clear(spaces[d].reduce);
    }
    flint_free(spaces);key_tables_clear();
    for(slong i=0;i<input_count;i++)flint_free(input[i].coeff);flint_free(input);
    flint_free(degrees);flint_free(prefix_degree);flint_free(hilbert_dims);
    flint_set_num_threads(saved_threads);stats->total=seconds()-begin;
    return ok;
}

int dr_degree_quotient_resultant(unified_mpoly_struct *result,
    const unified_mpoly_struct *polys,slong nv,slong np,slong requested_workers,mq_quotient_stats *stats)
{
    if(fq_nmod_ctx_degree(polys[0].ctx)>1)
        return dr_extension_quotient_resultant(result,polys,nv,np,requested_workers,stats);

    int ok;mq_quotient_stats local;if(!stats)stats=&local;
    #pragma omp critical(drsolve_quotient_context)
    {ok=compute(result,polys,nv,np,requested_workers,stats);}
    return ok;
}

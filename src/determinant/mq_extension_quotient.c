/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_quotient.h"
#include "fq_poly_mat_det.h"
#define MQ_QUOTIENT_EXTENSION
#include "degree_quotient_extension_field.h"
#include "mq_quotient_impl.h"
#include "quotient_extension_interpolation.h"
#include "dixon_flint.h"
#include <stdarg.h>

static void quotient_info(const char *fmt, ...)
{
    if (g_dixon_verbose_level < 1) return;
    va_list args;
    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    fflush(stdout);
}

static int quotient_compute(unified_mpoly_struct *result,
    const unified_mpoly_struct *polys,slong nv,slong np,slong workers,
    mq_quotient_stats *stats)
{
    double begin=now();memset(stats,0,sizeof(*stats));
    stats->reason="requires a square quadratic system over an extension field with one parameter";
    if(np!=1 || nv<1 || fq_nmod_ctx_degree(polys[0].ctx)<=1 || g_field_equation_reduction || g_field_equation_final_only)return 0;
    stats->reason="monomial keys or subset masks exceed machine-word capacity";
    if(nv>=FLINT_BITS-2 || !quotient_dimensions(nv+1))return -1;
    threads=FLINT_MIN(N+4,FLINT_MAX(1,workers));
    prime=fq_nmod_ctx_modulus(polys[0].ctx)->mod.n;
    stats->reason="field too small for distinct interpolation points";
    if(calculate_field_size(prime,fq_nmod_ctx_degree(polys[0].ctx))<(ulong)(N+4))return 0;
    nmod_init(&mod,prime);normal_time=map_time=dp_time=0;
    recurrence_kernel=4;verify_recurrence=0;frames_enabled=0;
    const fq_nmod_ctx_struct *ctx=polys[0].ctx;
    dq_field_setup(ctx);
    slong w=N+1;dq_coeff *input=dq_vec_new(N*w*w);
    int ok=0;
    stats->reason="input is not a quadratic system";
    for(slong f=0;f<N;f++) {
        if(polys[f].nvars!=nv || polys[f].npars!=1)goto input_done;
        for(slong t=0;t<dr_mpoly_length(polys+f);t++) {
            DR_MPOLY_TERM(term,polys+f,t);
            slong indices[2]={N,N},degree=0;
            for(slong v=0;v<N;v++) {
                slong e=v<nv?term.var_exp[v]:term.par_exp[0];
                if(e<0 || e>2 || degree+e>2)goto input_done;
                while(e--)indices[degree++]=v;
            }
            fq_nmod_to_field_elem(input+(f*w+indices[0])*w+indices[1],term.coeff,&extension_field);
        }
    }
    {
        quotient_info("  Extension field; %s; shared leading-frame MQ kernel\n",field_ctx_get_description(&extension_field));
        if(mq_byte_enabled())quotient_info("  Compact byte minors; SIMD field accumulation; layer-local maps\n");
        quotient_info("\nStep 1: Prepare quotient bases for the Dixon Bezoutian\n");
        quotient_info("  Quadratic equations: %ld; eliminated variables: %ld; parameter: 1\n",N,nv);
        quotient_info("  Factor leading Macaulay blocks through degree %ld; check quotient dimensions.\n",L);
        int saved_threads=flint_get_num_threads();flint_set_num_threads(1);
        dq_coeff zero;dq_init(zero);field_set_zero(&zero,DQ_FIELD);
        dq_coeff *q=dq_vec_new(N*N*N);specialize(q,input,zero);dq_clear(zero);
        /* Shared factorization may use workers; sample teams use serial FLINT. */
        flint_set_num_threads(FLINT_MIN(threads,8));
        double t=now();int admitted=leading_init(q);stats->leading=now()-t;
        flint_set_num_threads(1);
        stats->reason="leading Macaulay blocks do not have the required ranks";
        if(!admitted){leading_clear();dq_vec_free(q,N*N*N);flint_set_num_threads(saved_threads);goto input_done;}
        frames_enabled=1;
        slong h=choose(N,M),samples=N+2;
        stats->core_size=h;stats->samples=samples;
        uint64_t *basis=flint_malloc(h*sizeof(uint64_t)),*basis_y=flint_malloc(h*sizeof(uint64_t));
        leading_basis(basis,M);leading_basis(basis_y,K);
        quotient_info("  Row basis: degree %ld, dimension C(%ld,%ld) = %ld\n",M,N,M,h);
        quotient_info("  Column basis: degree %ld, dimension C(%ld,%ld) = %ld\n",K,N,K,choose(N,K));
        if(g_dixon_verbose_level>=2)for(slong d=2;d<=L;d++)
            quotient_info("  Degree %ld leading block: %ld x %ld, rank %ld\n",
                d,N*frames[d].lower_s,frames[d].sx,frames[d].rank);
        dixon_maybe_print_step_time("Step 1",stats->leading);
        dq_vec_free(q,N*N*N);
        slong points=samples+2;
        dq_coeff *xs=dq_vec_new(points),*ys=dq_vec_new(h*h*points);
        fq_nmod_t node,value,check;fq_nmod_init(node,ctx);fq_nmod_init(value,ctx);fq_nmod_init(check,ctx);
        for(slong i=0;i<points;i++) {
            ulong digits=i;fq_nmod_zero(node,ctx);
            for(slong j=0;digits;j++,digits/=prime)nmod_poly_set_coeff_ui(node,j,digits%prime);
            fq_nmod_to_field_elem(xs+i,node,&extension_field);
        }
        quotient_info("\nStep 2: Construct Dixon Bezoutian cores in quotient space\n");
        quotient_info("  Reuse the decomposition for normal forms and multiplication maps at each point.\n");
        quotient_info("  Build each reduced Dixon Bezoutian by minor recursion in the two quotient bases.\n");
        quotient_info("  Core: %ld x %ld, evaluation points: %ld; validation points: 2\n",h,h,samples);
        t=now();
        #pragma omp parallel for schedule(dynamic) num_threads(threads) if(threads>1)
        for(slong sample=0;sample<points;sample++) {
            dq_mat_t core;dq_mat_init(core,h,h,prime);
            evaluate(core,input,xs[sample],basis,basis_y,0);
            for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)dq_set(ys[(i*h+j)*points+sample],dq_mat_entry(core,i,j));
            dq_mat_clear(core);
        }
        stats->evaluation=now()-t;
        dixon_maybe_print_step_time("Step 2",stats->evaluation);
        quotient_info("\nStep 3: Interpolate the Dixon Bezoutian core matrix\n");
        quotient_info("  Interpolate %ld x %ld entries from %ld points; check entry degree bounds.\n",h,h,samples);
        t=now();
        fq_nmod_poly_mat_t polynomial;fq_nmod_poly_mat_init(polynomial,h,h,ctx);
        fq_nmod_poly_t det;fq_nmod_poly_init(det,ctx);ok=1;
        fq_nmod_poly_struct *weights=interpolation_weights(xs,samples,ctx);
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
            fq_nmod_poly_struct *entry=fq_nmod_poly_mat_entry(polynomial,i,j);
            for(slong sample=0;sample<samples;sample++) {
                field_elem_to_fq_nmod(value,ys+(i*h+j)*points+sample,&extension_field);
                fq_nmod_poly_scalar_addmul_fq_nmod(entry,weights+sample,value,ctx);
            }
            for(slong sample=samples;sample<points;sample++) {
                field_elem_to_fq_nmod(node,xs+sample,&extension_field);
                field_elem_to_fq_nmod(value,ys+(i*h+j)*points+sample,&extension_field);
                fq_nmod_poly_evaluate_fq_nmod(check,entry,node,ctx);
                if(!fq_nmod_equal(check,value,ctx))ok=0;
            }
            slong wi=M-(basis[i]>>(key_bits*(N-1))),wj=K-(basis_y[j]>>(key_bits*(N-1)));
            if(fq_nmod_poly_degree(entry,ctx)>N+1-wi-wj)ok=0;
        }
        stats->interpolation=now()-t;
        dixon_maybe_print_step_time("Step 3",stats->interpolation);
        int determinant_attempted=ok;
        if(ok) {
            quotient_info("\nStep 4: Compute Dixon resultant from the core determinant\n");
            quotient_info("  Determinant method: unified extension-field backend\n");
            quotient_info("  Expected resultant degree: %ld\n",1L<<N);
        }
        t=now();
        stats->reason="core degree or determinant check failed";
        if(ok)fq_nmod_poly_mat_det_iter(det,polynomial,ctx);
        if(ok)ok=fq_nmod_poly_degree(det,ctx)==(1L<<N);
        stats->determinant=now()-t;
        if(ok)quotient_info("  Resultant degree: %ld; normalize to monic.\n",fq_nmod_poly_degree(det,ctx));
        if(determinant_attempted)
            dixon_maybe_print_step_time("Step 4",stats->determinant);
        if(ok) {
            fq_nmod_poly_make_monic(det,det,ctx);dr_mpoly_init(result,0,1,ctx);
            for(slong i=0;i<fq_nmod_poly_length(det,ctx);i++) {
                fq_nmod_poly_get_coeff(value,det,i,ctx);
                if(!fq_nmod_is_zero(value,ctx))dr_mpoly_add_term_fast(result,NULL,&i,value);
            }
            dr_mpoly_normalize(result);stats->reason=NULL;
        }
        for(slong i=0;i<samples;i++)fq_nmod_poly_clear(weights+i,ctx);flint_free(weights);
        fq_nmod_poly_clear(det,ctx);fq_nmod_poly_mat_clear(polynomial,ctx);
        fq_nmod_clear(node,ctx);fq_nmod_clear(value,ctx);fq_nmod_clear(check,ctx);
        flint_free(basis);flint_free(basis_y);dq_vec_free(xs,points);dq_vec_free(ys,h*h*points);
        leading_clear();frames_enabled=0;flint_set_num_threads(saved_threads);
    }
input_done:
    dq_vec_free(input,N*w*w);field_ctx_clear(&extension_field);stats->total=now()-begin;return ok;
}

int dr_mq_extension_quotient_resultant(unified_mpoly_struct *result,
    const unified_mpoly_struct *polys,slong nv,slong np,slong workers,mq_quotient_stats *stats)
{
    int ok;mq_quotient_stats local;if(!stats)stats=&local;
    #pragma omp critical(drsolve_quotient_context)
    {ok=quotient_compute(result,polys,nv,np,workers,stats);}
    /* A packed key limit selects the dynamic-key kernel, never rejects a
     * system that the general representation can store. */
    if(ok<0)return dr_extension_quotient_resultant(result,polys,nv,np,workers,stats);
    return ok;
}

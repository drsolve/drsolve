/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_quotient.h"
#include "mq_quotient_impl.h"
#include "dixon_flint.h"
#include <stdarg.h>
#include <time.h>

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
    stats->reason="requires a square quadratic system over a prime field with one parameter";
    if(np!=1 || nv<1 || fq_nmod_ctx_degree(polys[0].ctx)!=1 ||
       polys[0].field_id!=FIELD_ID_NMOD || g_field_equation_reduction || g_field_equation_final_only)return 0;
    stats->reason="monomial keys or subset masks exceed machine-word capacity";
    if(nv>=FLINT_BITS-2 || !quotient_dimensions(nv+1))return 0;
    threads=FLINT_MAX(1,workers);
    prime=fq_nmod_ctx_modulus(polys[0].ctx)->mod.n;
    stats->reason="field too small for distinct interpolation points";
    if(prime<=(ulong)(N+1))return 0;
    nmod_init(&mod,prime);normal_time=map_time=dp_time=0;
    recurrence_kernel=packed_admitted()?5:1;verify_recurrence=0;frames_enabled=0;
    slong w=N+1;ulong *input=flint_calloc(N*w*w,sizeof(ulong));
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
            input[(f*w+indices[0])*w+indices[1]]=nmod_poly_get_coeff_ui(term.coeff,0);
        }
    }
    {
        quotient_info("\nStep 1: Prepare quotient bases for the Dixon Bezoutian\n");
        quotient_info("  Quadratic equations: %ld; eliminated variables: %ld; parameter: 1\n",N,nv);
        quotient_info("  Factor leading Macaulay blocks through degree %ld; check quotient dimensions.\n",L);
        int saved_threads=flint_get_num_threads();flint_set_num_threads(1);
        ulong *q=flint_calloc(N*N*N,sizeof(ulong));specialize(q,input,0);
        /* Shared factorization may use workers; sample teams use serial FLINT. */
        flint_set_num_threads(FLINT_MIN(threads,8));
        clock_t cpu_start=clock();
        double t=now();int admitted=leading_init(q);stats->leading=now()-t;
        double cpu_elapsed=(double)(clock()-cpu_start)/CLOCKS_PER_SEC;
        flint_set_num_threads(1);
        stats->reason="leading Macaulay blocks do not have the required ranks";
        if(!admitted){leading_clear();flint_free(q);flint_set_num_threads(saved_threads);goto input_done;}
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
        dixon_maybe_print_parallel_step_time("Step 1",cpu_elapsed,stats->leading);
        flint_free(q);
        ulong *xs=flint_malloc(samples*sizeof(ulong)),*ys=flint_malloc(h*h*samples*sizeof(ulong));
        quotient_info("\nStep 2: Construct Dixon Bezoutian cores in quotient space\n");
        quotient_info("  Reuse the decomposition for normal forms and multiplication maps at each point.\n");
        quotient_info("  Build each reduced Dixon Bezoutian by minor recursion in the two quotient bases.\n");
        quotient_info("  Core: %ld x %ld, evaluation points: %ld\n",h,h,samples);
        cpu_start=clock();t=now();
        #pragma omp parallel for schedule(dynamic) num_threads(threads) if(threads>1)
        for(slong sample=0;sample<samples;sample++) {
            nmod_mat_t core;nmod_mat_init(core,h,h,prime);
            evaluate(core,input,sample,basis,basis_y,0);xs[sample]=sample;
            for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)ys[(i*h+j)*samples+sample]=nmod_mat_entry(core,i,j);
            nmod_mat_clear(core);
        }
        stats->evaluation=now()-t;
        cpu_elapsed=(double)(clock()-cpu_start)/CLOCKS_PER_SEC;
        dixon_maybe_print_parallel_step_time("Step 2",cpu_elapsed,stats->evaluation);
        quotient_info("\nStep 3: Interpolate the Dixon Bezoutian core matrix\n");
        quotient_info("  Interpolate %ld x %ld entries from %ld points; check entry degree bounds.\n",h,h,samples);
        cpu_start=clock();t=now();
        nmod_poly_mat_t polynomial;nmod_poly_mat_init(polynomial,h,h,prime);
        nmod_poly_t det;nmod_poly_init(det,prime);ok=1;
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
            nmod_poly_struct *entry=nmod_poly_mat_entry(polynomial,i,j);
            nmod_poly_interpolate_nmod_vec(entry,xs,ys+(i*h+j)*samples,samples);
            slong wi=M-(basis[i]>>(key_bits*(N-1))),wj=K-(basis_y[j]>>(key_bits*(N-1)));
            if(nmod_poly_degree(entry)>N+1-wi-wj)ok=0;
        }
        stats->interpolation=now()-t;
        cpu_elapsed=(double)(clock()-cpu_start)/CLOCKS_PER_SEC;
        dixon_maybe_print_parallel_step_time("Step 3",cpu_elapsed,stats->interpolation);
        int determinant_attempted=ok;
        if(ok) {
            quotient_info("\nStep 4: Compute Dixon resultant from the core determinant\n");
            quotient_info("  Determinant method: fast HNF\n");
            quotient_info("  Expected resultant degree: %ld\n",1L<<N);
        }
        cpu_start=clock();t=now();
        stats->reason="core degree or determinant check failed";
        if(ok)ok=nmod_poly_mat_det_hnf(det,polynomial);
        if(ok)ok=nmod_poly_degree(det)==(1L<<N);
        stats->determinant=now()-t;
        cpu_elapsed=(double)(clock()-cpu_start)/CLOCKS_PER_SEC;
        if(ok)quotient_info("  Resultant degree: %ld; normalize to monic.\n",nmod_poly_degree(det));
        if(determinant_attempted)
            dixon_maybe_print_parallel_step_time("Step 4",cpu_elapsed,stats->determinant);
        if(ok) {
            nmod_poly_make_monic(det,det);dr_mpoly_init(result,0,1,polys[0].ctx);
            for(slong i=0;i<det->length;i++)if(det->coeffs[i]) {
                ulong e=(ulong)i;
                nmod_mpoly_push_term_ui_ui(GET_NMOD_POLY(result),det->coeffs[i],&e,GET_NMOD_CTX(result->ctx_ptr));
            }
            nmod_mpoly_sort_terms(GET_NMOD_POLY(result),GET_NMOD_CTX(result->ctx_ptr));
            dr_mpoly_normalize(result);stats->reason=NULL;
        }
        nmod_poly_clear(det);nmod_poly_mat_clear(polynomial);
        flint_free(basis);flint_free(basis_y);flint_free(xs);flint_free(ys);
        leading_clear();frames_enabled=0;flint_set_num_threads(saved_threads);
    }
input_done:
    flint_free(input);stats->total=now()-begin;return ok;
}

int dr_mq_quotient_resultant(unified_mpoly_struct *result,
    const unified_mpoly_struct *polys,slong nv,slong np,slong workers,mq_quotient_stats *stats)
{
    if(fq_nmod_ctx_degree(polys[0].ctx)>1)
        return dr_mq_extension_quotient_resultant(result,polys,nv,np,workers,stats);

    int ok;mq_quotient_stats local;if(!stats)stats=&local;
    /* The recurrence's shared context is immutable inside the sample team. */
    #pragma omp critical(drsolve_quotient_context)
    {ok=quotient_compute(result,polys,nv,np,workers,stats);}
    return ok;
}

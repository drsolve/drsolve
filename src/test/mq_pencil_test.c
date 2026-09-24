/* SPDX-License-Identifier: GPL-2.0-or-later */
#define main mq_filter_regression_main
#include "dixon_mq_filter_test.c"
#undef main
#include "mq_pencil_det.h"
#include <time.h>

static double pencil_test_time(void)
{
    struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;
}

int main(int argc,char **argv)
{
    slong n=argc>1?atol(argv[1]):4;
    ulong prime=argc>2?strtoul(argv[2],NULL,10):65537;
    int threads=argc>3?atoi(argv[3]):1;
    if(n<2 || n>8 || !n_is_prime(prime) || threads<1 || threads>64)return 2;
    omp_set_num_threads(threads);flint_set_num_threads(threads);
    g_dixon_verbose_level=0;g_dixon_det_cache_limit=100000;
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,132,941);
    fq_nmod_ctx_t fq;fq_nmod_ctx_init_ui(fq,prime,1,"a");
    fq_mvpoly_t *p=random_mq(n-1,fq,rng),**m,**a,full,actual;
    build_fq_cancellation_matrix_mvpoly(&m,p,n-1,1);
    perform_fq_matrix_row_operations_mvpoly(&a,&m,n-1,1);
    double start=pencil_test_time();
    compute_fq_cancel_matrix_det(&full,a,n-1,1,DET_METHOD_RECURSIVE);
    double dp=pencil_test_time()-start;
    nmod_mpoly_ctx_t ctx;nmod_mpoly_ctx_init(ctx,2*n-1,ORD_LEX,prime);
    nmod_mpoly_struct *saved=flint_malloc(n*n*sizeof(*saved));
    for(slong i=0;i<n;i++)for(slong j=0;j<n;j++) {
        nmod_mpoly_init(saved+i*n+j,ctx);
        fq_mvpoly_to_nmod_mpoly(saved+i*n+j,&a[i][j],ctx);
    }
    mq_pencil_stats stats;
    if(prime<=(ulong)n-1) {
        fq_mvpoly_copy(&actual,&full);fq_monomial_t *ptr=actual.terms;
        assert(!compute_fq_det_mq_pencil(&actual,a,n,&stats));
        assert(actual.terms==ptr && actual.nterms==full.nterms);
        fq_mvpoly_clear(&actual);
        printf("{\"n\":%ld,\"q\":%lu,\"characteristic_rejected\":true,\"preserved\":true}\n",n,prime);
    } else {
        int ok=compute_fq_det_mq_pencil(&actual,a,n,&stats);
        if(!ok) {
            assert(strstr(stats.reason,"rank deficient"));
            printf("{\"n\":%ld,\"q\":%lu,\"rank_rejected\":true}\n",n,prime);
        } else {
            nmod_mpoly_t x,y;nmod_mpoly_init(x,ctx);nmod_mpoly_init(y,ctx);
            fq_mvpoly_to_nmod_mpoly(x,&actual,ctx);fq_mvpoly_to_nmod_mpoly(y,&full,ctx);
            assert(nmod_mpoly_equal(x,y,ctx));
            for(slong i=0;i<n;i++)for(slong j=0;j<n;j++) {
                fq_mvpoly_to_nmod_mpoly(x,&a[i][j],ctx);
                assert(nmod_mpoly_equal(x,saved+i*n+j,ctx));
            }
            printf("{\"n\":%ld,\"q\":%lu,\"threads\":%ld,\"terms\":%ld,\"full_dp\":%.9f,\"normalization\":%.9f,\"recurrence\":%.9f,\"assembly\":%.9f,\"total\":%.9f,\"peak_terms\":%ld,\"equal\":true}\n",
                n,prime,stats.threads,actual.nterms,dp,stats.normalization,stats.recurrence,
                stats.assembly,stats.total,stats.peak_terms);
            /* Non-downward-closed, rectangular targets: intermediate divisors
             * must survive even though they are not final requested terms. */
            slong rows[24]={0},cols[16]={0};
            rows[0]=2;rows[n-1+1]=1;rows[2*(n-1)+n-2]=2;
            cols[0]=1;cols[2*(n-1)-1]=2;
            fq_mvpoly_t projected,expected;
            assert(compute_fq_det_mq_pencil_projected(&projected,a,n,rows,3,cols,2,&stats));
            assert(fq_mq_project_full(&expected,&full,rows,3,cols,2));
            fq_mvpoly_to_nmod_mpoly(x,&projected,ctx);fq_mvpoly_to_nmod_mpoly(y,&expected,ctx);
            assert(nmod_mpoly_equal(x,y,ctx));
            fq_mvpoly_clear(&projected);fq_mvpoly_clear(&expected);
            if(n>=5) {
                /* Larger irregular ideals exercise packed closure lookups. */
                slong wide_rows[128]={0},wide_cols[128]={0};
                for(slong i=0;i<16;i++) {
                    for(slong d=0;d<i%n;d++)wide_rows[i*(n-1)+n_randint(rng,n-1)]++;
                    for(slong d=0;d<(3*i+1)%n;d++)wide_cols[i*(n-1)+n_randint(rng,n-1)]++;
                }
                assert(compute_fq_det_mq_pencil_projected(&projected,a,n,wide_rows,16,wide_cols,16,&stats));
                assert(fq_mq_project_full(&expected,&full,wide_rows,16,wide_cols,16));
                fq_mvpoly_to_nmod_mpoly(x,&projected,ctx);fq_mvpoly_to_nmod_mpoly(y,&expected,ctx);
                assert(nmod_mpoly_equal(x,y,ctx));
                fq_mvpoly_clear(&projected);fq_mvpoly_clear(&expected);
            }
            rows[0]=-1;
            fq_monomial_t *preserved=actual.terms;
            assert(!compute_fq_det_mq_pencil_projected(&actual,a,n,rows,3,cols,2,&stats));
            assert(actual.terms==preserved);
            assert(!compute_fq_det_mq_pencil_projected(&actual,a,n,rows,0,cols,2,&stats));
            assert(actual.terms==preserved && stats.reason);
            fq_mvpoly_clear(&actual);
            nmod_mpoly_clear(x,ctx);nmod_mpoly_clear(y,ctx);
            /* Force C=[I|0], so only the LAST free-column choice works.
             * For even n this also forces an odd column permutation. */
            if(n<=5) {
                slong exps[64]={0},par=1;
                fq_nmod_t one;fq_nmod_init(one,fq);fq_nmod_one(one,fq);
                for(slong i=1;i<n;i++)for(slong j=0;j<n;j++) {
                    for(slong k=0;k<a[i][j].nterms;k++)
                        if(a[i][j].terms[k].par_exp[0])fq_nmod_zero(a[i][j].terms[k].coeff,fq);
                    if(j==i-1)fq_mvpoly_add_term_fast(&a[i][j],exps,&par,one);
                }
                fq_nmod_clear(one,fq);
                fq_mvpoly_t expected;
                compute_fq_cancel_matrix_det(&expected,a,n-1,1,DET_METHOD_RECURSIVE);
                assert(compute_fq_det_mq_pencil(&actual,a,n,&stats));
                nmod_mpoly_init(x,ctx);nmod_mpoly_init(y,ctx);
                fq_mvpoly_to_nmod_mpoly(x,&actual,ctx);fq_mvpoly_to_nmod_mpoly(y,&expected,ctx);
                assert(nmod_mpoly_equal(x,y,ctx));
                nmod_mpoly_clear(x,ctx);nmod_mpoly_clear(y,ctx);
                fq_mvpoly_clear(&expected);fq_mvpoly_clear(&actual);
            }
            /* Singular complete determinant, but the parameter normalization
             * remains valid: zero must be a successful exact result. */
            if(n<=5) {
                for(slong j=0;j<n;j++)for(slong k=0;k<a[0][j].nterms;k++)
                    fq_nmod_zero(a[0][j].terms[k].coeff,fq);
                assert(compute_fq_det_mq_pencil(&actual,a,n,&stats));
                assert(actual.nterms==0);fq_mvpoly_clear(&actual);
            }
        }
        /* Rank-deficient parameter coefficient matrix rejects without writing
         * to the caller's already initialized polynomial. */
        for(slong i=1;i<n;i++)for(slong j=0;j<n;j++)for(slong k=0;k<a[i][j].nterms;k++)
            if(a[i][j].terms[k].par_exp[0])fq_nmod_zero(a[i][j].terms[k].coeff,fq);
        fq_mvpoly_copy(&actual,&full);fq_monomial_t *ptr=actual.terms;
        assert(!compute_fq_det_mq_pencil(&actual,a,n,&stats));
        assert(strstr(stats.reason,"rank deficient"));
        assert(actual.terms==ptr && actual.nterms==full.nterms);
        fq_mvpoly_clear(&actual);
    }
    if(n==4 && prime==65537) {
        /* Exercise the memory guard without performing a large computation. */
        fq_mvpoly_t **large=flint_malloc(10*sizeof(*large));
        for(slong i=0;i<10;i++) {
            large[i]=flint_malloc(10*sizeof(**large));
            for(slong j=0;j<10;j++)fq_mvpoly_init(&large[i][j],18,1,fq);
        }
        fq_mvpoly_copy(&actual,&full);fq_monomial_t *ptr=actual.terms;
        assert(!compute_fq_det_mq_pencil(&actual,large,10,&stats));
        assert(strstr(stats.reason,"workspace") && actual.terms==ptr);
        fq_mvpoly_clear(&actual);
        for(slong i=0;i<10;i++) {
            for(slong j=0;j<10;j++)fq_mvpoly_clear(&large[i][j]);
            flint_free(large[i]);
        }
        flint_free(large);
    }
    for(slong i=0;i<n*n;i++)nmod_mpoly_clear(saved+i,ctx);
    flint_free(saved);nmod_mpoly_ctx_clear(ctx);
    fq_mvpoly_clear(&full);clear_input(p,m,a,n-1);fq_nmod_ctx_clear(fq);flint_rand_clear(rng);
    flint_cleanup_master();return 0;
}

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
    unified_mpoly_struct *p = random_mq(n - 1, fq, rng), **m, **a, full = {0}, actual = {0};
    build_fq_cancellation_matrix(&m, p, n - 1, 1);
    perform_fq_matrix_row_operations(&a, &m, n - 1, 1);
    double start=pencil_test_time();
    compute_fq_cancel_matrix_det(&full,a,n-1,1,DET_METHOD_RECURSIVE);
    double dp=pencil_test_time()-start;
    nmod_mpoly_ctx_t ctx;nmod_mpoly_ctx_init(ctx,2*n-1,ORD_LEX,prime);
    nmod_mpoly_struct *saved=flint_malloc(n*n*sizeof(*saved));
    for(slong i=0;i<n;i++)for(slong j=0;j<n;j++) {
        nmod_mpoly_init(saved+i*n+j,ctx);
        dr_mpoly_to_nmod_mpoly(saved + i * n + j, &a[i][j], ctx);
    }
    mq_pencil_stats stats;
    if(prime<=(ulong)n-1) {
        dr_mpoly_copy(&actual, &full);
        const void *ptr = actual.data.nmod_poly.coeffs;
        assert(!compute_fq_det_mq_pencil(&actual,a,n,&stats));
        assert(actual.data.nmod_poly.coeffs == ptr &&
               dr_mpoly_length(&(actual)) == dr_mpoly_length(&(full)));
        dr_mpoly_clear(&actual);
        printf("{\"n\":%ld,\"q\":%lu,\"characteristic_rejected\":true,\"preserved\":true}\n",n,prime);
    } else {
        int ok=compute_fq_det_mq_pencil(&actual,a,n,&stats);
        if(!ok) {
            assert(strstr(stats.reason,"rank deficient"));
            printf("{\"n\":%ld,\"q\":%lu,\"rank_rejected\":true}\n",n,prime);
        } else {
            nmod_mpoly_t x,y;nmod_mpoly_init(x,ctx);nmod_mpoly_init(y,ctx);
            dr_mpoly_to_nmod_mpoly(x, &actual, ctx);
            dr_mpoly_to_nmod_mpoly(y, &full, ctx);
            assert(nmod_mpoly_equal(x,y,ctx));
            for(slong i=0;i<n;i++)for(slong j=0;j<n;j++) {
                    dr_mpoly_to_nmod_mpoly(x, &a[i][j], ctx);
                    assert(nmod_mpoly_equal(x, saved + i * n + j, ctx));
            }
            printf("{\"n\":%ld,\"q\":%lu,\"threads\":%ld,\"terms\":%ld,\"full_dp\":%.9f,"
                   "\"normalization\":%.9f,\"recurrence\":%.9f,\"assembly\":%.9f,\"total\":%.9f,"
                   "\"peak_terms\":%ld,\"equal\":true}\n",
                   n, prime, stats.threads, dr_mpoly_length(&(actual)), dp, stats.normalization,
                   stats.recurrence, stats.assembly, stats.total, stats.peak_terms);
            /* Non-downward-closed, rectangular targets: intermediate divisors
             * must survive even though they are not final requested terms. */
            slong rows[24]={0},cols[16]={0};
            rows[0]=2;rows[n-1+1]=1;rows[2*(n-1)+n-2]=2;
            cols[0]=1;cols[2*(n-1)-1]=2;
            unified_mpoly_struct projected = {0}, expected = {0};
            assert(compute_fq_det_mq_pencil_projected(&projected,a,n,rows,3,cols,2,&stats));
            assert(fq_mq_project_full(&expected,&full,rows,3,cols,2));
            dr_mpoly_to_nmod_mpoly(x, &projected, ctx);
            dr_mpoly_to_nmod_mpoly(y, &expected, ctx);
            assert(nmod_mpoly_equal(x,y,ctx));
            dr_mpoly_clear(&projected);
            dr_mpoly_clear(&expected);
            if(n>=5) {
                /* Larger irregular ideals exercise packed closure lookups. */
                slong wide_rows[128]={0},wide_cols[128]={0};
                for(slong i=0;i<16;i++) {
                    for(slong d=0;d<i%n;d++)wide_rows[i*(n-1)+n_randint(rng,n-1)]++;
                    for(slong d=0;d<(3*i+1)%n;d++)wide_cols[i*(n-1)+n_randint(rng,n-1)]++;
                }
                assert(compute_fq_det_mq_pencil_projected(&projected,a,n,wide_rows,16,wide_cols,16,&stats));
                assert(fq_mq_project_full(&expected,&full,wide_rows,16,wide_cols,16));
                dr_mpoly_to_nmod_mpoly(x, &projected, ctx);
                dr_mpoly_to_nmod_mpoly(y, &expected, ctx);
                assert(nmod_mpoly_equal(x,y,ctx));
                dr_mpoly_clear(&projected);
                dr_mpoly_clear(&expected);
            }
            rows[0]=-1;
            const void *preserved = actual.data.nmod_poly.coeffs;
            assert(!compute_fq_det_mq_pencil_projected(&actual,a,n,rows,3,cols,2,&stats));
            assert(actual.data.nmod_poly.coeffs == preserved);
            assert(!compute_fq_det_mq_pencil_projected(&actual,a,n,rows,0,cols,2,&stats));
            assert(actual.data.nmod_poly.coeffs == preserved && stats.reason);
            dr_mpoly_clear(&actual);
            nmod_mpoly_clear(x,ctx);nmod_mpoly_clear(y,ctx);
            /* Force C=[I|0], so only the LAST free-column choice works.
             * For even n this also forces an odd column permutation. */
            if(n<=5) {
                slong exps[64]={0},par=1;
                fq_nmod_t one;fq_nmod_init(one,fq);fq_nmod_one(one,fq);
                for(slong i=1;i<n;i++)for(slong j=0;j<n;j++) {
                        for (slong k = 0; k < dr_mpoly_length(&(a[i][j])); k++) {
                            DR_MPOLY_TERM(term_1, &(a[i][j]), k);
                            if (term_1.par_exp[0])
                                nmod_mpoly_set_term_coeff_ui(GET_NMOD_POLY(&a[i][j]), k, 0,
                                                             GET_NMOD_CTX((&a[i][j])->ctx_ptr));
                        }
                        if (j == i - 1)
                            dr_mpoly_add_term_fast(&a[i][j], exps, &par, one);
                }
                fq_nmod_clear(one,fq);
                unified_mpoly_struct expected = {0};
                compute_fq_cancel_matrix_det(&expected,a,n-1,1,DET_METHOD_RECURSIVE);
                assert(compute_fq_det_mq_pencil(&actual,a,n,&stats));
                nmod_mpoly_init(x,ctx);nmod_mpoly_init(y,ctx);
                dr_mpoly_to_nmod_mpoly(x, &actual, ctx);
                dr_mpoly_to_nmod_mpoly(y, &expected, ctx);
                assert(nmod_mpoly_equal(x,y,ctx));
                nmod_mpoly_clear(x,ctx);nmod_mpoly_clear(y,ctx);
                dr_mpoly_clear(&expected);
                dr_mpoly_clear(&actual);
            }
            /* Singular complete determinant, but the parameter normalization
             * remains valid: zero must be a successful exact result. */
            if(n<=5) {
                for (slong j = 0; j < n; j++)
                    unified_mpoly_zero(&a[0][j]);
                assert(compute_fq_det_mq_pencil(&actual,a,n,&stats));
                assert(dr_mpoly_length(&(actual)) == 0);
                dr_mpoly_clear(&actual);
            }
        }
        /* Rank-deficient parameter coefficient matrix rejects without writing
         * to the caller's already initialized polynomial. */
        for (slong i = 1; i < n; i++)
            for (slong j = 0; j < n; j++)
                for (slong k = 0; k < dr_mpoly_length(&(a[i][j])); k++) {
                    DR_MPOLY_TERM(term_3, &(a[i][j]), k);
                    if (term_3.par_exp[0]) {
                        nmod_mpoly_set_term_coeff_ui(GET_NMOD_POLY(&a[i][j]), k, 0,
                                                     GET_NMOD_CTX((&a[i][j])->ctx_ptr));
                        a[i][j].canonical = 0;
                    }
                }
        dr_mpoly_copy(&actual, &full);
        const void *ptr = actual.data.nmod_poly.coeffs;
        assert(!compute_fq_det_mq_pencil(&actual,a,n,&stats));
        assert(strstr(stats.reason,"rank deficient"));
        assert(actual.data.nmod_poly.coeffs == ptr &&
               dr_mpoly_length(&(actual)) == dr_mpoly_length(&(full)));
        dr_mpoly_clear(&actual);
    }
    if(n==4 && prime==65537) {
        /* n=10 must reach both recurrence entry points without a dense MQ
         * allocation. diag(1,t,...,t) has the exact determinant t^9. n=11
         * still rejects before touching an initialized caller output. */
        for(slong size=10;size<=11;size++) {
            unified_mpoly_struct **large = flint_calloc(1, size * sizeof(*large));
            for(slong i=0;i<size;i++) {
                large[i] = flint_calloc(1, size * sizeof(**large));
                for (slong j = 0; j < size; j++)
                    dr_mpoly_init(&large[i][j], 2 * (size - 1), 1, fq);
            }
            if(size==10) {
                slong exps[20]={0},par=0;
                fq_nmod_t one;fq_nmod_init(one,fq);fq_nmod_one(one,fq);
                for(slong i=0;i<size;i++) {
                    par=i?1:0;
                    dr_mpoly_add_term_fast(&large[i][i], exps, &par, one);
                }
                fq_nmod_clear(one,fq);
                for(int projected=0;projected<2;projected++) {
                    assert(projected
                        ? compute_fq_det_mq_pencil_projected(&actual,large,size,exps,1,exps,1,&stats)
                        : compute_fq_det_mq_pencil(&actual,large,size,&stats));
                    DR_MPOLY_TERM(term_4, &(actual), 0);
                    DR_MPOLY_TERM(term_5, &(actual), 0);
                    assert(dr_mpoly_length(&(actual)) == 1 && term_5.par_exp[0] == 9);
                    DR_MPOLY_TERM(term_6, &(actual), 0);
                    assert(fq_nmod_is_one(term_6.coeff, fq));
                    for (slong v = 0; v < 18; v++) {
                        DR_MPOLY_TERM(term_7, &(actual), 0);
                        assert(term_7.var_exp[v] == 0);
                    }
                    dr_mpoly_clear(&actual);
                }
                puts("n=10 full/projected pencil admission and t^9 determinant passed");
            } else {
                dr_mpoly_copy(&actual, &full);
                const void *ptr = actual.data.nmod_poly.coeffs;
                assert(!compute_fq_det_mq_pencil(&actual,large,size,&stats));
                assert(strstr(stats.reason, "workspace") && actual.data.nmod_poly.coeffs == ptr);
                dr_mpoly_clear(&actual);
            }
            for(slong i=0;i<size;i++) {
                for (slong j = 0; j < size; j++)
                    dr_mpoly_clear(&large[i][j]);
                flint_free(large[i]);
            }
            flint_free(large);
        }
    }
    for(slong i=0;i<n*n;i++)nmod_mpoly_clear(saved+i,ctx);
    flint_free(saved);nmod_mpoly_ctx_clear(ctx);
    dr_mpoly_clear(&full);
    clear_input(p, m, a, n - 1);
    fq_nmod_ctx_clear(fq);
    flint_rand_clear(rng);
    flint_cleanup_master();return 0;
}

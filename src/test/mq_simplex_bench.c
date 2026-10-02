/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Exact full-polynomial comparison and stage timings for the production engine.
 * n is the equation count. Timings exclude common cancellation row operations. */
#define main mq_filter_regression_main
#include "dixon_mq_filter_test.c"
#undef main
#include "mq_simplex_det.h"
#include <time.h>

static double seconds(void)
{
    struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);
    return t.tv_sec+1e-9*t.tv_nsec;
}

int main(int argc,char **argv)
{
    slong n=argc>1?atol(argv[1]):5;
    ulong prime=argc>2?strtoul(argv[2],NULL,10):65537;
    int threads=argc>3?atoi(argv[3]):1;
    if(n<3 || n>8 || prime<2 || !n_is_prime(prime) || threads<1 || threads>64)return 2;
    omp_set_num_threads(threads);flint_set_num_threads(threads);
    g_dixon_verbose_level=0;g_dixon_det_cache_limit=100000;
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,132,941);
    fq_nmod_ctx_t fq;fq_nmod_ctx_init_ui(fq,prime,1,"a");
    unified_mpoly_struct *p = random_mq(n - 1, fq, rng), **m, **a, full = {0}, projected = {0},
                         actual = {0};
    build_fq_cancellation_matrix(&m, p, n - 1, 1);
    perform_fq_matrix_row_operations(&a, &m, n - 1, 1);
    double start=seconds();
    compute_fq_cancel_matrix_det(&full,a,n-1,1,DET_METHOD_RECURSIVE);
    double dp=seconds()-start;
    start=seconds();
    int success=dixon_try_mq_projection(&projected,a,p,n-1,1,DET_METHOD_RECURSIVE);
    double proj=seconds()-start;
    if (success)
        dr_mpoly_clear(&projected);
    mq_simplex_stats stats;
    if(prime<=(ulong)n+1) {
        /* Rejection must leave an already initialized caller output intact. */
        dr_mpoly_copy(&actual, &full);
        const void *saved = actual.data.nmod_poly.coeffs;
        assert(!compute_fq_det_mq_simplex(&actual,a,n,&stats));
        assert(actual.data.nmod_poly.coeffs == saved &&
               dr_mpoly_length(&(actual)) == dr_mpoly_length(&(full)));
        dr_mpoly_clear(&actual);
        printf("{\"n\":%ld,\"q\":%lu,\"rejected\":true,\"preserved\":true}\n",n,prime);
    } else {
        assert(compute_fq_det_mq_simplex(&actual,a,n,&stats));
        nmod_mpoly_ctx_t ctx;nmod_mpoly_ctx_init(ctx,2*n-1,ORD_LEX,prime);
        nmod_mpoly_t x,y;nmod_mpoly_init(x,ctx);nmod_mpoly_init(y,ctx);
        dr_mpoly_to_nmod_mpoly(x, &actual, ctx);
        dr_mpoly_to_nmod_mpoly(y, &full, ctx);
        assert(nmod_mpoly_equal(x,y,ctx));
        /* Cached projection and repair path must produce the same selected
         * coefficients as ordinary DP, including rectangular target sets. */
        slong rows[32]={0},cols[32]={0}; rows[n-1]=2;cols[2*(n-1)-1]=2;
        unified_mpoly_struct cached = {0}, independent = {0};
        assert(fq_mq_project_full(&cached,&actual,rows,2,cols,2));
        assert(compute_fq_det_mq_projected(&independent,a,n,rows,cols,2));
        dr_mpoly_to_nmod_mpoly(x, &cached, ctx);
        dr_mpoly_to_nmod_mpoly(y, &independent, ctx);
        assert(nmod_mpoly_equal(x,y,ctx));
        dr_mpoly_clear(&cached);
        dr_mpoly_clear(&independent);
        printf("{\"n\":%ld,\"q\":%lu,\"threads\":%ld,\"points\":%ld,\"terms\":%ld,\"full_dp\":%.9f,"
               "\"projected_dp\":%.9f,\"projection_success\":%d,\"setup\":%.9f,\"entry_eval\":%.9f,"
               "\"determinants\":%.9f,\"interpolation\":%.9f,\"packing\":%.9f,\"total\":%.9f,"
               "\"equal\":true}\n",
               n, prime, stats.threads, stats.points, dr_mpoly_length(&(actual)), dp, proj, success,
               stats.setup, stats.entry_eval, stats.determinants, stats.interpolation,
               stats.packing, stats.total);
        nmod_mpoly_clear(x,ctx);nmod_mpoly_clear(y,ctx);nmod_mpoly_ctx_clear(ctx);
        dr_mpoly_clear(&actual);
    }
    dr_mpoly_clear(&full);
    clear_input(p, m, a, n - 1);
    fq_nmod_ctx_clear(fq);
    flint_rand_clear(rng);
    flint_cleanup_master();return 0;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
#define main mq_layout_fixture_main
#include "dixon_mq_filter_test.c"
#undef main
#include <sys/resource.h>
#include <stdint.h>

void mq_layout_configure(int enabled, int profile, int rotate);
int mq_layout_shared_calls(void);
double mq_layout_plan_seconds(void);
double mq_layout_arithmetic_seconds(void);

static int layout_compute(fq_mvpoly_t *out, fq_mvpoly_t **matrix,
                           fq_mvpoly_t *polys, slong n, int projected)
{
    if (projected && dixon_try_mq_projection(out, matrix, polys, n-1, 1, DET_METHOD_RECURSIVE)) return 1;
    compute_fq_cancel_matrix_det(out, matrix, n-1, 1, DET_METHOD_RECURSIVE);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc != 10) {
        fprintf(stderr, "usage: %s n threads shared projected quadratic-first profile audit seed prime\n", argv[0]);
        return 2;
    }
    slong n = atol(argv[1]); int threads = atoi(argv[2]), shared = atoi(argv[3]);
    int projected = atoi(argv[4]), rotate = atoi(argv[5]), profile = atoi(argv[6]), audit = atoi(argv[7]);
    ulong seed = strtoul(argv[8], NULL, 10), prime = strtoul(argv[9], NULL, 10);
    if (n < 4 || n > 8 || threads < 1 || threads > 16 || !n_is_prime(prime)) return 2;
    omp_set_num_threads(threads); flint_set_num_threads(threads);
    g_dixon_verbose_level = 0; g_dixon_det_cache_limit = 100000;
    flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng, seed, 941);
    fq_nmod_ctx_t fq; fq_nmod_ctx_init_ui(fq, prime, 1, "a");
    fq_mvpoly_t *p = random_mq(n-1, fq, rng), **m, **a, actual, expected;
    build_fq_cancellation_matrix_mvpoly(&m, p, n-1, 1);
    perform_fq_matrix_row_operations_mvpoly(&a, &m, n-1, 1);
    mq_layout_configure(shared, profile, rotate);
    double start = omp_get_wtime();
    int verified = layout_compute(&actual, a, p, n, projected);
    double seconds = omp_get_wtime()-start;
    int shared_calls = mq_layout_shared_calls();
    if (shared >= 4 && !rotate && prime == 65537) assert(shared_calls > 0);
    double plan_seconds = mq_layout_plan_seconds(), arithmetic_seconds = mq_layout_arithmetic_seconds();
    struct rusage usage; getrusage(RUSAGE_SELF, &usage);
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx, 2*n-1, ORD_LEX, prime);
    nmod_mpoly_t x, y; nmod_mpoly_init(x, ctx); nmod_mpoly_init(y, ctx);
    fq_mvpoly_to_nmod_mpoly(x, &actual, ctx);
    assert(nmod_mpoly_is_canonical(x, ctx));
    if (audit) {
        mq_layout_configure(0, 0, 0);
        int reference_verified = layout_compute(&expected, a, p, n, projected);
        assert(reference_verified == verified);
        fq_mvpoly_to_nmod_mpoly(y, &expected, ctx);
        assert(nmod_mpoly_equal(x, y, ctx)); fq_mvpoly_clear(&expected);
    }
    /* Stable, canonical coefficient fingerprint for repeated timing runs.
     * The audit above uses exact coefficient equality, not this fingerprint. */
    uint64_t hash = UINT64_C(14695981039346656037); ulong exp[16];
    for (slong t = 0; t < x->length; t++) {
        nmod_mpoly_get_term_exp_ui(exp, x, t, ctx);
        hash = (hash ^ x->coeffs[t])*UINT64_C(1099511628211);
        for (slong v = 0; v < 2*n-1; v++) hash = (hash ^ exp[v])*UINT64_C(1099511628211);
    }
    printf("{\"kind\":\"run\",\"n\":%ld,\"q\":%lu,\"seed\":%lu,\"threads\":%d,\"shared\":%d,\"projected\":%d,\"verified\":%d,\"quadratic_first\":%d,\"profile\":%d,\"audit\":%d,\"seconds\":%.9f,\"peak_rss_kib_before_audit\":%ld,\"terms\":%ld,\"fingerprint\":\"%016llx\",\"plan_seconds\":%.9f,\"arithmetic_seconds\":%.9f}\n",
        n, prime, seed, threads, shared, projected, verified, rotate, profile, audit, seconds,
        usage.ru_maxrss, actual.nterms, (unsigned long long)hash, plan_seconds, arithmetic_seconds);
    nmod_mpoly_clear(x, ctx); nmod_mpoly_clear(y, ctx); nmod_mpoly_ctx_clear(ctx);
    fq_mvpoly_clear(&actual); clear_input(p, m, a, n-1);
    fq_nmod_ctx_clear(fq); flint_rand_clear(rng); flint_cleanup_master(); return 0;
}

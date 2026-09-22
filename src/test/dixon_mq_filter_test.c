/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Private integration helpers are exercised without extending the public ABI. */
#include "../dixon/dixon_flint.c"
#include <assert.h>

static fq_mvpoly_t *random_mq(slong n, const fq_nmod_ctx_t ctx, flint_rand_t state)
{
    fq_mvpoly_t *p = flint_malloc((size_t) (n + 1) * sizeof(*p));
    slong *e = flint_calloc((size_t) n + 1, sizeof(slong));
    fq_nmod_t c;
    fq_nmod_init(c, ctx);
    ulong prime = fq_nmod_ctx_modulus(ctx)->mod.n;
    for (slong k = 0; k <= n; k++) {
        fq_mvpoly_init(p + k, n, 1, ctx);
        for (slong i = -1; i <= n; i++) {
            for (slong j = i; j <= n; j++) {
                if (i >= 0) e[i]++;
                if (j >= 0) e[j]++;
                ulong coeff = n_randint(state, prime);
                /* Ensure elimination degree two even over F_2. */
                if (i == 0 && j == 0) coeff = 1;
                fq_nmod_set_ui(c, coeff, ctx);
                if (coeff) fq_mvpoly_add_term_fast(p + k, e, e + n, c);
                if (i >= 0) e[i]--;
                if (j >= 0) e[j]--;
            }
        }
    }
    fq_nmod_clear(c, ctx); flint_free(e);
    return p;
}

static void clear_input(fq_mvpoly_t *p, fq_mvpoly_t **m, fq_mvpoly_t **a, slong n)
{
    for (slong i = 0; i <= n; i++) {
        fq_mvpoly_clear(p + i);
        for (slong j = 0; j <= n; j++) {
            fq_mvpoly_clear(&m[i][j]); fq_mvpoly_clear(&a[i][j]);
        }
        flint_free(m[i]); flint_free(a[i]);
    }
    flint_free(p); flint_free(m); flint_free(a);
}

static int in_targets(const slong *exp, const slong *targets, slong count, slong n)
{
    for (slong i = 0; i < count; i++)
        if (!memcmp(exp, targets + i * n, (size_t) n * sizeof(slong))) return 1;
    return 0;
}

/* Compare every coefficient, including all parameter powers and coefficients
 * missing from the projected result, with an independently computed full D. */
static void check_block(const fq_mvpoly_t *actual, const fq_mvpoly_t *full,
                         const slong *rows, const slong *cols, slong count, slong n)
{
    fq_mvpoly_t expected;
    fq_mvpoly_init(&expected, 2 * n, 1, full->ctx);
    for (slong t = 0; t < full->nterms; t++) {
        const fq_monomial_t *term = &full->terms[t];
        if (in_targets(term->var_exp, rows, count, n) &&
            in_targets(term->var_exp + n, cols, count, n))
            fq_mvpoly_add_term_fast(&expected, term->var_exp, term->par_exp, term->coeff);
    }
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, 2 * n + 1, ORD_LEX, fq_nmod_ctx_modulus(full->ctx)->mod.n);
    nmod_mpoly_t x, y;
    nmod_mpoly_init(x, ctx); nmod_mpoly_init(y, ctx);
    fq_mvpoly_to_nmod_mpoly(x, actual, ctx); fq_mvpoly_to_nmod_mpoly(y, &expected, ctx);
    assert(nmod_mpoly_equal(x, y, ctx));
    nmod_mpoly_clear(x, ctx); nmod_mpoly_clear(y, ctx); nmod_mpoly_ctx_clear(ctx);
    fq_mvpoly_clear(&expected);
}

static void check_predicted_block(const fq_mvpoly_t *actual, const fq_mvpoly_t *full, slong n)
{
    /* Every target row/column occurs after successful full-rank verification. */
    slong *rows = flint_malloc((size_t) actual->nterms * n * sizeof(slong));
    slong *cols = flint_malloc((size_t) actual->nterms * n * sizeof(slong));
    for (slong t = 0; t < actual->nterms; t++) {
        memcpy(rows + t * n, actual->terms[t].var_exp, (size_t) n * sizeof(slong));
        memcpy(cols + t * n, actual->terms[t].var_exp + n, (size_t) n * sizeof(slong));
    }
    check_block(actual, full, rows, cols, actual->nterms, n);
    flint_free(rows); flint_free(cols);
}

static void check_coefficients(void)
{
    ulong primes[] = {2, 3, 101, 65537};
    slong limits[] = {0, 3, 100000};
    flint_rand_t state;
    flint_rand_init(state); flint_rand_set_seed(state, 132, 941);
    int successes = 0;
    for (slong pi = 0; pi < 4; pi++) for (slong n = 2; n <= 5; n++) {
        fq_nmod_ctx_t ctx;
        fq_nmod_ctx_init_ui(ctx, primes[pi], 1, "a");
        fq_mvpoly_t *p = random_mq(n, ctx, state), **m, **a, full, actual;
        build_fq_cancellation_matrix_mvpoly(&m, p, n, 1);
        perform_fq_matrix_row_operations_mvpoly(&a, &m, n, 1);
        g_dixon_det_cache_limit = 100000;
        compute_fq_cancel_matrix_det(&full, a, n, 1, DET_METHOD_RECURSIVE);
        slong *rows = flint_calloc((size_t) 3 * n, sizeof(slong));
        slong *cols = flint_calloc((size_t) 3 * n, sizeof(slong));
        rows[n] = cols[2 * n - 1] = 2;
        rows[2 * n + 1] = cols[3 * n - 2] = 1;
        /* Non-downward-closed targets force retention of intermediate terms. */
        for (int parallel = 0; parallel <= 1; parallel++) {
            omp_set_num_threads(parallel ? 4 : 1);
            for (int k = 0; k < 3; k++) {
                g_dixon_det_cache_limit = limits[k];
                assert(compute_fq_det_mq_projected(&actual, a, n + 1, rows, cols, 3));
                check_block(&actual, &full, rows, cols, 3, n);
                fq_mvpoly_clear(&actual);
            }
        }
        g_dixon_det_cache_limit = 100000;
        int ok = dixon_try_mq_projection(&actual, a, p, n, 1, DET_METHOD_RECURSIVE);
        if (ok) {
            check_predicted_block(&actual, &full, n);
            successes++;
            fq_mvpoly_clear(&actual);
        }
        printf("MQ projection p=%lu n=%ld: coefficients exact, candidate=%s\n",
               primes[pi], n, ok ? "verified" : "fallback/ineligible");
        assert(!dixon_try_mq_projection(&actual, a, p, n, 1, DET_METHOD_KRONECKER));
        g_dixon_mq_step1_filter = 0;
        assert(!dixon_try_mq_projection(&actual, a, p, n, 1, DET_METHOD_RECURSIVE));
        g_dixon_mq_step1_filter = 1;
        flint_free(rows); flint_free(cols);
        fq_mvpoly_clear(&full); clear_input(p, m, a, n); fq_nmod_ctx_clear(ctx);
    }
    assert(successes > 0);
    flint_rand_clear(state);
}

static void assert_same_result(const fq_mvpoly_t *a, const fq_mvpoly_t *b)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, 1, ORD_LEX, fq_nmod_ctx_modulus(a->ctx)->mod.n);
    nmod_mpoly_t x, y;
    nmod_mpoly_init(x, ctx); nmod_mpoly_init(y, ctx);
    fq_mvpoly_to_nmod_mpoly(x, a, ctx); fq_mvpoly_to_nmod_mpoly(y, b, ctx);
    assert(nmod_mpoly_equal(x, y, ctx));
    nmod_mpoly_clear(x, ctx); nmod_mpoly_clear(y, ctx); nmod_mpoly_ctx_clear(ctx);
}

static void check_entry_points(int homogeneous)
{
    const slong n = 3;
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx, 65537, 1, "a");
    flint_rand_t state;
    flint_rand_init(state); flint_rand_set_seed(state, 172, 591);
    fq_mvpoly_t *p = random_mq(n, ctx, state), **m, **a, projected;
    if (homogeneous) {
        for (slong i = 0; i <= n; i++) {
            fq_mvpoly_t h;
            fq_mvpoly_init(&h, n, 1, ctx);
            for (slong t = 0; t < p[i].nterms; t++) {
                fq_monomial_t *term = &p[i].terms[t];
                slong d = term->par_exp[0];
                for (slong v = 0; v < n; v++) d += term->var_exp[v];
                if (d == 2) fq_mvpoly_add_term_fast(&h, term->var_exp, term->par_exp, term->coeff);
            }
            fq_mvpoly_clear(p + i); p[i] = h;
        }
    }
    if (homogeneous == 2) {
        /* Replace t by t-1: point 1 is singular, but point 0 can certify.
         * This exercises retry/zero handling after changing point order. */
        for (slong i = 0; i <= n; i++) {
            fq_mvpoly_t shifted;
            fq_mvpoly_init(&shifted, n, 1, ctx);
            fq_nmod_t c;
            fq_nmod_init(c, ctx);
            for (slong t = 0; t < p[i].nterms; t++) {
                const fq_monomial_t *term = &p[i].terms[t];
                slong power = term->par_exp[0];
                for (slong d = 0; d <= power; d++) {
                    fq_nmod_set(c, term->coeff, ctx);
                    if (power == 2 && d == 1) fq_nmod_add(c, c, c, ctx);
                    if ((power - d) & 1) fq_nmod_neg(c, c, ctx);
                    fq_mvpoly_add_term(&shifted, term->var_exp, &d, c);
                }
            }
            fq_nmod_clear(c, ctx);
            fq_mvpoly_clear(p + i); p[i] = shifted;
        }
    }
    build_fq_cancellation_matrix_mvpoly(&m, p, n, 1);
    perform_fq_matrix_row_operations_mvpoly(&a, &m, n, 1);
    assert(dixon_try_mq_projection(&projected, a, p, n, 1, DET_METHOD_RECURSIVE));
    /* The public (unverified) extractor independently selects/evaluates the
     * already complete target block. Restore its extracted t-content. */
    fq_mvpoly_t **cm = NULL, expected, actual, named, forced;
    fq_nmod_poly_mat_t pm;
    slong *ri = flint_malloc((size_t) projected.nterms * sizeof(slong));
    slong *ci = flint_malloc((size_t) projected.nterms * sizeof(slong));
    slong size, content;
    long degrees[] = {2, 2, 2, 2};
    extract_fq_coefficient_matrix_from_dixon(&cm, &pm, ri, ci, &size, &content,
                              &projected, n, 1, NULL, NULL, NULL, degrees, n + 1);
    assert(size > 0);
    if (homogeneous == 1) assert(content > 0);
    fq_nmod_poly_t det;
    fq_nmod_poly_init(det, ctx); fq_nmod_poly_mat_det_iter(det, pm, ctx);
    fq_mvpoly_init(&expected, 0, 1, ctx);
    fq_nmod_t coeff;
    fq_nmod_init(coeff, ctx);
    for (slong k = 0; k <= fq_nmod_poly_degree(det, ctx); k++) {
        fq_nmod_poly_get_coeff(coeff, det, k, ctx);
        slong power = k + content;
        if (!fq_nmod_is_zero(coeff, ctx)) fq_mvpoly_add_term_fast(&expected, NULL, &power, coeff);
    }
    fq_mvpoly_make_monic(&expected);
    fq_dixon_resultant(&actual, p, n, 1);
    char *vars[] = {"x", "y", "z"}, *pars[] = {"t"};
    fq_dixon_resultant_with_names(&named, p, n, 1, vars, pars, "a");
    dixon_global_method_step4 = DET_METHOD_KRONECKER;
    fq_dixon_resultant(&forced, p, n, 1);
    dixon_global_method_step4 = -1;
    assert_same_result(&actual, &expected); assert_same_result(&named, &expected);
    assert_same_result(&forced, &expected);
    fq_mvpoly_clear(&forced); fq_mvpoly_clear(&named); fq_mvpoly_clear(&actual);
    fq_mvpoly_clear(&expected); fq_mvpoly_clear(&projected);
    fq_nmod_clear(coeff, ctx); fq_nmod_poly_clear(det, ctx); fq_nmod_poly_mat_clear(pm, ctx);
    flint_free(ri); flint_free(ci);
    clear_input(p, m, a, n); flint_rand_clear(state); fq_nmod_ctx_clear(ctx);
}

static void check_fallback_and_gates(void)
{
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx, 101, 1, "a");
    flint_rand_t state;
    flint_rand_init(state); flint_rand_set_seed(state, 771, 882);
    slong n = 3;
    fq_mvpoly_t *p = random_mq(n, ctx, state), **m, **a, result, old;
    /* An identically singular candidate MUST recompute the full polynomial. */
    fq_mvpoly_clear(p + 1); fq_mvpoly_copy(p + 1, p);
    build_fq_cancellation_matrix_mvpoly(&m, p, n, 1);
    perform_fq_matrix_row_operations_mvpoly(&a, &m, n, 1);
    assert(!dixon_try_mq_projection(&result, a, p, n, 1, DET_METHOD_RECURSIVE));
    fq_dixon_resultant(&result, p, n, 1);
    g_dixon_mq_step1_filter = 0;
    fq_dixon_resultant(&old, p, n, 1);
    g_dixon_mq_step1_filter = 1;
    assert_same_result(&result, &old); assert(result.nterms == 0);
    fq_mvpoly_clear(&result); fq_mvpoly_clear(&old);
    /* Parameter degree is part of the MQ gate. */
    slong exp[] = {0, 0, 0}, pe = 3;
    fq_nmod_t c; fq_nmod_init(c, ctx); fq_nmod_one(c, ctx);
    fq_mvpoly_add_term_fast(p, exp, &pe, c);
    assert(!dixon_try_mq_projection(&result, a, p, n, 1, DET_METHOD_RECURSIVE));
    slong targets[] = {0, 0, 0};
    a[0][0].terms[0].par_exp[0] = 3;
    assert(!compute_fq_det_mq_projected(&result, a, n + 1, targets, targets, 1));
    fq_nmod_clear(c, ctx); clear_input(p, m, a, n);
    fq_nmod_ctx_clear(ctx);
    /* Extension fields and multiple parameters keep their existing backend. */
    fq_nmod_ctx_init_ui(ctx, 3, 2, "a");
    p = random_mq(n, ctx, state);
    build_fq_cancellation_matrix_mvpoly(&m, p, n, 1);
    perform_fq_matrix_row_operations_mvpoly(&a, &m, n, 1);
    assert(!dixon_try_mq_projection(&result, a, p, n, 1, DET_METHOD_RECURSIVE));
    assert(!compute_fq_det_mq_projected(&result, a, n + 1, targets, targets, 1));
    assert(!dixon_try_mq_projection(&result, a, p, n, 2, DET_METHOD_RECURSIVE));
    assert(!dixon_try_mq_projection(&result, a, p, n, 0, DET_METHOD_RECURSIVE));
    clear_input(p, m, a, n); fq_nmod_ctx_clear(ctx); flint_rand_clear(state);
}

/* Reproducible timing probe; full is determinant only, projected includes
 * candidate construction AND rank verification. No speed assertion. */
static void benchmark_projection(void)
{
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx, 65537, 1, "a");
    flint_rand_t state;
    flint_rand_init(state); flint_rand_set_seed(state, 12, 53);
    g_dixon_det_cache_limit = 100000;
    for (slong n = 3; n <= 7; n++) {
        fq_mvpoly_t *p = random_mq(n, ctx, state), **m, **a, full, projected;
        build_fq_cancellation_matrix_mvpoly(&m, p, n, 1);
        perform_fq_matrix_row_operations_mvpoly(&a, &m, n, 1);
        double start = get_wall_time();
        compute_fq_cancel_matrix_det(&full, a, n, 1, DET_METHOD_RECURSIVE);
        double full_time = get_wall_time() - start;
        slong full_terms = full.nterms;
        fq_mvpoly_clear(&full);
        start = get_wall_time();
        int ok = dixon_try_mq_projection(&projected, a, p, n, 1, DET_METHOD_RECURSIVE);
        double projected_time = get_wall_time() - start;
        printf("n=%ld full_terms=%ld projected_terms=%ld full=%.4fs projected+verify=%.4fs accepted=%d\n",
               n, full_terms, ok ? projected.nterms : 0, full_time, projected_time, ok);
        if (ok) fq_mvpoly_clear(&projected);
        clear_input(p, m, a, n);
    }
    flint_rand_clear(state); fq_nmod_ctx_clear(ctx);
}

int main(int argc, char **argv)
{
    g_dixon_verbose_level = 0;
    omp_set_num_threads(4);
    assert(g_dixon_mq_step1_filter == 1);
    unsetenv("DRSOLVE_PREDICT_MAXRANK");
    if (argc == 2 && strcmp(argv[1], "--bench") == 0) {
        benchmark_projection();
        flint_cleanup_master();
        return 0;
    }
    check_coefficients();
    check_entry_points(0);
    check_entry_points(1);
    check_entry_points(2);
    check_fallback_and_gates();
    puts("MQ projection coefficient, integration and fallback tests passed");
    flint_cleanup_master();
    return 0;
}

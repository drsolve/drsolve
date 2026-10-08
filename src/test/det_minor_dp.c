/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_mpoly_det.h"
#include "fq_mpoly_mat_det.h"
#include <flint/nmod_mat.h>

/* Observe the actual determinant backends in the instrumented test build.
 * Normal library builds do not contain these callbacks. */
static struct {
    int enabled, parallel, begins, root_products, wrong_team;
    slong size;
} scheduling_probe;

void drsolve_det_test_event(int event, slong size)
{
    if (!scheduling_probe.enabled || size != scheduling_probe.size) return;
    if (event == 0) {
#ifdef _OPENMP
        #pragma omp atomic update
#endif
        scheduling_probe.begins++;
    } else {
#ifdef _OPENMP
        #pragma omp atomic update
#endif
        scheduling_probe.root_products++;
#ifdef _OPENMP
        if ((omp_get_num_threads() > 1) != scheduling_probe.parallel) {
            #pragma omp atomic update
            scheduling_probe.wrong_team++;
        }
#endif
    }
}

static void start_scheduling_probe(slong size, int parallel)
{
    memset(&scheduling_probe, 0, sizeof(scheduling_probe));
    scheduling_probe.enabled = 1;
    scheduling_probe.size = size;
    scheduling_probe.parallel = parallel;
}

static void finish_scheduling_probe(int begins, int products)
{
    scheduling_probe.enabled = 0;
    if (scheduling_probe.begins != begins ||
        scheduling_probe.root_products != products || scheduling_probe.wrong_team) {
        fprintf(stderr, "minor DP scheduling regression: n=%ld parallel=%d begins=%d products=%d wrong_team=%d\n",
                scheduling_probe.size, scheduling_probe.parallel, scheduling_probe.begins,
                scheduling_probe.root_products, scheduling_probe.wrong_team);
        abort();
    }
}

static void require_equal(unified_mpoly_t actual, unified_mpoly_t expected,
                          slong n, slong limit, int parallel)
{
    if (!unified_mpoly_equal(actual, expected)) {
        fprintf(stderr, "minor DP mismatch: field=%d n=%ld limit=%ld parallel=%d\n",
                actual->field_id, n, limit, parallel);
        abort();
    }
}

/* Exercise the public prime-field route, which bypasses unified_mpoly. */
static void check_direct_nmod(unified_mpoly_t **matrix, slong n,
                              unified_mpoly_t expected, const fq_nmod_ctx_t fq)
{
    unified_mpoly_ctx_t ctx = expected->ctx_ptr;
    unified_mpoly_struct **input = calloc(1, (size_t)n * sizeof(*input));
    unified_mpoly_struct result = {0};
    nmod_mpoly_t actual;
    nmod_mpoly_init(actual, GET_NMOD_CTX(ctx));
    for (slong i = 0; i < n; i++) {
        input[i] = calloc(1, (size_t)n * sizeof(**input));
        for (slong j = 0; j < n; j++)
            nmod_mpoly_to_dr_mpoly(&input[i][j], GET_NMOD_POLY(matrix[i][j]), 2, 0,
                                   GET_NMOD_CTX(ctx), fq);
    }
    compute_fq_det_unified_interface(&result, input, n);
    dr_mpoly_to_nmod_mpoly(actual, &result, GET_NMOD_CTX(ctx));
    if (!nmod_mpoly_equal(actual, GET_NMOD_POLY(expected), GET_NMOD_CTX(ctx))) abort();
    dr_mpoly_clear(&result);
    nmod_mpoly_clear(actual, GET_NMOD_CTX(ctx));
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++)
            dr_mpoly_clear(&input[i][j]);
        free(input[i]);
    }
    free(input);
}

static void check_root_scheduling(void)
{
    fq_nmod_ctx_t fq;
    field_ctx_t field;
    fq_nmod_ctx_init_ui(fq, 101, 1, "a");
    field_ctx_init_enhanced(&field, fq, 0);
    unified_mpoly_ctx_t ctx = unified_mpoly_ctx_init(2, ORD_LEX, &field);
    nmod_mpoly_t factor;
    nmod_mpoly_init(factor, GET_NMOD_CTX(ctx));
    const char *vars[] = {"x", "y"};
    if (nmod_mpoly_set_str_pretty(factor, "x+y+1", vars, GET_NMOD_CTX(ctx))) abort();

    for (slong n = 4; n <= 7; n++) {
        nmod_mat_t numeric;
        nmod_mat_init(numeric, n, n, 101);
        unified_mpoly_t **matrix = unified_mpoly_mat_init(n, n, ctx);
        unified_mpoly_t expected = unified_mpoly_init(ctx);
        unified_mpoly_t actual = unified_mpoly_init(ctx);
        unified_mpoly_t entry = unified_mpoly_init(ctx);
        for (slong i = 0; i < n; i++) {
            for (slong j = 0; j < n; j++) {
                ulong value = 1;
                for (slong e = 0; e < (i == 0 ? n - 1 : i - 1); e++)
                    value = (value * (j + 1)) % 101;
                nmod_mat_entry(numeric, i, j) = value;
                nmod_mpoly_scalar_mul_ui(GET_NMOD_POLY(matrix[i][j]), factor, value, GET_NMOD_CTX(ctx));
            }
        }
        /* Vandermonde minors are nonzero. Both variables ensure that the
         * public nmod route exercises minor DP rather than univariate code. */
        if (!nmod_mpoly_pow_ui(GET_NMOD_POLY(expected), factor, n, GET_NMOD_CTX(ctx))) abort();
        nmod_mpoly_scalar_mul_ui(GET_NMOD_POLY(expected), GET_NMOD_POLY(expected),
                                  nmod_mat_det(numeric), GET_NMOD_CTX(ctx));
        for (int parallel = 0; parallel <= 1; parallel++) {
#ifdef _OPENMP
            omp_set_num_threads(parallel ? 4 : 1);
#endif
            g_dixon_det_cache_limit = n == 6 ? 35 : 1024;
            for (int direct = 0; direct <= 1; direct++) {
                start_scheduling_probe(n, parallel);
                if (direct) check_direct_nmod(matrix, n, expected, fq);
                else {
                    compute_unified_mpoly_det(actual, matrix, n, ctx, parallel);
                    require_equal(actual, expected, n, g_dixon_det_cache_limit, parallel);
                }
                finish_scheduling_probe(1, n);
            }
        }
        if (n == 6) {
#ifdef _OPENMP
            omp_set_num_threads(1);
#endif
            /* A 34-entry budget requires six n=5 DP children; zero disables
             * DP. Check both real backends and the read-only matrix views. */
            for (int direct = 0; direct <= 1; direct++) {
                for (int disabled = 0; disabled <= 1; disabled++) {
                    g_dixon_det_cache_limit = disabled ? 0 : 34;
                    start_scheduling_probe(5, 0);
                    if (direct) check_direct_nmod(matrix, n, expected, fq);
                    else {
                        compute_unified_mpoly_det(actual, matrix, n, ctx, 0);
                        require_equal(actual, expected, n, g_dixon_det_cache_limit, 0);
                    }
                    finish_scheduling_probe(disabled ? 0 : 6, disabled ? 0 : 30);
                }
            }
            for (slong i = 0; i < n; i++) {
                for (slong j = 0; j < n; j++) {
                    nmod_mpoly_scalar_mul_ui(GET_NMOD_POLY(entry), factor,
                                              nmod_mat_entry(numeric, i, j), GET_NMOD_CTX(ctx));
                    require_equal(matrix[i][j], entry, n, 34, 0);
                }
            }
        }
        unified_mpoly_clear(entry);
        unified_mpoly_clear(expected);
        unified_mpoly_clear(actual);
        unified_mpoly_mat_clear(matrix, n, n);
        nmod_mat_clear(numeric);
    }
    nmod_mpoly_clear(factor, GET_NMOD_CTX(ctx));
    unified_mpoly_ctx_clear(ctx);
    field_ctx_clear(&field);
    fq_nmod_ctx_clear(fq);
#ifdef _OPENMP
    omp_set_num_threads(4);
#endif
}

static void check_field(ulong prime, slong degree, ulong zech_limit)
{
    fq_nmod_ctx_t fq;
    field_ctx_t field;
    flint_rand_t state;
    fq_nmod_ctx_init_ui(fq, prime, degree, "a");
    field_ctx_init_enhanced(&field, fq, zech_limit);
    unified_mpoly_ctx_t ctx = unified_mpoly_ctx_init(2, ORD_LEX, &field);
    flint_rand_init(state);
    flint_rand_set_seed(state, 123, 456);
    for (slong n = 0; n <= 7; n++) {
        unified_mpoly_t **matrix = unified_mpoly_mat_init(n, n, ctx);
        unified_mpoly_t expected = unified_mpoly_init(ctx);
        unified_mpoly_t actual = unified_mpoly_init(ctx);
        for (int variant = 0; variant < 4; variant++) {
            for (slong i = 0; i < n; i++) {
                for (slong j = 0; j < n; j++) {
                    unified_mpoly_randtest(matrix[i][j], state, 2, 2);
                    if ((variant == 1 && (i + j) % 3 == 0) ||
                        (variant == 2 && i > j)) unified_mpoly_zero(matrix[i][j]);
                }
            }
            if (variant == 3 && n > 1)
                for (slong j = 0; j < n; j++) unified_mpoly_set(matrix[1][j], matrix[0][j]);
            compute_unified_mpoly_det_recursive(expected, matrix, n, ctx);
            /* n=6 needs exactly 35 entries for two adjacent layers. */
            const slong limits[] = {0, 1, 34, 35, 1024};
            for (size_t l = 0; l < sizeof(limits) / sizeof(limits[0]); l++) {
                g_dixon_det_cache_limit = limits[l];
                for (int parallel = 0; parallel <= 1; parallel++) {
                    compute_unified_mpoly_det(actual, matrix, n, ctx, parallel);
                    require_equal(actual, expected, n, limits[l], parallel);
                }
                if (field.field_id == FIELD_ID_NMOD && n >= 4)
                    check_direct_nmod(matrix, n, expected, fq);
            }
        }
        unified_mpoly_clear(expected);
        unified_mpoly_clear(actual);
        unified_mpoly_mat_clear(matrix, n, n);
    }

    if (field.field_id == FIELD_ID_NMOD) {
        /* Larger layers, checked against independent numeric elimination. */
        slong n = 11;
        nmod_mat_t numeric;
        nmod_mat_init(numeric, n, n, prime);
        nmod_mat_randtest(numeric, state);
        unified_mpoly_t **matrix = unified_mpoly_mat_init(n, n, ctx);
        unified_mpoly_t expected = unified_mpoly_init(ctx);
        unified_mpoly_t actual = unified_mpoly_init(ctx);
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++)
                nmod_mpoly_set_ui(GET_NMOD_POLY(matrix[i][j]),
                                  nmod_mat_entry(numeric, i, j), GET_NMOD_CTX(ctx));
        nmod_mpoly_set_ui(GET_NMOD_POLY(expected), nmod_mat_det(numeric), GET_NMOD_CTX(ctx));
        for (slong limit = 923; limit <= 924; limit++) {
            g_dixon_det_cache_limit = limit;
            for (int parallel = 0; parallel <= 1; parallel++) {
                compute_unified_mpoly_det(actual, matrix, n, ctx, parallel);
                require_equal(actual, expected, n, limit, parallel);
            }
        }
        unified_mpoly_clear(expected);
        unified_mpoly_clear(actual);
        unified_mpoly_mat_clear(matrix, n, n);
        nmod_mat_clear(numeric);
    }
    printf("minor DP tests passed: field=%d\n", field.field_id);
    flint_rand_clear(state);
    unified_mpoly_ctx_clear(ctx);
    field_ctx_clear(&field);
    fq_nmod_ctx_clear(fq);
}

/* Exercise native-byte merges across packing widths and monomial orders. */
static void check_byte_packing(slong degree)
{
    fq_nmod_ctx_t fq;
    fq_nmod_ctx_init_ui(fq, 2, degree, "a");
    field_ctx_t field;
    field_ctx_init(&field, fq);
    flint_rand_t random;
    flint_rand_init(random); flint_rand_set_seed(random, degree, 89);
    ordering_t orders[] = {ORD_LEX, ORD_DEGLEX, ORD_DEGREVLEX};
    for (int order = 0; order < 3; ++order) {
        unified_mpoly_ctx_t ctx = unified_mpoly_ctx_init(3, orders[order], &field);
        unified_mpoly_t **m = unified_mpoly_mat_init(4, 4, ctx);
        unified_mpoly_t **copy = unified_mpoly_mat_init(4, 4, ctx);
        unified_mpoly_t expected = unified_mpoly_init(ctx), actual = unified_mpoly_init(ctx);
        for (int trial = 0; trial < 4; ++trial) {
            for (slong i = 0; i < 4; ++i)
                for (slong j = 0; j < 4; ++j) {
                    fq_nmod_mpoly_randtest_bits(GET_FQ_POLY(m[i][j]), random, 2,
                        (i+j) % 3 ? 3 : trial % 2 ? 80 : 15, GET_FQ_CTX(ctx));
                    if (trial > 1 && i == 1) unified_mpoly_set(m[i][j], m[0][j]);
                    unified_mpoly_set(copy[i][j], m[i][j]);
                }
            compute_unified_mpoly_det_recursive(expected, m, 4, ctx);
            g_dixon_det_cache_limit = 1024;
            for (int parallel = 0; parallel <= 1; ++parallel) {
                compute_unified_mpoly_det(actual, m, 4, ctx, parallel);
                require_equal(actual, expected, 4, 1024, parallel);
                for (slong i = 0; i < 4; ++i)
                    for (slong j = 0; j < 4; ++j)
                        require_equal(m[i][j], copy[i][j], 4, 1024, parallel);
            }
        }
        unified_mpoly_clear(expected); unified_mpoly_clear(actual);
        unified_mpoly_mat_clear(m, 4, 4); unified_mpoly_mat_clear(copy, 4, 4);
        unified_mpoly_ctx_clear(ctx);
    }
    flint_rand_clear(random); field_ctx_clear(&field); fq_nmod_ctx_clear(fq);
}

int main(void)
{
    g_dixon_verbose_level = 0;
    g_field_equation_reduction = 0;
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(4);
    omp_set_max_active_levels(1);
#endif
    check_root_scheduling();
    check_field(101, 1, 0);
    check_field(2, 3, 1024);
    check_field(3, 2, 0);
    check_field(2, 1, 0);
    check_field(2, 4, 0);
    check_field(2, 8, 0);
    check_byte_packing(4);
    check_byte_packing(8);
    return 0;
}

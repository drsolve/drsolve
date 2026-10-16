/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "fq_sparse_interpolation.h"

static void test_random_polynomial(void)
{
    printf("\n========== Testing Sparse.pdf Derivative Interpolation ==========\n");

    slong n = 3;
    slong T = 8000;
    slong d = 1000;
    mp_limb_t p = n_nextprime(2 * (n + 2) * T * T * d, 1);

    printf("Parameters: n=%ld, T=%ld, d=%ld\n", n, T, d);
    printf("Using prime: %lu\n", p);

    nmod_t mod;
    nmod_init(&mod, p);

    char** vars = (char**) malloc((size_t) n * sizeof(char*));
    for (slong i = 0; i < n; i++) {
        vars[i] = (char*) malloc(10);
        sprintf(vars[i], "x%ld", i);
    }

    nmod_mpoly_ctx_t mctx;
    nmod_mpoly_ctx_init(mctx, n, ORD_LEX, p);

    nmod_mpoly_t f, g;
    nmod_mpoly_init(f, mctx);
    nmod_mpoly_init(g, mctx);

    myrandpoly(f, n, T, d, mod, mctx);

    slong* max_degs = (slong*) malloc((size_t) n * sizeof(slong));
    poly_max_degrees_per_var(max_degs, f, mctx);

    if (!DerivativeSparseInterpolatePolynomial(g, f, n, mod, max_degs, mctx)) {
        printf("Derivative sparse interpolation failed.\n");
    }

    nmod_mpoly_sort_terms(f, mctx);
    nmod_mpoly_sort_terms(g, mctx);

    printf("Recovered polynomial with %ld terms\n", nmod_mpoly_length(g, mctx));
    if (nmod_mpoly_equal(f, g, mctx)) {
        printf("Success: Polynomials match!\n");
    } else {
        printf("Error: Polynomials do not match!\n");
    }

    free(max_degs);
    nmod_mpoly_clear(f, mctx);
    nmod_mpoly_clear(g, mctx);
    nmod_mpoly_ctx_clear(mctx);

    for (slong i = 0; i < n; i++) {
        free(vars[i]);
    }
    free(vars);
}

int test_sparse_interpolation(void)
{
    printf("Code version: probe-based Sparse.pdf interpolation without old method\n");

    printf("========== Testing Polynomial Matrix Determinant ==========\n");

    mp_limb_t p = 65537;
    slong n = 3;
    slong k = 2;

    char** vars = (char**) malloc((size_t) n * sizeof(char*));
    for (slong i = 0; i < n; i++) {
        vars[i] = (char*) malloc(10);
        sprintf(vars[i], "x%ld", i);
    }

    nmod_mpoly_ctx_t mctx;
    nmod_mpoly_ctx_init(mctx, n, ORD_LEX, p);

    poly_mat_t A;
    poly_mat_init(&A, k, k, mctx);

    nmod_mpoly_t entry;
    nmod_mpoly_init(entry, mctx);
    ulong* exp = (ulong*) calloc((size_t) n, sizeof(ulong));

    nmod_mpoly_zero(entry, mctx);
    exp[0] = 0; exp[1] = 1; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 3, exp, mctx);
    exp[0] = 3; exp[1] = 1; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 22, exp, mctx);
    poly_mat_entry_set(&A, 0, 0, entry, mctx);

    nmod_mpoly_zero(entry, mctx);
    exp[0] = 1; exp[1] = 2; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 64, exp, mctx);
    poly_mat_entry_set(&A, 0, 1, entry, mctx);

    nmod_mpoly_zero(entry, mctx);
    exp[0] = 2; exp[1] = 0; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 61, exp, mctx);
    exp[0] = 2; exp[1] = 1; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 91, exp, mctx);
    poly_mat_entry_set(&A, 1, 0, entry, mctx);

    nmod_mpoly_zero(entry, mctx);
    exp[0] = 1; exp[1] = 3; exp[2] = 0;
    nmod_mpoly_push_term_ui_ui(entry, 87, exp, mctx);
    exp[0] = 0; exp[1] = 2; exp[2] = 4;
    nmod_mpoly_push_term_ui_ui(entry, 26, exp, mctx);
    exp[0] = 0; exp[1] = 3; exp[2] = 2;
    nmod_mpoly_push_term_ui_ui(entry, 89, exp, mctx);
    poly_mat_entry_set(&A, 1, 1, entry, mctx);

    free(exp);
    nmod_mpoly_clear(entry, mctx);

    nmod_mpoly_t det_poly, actual_det;
    nmod_mpoly_init(det_poly, mctx);
    nmod_mpoly_init(actual_det, mctx);

    poly_mat_det(actual_det, &A, mctx);
    ComputePolyMatrixDet(det_poly, &A, n, p, mctx);

    nmod_mpoly_sort_terms(det_poly, mctx);
    nmod_mpoly_sort_terms(actual_det, mctx);

    printf("\n========== FINAL RESULTS ==========\n");
    printf("Computed determinant:\n");
    nmod_mpoly_print_pretty(det_poly, (const char**) vars, mctx);
    printf("\n\nActual determinant:\n");
    nmod_mpoly_print_pretty(actual_det, (const char**) vars, mctx);
    printf("\n");

    if (nmod_mpoly_equal(det_poly, actual_det, mctx)) {
        printf("\nSuccess: Determinants match!\n");
    } else {
        printf("\nError: Determinants do not match!\n");
    }

    nmod_mpoly_clear(det_poly, mctx);
    nmod_mpoly_clear(actual_det, mctx);
    poly_mat_clear(&A, mctx);
    nmod_mpoly_ctx_clear(mctx);

    for (slong i = 0; i < n; i++) {
        free(vars[i]);
    }
    free(vars);

    test_random_polynomial();

    fq_sparse_interpolation_cleanup();
    flint_cleanup_master();
    return 0;
}

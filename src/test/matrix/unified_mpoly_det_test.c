/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_mpoly_det.h"

void test_unified_mpoly_det(void) {
    printf("\n=== Testing Unified Polynomial Matrix Determinant ===\n");
    
    /* Initialize field context for GF(2^8) */
    fq_nmod_ctx_t fq_ctx;
    nmod_poly_t mod;
    nmod_poly_init(mod, 2);
    
    /* Set modulus for GF(2^8): x^8 + x^4 + x^3 + x^2 + 1 */
    nmod_poly_set_coeff_ui(mod, 0, 1);
    nmod_poly_set_coeff_ui(mod, 2, 1);
    nmod_poly_set_coeff_ui(mod, 3, 1);
    nmod_poly_set_coeff_ui(mod, 4, 1);
    nmod_poly_set_coeff_ui(mod, 8, 1);
    
    fq_nmod_ctx_init_modulus(fq_ctx, mod, "a");
    nmod_poly_clear(mod);
    
    /* Create field context */
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, fq_ctx);
    
    /* Create multivariate context for 2 variables */
    unified_mpoly_ctx_t ctx = unified_mpoly_ctx_init(2, ORD_LEX, &field_ctx);
    
    /* Test 3x3 determinant */
    slong size = 3;
    unified_mpoly_t **mat = unified_mpoly_mat_init(size, size, ctx);
    
    /* Set up a simple test matrix:
     * [x+1,  y,   1]
     * [y,    x,   0]
     * [1,    0,   x+y]
     */
    field_elem_u one;
    field_set_one(&one, field_ctx.field_id, (void*)field_ctx.ctx.fq_ctx);
    
    ulong exp[2];
    
    /* mat[0][0] = x + 1 */
    exp[0] = 1; exp[1] = 0;  /* x */
    unified_mpoly_set_coeff_ui(mat[0][0], &one, exp);
    exp[0] = 0; exp[1] = 0;  /* 1 */
    unified_mpoly_set_coeff_ui(mat[0][0], &one, exp);
    
    /* mat[0][1] = y */
    exp[0] = 0; exp[1] = 1;  /* y */
    unified_mpoly_set_coeff_ui(mat[0][1], &one, exp);
    
    /* mat[0][2] = 1 */
    exp[0] = 0; exp[1] = 0;  /* 1 */
    unified_mpoly_set_coeff_ui(mat[0][2], &one, exp);
    
    /* mat[1][0] = y */
    exp[0] = 0; exp[1] = 1;  /* y */
    unified_mpoly_set_coeff_ui(mat[1][0], &one, exp);
    
    /* mat[1][1] = x */
    exp[0] = 1; exp[1] = 0;  /* x */
    unified_mpoly_set_coeff_ui(mat[1][1], &one, exp);
    
    /* mat[1][2] = 0 */
    unified_mpoly_zero(mat[1][2]);
    
    /* mat[2][0] = 1 */
    exp[0] = 0; exp[1] = 0;  /* 1 */
    unified_mpoly_set_coeff_ui(mat[2][0], &one, exp);
    
    /* mat[2][1] = 0 */
    unified_mpoly_zero(mat[2][1]);
    
    /* mat[2][2] = x + y */
    exp[0] = 1; exp[1] = 0;  /* x */
    unified_mpoly_set_coeff_ui(mat[2][2], &one, exp);
    exp[0] = 0; exp[1] = 1;  /* y */
    unified_mpoly_set_coeff_ui(mat[2][2], &one, exp);
    
    /* Print the matrix */
    const char *vars[] = {"x", "y"};
    printf("Test matrix:\n");
    unified_mpoly_mat_print_pretty(mat, size, size, vars);
    
    /* Compute determinant */
    unified_mpoly_t det = unified_mpoly_init(ctx);
    
    /* Test sequential computation */
    printf("\nSequential computation:\n");
    clock_t start = clock();
    compute_unified_mpoly_det(det, mat, size, ctx, 0);
    clock_t end = clock();
    double seq_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Determinant = ");
    unified_mpoly_print_pretty(det, vars);
    printf("\n");
    printf("Time: %.6f seconds\n", seq_time);
    
    /* Test parallel computation (if available) */
    #ifdef _OPENMP
    printf("\nParallel computation:\n");
    unified_mpoly_t det_parallel = unified_mpoly_init(ctx);
    
    start = clock();
    compute_unified_mpoly_det(det_parallel, mat, size, ctx, 1);
    end = clock();
    double par_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Determinant = ");
    unified_mpoly_print_pretty(det_parallel, vars);
    printf("\n");
    printf("Time: %.6f seconds\n", par_time);
    
    /* Verify results match */
    unified_mpoly_sub(det_parallel, det_parallel, det);
    if (unified_mpoly_is_zero(det_parallel)) {
        printf("Results match!\n");
    } else {
        printf("ERROR: Results don't match!\n");
    }
    
    unified_mpoly_clear(det_parallel);
    #endif
    
    /* Cleanup */
    unified_mpoly_clear(det);
    unified_mpoly_mat_clear(mat, size, size);
    unified_mpoly_ctx_clear(ctx);
    field_ctx_clear(&field_ctx);
    fq_nmod_ctx_clear(fq_ctx);
    
    /* Test with larger matrix */
    printf("\n=== Testing with larger 5x5 matrix ===\n");
    
    /* Reinitialize for prime field test */
    fq_nmod_ctx_init_ui(fq_ctx, 101, 1, "x");
    field_ctx_init(&field_ctx, fq_ctx);
    ctx = unified_mpoly_ctx_init(2, ORD_LEX, &field_ctx);
    
    size = 5;
    mat = unified_mpoly_mat_init(size, size, ctx);
    
    /* Create a random-ish matrix */
    flint_rand_t state;
    flint_rand_init(state);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            /* Add some random terms */
            field_elem_u coeff;
            for (int k = 0; k < 3; k++) {
                exp[0] = n_randint(state, 3);
                exp[1] = n_randint(state, 3);
                coeff.nmod = n_randint(state, 100) + 1;
                unified_mpoly_set_coeff_ui(mat[i][j], &coeff, exp);
            }
        }
    }
    
    printf("Computing 5x5 determinant...\n");
    det = unified_mpoly_init(ctx);
    
    start = clock();
    compute_unified_mpoly_det(det, mat, size, ctx, 1);  /* Use parallel if available */
    end = clock();
    double time_5x5 = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    printf("Determinant computed in %.6f seconds\n", time_5x5);
    printf("Result has %ld terms\n", unified_mpoly_length(det));
    
    /* Cleanup */
    flint_rand_clear(state);
    unified_mpoly_clear(det);
    unified_mpoly_mat_clear(mat, size, size);
    unified_mpoly_ctx_clear(ctx);
    field_ctx_clear(&field_ctx);
    fq_nmod_ctx_clear(fq_ctx);
}

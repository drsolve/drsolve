/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Method reference for the stored-minor expansion path: Gentleman and
 * Johnson. https://dl.acm.org/doi/10.1145/355694.355696
 */
/* unified_mpoly_det.c - Implementation of unified polynomial matrix determinant computation */

#include "unified_mpoly_det.h"

static void debug_print_unified_mpoly_generic(const unified_mpoly_t poly)
{
    slong nvars = (poly != NULL && poly->ctx_ptr != NULL) ? poly->ctx_ptr->nvars : 0;
    if (nvars <= 0) {
        const char *fallback[] = {"x"};
        unified_mpoly_print_pretty(poly, fallback);
        return;
    }

    char **vars = (char **) malloc((size_t) nvars * sizeof(char *));
    if (vars == NULL) {
        printf("<debug-print allocation failed>");
        return;
    }
    for (slong i = 0; i < nvars; i++) {
        vars[i] = (char *) malloc(24);
        if (vars[i] == NULL) {
            for (slong j = 0; j < i; j++) free(vars[j]);
            free(vars);
            printf("<debug-print allocation failed>");
            return;
        }
        snprintf(vars[i], 24, "x%ld", i + 1);
    }

    unified_mpoly_print_pretty(poly, (const char **) vars);

    for (slong i = 0; i < nvars; i++) free(vars[i]);
    free(vars);
}

static int debug_perm_sign(const slong *perm, slong size)
{
    int sign = 1;
    for (slong i = 0; i < size; i++) {
        for (slong j = i + 1; j < size; j++) {
            if (perm[i] > perm[j]) sign = -sign;
        }
    }
    return sign;
}

static slong debug_small_factorial(slong n)
{
    slong r = 1;
    for (slong i = 2; i <= n; i++) r *= i;
    return r;
}

static void debug_dump_leibniz_paths_dfs(unified_mpoly_t **mpoly_matrix,
                                         slong size,
                                         unified_mpoly_ctx_t ctx,
                                         slong row,
                                         slong *perm,
                                         int *used,
                                         slong *path_index)
{
    if (row == size) {
        unified_mpoly_t prefix = unified_mpoly_init(ctx);
        unified_mpoly_t full = unified_mpoly_init(ctx);
        slong prefix_terms;
        slong last_terms;
        int sign = debug_perm_sign(perm, size);

        unified_mpoly_one(prefix);
        for (slong i = 0; i < size - 1; i++) {
            unified_mpoly_mul(prefix, prefix, mpoly_matrix[i][perm[i]]);
        }
        prefix_terms = unified_mpoly_length(prefix);
        last_terms = unified_mpoly_length(mpoly_matrix[size - 1][perm[size - 1]]);

        unified_mpoly_mul(full, prefix, mpoly_matrix[size - 1][perm[size - 1]]);
        if (sign < 0) {
            unified_mpoly_neg(full, full);
        }

        printf("[direct-det path %ld/%ld] perm=[",
               *path_index + 1, debug_small_factorial(size));
        for (slong i = 0; i < size; i++) {
            printf("%ld%s", perm[i], (i + 1 < size) ? "," : "");
        }
        printf("] sign=%c\n", sign > 0 ? '+' : '-');
        printf("  final multiply operand term counts: left=%ld, right=%ld\n",
               prefix_terms, last_terms);
        printf("  path result before summation: terms=%ld\n",
               unified_mpoly_length(full));

        unified_mpoly_clear(prefix);
        unified_mpoly_clear(full);
        (*path_index)++;
        return;
    }

    for (slong col = 0; col < size; col++) {
        if (used[col]) continue;
        used[col] = 1;
        perm[row] = col;
        debug_dump_leibniz_paths_dfs(mpoly_matrix, size, ctx, row + 1, perm, used, path_index);
        used[col] = 0;
    }
}

static void debug_dump_leibniz_paths(unified_mpoly_t **mpoly_matrix,
                                     slong size,
                                     unified_mpoly_ctx_t ctx)
{
    slong *perm = NULL;
    int *used = NULL;
    slong path_index = 0;

    if (g_dixon_verbose_level < 3 || size > 5 || size <= 0) {
        return;
    }

    perm = (slong *) calloc((size_t) size, sizeof(slong));
    used = (int *) calloc((size_t) size, sizeof(int));
    if (perm == NULL || used == NULL) {
        free(perm);
        free(used);
        return;
    }

    printf("\n=== Direct determinant Leibniz path dump (size=%ld) ===\n", size);
    debug_dump_leibniz_paths_dfs(mpoly_matrix, size, ctx, 0, perm, used, &path_index);
    printf("=== End direct determinant Leibniz path dump ===\n");

    free(perm);
    free(used);
}

static slong debug_binomial_slong(slong n, slong k)
{
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    slong r = 1;
    for (slong i = 1; i <= k; i++) {
        r = (r * (n - k + i)) / i;
    }
    return r;
}

static void debug_print_balanced_split_summary(slong size)
{
    slong left = size / 2;
    slong subsets = debug_binomial_slong(size, left);
    if (g_dixon_verbose_level >= 2) {
        printf("  experimental balanced split: top-level left-size=%ld, right-size=%ld, column-subsets=%ld\n",
               left, size - left, subsets);
    }
}

static int debug_subset_sign_from_cols(const int *cols, slong left_size)
{
    slong sum = 0;
    for (slong i = 0; i < left_size; i++) sum += cols[i];
    return (sum & 1) ? -1 : 1;
}

static void build_complement_cols(const int *cols,
                                  slong left_size,
                                  slong size,
                                  int *right_cols)
{
    slong p = 0;
    slong q = 0;
    for (slong c = 0; c < size; c++) {
        if (p < left_size && cols[p] == c) {
            p++;
        } else {
            right_cols[q++] = (int) c;
        }
    }
}

static void fill_submatrix_from_rows_cols(unified_mpoly_t **dst,
                                          unified_mpoly_t **src,
                                          slong row_start,
                                          slong row_count,
                                          const int *cols,
                                          slong col_count)
{
    for (slong i = 0; i < row_count; i++) {
        for (slong j = 0; j < col_count; j++) {
            unified_mpoly_set(dst[i][j], src[row_start + i][cols[j]]);
        }
    }
}

static int next_k_subset_inplace(int *subset, slong k, slong n)
{
    for (slong i = k - 1; i >= 0; i--) {
        if (subset[i] < (int) (n - k + i)) {
            subset[i]++;
            for (slong j = i + 1; j < k; j++) {
                subset[j] = subset[j - 1] + 1;
            }
            return 1;
        }
    }
    return 0;
}

static void compute_unified_mpoly_det_balanced_split_top_level(unified_mpoly_t det_result,
                                                               unified_mpoly_t **mpoly_matrix,
                                                               slong size,
                                                               unified_mpoly_ctx_t ctx)
{
    slong left_size = size / 2;
    slong right_size = size - left_size;
    unified_mpoly_t accum = unified_mpoly_init(ctx);
    unified_mpoly_t term = unified_mpoly_init(ctx);
    unified_mpoly_t left_det = unified_mpoly_init(ctx);
    unified_mpoly_t right_det = unified_mpoly_init(ctx);
    unified_mpoly_t temp = unified_mpoly_init(ctx);
    unified_mpoly_t **left_mat = NULL;
    unified_mpoly_t **right_mat = NULL;
    int *left_cols = NULL;
    int *right_cols = NULL;

    unified_mpoly_zero(accum);

    if (size <= 1) {
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        unified_mpoly_clear(accum);
        unified_mpoly_clear(term);
        unified_mpoly_clear(left_det);
        unified_mpoly_clear(right_det);
        unified_mpoly_clear(temp);
        return;
    }

    left_mat = unified_mpoly_mat_init(left_size, left_size, ctx);
    right_mat = unified_mpoly_mat_init(right_size, right_size, ctx);
    left_cols = (int *) malloc((size_t) left_size * sizeof(int));
    right_cols = (int *) malloc((size_t) right_size * sizeof(int));
    if (left_mat == NULL || right_mat == NULL || left_cols == NULL || right_cols == NULL) {
        if (g_dixon_verbose_level >= 1) {
            printf("  balanced split experimental: allocation failed, falling back to minor expansion\n");
        }
        if (left_mat != NULL) unified_mpoly_mat_clear(left_mat, left_size, left_size);
        if (right_mat != NULL) unified_mpoly_mat_clear(right_mat, right_size, right_size);
        free(left_cols);
        free(right_cols);
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        unified_mpoly_clear(accum);
        unified_mpoly_clear(term);
        unified_mpoly_clear(left_det);
        unified_mpoly_clear(right_det);
        unified_mpoly_clear(temp);
        return;
    }

    for (slong i = 0; i < left_size; i++) left_cols[i] = (int) i;

    do {
        int sign;
        build_complement_cols(left_cols, left_size, size, right_cols);
        fill_submatrix_from_rows_cols(left_mat, mpoly_matrix, 0, left_size, left_cols, left_size);
        fill_submatrix_from_rows_cols(right_mat, mpoly_matrix, left_size, right_size, right_cols, right_size);

        compute_unified_mpoly_det_recursive(left_det, left_mat, left_size, ctx);
        compute_unified_mpoly_det_recursive(right_det, right_mat, right_size, ctx);
        unified_mpoly_mul(term, left_det, right_det);

        sign = debug_subset_sign_from_cols(left_cols, left_size);
        if (sign > 0) {
            unified_mpoly_add(temp, accum, term);
        } else {
            unified_mpoly_sub(temp, accum, term);
        }
        unified_mpoly_set(accum, temp);
    } while (next_k_subset_inplace(left_cols, left_size, size));

    unified_mpoly_set(det_result, accum);

    unified_mpoly_mat_clear(left_mat, left_size, left_size);
    unified_mpoly_mat_clear(right_mat, right_size, right_size);
    free(left_cols);
    free(right_cols);
    unified_mpoly_clear(accum);
    unified_mpoly_clear(term);
    unified_mpoly_clear(left_det);
    unified_mpoly_clear(right_det);
    unified_mpoly_clear(temp);
}

/* ============================================================================
   RECURSIVE DETERMINANT COMPUTATION WITH UNIFIED INTERFACE
   ============================================================================ */

/* Hand-optimized 3x3 determinant for unified_mpoly */
static void compute_det_3x3_unified(unified_mpoly_t det, 
                                   unified_mpoly_t **m,
                                   unified_mpoly_ctx_t ctx) {
    unified_mpoly_t t1, t2, t3, t4, t5, t6, sum;
    /*
    const char *vars[] = {"x", "y", "z", "u", "v", "w"};
    for (int i = 0; i < 3; i++) {
        printf("[ ");
        for (int j = 0; j < 3; j++) {
            unified_mpoly_print_pretty(m[i][j], vars);
            if (j < 2) printf(",\n  ");
        }
        printf(" ]");
        if (i < 2) printf(",\n");
        printf("\n");
    }
    */
    
    /* Initialize temporaries */
    t1 = unified_mpoly_init(ctx);
    t2 = unified_mpoly_init(ctx);
    t3 = unified_mpoly_init(ctx);
    t4 = unified_mpoly_init(ctx);
    t5 = unified_mpoly_init(ctx);
    t6 = unified_mpoly_init(ctx);
    sum = unified_mpoly_init(ctx);
    
    /* Compute 6 products in parallel if beneficial */
    #ifdef _OPENMP
    #pragma omp parallel sections if(omp_get_max_threads() > 2)
    {
        #pragma omp section
        {
            unified_mpoly_mul(t1, m[1][1], m[2][2]);
            unified_mpoly_mul(t1, m[0][0], t1);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t2, m[1][2], m[2][0]);
            unified_mpoly_mul(t2, m[0][1], t2);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t3, m[1][0], m[2][1]);
            unified_mpoly_mul(t3, m[0][2], t3);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t4, m[1][0], m[2][2]);
            unified_mpoly_mul(t4, m[0][1], t4);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t5, m[1][1], m[2][0]);
            unified_mpoly_mul(t5, m[0][2], t5);
        }
        #pragma omp section
        {
            unified_mpoly_mul(t6, m[1][2], m[2][1]);
            unified_mpoly_mul(t6, m[0][0], t6);
        }
    }
    #else
    /* Sequential computation when OpenMP is not available */
    unified_mpoly_mul(t1, m[1][1], m[2][2]);
    unified_mpoly_mul(t1, m[0][0], t1);
    
    unified_mpoly_mul(t2, m[1][2], m[2][0]);
    unified_mpoly_mul(t2, m[0][1], t2);
    
    unified_mpoly_mul(t3, m[1][0], m[2][1]);
    unified_mpoly_mul(t3, m[0][2], t3);
    
    unified_mpoly_mul(t4, m[1][0], m[2][2]);
    unified_mpoly_mul(t4, m[0][1], t4);
    
    unified_mpoly_mul(t5, m[1][1], m[2][0]);
    unified_mpoly_mul(t5, m[0][2], t5);
    
    unified_mpoly_mul(t6, m[1][2], m[2][1]);
    unified_mpoly_mul(t6, m[0][0], t6);
    #endif
    
    /* Sum with signs: det = t1 + t2 + t3 - t4 - t5 - t6 */
    unified_mpoly_add(sum, t1, t2);
    unified_mpoly_add(sum, sum, t3);
    unified_mpoly_sub(sum, sum, t4);
    unified_mpoly_sub(sum, sum, t5);
    unified_mpoly_sub(det, sum, t6);
    
    /* Cleanup */
    unified_mpoly_clear(t1);
    unified_mpoly_clear(t2);
    unified_mpoly_clear(t3);
    unified_mpoly_clear(t4);
    unified_mpoly_clear(t5);
    unified_mpoly_clear(t6);
    unified_mpoly_clear(sum);
}

static int unified_mpoly_find_nonzero_pivot(unified_mpoly_t **mat,
                                            slong size,
                                            slong start)
{
    for (slong i = start; i < size; i++) {
        for (slong j = start; j < size; j++) {
            if (!unified_mpoly_is_zero(mat[i][j])) {
                return (int) i;
            }
        }
    }
    return -1;
}

static int unified_mpoly_find_nonzero_pivot_col(unified_mpoly_t **mat,
                                                slong size,
                                                slong row,
                                                slong start_col)
{
    for (slong j = start_col; j < size; j++) {
        if (!unified_mpoly_is_zero(mat[row][j])) {
            return (int) j;
        }
    }
    return -1;
}

int compute_unified_mpoly_det_bareiss(unified_mpoly_t det_result,
                                     unified_mpoly_t **mpoly_matrix,
                                     slong size,
                                     unified_mpoly_ctx_t ctx) {
    if (size <= 0) {
        unified_mpoly_one(det_result);
        return 1;
    }

    if (size == 1) {
        unified_mpoly_set(det_result, mpoly_matrix[0][0]);
        return 1;
    }

    if (size == 2) {
        unified_mpoly_t ad, bc;
        ad = unified_mpoly_init(ctx);
        bc = unified_mpoly_init(ctx);
        unified_mpoly_mul(ad, mpoly_matrix[0][0], mpoly_matrix[1][1]);
        unified_mpoly_mul(bc, mpoly_matrix[0][1], mpoly_matrix[1][0]);
        unified_mpoly_sub(det_result, ad, bc);
        unified_mpoly_clear(ad);
        unified_mpoly_clear(bc);
        return 1;
    }

    unified_mpoly_t **work = unified_mpoly_mat_init(size, size, ctx);
    unified_mpoly_t prev_pivot = unified_mpoly_init(ctx);
    unified_mpoly_t pivot = unified_mpoly_init(ctx);
    unified_mpoly_t zero_poly = unified_mpoly_init(ctx);
    int sign = 1;
    int ok = 1;

    unified_mpoly_zero(zero_poly);
    unified_mpoly_one(prev_pivot);

    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            unified_mpoly_set(work[i][j], mpoly_matrix[i][j]);
        }
    }

    for (slong k = 0; k < size - 1; k++) {
        slong pivot_row = k;
        slong pivot_col = k;

        if (unified_mpoly_is_zero(work[k][k])) {
            int row_found = unified_mpoly_find_nonzero_pivot(work, size, k);
            if (row_found < 0) {
                unified_mpoly_zero(det_result);
                goto cleanup;
            }
            pivot_row = (slong) row_found;
            pivot_col = (slong) unified_mpoly_find_nonzero_pivot_col(work, size, pivot_row, k);
            if (pivot_col < 0) {
                unified_mpoly_zero(det_result);
                goto cleanup;
            }
            if (pivot_row != k) {
                unified_mpoly_t *tmp_row = work[k];
                work[k] = work[pivot_row];
                work[pivot_row] = tmp_row;
                sign = -sign;
            }
            if (pivot_col != k) {
                for (slong i = 0; i < size; i++) {
                    unified_mpoly_swap(work[i][k], work[i][pivot_col]);
                }
                sign = -sign;
            }
        }

        unified_mpoly_set(pivot, work[k][k]);
        if (unified_mpoly_is_zero(pivot)) {
            unified_mpoly_zero(det_result);
            goto cleanup;
        }

        if (unified_mpoly_is_zero(prev_pivot) && k > 0) {
            unified_mpoly_zero(det_result);
            ok = 0;
            goto cleanup;
        }

        #ifdef _OPENMP
        if ((size - (k + 1)) >= 2 && omp_get_max_threads() > 1) {
            volatile int parallel_failed = 0;
            #pragma omp parallel for collapse(2) schedule(static)
            for (slong i = k + 1; i < size; i++) {
                for (slong j = k + 1; j < size; j++) {
                    unified_mpoly_t t1 = unified_mpoly_init(ctx);
                    unified_mpoly_t t2 = unified_mpoly_init(ctx);
                    unified_mpoly_t numer = unified_mpoly_init(ctx);
                    unified_mpoly_t quotient = unified_mpoly_init(ctx);

                    unified_mpoly_mul(t1, work[i][j], work[k][k]);
                    unified_mpoly_mul(t2, work[i][k], work[k][j]);
                    unified_mpoly_sub(numer, t1, t2);

                    if (k == 0) {
                        unified_mpoly_set(work[i][j], numer);
                    } else {
                        if (!unified_mpoly_divexact(quotient, numer, prev_pivot)) {
                            parallel_failed = 1;
                        } else {
                            unified_mpoly_set(work[i][j], quotient);
                        }
                    }

                    unified_mpoly_clear(t1);
                    unified_mpoly_clear(t2);
                    unified_mpoly_clear(numer);
                    unified_mpoly_clear(quotient);
                }
            }
            if (parallel_failed) {
                ok = 0;
                goto cleanup;
            }
        } else
        #endif
        {
            unified_mpoly_t t1 = unified_mpoly_init(ctx);
            unified_mpoly_t t2 = unified_mpoly_init(ctx);
            unified_mpoly_t numer = unified_mpoly_init(ctx);
            unified_mpoly_t quotient = unified_mpoly_init(ctx);

            for (slong i = k + 1; i < size; i++) {
                for (slong j = k + 1; j < size; j++) {
                    unified_mpoly_mul(t1, work[i][j], work[k][k]);
                    unified_mpoly_mul(t2, work[i][k], work[k][j]);
                    unified_mpoly_sub(numer, t1, t2);

                    if (k == 0) {
                        unified_mpoly_set(work[i][j], numer);
                    } else {
                        if (!unified_mpoly_divexact(quotient, numer, prev_pivot)) {
                            ok = 0;
                            unified_mpoly_clear(t1);
                            unified_mpoly_clear(t2);
                            unified_mpoly_clear(numer);
                            unified_mpoly_clear(quotient);
                            goto cleanup;
                        }
                        unified_mpoly_set(work[i][j], quotient);
                    }
                }
            }

            unified_mpoly_clear(t1);
            unified_mpoly_clear(t2);
            unified_mpoly_clear(numer);
            unified_mpoly_clear(quotient);
        }

        for (slong i = k + 1; i < size; i++) {
            unified_mpoly_zero(work[i][k]);
            unified_mpoly_zero(work[k][i]);
        }

        unified_mpoly_set(prev_pivot, pivot);
    }

    unified_mpoly_set(det_result, work[size - 1][size - 1]);
    if (sign < 0) {
        unified_mpoly_sub(det_result, zero_poly, det_result);
    }

cleanup:
    unified_mpoly_clear(prev_pivot);
    unified_mpoly_clear(pivot);
    unified_mpoly_clear(zero_poly);
    unified_mpoly_mat_clear(work, size, size);
    return ok;
}

/* Recursive determinant computation using unified_mpoly interface */
void compute_unified_mpoly_det_recursive(unified_mpoly_t det_result, 
                                        unified_mpoly_t **mpoly_matrix, 
                                        slong size, 
                                        unified_mpoly_ctx_t ctx) {
    
    
    if (size <= 0) {
        unified_mpoly_one(det_result);
        return;
    }
    
    if (size == 1) {
        unified_mpoly_set(det_result, mpoly_matrix[0][0]);
        return;
    }
    
    if (size == 2) {
        unified_mpoly_t ad, bc;
        ad = unified_mpoly_init(ctx);
        bc = unified_mpoly_init(ctx);
        
        /* det = a*d - b*c */
        unified_mpoly_mul(ad, mpoly_matrix[0][0], mpoly_matrix[1][1]);
        unified_mpoly_mul(bc, mpoly_matrix[0][1], mpoly_matrix[1][0]);
        unified_mpoly_sub(det_result, ad, bc);
        
        unified_mpoly_clear(ad);
        unified_mpoly_clear(bc);
        return;
    }
    
    if (size == 3) {
        compute_det_3x3_unified(det_result, mpoly_matrix, ctx);
        return;
    }
    
    /* General case: Laplace expansion along first row */
    unified_mpoly_zero(det_result);
    
    unified_mpoly_t temp_result, cofactor, subdet;
    temp_result = unified_mpoly_init(ctx);
    cofactor = unified_mpoly_init(ctx);
    subdet = unified_mpoly_init(ctx);
    
    /* Allocate submatrix */
    unified_mpoly_t **submatrix = (unified_mpoly_t**) malloc((size-1) * sizeof(unified_mpoly_t*));
    for (slong i = 0; i < size-1; i++) {
        submatrix[i] = (unified_mpoly_t*) malloc((size-1) * sizeof(unified_mpoly_t));
        for (slong j = 0; j < size-1; j++) {
            submatrix[i][j] = unified_mpoly_init(ctx);
        }
    }
    
    /* Laplace expansion along the first row */
    for (slong col = 0; col < size; col++) {
        /* Skip if the element is zero */
        int is_zero = unified_mpoly_is_zero(mpoly_matrix[0][col]);
        
        if (is_zero) {
            continue;
        }
        
        /* Build submatrix by removing row 0 and column col */
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    unified_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j]);
                    sub_j++;
                }
            }
        }
        
        /* Recursive computation of minor */
        compute_unified_mpoly_det_recursive(subdet, submatrix, size-1, ctx);
        
        /* Compute cofactor = (-1)^col * matrix[0][col] * subdet */
        unified_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet);
        
        /* Add or subtract based on sign */
        if (col % 2 == 0) {
            unified_mpoly_add(temp_result, det_result, cofactor);
        } else {
            unified_mpoly_sub(temp_result, det_result, cofactor);
        }
        unified_mpoly_set(det_result, temp_result);
    }
    
    /* Cleanup submatrix */
    for (slong i = 0; i < size-1; i++) {
        for (slong j = 0; j < size-1; j++) {
            unified_mpoly_clear(submatrix[i][j]);
        }
        free(submatrix[i]);
    }
    free(submatrix);
    
    unified_mpoly_clear(temp_result);
    unified_mpoly_clear(cofactor);
    unified_mpoly_clear(subdet);

}

/* Layered minor DP, using colex column-subset ranks and adjacent layers.
 * The final expansion is parallelized by cofactor, followed by a tree sum.
 * The entry limit excludes arithmetic temporaries and matrix views.
 */
#ifdef DRSOLVE_DET_TESTING
/* Test-only observation; normal builds contain no callbacks. */
void drsolve_det_test_event(int event, slong size);
#endif
static int compute_unified_mpoly_det_layered_dp(unified_mpoly_t result, unified_mpoly_t **matrix,
                       slong size, unified_mpoly_ctx_t ctx, int use_parallel, slong limit)
{
    ulong choose[FLINT_BITS][FLINT_BITS] = {{0}};
    unified_mpoly_t *previous = NULL, *current = NULL;
    slong previous_count = 0, current_count = 0;
    slong peak = 0;
    int ok = 0;

    if (size <= 3 || size >= FLINT_BITS || limit <= 0) return 0;

    /* Saturation avoids overflow even when the requested matrix is far too
     * large. All ranks used after this preflight fit within the entry limit. */
    for (slong n = 0; n <= size; n++) {
        choose[n][0] = 1;
        for (slong k = 1; k <= n; k++) {
            ulong a = choose[n - 1][k - 1], b = choose[n - 1][k];
            choose[n][k] = (a > (ulong) limit || b > (ulong) limit - a)
                            ? (ulong) limit + 1 : a + b;
        }
    }
    for (slong k = 1; k <= size; k++) {
        ulong count = k == size ? (ulong) size : choose[size][k];
        ulong prev = k == 1 ? 0 : choose[size][k - 1];
        if (count > (ulong) limit || prev > (ulong) limit - count)
            return 0;
        if (count > (size_t) -1 / sizeof(unified_mpoly_t)) return 0;
        if ((slong) (count + prev) > peak) peak = (slong) (count + prev);
    }

    /* Keep cheap sparse expansions on the demand-driven path. The product
     * of row nonzero counts bounds the number of recursive branches; stop
     * counting once it already exceeds the full DP multiplication count. */
    {
        double dp_work = 0, recursive_work = 0, branches = 1;
        for (slong k = 2; k <= size; k++) dp_work += k * (double) choose[size][k];
        for (slong row = 0; row < size - 3; row++) {
            slong nonzero = 0;
            for (slong col = 0; col < size; col++)
                if (!unified_mpoly_is_zero(matrix[row][col])) nonzero++;
            branches *= FLINT_MIN(nonzero, size - row);
            recursive_work += branches;
            if (recursive_work >= dp_work) break;
        }
        if (recursive_work + 12 * branches < dp_work) return 0;
    }

    if (g_dixon_verbose_level >= 2)
        printf("  determinant layered DP (%s): size=%ld, peak entries=%ld, limit=%ld, parallel=%s\n",
               "unified", size, peak, limit, use_parallel ? "yes" : "no");

#ifdef DRSOLVE_DET_TESTING
    drsolve_det_test_event(0, size);
#endif

    for (slong k = 1; k <= size; k++) {
        slong count = k == size ? size : (slong) choose[size][k];
        int failed = 0;
        current = malloc((size_t) count * sizeof(unified_mpoly_t));
        if (current == NULL) goto cleanup;
        for (current_count = 0; current_count < count; current_count++) {
            if (!((current[current_count] = unified_mpoly_init(ctx)) != NULL)) goto cleanup;
        }

        if (k == 1) {
            for (slong i = 0; i < count; i++)
                unified_mpoly_set(current[i], matrix[size - 1][i]);
        } else {
            /* One writer per output, immutable previous layer, implicit
             * barrier before freeing it. No hash locks or duplicate work. */
#ifdef _OPENMP
            #pragma omp parallel if(use_parallel && count > 1 && !omp_in_parallel()) num_threads(FLINT_MIN(count, omp_get_max_threads())) reduction(|:failed)
#endif
            {
                unified_mpoly_t product, sum;
                int product_ok = ((product = unified_mpoly_init(ctx)) != NULL);
                int sum_ok = ((sum = unified_mpoly_init(ctx)) != NULL);
                if (!product_ok || !sum_ok) failed = 1;
#ifdef _OPENMP
                #pragma omp for schedule(dynamic, 1)
#endif
                for (slong index = 0; index < count; index++) {
                    slong cols[FLINT_BITS];
                    ulong prefix[FLINT_BITS], suffix[FLINT_BITS];
                    ulong rank = (ulong) index;
                    slong col = size - 1;
                    if (!product_ok || !sum_ok) continue;

                    /* There is only one root state. Parallelize its expansion
                     * terms instead of serializing all n large products in
                     * one worker. The colex rank of the (n-1)-subset missing
                     * column index is n-1-index. For n>=4, these n outputs and
                     * n inputs fit below the middle-layer entry peak. */
                    if (k == size) {
                        if (unified_mpoly_is_zero(matrix[0][index]) ||
                            unified_mpoly_is_zero(previous[size - 1 - index])) continue;
#ifdef DRSOLVE_DET_TESTING
                        drsolve_det_test_event(1, size);
#endif
                        unified_mpoly_mul(current[index], matrix[0][index], previous[size - 1 - index]);
                        if (index & 1) {
                            unified_mpoly_zero(product);
                            unified_mpoly_sub(current[index], product, current[index]);
                        }
                        continue;
                    }

                    /* Unrank sum_i C(cols[i], i+1) in O(n) time. */
                    for (slong j = k; j > 0; j--) {
                        while (choose[col][j] > rank) col--;
                        cols[j - 1] = col;
                        rank -= choose[col][j];
                        col--;
                    }
                    prefix[0] = 0;
                    for (slong j = 0; j < k; j++)
                        prefix[j + 1] = prefix[j] + choose[cols[j]][j + 1];
                    suffix[k] = 0;
                    for (slong j = k - 1; j > 0; j--)
                        suffix[j] = suffix[j + 1] + choose[cols[j]][j];

                    unified_mpoly_zero(current[index]);
                    for (slong j = 0; j < k; j++) {
                        ulong child = prefix[j] + suffix[j + 1];
                        if (unified_mpoly_is_zero(matrix[size - k][cols[j]]) ||
                            unified_mpoly_is_zero(previous[child])) continue;
                        unified_mpoly_mul(product, matrix[size - k][cols[j]], previous[child]);
                        if (j & 1)
                            unified_mpoly_sub(sum, current[index], product);
                        else
                            unified_mpoly_add(sum, current[index], product);
                        unified_mpoly_swap(current[index], sum);
                    }
                }
                /* The worksharing barrier above publishes all root terms.
                 * A balanced reduction keeps large sums parallel and avoids
                 * repeatedly merging one term into an ever-growing prefix. */
                if (k == size) {
                    for (slong stride = 1; stride < count; stride *= 2) {
#ifdef _OPENMP
                        #pragma omp for schedule(static)
#endif
                        for (slong left = 0; left < count; left += 2 * stride) {
                            if (!product_ok || !sum_ok || left + stride >= count) continue;
                            unified_mpoly_add(sum, current[left], current[left + stride]);
                            unified_mpoly_swap(current[left], sum);
                        }
                    }
                }
                if (product_ok) unified_mpoly_clear(product);
                if (sum_ok) unified_mpoly_clear(sum);
            }
        }
        if (failed) goto cleanup;
        for (slong i = 0; i < previous_count; i++) unified_mpoly_clear(previous[i]);
        free(previous);
        previous = current;
        previous_count = current_count;
        current = NULL;
        current_count = 0;
    }
    unified_mpoly_swap(result, previous[0]);
    ok = 1;

cleanup:
    for (slong i = 0; i < current_count; i++) unified_mpoly_clear(current[i]);
    for (slong i = 0; i < previous_count; i++) unified_mpoly_clear(previous[i]);
    free(current);
    free(previous);
    return ok;
}

/* Method 0 backend. If a complete layer does not fit, expand one row and
 * retry DP on each child. Siblings execute sequentially, so their DP storage
 * never multiplies the entry budget; each child may still use layer threads.
 * The submatrix is a shallow, read-only view of the input polynomials. */
static void compute_unified_mpoly_det_minor(unified_mpoly_t result, unified_mpoly_t **matrix,
                         slong size, unified_mpoly_ctx_t ctx, int use_parallel, slong limit)
{
    unified_mpoly_t **rows = NULL, *entries = NULL;
    unified_mpoly_t accum, child, product, sum;
    int accum_ok, child_ok, product_ok, sum_ok;
    size_t width;

    if (size <= 3 || limit <= 0) {
        compute_unified_mpoly_det_recursive(result, matrix, size, ctx);
        return;
    }
    if (compute_unified_mpoly_det_layered_dp(result, matrix, size, ctx, use_parallel, limit)) return;

    width = (size_t) (size - 1);
    if (width > (size_t) -1 / sizeof(*rows) ||
        width > (size_t) -1 / sizeof(*entries) / width) {
        compute_unified_mpoly_det_recursive(result, matrix, size, ctx);
        return;
    }
    rows = malloc(width * sizeof(*rows));
    entries = malloc(width * width * sizeof(*entries));
    if (rows == NULL || entries == NULL) {
        free(rows);
        free(entries);
        compute_unified_mpoly_det_recursive(result, matrix, size, ctx);
        return;
    }
    for (size_t i = 0; i < width; i++) rows[i] = entries + i * width;
    accum_ok = ((accum = unified_mpoly_init(ctx)) != NULL);
    child_ok = ((child = unified_mpoly_init(ctx)) != NULL);
    product_ok = ((product = unified_mpoly_init(ctx)) != NULL);
    sum_ok = ((sum = unified_mpoly_init(ctx)) != NULL);
    if (accum_ok && child_ok && product_ok && sum_ok) {
        unified_mpoly_zero(accum);
        for (slong col = 0; col < size; col++) {
            if (unified_mpoly_is_zero(matrix[0][col])) continue;
            for (slong i = 1; i < size; i++) {
                slong dst = 0;
                for (slong j = 0; j < size; j++) {
                    if (j == col) continue;
                    memcpy(&rows[i - 1][dst++], &matrix[i][j], sizeof(*entries));
                }
            }
            compute_unified_mpoly_det_minor(child, rows, size - 1, ctx, use_parallel, limit);
            if (unified_mpoly_is_zero(child)) continue;
            unified_mpoly_mul(product, matrix[0][col], child);
            if (col & 1)
                unified_mpoly_sub(sum, accum, product);
            else
                unified_mpoly_add(sum, accum, product);
            unified_mpoly_swap(accum, sum);
        }
        unified_mpoly_swap(result, accum);
    } else {
        compute_unified_mpoly_det_recursive(result, matrix, size, ctx);
    }
    if (accum_ok) unified_mpoly_clear(accum);
    if (child_ok) unified_mpoly_clear(child);
    if (product_ok) unified_mpoly_clear(product);
    if (sum_ok) unified_mpoly_clear(sum);
    free(entries);
    free(rows);
}

/* ============================================================================
   PARALLEL DETERMINANT COMPUTATION WITH NESTED PARALLELISM
   ============================================================================ */

/* Initialize nested parallelism */
void init_nested_parallelism(void) {
    #ifdef _OPENMP
    omp_set_nested(1);  /* Enable nested parallelism */
    omp_set_max_active_levels(2);  /* Allow 2 levels of parallelism */
    #endif
}

/* Parallel determinant computation with proper nested parallelism */
void compute_unified_mpoly_det_parallel(unified_mpoly_t det_result, 
                                       unified_mpoly_t **mpoly_matrix, 
                                       slong size, 
                                       unified_mpoly_ctx_t ctx,
                                       slong depth) {
    /* For deep recursion or small matrices, use sequential */
    if (size < PARALLEL_THRESHOLD || depth >= MAX_PARALLEL_DEPTH) {
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        return;
    }
    
    if (size <= 3) {
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        return;
    }
    
    unified_mpoly_zero(det_result);
    
    /* Count non-zero entries in first row */
    slong nonzero_count = 0;
    for (slong col = 0; col < size; col++) {
        if (!unified_mpoly_is_zero(mpoly_matrix[0][col])) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count < 2) {
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        return;
    }
    
    /* Allocate space for partial results */
    unified_mpoly_t *partial_results = (unified_mpoly_t*) malloc(size * sizeof(unified_mpoly_t));
    for (slong i = 0; i < size; i++) {
        partial_results[i] = unified_mpoly_init(ctx);
        unified_mpoly_zero(partial_results[i]);
    }
    
    /* Determine parallelism strategy based on depth */
    #ifdef _OPENMP
    if (depth == 0) {
        /* First level: use parallel for with nested parallelism enabled */
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, omp_get_max_threads()))
        for (slong col = 0; col < size; col++) {
            if (unified_mpoly_is_zero(mpoly_matrix[0][col])) {
                continue;
            }
            
            unified_mpoly_t cofactor = unified_mpoly_init(ctx);
            unified_mpoly_t subdet = unified_mpoly_init(ctx);
            
            /* Create submatrix */
            unified_mpoly_t **submatrix = (unified_mpoly_t**) malloc((size-1) * sizeof(unified_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (unified_mpoly_t*) malloc((size-1) * sizeof(unified_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    submatrix[i][j] = unified_mpoly_init(ctx);
                }
            }
            
            /* Fill submatrix */
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        unified_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j]);
                        sub_j++;
                    }
                }
            }
            
            /* Recursive call - this will use nested parallelism at depth 1 */
            compute_unified_mpoly_det_parallel(subdet, submatrix, size-1, ctx, depth+1);
            
            /* Compute cofactor */
            unified_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet);
            
            /* Store with sign */
            if (col % 2 == 0) {
                unified_mpoly_set(partial_results[col], cofactor);
            } else {
                unified_mpoly_neg(partial_results[col], cofactor);
            }
            
            /* Cleanup */
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    unified_mpoly_clear(submatrix[i][j]);
                }
                free(submatrix[i]);
            }
            free(submatrix);
            
            unified_mpoly_clear(cofactor);
            unified_mpoly_clear(subdet);
        }
    } else if (depth == 1 && size >= PARALLEL_THRESHOLD) {
        /* Second level: also use parallel for, but with fewer threads */
        slong max_threads_level2 = FLINT_MAX(1, omp_get_max_threads() / size);
        
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, max_threads_level2)) if(nonzero_count >= 3)
        for (slong col = 0; col < size; col++) {
            if (unified_mpoly_is_zero(mpoly_matrix[0][col])) {
                continue;
            }
            
            unified_mpoly_t cofactor = unified_mpoly_init(ctx);
            unified_mpoly_t subdet = unified_mpoly_init(ctx);
            
            /* Create submatrix */
            unified_mpoly_t **submatrix = (unified_mpoly_t**) malloc((size-1) * sizeof(unified_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (unified_mpoly_t*) malloc((size-1) * sizeof(unified_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    submatrix[i][j] = unified_mpoly_init(ctx);
                }
            }
            
            /* Fill submatrix */
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        unified_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j]);
                        sub_j++;
                    }
                }
            }
            
            /* At depth > 1, use sequential computation */
            compute_unified_mpoly_det_recursive(subdet, submatrix, size-1, ctx);
            
            /* Compute cofactor */
            unified_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet);
            
            /* Store with sign */
            if (col % 2 == 0) {
                unified_mpoly_set(partial_results[col], cofactor);
            } else {
                unified_mpoly_neg(partial_results[col], cofactor);
            }
            
            /* Cleanup */
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    unified_mpoly_clear(submatrix[i][j]);
                }
                free(submatrix[i]);
            }
            free(submatrix);
            
            unified_mpoly_clear(cofactor);
            unified_mpoly_clear(subdet);
        }
    } else
    #endif
    {
        /* Sequential fallback for deeper levels or when OpenMP is not available */
        for (slong col = 0; col < size; col++) {
            if (unified_mpoly_is_zero(mpoly_matrix[0][col])) {
                continue;
            }
            
            unified_mpoly_t cofactor = unified_mpoly_init(ctx);
            unified_mpoly_t subdet = unified_mpoly_init(ctx);
            
            /* Create and fill submatrix */
            unified_mpoly_t **submatrix = (unified_mpoly_t**) malloc((size-1) * sizeof(unified_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (unified_mpoly_t*) malloc((size-1) * sizeof(unified_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    submatrix[i][j] = unified_mpoly_init(ctx);
                }
            }
            
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        unified_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j]);
                        sub_j++;
                    }
                }
            }
            
            /* Sequential computation */
            compute_unified_mpoly_det_recursive(subdet, submatrix, size-1, ctx);
            
            /* Compute cofactor and store */
            unified_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet);
            if (col % 2 == 0) {
                unified_mpoly_set(partial_results[col], cofactor);
            } else {
                unified_mpoly_neg(partial_results[col], cofactor);
            }
            
            /* Cleanup */
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    unified_mpoly_clear(submatrix[i][j]);
                }
                free(submatrix[i]);
            }
            free(submatrix);
            
            unified_mpoly_clear(cofactor);
            unified_mpoly_clear(subdet);
        }
    }
    
    /* Sum results (sequential to avoid race conditions) */
    unified_mpoly_t temp_sum = unified_mpoly_init(ctx);
    
    for (slong col = 0; col < size; col++) {
        if (!unified_mpoly_is_zero(partial_results[col])) {
            unified_mpoly_add(temp_sum, det_result, partial_results[col]);
            unified_mpoly_set(det_result, temp_sum);
        }
        unified_mpoly_clear(partial_results[col]);
    }
    
    unified_mpoly_clear(temp_sum);
    free(partial_results);
}

/* ============================================================================
   CONVENIENCE WRAPPER FUNCTIONS
   ============================================================================ */

/* Main entry point for determinant computation */
void compute_unified_mpoly_det(unified_mpoly_t det_result,
                              unified_mpoly_t **mpoly_matrix,
                              slong size,
                              unified_mpoly_ctx_t ctx,
                              int use_parallel) {
    compute_unified_mpoly_det_with_method(det_result, mpoly_matrix, size, ctx,
                                          use_parallel, 0);
}

void compute_unified_mpoly_det_with_method(unified_mpoly_t det_result,
                                          unified_mpoly_t **mpoly_matrix,
                                          slong size,
                                          unified_mpoly_ctx_t ctx,
                                          int use_parallel,
                                          int method) {
    if (method != DET_METHOD_KRONECKER_NMOD &&
        g_dixon_verbose_level >= 3 && size <= 3) {
        debug_dump_leibniz_paths(mpoly_matrix, size, ctx);
    }

    if (method == DET_METHOD_KRONECKER_NMOD) {
        if (!compute_unified_mpoly_det_bareiss(det_result, mpoly_matrix, size, ctx)) {
            compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        }
        return;
    }

    if (method == DET_METHOD_BALANCED_SPLIT) {
        if (g_dixon_verbose_level >= 2) {
            printf("  experimental determinant method: balanced split Laplace skeleton\n");
            debug_print_balanced_split_summary(size);
            printf("  current implementation note: one top-level balanced split; subdeterminants use the existing recursive backend\n");
        }
        compute_unified_mpoly_det_balanced_split_top_level(det_result, mpoly_matrix, size, ctx);
        return;
    }

    if (method == DET_METHOD_RECURSIVE) {
        compute_unified_mpoly_det_minor(det_result, mpoly_matrix, size, ctx,
                                        use_parallel, g_dixon_det_cache_limit);
        return;
    }
    
    /* Special case: 3x3 matrix with parallel computation available */
    #ifdef _OPENMP
    if (size == 3 && use_parallel && omp_get_max_threads() > 1) {
        compute_det_3x3_unified(det_result, mpoly_matrix, ctx);
        return;
    }
    #endif
    
    if (use_parallel && size >= PARALLEL_THRESHOLD) {
        #ifdef _OPENMP
        static int nested_init = 0;
        if (!nested_init) {
            init_nested_parallelism();
            nested_init = 1;
        }
        if (g_dixon_verbose_level >= 2) {
            printf("  determinant path: parallel Laplace without shared memo (size=%ld)\n", size);
        }
        compute_unified_mpoly_det_parallel(det_result, mpoly_matrix, size, ctx, 0);
        #else
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
        #endif
    } else {
        compute_unified_mpoly_det_recursive(det_result, mpoly_matrix, size, ctx);
    }
}

/* ============================================================================
   HELPER FUNCTIONS FOR MATRIX ALLOCATION/DEALLOCATION
   ============================================================================ */

/* Allocate a unified_mpoly matrix */
unified_mpoly_t** unified_mpoly_mat_init(slong rows, slong cols, unified_mpoly_ctx_t ctx) {
    unified_mpoly_t **mat = (unified_mpoly_t**) malloc(rows * sizeof(unified_mpoly_t*));
    for (slong i = 0; i < rows; i++) {
        mat[i] = (unified_mpoly_t*) malloc(cols * sizeof(unified_mpoly_t));
        for (slong j = 0; j < cols; j++) {
            mat[i][j] = unified_mpoly_init(ctx);
        }
    }
    return mat;
}

/* Free a unified_mpoly matrix */
void unified_mpoly_mat_clear(unified_mpoly_t **mat, slong rows, slong cols) {
    for (slong i = 0; i < rows; i++) {
        for (slong j = 0; j < cols; j++) {
            unified_mpoly_clear(mat[i][j]);
        }
        free(mat[i]);
    }
    free(mat);
}

/* Print a unified_mpoly matrix */
void unified_mpoly_mat_print_pretty(unified_mpoly_t **mat, slong rows, slong cols,
                                   const char **vars) {
    printf("[\n");
    for (slong i = 0; i < rows; i++) {
        printf("  [");
        for (slong j = 0; j < cols; j++) {
            unified_mpoly_print_pretty(mat[i][j], vars);
            if (j < cols - 1) printf(", ");
        }
        printf("]");
        if (i < rows - 1) printf(",");
        printf("\n");
    }
    printf("]\n");
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Method reference for the stored-minor expansion path: Gentleman and
 * Johnson. https://dl.acm.org/doi/10.1145/355694.355696
 */
/*
 * Implementation of optimized polynomial matrix determinant computation
 * Contains all algorithm implementations for various determinant methods
 */

#include "fq_mpoly_mat_det.h"
extern int g_field_equation_reduction;
extern int g_dixon_debug_mode;
extern int g_dixon_verbose_level;
extern slong g_dixon_det_cache_limit;
static inline ulong reduce_exp_field_ui(ulong e, ulong q);
static ulong field_size_q_from_fq_ctx(const fq_nmod_ctx_t ctx);
static void fq_nmod_poly_reduce_field_equation_inplace(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx);
static void fq_nmod_mpoly_reduce_field_equation_inplace(fq_nmod_mpoly_t poly, const fq_nmod_mpoly_ctx_t ctx);
// ============= Timing Utilities Implementation =============
double get_cpu_time(void) {
    return ((double)clock()) / CLOCKS_PER_SEC;
}

timing_info_t start_timing(void) {
    timing_info_t t;
    t.wall_time = get_wall_time();
    t.cpu_time = get_cpu_time();
    return t;
}

timing_info_t end_timing(timing_info_t start) {
    timing_info_t elapsed;
    elapsed.wall_time = get_wall_time() - start.wall_time;
    elapsed.cpu_time = get_cpu_time() - start.cpu_time;
    return elapsed;
}

void print_timing(const char* label, timing_info_t elapsed) {
    DET_PRINT("%s: Wall time: %.6f s, CPU time: %.6f s", 
              label, elapsed.wall_time, elapsed.cpu_time);
    if (elapsed.wall_time > 0) {
        DET_PRINT(" (CPU efficiency: %.1f%%)\n", 
                  (elapsed.cpu_time / elapsed.wall_time) * 100.0);
    } else {
        DET_PRINT("\n");
    }
}

// ============= Polynomial Operation Optimizations Implementation =============

// Custom multiplication for dense polynomials
static inline void poly_mul_dense_optimized(fq_nmod_mpoly_t c, 
                                           const fq_nmod_mpoly_t a,
                                           const fq_nmod_mpoly_t b,
                                           const fq_nmod_mpoly_ctx_t ctx) {
    slong alen = fq_nmod_mpoly_length(a, ctx);
    slong blen = fq_nmod_mpoly_length(b, ctx);
    
    // For very dense polynomials, try different multiplication algorithms
    if (alen > 100 && blen > 100) {
        // Try Johnson's multiplication for dense polynomials
        fq_nmod_mpoly_mul_johnson(c, a, b, ctx);
    } else {
        // Default multiplication
        fq_nmod_mpoly_mul(c, a, b, ctx);
    }
    fq_nmod_mpoly_reduce_field_equation_inplace(c, ctx);
}

// ============= Conversion Functions for Polynomial Recursive Implementation =============

static slong kronecker_max_univariate_degree_poly_array(fq_nmod_poly_t **poly_matrix,
                                                        slong rows,
                                                        slong cols,
                                                        const fq_nmod_ctx_t ctx) {
    slong max_degree = -1;
    for (slong i = 0; i < rows; i++) {
        for (slong j = 0; j < cols; j++) {
            slong deg = fq_nmod_poly_degree(poly_matrix[i][j], ctx);
            if (deg > max_degree) {
                max_degree = deg;
            }
        }
    }
    return max_degree;
}

static slong kronecker_max_univariate_degree_poly_mat(const fq_nmod_poly_mat_t mat,
                                                      const fq_nmod_ctx_t ctx) {
    slong max_degree = -1;
    slong rows = mat->r;
    slong cols = mat->c;
    for (slong i = 0; i < rows; i++) {
        for (slong j = 0; j < cols; j++) {
            slong deg = fq_nmod_poly_degree(fq_nmod_poly_mat_entry(mat, i, j), ctx);
            if (deg > max_degree) {
                max_degree = deg;
            }
        }
    }
    return max_degree;
}

static int compare_long_desc_local(const void *a, const void *b) {
    long av = *(const long *) a;
    long bv = *(const long *) b;
    if (av < bv) return 1;
    if (av > bv) return -1;
    return 0;
}

static inline ulong addmod_ui_local(ulong a, ulong b, ulong mod) {
    ulong sum = a + b;
    if (sum >= mod || sum < a) {
        sum -= mod;
    }
    return sum;
}

static void kronecker_print_slong_array(const char *label,
                                        const slong *values,
                                        slong length) {
    DET_PRINT("[Kronecker] %s:", label);
    for (slong i = 0; i < length; i++) {
        DET_PRINT(" [%ld]=%ld", i, values[i]);
    }
    DET_PRINT("\n");
}

// Convert fq_mvpoly to fq_nmod_poly for a specific variable
void mvpoly_to_fq_nmod_poly(fq_nmod_poly_t poly, const fq_mvpoly_t *mvpoly, 
                           slong var_index, const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(poly, ctx);
    
    slong nvars = mvpoly->nvars;
    
    for (slong t = 0; t < mvpoly->nterms; t++) {
        slong degree = 0;
        int other_vars_zero = 1;
        
        // Check if this term has non-zero exponents in other variables
        if (mvpoly->terms[t].var_exp) {
            for (slong v = 0; v < nvars; v++) {
                if (v == var_index && var_index < nvars) {
                    degree = mvpoly->terms[t].var_exp[v];
                } else if (mvpoly->terms[t].var_exp[v] > 0) {
                    other_vars_zero = 0;
                    break;
                }
            }
        }
        
        if (mvpoly->terms[t].par_exp && other_vars_zero) {
            for (slong p = 0; p < mvpoly->npars; p++) {
                slong idx = nvars + p;
                if (idx == var_index) {
                    degree = mvpoly->terms[t].par_exp[p];
                } else if (mvpoly->terms[t].par_exp[p] > 0) {
                    other_vars_zero = 0;
                    break;
                }
            }
        }
        
        // Only include terms where all other variables have zero exponent
        if (other_vars_zero) {
            fq_nmod_t existing;
            fq_nmod_init(existing, ctx);
            fq_nmod_poly_get_coeff(existing, poly, degree, ctx);
            fq_nmod_add(existing, existing, mvpoly->terms[t].coeff, ctx);
            fq_nmod_poly_set_coeff(poly, degree, existing, ctx);
            fq_nmod_clear(existing, ctx);
        }
    }
}

// Convert fq_nmod_poly back to fq_mvpoly
void fq_nmod_poly_to_mvpoly(fq_mvpoly_t *mvpoly, const fq_nmod_poly_t poly,
                           slong var_index, slong nvars, slong npars,
                           const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(mvpoly, nvars, npars, ctx);
    
    slong degree = fq_nmod_poly_degree(poly, ctx);
    if (degree < 0) return;
    
    for (slong d = 0; d <= degree; d++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, poly, d, ctx);
        
        if (!fq_nmod_is_zero(coeff, ctx)) {
            slong *var_exp = NULL;
            slong *par_exp = NULL;
            
            if (var_index < nvars && nvars > 0) {
                var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                var_exp[var_index] = d;
            } else if (var_index >= nvars && npars > 0) {
                par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                par_exp[var_index - nvars] = d;
            }
            
            fq_mvpoly_add_term(mvpoly, var_exp, par_exp, coeff);
            
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
        
        fq_nmod_clear(coeff, ctx);
    }
}

// ============= Polynomial Recursive Determinant Implementation =============

// Recursive determinant computation using fq_nmod_poly operations
void compute_det_poly_recursive_helper(fq_nmod_poly_t det, 
                                      fq_nmod_poly_t **matrix,
                                      slong size, 
                                      const fq_nmod_ctx_t ctx) {
    if (size == 0) {
        fq_nmod_poly_one(det, ctx);
        return;
    }
    
    if (size == 1) {
        fq_nmod_poly_set(det, matrix[0][0], ctx);
        return;
    }
    
    if (size == 2) {
        fq_nmod_poly_t ad, bc;
        fq_nmod_poly_init(ad, ctx);
        fq_nmod_poly_init(bc, ctx);
        
        fq_nmod_poly_mul(ad, matrix[0][0], matrix[1][1], ctx);
        fq_nmod_poly_reduce_field_equation_inplace(ad, ctx);
        fq_nmod_poly_mul(bc, matrix[0][1], matrix[1][0], ctx);
        fq_nmod_poly_reduce_field_equation_inplace(bc, ctx);
        fq_nmod_poly_sub(det, ad, bc, ctx);
        
        fq_nmod_poly_clear(ad, ctx);
        fq_nmod_poly_clear(bc, ctx);
        return;
    }
    
    // General case: Laplace expansion
    fq_nmod_poly_zero(det, ctx);
    
    fq_nmod_poly_t cofactor, subdet;
    fq_nmod_poly_init(cofactor, ctx);
    fq_nmod_poly_init(subdet, ctx);
    
    // Allocate submatrix
    fq_nmod_poly_t **submatrix = (fq_nmod_poly_t**) malloc((size-1) * sizeof(fq_nmod_poly_t*));
    for (slong i = 0; i < size-1; i++) {
        submatrix[i] = (fq_nmod_poly_t*) malloc((size-1) * sizeof(fq_nmod_poly_t));
        for (slong j = 0; j < size-1; j++) {
            fq_nmod_poly_init(submatrix[i][j], ctx);
        }
    }
    
    for (slong col = 0; col < size; col++) {
        if (fq_nmod_poly_is_zero(matrix[0][col], ctx)) {
            continue;
        }
        
        // Build submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    fq_nmod_poly_set(submatrix[i-1][sub_j], matrix[i][j], ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive call
        compute_det_poly_recursive_helper(subdet, submatrix, size-1, ctx);
        
        // Multiply by matrix element
        fq_nmod_poly_mul(cofactor, matrix[0][col], subdet, ctx);
        fq_nmod_poly_reduce_field_equation_inplace(cofactor, ctx);
        
        // Add or subtract based on sign
        if (col % 2 == 0) {
            fq_nmod_poly_add(det, det, cofactor, ctx);
        } else {
            fq_nmod_poly_sub(det, det, cofactor, ctx);
        }
    }
    
    // Cleanup
    for (slong i = 0; i < size-1; i++) {
        for (slong j = 0; j < size-1; j++) {
            fq_nmod_poly_clear(submatrix[i][j], ctx);
        }
        free(submatrix[i]);
    }
    free(submatrix);
    
    fq_nmod_poly_clear(cofactor, ctx);
    fq_nmod_poly_clear(subdet, ctx);
}

static slong kronecker_max_univariate_degree_nmod_poly_array(nmod_poly_t **poly_matrix,
                                                             slong rows,
                                                             slong cols) {
    slong max_degree = -1;
    for (slong i = 0; i < rows; i++) {
        for (slong j = 0; j < cols; j++) {
            slong deg = nmod_poly_degree(poly_matrix[i][j]);
            if (deg > max_degree) {
                max_degree = deg;
            }
        }
    }
    return max_degree;
}

static void mvpoly_to_univariate_kronecker_nmod(nmod_poly_t uni_poly,
                                                const fq_mvpoly_t *mv_poly,
                                                const slong *substitution_powers,
                                                ulong prime) {
    nmod_poly_zero(uni_poly);

    if (mv_poly->nterms == 0) return;

    for (slong t = 0; t < mv_poly->nterms; t++) {
        slong uni_exp = 0;
        ulong coeff_val;
        ulong existing;

        if (mv_poly->terms[t].var_exp) {
            for (slong v = 0; v < mv_poly->nvars; v++) {
                uni_exp += mv_poly->terms[t].var_exp[v] * substitution_powers[v];
            }
        }

        if (mv_poly->terms[t].par_exp) {
            for (slong p = 0; p < mv_poly->npars; p++) {
                uni_exp += mv_poly->terms[t].par_exp[p] *
                          substitution_powers[mv_poly->nvars + p];
            }
        }

        coeff_val = nmod_poly_get_coeff_ui(mv_poly->terms[t].coeff, 0);
        if (coeff_val == 0) {
            continue;
        }

        existing = nmod_poly_get_coeff_ui(uni_poly, uni_exp);
        nmod_poly_set_coeff_ui(uni_poly, uni_exp,
                               addmod_ui_local(existing, coeff_val, prime));
    }
}

static void nmod_poly_to_fq_nmod_poly_prime_field(fq_nmod_poly_t dst,
                                                  const nmod_poly_t src,
                                                  const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(dst, ctx);

    for (slong i = 0; i < nmod_poly_length(src); i++) {
        ulong coeff_val = nmod_poly_get_coeff_ui(src, i);
        if (coeff_val != 0) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_set_ui(coeff, coeff_val, ctx);
            fq_nmod_poly_set_coeff(dst, i, coeff, ctx);
            fq_nmod_clear(coeff, ctx);
        }
    }
}

static void compute_det_nmod_poly_recursive_helper(nmod_poly_t det,
                                                   nmod_poly_t **matrix,
                                                   slong size) {
    if (size == 0) {
        nmod_poly_one(det);
        return;
    }

    if (size == 1) {
        nmod_poly_set(det, matrix[0][0]);
        return;
    }

    if (size == 2) {
        nmod_poly_t ad, bc;
        nmod_poly_init(ad, det->mod.n);
        nmod_poly_init(bc, det->mod.n);

        nmod_poly_mul(ad, matrix[0][0], matrix[1][1]);
        nmod_poly_mul(bc, matrix[0][1], matrix[1][0]);
        nmod_poly_sub(det, ad, bc);

        nmod_poly_clear(ad);
        nmod_poly_clear(bc);
        return;
    }

    nmod_poly_zero(det);

    nmod_poly_t cofactor, subdet;
    nmod_poly_init(cofactor, det->mod.n);
    nmod_poly_init(subdet, det->mod.n);

    nmod_poly_t **submatrix = (nmod_poly_t**) malloc((size_t) (size - 1) * sizeof(nmod_poly_t*));
    for (slong i = 0; i < size - 1; i++) {
        submatrix[i] = (nmod_poly_t*) malloc((size_t) (size - 1) * sizeof(nmod_poly_t));
        for (slong j = 0; j < size - 1; j++) {
            nmod_poly_init(submatrix[i][j], det->mod.n);
        }
    }

    for (slong col = 0; col < size; col++) {
        if (nmod_poly_is_zero(matrix[0][col])) {
            continue;
        }

        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    nmod_poly_set(submatrix[i - 1][sub_j], matrix[i][j]);
                    sub_j++;
                }
            }
        }

        compute_det_nmod_poly_recursive_helper(subdet, submatrix, size - 1);
        nmod_poly_mul(cofactor, matrix[0][col], subdet);

        if ((col & 1) == 0) {
            nmod_poly_add(det, det, cofactor);
        } else {
            nmod_poly_sub(det, det, cofactor);
        }
    }

    for (slong i = 0; i < size - 1; i++) {
        for (slong j = 0; j < size - 1; j++) {
            nmod_poly_clear(submatrix[i][j]);
        }
        free(submatrix[i]);
    }
    free(submatrix);

    nmod_poly_clear(cofactor);
    nmod_poly_clear(subdet);
}

// Main function for polynomial recursive algorithm
void compute_fq_det_poly_recursive(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    
    DET_PRINT("Computing %ldx%ld determinant via polynomial recursive with Kronecker\n", size, size);
    DET_PRINT("Variables: %ld, Parameters: %ld\n", nvars, npars);
    
    // Special case: if already univariate, no Kronecker needed
    if (total_vars == 1) {
        DET_PRINT("Already univariate, using direct polynomial recursive\n");
        
        // Convert to fq_nmod_poly format
        fq_nmod_poly_t **poly_matrix = (fq_nmod_poly_t**) malloc(size * sizeof(fq_nmod_poly_t*));
        for (slong i = 0; i < size; i++) {
            poly_matrix[i] = (fq_nmod_poly_t*) malloc(size * sizeof(fq_nmod_poly_t));
            for (slong j = 0; j < size; j++) {
                fq_nmod_poly_init(poly_matrix[i][j], ctx);
                
                // Convert mvpoly to poly (direct conversion for univariate)
                for (slong t = 0; t < matrix[i][j].nterms; t++) {
                    slong degree = 0;
                    if (matrix[i][j].terms[t].var_exp && nvars > 0) {
                        degree = matrix[i][j].terms[t].var_exp[0];
                    } else if (matrix[i][j].terms[t].par_exp && npars > 0) {
                        degree = matrix[i][j].terms[t].par_exp[0];
                    }
                    fq_nmod_poly_set_coeff(poly_matrix[i][j], degree, 
                                          matrix[i][j].terms[t].coeff, ctx);
                }
            }
        }
        
        // Compute determinant
        fq_nmod_poly_t det_poly;
        fq_nmod_poly_init(det_poly, ctx);
        compute_det_poly_recursive_helper(det_poly, poly_matrix, size, ctx);
        
        // Convert back to mvpoly
        fq_mvpoly_init(result, nvars, npars, ctx);
        slong degree = fq_nmod_poly_degree(det_poly, ctx);
        for (slong d = 0; d <= degree; d++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_poly_get_coeff(coeff, det_poly, d, ctx);
            
            if (!fq_nmod_is_zero(coeff, ctx)) {
                if (nvars > 0) {
                    slong *var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                    var_exp[0] = d;
                    fq_mvpoly_add_term(result, var_exp, NULL, coeff);
                    flint_free(var_exp);
                } else {
                    slong *par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                    par_exp[0] = d;
                    fq_mvpoly_add_term(result, NULL, par_exp, coeff);
                    flint_free(par_exp);
                }
            }
            
            fq_nmod_clear(coeff, ctx);
        }
        
        // Cleanup
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_nmod_poly_clear(poly_matrix[i][j], ctx);
            }
            free(poly_matrix[i]);
        }
        free(poly_matrix);
        fq_nmod_poly_clear(det_poly, ctx);
        
        timing_info_t total_elapsed = end_timing(total_start);
        print_timing("Total poly recursive (univariate)", total_elapsed);
        return;
    }
    
    // Multivariate case: use Kronecker+HNF
    
    // Step 1: Compute variable bounds (same as in Kronecker algorithm)
    timing_info_t bounds_start = start_timing();
    slong *var_bounds = (slong*) malloc(total_vars * sizeof(slong));
    compute_kronecker_bounds(var_bounds, matrix, size, nvars, npars);
    timing_info_t bounds_elapsed = end_timing(bounds_start);
    kronecker_print_slong_array("var_bounds", var_bounds, total_vars);
    
    // Step 2: Compute substitution powers
    slong *substitution_powers = (slong*) malloc(total_vars * sizeof(slong));
    substitution_powers[0] = 1;
    for (slong v = 1; v < total_vars; v++) {
        substitution_powers[v] = substitution_powers[v-1] * var_bounds[v-1];
    }
    kronecker_print_slong_array("substitution_powers", substitution_powers, total_vars);
    
    DET_PRINT("Substitution powers: ");
    for (slong v = 0; v < total_vars; v++) {
        DET_PRINT("%ld ", substitution_powers[v]);
    }
    DET_PRINT("\n");
    
    // Step 3: Convert matrix to univariate using Kronecker
    timing_info_t convert_start = start_timing();
    fq_nmod_poly_t **poly_matrix = (fq_nmod_poly_t**) malloc(size * sizeof(fq_nmod_poly_t*));
    for (slong i = 0; i < size; i++) {
        poly_matrix[i] = (fq_nmod_poly_t*) malloc(size * sizeof(fq_nmod_poly_t));
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_init(poly_matrix[i][j], ctx);
            mvpoly_to_univariate_kronecker(poly_matrix[i][j], &matrix[i][j], 
                                          substitution_powers, ctx);
        }
    }
    timing_info_t convert_elapsed = end_timing(convert_start);
    DET_PRINT("[Kronecker] Maximum univariate degree after substitution: %ld\n",
           kronecker_max_univariate_degree_poly_array(poly_matrix, size, size, ctx));
    
    // Step 4: Compute determinant using recursive algorithm
    timing_info_t det_start = start_timing();
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    
    compute_det_poly_recursive_helper(det_poly, poly_matrix, size, ctx);
    
    timing_info_t det_elapsed = end_timing(det_start);
    DET_PRINT("Univariate determinant degree: %ld\n", fq_nmod_poly_degree(det_poly, ctx));
    
    // Step 5: Convert back to multivariate
    timing_info_t back_start = start_timing();
    univariate_to_mvpoly_kronecker(result, det_poly, substitution_powers, 
                                  var_bounds, nvars, npars, ctx);
    timing_info_t back_elapsed = end_timing(back_start);
    
    // Cleanup
    free(var_bounds);
    free(substitution_powers);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_clear(poly_matrix[i][j], ctx);
        }
        free(poly_matrix[i]);
    }
    free(poly_matrix);
    fq_nmod_poly_clear(det_poly, ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    
    if (g_dixon_verbose_level >= 3) {
        printf("\n=== Polynomial Recursive Time Statistics ===\n");
        print_timing("Compute bounds", bounds_elapsed);
        print_timing("Convert to univariate", convert_elapsed);
        print_timing("Recursive determinant", det_elapsed);
        print_timing("Convert back", back_elapsed);
        print_timing("Total poly recursive", total_elapsed);
        printf("Final result: %ld terms\n", result->nterms);
        printf("============================================\n");
    }
}

// ============= Kronecker+HNF Implementation =============

// Compute the Kronecker bound for a multivariate polynomial matrix
void compute_kronecker_bounds(slong *var_bounds, fq_mvpoly_t **matrix, 
                             slong size, slong nvars, slong npars) {
    slong total_vars = nvars + npars;
    slong *entry_max = NULL;
    int is_step1_structured = (size > 0 && nvars > 0 && (nvars % 2) == 0 &&
                               size == nvars / 2 + 1);
    
    // Initialize bounds
    for (slong v = 0; v < total_vars; v++) {
        var_bounds[v] = 0;
    }

    if (total_vars > 0) {
        entry_max = (slong *) flint_calloc(total_vars, sizeof(slong));
    }
    
    // Find maximum degree in each variable across all matrix entries
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *poly = &matrix[i][j];
            
            for (slong t = 0; t < poly->nterms; t++) {
                // Check variable degrees
                if (poly->terms[t].var_exp && nvars > 0) {
                    for (slong v = 0; v < nvars && v < poly->nvars; v++) {
                        if (poly->terms[t].var_exp[v] > entry_max[v]) {
                            entry_max[v] = poly->terms[t].var_exp[v];
                        }
                    }
                }
                
                // Check parameter degrees
                if (poly->terms[t].par_exp && npars > 0) {
                    for (slong p = 0; p < npars && p < poly->npars; p++) {
                        if (poly->terms[t].par_exp[p] > entry_max[nvars + p]) {
                            entry_max[nvars + p] = poly->terms[t].par_exp[p];
                        }
                    }
                }
            }
        }
    }

    if (is_step1_structured) {
        slong elim_vars = nvars / 2;
        long *col_degrees = (long *) flint_calloc(size, sizeof(long));
        long *sorted = (long *) flint_calloc(size, sizeof(long));
        long *prefix = (long *) flint_calloc(size + 1, sizeof(long));

        for (slong col = 0; col < size; col++) {
            long col_max_total = 0;
            for (slong row = 0; row < size; row++) {
                fq_mvpoly_t *poly = &matrix[row][col];
                for (slong t = 0; t < poly->nterms; t++) {
                    long total_deg = 0;
                    if (poly->terms[t].var_exp) {
                        for (slong v = 0; v < poly->nvars; v++) {
                            total_deg += poly->terms[t].var_exp[v];
                        }
                    }
                    if (poly->terms[t].par_exp) {
                        for (slong p = 0; p < poly->npars; p++) {
                            total_deg += poly->terms[t].par_exp[p];
                        }
                    }
                    if (total_deg > col_max_total) {
                        col_max_total = total_deg;
                    }
                }
            }
            col_degrees[col] = col_max_total;
            sorted[col] = col_max_total;
        }

        qsort(sorted, (size_t) size, sizeof(long), compare_long_desc_local);
        for (slong i = 0; i < size; i++) {
            prefix[i + 1] = prefix[i] + sorted[i];
        }

        for (slong i = 0; i < elim_vars; i++) {
            slong orig_count = i + 1;
            slong dual_count = elim_vars - i;
            slong orig_det_bound = prefix[orig_count] - orig_count;
            slong dual_det_bound = prefix[dual_count] - dual_count;
            if (orig_det_bound < 0) orig_det_bound = 0;
            if (dual_det_bound < 0) dual_det_bound = 0;

            var_bounds[i] = FLINT_MAX(entry_max[i], orig_det_bound) + 1;
            var_bounds[elim_vars + i] = FLINT_MAX(entry_max[elim_vars + i], dual_det_bound) + 1;
        }

        for (slong p = 0; p < npars; p++) {
            slong param_det_bound = prefix[size] - elim_vars;
            if (param_det_bound < 0) param_det_bound = 0;
            var_bounds[nvars + p] = FLINT_MAX(entry_max[nvars + p], param_det_bound) + 1;
        }

        flint_free(col_degrees);
        flint_free(sorted);
        flint_free(prefix);
        if (entry_max) flint_free(entry_max);
        return;
    }
    
    // Compute degree bound for determinant (sum of row maximums)
    for (slong v = 0; v < total_vars; v++) {
        slong det_bound = 0;
        
        for (slong row = 0; row < size; row++) {
            slong row_max = 0;
            
            for (slong col = 0; col < size; col++) {
                fq_mvpoly_t *poly = &matrix[row][col];
                
                for (slong t = 0; t < poly->nterms; t++) {
                    slong deg = 0;
                    
                    if (v < nvars && poly->terms[t].var_exp && v < poly->nvars) {
                        deg = poly->terms[t].var_exp[v];
                    } else if (v >= nvars && poly->terms[t].par_exp && 
                              v - nvars < poly->npars) {
                        deg = poly->terms[t].par_exp[v - nvars];
                    }
                    
                    if (deg > row_max) row_max = deg;
                }
            }
            det_bound += row_max;
        }
        
        var_bounds[v] = FLINT_MAX(entry_max ? entry_max[v] : 0, det_bound) + 1;  // Add 1 for safety
    }

    if (entry_max) flint_free(entry_max);
}

// Convert multivariate polynomial to univariate using Kronecker+HNF
void mvpoly_to_univariate_kronecker(fq_nmod_poly_t uni_poly,
                                   const fq_mvpoly_t *mv_poly,
                                   const slong *substitution_powers,
                                   const fq_nmod_ctx_t ctx) {
    fq_nmod_poly_zero(uni_poly, ctx);
    
    if (mv_poly->nterms == 0) return;
    
    slong total_vars = mv_poly->nvars + mv_poly->npars;
    
    for (slong t = 0; t < mv_poly->nterms; t++) {
        slong uni_exp = 0;
        
        // Compute univariate exponent: sum of var_exp[i] * substitution_powers[i]
        if (mv_poly->terms[t].var_exp) {
            for (slong v = 0; v < mv_poly->nvars; v++) {
                uni_exp += mv_poly->terms[t].var_exp[v] * substitution_powers[v];
            }
        }
        
        if (mv_poly->terms[t].par_exp) {
            for (slong p = 0; p < mv_poly->npars; p++) {
                uni_exp += mv_poly->terms[t].par_exp[p] * 
                          substitution_powers[mv_poly->nvars + p];
            }
        }
        
        // Add coefficient at computed degree
        fq_nmod_t existing;
        fq_nmod_init(existing, ctx);
        fq_nmod_poly_get_coeff(existing, uni_poly, uni_exp, ctx);
        fq_nmod_add(existing, existing, mv_poly->terms[t].coeff, ctx);
        fq_nmod_poly_set_coeff(uni_poly, uni_exp, existing, ctx);
        fq_nmod_clear(existing, ctx);
    }
}

// Convert univariate polynomial back to multivariate
void univariate_to_mvpoly_kronecker(fq_mvpoly_t *mv_poly,
                                   const fq_nmod_poly_t uni_poly,
                                   const slong *substitution_powers,
                                   const slong *var_bounds,
                                   slong nvars, slong npars,
                                   const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(mv_poly, nvars, npars, ctx);
    
    slong degree = fq_nmod_poly_degree(uni_poly, ctx);
    if (degree < 0) return;
    
    slong total_vars = nvars + npars;
    
    for (slong d = 0; d <= degree; d++) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        fq_nmod_poly_get_coeff(coeff, uni_poly, d, ctx);
        
        if (!fq_nmod_is_zero(coeff, ctx)) {
            // Decompose d into multivariate exponents
            slong *var_exp = NULL;
            slong *par_exp = NULL;
            
            if (nvars > 0) {
                var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            }
            if (npars > 0) {
                par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            }
            
            slong remaining = d;
            
            // Extract exponents in reverse order (largest substitution power first)
            for (slong v = total_vars - 1; v >= 0; v--) {
                slong exp = remaining / substitution_powers[v];
                remaining = remaining % substitution_powers[v];
                
                if (v < nvars && var_exp) {
                    var_exp[v] = exp;
                } else if (v >= nvars && par_exp) {
                    par_exp[v - nvars] = exp;
                }
            }
            
            fq_mvpoly_add_term(mv_poly, var_exp, par_exp, coeff);
            
            if (var_exp) flint_free(var_exp);
            if (par_exp) flint_free(par_exp);
        }
        
        fq_nmod_clear(coeff, ctx);
    }
}

// Compute determinant using Kronecker+HNF
void compute_fq_det_kronecker(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    
    DET_PRINT("Computing %ldx%ld determinant via Kronecker+HNF\n", size, size);
    DET_PRINT("Variables: %ld, Parameters: %ld\n", nvars, npars);
    
    // Step 1: Compute variable bounds
    timing_info_t bounds_start = start_timing();
    slong *var_bounds = (slong*) malloc(total_vars * sizeof(slong));
    compute_kronecker_bounds(var_bounds, matrix, size, nvars, npars);
    timing_info_t bounds_elapsed = end_timing(bounds_start);
    kronecker_print_slong_array("var_bounds", var_bounds, total_vars);
    
    DET_PRINT("Variable bounds: ");
    for (slong v = 0; v < total_vars; v++) {
        DET_PRINT("%ld ", var_bounds[v]);
    }
    DET_PRINT("\n");
    
    // Step 2: Compute substitution powers
    slong *substitution_powers = (slong*) malloc(total_vars * sizeof(slong));
    substitution_powers[0] = 1;
    for (slong v = 1; v < total_vars; v++) {
        substitution_powers[v] = substitution_powers[v-1] * var_bounds[v-1];
    }
    kronecker_print_slong_array("substitution_powers", substitution_powers, total_vars);
    
    DET_PRINT("Substitution powers: ");
    for (slong v = 0; v < total_vars; v++) {
        DET_PRINT("%ld ", substitution_powers[v]);
    }
    DET_PRINT("\n");
    
    // Step 3: Convert matrix to univariate
    timing_info_t convert_start = start_timing();
    fq_nmod_poly_mat_t uni_mat;
    fq_nmod_poly_mat_init(uni_mat, size, size, ctx);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            mvpoly_to_univariate_kronecker(fq_nmod_poly_mat_entry(uni_mat, i, j),
                                          &matrix[i][j], substitution_powers, ctx);
        }
    }
    timing_info_t convert_elapsed = end_timing(convert_start);
    DET_PRINT("[Kronecker] Maximum univariate degree after substitution: %ld\n",
           kronecker_max_univariate_degree_poly_mat(uni_mat, ctx));
    
    // Step 4: Compute univariate determinant
    timing_info_t det_start = start_timing();
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    
    // Use the optimized univariate determinant function
    fq_nmod_poly_mat_det_iter(det_poly, uni_mat, ctx);
    
    timing_info_t det_elapsed = end_timing(det_start);
    DET_PRINT("Univariate determinant degree: %ld\n", fq_nmod_poly_degree(det_poly, ctx));
    
    // Step 5: Convert back to multivariate
    timing_info_t back_start = start_timing();
    univariate_to_mvpoly_kronecker(result, det_poly, substitution_powers, 
                                  var_bounds, nvars, npars, ctx);
    timing_info_t back_elapsed = end_timing(back_start);
    
    // Cleanup
    free(var_bounds);
    free(substitution_powers);
    fq_nmod_poly_clear(det_poly, ctx);
    fq_nmod_poly_mat_clear(uni_mat, ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    
    
    if (g_dixon_verbose_level >= 3) {
        printf("\n=== Kronecker+HNF Time Statistics ===\n");
        print_timing("Compute bounds", bounds_elapsed);
        print_timing("Convert to univariate", convert_elapsed);
        print_timing("Univariate determinant", det_elapsed);
        print_timing("Convert back", back_elapsed);
        print_timing("Total Kronecker", total_elapsed);
        printf("Final result: %ld terms\n", result->nterms);
        printf("==============================================\n");
    }
    
}

// ============= Prime Field Conversion Functions Implementation =============

void fq_mvpoly_to_nmod_mpoly(nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                             nmod_mpoly_ctx_t mpoly_ctx) {
    nmod_mpoly_zero(mpoly, mpoly_ctx);
    
    if (poly->nterms == 0) return;
    
    slong total_vars = poly->nvars + poly->npars;
    
    // Pre-allocate space for better performance
    nmod_mpoly_fit_length(mpoly, poly->nterms, mpoly_ctx);
    
    for (slong i = 0; i < poly->nterms; i++) {
        ulong *exps = (ulong*) flint_calloc(total_vars, sizeof(ulong));
        
        if (poly->terms[i].var_exp && poly->nvars > 0) {
            for (slong j = 0; j < poly->nvars; j++) {
                exps[j] = (ulong)poly->terms[i].var_exp[j];
            }
        }
        
        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                exps[poly->nvars + j] = (ulong)poly->terms[i].par_exp[j];
            }
        }
        
        // For prime fields, extract the coefficient as ulong
        ulong coeff_val = nmod_poly_get_coeff_ui(poly->terms[i].coeff, 0);
        nmod_mpoly_push_term_ui_ui(mpoly, coeff_val, exps, mpoly_ctx);
        flint_free(exps);
    }
    
    nmod_mpoly_sort_terms(mpoly, mpoly_ctx);
    nmod_mpoly_combine_like_terms(mpoly, mpoly_ctx);
}

void nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *result, const nmod_mpoly_t poly,
                             slong nvars, slong npars,
                             const nmod_mpoly_ctx_t mpoly_ctx,
                             const fq_nmod_ctx_t field_ctx) {
    fq_mvpoly_init(result, nvars, npars, field_ctx);
    
    slong nterms = nmod_mpoly_length(poly, mpoly_ctx);
    if (nterms == 0) return;
    
    slong total_vars = nmod_mpoly_ctx_nvars(mpoly_ctx);
    
    if (result->alloc < nterms) {
        result->alloc = nterms;
        result->terms = (fq_monomial_t*) flint_realloc(result->terms, 
                                                        result->alloc * sizeof(fq_monomial_t));
    }
    
    // Batch allocate exponent buffer
    ulong *exp_buffer = (ulong*) flint_malloc(total_vars * sizeof(ulong));
    
    // Batch convert with proper initialization
    for (slong i = 0; i < nterms; i++) {
        // Get coefficient
        mp_limb_t coeff_limb = nmod_mpoly_get_term_coeff_ui(poly, i, mpoly_ctx);
        
        // Initialize coefficient (important!)
        fq_nmod_init(result->terms[i].coeff, field_ctx);
        fq_nmod_set_ui(result->terms[i].coeff, coeff_limb, field_ctx);
        
        // Get exponents
        nmod_mpoly_get_term_exp_ui(exp_buffer, poly, i, mpoly_ctx);
        
        // Allocate and set variable exponents
        if (nvars > 0) {
            result->terms[i].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars && j < total_vars; j++) {
                result->terms[i].var_exp[j] = (slong)exp_buffer[j];
            }
        } else {
            result->terms[i].var_exp = NULL;
        }
        
        // Allocate and set parameter exponents
        if (npars > 0 && total_vars > nvars) {
            result->terms[i].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars && (nvars + j) < total_vars; j++) {
                result->terms[i].par_exp[j] = (slong)exp_buffer[nvars + j];
            }
        } else {
            result->terms[i].par_exp = NULL;
        }
    }
    
    // Set term count
    result->nterms = nterms;
    
    // Cleanup
    flint_free(exp_buffer);
}

void fq_matrix_mvpoly_to_nmod_mpoly(nmod_mpoly_t **mpoly_matrix, 
                                   fq_mvpoly_t **mvpoly_matrix, 
                                   slong size, 
                                   nmod_mpoly_ctx_t mpoly_ctx) {
    DET_PRINT("Converting %ld x %ld matrix to nmod_mpoly\n", size, size);
    
    timing_info_t start = start_timing();
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            nmod_mpoly_init(mpoly_matrix[i][j], mpoly_ctx);
            fq_mvpoly_to_nmod_mpoly(mpoly_matrix[i][j], &mvpoly_matrix[i][j], mpoly_ctx);
        }
    }
    
    timing_info_t elapsed = end_timing(start);
    print_timing("Matrix conversion to nmod_mpoly", elapsed);
}

// ============= Prime Field Determinant Computation Implementation =============

// Hand-optimized 3x3 determinant for nmod_mpoly
void compute_det_3x3_nmod_optimized(nmod_mpoly_t det, 
                                   nmod_mpoly_t **m,
                                   nmod_mpoly_ctx_t ctx) {
    nmod_mpoly_t t1, t2, t3, t4, t5, t6, sum;
    
    // Initialize temporaries
    nmod_mpoly_init(t1, ctx);
    nmod_mpoly_init(t2, ctx);
    nmod_mpoly_init(t3, ctx);
    nmod_mpoly_init(t4, ctx);
    nmod_mpoly_init(t5, ctx);
    nmod_mpoly_init(t6, ctx);
    nmod_mpoly_init(sum, ctx);
    
    // Compute 6 products in parallel if beneficial
    #pragma omp parallel sections if(omp_get_max_threads() > 2)
    {
        #pragma omp section
        {
            nmod_mpoly_mul(t1, m[1][1], m[2][2], ctx);
            nmod_mpoly_mul(t1, m[0][0], t1, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t2, m[1][2], m[2][0], ctx);
            nmod_mpoly_mul(t2, m[0][1], t2, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t3, m[1][0], m[2][1], ctx);
            nmod_mpoly_mul(t3, m[0][2], t3, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t4, m[1][0], m[2][2], ctx);
            nmod_mpoly_mul(t4, m[0][1], t4, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t5, m[1][1], m[2][0], ctx);
            nmod_mpoly_mul(t5, m[0][2], t5, ctx);
        }
        #pragma omp section
        {
            nmod_mpoly_mul(t6, m[1][2], m[2][1], ctx);
            nmod_mpoly_mul(t6, m[0][0], t6, ctx);
        }
    }
    
    // Sum with signs
    nmod_mpoly_add(sum, t1, t2, ctx);
    nmod_mpoly_add(sum, sum, t3, ctx);
    nmod_mpoly_sub(sum, sum, t4, ctx);
    nmod_mpoly_sub(sum, sum, t5, ctx);
    nmod_mpoly_sub(det, sum, t6, ctx);
    
    // Cleanup
    nmod_mpoly_clear(t1, ctx);
    nmod_mpoly_clear(t2, ctx);
    nmod_mpoly_clear(t3, ctx);
    nmod_mpoly_clear(t4, ctx);
    nmod_mpoly_clear(t5, ctx);
    nmod_mpoly_clear(t6, ctx);
    nmod_mpoly_clear(sum, ctx);
}

// Recursive determinant for nmod_mpoly
void compute_nmod_mpoly_det_recursive(nmod_mpoly_t det_result, 
                                     nmod_mpoly_t **mpoly_matrix, 
                                     slong size, 
                                     nmod_mpoly_ctx_t mpoly_ctx) {
    if (size <= 0) {
        nmod_mpoly_one(det_result, mpoly_ctx);
        return;
    }
    
    if (size == 1) {
        nmod_mpoly_set(det_result, mpoly_matrix[0][0], mpoly_ctx);
        return;
    }
    
    if (size == 2) {
        nmod_mpoly_t ad, bc;
        nmod_mpoly_init(ad, mpoly_ctx);
        nmod_mpoly_init(bc, mpoly_ctx);
        
        nmod_mpoly_mul(ad, mpoly_matrix[0][0], mpoly_matrix[1][1], mpoly_ctx);
        nmod_mpoly_mul(bc, mpoly_matrix[0][1], mpoly_matrix[1][0], mpoly_ctx);
        nmod_mpoly_sub(det_result, ad, bc, mpoly_ctx);
        
        nmod_mpoly_clear(ad, mpoly_ctx);
        nmod_mpoly_clear(bc, mpoly_ctx);
        return;
    }
    
    if (size == 3) {
        compute_det_3x3_nmod_optimized(det_result, mpoly_matrix, mpoly_ctx);
        return;
    }
    
    // General case: Laplace expansion
    nmod_mpoly_zero(det_result, mpoly_ctx);
    
    nmod_mpoly_t temp_result, cofactor, subdet;
    nmod_mpoly_init(temp_result, mpoly_ctx);
    nmod_mpoly_init(cofactor, mpoly_ctx);
    nmod_mpoly_init(subdet, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            continue;
        }
        
        // Create submatrix
        nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
            for (slong j = 0; j < size-1; j++) {
                nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
            }
        }
        
        // Fill submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive computation
        compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
        
        // Compute cofactor
        nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
        
        // Add/subtract to result
        if (col % 2 == 0) {
            nmod_mpoly_add(temp_result, det_result, cofactor, mpoly_ctx);
        } else {
            nmod_mpoly_sub(temp_result, det_result, cofactor, mpoly_ctx);
        }
        nmod_mpoly_set(det_result, temp_result, mpoly_ctx);
        
        // Cleanup submatrix
        for (slong i = 0; i < size-1; i++) {
            for (slong j = 0; j < size-1; j++) {
                nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
    
    nmod_mpoly_clear(temp_result, mpoly_ctx);
    nmod_mpoly_clear(cofactor, mpoly_ctx);
    nmod_mpoly_clear(subdet, mpoly_ctx);
}

// Parallel determinant computation for nmod_mpoly with proper nested parallelism
void compute_nmod_mpoly_det_parallel_optimized(nmod_mpoly_t det_result, 
                                              nmod_mpoly_t **mpoly_matrix, 
                                              slong size, 
                                              nmod_mpoly_ctx_t mpoly_ctx,
                                              slong depth) {
    // For deep recursion or small matrices, use sequential
    if (size < PARALLEL_THRESHOLD || depth >= MAX_PARALLEL_DEPTH) {
        compute_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    DET_PRINT("Parallel nmod computation for %ld x %ld matrix (depth %ld)\n", size, size, depth);
    
    if (size <= 3) {
        compute_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    nmod_mpoly_zero(det_result, mpoly_ctx);
    
    // Count non-zero entries in first row
    slong nonzero_count = 0;
    for (slong col = 0; col < size; col++) {
        if (!nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count < 2) {
        compute_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    // Allocate space for partial results
    nmod_mpoly_t *partial_results = (nmod_mpoly_t*) flint_malloc(size * sizeof(nmod_mpoly_t));
    for (slong i = 0; i < size; i++) {
        nmod_mpoly_init(partial_results[i], mpoly_ctx);
        nmod_mpoly_zero(partial_results[i], mpoly_ctx);
    }
    
    // Determine parallelism strategy based on depth
    if (depth == 0) {
        // First level: use parallel for with nested parallelism enabled
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, omp_get_max_threads()))
        for (slong col = 0; col < size; col++) {
            if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            nmod_mpoly_t cofactor, subdet;
            nmod_mpoly_init(cofactor, mpoly_ctx);
            nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create submatrix
            nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            // Fill submatrix
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Recursive call - this will use nested parallelism at depth 1
            compute_nmod_mpoly_det_parallel_optimized(subdet, submatrix, size-1, mpoly_ctx, depth+1);
            
            // Compute cofactor
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
            // Store with sign
            if (col % 2 == 0) {
                nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            nmod_mpoly_clear(cofactor, mpoly_ctx);
            nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    } else if (depth == 1 && size >= PARALLEL_THRESHOLD) {
        // Second level: also use parallel for, but with fewer threads
        slong max_threads_level2 = FLINT_MAX(1, omp_get_max_threads() / size);
        
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, max_threads_level2)) if(nonzero_count >= 3)
        for (slong col = 0; col < size; col++) {
            if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            nmod_mpoly_t cofactor, subdet;
            nmod_mpoly_init(cofactor, mpoly_ctx);
            nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create submatrix
            nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            // Fill submatrix
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // At depth > 1, use sequential computation
            compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
            // Store with sign
            if (col % 2 == 0) {
                nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            nmod_mpoly_clear(cofactor, mpoly_ctx);
            nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    } else {
        // Sequential fallback for deeper levels
        for (slong col = 0; col < size; col++) {
            if (nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            nmod_mpoly_t cofactor, subdet;
            nmod_mpoly_init(cofactor, mpoly_ctx);
            nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create and fill submatrix
            nmod_mpoly_t **submatrix = (nmod_mpoly_t**) flint_malloc((size-1) * sizeof(nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (nmod_mpoly_t*) flint_malloc((size-1) * sizeof(nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Sequential computation
            compute_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor and store
            nmod_mpoly_mul(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            if (col % 2 == 0) {
                nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            nmod_mpoly_clear(cofactor, mpoly_ctx);
            nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    }
    
    // Sum results (sequential to avoid race conditions)
    nmod_mpoly_t temp_sum;
    nmod_mpoly_init(temp_sum, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (!nmod_mpoly_is_zero(partial_results[col], mpoly_ctx)) {
            nmod_mpoly_add(temp_sum, det_result, partial_results[col], mpoly_ctx);
            nmod_mpoly_set(det_result, temp_sum, mpoly_ctx);
        }
        nmod_mpoly_clear(partial_results[col], mpoly_ctx);
    }
    
    nmod_mpoly_clear(temp_sum, mpoly_ctx);
    flint_free(partial_results);
}

/* Separate x/y order ideals: down(A x B) = down(A) x down(B).
 * Packed keys are only used when all coordinates fit in one limb. Zero is
 * represented by key 1, leaving 0 as the empty hash bucket. The tables are
 * immutable during the parallel minor DP. Parameters are never truncated. */
typedef struct {
    ulong *keys;
    slong alloc, count;
} mq_monom_set;

typedef struct {
    slong nvars;
    unsigned bits;
    ulong digit_mask;
    mq_monom_set rows, cols, row_targets, col_targets;
    slong safe_linear_layers;
    flint_bitcnt_t packed_bits;
    mq_monom_set packed[4];
} mq_det_filter;

#define MQ_FILTER_MAX_MONOMS (1L << 20)

static ulong mq_monom_hash(ulong key)
{
    /* Native FLINT keys can have their only nonzero coordinate near bit 56.
     * Mix high bits before masking a small power-of-two table. */
#if FLINT_BITS == 64
    key ^= key >> 33;
    key *= UWORD(0xff51afd7ed558ccd);
    key ^= key >> 33;
    key *= UWORD(0xc4ceb9fe1a85ec53);
    return key ^ (key >> 33);
#else
    key ^= key >> 16;
    key *= UWORD(0x85ebca6b);
    key ^= key >> 13;
    key *= UWORD(0xc2b2ae35);
    return key ^ (key >> 16);
#endif
}

static int mq_monom_contains(const mq_monom_set *set, ulong code)
{
    ulong key = code + 1, pos = mq_monom_hash(key) & (set->alloc - 1);
    while (set->keys[pos]) {
        if (set->keys[pos] == key) return 1;
        pos = (pos + 1) & (set->alloc - 1);
    }
    return 0;
}

/* Return 1 for a new key, 0 for a duplicate, -1 at the memory budget. */
static int mq_monom_insert(mq_monom_set *set, ulong code)
{
    if (mq_monom_contains(set, code)) return 0;
    if (set->count >= MQ_FILTER_MAX_MONOMS) return -1;
    if (2 * set->count >= set->alloc) {
        slong old_alloc = set->alloc;
        ulong *old = set->keys;
        set->alloc *= 2;
        set->keys = flint_calloc((size_t) set->alloc, sizeof(ulong));
        for (slong i = 0; i < old_alloc; i++) if (old[i]) {
            ulong pos = mq_monom_hash(old[i]) & (set->alloc - 1);
            while (set->keys[pos]) pos = (pos + 1) & (set->alloc - 1);
            set->keys[pos] = old[i];
        }
        flint_free(old);
    }
    ulong key = code + 1, pos = mq_monom_hash(key) & (set->alloc - 1);
    while (set->keys[pos]) pos = (pos + 1) & (set->alloc - 1);
    set->keys[pos] = key;
    set->count++;
    return 1;
}

static int mq_monom_close(mq_monom_set *set, ulong code, const mq_det_filter *f)
{
    int inserted = mq_monom_insert(set, code);
    if (inserted <= 0) return inserted == 0;
    for (slong v = 0; v < f->nvars; v++) {
        unsigned shift = v * f->bits;
        if (((code >> shift) & f->digit_mask) &&
            !mq_monom_close(set, code - (UWORD(1) << shift), f)) return 0;
    }
    return 1;
}

static void mq_filter_clear(mq_det_filter *f)
{
    flint_free(f->rows.keys); flint_free(f->cols.keys);
    flint_free(f->row_targets.keys); flint_free(f->col_targets.keys);
    for (int i = 0; i < 4; i++) flint_free(f->packed[i].keys);
}

static int mq_filter_init(mq_det_filter *f, slong nvars,
                         const slong *rows, const slong *cols, slong count)
{
    ulong largest = 0;
    memset(f, 0, sizeof(*f));
    if (nvars <= 0 || count <= 0) return 0;
    for (slong i = 0; i < count * nvars; i++) {
        if (rows[i] < 0 || cols[i] < 0) return 0;
        largest = FLINT_MAX(largest, (ulong) FLINT_MAX(rows[i], cols[i]));
    }
    f->nvars = nvars; f->bits = 1;
    while (largest >>= 1) f->bits++;
    if (nvars > (FLINT_BITS - 1) / f->bits) return 0;
    f->digit_mask = (UWORD(1) << f->bits) - 1;
    mq_monom_set *sets[] = {&f->rows, &f->cols, &f->row_targets, &f->col_targets};
    for (int i = 0; i < 4; i++) {
        sets[i]->alloc = 16;
        sets[i]->keys = flint_calloc(16, sizeof(ulong));
    }
    for (slong i = 0; i < count; i++) {
        ulong r = 0, c = 0;
        for (slong v = 0; v < nvars; v++) {
            r |= (ulong) rows[i * nvars + v] << (v * f->bits);
            c |= (ulong) cols[i * nvars + v] << (v * f->bits);
        }
        if (mq_monom_insert(&f->row_targets, r) < 0 ||
            mq_monom_insert(&f->col_targets, c) < 0 ||
            !mq_monom_close(&f->rows, r, f) || !mq_monom_close(&f->cols, c, f)) {
            mq_filter_clear(f);
            return 0;
        }
    }
    return 1;
}

/* Every linear difference row contributes either zero or one axis variable.
 * Ignore coefficients and column exclusivity, and include zero even when it
 * is absent: the resulting support is a safe upper bound for EVERY minor of
 * these trailing rows. Stop at the first term outside the order ideal. Since
 * zero was included, this bound only grows with the number of rows. */
static slong mq_safe_axis_layers(const mq_det_filter *f, const mq_monom_set *closure,
                                 fq_mvpoly_t **matrix, slong size, slong axis)
{
    mq_monom_set previous = {flint_calloc(16, sizeof(ulong)), 16, 0};
    mq_monom_insert(&previous, 0);
    slong safe = 0;
    for (slong row = size - 1; row >= 1; row--) {
        ulong active = 0;
        for (slong col = 0; col < size; col++) {
            const fq_mvpoly_t *p = &matrix[row][col];
            for (slong t = 0; t < p->nterms; t++) if (p->terms[t].var_exp)
                for (slong v = 0; v < f->nvars; v++)
                    if (p->terms[t].var_exp[axis * f->nvars + v]) active |= UWORD(1) << v;
        }
        mq_monom_set next = {flint_calloc(16, sizeof(ulong)), 16, 0};
        int inside = 1;
        for (slong i = 0; inside && i < previous.alloc; i++) if (previous.keys[i]) {
            ulong code = previous.keys[i] - 1;
            if (mq_monom_insert(&next, code) < 0) { inside = 0; break; }
            for (slong v = 0; v < f->nvars; v++) if (active & (UWORD(1) << v)) {
                unsigned shift = v * f->bits;
                /* Do not let a packed-coordinate carry alias another monomial. */
                if (((code >> shift) & f->digit_mask) == f->digit_mask) { inside = 0; break; }
                ulong added = code + (UWORD(1) << shift);
                if (!mq_monom_contains(closure, added) || mq_monom_insert(&next, added) < 0) {
                    inside = 0; break;
                }
            }
        }
        flint_free(previous.keys);
        previous = next;
        if (!inside) break;
        safe++;
    }
    flint_free(previous.keys);
    return safe;
}

/* For short axes use the SAME packing as FLINT's lex exponents. An axis can
 * then be extracted with at most two word loads, without visiting variables.
 * Build these read-only tables once, before entering any parallel region. */
static void mq_filter_prepare_packed(mq_det_filter *f, const nmod_mpoly_ctx_t ctx,
                                     slong degree_bound)
{
    flint_bitcnt_t bits = mpoly_fix_bits(1 + FLINT_BIT_COUNT((ulong) degree_bound), ctx->minfo);
    if (ctx->minfo->ord != ORD_LEX || bits >= FLINT_BITS ||
        f->nvars > (FLINT_BITS - 1) / bits) return;
    f->packed_bits = bits;
    mq_monom_set *source[] = {&f->rows, &f->cols, &f->row_targets, &f->col_targets};
    for (int s = 0; s < 4; s++) {
        mq_monom_set *dest = &f->packed[s];
        dest->alloc = source[s]->alloc;
        dest->keys = flint_calloc((size_t) dest->alloc, sizeof(ulong));
        for (slong i = 0; i < source[s]->alloc; i++) if (source[s]->keys[i]) {
            ulong code = source[s]->keys[i] - 1, native = 0;
            for (slong v = 0; v < f->nvars; v++)
                native |= ((code >> (v * f->bits)) & f->digit_mask)
                            << ((f->nvars - 1 - v) * bits);
            mq_monom_insert(dest, native);
        }
    }
}

static int mq_filter_accepts(const mq_det_filter *f, const ulong *exp, int target)
{
    ulong r = 0, c = 0;
    for (slong v = 0; v < f->nvars; v++) {
        if (exp[v] > f->digit_mask || exp[f->nvars + v] > f->digit_mask) return 0;
        r |= exp[v] << (v * f->bits);
        c |= exp[f->nvars + v] << (v * f->bits);
    }
    return mq_monom_contains(target ? &f->row_targets : &f->rows, r) &&
           mq_monom_contains(target ? &f->col_targets : &f->cols, c);
}

/* A lex axis may straddle a word boundary with unused padding between
 * words (e.g. seven 9-bit fields per limb). Drop that padding when joining
 * the two pieces. The prepared widths are nonzero and below FLINT_BITS. */
static inline ulong mq_packed_axis(const ulong *exp, slong word, slong shift,
                                   slong low_bits, slong width)
{
    ulong code = (exp[word] >> shift) & (UWORD_MAX >> (FLINT_BITS - low_bits));
    if (width > low_bits) code |= exp[word + 1] << low_bits;
    return code & (UWORD_MAX >> (FLINT_BITS - width));
}

/* Compact canonical FLINT storage in place. Native lex axes need no exponent
 * buffer. Larger axes use cached field locations; multiword fields retain
 * the generic FLINT unpacker. All scratch is bounded stack storage. */
static void mq_filter_poly(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx,
                           const mq_det_filter *f, int target)
{
    slong out = 0, words = mpoly_words_per_exp(poly->bits, ctx->minfo);
    slong offset[2 * FLINT_BITS], shift[2 * FLINT_BITS];
    ulong exp[2 * FLINT_BITS + 1];
    int native = f->packed_bits && poly->bits == f->packed_bits && ctx->minfo->ord == ORD_LEX;
    int small_fields = poly->bits <= FLINT_BITS;
    slong axis_offset[2], axis_shift[2], low_bits[2];
    slong width = f->nvars * poly->bits;
    if (native) {
        for (slong a = 0; a < 2; a++) {
            mpoly_gen_offset_shift_sp(axis_offset + a, axis_shift + a,
                                      (a + 1) * f->nvars - 1, poly->bits, ctx->minfo);
            low_bits[a] = FLINT_MIN(width, ((FLINT_BITS - axis_shift[a]) / poly->bits) * poly->bits);
        }
    }
    ulong field_mask = small_fields ? UWORD_MAX >> (FLINT_BITS - poly->bits) : 0;
    if (!native && small_fields)
        for (slong v = 0; v < 2 * f->nvars; v++)
            mpoly_gen_offset_shift_sp(offset + v, shift + v, v, poly->bits, ctx->minfo);
    const mq_monom_set *rows = target ? &f->row_targets : &f->rows;
    const mq_monom_set *cols = target ? &f->col_targets : &f->cols;
    for (slong i = 0; i < poly->length; i++) {
        const ulong *packed = poly->exps + i * words;
        int keep;
        if (native) {
            ulong r = mq_packed_axis(packed, axis_offset[0], axis_shift[0], low_bits[0], width);
            ulong c = mq_packed_axis(packed, axis_offset[1], axis_shift[1], low_bits[1], width);
            keep = mq_monom_contains(&f->packed[target ? 2 : 0], r) &&
                   mq_monom_contains(&f->packed[target ? 3 : 1], c);
        } else if (small_fields) {
            ulong r = 0, c = 0;
            keep = 1;
            for (slong v = 0; v < f->nvars; v++) {
                ulong e = (packed[offset[v]] >> shift[v]) & field_mask;
                if (e > f->digit_mask) { keep = 0; break; }
                r |= e << (v * f->bits);
            }
            if (keep) keep = mq_monom_contains(rows, r);
            for (slong v = 0; keep && v < f->nvars; v++) {
                ulong e = (packed[offset[f->nvars + v]] >> shift[f->nvars + v]) & field_mask;
                if (e > f->digit_mask) { keep = 0; break; }
                c |= e << (v * f->bits);
            }
            if (keep) keep = mq_monom_contains(cols, c);
        } else {
            nmod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
            keep = mq_filter_accepts(f, exp, target);
        }
        if (!keep) continue;
        if (out != i) {
            poly->coeffs[out] = poly->coeffs[i];
            memcpy(poly->exps + out * words, packed, (size_t) words * sizeof(ulong));
        }
        out++;
    }
    _nmod_mpoly_set_length(poly, out, ctx);
}

/* Retain FLINT's packed-exponent multiplication and compact immediately,
 * before this product enters a sum or feeds the next DP layer. Generating
 * pairs with push_term/sort is substantially slower on dense random MQ. */
static void mq_filtered_mul(nmod_mpoly_t out, const nmod_mpoly_t a,
                            const nmod_mpoly_t b, const nmod_mpoly_ctx_t ctx,
                            const mq_det_filter *f)
{
    nmod_mpoly_mul(out, a, b, ctx);
    if (f) mq_filter_poly(out, ctx, f, 0);
}

/* Layered minor DP, using colex column-subset ranks and adjacent layers.
 * The final expansion is parallelized by cofactor, followed by a tree sum.
 * The entry limit excludes arithmetic temporaries and matrix views.
 */
#ifdef DRSOLVE_DET_TESTING
/* Test-only observation; normal builds contain no callbacks. */
void drsolve_det_test_event(int event, slong size);
#endif
static int compute_nmod_mpoly_det_layered_dp(nmod_mpoly_t result, nmod_mpoly_t **matrix,
                       slong size, nmod_mpoly_ctx_t ctx, int use_parallel, slong limit,
                       const mq_det_filter *filter)
{
    ulong choose[FLINT_BITS][FLINT_BITS] = {{0}};
    nmod_mpoly_t *previous = NULL, *current = NULL;
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
        if (count > (size_t) -1 / sizeof(nmod_mpoly_t)) return 0;
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
                if (!nmod_mpoly_is_zero(matrix[row][col], ctx)) nonzero++;
            branches *= FLINT_MIN(nonzero, size - row);
            recursive_work += branches;
            if (recursive_work >= dp_work) break;
        }
        if (recursive_work + 12 * branches < dp_work) return 0;
    }

    if (g_dixon_verbose_level >= 2)
        printf("  determinant layered DP (%s): size=%ld, peak entries=%ld, limit=%ld, parallel=%s\n",
               "nmod", size, peak, limit, use_parallel ? "yes" : "no");

#ifdef DRSOLVE_DET_TESTING
    drsolve_det_test_event(0, size);
#endif

    for (slong k = 1; k <= size; k++) {
        slong count = k == size ? size : (slong) choose[size][k];
        const mq_det_filter *layer_filter = filter && k > filter->safe_linear_layers ? filter : NULL;
        int failed = 0;
        current = malloc((size_t) count * sizeof(nmod_mpoly_t));
        if (current == NULL) goto cleanup;
        for (current_count = 0; current_count < count; current_count++) {
            if (!(nmod_mpoly_init(current[current_count], ctx), 1)) goto cleanup;
        }

        if (k == 1) {
            for (slong i = 0; i < count; i++)
                nmod_mpoly_set(current[i], matrix[size - 1][i], ctx);
        } else {
            /* One writer per output, immutable previous layer, implicit
             * barrier before freeing it. No hash locks or duplicate work. */
#ifdef _OPENMP
            #pragma omp parallel if(use_parallel && count > 1 && !omp_in_parallel()) num_threads(FLINT_MIN(count, omp_get_max_threads())) reduction(|:failed)
#endif
            {
                nmod_mpoly_t product, sum;
                int product_ok = (nmod_mpoly_init(product, ctx), 1);
                int sum_ok = (nmod_mpoly_init(sum, ctx), 1);
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
                        if (nmod_mpoly_is_zero(matrix[0][index], ctx) ||
                            nmod_mpoly_is_zero(previous[size - 1 - index], ctx)) continue;
#ifdef DRSOLVE_DET_TESTING
                        drsolve_det_test_event(1, size);
#endif
                        mq_filtered_mul(current[index], matrix[0][index], previous[size - 1 - index], ctx, layer_filter);
                        if (index & 1) {
                            nmod_mpoly_zero(product, ctx);
                            nmod_mpoly_sub(current[index], product, current[index], ctx);
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

                    nmod_mpoly_zero(current[index], ctx);
                    for (slong j = 0; j < k; j++) {
                        ulong child = prefix[j] + suffix[j + 1];
                        if (nmod_mpoly_is_zero(matrix[size - k][cols[j]], ctx) ||
                            nmod_mpoly_is_zero(previous[child], ctx)) continue;
                        mq_filtered_mul(product, matrix[size - k][cols[j]], previous[child], ctx, layer_filter);
                        if (j & 1)
                            nmod_mpoly_sub(sum, current[index], product, ctx);
                        else
                            nmod_mpoly_add(sum, current[index], product, ctx);
                        nmod_mpoly_swap(current[index], sum, ctx);
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
                            nmod_mpoly_add(sum, current[left], current[left + stride], ctx);
                            nmod_mpoly_swap(current[left], sum, ctx);
                        }
                    }
                }
                if (product_ok) nmod_mpoly_clear(product, ctx);
                if (sum_ok) nmod_mpoly_clear(sum, ctx);
            }
        }
        if (failed) goto cleanup;
        if (filter && g_dixon_verbose_level >= 3) {
            slong terms = 0;
            /* At the root only current[0] remains live after tree reduction. */
            slong live = k == size ? 1 : count;
            for (slong i = 0; i < live; i++) terms += nmod_mpoly_length(current[i], ctx);
            printf("  MQ minor DP layer %ld/%ld: %ld live minors, %ld retained terms\n",
                   k, size, live, terms);
        }
        for (slong i = 0; i < previous_count; i++) nmod_mpoly_clear(previous[i], ctx);
        free(previous);
        previous = current;
        previous_count = current_count;
        current = NULL;
        current_count = 0;
    }
    nmod_mpoly_swap(result, previous[0], ctx);
    ok = 1;

cleanup:
    for (slong i = 0; i < current_count; i++) nmod_mpoly_clear(current[i], ctx);
    for (slong i = 0; i < previous_count; i++) nmod_mpoly_clear(previous[i], ctx);
    free(current);
    free(previous);
    return ok;
}

/* Method 0 backend. If a complete layer does not fit, expand one row and
 * retry DP on each child. Siblings execute sequentially, so their DP storage
 * never multiplies the entry budget; each child may still use layer threads.
 * The submatrix is a shallow, read-only view of the input polynomials. */
static void compute_nmod_mpoly_det_minor(nmod_mpoly_t result, nmod_mpoly_t **matrix,
                         slong size, nmod_mpoly_ctx_t ctx, int use_parallel, slong limit,
                       const mq_det_filter *filter)
{
    nmod_mpoly_t **rows = NULL, *entries = NULL;
    nmod_mpoly_t accum, child, product, sum;
    int accum_ok, child_ok, product_ok, sum_ok;
    size_t width;

    if (size <= 1 || ((!filter || size <= filter->safe_linear_layers) && (size <= 3 || limit <= 0))) {
        compute_nmod_mpoly_det_recursive(result, matrix, size, ctx);
        return;
    }
    if (compute_nmod_mpoly_det_layered_dp(result, matrix, size, ctx, use_parallel, limit, filter)) return;

    width = (size_t) (size - 1);
    if (width > (size_t) -1 / sizeof(*rows) ||
        width > (size_t) -1 / sizeof(*entries) / width) {
        compute_nmod_mpoly_det_recursive(result, matrix, size, ctx);
        return;
    }
    rows = malloc(width * sizeof(*rows));
    entries = malloc(width * width * sizeof(*entries));
    if (rows == NULL || entries == NULL) {
        free(rows);
        free(entries);
        compute_nmod_mpoly_det_recursive(result, matrix, size, ctx);
        return;
    }
    for (size_t i = 0; i < width; i++) rows[i] = entries + i * width;
    accum_ok = (nmod_mpoly_init(accum, ctx), 1);
    child_ok = (nmod_mpoly_init(child, ctx), 1);
    product_ok = (nmod_mpoly_init(product, ctx), 1);
    sum_ok = (nmod_mpoly_init(sum, ctx), 1);
    if (accum_ok && child_ok && product_ok && sum_ok) {
        nmod_mpoly_zero(accum, ctx);
        for (slong col = 0; col < size; col++) {
            if (nmod_mpoly_is_zero(matrix[0][col], ctx)) continue;
            for (slong i = 1; i < size; i++) {
                slong dst = 0;
                for (slong j = 0; j < size; j++) {
                    if (j == col) continue;
                    memcpy(&rows[i - 1][dst++], &matrix[i][j], sizeof(*entries));
                }
            }
            compute_nmod_mpoly_det_minor(child, rows, size - 1, ctx, use_parallel, limit, filter);
            if (nmod_mpoly_is_zero(child, ctx)) continue;
            mq_filtered_mul(product, matrix[0][col], child, ctx,
                            filter && size > filter->safe_linear_layers ? filter : NULL);
            if (col & 1)
                nmod_mpoly_sub(sum, accum, product, ctx);
            else
                nmod_mpoly_add(sum, accum, product, ctx);
            nmod_mpoly_swap(accum, sum, ctx);
        }
        nmod_mpoly_swap(result, accum, ctx);
    } else {
        compute_nmod_mpoly_det_recursive(result, matrix, size, ctx);
    }
    if (accum_ok) nmod_mpoly_clear(accum, ctx);
    if (child_ok) nmod_mpoly_clear(child, ctx);
    if (product_ok) nmod_mpoly_clear(product, ctx);
    if (sum_ok) nmod_mpoly_clear(sum, ctx);
    free(entries);
    free(rows);
}

static void compute_fq_det_nmod_minor_direct(fq_mvpoly_t *result,
                                              fq_mvpoly_t **matrix,
                                              slong size)
{
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    mp_limb_t p = fq_nmod_ctx_modulus(ctx)->mod.n;
    int use_parallel = (size >= PARALLEL_THRESHOLD && omp_get_max_threads() > 1);

    nmod_mpoly_ctx_t nmod_ctx;
    nmod_mpoly_ctx_init(nmod_ctx, total_vars, ORD_LEX, p);

    nmod_mpoly_t **nmod_matrix = (nmod_mpoly_t **) flint_malloc((size_t) size * sizeof(nmod_mpoly_t *));
    for (slong i = 0; i < size; i++) {
        nmod_matrix[i] = (nmod_mpoly_t *) flint_malloc((size_t) size * sizeof(nmod_mpoly_t));
    }
    fq_matrix_mvpoly_to_nmod_mpoly(nmod_matrix, matrix, size, nmod_ctx);

    nmod_mpoly_t det_nmod;
    nmod_mpoly_init(det_nmod, nmod_ctx);
    compute_nmod_mpoly_det_minor(det_nmod, nmod_matrix, size, nmod_ctx,
                                  use_parallel, g_dixon_det_cache_limit, NULL);

    fq_mvpoly_clear(result);
    nmod_mpoly_to_fq_mvpoly(result, det_nmod, nvars, npars, nmod_ctx, ctx);

    nmod_mpoly_clear(det_nmod, nmod_ctx);
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            nmod_mpoly_clear(nmod_matrix[i][j], nmod_ctx);
        }
        flint_free(nmod_matrix[i]);
    }
    flint_free(nmod_matrix);
    nmod_mpoly_ctx_clear(nmod_ctx);
}

/* Compute exactly the requested coefficient block, without claiming anything
 * about its rank. On ineligibility/budget failure, leave result untouched.
 * rows/cols are count contiguous exponent vectors of length size-1. */
int compute_fq_det_mq_projected(fq_mvpoly_t *result, fq_mvpoly_t **matrix,
                              slong size, const slong *rows,
                              const slong *cols, slong count)
{
    if (size < 2 || size >= FLINT_BITS || count <= 0 || count > MQ_FILTER_MAX_MONOMS ||
        !rows || !cols || !is_prime_field(matrix[0][0].ctx) ||
        matrix[0][0].nvars != 2 * (size - 1) || matrix[0][0].npars != 1) return 0;
    /* This backend intentionally accepts only divided-difference MQ matrices.
     * Check all entries, including the parameter degree, before packing. */
    for (slong i = 0; i < size; i++) for (slong j = 0; j < size; j++) {
        const fq_mvpoly_t *p = &matrix[i][j];
        if (p->nvars != 2 * (size - 1) || p->npars != 1) return 0;
        for (slong k = 0; k < p->nterms; k++) {
            slong degree = p->terms[k].par_exp ? p->terms[k].par_exp[0] : 0;
            if (p->terms[k].var_exp)
                for (slong v = 0; v < p->nvars; v++) degree += p->terms[k].var_exp[v];
            if (degree > (i == 0 ? 2 : 1)) return 0;
        }
    }
    for (slong i = 0; i < count * (size - 1); i++)
        if (rows[i] < 0 || cols[i] < 0 || rows[i] > size + 1 || cols[i] > size + 1) return 0;
    mq_det_filter filter;
    if (!mq_filter_init(&filter, size - 1, rows, cols, count)) return 0;
    const fq_nmod_ctx_struct *fq = matrix[0][0].ctx;
    slong nv = matrix[0][0].nvars;
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, nv + 1, ORD_LEX, fq_nmod_ctx_modulus(fq)->mod.n);
    slong row_safe = mq_safe_axis_layers(&filter, &filter.rows, matrix, size, 0);
    slong col_safe = mq_safe_axis_layers(&filter, &filter.cols, matrix, size, 1);
    filter.safe_linear_layers = FLINT_MIN(row_safe, col_safe);
    mq_filter_prepare_packed(&filter, ctx, size + 1);
    if (g_dixon_verbose_level >= 2)
        printf("  MQ filter: skipping %ld certified linear layers; packed axes=%s\n",
               filter.safe_linear_layers, filter.packed_bits ? "yes" : "no");
    nmod_mpoly_t **m = flint_malloc((size_t) size * sizeof(*m));
    for (slong i = 0; i < size; i++) m[i] = flint_malloc((size_t) size * sizeof(**m));
    fq_matrix_mvpoly_to_nmod_mpoly(m, matrix, size, ctx);
    nmod_mpoly_t det;
    nmod_mpoly_init(det, ctx);
    compute_nmod_mpoly_det_minor(det, m, size, ctx,
        size >= PARALLEL_THRESHOLD && omp_get_max_threads() > 1,
        g_dixon_det_cache_limit, &filter);
    mq_filter_poly(det, ctx, &filter, 1);
    nmod_mpoly_to_fq_mvpoly(result, det, nv, 1, ctx, fq);
    if (g_dixon_verbose_level >= 2)
        printf("  MQ projected minor DP: targets=%ld x %ld, closures=%ld x %ld, output=%ld terms\n",
               filter.row_targets.count, filter.col_targets.count,
               filter.rows.count, filter.cols.count, result->nterms);
    nmod_mpoly_clear(det, ctx);
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) nmod_mpoly_clear(m[i][j], ctx);
        flint_free(m[i]);
    }
    flint_free(m); nmod_mpoly_ctx_clear(ctx); mq_filter_clear(&filter);
    return 1;
}

// ============= Univariate Optimization Implementation =============

int is_univariate_matrix(fq_mvpoly_t **matrix, slong size) {
    if (size == 0) return 0;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    return (nvars == 1 && npars == 0);
}

void compute_fq_det_univariate_optimized(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    DET_PRINT("Using univariate polynomial matrix optimization for %ldx%ld matrix\n", size, size);
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    fq_mvpoly_init(result, 1, 0, ctx);
    
    timing_info_t start = start_timing();
    
    fq_nmod_poly_mat_t poly_mat;
    fq_nmod_poly_mat_init(poly_mat, size, size, ctx);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(poly_mat, i, j);
            fq_nmod_poly_zero(entry, ctx);
            
            for (slong k = 0; k < matrix[i][j].nterms; k++) {
                fq_monomial_t *term = &matrix[i][j].terms[k];
                slong degree = 0;
                if (term->var_exp && matrix[i][j].nvars > 0) {
                    degree = term->var_exp[0];
                }
                fq_nmod_poly_set_coeff(entry, degree, term->coeff, ctx);
            }
        }
    }
    
    fq_nmod_poly_t det_poly;
    fq_nmod_poly_init(det_poly, ctx);
    
    fq_nmod_poly_mat_det_iter(det_poly, poly_mat, ctx);
    
    timing_info_t conv_elapsed = end_timing(start);
    print_timing("Univariate matrix determinant", conv_elapsed);
    
    slong degree = fq_nmod_poly_degree(det_poly, ctx);
    if (degree >= 0) {
        for (slong d = 0; d <= degree; d++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, ctx);
            fq_nmod_poly_get_coeff(coeff, det_poly, d, ctx);
            
            if (!fq_nmod_is_zero(coeff, ctx)) {
                slong *var_exp = (slong*) flint_calloc(1, sizeof(slong));
                var_exp[0] = d;
                fq_mvpoly_add_term(result, var_exp, NULL, coeff);
                flint_free(var_exp);
            }
            
            fq_nmod_clear(coeff, ctx);
        }
    }
    
    fq_nmod_poly_clear(det_poly, ctx);
    fq_nmod_poly_mat_clear(poly_mat, ctx);
}

// ============= Conversion Functions Implementation =============

void fq_mvpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t mpoly, const fq_mvpoly_t *poly, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx) {
    fq_nmod_mpoly_zero(mpoly, mpoly_ctx);
    
    if (poly->nterms == 0) return;
    
    slong total_vars = poly->nvars + poly->npars;
    
    // Pre-allocate space for better performance
    fq_nmod_mpoly_fit_length(mpoly, poly->nterms, mpoly_ctx);
    
    for (slong i = 0; i < poly->nterms; i++) {
        ulong *exps = (ulong*) flint_calloc(total_vars, sizeof(ulong));
        
        if (poly->terms[i].var_exp && poly->nvars > 0) {
            for (slong j = 0; j < poly->nvars; j++) {
                exps[j] = (ulong)poly->terms[i].var_exp[j];
            }
        }
        
        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                exps[poly->nvars + j] = (ulong)poly->terms[i].par_exp[j];
            }
        }
        
        fq_nmod_mpoly_push_term_fq_nmod_ui(mpoly, poly->terms[i].coeff, exps, mpoly_ctx);
        flint_free(exps);
    }
    
    fq_nmod_mpoly_sort_terms(mpoly, mpoly_ctx);
    fq_nmod_mpoly_combine_like_terms(mpoly, mpoly_ctx);
}

void fq_nmod_mpoly_to_fq_mvpoly(fq_mvpoly_t *poly, const fq_nmod_mpoly_t mpoly,
                               slong nvars, slong npars, 
                               fq_nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx) {
    fq_mvpoly_init(poly, nvars, npars, ctx);
    slong nterms = fq_nmod_mpoly_length(mpoly, mpoly_ctx);
    if (nterms == 0) return;
    
    slong total_vars = fq_nmod_mpoly_ctx_nvars(mpoly_ctx);
    
    // Pre-allocate the terms array
    if (poly->alloc < nterms) {
        // First, clear any existing terms
        for (slong i = 0; i < poly->nterms; i++) {
            fq_nmod_clear(poly->terms[i].coeff, ctx);
            if (poly->terms[i].var_exp) flint_free(poly->terms[i].var_exp);
            if (poly->terms[i].par_exp) flint_free(poly->terms[i].par_exp);
        }
        
        poly->alloc = nterms;
        poly->terms = (fq_monomial_t*) flint_realloc(poly->terms, 
                                                      poly->alloc * sizeof(fq_monomial_t));
        
        // Initialize all term structures
        for (slong i = 0; i < poly->alloc; i++) {
            // Zero out the structure first
            memset(&poly->terms[i], 0, sizeof(fq_monomial_t));
        }
    }
   
    // Allocate a single exponent buffer for reading
    ulong *exp_buffer = (ulong*) flint_malloc(total_vars * sizeof(ulong));
    
    // Process all terms - but allocate individually for compatibility
    for (slong i = 0; i < nterms; i++) {
        // Initialize coefficient
        fq_nmod_init(poly->terms[i].coeff, ctx);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(poly->terms[i].coeff, mpoly, i, mpoly_ctx);
        // Get exponents for this term
        fq_nmod_mpoly_get_term_exp_ui(exp_buffer, mpoly, i, mpoly_ctx);
       
        // Allocate and set variable exponents
        if (nvars > 0) {
            poly->terms[i].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
            for (slong j = 0; j < nvars && j < total_vars; j++) {
                poly->terms[i].var_exp[j] = (slong)exp_buffer[j];
            }
        } else {
            poly->terms[i].var_exp = NULL;
        }
        // Allocate and set parameter exponents
        if (npars > 0 && total_vars > nvars) {
            poly->terms[i].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
            for (slong j = 0; j < npars && (nvars + j) < total_vars; j++) {
                poly->terms[i].par_exp[j] = (slong)exp_buffer[nvars + j];
            }
        } else {
            poly->terms[i].par_exp = NULL;
        }
    }

// Set the number of terms
    poly->nterms = nterms;
    if (g_field_equation_reduction) {
        fq_mvpoly_reduce_field_equation(poly);
    }
    
    // Cleanup
    flint_free(exp_buffer);
}

void fq_matrix_mvpoly_to_mpoly(fq_nmod_mpoly_t **mpoly_matrix, 
                              fq_mvpoly_t **mvpoly_matrix, 
                              slong size, 
                              fq_nmod_mpoly_ctx_t mpoly_ctx) {
    DET_PRINT("Converting %ld x %ld matrix\n", size, size);
    
    timing_info_t start = start_timing();
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_nmod_mpoly_init(mpoly_matrix[i][j], mpoly_ctx);
            fq_mvpoly_to_fq_nmod_mpoly(mpoly_matrix[i][j], &mvpoly_matrix[i][j], mpoly_ctx);
        }
    }
    
    timing_info_t elapsed = end_timing(start);
    print_timing("Matrix conversion", elapsed);
}

// ============= Optimized Determinant Computation Implementation =============

// Hand-optimized 3x3 determinant
void compute_det_3x3_optimized(fq_nmod_mpoly_t det, 
                              fq_nmod_mpoly_t **m,
                              fq_nmod_mpoly_ctx_t ctx) {
    fq_nmod_mpoly_t t1, t2, t3, t4, t5, t6, sum;
    
    // Initialize temporaries
    fq_nmod_mpoly_init(t1, ctx);
    fq_nmod_mpoly_init(t2, ctx);
    fq_nmod_mpoly_init(t3, ctx);
    fq_nmod_mpoly_init(t4, ctx);
    fq_nmod_mpoly_init(t5, ctx);
    fq_nmod_mpoly_init(t6, ctx);
    fq_nmod_mpoly_init(sum, ctx);
    
    // Compute 6 products in parallel if beneficial
    #pragma omp parallel sections if(omp_get_max_threads() > 2)
    {
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t1, m[1][1], m[2][2], ctx);
            poly_mul_dense_optimized(t1, m[0][0], t1, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t2, m[1][2], m[2][0], ctx);
            poly_mul_dense_optimized(t2, m[0][1], t2, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t3, m[1][0], m[2][1], ctx);
            poly_mul_dense_optimized(t3, m[0][2], t3, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t4, m[1][0], m[2][2], ctx);
            poly_mul_dense_optimized(t4, m[0][1], t4, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t5, m[1][1], m[2][0], ctx);
            poly_mul_dense_optimized(t5, m[0][2], t5, ctx);
        }
        #pragma omp section
        {
            fq_nmod_mpoly_mul(t6, m[1][2], m[2][1], ctx);
            poly_mul_dense_optimized(t6, m[0][0], t6, ctx);
        }
    }
    
    // Sum with signs
    fq_nmod_mpoly_add(sum, t1, t2, ctx);
    fq_nmod_mpoly_add(sum, sum, t3, ctx);
    fq_nmod_mpoly_sub(sum, sum, t4, ctx);
    fq_nmod_mpoly_sub(sum, sum, t5, ctx);
    fq_nmod_mpoly_sub(det, sum, t6, ctx);
    
    // Cleanup
    fq_nmod_mpoly_clear(t1, ctx);
    fq_nmod_mpoly_clear(t2, ctx);
    fq_nmod_mpoly_clear(t3, ctx);
    fq_nmod_mpoly_clear(t4, ctx);
    fq_nmod_mpoly_clear(t5, ctx);
    fq_nmod_mpoly_clear(t6, ctx);
    fq_nmod_mpoly_clear(sum, ctx);
}

// Recursive determinant with optimizations
void compute_fq_nmod_mpoly_det_recursive(fq_nmod_mpoly_t det_result, 
                                        fq_nmod_mpoly_t **mpoly_matrix, 
                                        slong size, 
                                        fq_nmod_mpoly_ctx_t mpoly_ctx) {
    if (size <= 0) {
        fq_nmod_mpoly_one(det_result, mpoly_ctx);
        return;
    }
    
    if (size == 1) {
        fq_nmod_mpoly_set(det_result, mpoly_matrix[0][0], mpoly_ctx);
        return;
    }
    
    if (size == 2) {
        fq_nmod_mpoly_t ad, bc;
        fq_nmod_mpoly_init(ad, mpoly_ctx);
        fq_nmod_mpoly_init(bc, mpoly_ctx);
        
        poly_mul_dense_optimized(ad, mpoly_matrix[0][0], mpoly_matrix[1][1], mpoly_ctx);
        poly_mul_dense_optimized(bc, mpoly_matrix[0][1], mpoly_matrix[1][0], mpoly_ctx);
        fq_nmod_mpoly_sub(det_result, ad, bc, mpoly_ctx);
        
        fq_nmod_mpoly_clear(ad, mpoly_ctx);
        fq_nmod_mpoly_clear(bc, mpoly_ctx);
        return;
    }
    
    if (size == 3) {
        compute_det_3x3_optimized(det_result, mpoly_matrix, mpoly_ctx);
        return;
    }
    
    // General case: Laplace expansion
    fq_nmod_mpoly_zero(det_result, mpoly_ctx);
    
    fq_nmod_mpoly_t temp_result, cofactor, subdet;
    fq_nmod_mpoly_init(temp_result, mpoly_ctx);
    fq_nmod_mpoly_init(cofactor, mpoly_ctx);
    fq_nmod_mpoly_init(subdet, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            continue;
        }
        
        // Create submatrix
        fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
        for (slong i = 0; i < size-1; i++) {
            submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
            }
        }
        
        // Fill submatrix
        for (slong i = 1; i < size; i++) {
            slong sub_j = 0;
            for (slong j = 0; j < size; j++) {
                if (j != col) {
                    fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                    sub_j++;
                }
            }
        }
        
        // Recursive computation
        compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
        
        // Compute cofactor
        poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
        
        // Add/subtract to result
        if (col % 2 == 0) {
            fq_nmod_mpoly_add(temp_result, det_result, cofactor, mpoly_ctx);
        } else {
            fq_nmod_mpoly_sub(temp_result, det_result, cofactor, mpoly_ctx);
        }
        fq_nmod_mpoly_set(det_result, temp_result, mpoly_ctx);
        
        // Cleanup submatrix
        for (slong i = 0; i < size-1; i++) {
            for (slong j = 0; j < size-1; j++) {
                fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
            }
            flint_free(submatrix[i]);
        }
        flint_free(submatrix);
    }
    
    fq_nmod_mpoly_clear(temp_result, mpoly_ctx);
    fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
    fq_nmod_mpoly_clear(subdet, mpoly_ctx);
}

void compute_fq_nmod_mpoly_det_parallel_optimized(fq_nmod_mpoly_t det_result, 
                                                  fq_nmod_mpoly_t **mpoly_matrix, 
                                                  slong size, 
                                                  fq_nmod_mpoly_ctx_t mpoly_ctx,
                                                  slong depth) {
    if (size < PARALLEL_THRESHOLD || depth >= MAX_PARALLEL_DEPTH) {
        compute_fq_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    // DET_PRINT("Parallel computation for %ld x %ld matrix (depth %ld)\n", size, size, depth);
    
    if (size <= 3) {
        compute_fq_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    fq_nmod_mpoly_zero(det_result, mpoly_ctx);
    
    // Count non-zero entries in first row
    slong nonzero_count = 0;
    for (slong col = 0; col < size; col++) {
        if (!fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
            nonzero_count++;
        }
    }
    
    if (nonzero_count < 2) {
        compute_fq_nmod_mpoly_det_recursive(det_result, mpoly_matrix, size, mpoly_ctx);
        return;
    }
    
    // Allocate space for partial results
    fq_nmod_mpoly_t *partial_results = (fq_nmod_mpoly_t*) flint_malloc(size * sizeof(fq_nmod_mpoly_t));
    for (slong i = 0; i < size; i++) {
        fq_nmod_mpoly_init(partial_results[i], mpoly_ctx);
        fq_nmod_mpoly_zero(partial_results[i], mpoly_ctx);
    }
    
    // Determine parallelism strategy based on depth
    if (depth == 0) {
        // First level: use parallel for with nested parallelism enabled
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, omp_get_max_threads()))
        for (slong col = 0; col < size; col++) {
            if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            fq_nmod_mpoly_t cofactor, subdet;
            fq_nmod_mpoly_init(cofactor, mpoly_ctx);
            fq_nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create submatrix
            fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            // Fill submatrix
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Recursive call - this will use nested parallelism at depth 1
            compute_fq_nmod_mpoly_det_parallel_optimized(subdet, submatrix, size-1, mpoly_ctx, depth+1);
            
            // Compute cofactor
            poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
            // Store with sign
            if (col % 2 == 0) {
                fq_nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                fq_nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
            fq_nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    } else if (depth == 1 && size >= PARALLEL_THRESHOLD) {
        // Second level: also use parallel for, but with fewer threads
        slong max_threads_level2 = FLINT_MAX(1, omp_get_max_threads() / size);
        
        #pragma omp parallel for schedule(static) num_threads(FLINT_MIN(nonzero_count, max_threads_level2)) if(nonzero_count >= 3)
        for (slong col = 0; col < size; col++) {
            if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            fq_nmod_mpoly_t cofactor, subdet;
            fq_nmod_mpoly_init(cofactor, mpoly_ctx);
            fq_nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create submatrix
            fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            // Fill submatrix
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // At depth > 1, use sequential computation
            compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor
            poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            
            // Store with sign
            if (col % 2 == 0) {
                fq_nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                fq_nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
            fq_nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    } else {
        // Sequential fallback for deeper levels or small matrices
        for (slong col = 0; col < size; col++) {
            if (fq_nmod_mpoly_is_zero(mpoly_matrix[0][col], mpoly_ctx)) {
                continue;
            }
            
            fq_nmod_mpoly_t cofactor, subdet;
            fq_nmod_mpoly_init(cofactor, mpoly_ctx);
            fq_nmod_mpoly_init(subdet, mpoly_ctx);
            
            // Create and fill submatrix
            fq_nmod_mpoly_t **submatrix = (fq_nmod_mpoly_t**) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t*));
            for (slong i = 0; i < size-1; i++) {
                submatrix[i] = (fq_nmod_mpoly_t*) flint_malloc((size-1) * sizeof(fq_nmod_mpoly_t));
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_init(submatrix[i][j], mpoly_ctx);
                }
            }
            
            for (slong i = 1; i < size; i++) {
                slong sub_j = 0;
                for (slong j = 0; j < size; j++) {
                    if (j != col) {
                        fq_nmod_mpoly_set(submatrix[i-1][sub_j], mpoly_matrix[i][j], mpoly_ctx);
                        sub_j++;
                    }
                }
            }
            
            // Sequential computation
            compute_fq_nmod_mpoly_det_recursive(subdet, submatrix, size-1, mpoly_ctx);
            
            // Compute cofactor and store
            poly_mul_dense_optimized(cofactor, mpoly_matrix[0][col], subdet, mpoly_ctx);
            if (col % 2 == 0) {
                fq_nmod_mpoly_set(partial_results[col], cofactor, mpoly_ctx);
            } else {
                fq_nmod_mpoly_neg(partial_results[col], cofactor, mpoly_ctx);
            }
            
            // Cleanup
            for (slong i = 0; i < size-1; i++) {
                for (slong j = 0; j < size-1; j++) {
                    fq_nmod_mpoly_clear(submatrix[i][j], mpoly_ctx);
                }
                flint_free(submatrix[i]);
            }
            flint_free(submatrix);
            
            fq_nmod_mpoly_clear(cofactor, mpoly_ctx);
            fq_nmod_mpoly_clear(subdet, mpoly_ctx);
        }
    }
    
    // Sum results (sequential to avoid race conditions)
    fq_nmod_mpoly_t temp_sum;
    fq_nmod_mpoly_init(temp_sum, mpoly_ctx);
    
    for (slong col = 0; col < size; col++) {
        if (!fq_nmod_mpoly_is_zero(partial_results[col], mpoly_ctx)) {
            fq_nmod_mpoly_add(temp_sum, det_result, partial_results[col], mpoly_ctx);
            fq_nmod_mpoly_set(det_result, temp_sum, mpoly_ctx);
        }
        fq_nmod_mpoly_clear(partial_results[col], mpoly_ctx);
    }
    
    fq_nmod_mpoly_clear(temp_sum, mpoly_ctx);
    flint_free(partial_results);
}

void compute_fq_det_huang_interpolation(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    /*
    // Check if we're in a prime field
    if (!is_prime_field(ctx)) {
        printf("ERROR: sparse interpolation requires prime field\n");
        compute_fq_det_recursive(result, matrix, size);
        return;
    }
    */
    DET_PRINT("Computing %ldx%ld determinant via sparse interpolation\n", size, size);
    DET_PRINT("Variables: %ld, Parameters: %ld\n", nvars, npars);
    
    // Get the prime
    mp_limb_t p = fq_nmod_ctx_modulus(ctx)->mod.n;
    
    // Create nmod_mpoly context - FIX: Add modulus parameter
    nmod_mpoly_ctx_t nmod_ctx;
    nmod_mpoly_ctx_init(nmod_ctx, nvars + npars, ORD_LEX, p);  // Added p as 4th parameter
    
    // Convert matrix to poly_mat_t format for huang.h
    poly_mat_t huang_mat;
    poly_mat_init(&huang_mat, size, size, nmod_ctx);
    
    // Convert each entry
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            nmod_mpoly_t temp;
            nmod_mpoly_init(temp, nmod_ctx);
            
            // Convert fq_mvpoly to nmod_mpoly
            fq_mvpoly_to_nmod_mpoly(temp, &matrix[i][j], nmod_ctx);
            poly_mat_entry_set(&huang_mat, i, j, temp, nmod_ctx);
            
            nmod_mpoly_clear(temp, nmod_ctx);
        }
    }
    
    // Call Sparse.pdf-style determinant probing
    nmod_mpoly_t det_nmod;
    nmod_mpoly_init(det_nmod, nmod_ctx);
    
    timing_info_t huang_start = start_timing();
    ComputePolyMatrixDet(det_nmod, &huang_mat, nvars + npars, p, nmod_ctx);
    timing_info_t huang_elapsed = end_timing(huang_start);
    print_timing("sparse interpolation", huang_elapsed);
    
    // Convert result back to fq_mvpoly
    nmod_mpoly_to_fq_mvpoly(result, det_nmod, nvars, npars, nmod_ctx, ctx);
    
    DET_PRINT("Final result: %ld terms\n", result->nterms);
    
    // Cleanup
    poly_mat_clear(&huang_mat, nmod_ctx);
    nmod_mpoly_clear(det_nmod, nmod_ctx);
    nmod_mpoly_ctx_clear(nmod_ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    print_timing("Total sparse interpolation method", total_elapsed);
}

/* Main function extracted from the #else branch */
static void compute_fq_det_unified_interface_impl(fq_mvpoly_t *result,
                                                  fq_mvpoly_t **matrix,
                                                  slong size,
                                                  int method) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }

    timing_info_t total_start = start_timing();
    
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    slong total_vars = nvars + npars;
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    slong max_threads = omp_get_max_threads();
    
    DET_PRINT("Computing %ldx%ld determinant%s (OpenMP: %ld threads available)\n",
              size, size,
              (method == DET_METHOD_KRONECKER_NMOD) ? " via Bareiss" :
              (method == DET_METHOD_BALANCED_SPLIT) ? " via experimental balanced split" : "",
              max_threads);

    fq_mvpoly_init(result, nvars, npars, ctx);

    // Check for univariate optimization
    if (is_univariate_matrix(matrix, size) && size >= UNIVARIATE_THRESHOLD) {
        DET_PRINT("Detected univariate matrix, using specialized optimization\n");
        compute_fq_det_univariate_optimized(result, matrix, size);
        
        timing_info_t total_elapsed = end_timing(total_start);
        //print_timing("Total univariate computation", total_elapsed);
        return;
    }

    if (is_prime_field(ctx) && method == 0) {
        if (g_dixon_verbose_level >= 2) {
            printf("  prime-field recursive determinant: using direct nmod_mpoly minor DP path\n");
        }
        compute_fq_det_nmod_minor_direct(result, matrix, size);
        timing_info_t total_elapsed = end_timing(total_start);
        (void) total_elapsed;
        return;
    }

    // ===== USE UNIFIED INTERFACE =====
    DET_PRINT("Using unified multivariate polynomial interface\n");
    
    // DEBUG: Print original matrix
    //debug_print_fq_mvpoly_matrix(matrix, size, "ORIGINAL");
    // Step 1: Create field context wrapper
    field_ctx_t field_ctx;
    field_ctx_init(&field_ctx, ctx);  // Properly initialize the field context

    // Verify the field type detection
    mp_limb_t p = fq_nmod_ctx_modulus(ctx)->mod.n;
    slong degree = fq_nmod_ctx_degree(ctx);
    DET_PRINT("Field: p=%lu, degree=%ld, detected type=%d\n", p, degree, field_ctx.field_id);
    
    // Debug: print field context details
    //printf("Field context details:\n");
    //printf("  field_id: %d\n", field_ctx.field_id);
    //printf("  modulus: ");
    //fq_nmod_ctx_modulus_print_pretty(ctx, "t");
    //printf("\n");
    //printf("  degree: %ld\n", degree);
    //printf("  characteristic: %lu\n", p);
    // Step 2: Create unified multivariate context
    unified_mpoly_ctx_t unified_ctx = unified_mpoly_ctx_init(total_vars, ORD_LEX, &field_ctx);
    if (!unified_ctx) {
        printf("ERROR: Failed to create unified context\n");
        field_ctx_clear(&field_ctx);
        return;
    }
    // Step 3: Allocate unified polynomial matrix
    unified_mpoly_t **unified_matrix = (unified_mpoly_t**) malloc(size * sizeof(unified_mpoly_t*));
    for (slong i = 0; i < size; i++) {
        unified_matrix[i] = (unified_mpoly_t*) malloc(size * sizeof(unified_mpoly_t));
        for (slong j = 0; j < size; j++) {
            unified_matrix[i][j] = unified_mpoly_init(unified_ctx);
            if (!unified_matrix[i][j]) {
                printf("ERROR: Failed to initialize unified polynomial at (%ld,%ld)\n", i, j);
                for (slong ii = 0; ii <= i; ii++) {
                    slong limit = (ii == i) ? j : size;
                    for (slong jj = 0; jj < limit; jj++) {
                        unified_mpoly_clear(unified_matrix[ii][jj]);
                    }
                    free(unified_matrix[ii]);
                }
                free(unified_matrix);
                unified_mpoly_ctx_clear(unified_ctx);
                field_ctx_clear(&field_ctx);
                return;
            }
        }
    }
    // Step 4: Convert fq_mvpoly matrix to unified_mpoly matrix - Fast batch version
    timing_info_t convert_start = start_timing();
    DET_PRINT("Converting %ldx%ld matrix to unified format\n", size, size);
    
    // Create mpoly context for intermediate conversion
    fq_nmod_mpoly_ctx_t mpoly_ctx;
    fq_nmod_mpoly_ctx_init(mpoly_ctx, total_vars, ORD_LEX, ctx);
    
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            unified_mpoly_zero(unified_matrix[i][j]);
            
            fq_mvpoly_t *src_poly = &matrix[i][j];
            
            if (src_poly->nterms == 0) {
                continue; // Skip empty polynomials
            }
            
            // Method 1: Use push_term interface for batch processing
            fq_nmod_mpoly_t temp_poly;
            fq_nmod_mpoly_init(temp_poly, mpoly_ctx);
            
            // Push all terms at once
            for (slong t = 0; t < src_poly->nterms; t++) {
                if (fq_nmod_is_zero(src_poly->terms[t].coeff, ctx)) {
                    continue;
                }
                
                // Build combined exponent vector
                ulong *exp = (ulong*) flint_calloc(total_vars, sizeof(ulong));
                
                // Copy variable exponents
                if (src_poly->terms[t].var_exp && nvars > 0) {
                    for (slong v = 0; v < nvars; v++) {
                        exp[v] = (ulong)src_poly->terms[t].var_exp[v];
                    }
                }
                
                // Copy parameter exponents  
                if (src_poly->terms[t].par_exp && npars > 0) {
                    for (slong p = 0; p < npars; p++) {
                        exp[nvars + p] = (ulong)src_poly->terms[t].par_exp[p];
                    }
                }
                
                // Push term directly - much faster than repeated add
                fq_nmod_mpoly_push_term_fq_nmod_ui(temp_poly, src_poly->terms[t].coeff, exp, mpoly_ctx);
                
                flint_free(exp);
            }
            
            // Sort and combine like terms once at the end
            fq_nmod_mpoly_sort_terms(temp_poly, mpoly_ctx);
            fq_nmod_mpoly_combine_like_terms(temp_poly, mpoly_ctx);
            
            // Convert from fq_nmod_mpoly to unified_mpoly
            slong temp_length = fq_nmod_mpoly_length(temp_poly, mpoly_ctx);
            
            for (slong k = 0; k < temp_length; k++) {
                // Get coefficient
                fq_nmod_t coeff;
                fq_nmod_init(coeff, ctx);
                fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, temp_poly, k, mpoly_ctx);
                
                // Get exponent vector
                ulong *exp = (ulong*) flint_calloc(total_vars, sizeof(ulong));
                fq_nmod_mpoly_get_term_exp_ui(exp, temp_poly, k, mpoly_ctx);
                
                // Convert coefficient to field element
                field_elem_u field_coeff;
                void *ctx_ptr = (field_ctx.field_id == FIELD_ID_NMOD) ?
                               (void*)&field_ctx.ctx.nmod_ctx :
                               (field_ctx.field_id == FIELD_ID_FQ_ZECH) ?
                               (void*)field_ctx.ctx.zech_ctx :
                               (void*)field_ctx.ctx.fq_ctx;
                field_init_elem(&field_coeff, field_ctx.field_id, ctx_ptr);
                fq_nmod_to_field_elem(&field_coeff, coeff, &field_ctx);
                
                // Set in unified polynomial
                unified_mpoly_set_coeff_ui(unified_matrix[i][j], &field_coeff, exp);
                
                // Cleanup
                field_clear_elem(&field_coeff, field_ctx.field_id, ctx_ptr);
                fq_nmod_clear(coeff, ctx);
                flint_free(exp);
            }
            
            fq_nmod_mpoly_clear(temp_poly, mpoly_ctx);
        }
    }
    
    fq_nmod_mpoly_ctx_clear(mpoly_ctx);
    timing_info_t convert_elapsed = end_timing(convert_start);
    // DEBUG: Print converted matrix
    //debug_print_unified_matrix(unified_matrix, size, "CONVERTED");
    // Step 5: Enable optimizations if applicable
    /*
    if (field_ctx.field_id == FIELD_ID_GF28) {
        unified_mpoly_enable_optimizations(FIELD_ID_GF28, 1);
        DET_PRINT("Enabled GF(2^8) optimizations\n");
    } else if (field_ctx.field_id == FIELD_ID_GF2128) {
        unified_mpoly_enable_optimizations(FIELD_ID_GF2128, 1);
        DET_PRINT("Enabled GF(2^128) optimizations\n");
    }
    */
    // Step 6: Compute determinant using unified interface
    unified_mpoly_t det_unified = unified_mpoly_init(unified_ctx);
    if (!det_unified) {
        printf("ERROR: Failed to initialize determinant polynomial\n");
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                unified_mpoly_clear(unified_matrix[i][j]);
            }
            free(unified_matrix[i]);
        }
        free(unified_matrix);
        unified_mpoly_ctx_clear(unified_ctx);
        field_ctx_clear(&field_ctx);
        return;
    }

    timing_info_t det_start = start_timing();
    int use_parallel = (size >= PARALLEL_THRESHOLD && max_threads > 1);
    compute_unified_mpoly_det_with_method(det_unified, unified_matrix, size, unified_ctx,
                                          use_parallel, method);
    timing_info_t det_elapsed = end_timing(det_start);
    //print_timing("Determinant computation (unified)", det_elapsed);

    DET_PRINT("Unified determinant has %ld terms\n", unified_mpoly_length(det_unified));
    // Step 7: Convert result back to fq_mvpoly (with debugging)
    timing_info_t result_start = start_timing();
    fq_mvpoly_clear(result);  // Clear the initialization from the beginning
    fq_mvpoly_init(result, nvars, npars, ctx);
    
   // printf("\n--- Converting result back to fq_mvpoly ---\n");
    //printf("Field type for result conversion: %d\n", field_ctx.field_id);

    // Convert based on field type
    if (field_ctx.field_id == FIELD_ID_NMOD || is_prime_field(ctx)) {
        //printf("Using prime field result conversion\n");
        // For prime fields, convert from nmod_mpoly
        nmod_mpoly_struct *nmod_poly = GET_NMOD_POLY(det_unified);
        nmod_mpoly_ctx_struct *nmod_ctx = &(unified_ctx->ctx.nmod_ctx);

        // Use the existing conversion function
        slong nterms = nmod_mpoly_length(nmod_poly, nmod_ctx);
        DET_PRINT("Converting nmod_mpoly with %ld terms\n", nterms);

        // Pre-allocate the result polynomial to avoid reallocations
        if (result->alloc < nterms) {
            result->alloc = nterms + nterms/10; // Add 10% extra space
            result->terms = (fq_monomial_t*) flint_realloc(result->terms, 
                                                            result->alloc * sizeof(fq_monomial_t));
        }

        // Allocate a temporary exponent array once
        ulong *exp_ui = (ulong*) flint_malloc(total_vars * sizeof(ulong));

        // Batch convert all terms
        result->nterms = 0;
        for (slong i = 0; i < nterms; i++) {
            // Get coefficient
            mp_limb_t coeff_ui = nmod_mpoly_get_term_coeff_ui(nmod_poly, i, nmod_ctx);
            //printf("Result term %ld: coeff_ui = %lu\n", i, coeff_ui);

            // Get exponents
            nmod_mpoly_get_term_exp_ui(exp_ui, nmod_poly, i, nmod_ctx);
            /*
            printf("  exponents: ");
            for (slong k = 0; k < total_vars; k++) {
                printf("%lu ", exp_ui[k]);
            }
            printf("\n");
            */
            // Directly set the term without using fq_mvpoly_add_term
            fq_nmod_init(result->terms[result->nterms].coeff, ctx);
            fq_nmod_set_ui(result->terms[result->nterms].coeff, coeff_ui, ctx);

            // Split exponents
            if (nvars > 0) {
                result->terms[result->nterms].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                for (slong v = 0; v < nvars; v++) {
                    result->terms[result->nterms].var_exp[v] = (slong)exp_ui[v];
                }
            } else {
                result->terms[result->nterms].var_exp = NULL;
            }

            if (npars > 0) {
                result->terms[result->nterms].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                for (slong p = 0; p < npars; p++) {
                    result->terms[result->nterms].par_exp[p] = (slong)exp_ui[nvars + p];
                }
            } else {
                result->terms[result->nterms].par_exp = NULL;
            }

            result->nterms++;

            // Progress indicator for large conversions
            if (i > 0 && i % 10000 == 0) {
                DET_PRINT("Converted %ld/%ld terms...\n", i, nterms);
            }
        }

        flint_free(exp_ui);
        DET_PRINT("Conversion complete: %ld terms\n", result->nterms);

    } else if (field_ctx.field_id == FIELD_ID_FQ_ZECH) {
        //printf("Using Zech field result conversion\n");
        // For Zech logarithm fields, convert from fq_zech_mpoly
        fq_zech_mpoly_struct *zech_poly = GET_ZECH_POLY(det_unified);
        fq_zech_mpoly_ctx_struct *zech_ctx = &(unified_ctx->ctx.zech_ctx);

        slong nterms = fq_zech_mpoly_length(zech_poly, zech_ctx);
        DET_PRINT("Converting fq_zech_mpoly with %ld terms\n", nterms);

        if (nterms > 0) {
            // Pre-allocate the result polynomial
            if (result->alloc < nterms) {
                result->alloc = nterms + nterms/10; // Add 10% extra space
                result->terms = (fq_monomial_t*) flint_realloc(result->terms, 
                                                                result->alloc * sizeof(fq_monomial_t));
            }

            // Allocate temporary storage
            ulong *exp_ui = (ulong*) flint_malloc(total_vars * sizeof(ulong));
            fq_zech_t zech_coeff;
            fq_zech_init(zech_coeff, field_ctx.ctx.zech_ctx);

            // Convert each term
            result->nterms = 0;
            for (slong i = 0; i < nterms; i++) {
                // Get coefficient from Zech representation
                fq_zech_mpoly_get_term_coeff_fq_zech(zech_coeff, zech_poly, i, zech_ctx);

                // Get exponents
                fq_zech_mpoly_get_term_exp_ui(exp_ui, zech_poly, i, zech_ctx);

                // Initialize result coefficient
                fq_nmod_init(result->terms[result->nterms].coeff, ctx);

                // Convert Zech coefficient to fq_nmod
                fq_zech_get_fq_nmod(result->terms[result->nterms].coeff, zech_coeff, field_ctx.ctx.zech_ctx);

                // Split exponents
                if (nvars > 0) {
                    result->terms[result->nterms].var_exp = (slong*) flint_calloc(nvars, sizeof(slong));
                    for (slong v = 0; v < nvars; v++) {
                        result->terms[result->nterms].var_exp[v] = (slong)exp_ui[v];
                    }
                } else {
                    result->terms[result->nterms].var_exp = NULL;
                }

                if (npars > 0) {
                    result->terms[result->nterms].par_exp = (slong*) flint_calloc(npars, sizeof(slong));
                    for (slong p = 0; p < npars; p++) {
                        result->terms[result->nterms].par_exp[p] = (slong)exp_ui[nvars + p];
                    }
                } else {
                    result->terms[result->nterms].par_exp = NULL;
                }

                result->nterms++;

                // Progress indicator for large conversions
                if (i > 0 && i % 10000 == 0) {
                    DET_PRINT("Converted %ld/%ld terms from Zech...\n", i, nterms);
                }
            }

            // Cleanup
            flint_free(exp_ui);
            fq_zech_clear(zech_coeff, field_ctx.ctx.zech_ctx);
            DET_PRINT("Zech conversion complete: %ld terms\n", result->nterms);
        }

    } else {
        //printf("Using extension field result conversion\n");
        // For extension fields, convert from fq_nmod_mpoly
        fq_nmod_mpoly_struct *fq_poly = GET_FQ_POLY(det_unified);
        fq_nmod_mpoly_ctx_struct *fq_ctx = &(unified_ctx->ctx.fq_ctx);

        // Use the existing conversion function if available
        slong nterms = fq_nmod_mpoly_length(fq_poly, fq_ctx);
        DET_PRINT("Converting fq_nmod_mpoly with %ld terms\n", nterms);

        if (nterms > 0) {
            // Use existing conversion function
            fq_mvpoly_clear(result);
            fq_nmod_mpoly_to_fq_mvpoly(result, fq_poly, nvars, npars, fq_ctx, ctx);
        }
    }
    
    timing_info_t result_elapsed = end_timing(result_start);
    //print_timing("Result conversion from unified format", result_elapsed);

    DET_PRINT("Final result: %ld terms\n", result->nterms);

    // DEBUG: Print final result
    //debug_print_fq_mvpoly_matrix(&result, 1, "FINAL RESULT");
    // Step 8: Cleanup
/*
    // Disable optimizations
    if (field_ctx.field_id == FIELD_ID_GF28) {
        unified_mpoly_enable_optimizations(FIELD_ID_GF28, 0);
    } else if (field_ctx.field_id == FIELD_ID_GF2128) {
        unified_mpoly_enable_optimizations(FIELD_ID_GF2128, 0);
    }
*/
    // Free unified matrix
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            unified_mpoly_clear(unified_matrix[i][j]);
        }
        free(unified_matrix[i]);
    }
    free(unified_matrix);

    // Clear unified polynomial and context
    unified_mpoly_clear(det_unified);
    unified_mpoly_ctx_clear(unified_ctx);
    field_ctx_clear(&field_ctx);
    
    timing_info_t total_elapsed = end_timing(total_start);
    (void) total_elapsed;
}

void compute_fq_det_unified_interface(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    compute_fq_det_unified_interface_impl(result, matrix, size, 0);
}

void compute_fq_det_bareiss(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    compute_fq_det_unified_interface_impl(result, matrix, size, 4);
}

void compute_fq_det_balanced_split_experimental(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    compute_fq_det_unified_interface_impl(result, matrix, size, 6);
}
// ============= Main Interface with Algorithm Selection Implementation =============

void compute_fq_det_recursive_flint(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    if (size <= 0) {
        fq_mvpoly_init(result, matrix[0][0].nvars, matrix[0][0].npars, matrix[0][0].ctx);
        return;
    }
    
    timing_info_t total_start = start_timing();
    
    slong nvars = matrix[0][0].nvars;
    slong npars = matrix[0][0].npars;
    const fq_nmod_ctx_struct *ctx = matrix[0][0].ctx;
    
    // Choose algorithm based on configuration
    #if DET_ALGORITHM == DET_ALGORITHM_INTERPOLATION
    {
        printf("Using multivariate interpolation algorithm\n");
        
        // Include the interpolation header if not already included
        #ifndef FQ_NMOD_INTERPOLATION_OPTIMIZED_H
        #include "fq_multivariate_interpolation.h"
        #endif
        slong total_vars = nvars + npars;
        slong *var_bounds = (slong*) malloc(total_vars * sizeof(slong));
        compute_kronecker_bounds(var_bounds, matrix, size, nvars, npars);
        // Use interpolation algorithm
        fq_compute_det_by_interpolation_optimized(result, matrix, size, 
                                                 nvars, npars, ctx, var_bounds);
        return;
    }
    #elif DET_ALGORITHM == DET_ALGORITHM_KRONECKER
    {
        printf("Using Kronecker+HNF algorithm\n");
        compute_fq_det_kronecker(result, matrix, size);
        return;
    }
    #elif DET_ALGORITHM == DET_ALGORITHM_POLY_RECURSIVE
    {
        printf("Using polynomial recursive algorithm\n");
        compute_fq_det_poly_recursive(result, matrix, size);
        return;
    }
    #elif DET_ALGORITHM == DET_ALGORITHM_HUANG
    {
        // Check if we're in a prime field
        if (!is_prime_field(ctx)) {
            // Fall back to recursive algorithm (use code below)
        } else {
            compute_fq_det_huang_interpolation(result, matrix, size);
            return;
        }
    }
    #else
    {
        compute_fq_det_unified_interface(result, matrix, size);
    }
    #endif
}

// Compatibility interface
void compute_fq_det_recursive(fq_mvpoly_t *result, fq_mvpoly_t **matrix, slong size) {
    compute_fq_det_recursive_flint(result, matrix, size);
}

static inline ulong reduce_exp_field_ui(ulong e, ulong q) {
    if (e == 0 || e < q) return e;
    return ((e - 1) % (q - 1)) + 1;
}

static ulong field_size_q_from_fq_ctx(const fq_nmod_ctx_t ctx) {
    mp_limb_t p = fq_nmod_ctx_prime(ctx);
    slong d = fq_nmod_ctx_degree(ctx);
    ulong q = 1;
    for (slong i = 0; i < d; i++) {
        if (q > WORD_MAX / p) return WORD_MAX;
        q *= p;
    }
    return q;
}

static void fq_nmod_poly_reduce_field_equation_inplace(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx) {
    if (!g_field_equation_reduction) return;
    slong deg = fq_nmod_poly_degree(poly, ctx);
    if (deg <= 0) return;
    ulong q = field_size_q_from_fq_ctx(ctx);
    if (q <= 1 || q == WORD_MAX) return;
    if ((ulong)deg < q) return;

    fq_nmod_poly_t reduced;
    fq_nmod_poly_init(reduced, ctx);
    fq_nmod_poly_zero(reduced, ctx);
    fq_nmod_t coeff, acc;
    fq_nmod_init(coeff, ctx);
    fq_nmod_init(acc, ctx);

    for (slong i = 0; i <= deg; i++) {
        fq_nmod_poly_get_coeff(coeff, poly, i, ctx);
        if (fq_nmod_is_zero(coeff, ctx)) continue;
        slong tgt = (i == 0 || (ulong)i < q) ? i : (slong)(((ulong)(i - 1) % (q - 1)) + 1);
        fq_nmod_poly_get_coeff(acc, reduced, tgt, ctx);
        fq_nmod_add(acc, acc, coeff, ctx);
        fq_nmod_poly_set_coeff(reduced, tgt, acc, ctx);
    }

    fq_nmod_poly_set(poly, reduced, ctx);
    fq_nmod_clear(coeff, ctx);
    fq_nmod_clear(acc, ctx);
    fq_nmod_poly_clear(reduced, ctx);
}

static void fq_nmod_mpoly_reduce_field_equation_inplace(fq_nmod_mpoly_t poly, const fq_nmod_mpoly_ctx_t ctx) {
    if (!g_field_equation_reduction) return;
    slong nterms = fq_nmod_mpoly_length(poly, ctx);
    if (nterms <= 0) return;
    ulong q = field_size_q_from_fq_ctx(ctx->fqctx);
    if (q <= 1 || q == WORD_MAX) return;

    slong nvars = fq_nmod_mpoly_ctx_nvars(ctx);
    ulong *exp = (ulong*) flint_malloc(nvars * sizeof(ulong));
    fq_nmod_t coeff;
    fq_nmod_init(coeff, ctx->fqctx);
    fq_nmod_mpoly_t reduced;
    fq_nmod_mpoly_init(reduced, ctx);

    for (slong i = 0; i < nterms; i++) {
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff, poly, i, ctx);
        fq_nmod_mpoly_get_term_exp_ui(exp, poly, i, ctx);
        for (slong k = 0; k < nvars; k++) exp[k] = reduce_exp_field_ui(exp[k], q);
        fq_nmod_mpoly_push_term_fq_nmod_ui(reduced, coeff, exp, ctx);
    }

    fq_nmod_mpoly_sort_terms(reduced, ctx);
    fq_nmod_mpoly_combine_like_terms(reduced, ctx);
    fq_nmod_mpoly_set(poly, reduced, ctx);
    fq_nmod_mpoly_clear(reduced, ctx);
    fq_nmod_clear(coeff, ctx->fqctx);
    flint_free(exp);
}

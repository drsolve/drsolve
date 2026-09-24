/* SPDX-License-Identifier: GPL-2.0-or-later */
/*
 * dixon_flint.c - Dixon Resultant Implementation for Finite Extension Fields
 *
 * This file contains the complete implementation of Dixon resultant computation
 * over finite fields using FLINT library for polynomial arithmetic and matrix operations.
 */

#include "dixon_flint.h"
#include "mq_poly_mat_det.h"
#include "mq_simplex_det.h"
#include "mq_pencil_det.h"

/* Internal row-basis state used only by the Dixon implementation. */
typedef struct {
    field_elem_u *reduced_rows;
    slong *pivot_cols;
    slong *selected_indices;
    slong current_rank;
    slong max_size;
    slong ncols;
    field_ctx_t *ctx;
    int initialized;
    field_elem_u *work_row;
    field_elem_u *temp_vars;
    int workspace_initialized;
} unified_row_basis_tracker_t;

static void find_pivot_rows_nmod_fixed(slong **selected_rows_out,
                                       slong *num_selected,
                                       const nmod_mat_t mat);
static void find_pivot_rows_simple(slong **selected_rows_out,
                                   slong *num_selected,
                                   const field_elem_u *unified_mat,
                                   slong nrows, slong ncols,
                                   field_ctx_t *ctx);

// Global method selection variable definitions
det_method_t dixon_global_method_step1 = -1;
det_method_t dixon_global_method_step4 = -1;
det_method_t dixon_global_method = -1; // deprecated compatibility alias
resultant_method_t g_resultant_method = RESULTANT_METHOD_DIXON;
int g_dixon_verbose_level = 1;
int g_dixon_debug_mode = 0;
int g_dixon_show_step_timing = 0;
int g_matrix_transpose_threshold = 10000;
rational_root_scan_mode_t g_rational_root_scan_mode = RATIONAL_ROOT_SCAN_AUTO;
int g_dixon_fast_use_ksy_precondition = 0;
slong g_dixon_fast_ksy_constant_col = 0;
int g_dixon_mq_step1_filter = 1;
int g_dixon_mq_step1_shared = 1;
int g_dixon_mq_step1_rank = 1;
int g_dixon_mq_step1_simplex = 0;
int g_dixon_mq_step1_pencil = 0;
int g_dixon_mq_step4_schur = 1;
int g_dixon_step3_second_verification = 0;
slong g_dixon_det_cache_limit = 1024;

static const char *dixon_det_method_name(det_method_t method)
{
    switch (method) {
        case DET_METHOD_RECURSIVE:
            return "minor expansion";
        case DET_METHOD_KRONECKER:
            return "fast HNF";
        case DET_METHOD_INTERPOLATION:
            return "interpolation";
        case DET_METHOD_HUANG:
            return "sparse interpolation";
        case DET_METHOD_KRONECKER_NMOD:
            return "Bareiss";
        case DET_METHOD_BALANCED_SPLIT:
            return "balanced split Laplace (experimental)";
        default:
            return "default";
    }
}

static int dixon_estimate_interpolation_points(slong *total_points_out,
                                               fq_mvpoly_t **matrix,
                                               slong size)
{
    slong actual_nvars;
    slong actual_npars;
    slong total_vars;
    slong *degree_bounds;
    ulong prime;
    slong field_degree;
    slong field_size = 1;
    int field_size_capped = 0;
    slong total_points = 1;

    if (!total_points_out || !matrix || size <= 0) return 0;

    actual_nvars = matrix[0][0].nvars;
    actual_npars = matrix[0][0].npars;
    total_vars = actual_nvars + actual_npars;
    if (total_vars <= 0) {
        *total_points_out = 1;
        return 1;
    }

    degree_bounds = (slong *) flint_malloc((size_t) total_vars * sizeof(slong));
    if (!degree_bounds) return 0;

    fq_compute_det_degree_bounds_optimized(degree_bounds, matrix, size, total_vars);

    prime = fq_nmod_ctx_prime(matrix[0][0].ctx);
    field_degree = fq_nmod_ctx_degree(matrix[0][0].ctx);
    for (slong i = 0; i < field_degree; i++) {
        if (prime != 0 && field_size > LONG_MAX / (slong) prime) {
            field_size = LONG_MAX;
            field_size_capped = 1;
            break;
        }
        field_size *= (slong) prime;
    }

    for (slong i = 0; i < total_vars; i++) {
        slong grid_size = degree_bounds[i] + 1;
        if (grid_size < 1) grid_size = 1;
        if (!field_size_capped && grid_size > field_size) {
            grid_size = field_size;
        }
        if (grid_size < 1) grid_size = 1;

        if (total_points > DIXON_INTERPOLATION_POINT_LIMIT / grid_size) {
            flint_free(degree_bounds);
            *total_points_out = DIXON_INTERPOLATION_POINT_LIMIT + 1;
            return 1;
        }
        total_points *= grid_size;
    }

    flint_free(degree_bounds);
    *total_points_out = total_points;
    return 1;
}

static void dixon_info_log(const char *fmt, ...)
{
    va_list args;

    if (g_dixon_verbose_level < 1) {
        return;
    }

    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}

static det_method_t dixon_resolve_step1_det_method(fq_mvpoly_t **modified_M_mvpoly,
                                                   slong nvars,
                                                   slong npars,
                                                   det_method_t requested_method,
                                                   int emit_warning)
{
    slong total_points = 0;

    if (requested_method != DET_METHOD_INTERPOLATION) return requested_method;
    if (nvars < 0 || npars < 0) return requested_method;

    if (!dixon_estimate_interpolation_points(&total_points, modified_M_mvpoly, nvars + 1) ||
        total_points > DIXON_INTERPOLATION_POINT_LIMIT) {
        if (emit_warning) {
            if (total_points > DIXON_INTERPOLATION_POINT_LIMIT) {
                dixon_info_log("  Warning: interpolation for Step 1 needs over %ld points; "
                               "falling back to minor expansion\n",
                               total_points);
            } else {
                dixon_info_log("  Warning: failed to estimate Step 1 interpolation size; "
                               "falling back to minor expansion\n");
            }
        }
        return DET_METHOD_RECURSIVE;
    }

    return requested_method;
}

int dixon_method_uses_parallel_timing(det_method_t method)
{
    return method == DET_METHOD_INTERPOLATION;
}

int dixon_get_effective_interpolation_threads(void)
{
    int threads = 1;

#ifdef _OPENMP
    if (g_interpolation_threads > 0) {
        return g_interpolation_threads;
    }
    threads = omp_get_max_threads();
    threads = (threads + 1) / 2;
#endif
    return threads > 0 ? threads : 1;
}

static int dixon_get_effective_parallel_threads(void)
{
    int threads = 1;

#ifdef _OPENMP
    threads = omp_get_max_threads();
#endif
    return threads > 0 ? threads : 1;
}


void dixon_maybe_print_step_time(const char *step_label, double wall_elapsed)
{
    if (g_dixon_verbose_level < 1 || !g_dixon_show_step_timing) {
        return;
    }
    dixon_info_log("%s time: %.3f seconds\n", step_label, wall_elapsed);
}

void dixon_maybe_print_parallel_step_time(const char *step_label,
                                          double cpu_elapsed,
                                          double wall_elapsed)
{
    if (g_dixon_verbose_level < 1 || !g_dixon_show_step_timing) {
        return;
    }

    dixon_info_log("%s time: CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n",
                   step_label, cpu_elapsed, wall_elapsed,
                   dixon_get_effective_parallel_threads());
}

static void dixon_maybe_print_step_detail_time(const char *label,
                                               clock_t cpu_start,
                                               double wall_start)
{
    if (g_dixon_verbose_level < 3 || !g_dixon_show_step_timing) {
        return;
    }

    dixon_info_log("  %s time: CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n",
                   label,
                   (double) (clock() - cpu_start) / CLOCKS_PER_SEC,
                   get_wall_time() - wall_start,
                   dixon_get_effective_parallel_threads());
}

void dixon_maybe_print_step_method_time(const char *step_label,
                                        det_method_t method,
                                        double cpu_elapsed,
                                        double wall_elapsed)
{
    int show_parallel_style = 0;

    if (g_dixon_verbose_level < 1 || !g_dixon_show_step_timing) {
        return;
    }

    show_parallel_style = dixon_method_uses_parallel_timing(method) ||
                          (strcmp(step_label, "Step 1") == 0 &&
                           method == DET_METHOD_RECURSIVE);

    if (show_parallel_style) {
        dixon_info_log("%s time: CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n",
                       step_label, cpu_elapsed, wall_elapsed,
                       dixon_get_effective_interpolation_threads());
    } else {
        dixon_info_log("%s time: %.3f seconds\n", step_label, wall_elapsed);
    }
}

static void dixon_debug_log(const char *fmt, ...)
{
    va_list args;

    if (!g_dixon_debug_mode) {
        return;
    }

    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
}

static ulong hash_monom_exponents(const slong *exp, slong nvars)
{
    ulong hash = 1469598103934665603UL;

    for (slong i = 0; i < nvars; i++) {
        hash ^= (ulong) exp[i] + 0x9e3779b97f4a7c15UL + (hash << 6) + (hash >> 2);
    }

    return hash;
}

static void free_monom_index(hash_entry_t **buckets, slong hash_size)
{
    if (buckets == NULL) {
        return;
    }

    for (slong i = 0; i < hash_size; i++) {
        hash_entry_t *entry = buckets[i];
        while (entry != NULL) {
            hash_entry_t *next = entry->next;
            flint_free(entry);
            entry = next;
        }
    }

    flint_free(buckets);
}

static slong lookup_monom_index(hash_entry_t **buckets,
                                slong hash_size,
                                const slong *exp,
                                slong nvars)
{
    if (buckets == NULL || exp == NULL || hash_size <= 0) {
        return -1;
    }

    ulong hash = hash_monom_exponents(exp, nvars) & (ulong) (hash_size - 1);
    hash_entry_t *entry = buckets[hash];

    while (entry != NULL) {
        if (memcmp(entry->exp, exp, (size_t) nvars * sizeof(slong)) == 0) {
            return entry->idx;
        }
        entry = entry->next;
    }

    return -1;
}

static int dixon_show_step_details(void)
{
    return g_dixon_verbose_level >= 2;
}

static int dixon_should_dump_small_matrix(slong nrows, slong ncols)
{
    return (g_dixon_verbose_level == 3 && nrows <= 3 && ncols <= 3) || (g_dixon_verbose_level >= 4 && nrows <= 100 && ncols <= 100);
}

static void dixon_print_basis_monomial(const monom_t *monom,
                                       slong nvars,
                                       char **var_names,
                                       int dual)
{
    static const char default_var_names[] = { 'x', 'y', 'z', 'w', 'v', 'u' };
    int printed = 0;

    if (!monom || !monom->exp) {
        printf("0");
        return;
    }

    for (slong i = 0; i < nvars; i++) {
        slong exp = monom->exp[i];
        if (exp <= 0) {
            continue;
        }

        if (printed) {
            printf("*");
        }

        if (dual) {
            printf("~");
        }

        if (var_names && var_names[i]) {
            printf("%s", var_names[i]);
        } else if (i < (slong)(sizeof(default_var_names) / sizeof(default_var_names[0]))) {
            printf("%c", default_var_names[i]);
        } else {
            printf("x_%ld", i);
        }

        if (exp > 1) {
            printf("^%ld", exp);
        }

        printed = 1;
    }

    if (!printed) {
        printf("1");
    }
}

static void dixon_print_small_named_dense_matrix(const char *title,
                                                 fq_mvpoly_t **matrix,
                                                 slong nrows,
                                                 slong ncols,
                                                 char **var_names,
                                                 char **par_names,
                                                 const char *gen_name,
                                                 int expanded_format)
{
    if (!dixon_should_dump_small_matrix(nrows, ncols)) {
        return;
    }

    dixon_info_log("  %s details:\n", title);
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            printf("      M[%ld][%ld] = ", i, j);
            fq_mvpoly_print_with_names(&matrix[i][j], "", var_names, par_names,
                                       gen_name, expanded_format);
        }
    }
}

static void dixon_print_small_sparse_matrix(const char *title,
                                            fq_mvpoly_t ***matrix,
                                            slong nrows,
                                            slong ncols,
                                            const monom_t *row_monoms,
                                            const monom_t *col_monoms,
                                            slong nvars,
                                            char **var_names,
                                            char **par_names,
                                            const char *gen_name)
{
    if (!dixon_should_dump_small_matrix(nrows, ncols)) {
        return;
    }

    dixon_info_log("  %s details:\n", title);
    dixon_info_log("    Row basis:\n");
    for (slong i = 0; i < nrows; i++) {
        printf("      r%ld = ", i);
        dixon_print_basis_monomial(&row_monoms[i], nvars, var_names, 0);
        printf("\n");
    }

    dixon_info_log("    Column basis:\n");
    for (slong j = 0; j < ncols; j++) {
        printf("      c%ld = ", j);
        dixon_print_basis_monomial(&col_monoms[j], nvars, var_names, 1);
        printf("\n");
    }

    dixon_info_log("    Entries:\n");
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            printf("      D[%ld][%ld] = ", i, j);
            if (matrix[i][j] != NULL) {
                fq_mvpoly_print_with_names(matrix[i][j], "", var_names, par_names,
                                           gen_name, 1);
            } else {
                printf("0\n");
            }
        }
    }
}

static void dixon_print_small_dense_submatrix(const char *title,
                                              fq_mvpoly_t **matrix,
                                              slong nrows,
                                              slong ncols,
                                              const slong *row_indices,
                                              const slong *col_indices,
                                              const monom_t *row_monoms,
                                              const monom_t *col_monoms,
                                              slong nvars,
                                              char **var_names,
                                              char **par_names,
                                              const char *gen_name)
{
    if (!dixon_should_dump_small_matrix(nrows, ncols)) {
        return;
    }

    dixon_info_log("  %s details:\n", title);
    dixon_info_log("    Selected row basis:\n");
    for (slong i = 0; i < nrows; i++) {
        printf("      r%ld <- row %ld = ", i, row_indices[i]);
        dixon_print_basis_monomial(&row_monoms[row_indices[i]], nvars, var_names, 0);
        printf("\n");
    }

    dixon_info_log("    Selected column basis:\n");
    for (slong j = 0; j < ncols; j++) {
        printf("      c%ld <- col %ld = ", j, col_indices[j]);
        dixon_print_basis_monomial(&col_monoms[col_indices[j]], nvars, var_names, 1);
        printf("\n");
    }

    dixon_info_log("    Entries:\n");
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            printf("      S[%ld][%ld] = ", i, j);
            fq_mvpoly_print_with_names(&matrix[i][j], "", var_names, par_names,
                                       gen_name, 1);
        }
    }
}

static void init_evaluation_parameters(fq_nmod_t *param_vals, slong npars,
                                      const fq_nmod_ctx_t ctx,
                                      slong attempt)
{
    for (slong i = 0; i < npars; i++) {
        fq_nmod_init(param_vals[i], ctx);
        fq_nmod_set_si(param_vals[i], 7 * (attempt + 1) * (i + 1) + 13, ctx);
    }
}

static void clear_evaluation_parameters(fq_nmod_t *param_vals, slong npars,
                                       const fq_nmod_ctx_t ctx)
{
    for (slong i = 0; i < npars; i++) {
        fq_nmod_clear(param_vals[i], ctx);
    }
    flint_free(param_vals);
}

static void init_extension_evaluation_parameters(fq_nmod_t *param_vals,
                                                slong npars,
                                                const fq_nmod_ctx_t ext_ctx,
                                                slong attempt)
{
    fq_nmod_t gen, constant;
    fq_nmod_init(gen, ext_ctx);
    fq_nmod_init(constant, ext_ctx);
    fq_nmod_gen(gen, ext_ctx);

    for (slong i = 0; i < npars; i++) {
        fq_nmod_init(param_vals[i], ext_ctx);
        fq_nmod_set(param_vals[i], gen, ext_ctx);
        fq_nmod_set_ui(constant, 7 * (attempt + 1) * (i + 1) + 13, ext_ctx);
        fq_nmod_add(param_vals[i], param_vals[i], constant, ext_ctx);
        if (fq_nmod_is_zero(param_vals[i], ext_ctx)) {
            fq_nmod_set(param_vals[i], gen, ext_ctx);
        }
    }

    fq_nmod_clear(gen, ext_ctx);
    fq_nmod_clear(constant, ext_ctx);
}

static void evaluate_fq_mvpoly_at_extension_params(fq_nmod_t result,
                                                  const fq_mvpoly_t *poly,
                                                  const fq_nmod_t *param_vals,
                                                  const fq_nmod_ctx_t ext_ctx)
{
    fq_nmod_t term_val, coeff_val;
    fq_nmod_init(term_val, ext_ctx);
    fq_nmod_init(coeff_val, ext_ctx);
    fq_nmod_zero(result, ext_ctx);

    if (poly == NULL || poly->nterms == 0) {
        fq_nmod_clear(term_val, ext_ctx);
        fq_nmod_clear(coeff_val, ext_ctx);
        return;
    }

    for (slong i = 0; i < poly->nterms; i++) {
        fq_nmod_set_ui(coeff_val, nmod_poly_get_coeff_ui(poly->terms[i].coeff, 0), ext_ctx);
        fq_nmod_set(term_val, coeff_val, ext_ctx);

        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                slong exp = poly->terms[i].par_exp[j];
                if (exp > 0) {
                    fq_nmod_t pow_val;
                    fq_nmod_init(pow_val, ext_ctx);
                    fq_nmod_pow_ui(pow_val, param_vals[j], exp, ext_ctx);
                    fq_nmod_mul(term_val, term_val, pow_val, ext_ctx);
                    fq_nmod_clear(pow_val, ext_ctx);
                }
            }
        }

        fq_nmod_add(result, result, term_val, ext_ctx);
    }

    fq_nmod_clear(term_val, ext_ctx);
    fq_nmod_clear(coeff_val, ext_ctx);
}

static slong evaluate_selected_submatrix_rank(fq_mvpoly_t ***full_matrix,
                                             const slong *row_indices,
                                             const slong *col_indices,
                                             slong size,
                                             const fq_nmod_t *param_vals,
                                             const fq_nmod_ctx_t ctx)
{
    fq_nmod_mat_t mat;
    fq_nmod_t value;
    fq_nmod_mat_init(mat, size, size, ctx);
#ifdef _OPENMP
    #pragma omp parallel
    {
        fq_nmod_t local_value;
        fq_nmod_init(local_value, ctx);
        #pragma omp for schedule(static)
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
                if (entry != NULL) {
                    evaluate_fq_mvpoly_at_params(local_value, entry, param_vals);
                    fq_nmod_set(fq_nmod_mat_entry(mat, i, j), local_value, ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ctx);
                }
            }
        }
        fq_nmod_clear(local_value, ctx);
    }
#else
    fq_nmod_init(value, ctx);
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
            if (entry != NULL) {
                evaluate_fq_mvpoly_at_params(value, entry, param_vals);
                fq_nmod_set(fq_nmod_mat_entry(mat, i, j), value, ctx);
            } else {
                fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ctx);
            }
        }
    }
    fq_nmod_clear(value, ctx);
#endif

    slong rank = fq_nmod_mat_rank(mat, ctx);
    fq_nmod_mat_clear(mat, ctx);
    return rank;
}

static slong evaluate_selected_rectangular_submatrix_rank(fq_mvpoly_t ***full_matrix,
                                                         const slong *row_indices,
                                                         slong num_rows,
                                                         const slong *col_indices,
                                                         slong num_cols,
                                                         const fq_nmod_t *param_vals,
                                                         const fq_nmod_ctx_t ctx)
{
    fq_nmod_mat_t mat;
    fq_nmod_t value;

    fq_nmod_mat_init(mat, num_rows, num_cols, ctx);
#ifdef _OPENMP
    #pragma omp parallel
    {
        fq_nmod_t local_value;
        fq_nmod_init(local_value, ctx);
        #pragma omp for schedule(static)
        for (slong i = 0; i < num_rows; i++) {
            for (slong j = 0; j < num_cols; j++) {
                fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
                if (entry != NULL) {
                    evaluate_fq_mvpoly_at_params(local_value, entry, param_vals);
                    fq_nmod_set(fq_nmod_mat_entry(mat, i, j), local_value, ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ctx);
                }
            }
        }
        fq_nmod_clear(local_value, ctx);
    }
#else
    fq_nmod_init(value, ctx);
    for (slong i = 0; i < num_rows; i++) {
        for (slong j = 0; j < num_cols; j++) {
            fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
            if (entry != NULL) {
                evaluate_fq_mvpoly_at_params(value, entry, param_vals);
                fq_nmod_set(fq_nmod_mat_entry(mat, i, j), value, ctx);
            } else {
                fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ctx);
            }
        }
    }
    fq_nmod_clear(value, ctx);
#endif

    slong rank = fq_nmod_mat_rank(mat, ctx);
    fq_nmod_mat_clear(mat, ctx);
    return rank;
}

static slong evaluate_selected_submatrix_rank_extension(fq_mvpoly_t ***full_matrix,
                                                       const slong *row_indices,
                                                       const slong *col_indices,
                                                       slong size,
                                                       const fq_nmod_t *param_vals,
                                                       const fq_nmod_ctx_t ext_ctx)
{
    fq_nmod_mat_t mat;
    fq_nmod_t value;
    fq_nmod_mat_init(mat, size, size, ext_ctx);
#ifdef _OPENMP
    #pragma omp parallel
    {
        fq_nmod_t local_value;
        fq_nmod_init(local_value, ext_ctx);
        #pragma omp for schedule(static)
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
                if (entry != NULL) {
                    evaluate_fq_mvpoly_at_extension_params(local_value, entry, param_vals, ext_ctx);
                    fq_nmod_set(fq_nmod_mat_entry(mat, i, j), local_value, ext_ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ext_ctx);
                }
            }
        }
        fq_nmod_clear(local_value, ext_ctx);
    }
#else
    fq_nmod_init(value, ext_ctx);
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
            if (entry != NULL) {
                evaluate_fq_mvpoly_at_extension_params(value, entry, param_vals, ext_ctx);
                fq_nmod_set(fq_nmod_mat_entry(mat, i, j), value, ext_ctx);
            } else {
                fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ext_ctx);
            }
        }
    }
    fq_nmod_clear(value, ext_ctx);
#endif

    slong rank = fq_nmod_mat_rank(mat, ext_ctx);
    fq_nmod_mat_clear(mat, ext_ctx);
    return rank;
}

static slong evaluate_selected_rectangular_submatrix_rank_extension(
    fq_mvpoly_t ***full_matrix,
    const slong *row_indices,
    slong num_rows,
    const slong *col_indices,
    slong num_cols,
    const fq_nmod_t *param_vals,
    const fq_nmod_ctx_t ext_ctx)
{
    fq_nmod_mat_t mat;
    fq_nmod_t value;

    fq_nmod_mat_init(mat, num_rows, num_cols, ext_ctx);
#ifdef _OPENMP
    #pragma omp parallel
    {
        fq_nmod_t local_value;
        fq_nmod_init(local_value, ext_ctx);
        #pragma omp for schedule(static)
        for (slong i = 0; i < num_rows; i++) {
            for (slong j = 0; j < num_cols; j++) {
                fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
                if (entry != NULL) {
                    evaluate_fq_mvpoly_at_extension_params(local_value, entry, param_vals, ext_ctx);
                    fq_nmod_set(fq_nmod_mat_entry(mat, i, j), local_value, ext_ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ext_ctx);
                }
            }
        }
        fq_nmod_clear(local_value, ext_ctx);
    }
#else
    fq_nmod_init(value, ext_ctx);
    for (slong i = 0; i < num_rows; i++) {
        for (slong j = 0; j < num_cols; j++) {
            fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
            if (entry != NULL) {
                evaluate_fq_mvpoly_at_extension_params(value, entry, param_vals, ext_ctx);
                fq_nmod_set(fq_nmod_mat_entry(mat, i, j), value, ext_ctx);
            } else {
                fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ext_ctx);
            }
        }
    }
    fq_nmod_clear(value, ext_ctx);
#endif

    slong rank = fq_nmod_mat_rank(mat, ext_ctx);
    fq_nmod_mat_clear(mat, ext_ctx);
    return rank;
}

static void evaluate_fq_mvpoly_at_nmod_params_direct(fq_nmod_t result,
                                                     const fq_mvpoly_t *poly,
                                                     const fq_nmod_t *param_vals,
                                                     const fq_nmod_ctx_t ctx)
{
    fq_nmod_t term_val, coeff_val, pow_val;

    fq_nmod_init(term_val, ctx);
    fq_nmod_init(coeff_val, ctx);
    fq_nmod_init(pow_val, ctx);
    fq_nmod_zero(result, ctx);

    if (poly == NULL || poly->nterms == 0) {
        fq_nmod_clear(term_val, ctx);
        fq_nmod_clear(coeff_val, ctx);
        fq_nmod_clear(pow_val, ctx);
        return;
    }

    for (slong i = 0; i < poly->nterms; i++) {
        fq_nmod_set_ui(coeff_val, nmod_poly_get_coeff_ui(poly->terms[i].coeff, 0), ctx);
        fq_nmod_set(term_val, coeff_val, ctx);

        if (poly->terms[i].par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                slong exp = poly->terms[i].par_exp[j];
                if (exp > 0) {
                    fq_nmod_pow_ui(pow_val, param_vals[j], exp, ctx);
                    fq_nmod_mul(term_val, term_val, pow_val, ctx);
                }
            }
        }

        fq_nmod_add(result, result, term_val, ctx);
    }

    fq_nmod_clear(term_val, ctx);
    fq_nmod_clear(coeff_val, ctx);
    fq_nmod_clear(pow_val, ctx);
}

static slong evaluate_selected_submatrix_rank_nmod(fq_mvpoly_t ***full_matrix,
                                                   const slong *row_indices,
                                                   const slong *col_indices,
                                                   slong size,
                                                   const fq_nmod_t *param_vals,
                                                   const fq_nmod_ctx_t ctx)
{
    fq_nmod_mat_t mat;
    fq_nmod_t value;

    fq_nmod_mat_init(mat, size, size, ctx);
    fq_nmod_init(value, ctx);

    #ifdef _OPENMP
    #pragma omp parallel
    {
        fq_nmod_t local_value;
        fq_nmod_init(local_value, ctx);
        #pragma omp for schedule(static)
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
                if (entry != NULL) {
                    evaluate_fq_mvpoly_at_nmod_params_direct(local_value, entry, param_vals, ctx);
                    fq_nmod_set(fq_nmod_mat_entry(mat, i, j), local_value, ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ctx);
                }
            }
        }
        fq_nmod_clear(local_value, ctx);
    }
    #else
    for (slong i = 0; i < size; i++) {
        for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
            if (entry != NULL) {
                evaluate_fq_mvpoly_at_nmod_params_direct(value, entry, param_vals, ctx);
                fq_nmod_set(fq_nmod_mat_entry(mat, i, j), value, ctx);
            } else {
                fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ctx);
            }
        }
    }
    #endif

    slong rank = fq_nmod_mat_rank(mat, ctx);
    fq_nmod_clear(value, ctx);
    fq_nmod_mat_clear(mat, ctx);
    return rank;
}

static slong evaluate_selected_rectangular_submatrix_rank_nmod(
    fq_mvpoly_t ***full_matrix,
    const slong *row_indices,
    slong num_rows,
    const slong *col_indices,
    slong num_cols,
    const fq_nmod_t *param_vals,
    const fq_nmod_ctx_t ctx)
{
    fq_nmod_mat_t mat;
    fq_nmod_t value;

    fq_nmod_mat_init(mat, num_rows, num_cols, ctx);
    fq_nmod_init(value, ctx);

    #ifdef _OPENMP
    #pragma omp parallel
    {
        fq_nmod_t local_value;
        fq_nmod_init(local_value, ctx);
        #pragma omp for schedule(static)
        for (slong i = 0; i < num_rows; i++) {
            for (slong j = 0; j < num_cols; j++) {
                fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
                if (entry != NULL) {
                    evaluate_fq_mvpoly_at_nmod_params_direct(local_value, entry, param_vals, ctx);
                    fq_nmod_set(fq_nmod_mat_entry(mat, i, j), local_value, ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ctx);
                }
            }
        }
        fq_nmod_clear(local_value, ctx);
    }
    #else
    for (slong i = 0; i < num_rows; i++) {
        for (slong j = 0; j < num_cols; j++) {
            fq_mvpoly_t *entry = full_matrix[row_indices[i]][col_indices[j]];
            if (entry != NULL) {
                evaluate_fq_mvpoly_at_nmod_params_direct(value, entry, param_vals, ctx);
                fq_nmod_set(fq_nmod_mat_entry(mat, i, j), value, ctx);
            } else {
                fq_nmod_zero(fq_nmod_mat_entry(mat, i, j), ctx);
            }
        }
    }
    #endif

    slong rank = fq_nmod_mat_rank(mat, ctx);
    fq_nmod_clear(value, ctx);
    fq_nmod_mat_clear(mat, ctx);
    return rank;
}

static int selected_rows_satisfy_ksy_precondition_nmod(
    fq_mvpoly_t ***full_matrix,
    const slong *row_indices,
    slong size,
    slong total_cols,
    slong ksy_constant_col,
    const fq_nmod_t *param_vals,
    const fq_nmod_ctx_t ctx)
{
    slong *reduced_cols;
    slong reduced_rank;

    if (ksy_constant_col < 0 || size <= 0 || total_cols <= 0) {
        return 1;
    }
    if (ksy_constant_col >= total_cols) {
        return 0;
    }
    if (size == 1) {
        return 1;
    }

    reduced_cols = (slong *) flint_malloc((size_t) (total_cols - 1) * sizeof(slong));
    for (slong j = 0, out = 0; j < total_cols; j++) {
        if (j == ksy_constant_col) {
            continue;
        }
        reduced_cols[out++] = j;
    }

    reduced_rank = evaluate_selected_rectangular_submatrix_rank_nmod(full_matrix,
                                                                     row_indices,
                                                                     size,
                                                                     reduced_cols,
                                                                     total_cols - 1,
                                                                     param_vals,
                                                                     ctx);
    flint_free(reduced_cols);
    return reduced_rank == size - 1;
}

static int selected_rows_satisfy_ksy_precondition(
    fq_mvpoly_t ***full_matrix,
    const slong *row_indices,
    slong size,
    slong total_cols,
    slong ksy_constant_col,
    const fq_nmod_t *param_vals,
    const fq_nmod_ctx_t ctx)
{
    slong *reduced_cols;
    slong reduced_rank;

    if (ksy_constant_col < 0 || size <= 0 || total_cols <= 0) {
        return 1;
    }

    if (ksy_constant_col >= total_cols) {
        return 0;
    }

    if (size == 1) {
        return 1;
    }

    reduced_cols = (slong *) flint_malloc((size_t) (total_cols - 1) * sizeof(slong));
    for (slong j = 0, out = 0; j < total_cols; j++) {
        if (j == ksy_constant_col) {
            continue;
        }
        reduced_cols[out++] = j;
    }

    reduced_rank = evaluate_selected_rectangular_submatrix_rank(full_matrix,
                                                                row_indices,
                                                                size,
                                                                reduced_cols,
                                                                total_cols - 1,
                                                                param_vals,
                                                                ctx);
    flint_free(reduced_cols);
    return reduced_rank == size - 1;
}

static int selected_rows_satisfy_ksy_precondition_extension(
    fq_mvpoly_t ***full_matrix,
    const slong *row_indices,
    slong size,
    slong total_cols,
    slong ksy_constant_col,
    const fq_nmod_t *param_vals,
    const fq_nmod_ctx_t ext_ctx)
{
    slong *reduced_cols;
    slong reduced_rank;

    if (ksy_constant_col < 0 || size <= 0 || total_cols <= 0) {
        return 1;
    }

    if (ksy_constant_col >= total_cols) {
        return 0;
    }

    if (size == 1) {
        return 1;
    }

    reduced_cols = (slong *) flint_malloc((size_t) (total_cols - 1) * sizeof(slong));
    for (slong j = 0, out = 0; j < total_cols; j++) {
        if (j == ksy_constant_col) {
            continue;
        }
        reduced_cols[out++] = j;
    }

    reduced_rank = evaluate_selected_rectangular_submatrix_rank_extension(full_matrix,
                                                                          row_indices,
                                                                          size,
                                                                          reduced_cols,
                                                                          total_cols - 1,
                                                                          param_vals,
                                                                          ext_ctx);
    flint_free(reduced_cols);
    return reduced_rank == size - 1;
}

// Build cancellation matrix in multivariate form
void build_fq_cancellation_matrix_mvpoly(fq_mvpoly_t ***M, fq_mvpoly_t *polys, 
                                        slong nvars, slong npars) {
    slong n = nvars + 1;
    
    // Allocate matrix space
    *M = (fq_mvpoly_t**) flint_malloc(n * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < n; i++) {
        (*M)[i] = (fq_mvpoly_t*) flint_malloc(n * sizeof(fq_mvpoly_t));
    }
    
    // Build matrix entries
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            fq_mvpoly_init(&(*M)[i][j], 2 * nvars, npars, polys[0].ctx);
            
            // Substitute variables according to row index
            for (slong t = 0; t < polys[j].nterms; t++) {
                slong *new_var_exp = (slong*) flint_calloc(2 * nvars, sizeof(slong));
                
                for (slong k = 0; k < nvars; k++) {
                    slong orig_exp = polys[j].terms[t].var_exp ? polys[j].terms[t].var_exp[k] : 0;
                    
                    if (k < i) {
                        // Use dual variable ~x_k
                        new_var_exp[nvars + k] = orig_exp;
                    } else {
                        // Use original variable x_k
                        new_var_exp[k] = orig_exp;
                    }
                }
                
                fq_mvpoly_add_term_fast(&(*M)[i][j], new_var_exp, polys[j].terms[t].par_exp, 
                                  polys[j].terms[t].coeff);
                flint_free(new_var_exp);
            }
        }
    }
}

void perform_fq_matrix_row_operations_mvpoly(fq_mvpoly_t ***new_matrix, fq_mvpoly_t ***original_matrix,
                                                   slong nvars, slong npars) {
    slong n = nvars + 1;
    
    *new_matrix = (fq_mvpoly_t**) flint_malloc(n * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < n; i++) {
        (*new_matrix)[i] = (fq_mvpoly_t*) flint_malloc(n * sizeof(fq_mvpoly_t));
    }
    
    for (slong j = 0; j < n; j++) {
        fq_mvpoly_copy(&(*new_matrix)[0][j], &(*original_matrix)[0][j]);
    }
    
    for (slong i = 1; i < n; i++) {
        for (slong j = 0; j < n; j++) {
            fq_mvpoly_t diff = {0};
            if (g_dixon_verbose_level > 3) {
                dixon_debug_log("    RowOp[%ld,%ld]: top terms=%ld, bottom terms=%ld\n",
                                i, j,
                                (*original_matrix)[i][j].nterms,
                                (*original_matrix)[i-1][j].nterms);
            }
            fq_mvpoly_sub(&diff, &(*original_matrix)[i][j], &(*original_matrix)[i-1][j]);
            if (g_dixon_verbose_level > 3) {
                dixon_debug_log("    RowOp[%ld,%ld]: diff terms=%ld\n", i, j, diff.nterms);
            }
      
            if (diff.nterms > 0) {
                divide_by_fq_linear_factor_flint(&(*new_matrix)[i][j], &diff, 
                                               i-1, 2*nvars, npars);
                if (g_dixon_verbose_level > 3) {
                    dixon_debug_log("    RowOp[%ld,%ld]: quotient terms=%ld\n",
                                    i, j, (*new_matrix)[i][j].nterms);
                }
            } else {
                fq_mvpoly_init(&(*new_matrix)[i][j], 2*nvars, npars, diff.ctx);
                if (g_dixon_verbose_level > 3) {
                    dixon_debug_log("    RowOp[%ld,%ld]: zero diff, initialized zero quotient\n",
                                    i, j);
                }
            }
             
            fq_mvpoly_clear(&diff);
        }
    }
}

// Compute Dixon resultant degree bound
slong compute_fq_dixon_resultant_degree_bound(fq_mvpoly_t *polys, slong npolys, slong nvars, slong npars) {
    slong degree_product = 1;
    
    for (slong i = 0; i < npolys; i++) {
        slong max_total_deg = 0;
        
        // Find maximum total degree of polynomial i
        for (slong t = 0; t < polys[i].nterms; t++) {
            slong total_deg = 0;
            
            // Sum variable degrees
            if (polys[i].terms[t].var_exp) {
                for (slong j = 0; j < nvars; j++) {
                    total_deg += polys[i].terms[t].var_exp[j];
                }
            }
            
            // Sum parameter degrees
            if (polys[i].terms[t].par_exp && npars > 0) {
                for (slong j = 0; j < npars; j++) {
                    total_deg += polys[i].terms[t].par_exp[j];
                }
            }
            
            if (total_deg > max_total_deg) {
                max_total_deg = total_deg;
            }
        }
        
        degree_product *= max_total_deg;
    }
    
    return degree_product + 1;
}

void compute_fq_coefficient_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **coeff_matrix,
                                       slong size, slong npars, const fq_nmod_ctx_t ctx,
                                       det_method_t method, slong res_deg_bound) {
    if (size == 0) {
        fq_mvpoly_init(result, 0, npars, ctx);
        return;
    }
    
    if (npars == 0) {
        fq_mvpoly_init(result, 0, npars, ctx);
        
        fq_nmod_mat_t scalar_mat;
        fq_nmod_mat_init(scalar_mat, size, size, ctx);
        
        for (slong i = 0; i < size; i++) {
            for (slong j = 0; j < size; j++) {
                if (coeff_matrix[i][j].nterms > 0) {
                    fq_nmod_set(fq_nmod_mat_entry(scalar_mat, i, j), 
                                coeff_matrix[i][j].terms[0].coeff, ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(scalar_mat, i, j), ctx);
                }
            }
        }
        
        dixon_debug_log("  Computing Resultant\n");
        clock_t cpu_start = clock();
        double wall_start = get_wall_time();
        
        fq_nmod_t det;
        fq_nmod_init(det, ctx);
        fq_nmod_mat_det(det, scalar_mat, ctx);
        
        clock_t cpu_end = clock();
        double wall_end = get_wall_time();
        double cpu_elapsed = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
        double wall_elapsed = wall_end - wall_start;
        
        int threads = 1;
        #ifdef _OPENMP
        threads = omp_get_max_threads();
        #endif
        //printf("CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", cpu_elapsed, wall_elapsed, threads);
        
        if (!fq_nmod_is_zero(det, ctx)) {
            fq_mvpoly_add_term_fast(result, NULL, NULL, det);
        }
        
        fq_nmod_clear(det, ctx);
        fq_nmod_mat_clear(scalar_mat, ctx);
        
    } else if (npars == 1) {
        clock_t cpu_start = clock();
        double wall_start = get_wall_time();
        
        if (method == DET_METHOD_INTERPOLATION) {
            dixon_debug_log("  Method: interpolation\n");
            
            fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                           0, npars, ctx, res_deg_bound);
        } else if (method == DET_METHOD_HUANG) {
            dixon_debug_log("  Method: sparse interpolation\n");

            compute_fq_det_huang_interpolation(result, coeff_matrix, size);
        } else if (method == DET_METHOD_KRONECKER_NMOD) {
            dixon_debug_log("  Method: Bareiss\n");
            compute_fq_det_bareiss(result, coeff_matrix, size);
        } else if (method == DET_METHOD_BALANCED_SPLIT) {
            dixon_debug_log("  Method: balanced split Laplace (experimental)\n");
            compute_fq_det_balanced_split_experimental(result, coeff_matrix, size);
        } else {
            fq_mvpoly_init(result, 0, npars, ctx);
            
            fq_nmod_poly_mat_t poly_mat;
            fq_nmod_poly_mat_init(poly_mat, size, size, ctx);
            
            for (slong i = 0; i < size; i++) {
                for (slong j = 0; j < size; j++) {
                    fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(poly_mat, i, j);
                    fq_nmod_poly_zero(entry, ctx);
                    
                    for (slong t = 0; t < coeff_matrix[i][j].nterms; t++) {
                        slong deg = coeff_matrix[i][j].terms[t].par_exp ? 
                                   coeff_matrix[i][j].terms[t].par_exp[0] : 0;
                        fq_nmod_poly_set_coeff(entry, deg, 
                                              coeff_matrix[i][j].terms[t].coeff, ctx);
                    }
                }
            }
            
            fq_nmod_poly_t det_poly;
            fq_nmod_poly_init(det_poly, ctx);
            
            dixon_debug_log("  Method: fast HNF\n");
            
            fq_nmod_poly_mat_det_iter(det_poly, poly_mat, ctx);
            
            slong det_deg = fq_nmod_poly_degree(det_poly, ctx);
            if (det_deg >= 0) {
                for (slong i = 0; i <= det_deg; i++) {
                    fq_nmod_t coeff;
                    fq_nmod_init(coeff, ctx);
                    fq_nmod_poly_get_coeff(coeff, det_poly, i, ctx);
                    if (!fq_nmod_is_zero(coeff, ctx)) {
                        slong par_exp[1] = {i};
                        fq_mvpoly_add_term_fast(result, NULL, par_exp, coeff);
                    }
                    fq_nmod_clear(coeff, ctx);
                }
            }
            
            fq_nmod_poly_clear(det_poly, ctx);
            fq_nmod_poly_mat_clear(poly_mat, ctx);
        }
        
        clock_t cpu_end = clock();
        double wall_end = get_wall_time();
        double cpu_elapsed = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
        double wall_elapsed = wall_end - wall_start;
        
        int threads = 1;
        #ifdef _OPENMP
        threads = omp_get_max_threads();
        #endif
        // printf("CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", cpu_elapsed, wall_elapsed, threads);
        
    } else {
        clock_t cpu_start = clock();
        double wall_start = get_wall_time();
        switch (method) {
            case DET_METHOD_INTERPOLATION:
                dixon_debug_log("  Method: interpolation\n");
                
                fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                               0, npars, ctx, res_deg_bound);
                break;
                
            case DET_METHOD_RECURSIVE:
                dixon_debug_log("  Method: minor expansion\n");
                
                compute_fq_det_recursive(result, coeff_matrix, size);
                break;
                
            case DET_METHOD_KRONECKER:
                dixon_debug_log("  Method: fast HNF\n");
                
                compute_fq_det_kronecker(result, coeff_matrix, size);
                break;

            case DET_METHOD_KRONECKER_NMOD:
                dixon_debug_log("  Method: Bareiss\n");
                compute_fq_det_bareiss(result, coeff_matrix, size);
                break;
            case DET_METHOD_BALANCED_SPLIT:
                dixon_debug_log("  Method: balanced split Laplace (experimental)\n");
                compute_fq_det_balanced_split_experimental(result, coeff_matrix, size);
                break;

            case DET_METHOD_HUANG:
                dixon_debug_log("  Method: sparse interpolation\n");
                
                compute_fq_det_huang_interpolation(result, coeff_matrix, size);
                break;
                
            default:
                dixon_debug_log("  Method: interpolation (default)\n");
                fq_compute_det_by_interpolation(result, coeff_matrix, size,
                                               0, npars, ctx, res_deg_bound);
                break;
        }
        clock_t cpu_end = clock();
        double wall_end = get_wall_time();
        double cpu_elapsed = (double)(cpu_end - cpu_start) / CLOCKS_PER_SEC;
        double wall_elapsed = wall_end - wall_start;
        
        int threads = 1;
        #ifdef _OPENMP
        threads = omp_get_max_threads();
        #endif
        // printf("CPU time: %.3f seconds | Wall time: %.3f seconds | Threads: %d\n", cpu_elapsed, wall_elapsed, threads);
    }

    if (g_field_equation_reduction || g_field_equation_final_only) {
        fq_mvpoly_reduce_field_equation(result);
    }
}

/* Cheap first deflation pass for the Step 4 coefficient matrix.
 * For the single remaining parameter x, extract only row/column powers of x.
 * The returned exponent is a factor of the determinant and must be retained
 * by the caller as metadata; it is not discarded mathematically. */
static slong fq_mvpoly_x_valuation(const fq_mvpoly_t *poly)
{
    slong min_exp = LONG_MAX;

    if (!poly || poly->nterms == 0 || poly->npars != 1)
        return 0;

    for (slong t = 0; t < poly->nterms; t++) {
        slong exp = poly->terms[t].par_exp ? poly->terms[t].par_exp[0] : 0;
        if (exp < min_exp)
            min_exp = exp;
    }

    return min_exp == LONG_MAX ? 0 : min_exp;
}

static void fq_mvpoly_divide_by_x_power_inplace(fq_mvpoly_t *poly, slong power)
{
    if (!poly || power <= 0 || poly->npars != 1)
        return;

    for (slong t = 0; t < poly->nterms; t++) {
        if (poly->terms[t].par_exp)
            poly->terms[t].par_exp[0] -= power;
    }
}

static void fq_mvpoly_multiply_by_x_power_inplace(fq_mvpoly_t *poly, slong power)
{
    if (!poly || power <= 0 || poly->npars != 1)
        return;

    for (slong t = 0; t < poly->nterms; t++) {
        if (poly->terms[t].par_exp)
            poly->terms[t].par_exp[0] += power;
    }
}

static slong extract_fq_matrix_x_content(fq_mvpoly_t **matrix, slong size,
                                         slong npars)
{
    slong extracted = 0;

    if (!matrix || size <= 0 || npars != 1)
        return 0;

    /* Row content. */
    for (slong i = 0; i < size; i++) {
        slong row_power = 0;
        int have_nonzero = 0;

        for (slong j = 0; j < size; j++) {
            if (matrix[i][j].nterms > 0) {
                slong power = fq_mvpoly_x_valuation(&matrix[i][j]);
                if (!have_nonzero || power < row_power)
                    row_power = power;
                have_nonzero = 1;
            }
        }

        if (have_nonzero && row_power > 0) {
            for (slong j = 0; j < size; j++)
                fq_mvpoly_divide_by_x_power_inplace(&matrix[i][j], row_power);
            extracted += row_power;
        }
    }

    /* Column content after row normalization. */
    for (slong j = 0; j < size; j++) {
        slong col_power = 0;
        int have_nonzero = 0;

        for (slong i = 0; i < size; i++) {
            if (matrix[i][j].nterms > 0) {
                slong power = fq_mvpoly_x_valuation(&matrix[i][j]);
                if (!have_nonzero || power < col_power)
                    col_power = power;
                have_nonzero = 1;
            }
        }

        if (have_nonzero && col_power > 0) {
            for (slong i = 0; i < size; i++)
                fq_mvpoly_divide_by_x_power_inplace(&matrix[i][j], col_power);
            extracted += col_power;
        }
    }

    return extracted;
}

/* Same deflation pass for the full sparse coefficient matrix, before
 * maximal-rank row/column selection. */
static slong extract_fq_full_matrix_x_content(fq_mvpoly_t ***matrix,
                                              slong nrows, slong ncols,
                                              slong npars,
                                              slong *row_powers,
                                              slong *col_powers)
{
    slong extracted = 0;

    if (!matrix || nrows <= 0 || ncols <= 0 || npars != 1)
        return 0;

    for (slong i = 0; i < nrows; i++) {
        slong row_power = 0;
        int have_nonzero = 0;

        for (slong j = 0; j < ncols; j++) {
            if (matrix[i][j] && matrix[i][j]->nterms > 0) {
                slong power = fq_mvpoly_x_valuation(matrix[i][j]);
                if (!have_nonzero || power < row_power)
                    row_power = power;
                have_nonzero = 1;
            }
        }

        if (have_nonzero && row_power > 0) {
            for (slong j = 0; j < ncols; j++)
                fq_mvpoly_divide_by_x_power_inplace(matrix[i][j], row_power);
            if (row_powers)
                row_powers[i] += row_power;
            extracted += row_power;
        }
    }

    for (slong j = 0; j < ncols; j++) {
        slong col_power = 0;
        int have_nonzero = 0;

        for (slong i = 0; i < nrows; i++) {
            if (matrix[i][j] && matrix[i][j]->nterms > 0) {
                slong power = fq_mvpoly_x_valuation(matrix[i][j]);
                if (!have_nonzero || power < col_power)
                    col_power = power;
                have_nonzero = 1;
            }
        }

        if (have_nonzero && col_power > 0) {
            for (slong i = 0; i < nrows; i++)
                fq_mvpoly_divide_by_x_power_inplace(matrix[i][j], col_power);
            if (col_powers)
                col_powers[j] += col_power;
            extracted += col_power;
        }
    }

    return extracted;
}

// Extended tracker structure with pre-allocated workspace
// Initialize optimized tracker
static void unified_row_basis_tracker_init(unified_row_basis_tracker_t *tracker, 
                                            slong max_size, slong ncols, 
                                            field_ctx_t *ctx) {
    tracker->max_size = max_size;
    tracker->ncols = ncols;
    tracker->ctx = ctx;
    tracker->current_rank = 0;
    tracker->initialized = 1;
    tracker->workspace_initialized = 0;
    
    void *ctx_ptr = (ctx->field_id == FIELD_ID_NMOD) ?
                   (void*)&ctx->ctx.nmod_ctx :
                   // (ctx->field_id == FIELD_ID_FQ_ZECH) ?
                   // (void*)ctx->ctx.zech_ctx :
                   (void*)ctx->ctx.fq_ctx;
    
    // Allocate main storage
    tracker->reduced_rows = (field_elem_u*) flint_calloc(max_size * ncols, sizeof(field_elem_u));
    tracker->pivot_cols = (slong*) flint_calloc(max_size, sizeof(slong));
    tracker->selected_indices = (slong*) flint_calloc(max_size, sizeof(slong));
    
    // Pre-allocate workspace
    tracker->work_row = (field_elem_u*) flint_malloc(ncols * sizeof(field_elem_u));
    tracker->temp_vars = (field_elem_u*) flint_malloc(4 * sizeof(field_elem_u)); // factor, temp, pivot_val, neg_temp
    
    // Initialize all field elements
    for (slong i = 0; i < max_size * ncols; i++) {
        field_init_elem(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
        field_set_zero(&tracker->reduced_rows[i], ctx->field_id, ctx_ptr);
    }
    
    // Initialize workspace
    for (slong j = 0; j < ncols; j++) {
        field_init_elem(&tracker->work_row[j], ctx->field_id, ctx_ptr);
    }
    for (slong i = 0; i < 4; i++) {
        field_init_elem(&tracker->temp_vars[i], ctx->field_id, ctx_ptr);
    }
    tracker->workspace_initialized = 1;
    
    // Initialize pivot columns to -1
    for (slong i = 0; i < max_size; i++) {
        tracker->pivot_cols[i] = -1;
        tracker->selected_indices[i] = -1;
    }
}

// Clear optimized tracker
static void unified_row_basis_tracker_clear(unified_row_basis_tracker_t *tracker) {
    if (!tracker->initialized) return;
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ?
                   (void*)&tracker->ctx->ctx.nmod_ctx :
                   // (tracker->ctx->field_id == FIELD_ID_FQ_ZECH) ?
                   // (void*)tracker->ctx->ctx.zech_ctx :
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // Clear main storage
    if (tracker->reduced_rows) {
        for (slong i = 0; i < tracker->max_size * tracker->ncols; i++) {
            field_clear_elem(&tracker->reduced_rows[i], tracker->ctx->field_id, ctx_ptr);
        }
        flint_free(tracker->reduced_rows);
    }
    
    // Clear workspace
    if (tracker->workspace_initialized) {
        for (slong j = 0; j < tracker->ncols; j++) {
            field_clear_elem(&tracker->work_row[j], tracker->ctx->field_id, ctx_ptr);
        }
        for (slong i = 0; i < 4; i++) {
            field_clear_elem(&tracker->temp_vars[i], tracker->ctx->field_id, ctx_ptr);
        }
        flint_free(tracker->work_row);
        flint_free(tracker->temp_vars);
    }
    
    if (tracker->pivot_cols) flint_free(tracker->pivot_cols);
    if (tracker->selected_indices) flint_free(tracker->selected_indices);
    
    tracker->initialized = 0;
    tracker->workspace_initialized = 0;
}

// Optimized version of adding row to basis - mathematical logic completely unchanged
static int unified_try_add_row_to_basis(unified_row_basis_tracker_t *tracker, 
                                         const field_elem_u *unified_mat,
                                         slong new_row_idx, slong ncols) {
    if (!tracker->initialized || tracker->current_rank >= tracker->max_size) {
        return 0;
    }
    
    void *ctx_ptr = (tracker->ctx->field_id == FIELD_ID_NMOD) ?
                   (void*)&tracker->ctx->ctx.nmod_ctx :
                   // (tracker->ctx->field_id == FIELD_ID_FQ_ZECH) ?
                   // (void*)tracker->ctx->ctx.zech_ctx :
                   (void*)tracker->ctx->ctx.fq_ctx;
    
    // Use pre-allocated work row to avoid allocation each time
    field_elem_u *work_row = tracker->work_row;
    
    // Copy input row to work row (reuse allocated space)
    for (slong j = 0; j < ncols; j++) {
        field_set_elem(&work_row[j], &unified_mat[new_row_idx * ncols + j], 
                      tracker->ctx->field_id, ctx_ptr);
    }
    
    // Use pre-allocated temporary variables
    field_elem_u *factor = &tracker->temp_vars[0];
    field_elem_u *temp = &tracker->temp_vars[1];
    field_elem_u *pivot_val = &tracker->temp_vars[2];
    field_elem_u *neg_temp = &tracker->temp_vars[3];
    
    // Perform elimination for each existing basis vector (mathematical logic unchanged)
    for (slong i = 0; i < tracker->current_rank; i++) {
        slong pivot_col = tracker->pivot_cols[i];
        
        // Validate pivot column index
        if (pivot_col < 0 || pivot_col >= ncols) continue;
        
        // If work row is non-zero at pivot position, perform elimination
        if (!field_is_zero(&work_row[pivot_col], tracker->ctx->field_id, ctx_ptr)) {
            // Get pivot value of basis vector
            slong base_idx = i * ncols + pivot_col;
            field_set_elem(pivot_val, &tracker->reduced_rows[base_idx], 
                          tracker->ctx->field_id, ctx_ptr);
            
            // Calculate elimination factor = work_row[pivot_col] / pivot_val
            field_inv(temp, pivot_val, tracker->ctx->field_id, ctx_ptr);
            field_mul(factor, &work_row[pivot_col], temp, tracker->ctx->field_id, ctx_ptr);
            
            // Perform elimination: work_row -= factor * basis_row[i]
            // Optimization: pre-determine if neg_temp is needed to avoid conditional branches in inner loop
            int use_direct_add = (tracker->ctx->field_id >= FIELD_ID_GF28 && 
                                 tracker->ctx->field_id <= FIELD_ID_GF2128);
            
            for (slong j = 0; j < ncols; j++) {
                slong idx = i * ncols + j;
                field_mul(temp, factor, &tracker->reduced_rows[idx], 
                         tracker->ctx->field_id, ctx_ptr);
                
                if (use_direct_add) {
                    // For GF(2^n), subtraction equals addition
                    field_add(&work_row[j], &work_row[j], temp, 
                             tracker->ctx->field_id, ctx_ptr);
                } else {
                    // For other fields, use actual subtraction
                    field_neg(neg_temp, temp, tracker->ctx->field_id, ctx_ptr);
                    field_add(&work_row[j], &work_row[j], neg_temp, 
                             tracker->ctx->field_id, ctx_ptr);
                }
            }
        }
    }
    
    // Find first non-zero position (mathematical logic unchanged)
    slong first_nonzero = -1;
    for (slong j = 0; j < ncols; j++) {
        if (!field_is_zero(&work_row[j], tracker->ctx->field_id, ctx_ptr)) {
            first_nonzero = j;
            break;
        }
    }
    
    // If all zeros, then linearly dependent
    if (first_nonzero == -1) {
        return 0;  // Work row doesn't need cleanup because it's pre-allocated
    }
    
    // Normalize work row (make first non-zero element 1)
    field_inv(temp, &work_row[first_nonzero], tracker->ctx->field_id, ctx_ptr);
    
    // Store normalized row to basis
    slong base_row_start = tracker->current_rank * ncols;
    for (slong j = 0; j < ncols; j++) {
        slong idx = base_row_start + j;
        if (j < first_nonzero) {
            // Positions before pivot should be 0
            field_set_zero(&tracker->reduced_rows[idx], tracker->ctx->field_id, ctx_ptr);
        } else {
            // Normalize and store
            field_mul(&tracker->reduced_rows[idx], &work_row[j], temp, 
                     tracker->ctx->field_id, ctx_ptr);
        }
    }
    
    // Update tracking information (mathematical logic unchanged)
    tracker->pivot_cols[tracker->current_rank] = first_nonzero;
    tracker->selected_indices[tracker->current_rank] = new_row_idx;
    tracker->current_rank++;
    
    return 1;
}

typedef struct {
    mp_limb_t *reduced_rows;
    slong *pivot_cols;
    slong *selected_indices;
    mp_limb_t *work_row;
    slong current_rank;
    slong max_size;
    slong ncols;
    nmod_t mod;
    int initialized;
} nmod_row_basis_tracker_t;

static void nmod_row_basis_tracker_init(nmod_row_basis_tracker_t *tracker,
                                        slong max_size,
                                        slong ncols,
                                        const nmod_t *mod)
{
    tracker->max_size = max_size;
    tracker->ncols = ncols;
    tracker->current_rank = 0;
    tracker->initialized = 1;
    tracker->mod = *mod;
    tracker->reduced_rows = (mp_limb_t *) flint_calloc((size_t) max_size * ncols, sizeof(mp_limb_t));
    tracker->pivot_cols = (slong *) flint_malloc((size_t) max_size * sizeof(slong));
    tracker->selected_indices = (slong *) flint_malloc((size_t) max_size * sizeof(slong));
    tracker->work_row = (mp_limb_t *) flint_malloc((size_t) ncols * sizeof(mp_limb_t));

    for (slong i = 0; i < max_size; i++) {
        tracker->pivot_cols[i] = -1;
        tracker->selected_indices[i] = -1;
    }
}

static void nmod_row_basis_tracker_clear(nmod_row_basis_tracker_t *tracker)
{
    if (!tracker->initialized) {
        return;
    }

    flint_free(tracker->reduced_rows);
    flint_free(tracker->pivot_cols);
    flint_free(tracker->selected_indices);
    flint_free(tracker->work_row);
    tracker->initialized = 0;
}

static int nmod_try_add_row_to_basis(nmod_row_basis_tracker_t *tracker,
                                     const field_elem_u *unified_mat,
                                     slong new_row_idx)
{
    if (!tracker->initialized || tracker->current_rank >= tracker->max_size) {
        return 0;
    }

    const mp_limb_t mod_n = tracker->mod.n;
    const ulong mod_ninv = tracker->mod.ninv;
    mp_limb_t *work_row = tracker->work_row;
    for (slong j = 0; j < tracker->ncols; j++) {
        work_row[j] = unified_mat[new_row_idx * tracker->ncols + j].nmod;
    }

    for (slong i = 0; i < tracker->current_rank; i++) {
        slong pivot_col = tracker->pivot_cols[i];
        mp_limb_t factor;
        mp_limb_t neg_factor;
        mp_limb_t *basis_row;

        if (pivot_col < 0 || pivot_col >= tracker->ncols) {
            continue;
        }
        factor = work_row[pivot_col];
        if (factor == 0) {
            continue;
        }

        basis_row = tracker->reduced_rows + i * tracker->ncols;
        neg_factor = nmod_neg(factor, tracker->mod);
        _nmod_vec_scalar_addmul_nmod(work_row + pivot_col,
                                     basis_row + pivot_col,
                                     tracker->ncols - pivot_col,
                                     neg_factor,
                                     tracker->mod);
    }

    slong first_nonzero = -1;
    for (slong j = 0; j < tracker->ncols; j++) {
        if (work_row[j] != 0) {
            first_nonzero = j;
            break;
        }
    }
    if (first_nonzero < 0) {
        return 0;
    }

    mp_limb_t pivot = work_row[first_nonzero];
    mp_limb_t pivot_inv = n_invmod(pivot, mod_n);
    mp_limb_t *dest = tracker->reduced_rows + tracker->current_rank * tracker->ncols;

    for (slong j = 0; j < first_nonzero; j++) {
        dest[j] = 0;
    }
    for (slong j = first_nonzero; j < tracker->ncols; j++) {
        dest[j] = n_mulmod2_preinv(work_row[j], pivot_inv, mod_n, mod_ninv);
    }

    tracker->pivot_cols[tracker->current_rank] = first_nonzero;
    tracker->selected_indices[tracker->current_rank] = new_row_idx;
    tracker->current_rank++;
    return 1;
}

static int nmod_try_add_row_to_basis_raw(nmod_row_basis_tracker_t *tracker,
                                         const mp_limb_t *mat,
                                         slong new_row_idx)
{
    if (!tracker->initialized || tracker->current_rank >= tracker->max_size) {
        return 0;
    }

    const mp_limb_t mod_n = tracker->mod.n;
    const ulong mod_ninv = tracker->mod.ninv;
    mp_limb_t *work_row = tracker->work_row;
    for (slong j = 0; j < tracker->ncols; j++) {
        work_row[j] = mat[new_row_idx * tracker->ncols + j];
    }

    for (slong i = 0; i < tracker->current_rank; i++) {
        slong pivot_col = tracker->pivot_cols[i];
        mp_limb_t factor;
        mp_limb_t neg_factor;
        mp_limb_t *basis_row;

        if (pivot_col < 0 || pivot_col >= tracker->ncols) {
            continue;
        }
        factor = work_row[pivot_col];
        if (factor == 0) {
            continue;
        }

        basis_row = tracker->reduced_rows + i * tracker->ncols;
        neg_factor = nmod_neg(factor, tracker->mod);
        _nmod_vec_scalar_addmul_nmod(work_row + pivot_col,
                                     basis_row + pivot_col,
                                     tracker->ncols - pivot_col,
                                     neg_factor,
                                     tracker->mod);
    }

    slong first_nonzero = -1;
    for (slong j = 0; j < tracker->ncols; j++) {
        if (work_row[j] != 0) {
            first_nonzero = j;
            break;
        }
    }
    if (first_nonzero < 0) {
        return 0;
    }

    mp_limb_t pivot = work_row[first_nonzero];
    mp_limb_t pivot_inv = n_invmod(pivot, mod_n);
    mp_limb_t *dest = tracker->reduced_rows + tracker->current_rank * tracker->ncols;

    for (slong j = 0; j < first_nonzero; j++) {
        dest[j] = 0;
    }
    for (slong j = first_nonzero; j < tracker->ncols; j++) {
        dest[j] = n_mulmod2_preinv(work_row[j], pivot_inv, mod_n, mod_ninv);
    }

    tracker->pivot_cols[tracker->current_rank] = first_nonzero;
    tracker->selected_indices[tracker->current_rank] = new_row_idx;
    tracker->current_rank++;
    return 1;
}

static void find_pivot_rows_nmod_dense(slong **selected_rows_out,
                                       slong *num_selected,
                                       const mp_limb_t *mat,
                                       slong nrows,
                                       slong ncols,
                                       const nmod_t *mod)
{
    slong max_rank = FLINT_MIN(nrows, ncols);
    slong rank = 0;
    slong *selected_rows;
    slong *perm;
    nmod_mat_t work;

    if (nrows == 0 || ncols == 0) {
        *selected_rows_out = NULL;
        *num_selected = 0;
        return;
    }

    selected_rows = (slong *) flint_malloc((size_t) max_rank * sizeof(slong));
    perm = (slong *) flint_malloc((size_t) nrows * sizeof(slong));
    nmod_mat_init(work, nrows, ncols, mod->n);

    for (slong i = 0; i < nrows; i++) {
        perm[i] = i;
        for (slong j = 0; j < ncols; j++) {
            nmod_mat_entry(work, i, j) = mat[i * ncols + j];
        }
    }

    rank = nmod_mat_lu(perm, work, 0);
    for (slong i = 0; i < rank; i++) {
        selected_rows[i] = perm[i];
    }

    nmod_mat_clear(work);
    flint_free(perm);

    if (rank > 0) {
        *selected_rows_out = (slong *) flint_realloc(selected_rows, (size_t) rank * sizeof(slong));
    } else {
        flint_free(selected_rows);
        *selected_rows_out = NULL;
    }
    *num_selected = rank;
}

// Helper function to compute maximum degree (not total degree) of a polynomial
static slong compute_fq_polynomial_total_degree(fq_mvpoly_t *poly, slong npars) {
    if (poly == NULL || poly->nterms == 0) {
        return 0;
    }
    
    slong max_degree = 0;
    
    for (slong t = 0; t < poly->nterms; t++) {
        // Check variable degrees
        if (poly->terms[t].var_exp) {
            for (slong v = 0; v < poly->nvars; v++) {
                if (poly->terms[t].var_exp[v] > max_degree) {
                    max_degree = poly->terms[t].var_exp[v];
                }
            }
        }
        
        // Check parameter degrees
        if (poly->terms[t].par_exp && npars > 0) {
            for (slong p = 0; p < npars; p++) {
                if (poly->terms[t].par_exp[p] > max_degree) {
                    max_degree = poly->terms[t].par_exp[p];
                }
            }
        }
    }
    
    return max_degree;
}

static slong **build_entry_degree_cache(fq_mvpoly_t ***full_matrix,
                                       slong nrows,
                                       slong ncols,
                                       slong npars)
{
    slong **degree_cache = (slong **) flint_malloc((size_t) nrows * sizeof(slong *));

    for (slong i = 0; i < nrows; i++) {
        degree_cache[i] = (slong *) flint_malloc((size_t) ncols * sizeof(slong));
    }

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (slong i = 0; i < nrows; i++) {
        for (slong j = 0; j < ncols; j++) {
            degree_cache[i][j] = compute_fq_polynomial_total_degree(full_matrix[i][j], npars);
        }
    }

    return degree_cache;
}

static void free_entry_degree_cache(slong **degree_cache, slong nrows)
{
    if (degree_cache == NULL) {
        return;
    }

    for (slong i = 0; i < nrows; i++) {
        flint_free(degree_cache[i]);
    }
    flint_free(degree_cache);
}

static void print_fq_degree_vector_with_names(const char *label,
                                              const slong *degrees,
                                              slong count,
                                              char **names,
                                              int is_dual,
                                              int is_parameter) {
    static const char default_var_names[] = {'x', 'y', 'z', 'w', 'v', 'u'};
    static const char default_par_names[] = {'a', 'b', 'c', 'd'};

    printf("  %s:", label);
    if (count <= 0) {
        printf(" (none)\n");
        return;
    }

    for (slong i = 0; i < count; i++) {
        printf(" ");
        if (is_dual) {
            printf("~");
        }

        if (names && names[i]) {
            printf("%s:%ld", names[i], degrees[i]);
        } else if (is_parameter) {
            if (i < (slong) (sizeof(default_par_names) / sizeof(default_par_names[0]))) {
                printf("%c:%ld", default_par_names[i], degrees[i]);
            } else {
                printf("p_%ld:%ld", i, degrees[i]);
            }
        } else {
            if (i < (slong) (sizeof(default_var_names) / sizeof(default_var_names[0]))) {
                printf("%c:%ld", default_var_names[i], degrees[i]);
            } else {
                printf("x_%ld:%ld", i, degrees[i]);
            }
        }
    }
    printf("\n");
}

static void print_dixon_poly_actual_degrees(const fq_mvpoly_t *poly,
                                            slong nvars,
                                            slong npars,
                                            char **var_names,
                                            char **par_names) {
    slong *orig_deg = (slong*) flint_calloc(nvars > 0 ? nvars : 1, sizeof(slong));
    slong *dual_deg = (slong*) flint_calloc(nvars > 0 ? nvars : 1, sizeof(slong));
    slong *par_deg = (slong*) flint_calloc(npars > 0 ? npars : 1, sizeof(slong));

    if (poly != NULL && poly->nterms > 0) {
        for (slong t = 0; t < poly->nterms; t++) {
            if (poly->terms[t].var_exp) {
                for (slong j = 0; j < nvars && j < poly->nvars; j++) {
                    if (poly->terms[t].var_exp[j] > orig_deg[j]) {
                        orig_deg[j] = poly->terms[t].var_exp[j];
                    }
                }
                for (slong j = 0; j < nvars && (nvars + j) < poly->nvars; j++) {
                    if (poly->terms[t].var_exp[nvars + j] > dual_deg[j]) {
                        dual_deg[j] = poly->terms[t].var_exp[nvars + j];
                    }
                }
            }

            if (poly->terms[t].par_exp) {
                for (slong j = 0; j < npars && j < poly->npars; j++) {
                    if (poly->terms[t].par_exp[j] > par_deg[j]) {
                        par_deg[j] = poly->terms[t].par_exp[j];
                    }
                }
            }
        }
    }

    printf("  Dixon polynomial actual degrees:\n");
    print_fq_degree_vector_with_names("original vars", orig_deg, nvars, var_names, 0, 0);
    print_fq_degree_vector_with_names("dual vars", dual_deg, nvars, var_names, 1, 0);
    print_fq_degree_vector_with_names("parameter vars", par_deg, npars, par_names, 0, 1);

    flint_free(orig_deg);
    flint_free(dual_deg);
    flint_free(par_deg);
}

// Compute row maximum total degree
static slong compute_fq_row_max_total_degree(fq_mvpoly_t **matrix_row, slong ncols, slong npars) {
    slong max_degree = -1;
    
    for (slong j = 0; j < ncols; j++) {
        slong poly_deg = compute_fq_polynomial_total_degree(matrix_row[j], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// Compute column maximum total degree
static slong compute_fq_col_max_total_degree(fq_mvpoly_t ***matrix, slong col_idx, slong nrows, slong npars) {
    slong max_degree = -1;
    
    for (slong i = 0; i < nrows; i++) {
        slong poly_deg = compute_fq_polynomial_total_degree(matrix[i][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}


// Compute maximum total degree of a column in selected rows submatrix
static slong compute_fq_selected_rows_col_max_total_degree(fq_mvpoly_t ***full_matrix, 
                                                           slong *selected_rows, 
                                                           slong num_selected_rows,
                                                           slong col_idx, 
                                                           slong npars) {
    slong max_degree = -1;
    
    for (slong i = 0; i < num_selected_rows; i++) {
        slong row_idx = selected_rows[i];
        slong poly_deg = compute_fq_polynomial_total_degree(full_matrix[row_idx][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

// Compute maximum total degree of a row in selected columns submatrix
static slong compute_fq_selected_cols_row_max_total_degree(fq_mvpoly_t ***full_matrix, 
                                                           slong row_idx,
                                                           slong *selected_cols, 
                                                           slong num_selected_cols,
                                                           slong npars) {
    slong max_degree = -1;
    
    for (slong j = 0; j < num_selected_cols; j++) {
        slong col_idx = selected_cols[j];
        slong poly_deg = compute_fq_polynomial_total_degree(full_matrix[row_idx][col_idx], npars);
        if (poly_deg > max_degree) {
            max_degree = poly_deg;
        }
    }
    
    return max_degree;
}

/* Reorder an already selected square minor using the same entry-degree
   priorities as the degree-aware selector.  The row and column sets remain
   unchanged, so this cannot invalidate the verified rank certificate. */
static void reorder_fq_selected_minor_by_degree(fq_mvpoly_t ***full_matrix,
                                                slong *selected_rows,
                                                slong *selected_cols,
                                                slong size,
                                                slong npars)
{
    fq_index_degree_pair *row_degrees;
    fq_index_degree_pair *col_degrees;

    if (full_matrix == NULL || selected_rows == NULL || selected_cols == NULL ||
        size <= 1)
        return;

    row_degrees = (fq_index_degree_pair *)
        flint_malloc((size_t) size * sizeof(fq_index_degree_pair));
    col_degrees = (fq_index_degree_pair *)
        flint_malloc((size_t) size * sizeof(fq_index_degree_pair));

    for (slong i = 0; i < size; i++) {
        row_degrees[i].index = selected_rows[i];
        row_degrees[i].degree =
            compute_fq_selected_cols_row_max_total_degree(
                full_matrix, selected_rows[i], selected_cols, size, npars);
    }
    qsort(row_degrees, (size_t) size, sizeof(fq_index_degree_pair),
          compare_fq_degrees);
    for (slong i = 0; i < size; i++) selected_rows[i] = row_degrees[i].index;

    for (slong j = 0; j < size; j++) {
        col_degrees[j].index = selected_cols[j];
        col_degrees[j].degree =
            compute_fq_selected_rows_col_max_total_degree(
                full_matrix, selected_rows, size, selected_cols[j], npars);
    }
    qsort(col_degrees, (size_t) size, sizeof(fq_index_degree_pair),
          compare_fq_degrees);
    for (slong j = 0; j < size; j++) selected_cols[j] = col_degrees[j].index;

    dixon_debug_log("  Reordered predicted minor by entry degree "
                    "(rows %ld..%ld, columns %ld..%ld)\n",
                    row_degrees[0].degree, row_degrees[size - 1].degree,
                    col_degrees[0].degree, col_degrees[size - 1].degree);
    flint_free(row_degrees);
    flint_free(col_degrees);
}

// Check if two index arrays are identical
static int indices_equal(slong *indices1, slong *indices2, slong size) {
    for (slong i = 0; i < size; i++) {
        if (indices1[i] != indices2[i]) {
            return 0;
        }
    }
    return 1;
}

static int select_ksy_columns_from_transposed(slong **selected_cols_out,
                                              slong *num_selected_cols_out,
                                              int *ksy_condition_met,
                                              const field_elem_u *transposed,
                                              slong ncols,
                                              slong target_rank,
                                              const fq_index_degree_pair *col_degrees,
                                              slong ksy_constant_col,
                                              field_ctx_t *ctx)
{
    unified_row_basis_tracker_t nonconst_tracker;
    unified_row_basis_tracker_t verify_tracker;
    slong *candidate_cols;
    slong nonconst_rank;
    slong output_size = 0;
    int constant_independent = 0;

    *selected_cols_out = NULL;
    *num_selected_cols_out = 0;
    *ksy_condition_met = 0;

    if (ksy_constant_col < 0 || ksy_constant_col >= ncols || target_rank <= 0) {
        return 0;
    }

    candidate_cols = (slong *) flint_malloc((size_t) target_rank * sizeof(slong));
    unified_row_basis_tracker_init(&nonconst_tracker, target_rank, target_rank, ctx);

    for (slong j = 0; j < ncols && nonconst_tracker.current_rank < target_rank; j++) {
        slong col_idx = col_degrees[j].index;
        if (col_idx == ksy_constant_col) {
            continue;
        }
        unified_try_add_row_to_basis(&nonconst_tracker, transposed, col_idx, target_rank);
    }

    nonconst_rank = nonconst_tracker.current_rank;
    unified_row_basis_tracker_init(&verify_tracker,
                                   FLINT_MIN(nonconst_rank + 1, target_rank),
                                   target_rank,
                                   ctx);
    for (slong i = 0; i < nonconst_rank && verify_tracker.current_rank < nonconst_rank; i++) {
        unified_try_add_row_to_basis(&verify_tracker,
                                     transposed,
                                     nonconst_tracker.selected_indices[i],
                                     target_rank);
    }
    if (verify_tracker.current_rank == nonconst_rank) {
        constant_independent =
            unified_try_add_row_to_basis(&verify_tracker,
                                         transposed,
                                         ksy_constant_col,
                                         target_rank);
    }

    unified_row_basis_tracker_clear(&verify_tracker);
    if (constant_independent) {
        output_size = nonconst_rank + 1;
        candidate_cols[0] = ksy_constant_col;
        if (nonconst_rank > 0) {
            memcpy(candidate_cols + 1,
                   nonconst_tracker.selected_indices,
                   (size_t) nonconst_rank * sizeof(slong));
        }
        *ksy_condition_met = 1;
    } else {
        output_size = FLINT_MIN(nonconst_rank, target_rank - 1);
        if (output_size > 0) {
            memcpy(candidate_cols,
                   nonconst_tracker.selected_indices,
                   (size_t) output_size * sizeof(slong));
        }
    }

    unified_row_basis_tracker_clear(&nonconst_tracker);

    if (output_size <= 0) {
        flint_free(candidate_cols);
        return 0;
    }

    *selected_cols_out = candidate_cols;
    *num_selected_cols_out = output_size;
    return 1;
}

static int select_ksy_columns_from_transposed_nmod(slong **selected_cols_out,
                                                   slong *num_selected_cols_out,
                                                   int *ksy_condition_met,
                                                   const mp_limb_t *transposed,
                                                   slong ncols,
                                                   slong target_rank,
                                                   const fq_index_degree_pair *col_degrees,
                                                   slong ksy_constant_col,
                                                   const nmod_t *mod)
{
    nmod_row_basis_tracker_t nonconst_tracker;
    nmod_row_basis_tracker_t verify_tracker;
    slong *candidate_cols;
    slong nonconst_rank;
    slong output_size = 0;
    int constant_independent = 0;

    *selected_cols_out = NULL;
    *num_selected_cols_out = 0;
    *ksy_condition_met = 0;

    if (ksy_constant_col < 0 || ksy_constant_col >= ncols || target_rank <= 0) {
        return 0;
    }

    candidate_cols = (slong *) flint_malloc((size_t) target_rank * sizeof(slong));
    nmod_row_basis_tracker_init(&nonconst_tracker, target_rank, target_rank, mod);

    for (slong j = 0; j < ncols && nonconst_tracker.current_rank < target_rank; j++) {
        slong col_idx = col_degrees[j].index;
        if (col_idx == ksy_constant_col) {
            continue;
        }
        nmod_try_add_row_to_basis_raw(&nonconst_tracker, transposed, col_idx);
    }

    nonconst_rank = nonconst_tracker.current_rank;
    nmod_row_basis_tracker_init(&verify_tracker,
                                FLINT_MIN(nonconst_rank + 1, target_rank),
                                target_rank,
                                mod);
    for (slong i = 0; i < nonconst_rank && verify_tracker.current_rank < nonconst_rank; i++) {
        nmod_try_add_row_to_basis_raw(&verify_tracker,
                                      transposed,
                                      nonconst_tracker.selected_indices[i]);
    }
    if (verify_tracker.current_rank == nonconst_rank) {
        constant_independent = nmod_try_add_row_to_basis_raw(&verify_tracker,
                                                             transposed,
                                                             ksy_constant_col);
    }

    nmod_row_basis_tracker_clear(&verify_tracker);
    if (constant_independent) {
        output_size = nonconst_rank + 1;
        candidate_cols[0] = ksy_constant_col;
        if (nonconst_rank > 0) {
            memcpy(candidate_cols + 1,
                   nonconst_tracker.selected_indices,
                   (size_t) nonconst_rank * sizeof(slong));
        }
        *ksy_condition_met = 1;
    } else {
        output_size = FLINT_MIN(nonconst_rank, target_rank - 1);
        if (output_size > 0) {
            memcpy(candidate_cols,
                   nonconst_tracker.selected_indices,
                   (size_t) output_size * sizeof(slong));
        }
    }

    nmod_row_basis_tracker_clear(&nonconst_tracker);

    if (output_size <= 0) {
        flint_free(candidate_cols);
        return 0;
    }

    *selected_cols_out = candidate_cols;
    *num_selected_cols_out = output_size;
    return 1;
}
// Find pivot rows for maximal rank submatrix - ensure linear independence
// More efficient version: using incremental rank checking
static void find_pivot_rows_nmod_fixed(slong **selected_rows_out, slong *num_selected,
                                       const nmod_mat_t mat) {
    
    slong nrows = nmod_mat_nrows(mat);
    slong ncols = nmod_mat_ncols(mat);
    nmod_t mod = mat->mod;
    
    if (nrows == 0 || ncols == 0) {
        *selected_rows_out = NULL;
        *num_selected = 0;
        return;
    }
    
    slong min_dim = FLINT_MIN(nrows, ncols);
    slong *selected_rows = (slong*) flint_malloc(min_dim * sizeof(slong));
    slong rank = 0;
    
    // Create working copy
    nmod_mat_t A;
    nmod_mat_init(A, nrows, ncols, mod.n);
    nmod_mat_set(A, mat);
    
    // Row permutation tracking
    slong *P = (slong*) flint_malloc(nrows * sizeof(slong));
    for (slong i = 0; i < nrows; i++) {
        P[i] = i;
    }
    
    // PROPER LU DECOMPOSITION using FLINT-style algorithm
    slong current_row = 0;
    for (slong col = 0; col < ncols && current_row < nrows; col++) {
        // printf("%d %d\n",col,ncols);
        // Find pivot in current column
        slong pivot_row = -1;
        for (slong i = current_row; i < nrows; i++) {
            if (nmod_mat_entry(A, i, col) != 0) {
                pivot_row = i;
                break;
            }
        }
        
        // No pivot found in this column
        if (pivot_row == -1) {
            continue;
        }
        
        // Record this pivot row
        selected_rows[rank] = P[pivot_row];
        rank++;
        
        // Swap rows if needed - use FLINT's efficient method
        if (pivot_row != current_row) {
            // Swap permutation
            slong temp_idx = P[current_row];
            P[current_row] = P[pivot_row];
            P[pivot_row] = temp_idx;
            
            // Swap matrix rows efficiently
            //nn_ptr temp_row = A->rows[current_row];
            //A->rows[current_row] = A->rows[pivot_row];
            //A->rows[pivot_row] = temp_row;
            // Swap matrix rows using compatible API
			
            for (slong j = 0; j < ncols; j++) {
                mp_limb_t temp_val = nmod_mat_entry(A, current_row, j);
                nmod_mat_entry(A, current_row, j) = nmod_mat_entry(A, pivot_row, j);
                nmod_mat_entry(A, pivot_row, j) = temp_val;
            }
        }
        
        // Eliminate below pivot using FLINT vectorized operations
        mp_limb_t pivot = nmod_mat_entry(A, current_row, col);
        mp_limb_t pivot_inv = n_invmod(pivot, mod.n);
        
        for (slong i = current_row + 1; i < nrows; i++) {
            mp_limb_t factor = nmod_mat_entry(A, i, col);
            if (factor == 0) continue;
            
            factor = n_mulmod2_preinv(factor, pivot_inv, mod.n, mod.ninv);
            mp_limb_t neg_factor = nmod_neg(factor, mod);
            
            // Use FLINT's vectorized subtraction
            /*
            slong remaining_cols = ncols - col;
            if (remaining_cols > 0) {
                _nmod_vec_scalar_addmul_nmod(A->rows[i] + col,
                                           A->rows[current_row] + col,
                                           remaining_cols, neg_factor, mod);
            }
            */
            mp_limb_t *row_i = &nmod_mat_entry(A, i, 0);
            mp_limb_t *row_current = &nmod_mat_entry(A, current_row, 0);
            if (ncols - col > 0) {
                _nmod_vec_scalar_addmul_nmod(row_i + col, 
                                           row_current + col, 
                                           ncols - col, 
                                           neg_factor, mod);
            }
        }
        
        current_row++;
    }
    
    // Set output
    *num_selected = rank;
    if (rank > 0) {
        *selected_rows_out = (slong*) flint_realloc(selected_rows, rank * sizeof(slong));
    } else {
        flint_free(selected_rows);
        *selected_rows_out = NULL;
    }
    
    // Cleanup
    nmod_mat_clear(A);
    flint_free(P);
}
// ============================================================================
// Optimized find_pivot_rows_simple - directly calls nmod version
// ============================================================================
static void find_pivot_rows_simple(slong **selected_rows_out, slong *num_selected,
                                   const field_elem_u *unified_mat,
                                   slong nrows, slong ncols,
                                   field_ctx_t *ctx) {
    
    // ============================================================================
    // Prime field fast path: directly use nmod_fixed implementation (optimal performance)
    // ============================================================================
    if (ctx->field_id == FIELD_ID_NMOD) {
        //printf("Using direct nmod_fixed implementation for prime field (optimal performance)\n");
        clock_t start = clock();
        
        // Convert to nmod_mat format
        nmod_mat_t nmod_mat;
        nmod_mat_init(nmod_mat, nrows, ncols, ctx->ctx.nmod_ctx.n);
        
        // Efficient data copy (direct access to nmod field)
        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                nmod_mat_entry(nmod_mat, i, j) = unified_mat[i * ncols + j].nmod;
            }
        }
        
        // Directly use nmod_fixed version (avoid adaptive layer selection overhead)
        find_pivot_rows_nmod_fixed(selected_rows_out, num_selected, nmod_mat);
        
        clock_t end = clock();
        double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
        
        nmod_mat_clear(nmod_mat);
        //printf("Prime field computation completed in %.4f seconds using direct nmod path\n", elapsed);
        return;
    }
    
    // ============================================================================
    // Non-prime field: use unified interface implementation
    // ============================================================================
    //printf("Using unified interface for non-prime field computation\n");
    
    void *ctx_ptr = // (ctx->field_id == FIELD_ID_FQ_ZECH) ? 
                   // (void*)ctx->ctx.zech_ctx : 
                   (void*)ctx->ctx.fq_ctx;
    
    slong max_rank = FLINT_MIN(nrows, ncols);
    slong *selected_rows = (slong*) flint_malloc(max_rank * sizeof(slong));
    slong rank = 0;
    
    // Lightweight linear independence checker
    typedef struct {
        field_elem_u *rows;     
        slong *pivot_positions; 
        slong count;            
        slong ncols;
    } basis_tracker_t;
    
    basis_tracker_t tracker;
    tracker.rows = (field_elem_u*) flint_calloc(max_rank * ncols, sizeof(field_elem_u));
    tracker.pivot_positions = (slong*) flint_malloc(max_rank * sizeof(slong));
    tracker.count = 0;
    tracker.ncols = ncols;
    
    // Initialize all field elements
    for (slong i = 0; i < max_rank * ncols; i++) {
        field_init_elem(&tracker.rows[i], ctx->field_id, ctx_ptr);
    }
    
    clock_t start = clock();
    
    // Process row by row
    for (slong row = 0; row < nrows && rank < max_rank; row++) {
        // Create test vector
        field_elem_u *test_vec = (field_elem_u*) flint_malloc(ncols * sizeof(field_elem_u));
        for (slong j = 0; j < ncols; j++) {
            field_init_elem(&test_vec[j], ctx->field_id, ctx_ptr);
            field_set_elem(&test_vec[j], &unified_mat[row * ncols + j], ctx->field_id, ctx_ptr);
        }
        
        // Eliminate using existing basis vectors
        for (slong i = 0; i < tracker.count; i++) {
            slong pivot_col = tracker.pivot_positions[i];
            if (!field_is_zero(&test_vec[pivot_col], ctx->field_id, ctx_ptr)) {
                field_elem_u factor;
                field_init_elem(&factor, ctx->field_id, ctx_ptr);
                field_set_elem(&factor, &test_vec[pivot_col], ctx->field_id, ctx_ptr);
                
                for (slong j = 0; j < ncols; j++) {
                    field_elem_u temp;
                    field_init_elem(&temp, ctx->field_id, ctx_ptr);
                    field_mul(&temp, &factor, &tracker.rows[i * ncols + j], ctx->field_id, ctx_ptr);
                    field_sub(&test_vec[j], &test_vec[j], &temp, ctx->field_id, ctx_ptr);
                    field_clear_elem(&temp, ctx->field_id, ctx_ptr);
                }
                
                field_clear_elem(&factor, ctx->field_id, ctx_ptr);
            }
        }
        
        // Find first non-zero position
        slong pivot_pos = -1;
        for (slong j = 0; j < ncols; j++) {
            if (!field_is_zero(&test_vec[j], ctx->field_id, ctx_ptr)) {
                pivot_pos = j;
                break;
            }
        }
        
        if (pivot_pos >= 0) {
            // Linearly independent, add to basis
            selected_rows[rank] = row;
            tracker.pivot_positions[tracker.count] = pivot_pos;
            
            // Normalize and store
            field_elem_u inv;
            field_init_elem(&inv, ctx->field_id, ctx_ptr);
            field_inv(&inv, &test_vec[pivot_pos], ctx->field_id, ctx_ptr);
            
            for (slong j = 0; j < ncols; j++) {
                field_mul(&tracker.rows[tracker.count * ncols + j], &test_vec[j], &inv, ctx->field_id, ctx_ptr);
            }
            
            field_clear_elem(&inv, ctx->field_id, ctx_ptr);
            tracker.count++;
            rank++;
        }
        
        // Clean up test vector
        for (slong j = 0; j < ncols; j++) {
            field_clear_elem(&test_vec[j], ctx->field_id, ctx_ptr);
        }
        flint_free(test_vec);
    }
    
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;
    
    // Set output
    if (rank > 0) {
        *selected_rows_out = (slong*) flint_malloc(rank * sizeof(slong));
        memcpy(*selected_rows_out, selected_rows, rank * sizeof(slong));
    } else {
        *selected_rows_out = NULL;
    }
    *num_selected = rank;
    
    // Cleanup
    for (slong i = 0; i < max_rank * ncols; i++) {
        field_clear_elem(&tracker.rows[i], ctx->field_id, ctx_ptr);
    }
    flint_free(tracker.rows);
    flint_free(tracker.pivot_positions);
    flint_free(selected_rows);
    
    //printf("Non-prime field computation completed in %.4f seconds using unified interface\n", elapsed);
}


// Optimized find_fq_optimal_maximal_rank_submatrix
static void find_fq_optimal_maximal_rank_submatrix_nmod(fq_mvpoly_t ***full_matrix,
                                                        slong nrows,
                                                        slong ncols,
                                                        slong **row_indices_out,
                                                        slong **col_indices_out,
                                                        slong *num_rows,
                                                        slong *num_cols,
                                                        slong npars,
                                                        slong ksy_constant_col,
                                                        const fq_nmod_ctx_t ctx,
                                                        slong max_selection_attempts)
{
    slong accepted_size = 0;
    slong *accepted_rows = NULL;
    slong *accepted_cols = NULL;
    slong accepted_score = -1;
    int selection_stable = 0;
    const int require_ksy_precondition = (ksy_constant_col >= 0);
    nmod_t mod = ctx->mod;

    for (slong selection_attempt = 0;
         selection_attempt < max_selection_attempts && !selection_stable;
         selection_attempt++) {
        fq_nmod_t *param_vals = (fq_nmod_t *) flint_malloc((size_t) npars * sizeof(fq_nmod_t));
        mp_limb_t *unified_mat = (mp_limb_t *) flint_malloc((size_t) nrows * ncols * sizeof(mp_limb_t));
        slong *current_row_indices = NULL;
        slong *current_col_indices = NULL;
        slong *prev_row_indices = NULL;
        slong *prev_col_indices = NULL;
        slong current_size = 0;
        slong prev_size = 0;
        const slong MAX_ITERATIONS = 3;
        slong iteration = 0;
        int converged = 0;
        int ksy_condition_met = 0;
        clock_t degree_cache_cpu_start = clock();
        double degree_cache_wall_start = get_wall_time();
        slong **entry_degree_cache;

        init_evaluation_parameters(param_vals, npars, ctx, selection_attempt);
        entry_degree_cache = build_entry_degree_cache(full_matrix, nrows, ncols, npars);
        dixon_maybe_print_step_detail_time("Step 3 degree-cache build",
                                           degree_cache_cpu_start,
                                           degree_cache_wall_start);

        clock_t unified_mat_cpu_start = clock();
        double unified_mat_wall_start = get_wall_time();
#ifdef _OPENMP
        #pragma omp parallel
        {
            fq_nmod_t val;
            fq_nmod_init(val, ctx);
            #pragma omp for schedule(static)
            for (slong i = 0; i < nrows; i++) {
                for (slong j = 0; j < ncols; j++) {
                    slong idx = i * ncols + j;
                    if (full_matrix[i][j] == NULL) {
                        unified_mat[idx] = 0;
                    } else {
                        evaluate_fq_mvpoly_at_nmod_params_direct(val, full_matrix[i][j], param_vals, ctx);
                        unified_mat[idx] = nmod_poly_get_coeff_ui(val, 0);
                    }
                }
            }
            fq_nmod_clear(val, ctx);
        }
#else
        {
            fq_nmod_t val;
            fq_nmod_init(val, ctx);
            for (slong i = 0; i < nrows; i++) {
                for (slong j = 0; j < ncols; j++) {
                    slong idx = i * ncols + j;
                    if (full_matrix[i][j] == NULL) {
                        unified_mat[idx] = 0;
                    } else {
                        evaluate_fq_mvpoly_at_nmod_params_direct(val, full_matrix[i][j], param_vals, ctx);
                        unified_mat[idx] = nmod_poly_get_coeff_ui(val, 0);
                    }
                }
            }
            fq_nmod_clear(val, ctx);
        }
#endif
        dixon_maybe_print_step_detail_time("Step 3 unified matrix build",
                                           unified_mat_cpu_start,
                                           unified_mat_wall_start);

        while (iteration < MAX_ITERATIONS && !converged) {
            slong row_rank_selected = 0;

            if (iteration > 0) {
                prev_row_indices = (slong *) flint_malloc((size_t) current_size * sizeof(slong));
                prev_col_indices = (slong *) flint_malloc((size_t) current_size * sizeof(slong));
                memcpy(prev_row_indices, current_row_indices, (size_t) current_size * sizeof(slong));
                memcpy(prev_col_indices, current_col_indices, (size_t) current_size * sizeof(slong));
                prev_size = current_size;
            }

            if (iteration == 0) {
                slong *selected_rows = NULL;
                slong num_selected = 0;
                clock_t initial_row_select_cpu_start = clock();
                double initial_row_select_wall_start = get_wall_time();

                find_pivot_rows_nmod_dense(&selected_rows, &num_selected, unified_mat, nrows, ncols, &mod);
                dixon_maybe_print_step_detail_time("Step 3 initial row selection",
                                                   initial_row_select_cpu_start,
                                                   initial_row_select_wall_start);

                current_size = num_selected;
                row_rank_selected = num_selected;
                current_row_indices = (slong *) flint_malloc((size_t) current_size * sizeof(slong));
                memcpy(current_row_indices, selected_rows, (size_t) current_size * sizeof(slong));

                if (current_size > g_matrix_transpose_threshold) {
                    mp_limb_t *transposed_mat = (mp_limb_t *) flint_malloc((size_t) ncols * current_size * sizeof(mp_limb_t));
                    slong *selected_cols = NULL;
                    slong num_selected_cols = 0;
                    clock_t large_transpose_cpu_start = clock();
                    double large_transpose_wall_start = get_wall_time();

#ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
#endif
                    for (slong i = 0; i < current_size; i++) {
                        slong orig_row = current_row_indices[i];
                        for (slong j = 0; j < ncols; j++) {
                            transposed_mat[j * current_size + i] = unified_mat[orig_row * ncols + j];
                        }
                    }
                    dixon_maybe_print_step_detail_time("Step 3 large transpose build",
                                                       large_transpose_cpu_start,
                                                       large_transpose_wall_start);

                    clock_t initial_col_select_cpu_start = clock();
                    double initial_col_select_wall_start = get_wall_time();
                    find_pivot_rows_nmod_dense(&selected_cols, &num_selected_cols,
                                               transposed_mat, ncols, current_size, &mod);
                    dixon_maybe_print_step_detail_time("Step 3 initial column selection",
                                                       initial_col_select_cpu_start,
                                                       initial_col_select_wall_start);

                    current_col_indices = (slong *) flint_malloc((size_t) num_selected_cols * sizeof(slong));
                    memcpy(current_col_indices, selected_cols, (size_t) num_selected_cols * sizeof(slong));
                    current_size = FLINT_MIN(current_size, num_selected_cols);

                    flint_free(transposed_mat);
                    flint_free(selected_rows);
                    flint_free(selected_cols);

                    converged = 1;
                    break;
                }

                flint_free(selected_rows);
            } else {
                fq_index_degree_pair *row_degrees = (fq_index_degree_pair *) flint_malloc((size_t) nrows * sizeof(fq_index_degree_pair));
                mp_limb_t *col_submat = (mp_limb_t *) flint_malloc((size_t) nrows * current_size * sizeof(mp_limb_t));
                slong selected_rank = 0;
                slong *selected_rank_indices = NULL;
                clock_t row_degree_cpu_start = clock();
                double row_degree_wall_start = get_wall_time();

#ifdef _OPENMP
                #pragma omp parallel for schedule(static)
#endif
                for (slong i = 0; i < nrows; i++) {
                    row_degrees[i].index = i;
                    row_degrees[i].degree = -1;
                    for (slong j = 0; j < current_size; j++) {
                        slong col_idx = current_col_indices[j];
                        if (entry_degree_cache[i][col_idx] > row_degrees[i].degree) {
                            row_degrees[i].degree = entry_degree_cache[i][col_idx];
                        }
                    }
                }
                dixon_maybe_print_step_detail_time("Step 3 row-degree scan",
                                                   row_degree_cpu_start,
                                                   row_degree_wall_start);

                qsort(row_degrees, (size_t) nrows, sizeof(fq_index_degree_pair), compare_fq_degrees);

                clock_t col_submat_cpu_start = clock();
                double col_submat_wall_start = get_wall_time();
#ifdef _OPENMP
                #pragma omp parallel for schedule(static)
#endif
                for (slong i = 0; i < nrows; i++) {
                    for (slong j = 0; j < current_size; j++) {
                        col_submat[i * current_size + j] = unified_mat[i * ncols + current_col_indices[j]];
                    }
                }
                dixon_maybe_print_step_detail_time("Step 3 column submatrix build",
                                                   col_submat_cpu_start,
                                                   col_submat_wall_start);

                clock_t row_basis_cpu_start = clock();
                double row_basis_wall_start = get_wall_time();
                nmod_row_basis_tracker_t row_tracker_nmod;
                nmod_row_basis_tracker_init(&row_tracker_nmod, current_size, current_size, &mod);
                for (slong i = 0; i < nrows && row_tracker_nmod.current_rank < current_size; i++) {
                    slong row_idx = row_degrees[i].index;
                    nmod_try_add_row_to_basis_raw(&row_tracker_nmod, col_submat, row_idx);
                }
                selected_rank = row_tracker_nmod.current_rank;
                selected_rank_indices = (slong *) flint_malloc((size_t) selected_rank * sizeof(slong));
                memcpy(selected_rank_indices,
                       row_tracker_nmod.selected_indices,
                       (size_t) selected_rank * sizeof(slong));
                nmod_row_basis_tracker_clear(&row_tracker_nmod);
                dixon_maybe_print_step_detail_time("Step 3 row basis selection",
                                                   row_basis_cpu_start,
                                                   row_basis_wall_start);

                flint_free(current_row_indices);
                current_row_indices = selected_rank_indices;
                row_rank_selected = selected_rank;
                current_size = selected_rank;

                flint_free(col_submat);
                flint_free(row_degrees);
            }

            fq_index_degree_pair *col_degrees = (fq_index_degree_pair *) flint_malloc((size_t) ncols * sizeof(fq_index_degree_pair));
            mp_limb_t *transposed = (mp_limb_t *) flint_malloc((size_t) ncols * current_size * sizeof(mp_limb_t));
            slong transposed_rows = current_size;
            slong col_rank_selected = 0;
            slong *selected_cols = NULL;
            clock_t col_degree_cpu_start = clock();
            double col_degree_wall_start = get_wall_time();

#ifdef _OPENMP
            #pragma omp parallel for schedule(static)
#endif
            for (slong j = 0; j < ncols; j++) {
                col_degrees[j].index = j;
                col_degrees[j].degree = -1;
                for (slong i = 0; i < current_size; i++) {
                    slong row_idx = current_row_indices[i];
                    if (entry_degree_cache[row_idx][j] > col_degrees[j].degree) {
                        col_degrees[j].degree = entry_degree_cache[row_idx][j];
                    }
                }
            }
            dixon_maybe_print_step_detail_time("Step 3 column-degree scan",
                                               col_degree_cpu_start,
                                               col_degree_wall_start);

            qsort(col_degrees, (size_t) ncols, sizeof(fq_index_degree_pair), compare_fq_degrees);

            clock_t transposed_cpu_start = clock();
            double transposed_wall_start = get_wall_time();
#ifdef _OPENMP
            #pragma omp parallel for schedule(static)
#endif
            for (slong i = 0; i < current_size; i++) {
                for (slong j = 0; j < ncols; j++) {
                    transposed[j * current_size + i] = unified_mat[current_row_indices[i] * ncols + j];
                }
            }
            dixon_maybe_print_step_detail_time("Step 3 transposed build",
                                               transposed_cpu_start,
                                               transposed_wall_start);

            if (require_ksy_precondition && current_size > 0) {
                int helper_ksy_met = 0;
                if (select_ksy_columns_from_transposed_nmod(&selected_cols,
                                                            &col_rank_selected,
                                                            &helper_ksy_met,
                                                            transposed,
                                                            ncols,
                                                            current_size,
                                                            col_degrees,
                                                            ksy_constant_col,
                                                            &mod)) {
                    ksy_condition_met = helper_ksy_met;
                }
            }

            if (!ksy_condition_met) {
                clock_t col_basis_cpu_start = clock();
                double col_basis_wall_start = get_wall_time();
                nmod_row_basis_tracker_t col_tracker_nmod;
                nmod_row_basis_tracker_init(&col_tracker_nmod, ncols, current_size, &mod);
                for (slong j = 0; j < ncols && col_tracker_nmod.current_rank < current_size; j++) {
                    slong col_idx = col_degrees[j].index;
                    nmod_try_add_row_to_basis_raw(&col_tracker_nmod, transposed, col_idx);
                }
                col_rank_selected = col_tracker_nmod.current_rank;
                selected_cols = (slong *) flint_malloc((size_t) col_rank_selected * sizeof(slong));
                memcpy(selected_cols,
                       col_tracker_nmod.selected_indices,
                       (size_t) col_rank_selected * sizeof(slong));
                nmod_row_basis_tracker_clear(&col_tracker_nmod);
                dixon_maybe_print_step_detail_time("Step 3 column basis selection",
                                                   col_basis_cpu_start,
                                                   col_basis_wall_start);
            }

            if (iteration == 0) {
                current_col_indices = (slong *) flint_malloc((size_t) col_rank_selected * sizeof(slong));
            } else {
                flint_free(current_col_indices);
                current_col_indices = (slong *) flint_malloc((size_t) col_rank_selected * sizeof(slong));
            }
            memcpy(current_col_indices, selected_cols, (size_t) col_rank_selected * sizeof(slong));
            current_size = FLINT_MIN(current_size, col_rank_selected);
            flint_free(selected_cols);
            flint_free(transposed);
            flint_free(col_degrees);

            if (iteration > 0) {
                if (current_size == prev_size &&
                    current_size == FLINT_MIN(col_rank_selected, row_rank_selected) &&
                    indices_equal(current_row_indices, prev_row_indices, current_size) &&
                    indices_equal(current_col_indices, prev_col_indices, current_size)) {
                    converged = 1;
                }

                flint_free(prev_row_indices);
                flint_free(prev_col_indices);
                prev_row_indices = NULL;
                prev_col_indices = NULL;
            }

            iteration++;
            (void) transposed_rows;
        }

        free_entry_degree_cache(entry_degree_cache, nrows);

        slong final_rank = evaluate_selected_submatrix_rank_nmod(
            full_matrix, current_row_indices, current_col_indices, current_size, param_vals, ctx);
        slong verify_rank = current_size;
        int final_ksy_ok = !require_ksy_precondition;
        int verify_ksy_ok = !require_ksy_precondition;
        clock_t verification_cpu_start = clock();
        double verification_wall_start = get_wall_time();

        if (g_dixon_step3_second_verification) {
            fq_nmod_t *verify_vals = (fq_nmod_t *) flint_malloc((size_t) npars * sizeof(fq_nmod_t));
            init_evaluation_parameters(verify_vals, npars, ctx, selection_attempt + 1);
            verify_rank = evaluate_selected_submatrix_rank_nmod(
                full_matrix, current_row_indices, current_col_indices, current_size, verify_vals, ctx);
            if (require_ksy_precondition) {
                final_ksy_ok = selected_rows_satisfy_ksy_precondition_nmod(
                    full_matrix, current_row_indices, current_size, ncols, ksy_constant_col, param_vals, ctx);
                verify_ksy_ok = selected_rows_satisfy_ksy_precondition_nmod(
                    full_matrix, current_row_indices, current_size, ncols, ksy_constant_col, verify_vals, ctx);
            }
            clear_evaluation_parameters(verify_vals, npars, ctx);
        } else if (require_ksy_precondition) {
            final_ksy_ok = selected_rows_satisfy_ksy_precondition_nmod(
                full_matrix, current_row_indices, current_size, ncols, ksy_constant_col, param_vals, ctx);
            verify_ksy_ok = final_ksy_ok;
        }
        dixon_maybe_print_step_detail_time("Step 3 verification",
                                           verification_cpu_start,
                                           verification_wall_start);

        slong candidate_score = 0;
        if (final_rank == current_size) candidate_score++;
        if (g_dixon_step3_second_verification && verify_rank == current_size) candidate_score++;
        if (require_ksy_precondition) {
            if (final_ksy_ok) candidate_score += 2;
            if (g_dixon_step3_second_verification && verify_ksy_ok) candidate_score += 2;
        }

        if (current_size > accepted_size ||
            (current_size == accepted_size && candidate_score > accepted_score)) {
            if (accepted_rows) flint_free(accepted_rows);
            if (accepted_cols) flint_free(accepted_cols);
            accepted_rows = (slong *) flint_malloc((size_t) current_size * sizeof(slong));
            accepted_cols = (slong *) flint_malloc((size_t) current_size * sizeof(slong));
            memcpy(accepted_rows, current_row_indices, (size_t) current_size * sizeof(slong));
            memcpy(accepted_cols, current_col_indices, (size_t) current_size * sizeof(slong));
            accepted_size = current_size;
            accepted_score = candidate_score;
        }

        if (candidate_score == (require_ksy_precondition
                                ? (g_dixon_step3_second_verification ? 6 : 2)
                                : (g_dixon_step3_second_verification ? 2 : 1))) {
            selection_stable = 1;
        } else if (require_ksy_precondition && (!final_ksy_ok || !verify_ksy_ok)) {
            dixon_debug_log("  Candidate submatrix failed the KSY precondition; retrying with a new specialization.\n");
        } else {
            dixon_debug_log("  Submatrix verification failed; retrying with a new specialization.\n");
        }

        if (prev_row_indices) flint_free(prev_row_indices);
        if (prev_col_indices) flint_free(prev_col_indices);
        flint_free(current_row_indices);
        flint_free(current_col_indices);
        clear_evaluation_parameters(param_vals, npars, ctx);
        flint_free(unified_mat);
    }

    if (!accepted_rows || !accepted_cols) {
        *row_indices_out = NULL;
        *col_indices_out = NULL;
        *num_rows = 0;
        *num_cols = 0;
        return;
    }

    if (!selection_stable) {
        if (require_ksy_precondition) {
            dixon_info_log("Warning: using best available submatrix after %ld attempts (KSY precondition not certified for every test).\n",
                           max_selection_attempts);
        } else {
            dixon_info_log("Warning: using best available submatrix after %ld attempts.\n",
                           max_selection_attempts);
        }
    }

    *row_indices_out = accepted_rows;
    *col_indices_out = accepted_cols;
    *num_rows = accepted_size;
    *num_cols = accepted_size;
}

void find_fq_optimal_maximal_rank_submatrix(fq_mvpoly_t ***full_matrix,
                                           slong nrows, slong ncols,
                                           slong **row_indices_out, 
                                           slong **col_indices_out,
                                           slong *num_rows, slong *num_cols,
                                           slong npars,
                                           slong ksy_constant_col) {
    // Get context
    const fq_nmod_ctx_struct *ctx = NULL;
    for (slong i = 0; i < nrows && !ctx; i++) {
        for (slong j = 0; j < ncols && !ctx; j++) {
            if (full_matrix[i][j] != NULL) {
                ctx = full_matrix[i][j]->ctx;
                break;
            }
        }
    }
    
    const int use_extension_specialization = (fq_nmod_ctx_degree(ctx) == 1 && fq_nmod_ctx_prime(ctx) <= 2 && npars > 0);
    const slong MAX_SELECTION_ATTEMPTS = use_extension_specialization ? 3 : 2;

    field_ctx_t selection_ctx;
    fq_nmod_ctx_t extension_eval_ctx;
    if (use_extension_specialization) {
        fmpz_t selection_prime;
        fmpz_init_set_ui(selection_prime, fq_nmod_ctx_prime(ctx));
        fq_nmod_ctx_init(extension_eval_ctx, selection_prime, 4, "u");
        fmpz_clear(selection_prime);
        field_ctx_init(&selection_ctx, extension_eval_ctx);
    } else {
        /* Step 3 rank selection is faster with plain fq_nmod than Zech on medium extension fields. */
        field_ctx_init_enhanced(&selection_ctx, ctx, 0);
    }

    if (selection_ctx.field_id == FIELD_ID_NMOD) {
        find_fq_optimal_maximal_rank_submatrix_nmod(full_matrix, nrows, ncols,
                                                    row_indices_out, col_indices_out,
                                                    num_rows, num_cols,
                                                    npars, ksy_constant_col,
                                                    ctx, MAX_SELECTION_ATTEMPTS);
        field_ctx_clear(&selection_ctx);
        if (use_extension_specialization) {
            fq_nmod_ctx_clear(extension_eval_ctx);
        }
        return;
    }

    void *ctx_ptr = (selection_ctx.field_id == FIELD_ID_NMOD) ?
                   (void*)&selection_ctx.ctx.nmod_ctx :
                   // (selection_ctx.field_id == FIELD_ID_FQ_ZECH) ?
                   // (void*)selection_ctx.ctx.zech_ctx :
                   (void*)selection_ctx.ctx.fq_ctx;

    slong accepted_size = 0;
    slong *accepted_rows = NULL;
    slong *accepted_cols = NULL;
    slong accepted_score = -1;
    int selection_stable = 0;
    const int require_ksy_precondition = (ksy_constant_col >= 0);

    for (slong selection_attempt = 0;
         selection_attempt < MAX_SELECTION_ATTEMPTS && !selection_stable;
         selection_attempt++) {
        fq_nmod_t *param_vals = NULL;
        if (use_extension_specialization) {
            param_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
            init_extension_evaluation_parameters(param_vals, npars,
                                                 extension_eval_ctx,
                                                 selection_attempt);
        } else {
            param_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
            init_evaluation_parameters(param_vals, npars, ctx, selection_attempt);
        }

        field_elem_u *unified_mat = (field_elem_u*) flint_malloc(nrows * ncols * sizeof(field_elem_u));
        slong *current_row_indices = NULL;
        slong *current_col_indices = NULL;
        slong *prev_row_indices = NULL;
        slong *prev_col_indices = NULL;
        slong current_size = 0;
        slong prev_size = 0;
        const slong MAX_ITERATIONS = 10;
        slong iteration = 0;
        int converged = 0;
        int ksy_condition_met = 0;
        clock_t degree_cache_cpu_start = clock();
        double degree_cache_wall_start = get_wall_time();
        slong **entry_degree_cache = build_entry_degree_cache(full_matrix, nrows, ncols, npars);
        dixon_maybe_print_step_detail_time("Step 3 degree-cache build",
                                           degree_cache_cpu_start,
                                           degree_cache_wall_start);

        clock_t unified_mat_cpu_start = clock();
        double unified_mat_wall_start = get_wall_time();

#ifdef _OPENMP
        #pragma omp parallel
        {
            fq_nmod_t val;
            fq_nmod_init(val, use_extension_specialization ? extension_eval_ctx : ctx);
            #pragma omp for schedule(static)
            for (slong i = 0; i < nrows; i++) {
                for (slong j = 0; j < ncols; j++) {
                    slong idx = i * ncols + j;
                    field_init_elem(&unified_mat[idx], selection_ctx.field_id, ctx_ptr);

                    if (use_extension_specialization) {
                        if (full_matrix[i][j] == NULL) {
                            fq_nmod_zero(val, extension_eval_ctx);
                        } else {
                            evaluate_fq_mvpoly_at_extension_params(val, full_matrix[i][j],
                                                                   param_vals, extension_eval_ctx);
                        }
                    } else {
                        if (full_matrix[i][j] == NULL) {
                            fq_nmod_zero(val, ctx);
                        } else {
                            evaluate_fq_mvpoly_at_params(val, full_matrix[i][j], param_vals);
                        }
                    }
                    fq_nmod_to_field_elem(&unified_mat[idx], val, &selection_ctx);
                }
            }
            fq_nmod_clear(val, use_extension_specialization ? extension_eval_ctx : ctx);
        }
#else
        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                slong idx = i * ncols + j;
                field_init_elem(&unified_mat[idx], selection_ctx.field_id, ctx_ptr);

                if (use_extension_specialization) {
                    fq_nmod_t val;
                    fq_nmod_init(val, extension_eval_ctx);
                    if (full_matrix[i][j] == NULL) {
                        fq_nmod_zero(val, extension_eval_ctx);
                    } else {
                        evaluate_fq_mvpoly_at_extension_params(val, full_matrix[i][j],
                                                               param_vals, extension_eval_ctx);
                    }
                    fq_nmod_to_field_elem(&unified_mat[idx], val, &selection_ctx);
                    fq_nmod_clear(val, extension_eval_ctx);
                } else {
                    fq_nmod_t val;
                    fq_nmod_init(val, ctx);
                    if (full_matrix[i][j] == NULL) {
                        fq_nmod_zero(val, ctx);
                    } else {
                        evaluate_fq_mvpoly_at_params(val, full_matrix[i][j], param_vals);
                    }
                    fq_nmod_to_field_elem(&unified_mat[idx], val, &selection_ctx);
                    fq_nmod_clear(val, ctx);
                }
            }
        }
#endif
        dixon_maybe_print_step_detail_time("Step 3 unified matrix build",
                                           unified_mat_cpu_start,
                                           unified_mat_wall_start);

        clock_t iter_start = clock();
        clock_t initial_row_select_cpu_start = 0;
        double initial_row_select_wall_start = 0.0;

        while (iteration < MAX_ITERATIONS && !converged) {
            if (iteration > 0) {
                prev_row_indices = (slong*) flint_malloc(current_size * sizeof(slong));
                prev_col_indices = (slong*) flint_malloc(current_size * sizeof(slong));
                memcpy(prev_row_indices, current_row_indices, current_size * sizeof(slong));
                memcpy(prev_col_indices, current_col_indices, current_size * sizeof(slong));
                prev_size = current_size;
            }

            slong row_rank_selected = 0;

            if (iteration == 0) {
                slong *selected_rows = NULL;
                slong num_selected = 0;
                initial_row_select_cpu_start = clock();
                initial_row_select_wall_start = get_wall_time();

                find_pivot_rows_simple(&selected_rows, &num_selected,
                                      unified_mat, nrows, ncols, &selection_ctx);
                dixon_maybe_print_step_detail_time("Step 3 initial row selection",
                                                   initial_row_select_cpu_start,
                                                   initial_row_select_wall_start);

                current_size = num_selected;
                row_rank_selected = num_selected;
                current_row_indices = (slong*) flint_malloc(current_size * sizeof(slong));
                memcpy(current_row_indices, selected_rows, current_size * sizeof(slong));

                if (current_size > g_matrix_transpose_threshold) {
                    field_elem_u *transposed_mat = (field_elem_u*) flint_malloc(ncols * current_size * sizeof(field_elem_u));
                    clock_t large_transpose_cpu_start = clock();
                    double large_transpose_wall_start = get_wall_time();

                    for (slong i = 0; i < ncols * current_size; i++) {
                        field_init_elem(&transposed_mat[i], selection_ctx.field_id, ctx_ptr);
                    }

#ifdef _OPENMP
                    #pragma omp parallel for schedule(static)
#endif
                    for (slong i = 0; i < current_size; i++) {
                        slong orig_row = current_row_indices[i];
                        for (slong j = 0; j < ncols; j++) {
                            slong src_idx = orig_row * ncols + j;
                            slong dst_idx = j * current_size + i;
                            field_set_elem(&transposed_mat[dst_idx], &unified_mat[src_idx],
                                          selection_ctx.field_id, ctx_ptr);
                        }
                    }
                    dixon_maybe_print_step_detail_time("Step 3 large transpose build",
                                                       large_transpose_cpu_start,
                                                       large_transpose_wall_start);

                    slong *selected_cols = NULL;
                    slong num_selected_cols = 0;
                    clock_t initial_col_select_cpu_start = clock();
                    double initial_col_select_wall_start = get_wall_time();

                    find_pivot_rows_simple(&selected_cols, &num_selected_cols,
                                          transposed_mat, ncols, current_size, &selection_ctx);
                    dixon_maybe_print_step_detail_time("Step 3 initial column selection",
                                                       initial_col_select_cpu_start,
                                                       initial_col_select_wall_start);

                    current_col_indices = (slong*) flint_malloc(num_selected_cols * sizeof(slong));
                    memcpy(current_col_indices, selected_cols, num_selected_cols * sizeof(slong));
                    current_size = FLINT_MIN(current_size, num_selected_cols);

                    for (slong i = 0; i < ncols * num_selected; i++) {
                        field_clear_elem(&transposed_mat[i], selection_ctx.field_id, ctx_ptr);
                    }
                    flint_free(transposed_mat);
                    flint_free(selected_rows);
                    flint_free(selected_cols);

                    converged = 1;
                    break;
                } else {
                    flint_free(selected_rows);
                }
            } else {
                fq_index_degree_pair *row_degrees = (fq_index_degree_pair*) flint_malloc(nrows * sizeof(fq_index_degree_pair));
                clock_t row_degree_cpu_start = clock();
                double row_degree_wall_start = get_wall_time();

                #ifdef _OPENMP
                #pragma omp parallel for schedule(static)
                #endif
                for (slong i = 0; i < nrows; i++) {
                    row_degrees[i].index = i;
                    row_degrees[i].degree = -1;
                    for (slong j = 0; j < current_size; j++) {
                        slong col_idx = current_col_indices[j];
                        if (entry_degree_cache[i][col_idx] > row_degrees[i].degree) {
                            row_degrees[i].degree = entry_degree_cache[i][col_idx];
                        }
                    }
                }
                dixon_maybe_print_step_detail_time("Step 3 row-degree scan",
                                                   row_degree_cpu_start,
                                                   row_degree_wall_start);

                qsort(row_degrees, nrows, sizeof(fq_index_degree_pair), compare_fq_degrees);

                field_elem_u *col_submat = (field_elem_u*) flint_malloc(nrows * current_size * sizeof(field_elem_u));
                slong col_submat_cols = current_size;
                clock_t col_submat_cpu_start = clock();
                double col_submat_wall_start = get_wall_time();
                for (slong i = 0; i < nrows * current_size; i++) {
                    field_init_elem(&col_submat[i], selection_ctx.field_id, ctx_ptr);
                }

#ifdef _OPENMP
                #pragma omp parallel for schedule(static)
#endif
                for (slong i = 0; i < nrows; i++) {
                    for (slong j = 0; j < current_size; j++) {
                        slong src_idx = i * ncols + current_col_indices[j];
                        slong dst_idx = i * current_size + j;
                        field_set_elem(&col_submat[dst_idx], &unified_mat[src_idx],
                                      selection_ctx.field_id, ctx_ptr);
                    }
                }
                dixon_maybe_print_step_detail_time("Step 3 column submatrix build",
                                                   col_submat_cpu_start,
                                                   col_submat_wall_start);

                clock_t row_basis_cpu_start = clock();
                double row_basis_wall_start = get_wall_time();
                slong selected_rank = 0;
                slong *selected_rank_indices = NULL;

                if (selection_ctx.field_id == FIELD_ID_NMOD) {
                    nmod_row_basis_tracker_t row_tracker_nmod;
                    nmod_row_basis_tracker_init(&row_tracker_nmod, current_size, current_size,
                                                &selection_ctx.ctx.nmod_ctx);
                    for (slong i = 0; i < nrows && row_tracker_nmod.current_rank < current_size; i++) {
                        slong row_idx = row_degrees[i].index;
                        nmod_try_add_row_to_basis(&row_tracker_nmod, col_submat, row_idx);
                    }
                    selected_rank = row_tracker_nmod.current_rank;
                    selected_rank_indices = (slong*) flint_malloc(selected_rank * sizeof(slong));
                    memcpy(selected_rank_indices, row_tracker_nmod.selected_indices,
                           selected_rank * sizeof(slong));
                    nmod_row_basis_tracker_clear(&row_tracker_nmod);
                } else {
                    unified_row_basis_tracker_t row_tracker;
                    unified_row_basis_tracker_init(&row_tracker, current_size, current_size, &selection_ctx);
                    for (slong i = 0; i < nrows && row_tracker.current_rank < current_size; i++) {
                        slong row_idx = row_degrees[i].index;
                        unified_try_add_row_to_basis(&row_tracker, col_submat, row_idx, current_size);
                    }
                    selected_rank = row_tracker.current_rank;
                    selected_rank_indices = (slong*) flint_malloc(selected_rank * sizeof(slong));
                    memcpy(selected_rank_indices, row_tracker.selected_indices,
                           selected_rank * sizeof(slong));
                    unified_row_basis_tracker_clear(&row_tracker);
                }
                dixon_maybe_print_step_detail_time("Step 3 row basis selection",
                                                   row_basis_cpu_start,
                                                   row_basis_wall_start);

                flint_free(current_row_indices);
                current_row_indices = selected_rank_indices;
                row_rank_selected = selected_rank;
                current_size = selected_rank;

                for (slong i = 0; i < nrows * col_submat_cols; i++) {
                    field_clear_elem(&col_submat[i], selection_ctx.field_id, ctx_ptr);
                }
                flint_free(col_submat);
                flint_free(row_degrees);
            }

            fq_index_degree_pair *col_degrees = (fq_index_degree_pair*) flint_malloc(ncols * sizeof(fq_index_degree_pair));
            clock_t col_degree_cpu_start = clock();
            double col_degree_wall_start = get_wall_time();

            #ifdef _OPENMP
            #pragma omp parallel for schedule(static)
            #endif
            for (slong j = 0; j < ncols; j++) {
                col_degrees[j].index = j;
                col_degrees[j].degree = -1;
                for (slong i = 0; i < current_size; i++) {
                    slong row_idx = current_row_indices[i];
                    if (entry_degree_cache[row_idx][j] > col_degrees[j].degree) {
                        col_degrees[j].degree = entry_degree_cache[row_idx][j];
                    }
                }
            }
            dixon_maybe_print_step_detail_time("Step 3 column-degree scan",
                                               col_degree_cpu_start,
                                               col_degree_wall_start);

            qsort(col_degrees, ncols, sizeof(fq_index_degree_pair), compare_fq_degrees);

            field_elem_u *transposed = (field_elem_u*) flint_malloc(ncols * current_size * sizeof(field_elem_u));
            slong transposed_rows = current_size;
            clock_t transposed_cpu_start = clock();
            double transposed_wall_start = get_wall_time();
            for (slong i = 0; i < ncols * current_size; i++) {
                field_init_elem(&transposed[i], selection_ctx.field_id, ctx_ptr);
            }

#ifdef _OPENMP
            #pragma omp parallel for schedule(static)
#endif
            for (slong i = 0; i < current_size; i++) {
                for (slong j = 0; j < ncols; j++) {
                    slong src_idx = current_row_indices[i] * ncols + j;
                    slong dst_idx = j * current_size + i;
                    field_set_elem(&transposed[dst_idx], &unified_mat[src_idx],
                                  selection_ctx.field_id, ctx_ptr);
                }
            }
            dixon_maybe_print_step_detail_time("Step 3 transposed build",
                                               transposed_cpu_start,
                                               transposed_wall_start);

            slong col_rank_selected = 0;
            slong *selected_cols = NULL;

            if (require_ksy_precondition && current_size > 0) {
                int helper_ksy_met = 0;
                if (select_ksy_columns_from_transposed(&selected_cols,
                                                       &col_rank_selected,
                                                       &helper_ksy_met,
                                                       transposed,
                                                       ncols,
                                                       current_size,
                                                       col_degrees,
                                                       ksy_constant_col,
                                                       &selection_ctx)) {
                    ksy_condition_met = helper_ksy_met;
                }
            }

            if (!ksy_condition_met) {
                clock_t col_basis_cpu_start = clock();
                double col_basis_wall_start = get_wall_time();
                if (selection_ctx.field_id == FIELD_ID_NMOD) {
                    nmod_row_basis_tracker_t col_tracker_nmod;
                    nmod_row_basis_tracker_init(&col_tracker_nmod, ncols, current_size,
                                                &selection_ctx.ctx.nmod_ctx);
                    for (slong j = 0; j < ncols && col_tracker_nmod.current_rank < current_size; j++) {
                        slong col_idx = col_degrees[j].index;
                        nmod_try_add_row_to_basis(&col_tracker_nmod, transposed, col_idx);
                    }
                    col_rank_selected = col_tracker_nmod.current_rank;
                    selected_cols = (slong*) flint_malloc(col_rank_selected * sizeof(slong));
                    memcpy(selected_cols, col_tracker_nmod.selected_indices,
                           col_rank_selected * sizeof(slong));
                    nmod_row_basis_tracker_clear(&col_tracker_nmod);
                } else {
                    unified_row_basis_tracker_t col_tracker;
                    unified_row_basis_tracker_init(&col_tracker, ncols, current_size, &selection_ctx);
                    for (slong j = 0; j < ncols && col_tracker.current_rank < current_size; j++) {
                        slong col_idx = col_degrees[j].index;
                        unified_try_add_row_to_basis(&col_tracker, transposed, col_idx, current_size);
                    }
                    col_rank_selected = col_tracker.current_rank;
                    selected_cols = (slong*) flint_malloc(col_rank_selected * sizeof(slong));
                    memcpy(selected_cols, col_tracker.selected_indices,
                           col_rank_selected * sizeof(slong));
                    unified_row_basis_tracker_clear(&col_tracker);
                }
                dixon_maybe_print_step_detail_time("Step 3 column basis selection",
                                                   col_basis_cpu_start,
                                                   col_basis_wall_start);
            }

            if (iteration == 0) {
                current_col_indices = (slong*) flint_malloc(col_rank_selected * sizeof(slong));
            } else {
                flint_free(current_col_indices);
                current_col_indices = (slong*) flint_malloc(col_rank_selected * sizeof(slong));
            }
            memcpy(current_col_indices, selected_cols, col_rank_selected * sizeof(slong));
            current_size = FLINT_MIN(current_size, col_rank_selected);
            if (selected_cols != NULL) {
                flint_free(selected_cols);
            }

            for (slong i = 0; i < ncols * transposed_rows; i++) {
                field_clear_elem(&transposed[i], selection_ctx.field_id, ctx_ptr);
            }
            flint_free(transposed);
            flint_free(col_degrees);

            if (iteration > 0) {
                if (current_size == prev_size &&
                    current_size == FLINT_MIN(col_rank_selected, row_rank_selected) &&
                    indices_equal(current_row_indices, prev_row_indices, current_size) &&
                    indices_equal(current_col_indices, prev_col_indices, current_size)) {
                    converged = 1;
                }

                flint_free(prev_row_indices);
                flint_free(prev_col_indices);
                prev_row_indices = NULL;
                prev_col_indices = NULL;
            }

            iteration++;
        }

        free_entry_degree_cache(entry_degree_cache, nrows);

        clock_t iter_end = clock();
        (void) iter_start;
        (void) iter_end;

        slong final_rank;
        slong verify_rank = current_size;
        int final_ksy_ok = !require_ksy_precondition;
        int verify_ksy_ok = !require_ksy_precondition;
        clock_t verification_cpu_start = clock();
        double verification_wall_start = get_wall_time();
        if (selection_ctx.field_id == FIELD_ID_NMOD) {
            final_rank = evaluate_selected_submatrix_rank_nmod(
                full_matrix, current_row_indices, current_col_indices, current_size,
                param_vals, ctx);

            if (g_dixon_step3_second_verification) {
                fq_nmod_t *verify_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
                init_evaluation_parameters(verify_vals, npars, ctx, selection_attempt + 1);
                verify_rank = evaluate_selected_submatrix_rank_nmod(
                    full_matrix, current_row_indices, current_col_indices, current_size,
                    verify_vals, ctx);
                if (require_ksy_precondition) {
                    final_ksy_ok = selected_rows_satisfy_ksy_precondition_nmod(
                        full_matrix, current_row_indices, current_size, ncols,
                        ksy_constant_col, param_vals, ctx);
                    verify_ksy_ok = selected_rows_satisfy_ksy_precondition_nmod(
                        full_matrix, current_row_indices, current_size, ncols,
                        ksy_constant_col, verify_vals, ctx);
                }
                clear_evaluation_parameters(verify_vals, npars, ctx);
            } else if (require_ksy_precondition) {
                final_ksy_ok = selected_rows_satisfy_ksy_precondition_nmod(
                    full_matrix, current_row_indices, current_size, ncols,
                    ksy_constant_col, param_vals, ctx);
                verify_ksy_ok = final_ksy_ok;
            }
        } else if (use_extension_specialization) {
            final_rank = evaluate_selected_submatrix_rank_extension(
                full_matrix, current_row_indices, current_col_indices, current_size,
                param_vals, extension_eval_ctx);

            if (g_dixon_step3_second_verification) {
                fq_nmod_t *verify_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
                init_extension_evaluation_parameters(verify_vals, npars,
                                                     extension_eval_ctx,
                                                     selection_attempt + 1);
                verify_rank = evaluate_selected_submatrix_rank_extension(
                    full_matrix, current_row_indices, current_col_indices, current_size,
                    verify_vals, extension_eval_ctx);
                if (require_ksy_precondition) {
                    final_ksy_ok =
                        selected_rows_satisfy_ksy_precondition_extension(
                            full_matrix, current_row_indices, current_size, ncols,
                            ksy_constant_col, param_vals, extension_eval_ctx);
                    verify_ksy_ok =
                        selected_rows_satisfy_ksy_precondition_extension(
                            full_matrix, current_row_indices, current_size, ncols,
                            ksy_constant_col, verify_vals, extension_eval_ctx);
                }
                clear_evaluation_parameters(verify_vals, npars, extension_eval_ctx);
            } else if (require_ksy_precondition) {
                final_ksy_ok =
                    selected_rows_satisfy_ksy_precondition_extension(
                        full_matrix, current_row_indices, current_size, ncols,
                        ksy_constant_col, param_vals, extension_eval_ctx);
                verify_ksy_ok = final_ksy_ok;
            }
        } else {
            final_rank = evaluate_selected_submatrix_rank(
                full_matrix, current_row_indices, current_col_indices, current_size,
                param_vals, ctx);

            if (g_dixon_step3_second_verification) {
                fq_nmod_t *verify_vals = (fq_nmod_t*) flint_malloc(npars * sizeof(fq_nmod_t));
                init_evaluation_parameters(verify_vals, npars, ctx, selection_attempt + 1);
                verify_rank = evaluate_selected_submatrix_rank(
                    full_matrix, current_row_indices, current_col_indices, current_size,
                    verify_vals, ctx);
                if (require_ksy_precondition) {
                    final_ksy_ok = selected_rows_satisfy_ksy_precondition(
                        full_matrix, current_row_indices, current_size, ncols,
                        ksy_constant_col, param_vals, ctx);
                    verify_ksy_ok = selected_rows_satisfy_ksy_precondition(
                        full_matrix, current_row_indices, current_size, ncols,
                        ksy_constant_col, verify_vals, ctx);
                }
                clear_evaluation_parameters(verify_vals, npars, ctx);
            } else if (require_ksy_precondition) {
                final_ksy_ok = selected_rows_satisfy_ksy_precondition(
                    full_matrix, current_row_indices, current_size, ncols,
                    ksy_constant_col, param_vals, ctx);
                verify_ksy_ok = final_ksy_ok;
                verify_rank = final_rank;
            }
        }
        dixon_maybe_print_step_detail_time("Step 3 verification",
                                           verification_cpu_start,
                                           verification_wall_start);

        slong candidate_score = 0;
        if (final_rank == current_size) candidate_score++;
        if (g_dixon_step3_second_verification && verify_rank == current_size) candidate_score++;
        if (require_ksy_precondition) {
            if (final_ksy_ok) candidate_score += 2;
            if (g_dixon_step3_second_verification && verify_ksy_ok) candidate_score += 2;
        }

        if (current_size > accepted_size ||
            (current_size == accepted_size && candidate_score > accepted_score)) {
            if (accepted_rows) flint_free(accepted_rows);
            if (accepted_cols) flint_free(accepted_cols);
            accepted_rows = (slong*) flint_malloc(current_size * sizeof(slong));
            accepted_cols = (slong*) flint_malloc(current_size * sizeof(slong));
            memcpy(accepted_rows, current_row_indices, current_size * sizeof(slong));
            memcpy(accepted_cols, current_col_indices, current_size * sizeof(slong));
            accepted_size = current_size;
            accepted_score = candidate_score;
        }

        if (candidate_score == (require_ksy_precondition
                                ? (g_dixon_step3_second_verification ? 6 : 2)
                                : (g_dixon_step3_second_verification ? 2 : 1))) {
            selection_stable = 1;
        } else {
            if (require_ksy_precondition && (!final_ksy_ok || !verify_ksy_ok)) {
                dixon_debug_log("  Candidate submatrix failed the KSY precondition; retrying with a new specialization.\n");
            } else {
                dixon_debug_log("  Submatrix verification failed; retrying with a new specialization.\n");
            }
        }

        if (prev_row_indices) flint_free(prev_row_indices);
        if (prev_col_indices) flint_free(prev_col_indices);
        flint_free(current_row_indices);
        flint_free(current_col_indices);
        if (param_vals) {
            clear_evaluation_parameters(param_vals, npars,
                                       use_extension_specialization ? extension_eval_ctx : ctx);
        }
        for (slong i = 0; i < nrows * ncols; i++) {
            field_clear_elem(&unified_mat[i], selection_ctx.field_id, ctx_ptr);
        }
        flint_free(unified_mat);
    }

    if (!accepted_rows || !accepted_cols) {
        *row_indices_out = NULL;
        *col_indices_out = NULL;
        *num_rows = 0;
        *num_cols = 0;
    } else {
        if (!selection_stable) {
            if (require_ksy_precondition) {
                dixon_info_log("Warning: using best available submatrix after %ld attempts (KSY precondition not certified for every test).\n",
                               MAX_SELECTION_ATTEMPTS);
            } else {
                dixon_info_log("Warning: using best available submatrix after %ld attempts.\n",
                               MAX_SELECTION_ATTEMPTS);
            }
        }
        *row_indices_out = accepted_rows;
        *col_indices_out = accepted_cols;
        *num_rows = accepted_size;
        *num_cols = accepted_size;
    }

    field_ctx_clear(&selection_ctx);
    if (use_extension_specialization) {
        fq_nmod_ctx_clear(extension_eval_ctx);
    }
}

/* Assign IDs in first-occurrence order. Only unique supports consume storage;
 * exponent pointers are borrowed until collection is complete. */
static slong dixon_intern_monom(monom_t **monoms, slong *count, slong *capacity,
                                hash_entry_t ***index, slong *hash_size,
                                const slong *exp, slong nvars)
{
    ulong hash = hash_monom_exponents(exp, nvars);
    slong bucket = hash & (*hash_size - 1);
    for (hash_entry_t *e = (*index)[bucket]; e; e = e->next)
        if (memcmp(e->exp, exp, (size_t) nvars * sizeof(slong)) == 0)
            return e->idx;

    if (*count >= *hash_size / 2) {
        slong new_size = *hash_size * 2;
        hash_entry_t **grown = flint_calloc((size_t) new_size, sizeof(*grown));
        for (slong i = 0; i < *hash_size; i++) {
            hash_entry_t *e = (*index)[i];
            while (e) {
                hash_entry_t *next = e->next;
                slong h = hash_monom_exponents(e->exp, nvars) & (new_size - 1);
                e->next = grown[h];
                grown[h] = e;
                e = next;
            }
        }
        flint_free(*index);
        *index = grown;
        *hash_size = new_size;
        bucket = hash & (new_size - 1);
    }
    if (*count == *capacity) {
        *capacity = *capacity ? 2 * *capacity : 16;
        *monoms = flint_realloc(*monoms, (size_t) *capacity * sizeof(**monoms));
    }
    slong id = (*count)++;
    (*monoms)[id].exp = (slong *) exp;
    (*monoms)[id].idx = id;
    hash_entry_t *e = flint_malloc(sizeof(*e));
    e->exp = (slong *) exp;
    e->idx = id;
    e->next = (*index)[bucket];
    (*index)[bucket] = e;
    return id;
}

/* Keep the existing single-free ownership convention for returned monomials.
 * Rebind the retained hash nodes to this compact, owned exponent storage. */
static void dixon_pack_monoms(monom_t **monoms, slong count, slong nvars,
                              hash_entry_t **index, slong hash_size)
{
    if (count == 0) {
        flint_free(*monoms);
        *monoms = NULL;
        return;
    }
    monom_t *packed = flint_malloc((size_t) count *
                                  (sizeof(monom_t) + (size_t) nvars * sizeof(slong)));
    slong *exps = (slong *) (packed + count);
    for (slong i = 0; i < count; i++) {
        packed[i].idx = i;
        packed[i].exp = exps + i * nvars;
        memcpy(packed[i].exp, (*monoms)[i].exp, (size_t) nvars * sizeof(slong));
    }
    for (slong i = 0; i < hash_size; i++)
        for (hash_entry_t *e = index[i]; e; e = e->next)
            e->exp = packed[e->idx].exp;
    flint_free(*monoms);
    *monoms = packed;
}

static void collect_unique_monomials(
    monom_t **x_monoms_out, slong *nx_monoms_out,
    monom_t **dual_monoms_out, slong *ndual_monoms_out,
    hash_entry_t ***x_index, slong *x_hash_size,
    hash_entry_t ***dual_index, slong *dual_hash_size,
    slong *term_rows, slong *term_cols,
    const fq_mvpoly_t *dixon_poly,
    const slong *d0, const slong *d1, slong nvars)
{
    slong x_capacity = 0, dual_capacity = 0;
    *x_monoms_out = NULL;
    *dual_monoms_out = NULL;
    *nx_monoms_out = *ndual_monoms_out = 0;
    *x_hash_size = *dual_hash_size = 16;
    *x_index = flint_calloc(16, sizeof(hash_entry_t *));
    *dual_index = flint_calloc(16, sizeof(hash_entry_t *));
    for (slong t = 0; t < dixon_poly->nterms; t++) {
        const slong *exp = dixon_poly->terms[t].var_exp;
        term_rows[t] = term_cols[t] = -1;
        if (!exp) continue;
        int valid = 1;
        for (slong k = 0; k < nvars; k++)
            if (exp[k] >= d0[k] || exp[nvars + k] >= d1[k]) {
                valid = 0;
                break;
            }
        if (!valid) continue;
        term_rows[t] = dixon_intern_monom(x_monoms_out, nx_monoms_out,
                            &x_capacity, x_index, x_hash_size, exp, nvars);
        term_cols[t] = dixon_intern_monom(dual_monoms_out, ndual_monoms_out,
                            &dual_capacity, dual_index, dual_hash_size, exp + nvars, nvars);
    }
    dixon_pack_monoms(x_monoms_out, *nx_monoms_out, nvars, *x_index, *x_hash_size);
    dixon_pack_monoms(dual_monoms_out, *ndual_monoms_out, nvars, *dual_index, *dual_hash_size);
}

// Allocate single element on demand
fq_mvpoly_t* get_matrix_entry_lazy(fq_mvpoly_t ***matrix, slong i, slong j,
                                  slong npars, const fq_nmod_ctx_t ctx) {
    if (!matrix[i][j]) {
        matrix[i][j] = (fq_mvpoly_t*) flint_malloc(sizeof(fq_mvpoly_t));
        fq_mvpoly_init(matrix[i][j], 0, npars, ctx);
    }
    return matrix[i][j];
}
static void fill_coefficient_matrix_optimized(fq_mvpoly_t ***full_matrix,
                                      slong nx_monoms, slong ndual_monoms,
                                      const fq_mvpoly_t *dixon_poly, slong npars,
                                      const slong *term_rows, const slong *term_cols)
{
    if (dixon_poly->nterms <= 0) return;
    slong *row_offsets = flint_calloc((size_t) nx_monoms + 1, sizeof(slong));
    slong *row_write_offsets = flint_malloc((size_t) nx_monoms * sizeof(slong));
    for (slong t = 0; t < dixon_poly->nterms; t++)
        if (term_rows[t] >= 0) row_offsets[term_rows[t] + 1]++;
    for (slong row = 0; row < nx_monoms; row++)
        row_offsets[row + 1] += row_offsets[row];
    memcpy(row_write_offsets, row_offsets, (size_t) nx_monoms * sizeof(slong));
    slong *row_term_indices = flint_malloc((size_t) row_offsets[nx_monoms] * sizeof(slong));
    for (slong t = 0; t < dixon_poly->nterms; t++)
        if (term_rows[t] >= 0) row_term_indices[row_write_offsets[term_rows[t]]++] = t;

#ifdef _OPENMP
    #pragma omp parallel if(nx_monoms > 1)
#endif
    {
        slong *entry_counts = flint_calloc((size_t) ndual_monoms, sizeof(slong));
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 1)
#endif
        for (slong row = 0; row < nx_monoms; row++) {
            for (slong idx = row_offsets[row]; idx < row_offsets[row + 1]; idx++)
                entry_counts[term_cols[row_term_indices[idx]]]++;
            for (slong col = 0; col < ndual_monoms; col++) {
                slong count = entry_counts[col];
                if (!count) continue;
                fq_mvpoly_t *entry = flint_malloc(sizeof(*entry));
                entry->nvars = 0;
                entry->npars = npars;
                entry->nterms = 0;
                entry->alloc = count;
                entry->ctx = dixon_poly->ctx;
                /* add_term_fast initializes every occupied slot. */
                entry->terms = flint_malloc((size_t) count * sizeof(*entry->terms));
                full_matrix[row][col] = entry;
                entry_counts[col] = 0;
            }
            for (slong idx = row_offsets[row]; idx < row_offsets[row + 1]; idx++) {
                slong t = row_term_indices[idx];
                fq_mvpoly_add_term_fast(full_matrix[row][term_cols[t]], NULL,
                                        dixon_poly->terms[t].par_exp,
                                        dixon_poly->terms[t].coeff);
            }
        }
        flint_free(entry_counts);
    }
    flint_free(row_term_indices);
    flint_free(row_write_offsets);
    flint_free(row_offsets);
}
// Optimized version of find_fq_optimal_maximal_rank_submatrix
// ============ Extract coefficient matrix ============

/* The determinant degree is bounded by either sum of row maxima or sum of
   column maxima.  Only inspect the selected minor, never the full matrix. */
static slong dixon_minor_degree_bound(fq_mvpoly_t ***matrix,
                                     slong *rows, slong *cols, slong size)
{
    slong row_sum = 0, col_sum = 0;
    for (slong i = 0; i < size; i++) {
        row_sum += FLINT_MAX(0, compute_fq_selected_cols_row_max_total_degree(
                                   matrix, rows[i], cols, size, 1));
        col_sum += FLINT_MAX(0, compute_fq_selected_rows_col_max_total_degree(
                                   matrix, rows, size, cols[i], 1));
    }
    return FLINT_MIN(row_sum, col_sum);
}

/* Sparse, per-specialization value cache shared by Schur completion and basis
   exchange.  Zero values are cached too; untouched matrix entries cost no space. */
typedef struct {
    size_t key;
    mp_limb_t value;
} dixon_eval_slot_t;

typedef struct {
    fq_mvpoly_t ***matrix;
    slong ncols;
    fq_nmod_t *params;
    fq_nmod_t scratch;
    const fq_nmod_ctx_struct *ctx;
    dixon_eval_slot_t *slots;
    size_t capacity, count;
} dixon_eval_cache_t;

static size_t dixon_eval_hash(size_t key, size_t capacity)
{
    key ^= key >> 16;
    key *= (size_t) 0x45d9f3b;
    key ^= key >> 16;
    return key & (capacity - 1);
}

static void dixon_eval_cache_init(dixon_eval_cache_t *cache,
                                  fq_mvpoly_t ***matrix, slong ncols,
                                  fq_nmod_t *params, const fq_nmod_ctx_t ctx)
{
    cache->matrix = matrix;
    cache->ncols = ncols;
    cache->params = params;
    cache->ctx = ctx;
    cache->capacity = 1024;
    cache->count = 0;
    cache->slots = flint_calloc(cache->capacity, sizeof(*cache->slots));
    fq_nmod_init(cache->scratch, ctx);
}

static mp_limb_t dixon_eval_cached(dixon_eval_cache_t *cache, slong row, slong col)
{
    size_t key = (size_t) row * cache->ncols + col + 1;
    size_t pos = dixon_eval_hash(key, cache->capacity);
    while (cache->slots[pos].key && cache->slots[pos].key != key)
        pos = (pos + 1) & (cache->capacity - 1);
    if (cache->slots[pos].key) return cache->slots[pos].value;
    /* At most 16 MiB of slots; large streamed sweeps bypass this cache. */
    const size_t max_capacity = ((size_t) 16 << 20) / sizeof(*cache->slots);
    if (cache->count >= cache->capacity / 2 && cache->capacity < max_capacity) {
        size_t old_capacity = cache->capacity;
        dixon_eval_slot_t *old = cache->slots;
        cache->capacity *= 2;
        cache->slots = flint_calloc(cache->capacity, sizeof(*cache->slots));
        for (size_t i = 0; i < old_capacity; i++) if (old[i].key) {
            size_t p = dixon_eval_hash(old[i].key, cache->capacity);
            while (cache->slots[p].key) p = (p + 1) & (cache->capacity - 1);
            cache->slots[p] = old[i];
        }
        flint_free(old);
        pos = dixon_eval_hash(key, cache->capacity);
        while (cache->slots[pos].key) pos = (pos + 1) & (cache->capacity - 1);
    }
    fq_mvpoly_t *entry = cache->matrix[row][col];
    mp_limb_t value = 0;
    if (entry && entry->nterms) {
        evaluate_fq_mvpoly_at_params(cache->scratch, entry, cache->params);
        value = nmod_poly_get_coeff_ui(cache->scratch, 0);
    }
    if (cache->count < cache->capacity / 2) {
        cache->slots[pos].key = key;
        cache->slots[pos].value = value;
        cache->count++;
    }
    return value;
}

static void dixon_eval_cache_clear(dixon_eval_cache_t *cache)
{
    fq_nmod_clear(cache->scratch, cache->ctx);
    flint_free(cache->slots);
}

/* H = E_k ... E_1 H0 F_1 ... F_l.  Each E replaces one row, each F one
   column.  Store inverse updates I + e_p h (row) or I + h e_p^T (column).
   After r exchanges, rebase on the current r x r minor to bound storage and
   the cost of applying the update chain. */
typedef struct {
    nmod_mat_t lu, transposed, rhs, tmp, answer, updates;
    slong *perm, *positions;
    unsigned char *columns;
    slong size, count, rebases;
} dixon_exchange_solver_t;

static void dixon_exchange_solver_init(dixon_exchange_solver_t *solver,
                                       slong size, mp_limb_t prime)
{
    solver->size = size;
    solver->count = solver->rebases = 0;
    nmod_mat_init(solver->lu, size, size, prime);
    nmod_mat_init(solver->transposed, size, size, prime);
    nmod_mat_init(solver->rhs, size, 1, prime);
    nmod_mat_init(solver->tmp, size, 1, prime);
    nmod_mat_init(solver->answer, size, 1, prime);
    nmod_mat_init(solver->updates, size, size, prime);
    solver->perm = flint_malloc((size_t) size * sizeof(slong));
    solver->positions = flint_malloc((size_t) size * sizeof(slong));
    solver->columns = flint_malloc((size_t) size);
}

static void dixon_exchange_solver_clear(dixon_exchange_solver_t *solver)
{
    flint_free(solver->columns);
    flint_free(solver->positions);
    flint_free(solver->perm);
    nmod_mat_clear(solver->updates);
    nmod_mat_clear(solver->answer);
    nmod_mat_clear(solver->tmp);
    nmod_mat_clear(solver->rhs);
    nmod_mat_clear(solver->transposed);
    nmod_mat_clear(solver->lu);
}

/* Assemble a full LU from the existing core LU and only a delta x delta LU.
   [A B; V W] = [L 0; V U^-1 Ls] [U L^-1 B; 0 Us], with the Schur row
   permutation applied to the bottom rows.  No new size x size factorization. */
static void dixon_exchange_seed(dixon_exchange_solver_t *solver,
                                const nmod_mat_t lower, const nmod_mat_t upper,
                                const nmod_mat_t x, const nmod_mat_t v,
                                const nmod_mat_t schur,
                                const slong *chosen_r, const slong *chosen_c)
{
    slong s = upper->r, delta = solver->size - s;
    mp_limb_t prime = upper->mod.n;
    nmod_mat_t small, ut, vt, bottom, xs, top;
    nmod_mat_init(small, delta, delta, prime);
    nmod_mat_init(ut, s, s, prime);
    nmod_mat_init(vt, s, delta, prime);
    nmod_mat_init(bottom, s, delta, prime);
    nmod_mat_init(xs, s, delta, prime);
    nmod_mat_init(top, s, delta, prime);
    slong *sp = flint_malloc((size_t) delta * sizeof(slong));
    for (slong i = 0; i < delta; i++) {
        for (slong j = 0; j < delta; j++)
            nmod_mat_entry(small, i, j) = nmod_mat_entry(schur, chosen_r[i], chosen_c[j]);
        for (slong j = 0; j < s; j++) {
            nmod_mat_entry(vt, j, i) = nmod_mat_entry(v, chosen_r[i], j);
            nmod_mat_entry(xs, j, i) = nmod_mat_entry(x, j, chosen_c[i]);
        }
    }
    if (nmod_mat_lu(sp, small, 0) != delta)
        flint_throw(FLINT_ERROR, "Dixon Schur seed lost rank\n");
    nmod_mat_transpose(ut, upper);
    nmod_mat_solve_tril(bottom, ut, vt, 0);
    nmod_mat_mul(top, upper, xs);
    for (slong i = 0; i < s; i++) {
        solver->perm[i] = i;
        for (slong j = 0; j < s; j++)
            nmod_mat_entry(solver->lu, i, j) = i > j
                ? nmod_mat_entry(lower, i, j) : nmod_mat_entry(upper, i, j);
        for (slong j = 0; j < delta; j++)
            nmod_mat_entry(solver->lu, i, s + j) = nmod_mat_entry(top, i, j);
    }
    for (slong i = 0; i < delta; i++) {
        solver->perm[s + i] = s + sp[i];
        for (slong j = 0; j < s; j++)
            nmod_mat_entry(solver->lu, s + i, j) = nmod_mat_entry(bottom, j, sp[i]);
        for (slong j = 0; j < delta; j++)
            nmod_mat_entry(solver->lu, s + i, s + j) = nmod_mat_entry(small, i, j);
    }
    nmod_mat_transpose(solver->transposed, solver->lu);
    flint_free(sp);
    nmod_mat_clear(top); nmod_mat_clear(xs); nmod_mat_clear(bottom);
    nmod_mat_clear(vt); nmod_mat_clear(ut); nmod_mat_clear(small);
}

static void dixon_exchange_apply(dixon_exchange_solver_t *solver,
                                 mp_limb_t *vector, slong update, int scatter)
{
    slong size = solver->size, p = solver->positions[update];
    nmod_t mod = solver->lu->mod;
    if (scatter) {
        mp_limb_t factor = vector[p];
        for (slong i = 0; i < size; i++)
            vector[i] = nmod_add(vector[i], nmod_mul(factor,
                                nmod_mat_entry(solver->updates, update, i), mod), mod);
    } else {
        mp_limb_t sum = 0;
        for (slong i = 0; i < size; i++)
            sum = nmod_add(sum, nmod_mul(vector[i],
                           nmod_mat_entry(solver->updates, update, i), mod), mod);
        vector[p] = nmod_add(vector[p], sum, mod);
    }
}

/* columns=0 solves lambda H = a; columns=1 solves H lambda = a. */
static void dixon_exchange_solve(dixon_exchange_solver_t *solver,
                                 mp_limb_t *vector, int columns)
{
    for (slong k = solver->count; k-- > 0;)
        if (solver->columns[k] != columns)
            dixon_exchange_apply(solver, vector, k, 0);
    for (slong i = 0; i < solver->size; i++)
        nmod_mat_entry(solver->rhs, i, 0) = vector[columns ? solver->perm[i] : i];
    if (columns) {
        nmod_mat_solve_tril(solver->tmp, solver->lu, solver->rhs, 1);
        nmod_mat_solve_triu(solver->answer, solver->lu, solver->tmp, 0);
    } else {
        nmod_mat_solve_tril(solver->tmp, solver->transposed, solver->rhs, 0);
        nmod_mat_solve_triu(solver->answer, solver->transposed, solver->tmp, 1);
    }
    for (slong i = 0; i < solver->size; i++)
        vector[columns ? i : solver->perm[i]] = nmod_mat_entry(solver->answer, i, 0);
    for (slong k = 0; k < solver->count; k++)
        if (solver->columns[k] == columns)
            dixon_exchange_apply(solver, vector, k, 1);
}

static void dixon_exchange_record(dixon_exchange_solver_t *solver,
                                  const mp_limb_t *coeff, slong position, int columns,
                                  dixon_eval_cache_t *cache, slong *rows, slong *cols)
{
    slong k = solver->count++, size = solver->size;
    nmod_t mod = solver->lu->mod;
    mp_limb_t inv = nmod_inv(coeff[position], mod);
    solver->positions[k] = position;
    solver->columns[k] = columns;
    for (slong i = 0; i < size; i++)
        nmod_mat_entry(solver->updates, k, i) =
            nmod_mul(nmod_sub(i == position, coeff[i], mod), inv, mod);
    if (solver->count == size) {
        for (slong i = 0; i < size; i++)
            for (slong j = 0; j < size; j++)
                nmod_mat_entry(solver->lu, i, j) = dixon_eval_cached(cache, rows[i], cols[j]);
        if (nmod_mat_lu(solver->perm, solver->lu, 0) != size)
            flint_throw(FLINT_ERROR, "Dixon basis exchange lost rank\n");
        nmod_mat_transpose(solver->transposed, solver->lu);
        solver->count = 0;
        solver->rebases++;
    }
}

/* Incremental minimum-weight basis selection over ALL candidates.  Weights
   are exactly the original selector's maximum parameter degree on the fixed
   opposite set, with the same index tie-break.  Circuit exchange removes the
   worst eligible basis element.  Unlike a determinant-bound filter, this also
   permits equal-degree exchanges needed to reproduce the greedy basis. */
static slong dixon_exchange_axis(dixon_exchange_solver_t *solver,
                                 dixon_eval_cache_t *cache, slong nrows, slong ncols,
                                 slong *rows, slong *cols, int columns)
{
    slong size = solver->size, count = columns ? ncols : nrows, exchanges = 0;
    slong *indices = columns ? cols : rows;
    fq_index_degree_pair *weights = flint_malloc((size_t) count * sizeof(*weights));
    fq_index_degree_pair *order = flint_malloc((size_t) count * sizeof(*order));
    slong *present = flint_malloc((size_t) count * sizeof(*present));
    mp_limb_t *coeff = flint_malloc((size_t) size * sizeof(*coeff));
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(size >= 256)
#endif
    for (slong i = 0; i < count; i++) {
        weights[i].index = i;
        weights[i].degree = columns
            ? compute_fq_selected_rows_col_max_total_degree(cache->matrix, rows, size, i, 1)
            : compute_fq_selected_cols_row_max_total_degree(cache->matrix, i, cols, size, 1);
        present[i] = -1;
    }
    memcpy(order, weights, (size_t) count * sizeof(*order));
    qsort(order, (size_t) count, sizeof(*order), compare_fq_degrees);
    for (slong i = 0; i < size; i++) present[indices[i]] = i;
    for (slong k = 0; k < count; k++) {
        slong candidate = order[k].index, worst = 0;
        if (present[candidate] >= 0 || order[k].degree < 0) continue;
        for (slong j = 1; j < size; j++)
            if (compare_fq_degrees(&weights[indices[j]], &weights[indices[worst]]) > 0) worst = j;
        if (compare_fq_degrees(&weights[candidate], &weights[indices[worst]]) >= 0) break;
        for (slong j = 0; j < size; j++)
            coeff[j] = columns ? dixon_eval_cached(cache, rows[j], candidate)
                               : dixon_eval_cached(cache, candidate, cols[j]);
        dixon_exchange_solve(solver, coeff, columns);
        worst = -1;
        for (slong j = 0; j < size; j++) if (coeff[j] &&
            (worst < 0 || compare_fq_degrees(&weights[indices[j]], &weights[indices[worst]]) > 0))
            worst = j;
        if (worst < 0 || compare_fq_degrees(&weights[candidate], &weights[indices[worst]]) >= 0)
            continue;
        present[indices[worst]] = -1;
        indices[worst] = candidate;
        present[candidate] = worst;
        dixon_exchange_record(solver, coeff, worst, columns, cache, rows, cols);
        exchanges++;
    }
    flint_free(coeff); flint_free(present); flint_free(order); flint_free(weights);
    return exchanges;
}

static void dixon_refine_schur_minor(dixon_eval_cache_t *cache,
                                    slong nrows, slong ncols, slong *rows, slong *cols,
                                    const nmod_mat_t lower, const nmod_mat_t upper,
                                    const nmod_mat_t x, const nmod_mat_t v,
                                    const nmod_mat_t schur, slong delta,
                                    const slong *chosen_r, const slong *chosen_c)
{
    dixon_exchange_solver_t solver;
    slong size = upper->r + delta;
    slong before = dixon_minor_degree_bound(cache->matrix, rows, cols, size);
    clock_t cpu_start = clock();
    double wall_start = get_wall_time();
    dixon_exchange_solver_init(&solver, size, upper->mod.n);
    dixon_exchange_seed(&solver, lower, upper, x, v, schur, chosen_r, chosen_c);
    for (int pass = 0; pass < 3; pass++) {
        slong nr = dixon_exchange_axis(&solver, cache, nrows, ncols, rows, cols, 0);
        slong nc = dixon_exchange_axis(&solver, cache, nrows, ncols, rows, cols, 1);
        dixon_debug_log("  Degree-aware exchange pass %d: rows=%ld, columns=%ld\n", pass + 1, nr, nc);
        if (!nr && !nc) break;
    }
    dixon_debug_log("  Degree-aware exchange: degree bound %ld -> %ld, rebases=%ld\n",
                    before, dixon_minor_degree_bound(cache->matrix, rows, cols, size), solver.rebases);
    dixon_exchange_solver_clear(&solver);
    dixon_maybe_print_step_detail_time("Step 3 degree-aware basis exchange", cpu_start, wall_start);
}

/* Large minors use streamed greedy selection, avoiding one solve per
   candidate and the long elementary-update chain.  Each panel contains
   candidates in EXACT degree/index order as columns.  LU scans these columns
   from left to right, so its pivot columns are the original greedy basis.
   Only L and its row permutation survive between panels: the old U is not
   needed to reduce subsequent candidates modulo the accepted span. */
#define DIXON_DEGREE_BLOCK_THRESHOLD 256
#define DIXON_DEGREE_PANEL_SIZE 128

/* Few external directions favour reusing the Schur-seeded LU, even when
 * the minor itself is large. Streaming pays for building the entire basis
 * again on every axis/pass; incremental exchange only solves outsiders. */
static int dixon_use_schur_exchange(slong size, slong nrows, slong ncols)
{
    return size < DIXON_DEGREE_BLOCK_THRESHOLD ||
           (nrows - size <= 32 && ncols - size <= 32);
}

static void dixon_permute_basis_prefix(nmod_mat_t lower, slong *perm,
                                       slong rank, const slong *action)
{
    slong remaining = lower->r - rank;
    unsigned char *seen = flint_calloc((size_t) remaining, 1);
    mp_limb_t *saved = flint_malloc((size_t) FLINT_MAX(rank, 1) * sizeof(*saved));
    for (slong i = 0; i < remaining; i++) if (!seen[i]) {
        slong j = i, saved_perm = perm[rank + i];
        if (rank) memcpy(saved, nmod_mat_entry_ptr(lower, rank + i, 0), (size_t) rank * sizeof(*saved));
        while (action[j] != i) {
            slong src = action[j];
            if (rank) memcpy(nmod_mat_entry_ptr(lower, rank + j, 0),
                             nmod_mat_entry_ptr(lower, rank + src, 0), (size_t) rank * sizeof(*saved));
            perm[rank + j] = perm[rank + src];
            seen[j] = 1;
            j = src;
        }
        if (rank) memcpy(nmod_mat_entry_ptr(lower, rank + j, 0), saved, (size_t) rank * sizeof(*saved));
        perm[rank + j] = saved_perm;
        seen[j] = 1;
    }
    flint_free(saved);
    flint_free(seen);
}

static slong dixon_degree_stream_axis(fq_mvpoly_t ***matrix,
                                      slong nrows, slong ncols, slong size,
                                      slong *rows, slong *cols, int columns,
                                      fq_nmod_t *params, const fq_nmod_ctx_t ctx,
                                      slong panel_size)
{
    slong count = columns ? ncols : nrows;
    slong *indices = columns ? cols : rows;
    slong *opposite = columns ? rows : cols;
    fq_index_degree_pair *order = flint_malloc((size_t) count * sizeof(*order));
    unsigned char *present = flint_calloc((size_t) count, 1);
    slong *selected = flint_malloc((size_t) size * sizeof(*selected));
    slong *perm = flint_malloc((size_t) size * sizeof(*perm));
    slong *action = flint_malloc((size_t) size * sizeof(*action));
    double start = get_wall_time(), last_report = start;
    dixon_debug_log("  Degree-aware blocked %s: scanning %ld candidates, target=%ld, panel=%ld\n",
                    columns ? "columns" : "rows", count, size, panel_size);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (slong i = 0; i < count; i++) {
        order[i].index = i;
        order[i].degree = columns
            ? compute_fq_selected_rows_col_max_total_degree(matrix, rows, size, i, 1)
            : compute_fq_selected_cols_row_max_total_degree(matrix, i, cols, size, 1);
    }
    qsort(order, (size_t) count, sizeof(*order), compare_fq_degrees);
    for (slong i = 0; i < size; i++) { perm[i] = i; present[indices[i]] = 1; }
    nmod_mat_t lower;
    nmod_mat_init(lower, size, size, fq_nmod_ctx_prime(ctx));
    slong rank = 0, changed = 0, processed = 0;
    while (processed < count && order[processed].degree < 0) processed++;
    while (processed < count && rank < size) {
        slong width = FLINT_MIN(panel_size, count - processed);
        nmod_mat_t panel, residual;
        nmod_mat_init(panel, size, width, fq_nmod_ctx_prime(ctx));
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            fq_nmod_t value;
            fq_nmod_init(value, ctx);
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            for (slong i = 0; i < size; i++) {
                slong fixed = opposite[perm[i]];
                for (slong j = 0; j < width; j++) {
                    slong candidate = order[processed + j].index;
                    fq_mvpoly_t *entry = columns ? matrix[fixed][candidate] : matrix[candidate][fixed];
                    if (entry && entry->nterms) {
                        evaluate_fq_mvpoly_at_nmod_params_direct(value, entry, params, ctx);
                        nmod_mat_entry(panel, i, j) = nmod_poly_get_coeff_ui(value, 0);
                    }
                }
            }
            fq_nmod_clear(value, ctx);
        }
        nmod_mat_window_init(residual, panel, rank, 0, size, width);
        if (rank) {
            nmod_mat_t top, l11, l21, solved, product;
            nmod_mat_window_init(top, panel, 0, 0, rank, width);
            nmod_mat_window_init(l11, lower, 0, 0, rank, rank);
            nmod_mat_window_init(l21, lower, rank, 0, size, rank);
            nmod_mat_init(solved, rank, width, panel->mod.n);
            nmod_mat_init(product, size - rank, width, panel->mod.n);
            nmod_mat_solve_tril(solved, l11, top, 1);
            nmod_mat_mul(product, l21, solved);
            nmod_mat_sub(residual, residual, product);
            nmod_mat_clear(product); nmod_mat_clear(solved);
            nmod_mat_window_clear(l21); nmod_mat_window_clear(l11); nmod_mat_window_clear(top);
        }
        slong added = nmod_mat_lu(action, residual, 0), pivot = 0;
        if (added) {
            for (slong i = 0; i < added; i++) {
                while (pivot < width && !nmod_mat_entry(residual, i, pivot)) pivot++;
                if (pivot == width) flint_throw(FLINT_ERROR, "Dixon blocked selector: missing pivot\n");
                selected[rank + i] = order[processed + pivot++].index;
            }
            dixon_permute_basis_prefix(lower, perm, rank, action);
            for (slong i = 0; i < size - rank; i++)
                for (slong j = 0; j < FLINT_MIN(i, added); j++)
                    nmod_mat_entry(lower, rank + i, rank + j) = nmod_mat_entry(residual, i, j);
            rank += added;
        }
        nmod_mat_window_clear(residual);
        nmod_mat_clear(panel);
        processed += width;
        double now = get_wall_time();
        if (rank == size || now - last_report >= 1.0) {
            dixon_debug_log("    Degree-aware blocked %s: candidates=%ld/%ld, rank=%ld/%ld, wall=%.3fs\n",
                            columns ? "columns" : "rows", processed, count, rank, size, now - start);
            last_report = now;
        }
    }
    if (rank != size) flint_throw(FLINT_ERROR, "Dixon blocked selector lost the certified basis\n");
    for (slong i = 0; i < size; i++) if (!present[selected[i]]) changed++;
    memcpy(indices, selected, (size_t) size * sizeof(*indices));
    nmod_mat_clear(lower);
    flint_free(action); flint_free(perm); flint_free(selected); flint_free(present); flint_free(order);
    return changed;
}

static void dixon_refine_streamed_minor(fq_mvpoly_t ***matrix,
                                       slong nrows, slong ncols, slong size,
                                       slong *rows, slong *cols,
                                       fq_nmod_t *params, const fq_nmod_ctx_t ctx)
{
    clock_t cpu_start = clock();
    double wall_start = get_wall_time();
    dixon_debug_log("  Degree-aware refinement backend: blocked streaming (rank=%ld)\n", size);
    for (int pass = 0; pass < 3; pass++) {
        slong nr = dixon_degree_stream_axis(matrix, nrows, ncols, size, rows, cols, 0,
                                            params, ctx, DIXON_DEGREE_PANEL_SIZE);
        slong nc = dixon_degree_stream_axis(matrix, nrows, ncols, size, rows, cols, 1,
                                            params, ctx, DIXON_DEGREE_PANEL_SIZE);
        dixon_debug_log("  Degree-aware blocked pass %d: rows=%ld, columns=%ld\n", pass + 1, nr, nc);
        if (!nr && !nc) break;
    }
    dixon_maybe_print_step_detail_time("Step 3 blocked degree-aware selection", cpu_start, wall_start);
}

/* Reuse the rank-deficient candidate's packed LU.  The core consists of the
   first s permuted rows and U's pivot columns.  Its L is stored below the
   first s diagonal entries (not below the original pivot columns).
   Grow the reserve rectangle geometrically, caching V, A^-1 B and S.  Only
   the small Schur rectangle is refactored after growth.  A bounded unsuccessful
   search leaves the caller's indices untouched for the general fallback. */
static int dixon_repair_predicted_minor(fq_mvpoly_t ***matrix,
                                       slong nrows, slong ncols,
                                       slong *rows, slong *cols, slong size,
                                       const nmod_mat_t lu, const slong *perm, slong s,
                                       const fq_index_degree_pair *row_order,
                                       const fq_index_degree_pair *col_order,
                                       slong sigma, fq_nmod_t *params, const fq_nmod_ctx_t ctx)
{
    slong delta = size - s;
    if (s <= 0 || delta <= 0 || delta > 128) return 0;
    slong budget = FLINT_MIN(32 * delta, 128);
    slong max_r = FLINT_MIN(nrows - s, budget);
    slong max_c = FLINT_MIN(ncols - s, budget);
    if (max_r < delta || max_c < delta) return 0;
    int success = 0;
    int use_exchange = dixon_use_schur_exchange(size, nrows, ncols);
    slong *r0 = flint_malloc((size_t) s * sizeof(slong));
    slong *c0 = flint_malloc((size_t) s * sizeof(slong));
    slong *pivots = flint_malloc((size_t) s * sizeof(slong));
    slong *rx = flint_malloc((size_t) max_r * sizeof(slong));
    slong *cx = flint_malloc((size_t) max_c * sizeof(slong));
    unsigned char *used_r = flint_calloc((size_t) nrows, 1);
    unsigned char *used_c = flint_calloc((size_t) ncols, 1);
    slong next = 0;
    for (slong i = 0; i < s; i++) {
        r0[i] = rows[perm[i]];
        used_r[r0[i]] = 1;
        while (next < size && nmod_mat_entry(lu, i, next) == 0) next++;
        if (next == size) goto cleanup_indices;
        pivots[i] = next;
        c0[i] = cols[next++];
        used_c[c0[i]] = 1;
    }
    /* Admit complementary degree pairs together.  Two sorted cursors avoid
       scanning the full Cartesian product of unused rows and columns. */
    slong kr = 0, kc = 0, ri = 0, ci = ncols - 1;
    while (sigma >= 0 && ri < nrows && ci >= 0 && kr < max_r && kc < max_c) {
        if (used_r[row_order[ri].index]) { ri++; continue; }
        if (used_c[col_order[ci].index]) { ci--; continue; }
        slong sum = row_order[ri].degree + col_order[ci].degree;
        if (sum < sigma) { ri++; continue; }
        if (sum > sigma) { ci--; continue; }
        rx[kr++] = row_order[ri].index;
        cx[kc++] = col_order[ci].index;
        used_r[row_order[ri++].index] = 1;
        used_c[col_order[ci--].index] = 1;
    }
    /* Also retain dependent candidate directions: external rows or columns
       can make them useful.  Fill unmatched slots in ascending degree order. */
    for (slong i = 0; i < nrows && kr < max_r; i++)
        if (!used_r[row_order[i].index]) rx[kr++] = row_order[i].index;
    for (slong i = 0; i < ncols && kc < max_c; i++)
        if (!used_c[col_order[i].index]) cx[kc++] = col_order[i].index;

    mp_limb_t prime = fq_nmod_ctx_prime(ctx);
    nmod_mat_t lower, upper, x, v, schur;
    nmod_mat_init(lower, s, s, prime);
    nmod_mat_init(upper, s, s, prime);
    nmod_mat_init(x, s, max_c, prime);
    nmod_mat_init(v, max_r, s, prime);
    nmod_mat_init(schur, max_r, max_c, prime);
    for (slong i = 0; i < s; i++) {
        for (slong j = 0; j < i; j++)
            nmod_mat_entry(lower, i, j) = nmod_mat_entry(lu, i, j);
        for (slong j = i; j < s; j++)
            nmod_mat_entry(upper, i, j) = nmod_mat_entry(lu, i, pivots[j]);
    }
    dixon_eval_cache_t cache;
    dixon_eval_cache_init(&cache, matrix, ncols, params, ctx);
    slong old_r = 0, old_c = 0;
    slong nr = FLINT_MIN(2 * delta, max_r), nc = FLINT_MIN(2 * delta, max_c);
    for (;;) {
        if (nc > old_c) {
            nmod_mat_t b, y, new_x;
            nmod_mat_init(b, s, nc - old_c, prime);
            nmod_mat_init(y, s, nc - old_c, prime);
            nmod_mat_window_init(new_x, x, 0, old_c, s, nc);
            for (slong i = 0; i < s; i++)
                for (slong j = old_c; j < nc; j++)
                    nmod_mat_entry(b, i, j - old_c) =
                        dixon_eval_cached(&cache, r0[i], cx[j]);
            nmod_mat_solve_tril(y, lower, b, 1);
            nmod_mat_solve_triu(new_x, upper, y, 0);
            nmod_mat_window_clear(new_x);
            nmod_mat_clear(y);
            nmod_mat_clear(b);
        }
        for (slong i = old_r; i < nr; i++)
            for (slong j = 0; j < s; j++)
                nmod_mat_entry(v, i, j) = dixon_eval_cached(&cache, rx[i], c0[j]);
        /* Only the two new borders need evaluation and multiplication. */
        for (int border = 0; border < 2; border++) {
            slong rbegin = border ? old_r : 0, rend = border ? nr : old_r;
            slong cbegin = border ? 0 : old_c;
            if (rbegin == rend || cbegin == nc) continue;
            nmod_mat_t vv, xx, product;
            nmod_mat_window_init(vv, v, rbegin, 0, rend, s);
            nmod_mat_window_init(xx, x, 0, cbegin, s, nc);
            nmod_mat_init(product, rend - rbegin, nc - cbegin, prime);
            nmod_mat_mul(product, vv, xx);
            for (slong i = rbegin; i < rend; i++)
                for (slong j = cbegin; j < nc; j++)
                    nmod_mat_entry(schur, i, j) = nmod_sub(
                        dixon_eval_cached(&cache, rx[i], cx[j]),
                        nmod_mat_entry(product, i - rbegin, j - cbegin), schur->mod);
            nmod_mat_clear(product);
            nmod_mat_window_clear(xx);
            nmod_mat_window_clear(vv);
        }
        nmod_mat_t work, window;
        nmod_mat_window_init(window, schur, 0, 0, nr, nc);
        nmod_mat_init_set(work, window);
        nmod_mat_window_clear(window);
        slong *sp = flint_malloc((size_t) nr * sizeof(slong));
        slong sr = nmod_mat_lu(sp, work, 0);
        dixon_debug_log("  Schur reserve: %ld x %ld, rank=%ld, needed=%ld\n", nr, nc, sr, delta);
        if (sr >= delta) {
            slong *chosen_r = flint_malloc((size_t) delta * sizeof(slong));
            slong *chosen_c = flint_malloc((size_t) delta * sizeof(slong));
            next = 0;
            for (slong i = 0; i < delta; i++) {
                while (next < nc && nmod_mat_entry(work, i, next) == 0) next++;
                chosen_r[i] = sp[i];
                chosen_c[i] = next++;
            }
            memcpy(rows, r0, (size_t) s * sizeof(slong));
            memcpy(cols, c0, (size_t) s * sizeof(slong));
            for (slong i = 0; i < delta; i++) {
                rows[s + i] = rx[chosen_r[i]];
                cols[s + i] = cx[chosen_c[i]];
            }
            if (use_exchange)
                dixon_debug_log("  Degree-aware refinement backend: Schur-seeded exchange (rank=%ld, reserves=%ld/%ld)\n",
                                size, nrows - size, ncols - size);
            if (use_exchange)
                dixon_refine_schur_minor(&cache, nrows, ncols, rows, cols,
                                         lower, upper, x, v, schur, delta, chosen_r, chosen_c);
            flint_free(chosen_c);
            flint_free(chosen_r);
            success = 1;
        }
        flint_free(sp);
        nmod_mat_clear(work);
        if (success || (nr == max_r && nc == max_c)) break;
        old_r = nr;
        old_c = nc;
        nr = FLINT_MIN(2 * nr, max_r);
        nc = FLINT_MIN(2 * nc, max_c);
    }
    dixon_eval_cache_clear(&cache);
    nmod_mat_clear(schur);
    nmod_mat_clear(v);
    nmod_mat_clear(x);
    nmod_mat_clear(upper);
    nmod_mat_clear(lower);
    /* Release completion workspaces before allocating the streamed basis. */
    if (success && !use_exchange)
        dixon_refine_streamed_minor(matrix, nrows, ncols, size, rows, cols, params, ctx);
cleanup_indices:
    flint_free(used_c);
    flint_free(used_r);
    flint_free(cx);
    flint_free(rx);
    flint_free(pivots);
    flint_free(c0);
    flint_free(r0);
    return success;
}

static int dixon_build_predicted_mirror_indices(
    slong **rows_out, slong **cols_out, slong *size_out,
    const monom_t *row_monoms, slong nrows,
    const monom_t *col_monoms, slong ncols,
    hash_entry_t **row_index, slong row_hash_size,
    hash_entry_t **col_index, slong col_hash_size,
    slong nvars, const slong *H, slong h_len,
    slong sigma, slong dmin, slong predicted)
{
    unsigned char *standard = NULL, *selected = NULL;
    slong *rows = NULL, *cols = NULL, *degree = NULL, *tmp_exp = NULL;
    slong count = 0;

    if (!rows_out || !cols_out || !size_out || !H || predicted <= 0) return 0;
    standard = flint_calloc((size_t) ncols, 1);
    selected = flint_calloc((size_t) ncols, 1);
    degree = flint_malloc((size_t) ncols * sizeof(slong));
    tmp_exp = flint_malloc((size_t) nvars * sizeof(slong));
    if (!standard || !selected || !degree || !tmp_exp) goto fail;

    for (slong j = 0; j < ncols; j++) {
        degree[j] = 0;
        for (slong v = 0; v < nvars; v++) degree[j] += col_monoms[j].exp[v];
    }

    /* Monomial collection order is deterministic; within each degree its
       first H_e entries serve as the experimental almost-revlex staircase. */
    for (slong e = 0; e < h_len; e++) {
        slong need = H[e];
        for (slong j = 0; j < ncols && need > 0; j++) {
            if (degree[j] == e) { standard[j] = 1; need--; }
        }
        if (need != 0) goto fail;
    }

    for (slong j = 0; j < ncols; j++) {
        slong t = degree[j];
        if (standard[j] || (t >= dmin && 2 * t <= sigma && !standard[j]))
            selected[j] = 1;
    }
    for (slong j = 0; j < ncols; j++) {
        slong t = degree[j];
        if (!(t >= dmin && 2 * t <= sigma && !standard[j])) continue;
        slong large = sigma - t;
        if (large <= t) continue;
        memcpy(tmp_exp, col_monoms[j].exp, (size_t) nvars * sizeof(slong));
        tmp_exp[0] += large - t;
        slong shifted = lookup_monom_index(col_index, col_hash_size, tmp_exp, nvars);
        if (shifted < 0) goto fail;
        selected[shifted] = 1;
    }

    for (slong j = 0; j < ncols; j++) if (selected[j]) count++;
    if (count != predicted || count > nrows) goto fail;
    rows = flint_malloc((size_t) count * sizeof(slong));
    cols = flint_malloc((size_t) count * sizeof(slong));
    if (!rows || !cols) goto fail;
    count = 0;
    for (slong j = 0; j < ncols; j++) if (selected[j]) {
        for (slong v = 0; v < nvars; v++) tmp_exp[v] = col_monoms[j].exp[nvars - 1 - v];
        slong ri = lookup_monom_index(row_index, row_hash_size, tmp_exp, nvars);
        if (ri < 0) goto fail;
        rows[count] = ri; cols[count] = j; count++;
    }
    *rows_out = rows; *cols_out = cols; *size_out = count;
    flint_free(standard); flint_free(selected); flint_free(degree); flint_free(tmp_exp);
    return 1;
fail:
    flint_free(rows); flint_free(cols); flint_free(standard); flint_free(selected);
    flint_free(degree); flint_free(tmp_exp);
    return 0;
}

/* Enumerate the MQ dual support: each suffix has degree at most its
 * length. Rows are its reversal. This is a structural candidate universe,
 * independent of coefficient cancellations in any particular finite field. */
static int dixon_mq_support(slong *out, slong capacity, slong *count, slong *exp,
                            slong n, slong pos, slong remaining)
{
    if (pos == n) {
        if (*count >= capacity) return 0;
        memcpy(out + (*count)++ * n, exp, (size_t) n * sizeof(slong));
        return 1;
    }
    remaining = FLINT_MIN(remaining, n - pos);
    for (slong e = remaining; e >= 0; e--) {
        exp[pos] = e;
        if (!dixon_mq_support(out, capacity, count, exp, n, pos + 1, remaining - e)) return 0;
    }
    return 1;
}

/* Copy a computed coefficient rectangle into a complete local repair block. */
static void dixon_mq_copy_block(fq_mvpoly_t ***local, const fq_mvpoly_t *poly,
                              slong n, hash_entry_t **ri, slong rhs,
                              hash_entry_t **ci, slong chs,
                              const slong *rmap, const slong *cmap)
{
    for (slong t = 0; t < poly->nterms; t++) {
        const fq_monomial_t *term = &poly->terms[t];
        slong r = lookup_monom_index(ri, rhs, term->var_exp, n);
        slong c = lookup_monom_index(ci, chs, term->var_exp + n, n);
        FLINT_ASSERT(r >= 0 && c >= 0 && rmap[r] >= 0 && cmap[c] >= 0);
        fq_mvpoly_t *entry = get_matrix_entry_lazy(local, rmap[r], cmap[c], 1, poly->ctx);
        fq_mvpoly_add_term_fast(entry, NULL, term->par_exp, term->coeff);
    }
}

/* A low parameter-degree bound alone favours the extreme, very sparse
 * monomial layers. Select near the lost LU directions first, then include
 * complementary and neighbouring degree layers. Actual coefficient degrees
 * are still optimized by the existing refinement after rank completion. */
static void dixon_mq_choose_reserve(slong *indices, slong *map, slong count,
                                   slong size, slong extra, monom_t *monoms,
                                   slong n, const slong *lost,
                                   const slong *opposite_degrees,
                                   slong delta, slong sigma)
{
    for (slong k = 0; k < extra; k++) {
        slong anchor = k % delta;
        slong wanted = 0;
        const slong *reference = monoms[lost[anchor]].exp;
        for (slong v = 0; v < n; v++) wanted += reference[v];
        if (k >= (extra + 1) / 2) {
            slong round = (k - (extra + 1) / 2) / delta;
            slong offset = round == 0 ? 0 : (round & 1) ? (round + 1) / 2 : -round / 2;
            wanted = FLINT_MAX(0, FLINT_MIN(n, sigma - opposite_degrees[anchor] + offset));
        }
        slong best = -1, best_gap = WORD_MAX, best_distance = WORD_MAX, best_degree = -1;
        for (slong i = 0; i < count; i++) if (map[i] < 0) {
            slong degree = 0, distance = 0;
            for (slong v = 0; v < n; v++) {
                degree += monoms[i].exp[v];
                distance += FLINT_ABS(monoms[i].exp[v] - reference[v]);
            }
            slong gap = FLINT_ABS(degree - wanted);
            if (gap < best_gap || (gap == best_gap &&
                (distance < best_distance || (distance == best_distance && degree > best_degree)))) {
                best = i; best_gap = gap; best_distance = distance; best_degree = degree;
            }
        }
        FLINT_ASSERT(best >= 0);
        indices[size + k] = best; map[best] = size + k;
    }
}

/* Keep the exact candidate, project only two disjoint border strips, then
 * reuse Step 3's Schur completion and actual parameter-degree refinement.
 * Unknown entries are never presented as zeros to that machinery. On failure
 * result and the caller's chosen indices remain untouched. */
static int dixon_repair_mq_projection(fq_mvpoly_t *result, fq_mvpoly_t **matrix,
                                     slong n, monom_t *rm, slong nr,
                                     monom_t *cm, slong nc,
                                     hash_entry_t **ri, slong rhs,
                                     hash_entry_t **ci, slong chs,
                                     slong *rows, slong *cols, slong size,
                                     const nmod_mat_t lu, const slong *perm,
                                     slong s, ulong point, slong sigma,
                                     const fq_mvpoly_t *full)
{
    slong delta = size - s;
    if (s <= 0 || delta <= 0 || delta > 8) return 0;
    /* Border projections dominate the cost; use up to 32 useful directions
     * even for a one-dimensional deficit rather than eight extreme monomials. */
    slong budget = 32;
    slong lost[2][8], lost_degrees[2][8];
    for (slong i = 0; i < delta; i++) lost[0][i] = rows[perm[s + i]];
    slong pivot_row = 0, missing = 0;
    for (slong j = 0; j < size; j++) {
        if (pivot_row < s && nmod_mat_entry(lu, pivot_row, j) != 0) pivot_row++;
        else {
            if (missing >= delta) return 0;
            lost[1][missing++] = cols[j];
        }
    }
    if (missing != delta || pivot_row != s) return 0;
    for (slong i = 0; i < delta; i++) {
        lost_degrees[0][i] = lost_degrees[1][i] = 0;
        for (slong v = 0; v < n; v++) {
            lost_degrees[0][i] += rm[lost[0][i]].exp[v];
            lost_degrees[1][i] += cm[lost[1][i]].exp[v];
        }
    }
    slong counts[] = {nr, nc}, dims[2];
    monom_t *monoms[] = {rm, cm};
    slong *chosen[] = {rows, cols}, *indices[2], *targets[2], *maps[2];
    fq_index_degree_pair *orders[2];
    for (int axis = 0; axis < 2; axis++) {
        dims[axis] = FLINT_MIN(counts[axis], size + budget);
        indices[axis] = flint_malloc((size_t) dims[axis] * sizeof(slong));
        targets[axis] = flint_malloc((size_t) dims[axis] * n * sizeof(slong));
        maps[axis] = flint_malloc((size_t) counts[axis] * sizeof(slong));
        orders[axis] = flint_malloc((size_t) counts[axis] * sizeof(fq_index_degree_pair));
        for (slong i = 0; i < counts[axis]; i++) maps[axis][i] = -1;
        for (slong i = 0; i < size; i++) {
            indices[axis][i] = chosen[axis][i];
            maps[axis][chosen[axis][i]] = i;
        }
        dixon_mq_choose_reserve(indices[axis], maps[axis], counts[axis], size,
                                dims[axis] - size, monoms[axis], n, lost[axis],
                                lost_degrees[1 - axis], delta, sigma);
        if (g_dixon_verbose_level >= 2) {
            printf("  MQ repair %s: lost degrees", axis ? "columns" : "rows");
            for (slong i = 0; i < delta; i++) printf(" %ld", lost_degrees[axis][i]);
            printf("; reserve degrees");
            for (slong i = size; i < dims[axis]; i++) {
                slong degree = 0;
                for (slong v = 0; v < n; v++) degree += monoms[axis][indices[axis][i]].exp[v];
                printf(" %ld", degree);
            }
            printf("\n");
        }
        for (slong i = 0; i < dims[axis]; i++) {
            memcpy(targets[axis] + i * n, monoms[axis][indices[axis][i]].exp, (size_t) n * sizeof(slong));
            slong degree = 0;
            for (slong v = 0; v < n; v++) degree += targets[axis][i * n + v];
            orders[axis][i].index = i; orders[axis][i].degree = degree;
        }
        qsort(orders[axis], dims[axis], sizeof(fq_index_degree_pair), compare_fq_degrees);
    }
    fq_mvpoly_t ***local = flint_malloc((size_t) dims[0] * sizeof(*local));
    for (slong i = 0; i < dims[0]; i++)
        local[i] = flint_calloc((size_t) dims[1], sizeof(**local));
    int ok = 0;
    slong *lr = flint_malloc((size_t) size * sizeof(slong));
    slong *lc = flint_malloc((size_t) size * sizeof(slong));
    double phase_start = get_wall_time();
    dixon_mq_copy_block(local, result, n, ri, rhs, ci, chs, maps[0], maps[1]);
    dixon_debug_log("  MQ repair candidate copy: %.3fs\n", get_wall_time() - phase_start);
    dixon_info_log("  MQ Step 1 repair: deficit=%ld, adding %ld rows / %ld columns\n",
                   delta, dims[0] - size, dims[1] - size);
    for (int strip = 0; strip < 2; strip++) {
        slong rcount = strip == 0 ? dims[0] - size : size;
        slong ccount = strip == 0 ? dims[1] : dims[1] - size;
        if (!rcount || !ccount) continue;
        phase_start = get_wall_time();
        fq_mvpoly_t border;
        const slong *tr=targets[0]+(strip==0?size*n:0);
        const slong *tc=targets[1]+(strip==0?0:size*n);
        int border_ok=full ? fq_mq_project_full(&border,full,tr,rcount,tc,ccount)
                           : compute_fq_det_mq_projected_rect(&border,matrix,n+1,tr,rcount,tc,ccount);
        if (!border_ok) goto cleanup;
        dixon_debug_log("  MQ repair %s border DP: %.3fs\n",
                        strip == 0 ? "row" : "column", get_wall_time() - phase_start);
        phase_start = get_wall_time();
        dixon_mq_copy_block(local, &border, n, ri, rhs, ci, chs, maps[0], maps[1]);
        fq_mvpoly_clear(&border);
        dixon_debug_log("  MQ repair border insertion: %.3fs\n", get_wall_time() - phase_start);
    }
    for (slong i = 0; i < size; i++) lr[i] = lc[i] = i;
    fq_nmod_t params[1], value;
    fq_nmod_init(params[0], result->ctx); fq_nmod_init(value, result->ctx);
    fq_nmod_set_ui(params[0], point, result->ctx);
    phase_start = get_wall_time();
    ok = dixon_repair_predicted_minor(local, dims[0], dims[1], lr, lc, size,
                                     lu, perm, s, orders[0], orders[1], sigma,
                                     params, result->ctx);
    dixon_debug_log("  MQ repair Schur + degree selection: %.3fs\n", get_wall_time() - phase_start);
    if (ok) {
        phase_start = get_wall_time();
        nmod_mat_t check;
        nmod_mat_init(check, size, size, fq_nmod_ctx_prime(result->ctx));
        for (slong i = 0; i < size; i++) for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *entry = local[lr[i]][lc[j]];
            if (!entry) continue;
            evaluate_fq_mvpoly_at_params(value, entry, params);
            nmod_mat_entry(check, i, j) = nmod_poly_get_coeff_ui(value, 0);
        }
        ok = nmod_mat_rank(check) == size;
        nmod_mat_clear(check);
        dixon_debug_log("  MQ repair final rank verification: %.3fs\n", get_wall_time() - phase_start);
    }
    fq_nmod_clear(value, result->ctx); fq_nmod_clear(params[0], result->ctx);
    if (ok) {
        phase_start = get_wall_time();
        fq_mvpoly_t repaired;
        fq_mvpoly_init(&repaired, 2 * n, 1, result->ctx);
        slong *exp = flint_malloc((size_t) 2 * n * sizeof(slong));
        for (slong i = 0; i < size; i++) for (slong j = 0; j < size; j++) {
            fq_mvpoly_t *entry = local[lr[i]][lc[j]];
            if (!entry) continue;
            memcpy(exp, targets[0] + lr[i] * n, (size_t) n * sizeof(slong));
            memcpy(exp + n, targets[1] + lc[j] * n, (size_t) n * sizeof(slong));
            for (slong t = 0; t < entry->nterms; t++)
                fq_mvpoly_add_term_fast(&repaired, exp, entry->terms[t].par_exp, entry->terms[t].coeff);
        }
        flint_free(exp);
        fq_mvpoly_clear(result); *result = repaired;
        for (slong i = 0; i < size; i++) {
            rows[i] = indices[0][lr[i]]; cols[i] = indices[1][lc[i]];
        }
        dixon_debug_log("  MQ repair polynomial rebuild: %.3fs\n", get_wall_time() - phase_start);
    }
cleanup:
    for (slong i = 0; i < dims[0]; i++) {
        for (slong j = 0; j < dims[1]; j++) if (local[i][j]) {
            fq_mvpoly_clear(local[i][j]); flint_free(local[i][j]);
        }
        flint_free(local[i]);
    }
    flint_free(local); flint_free(lr); flint_free(lc);
    for (int axis = 0; axis < 2; axis++) {
        flint_free(indices[axis]); flint_free(targets[axis]);
        flint_free(maps[axis]); flint_free(orders[axis]);
    }
    return ok;
}

/* A successful return owns an initialized, projected Dixon polynomial.
 * Candidate verification follows the existing generic-rank heuristic: it
 * certifies non-singularity of this block, not an upper bound on full rank.
 * Deficient candidates first get a bounded, complete local Schur repair.
 * Only if that fails does the caller compute the FULL polynomial. */
static int dixon_mq_pencil_project(fq_mvpoly_t *result,fq_mvpoly_t **matrix,
    slong size,const slong *rows,slong nr,const slong *cols,slong nc)
{
    mq_pencil_stats stats={0};
    if(!compute_fq_det_mq_pencil_projected(result,matrix,size,rows,nr,cols,nc,&stats)) {
        dixon_info_log("  MQ Step 1 pencil: fallback (%s)\n",stats.reason?stats.reason:"unsupported projection");
        return compute_fq_det_mq_projected_rect(result,matrix,size,rows,nr,cols,nc);
    }
    dixon_info_log("  MQ Step 1 pencil: size=%ld, threads=%ld, total=%.6fs\n",stats.size,stats.threads,stats.total);
    dixon_info_log("    normalization=%.6fs, recurrence=%.6fs, assembly=%.6fs, peak matrix terms=%ld (closure projected)\n",
        stats.normalization,stats.recurrence,stats.assembly,stats.peak_terms);
    return 1;
}

static int dixon_try_mq_projection_from_full(fq_mvpoly_t *result, fq_mvpoly_t **matrix,
                                   const fq_mvpoly_t *polys, slong nvars,
                                   slong npars, det_method_t method, const fq_mvpoly_t *full)
{
    if (!g_dixon_mq_step1_filter ||
        method != DET_METHOD_RECURSIVE || npars != 1 || nvars < 2 ||
        nvars >= FLINT_BITS || fq_nmod_ctx_degree(polys[0].ctx) != 1) return 0;
    long *degrees = flint_malloc((size_t) (nvars + 1) * sizeof(long));
    int eligible = 1;
    for (slong i = 0; i <= nvars; i++) {
        degrees[i] = 0;
        for (slong t = 0; t < polys[i].nterms; t++) {
            slong d = 0;
            if (polys[i].terms[t].var_exp)
                for (slong v = 0; v < nvars; v++) d += polys[i].terms[t].var_exp[v];
            degrees[i] = FLINT_MAX(degrees[i], d);
            if (polys[i].terms[t].par_exp) d += polys[i].terms[t].par_exp[0];
            if (d > 2) eligible = 0;
        }
        if (degrees[i] != 2) eligible = 0;
    }
    slong *R = NULL, *H = NULL, rlen = 0, hlen = 0, sigma = 0, rank = 0;
    int model = eligible && dixon_rank_profile_from_degrees(&R, &rlen, &H, &hlen,
                                  &sigma, &rank, degrees, nvars + 1, nvars);
    flint_free(degrees);
    slong count = 0;
    for (slong i = 0; model && i < rlen; i++) {
        /* Bound support construction before allocating/enumerating. */
        if (R[i] > (1L << 18) - count) model = 0;
        else count += R[i];
    }
    flint_free(R);
    if (!model || rank <= 0 || rank >= count) { flint_free(H); return 0; }
    slong *exps = flint_malloc((size_t) count * nvars * sizeof(slong));
    slong *reverse = flint_malloc((size_t) count * nvars * sizeof(slong));
    slong *tmp = flint_calloc((size_t) nvars, sizeof(slong));
    slong actual = 0;
    if (!dixon_mq_support(exps, count, &actual, tmp, nvars, 0, nvars) || actual != count) {
        flint_free(tmp); flint_free(reverse); flint_free(exps); flint_free(H);
        return 0;
    }
    monom_t *rm = NULL, *cm = NULL;
    slong nr = 0, nc = 0, rcap = 0, ccap = 0, rhs = 16, chs = 16;
    hash_entry_t **ri = flint_calloc(16, sizeof(*ri));
    hash_entry_t **ci = flint_calloc(16, sizeof(*ci));
    /* Canonical MQ staircase: squarefree monomials first, then increasing
     * repeated-variable degree, breaking ties by descending lex order.
     * Unlike the legacy first-occurrence order this exists before Step 1.
     * Its candidate must pass verification; no almost-revlex theorem is
     * assumed for random inputs or small characteristic. */
    for (slong excess = 0; excess <= nvars; excess++) {
        for (slong j = 0; j < count; j++) {
            slong e = 0;
            for (slong v = 0; v < nvars; v++) e += FLINT_MAX(0, exps[j * nvars + v] - 1);
            if (e != excess) continue;
            for (slong v = 0; v < nvars; v++)
                reverse[j * nvars + v] = exps[j * nvars + nvars - 1 - v];
            dixon_intern_monom(&rm, &nr, &rcap, &ri, &rhs, reverse + j * nvars, nvars);
            dixon_intern_monom(&cm, &nc, &ccap, &ci, &chs, exps + j * nvars, nvars);
        }
    }
    slong *rows = NULL, *cols = NULL, size = 0;
    int ok = dixon_build_predicted_mirror_indices(&rows, &cols, &size,
                    rm, nr, cm, nc, ri, rhs, ci, chs, nvars, H, hlen, sigma, 2, rank);
    flint_free(H);
    slong *target_rows = NULL, *target_cols = NULL;
    int computed = 0;
    if (ok) {
        target_rows = flint_malloc((size_t) rank * nvars * sizeof(slong));
        target_cols = flint_malloc((size_t) rank * nvars * sizeof(slong));
        for (slong i = 0; i < rank; i++) {
            memcpy(target_rows + i * nvars, rm[rows[i]].exp, (size_t) nvars * sizeof(slong));
            memcpy(target_cols + i * nvars, cm[cols[i]].exp, (size_t) nvars * sizeof(slong));
        }
        computed = full ? fq_mq_project_full(result,full,target_rows,rank,target_cols,rank)
                        : g_dixon_mq_step1_pencil
                        ? dixon_mq_pencil_project(result,matrix,nvars+1,target_rows,rank,target_cols,rank)
                        : compute_fq_det_mq_projected(result,matrix,nvars+1,target_rows,target_cols,rank);
        ok = computed;
    }
    if (ok) {
        slong *rmap = flint_malloc((size_t) nr * sizeof(slong));
        slong *cmap = flint_malloc((size_t) nc * sizeof(slong));
        for (slong i = 0; i < nr; i++) rmap[i] = -1;
        for (slong j = 0; j < nc; j++) cmap[j] = -1;
        for (slong i = 0; i < rank; i++) { rmap[rows[i]] = i; cmap[cols[i]] = i; }
        ulong prime = fq_nmod_ctx_modulus(polys[0].ctx)->mod.n;
        nmod_mat_t values, best;
        nmod_mat_init(values, rank, rank, prime);
        int have_best = 0;
        slong best_s = 0;
        ulong best_point = 0;
        slong *perm = flint_malloc((size_t) rank * sizeof(slong));
        slong *best_perm = flint_malloc((size_t) rank * sizeof(slong));
        ok = 0;
        /* The same distinct points as before, but try 1 before 0 to avoid
         * a guaranteed rank drop when the candidate has parameter content.
         * MQ bounds the parameter degree by nvars+2. Cache powers once per
         * point; 0 and 1 require no modular multiplication per term. */
        ulong powers[FLINT_BITS + 2];
        for (ulong attempt = 0; attempt < FLINT_MIN(prime, UWORD(3)) && !ok; attempt++) {
            ulong point = attempt == 0 ? 1 : attempt == 1 ? 0 : 2;
            powers[0] = 1;
            for (slong d = 1; d <= nvars + 2; d++)
                powers[d] = nmod_mul(powers[d - 1], point, values->mod);
            nmod_mat_zero(values);
            for (slong t = 0; t < result->nterms; t++) {
                const fq_monomial_t *term = &result->terms[t];
                ulong power = term->par_exp ? term->par_exp[0] : 0;
                FLINT_ASSERT(power <= (ulong) nvars + 2);
                if (point == 0 && power != 0) continue;
                slong r = lookup_monom_index(ri, rhs, term->var_exp, nvars);
                slong c = lookup_monom_index(ci, chs, term->var_exp + nvars, nvars);
                FLINT_ASSERT(r >= 0 && c >= 0 && rmap[r] >= 0 && cmap[c] >= 0);
                ulong coeff = nmod_poly_get_coeff_ui(term->coeff, 0);
                if (point > 1 && power != 0)
                    coeff = nmod_mul(coeff, powers[power], values->mod);
                ulong *entry = nmod_mat_entry_ptr(values, rmap[r], cmap[c]);
                *entry = nmod_add(*entry, coeff, values->mod);
            }
            slong s = nmod_mat_lu(perm, values, 0);
            ok = s == rank;
            if (!ok && s > best_s && rank - s <= 8) {
                if (!have_best) { nmod_mat_init(best, rank, rank, prime); have_best = 1; }
                nmod_mat_set(best, values);
                memcpy(best_perm, perm, (size_t) rank * sizeof(slong));
                best_s = s; best_point = point;
            }
        }
        if (!ok && have_best) {
            ok = dixon_repair_mq_projection(result, matrix, nvars, rm, nr, cm, nc,
                     ri, rhs, ci, chs, rows, cols, rank, best, best_perm, best_s, best_point, sigma, full);
            if (ok) dixon_info_log("  MQ Step 1 Schur repair verified (degree-aware selection)\n");
        }
        if (have_best) nmod_mat_clear(best);
        flint_free(perm); flint_free(best_perm);
        nmod_mat_clear(values);
        flint_free(rmap); flint_free(cmap);
    }
    if (computed && !ok) fq_mvpoly_clear(result);
    if (computed)
        dixon_info_log("  MQ Step 1 projection: %ld x %ld candidate %s\n", rank, rank,
                       ok ? "verified" : full ? "failed; reusing full Dixon polynomial" : "failed; recomputing full Dixon polynomial");
    flint_free(target_rows); flint_free(target_cols);
    flint_free(rows); flint_free(cols);
    free_monom_index(ri, rhs); free_monom_index(ci, chs);
    flint_free(rm); flint_free(cm);
    flint_free(tmp); flint_free(reverse); flint_free(exps);
    return ok;
}

static int dixon_try_mq_projection(fq_mvpoly_t *result, fq_mvpoly_t **matrix,
                                   const fq_mvpoly_t *polys, slong nvars,
                                   slong npars, det_method_t method)
{
    return dixon_try_mq_projection_from_full(result,matrix,polys,nvars,npars,method,NULL);
}

/* A full experimental-backend polynomial is local to this call. Candidate extraction and
 * every repair border reuse it; failed verification never recomputes it. */
static int dixon_compute_step1(fq_mvpoly_t *result, fq_mvpoly_t **matrix,
                               const fq_mvpoly_t *polys, slong nvars,
                               slong npars, det_method_t method)
{
    if (g_dixon_mq_step1_pencil) {
        if(dixon_try_mq_projection(result,matrix,polys,nvars,npars,method))return 1;
        mq_pencil_stats stats={0};fq_mvpoly_t full;
        if (method==DET_METHOD_RECURSIVE && npars==1 &&
            compute_fq_det_mq_pencil(&full,matrix,nvars+1,&stats)) {
            dixon_info_log("  MQ Step 1 pencil: size=%ld, threads=%ld, total=%.6fs\n",stats.size,stats.threads,stats.total);
            dixon_info_log("    normalization=%.6fs, recurrence=%.6fs, assembly=%.6fs, peak matrix terms=%ld\n",
                stats.normalization,stats.recurrence,stats.assembly,stats.peak_terms);
            int verified=dixon_try_mq_projection_from_full(result,matrix,polys,nvars,npars,method,&full);
            if (verified)fq_mvpoly_clear(&full);else *result=full;
            return verified;
        }
        dixon_info_log("  MQ Step 1 pencil: fallback (%s)\n",stats.reason?stats.reason:"requires automatic/minor Step 1 and one parameter");
    } else if (g_dixon_mq_step1_simplex) {
        mq_simplex_stats stats={0};fq_mvpoly_t full;
        if (method==DET_METHOD_RECURSIVE && npars==1 &&
            compute_fq_det_mq_simplex(&full,matrix,nvars+1,&stats)) {
            dixon_info_log("  MQ Step 1 simplex: points=%ld, threads=%ld, total=%.6fs\n",stats.points,stats.threads,stats.total);
            dixon_info_log("    setup=%.6fs, entry evaluation=%.6fs, determinants=%.6fs, interpolation=%.6fs, packing=%.6fs (wall)\n",
                stats.setup,stats.entry_eval,stats.determinants,stats.interpolation,stats.packing);
            int verified=dixon_try_mq_projection_from_full(result,matrix,polys,nvars,npars,method,&full);
            if (verified) fq_mvpoly_clear(&full);
            else *result=full;
            return verified;
        }
        dixon_info_log("  MQ Step 1 simplex: fallback (%s)\n",stats.reason?stats.reason:"requires automatic/minor Step 1 and one parameter");
    }
    int verified=dixon_try_mq_projection(result,matrix,polys,nvars,npars,method);
    if (!verified)compute_fq_cancel_matrix_det(result,matrix,nvars,npars,method);
    return verified;
}

/* Metadata belongs to one selected matrix, never to a global monomial order.
 * In particular Step 3 may reorder or repair the candidate independently on
 * the two axes. The explicit permutations below retain its determinant sign. */
typedef struct {
    slong size, h, sigma;
    slong *rows, *cols, *rd, *cd;
    int odd;
} dixon_mq_step4_profile;

typedef struct {
    slong index, degree, excess, nvars;
    const slong *exp;
    int reverse;
} dixon_mq_step4_label;

static int dixon_mq_step4_label_cmp(const void *aa, const void *bb)
{
    const dixon_mq_step4_label *a = aa, *b = bb;
    if (a->degree != b->degree) return a->degree < b->degree ? -1 : 1;
    if (a->excess != b->excess) return a->excess < b->excess ? -1 : 1;
    for (slong v = 0; v < a->nvars; v++) {
        slong k = a->reverse ? a->nvars - 1 - v : v;
        if (a->exp[k] != b->exp[k]) return a->exp[k] > b->exp[k] ? -1 : 1;
    }
    return 0;
}

static int dixon_mq_step4_permutation_odd(const slong *p, slong n)
{
    unsigned char *seen = flint_calloc((size_t)n, 1);
    int odd = 0;
    for (slong i = 0; i < n; i++) if (!seen[i]) {
        slong len = 0, j = i;
        do { seen[j] = 1; len++; j = p[j]; } while (!seen[j]);
        odd ^= (len - 1) & 1;
    }
    flint_free(seen);
    return odd;
}

static void dixon_mq_step4_profile_clear(dixon_mq_step4_profile *p)
{
    flint_free(p->rows); flint_free(p->cols);
    flint_free(p->rd); flint_free(p->cd);
    memset(p, 0, sizeof(*p));
}

static int dixon_mq_step4_eligible(const fq_mvpoly_t *polys, slong m, slong npars)
{
    if (npars != 1 || m < 2 || fq_nmod_ctx_degree(polys[0].ctx) != 1) return 0;
    for (slong i = 0; i <= m; i++) {
        slong max_elim = 0;
        for (slong t = 0; t < polys[i].nterms; t++) {
            const fq_monomial_t *term = polys[i].terms + t;
            slong degree = 0;
            if (term->var_exp)
                for (slong v = 0; v < m; v++) degree += term->var_exp[v];
            max_elim = FLINT_MAX(max_elim, degree);
            if (term->par_exp) degree += term->par_exp[0];
            if (degree > 2) return 0;
        }
        if (max_elim != 2) return 0;
    }
    return 1;
}

static void dixon_mq_step4_prepare(dixon_mq_step4_profile *out,
                                  const monom_t *rm, const monom_t *cm,
                                  const slong *rows, const slong *cols,
                                  slong size, slong m, const long *degrees)
{
    slong *R = NULL, *H = NULL, rl = 0, hl = 0, sigma = 0, rho = 0, h = 0;
    if (!dixon_rank_profile_from_degrees(&R, &rl, &H, &hl, &sigma, &rho,
                                         degrees, m + 1, m)) return;
    flint_free(R);
    for (slong d = 0; d < hl; d++) h += H[d];
    if (size != rho || h <= 0 || h >= size) { flint_free(H); return; }
    out->rows = flint_malloc((size_t)size * sizeof(slong));
    out->cols = flint_malloc((size_t)size * sizeof(slong));
    out->rd = flint_malloc((size_t)size * sizeof(slong));
    out->cd = flint_malloc((size_t)size * sizeof(slong));
    dixon_mq_step4_label *labels = flint_malloc((size_t)size * sizeof(*labels));
    slong budget = (size - h) * sigma;
    int ok = 1;
    for (int axis = 0; axis < 2 && ok; axis++) {
        const monom_t *monoms = axis ? cm : rm;
        const slong *selected = axis ? cols : rows;
        slong *perm = axis ? out->cols : out->rows;
        slong *weight = axis ? out->cd : out->rd;
        for (slong i = 0; i < size; i++) {
            labels[i] = (dixon_mq_step4_label){i, 0, 0, m, monoms[selected[i]].exp, !axis};
            for (slong v = 0; v < m; v++) {
                labels[i].degree += labels[i].exp[v];
                labels[i].excess += FLINT_MAX(0, labels[i].exp[v] - 1);
            }
        }
        qsort(labels, (size_t)size, sizeof(*labels), dixon_mq_step4_label_cmp);
        slong *used = flint_calloc((size_t)hl, sizeof(slong));
        slong a = 0, b = h;
        for (slong i = 0; i < size; i++) {
            slong d = labels[i].degree;
            int core = d < hl && used[d] < H[d];
            if (core) used[d]++;
            slong k = core ? a++ : b++;
            if (k >= size) { ok = 0; break; }
            perm[k] = labels[i].index;
            weight[k] = d;
            if (!core) budget -= d;
        }
        for (slong d = 0; d < hl; d++) if (used[d] != H[d]) ok = 0;
        if (a != h || b != size) ok = 0;
        flint_free(used);
        if (ok) out->odd ^= dixon_mq_step4_permutation_odd(perm, size);
    }
    flint_free(labels); flint_free(H);
    if (!ok || budget != 0) { dixon_mq_step4_profile_clear(out); return; }
    out->size = size; out->h = h; out->sigma = sigma;
}

static int dixon_mq_step4_try(fq_nmod_poly_t det, const fq_nmod_poly_mat_t matrix,
                             const dixon_mq_step4_profile *p, const fq_nmod_ctx_t ctx)
{
    if (!p->size || p->size != matrix->r || matrix->r != matrix->c ||
        fq_nmod_ctx_degree(ctx) != 1) return 0;
    double start = get_wall_time();
    ulong prime = fq_nmod_ctx_prime(ctx), factor = 0;
    nmod_poly_mat_t B, core;
    nmod_poly_mat_init(B, p->size, p->size, prime);
    nmod_poly_mat_init(core, p->h, p->h, prime);
    for (slong i = 0; i < p->size; i++) for (slong j = 0; j < p->size; j++) {
        const fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(matrix, p->rows[i], p->cols[j]);
        nmod_poly_struct *target = nmod_poly_mat_entry(B, i, j);
        for (slong k = 0; k < entry->length; k++)
            nmod_poly_set_coeff_ui(target, k, nmod_poly_get_coeff_ui(entry->coeffs + k, 0));
    }
    int ok = nmod_poly_mat_mq_schur(core, &factor, B, p->rd, p->cd, p->h, p->sigma);
    nmod_poly_mat_clear(B);
    if (ok) {
        fq_nmod_poly_mat_t small;
        fq_nmod_poly_mat_init(small, p->h, p->h, ctx);
        fq_nmod_t coefficient;
        fq_nmod_init(coefficient, ctx);
        for (slong i = 0; i < p->h; i++) for (slong j = 0; j < p->h; j++) {
            const nmod_poly_struct *entry = nmod_poly_mat_entry(core, i, j);
            for (slong k = 0; k < entry->length; k++) {
                fq_nmod_set_ui(coefficient, nmod_poly_get_coeff_ui(entry, k), ctx);
                fq_nmod_poly_set_coeff(fq_nmod_poly_mat_entry(small, i, j), k, coefficient, ctx);
            }
        }
        dixon_info_log("  MQ Step 4 Schur: %ld -> %ld, compression %.3fs\n",
                       p->size, p->h, get_wall_time() - start);
        fq_nmod_poly_mat_det_iter(det, small, ctx);
        if (p->odd) factor = prime - factor;
        fq_nmod_set_ui(coefficient, factor, ctx);
        fq_nmod_poly_scalar_mul_fq_nmod(det, det, coefficient, ctx);
        fq_nmod_clear(coefficient, ctx);
        fq_nmod_poly_mat_clear(small, ctx);
    } else {
        dixon_info_log("  MQ Step 4 Schur: complement singular or degree check failed; using original determinant backend\n");
    }
    nmod_poly_mat_clear(core);
    return ok;
}

#include "dixon_projected_matrix.h"

static void extract_fq_coefficient_matrix_from_dixon_impl(fq_mvpoly_t ***coeff_matrix,
                                              fq_nmod_poly_mat_t *poly_matrix_out,
                                              slong *row_indices, slong *col_indices,
                                              slong *matrix_size,
                                              slong *extracted_x_power,
                                              const fq_mvpoly_t *dixon_poly,
                                              slong nvars, slong npars,
                                              char **var_names, char **par_names,
                                              const char *gen_name,
                                              const long *degrees, slong num_polys, int projected_verified,
                                              dixon_mq_step4_profile *mq_profile,
                                              nmod_poly_mat_t *prime_matrix_out, fq_mvpoly_t *consume) {
    dixon_info_log("\nStep 2: Construct Dixon matrix\n");
    if (extracted_x_power)
        *extracted_x_power = 0;
    clock_t step2_cpu_start = clock();
    double step2_wall_start = get_wall_time();
    double phase_start = step2_wall_start;
    
    slong *d0 = (slong*) flint_calloc(nvars, sizeof(slong));
    slong *d1 = (slong*) flint_calloc(nvars, sizeof(slong));
    dixon_debug_log("  Scanning Dixon polynomial degree bounds...\n");
    
    for (slong i = 0; i < dixon_poly->nterms; i++) {
        if (dixon_poly->terms[i].var_exp) {
            for (slong j = 0; j < nvars; j++) {
                if (dixon_poly->terms[i].var_exp[j] > d0[j]) {
                    d0[j] = dixon_poly->terms[i].var_exp[j];
                }
            }
            for (slong j = 0; j < nvars; j++) {
                if (dixon_poly->terms[i].var_exp[nvars + j] > d1[j]) {
                    d1[j] = dixon_poly->terms[i].var_exp[nvars + j];
                }
            }
        }
    }
    
    for (slong i = 0; i < nvars; i++) {
        d0[i]++;
        d1[i]++;
    }
    dixon_debug_log("  Degree bounds prepared in %.3f seconds\n",
                    get_wall_time() - phase_start);
    phase_start = get_wall_time();
    
    slong expected_rows = 1;
    slong expected_cols = 1;
    for (slong i = 0; i < nvars; i++) {
        expected_rows *= d0[i];
        expected_cols *= d1[i];
    }

    monom_t *x_monoms = NULL;
    monom_t *dual_monoms = NULL;
    slong nx_monoms = 0, ndual_monoms = 0;
    hash_entry_t **x_index = NULL;
    hash_entry_t **dual_index = NULL;
    slong x_hash_size = 0;
    slong dual_hash_size = 0;
    dixon_debug_log("  Collecting monomial supports...\n");
    
    slong *term_rows = flint_malloc((size_t) dixon_poly->nterms * sizeof(slong));
    slong *term_cols = flint_malloc((size_t) dixon_poly->nterms * sizeof(slong));
    collect_unique_monomials(&x_monoms, &nx_monoms,
                        &dual_monoms, &ndual_monoms,
                        &x_index, &x_hash_size, &dual_index, &dual_hash_size,
                        term_rows, term_cols, dixon_poly, d0, d1, nvars);
    dixon_debug_log("  Collected %ld row monomials and %ld dual monomials in %.3f seconds\n",
                    nx_monoms, ndual_monoms, get_wall_time() - phase_start);

    dixon_info_log("  Dixon matrix size: %ld x %ld\n", nx_monoms, ndual_monoms);
    
    if (nx_monoms == 0 || ndual_monoms == 0) {
        flint_free(term_rows);
        flint_free(term_cols);
        dixon_info_log("Warning: Empty coefficient matrix\n");
        *matrix_size = 0;
        dixon_maybe_print_step_time("Step 2", get_wall_time() - step2_wall_start);
        flint_free(d0);
        flint_free(d1);
        free_monom_index(x_index, x_hash_size);
        free_monom_index(dual_index, dual_hash_size);
        if (x_monoms) flint_free(x_monoms);
        if (dual_monoms) flint_free(dual_monoms);
        return;
    }

    int direct_projected = projected_verified && npars == 1 && (poly_matrix_out != NULL || prime_matrix_out != NULL)
                           && nx_monoms == ndual_monoms;
#ifdef DRSOLVE_MQ_STEP2_TEST
    if (dixon_step2_test_force_generic) direct_projected = 0;
#endif
    FLINT_ASSERT(!prime_matrix_out || direct_projected);
    if (direct_projected) {
        int own_indices = row_indices == NULL;
        if (own_indices) {
            row_indices = flint_malloc(nx_monoms*sizeof(slong));
            col_indices = flint_malloc(nx_monoms*sizeof(slong));
        }
        dixon_debug_log("  Filling verified projection directly into univariate matrix...\n");
        slong content = dixon_projected_poly_matrix(poly_matrix_out ? *poly_matrix_out : NULL,
                             prime_matrix_out, consume, row_indices,
                             col_indices, nx_monoms, dixon_poly, term_rows, term_cols);
        *coeff_matrix = NULL;
        *matrix_size = nx_monoms;
        if (extracted_x_power) *extracted_x_power = content;
        if (content > 0) {
            const char *v = (par_names && par_names[0]) ? par_names[0] : "x";
            dixon_info_log("  Pre-selection full-matrix %s-content: %s^%ld\n",v,v,content);
        }
        if (!prime_matrix_out) { flint_free(term_rows); flint_free(term_cols); }
        dixon_maybe_print_parallel_step_time("Step 2",
            (double)(clock()-step2_cpu_start)/CLOCKS_PER_SEC,get_wall_time()-step2_wall_start);
        clock_t direct_cpu = clock(); double direct_wall = get_wall_time();
        dixon_info_log("\nStep 3: Extract maximal-rank submatrix\n");
        dixon_debug_log("  Using MQ candidate verified in Step 1\n");
        dixon_info_log("  Submatrix size: %ld x %ld\n",nx_monoms,nx_monoms);
        if (mq_profile)
            dixon_mq_step4_prepare(mq_profile,x_monoms,dual_monoms,row_indices,
                                  col_indices,nx_monoms,nvars,degrees);
        if (own_indices) { flint_free(row_indices); flint_free(col_indices); }
        free_monom_index(x_index,x_hash_size); free_monom_index(dual_index,dual_hash_size);
        flint_free(x_monoms); flint_free(dual_monoms); flint_free(d0); flint_free(d1);
        dixon_maybe_print_parallel_step_time("Step 3",
            (double)(clock()-direct_cpu)/CLOCKS_PER_SEC,get_wall_time()-direct_wall);
        return;
    }

    fq_mvpoly_t ***full_matrix = (fq_mvpoly_t***) flint_malloc(nx_monoms * sizeof(fq_mvpoly_t**));
    slong *full_row_x_powers = NULL;
    slong *full_col_x_powers = NULL;
    for (slong i = 0; i < nx_monoms; i++) {
        full_matrix[i] = (fq_mvpoly_t**) flint_calloc(ndual_monoms, sizeof(fq_mvpoly_t*));
    }
    if (npars == 1) {
        full_row_x_powers = (slong *) flint_calloc((size_t) nx_monoms, sizeof(slong));
        full_col_x_powers = (slong *) flint_calloc((size_t) ndual_monoms, sizeof(slong));
    }

    phase_start = get_wall_time();
    dixon_debug_log("  Filling Dixon coefficient matrix...\n");
    fill_coefficient_matrix_optimized(full_matrix, nx_monoms, ndual_monoms,
                                     dixon_poly, npars, term_rows, term_cols);
    flint_free(term_rows);
    flint_free(term_cols);
    dixon_debug_log("  Coefficient matrix filled in %.3f seconds\n",
                    get_wall_time() - phase_start);
    slong preselection_x_power = extract_fq_full_matrix_x_content(full_matrix,
                                                                  nx_monoms,
                                                                  ndual_monoms,
                                                                  npars,
                                                                  full_row_x_powers,
                                                                  full_col_x_powers);
    if (preselection_x_power > 0) {
        const char *content_var = (par_names && par_names[0]) ? par_names[0] : "x";
        dixon_info_log("  Pre-selection full-matrix %s-content: %s^%ld\n",
                       content_var, content_var, preselection_x_power);
    }
    dixon_print_small_sparse_matrix("Dixon matrix", full_matrix, nx_monoms, ndual_monoms,
                                    x_monoms, dual_monoms, nvars,
                                    var_names, par_names, gen_name);
    dixon_maybe_print_parallel_step_time("Step 2",
                                         (double) (clock() - step2_cpu_start) / CLOCKS_PER_SEC,
                                         get_wall_time() - step2_wall_start);

    clock_t step3_cpu_start = clock();
    double step3_wall_start = get_wall_time();
    dixon_info_log("\nStep 3: Extract maximal-rank submatrix\n");
    slong *row_idx_array = NULL;
    slong *col_idx_array = NULL;
    slong num_rows, num_cols;
    
    if (projected_verified) {
        /* Step 1 emitted ONLY the verified target block. Every target row
         * and column occurs because its specialization is non-singular. */
        FLINT_ASSERT(nx_monoms == ndual_monoms);
        num_rows = num_cols = nx_monoms;
        row_idx_array = flint_malloc((size_t) num_rows * sizeof(slong));
        col_idx_array = flint_malloc((size_t) num_cols * sizeof(slong));
        for (slong i = 0; i < num_rows; i++) row_idx_array[i] = col_idx_array[i] = i;
        const char *reorder = getenv("DRSOLVE_PREDICT_REORDER");
        if (!reorder || strcmp(reorder, "0") != 0)
            reorder_fq_selected_minor_by_degree(full_matrix, row_idx_array,
                                               col_idx_array, num_rows, npars);
        dixon_debug_log("  Using MQ candidate verified in Step 1\n");
        goto coefficient_matrix_selected;
    }

    const char *predict_env = getenv("DRSOLVE_PREDICT_MAXRANK");
    int use_predicted_candidate =
        (npars == 1 && (predict_env == NULL || strcmp(predict_env, "0") != 0));
    if (npars == 0 || use_predicted_candidate) {
        dixon_debug_log("  Using rank-predicted scalar specialization (%ld parameter(s))...\n", npars);
        fq_nmod_mat_t eval_mat;
        int eval_mat_ready = 0;
        fq_nmod_t *eval_params = NULL;
        if (npars == 1) {
            eval_params = (fq_nmod_t *) flint_malloc(sizeof(fq_nmod_t));
            init_evaluation_parameters(eval_params, 1, dixon_poly->ctx, 0);
        }
        /* Keep the predicted npars=1 path candidate-only: the full scalar
           matrix is needed only by the non-predicted fallback. */
        if (!use_predicted_candidate) {
            fq_nmod_mat_init(eval_mat, nx_monoms, ndual_monoms, dixon_poly->ctx);
            eval_mat_ready = 1;
            clock_t scalar_eval_cpu_start = clock();
            double scalar_eval_wall_start = get_wall_time();
#ifdef _OPENMP
        #pragma omp parallel for schedule(static)
#endif
        for (slong i = 0; i < nx_monoms; i++) {
            for (slong j = 0; j < ndual_monoms; j++) {
                if (full_matrix[i][j] != NULL && full_matrix[i][j]->nterms > 0) {
                    if (npars == 1)
                        evaluate_fq_mvpoly_at_params(fq_nmod_mat_entry(eval_mat, i, j),
                                                     full_matrix[i][j], eval_params);
                    else
                        fq_nmod_set(fq_nmod_mat_entry(eval_mat, i, j),
                                   full_matrix[i][j]->terms[0].coeff, dixon_poly->ctx);
                } else {
                    fq_nmod_zero(fq_nmod_mat_entry(eval_mat, i, j), dixon_poly->ctx);
                }
            }
        }
            dixon_maybe_print_step_detail_time("Step 3 scalar matrix build",
                                               scalar_eval_cpu_start,
                                               scalar_eval_wall_start);
        }
        /* Heuristic rank candidate from the bihomogeneous support.  Rows are
         * ordered by decreasing x-degree and columns by increasing dual-degree
         * so that the candidate follows the Dixon anti-diagonal.  This is only
         * a fast path: sampled one-row/one-column extensions must not increase
         * rank, otherwise we fall back to the exact full-matrix rank. */
        clock_t prediction_cpu_start = clock();
        double prediction_wall_start = get_wall_time();
        slong sigma = 0;
        for (slong t = 0; t < dixon_poly->nterms; t++) {
            slong deg = 0;
            if (dixon_poly->terms[t].var_exp) {
                for (slong v = 0; v < 2 * nvars; v++)
                    deg += dixon_poly->terms[t].var_exp[v];
            }
            if (deg > sigma) sigma = deg;
        }
        fq_index_degree_pair *row_order =
            (fq_index_degree_pair *) flint_malloc((size_t) nx_monoms * sizeof(*row_order));
        fq_index_degree_pair *col_order =
            (fq_index_degree_pair *) flint_malloc((size_t) ndual_monoms * sizeof(*col_order));
        slong *row_counts = (slong *) flint_calloc((size_t) (sigma + 1), sizeof(slong));
        slong *col_counts = (slong *) flint_calloc((size_t) (sigma + 1), sizeof(slong));
        for (slong i = 0; i < nx_monoms; i++) {
            slong d = 0;
            for (slong v = 0; v < nvars; v++) d += x_monoms[i].exp[v];
            row_order[i].index = i; row_order[i].degree = d;
            if (d <= sigma) row_counts[d]++;
        }
        for (slong j = 0; j < ndual_monoms; j++) {
            slong d = 0;
            for (slong v = 0; v < nvars; v++) d += dual_monoms[j].exp[v];
            col_order[j].index = j; col_order[j].degree = d;
            if (d <= sigma) col_counts[d]++;
        }
        qsort(row_order, (size_t) nx_monoms, sizeof(*row_order), compare_fq_degrees);
        qsort(col_order, (size_t) ndual_monoms, sizeof(*col_order), compare_fq_degrees);
        slong predicted = 0;
        for (slong q = 0; q <= sigma; q++) {
            slong p = sigma - q;
            if (p <= sigma) predicted += FLINT_MIN(row_counts[p], col_counts[q]);
        }
        slong min_size = FLINT_MIN(nx_monoms, ndual_monoms);
        slong *model_H = NULL;
        slong model_h_len = 0;
        if (degrees != NULL && num_polys == nvars + 1) {
            slong *model_R = NULL;
            slong model_r_len = 0, model_sigma = 0, model_rank = 0;
            if (dixon_rank_profile_from_degrees(&model_R, &model_r_len,
                                                &model_H, &model_h_len,
                                                &model_sigma, &model_rank,
                                                degrees, num_polys, nvars)) {
                predicted = model_rank;
                sigma = model_sigma;
                dixon_debug_log("  Rank model: sigma=%ld, e0=%ld, predicted rank=%ld\n",
                                sigma, model_h_len, predicted);
            }
            flint_free(model_R);
        }
        predicted = FLINT_MIN(predicted, min_size);
        slong rank = 0;
        slong predicted_size = 0;
        long dmin = LONG_MAX;
        for (slong i = 0; i < num_polys; i++) if (degrees && degrees[i] < dmin) dmin = degrees[i];
        int candidate_ok =
            use_predicted_candidate && model_H != NULL && predicted > 0 && predicted < min_size &&
            dixon_build_predicted_mirror_indices(
                &row_idx_array, &col_idx_array, &predicted_size,
                x_monoms, nx_monoms, dual_monoms, ndual_monoms,
                x_index, x_hash_size, dual_index, dual_hash_size,
                nvars, model_H, model_h_len, sigma, (slong) dmin, predicted);
        if (use_predicted_candidate) {
            if (g_dixon_verbose_level >= 3)
                dixon_info_log("  Step 3 predicted rank: %ld\n", predicted);
            dixon_maybe_print_step_detail_time("Step 3 rank prediction",
                                               prediction_cpu_start,
                                               prediction_wall_start);
        }
        int candidate_constructed = candidate_ok;
        int schur_repaired = 0;
        clock_t verification_cpu_start = clock();
        double verification_wall_start = get_wall_time();
        if (candidate_ok) {
            if (fq_nmod_ctx_degree(dixon_poly->ctx) == 1) {
                if (g_dixon_verbose_level >= 3)
                    dixon_info_log("  Step 3 candidate verification backend: nmod_mat_lu\n");
                nmod_mat_t candidate;
                slong *perm = (slong *)
                    flint_malloc((size_t) predicted * sizeof(slong));
                nmod_mat_init(candidate, predicted, predicted,
                              fq_nmod_ctx_prime(dixon_poly->ctx));
                for (slong i = 0; i < predicted; i++) perm[i] = i;
#ifdef _OPENMP
                #pragma omp parallel
                {
                    fq_nmod_t value;
                    fq_nmod_init(value, dixon_poly->ctx);
                    #pragma omp for schedule(static)
                    for (slong i = 0; i < predicted; i++) {
                        for (slong j = 0; j < predicted; j++) {
                            fq_mvpoly_t *entry =
                                full_matrix[row_idx_array[i]][col_idx_array[j]];
                            if (entry != NULL && entry->nterms > 0) {
                                evaluate_fq_mvpoly_at_params(value, entry, eval_params);
                                nmod_mat_entry(candidate, i, j) =
                                    nmod_poly_get_coeff_ui(value, 0);
                            } else {
                                nmod_mat_entry(candidate, i, j) = 0;
                            }
                        }
                    }
                    fq_nmod_clear(value, dixon_poly->ctx);
                }
#else
                fq_nmod_t value;
                fq_nmod_init(value, dixon_poly->ctx);
                for (slong i = 0; i < predicted; i++) {
                    for (slong j = 0; j < predicted; j++) {
                        fq_mvpoly_t *entry =
                            full_matrix[row_idx_array[i]][col_idx_array[j]];
                        if (entry != NULL && entry->nterms > 0) {
                            evaluate_fq_mvpoly_at_params(value, entry, eval_params);
                            nmod_mat_entry(candidate, i, j) =
                                nmod_poly_get_coeff_ui(value, 0);
                        } else {
                            nmod_mat_entry(candidate, i, j) = 0;
                        }
                    }
                }
                fq_nmod_clear(value, dixon_poly->ctx);
#endif
                rank = nmod_mat_lu(perm, candidate, 0);
                if (rank < predicted && rank > 0 &&
                    dixon_repair_predicted_minor(full_matrix, nx_monoms, ndual_monoms,
                        row_idx_array, col_idx_array, predicted, candidate, perm, rank,
                        row_order, col_order, sigma, eval_params, dixon_poly->ctx)) {
                    num_rows = num_cols = predicted;
                    schur_repaired = 1;
                    dixon_debug_log("  Schur repair succeeded: base rank=%ld, added=%ld\n",
                                    rank, predicted - rank);
                }
                nmod_mat_clear(candidate);
                flint_free(perm);
            } else {
                if (g_dixon_verbose_level >= 3)
                    dixon_info_log("  Step 3 candidate verification backend: fq_nmod_mat_rank\n");
                fq_nmod_mat_t candidate;
                fq_nmod_mat_init(candidate, predicted, predicted, dixon_poly->ctx);
                for (slong i = 0; i < predicted; i++) {
                    for (slong j = 0; j < predicted; j++) {
                        fq_mvpoly_t *entry =
                            full_matrix[row_idx_array[i]][col_idx_array[j]];
                        if (entry != NULL && entry->nterms > 0)
                            evaluate_fq_mvpoly_at_params(
                                fq_nmod_mat_entry(candidate, i, j), entry, eval_params);
                        else
                            fq_nmod_zero(fq_nmod_mat_entry(candidate, i, j),
                                         dixon_poly->ctx);
                    }
                }
                rank = fq_nmod_mat_rank(candidate, dixon_poly->ctx);
                fq_nmod_mat_clear(candidate, dixon_poly->ctx);
            }
            if (!schur_repaired) {
                if (rank == predicted) num_rows = num_cols = predicted;
                else candidate_ok = 0;
            }
        }
        if (use_predicted_candidate) {
            const char *verification_label = candidate_ok
                ? "Step 3 predicted candidate verification (success)"
                : candidate_constructed
                    ? "Step 3 predicted candidate verification (rank mismatch)"
                    : "Step 3 predicted candidate verification (not constructed)";
            dixon_maybe_print_step_detail_time(verification_label,
                                               verification_cpu_start,
                                               verification_wall_start);
        }
        if (!candidate_ok) {
            flint_free(row_idx_array); flint_free(col_idx_array);
            row_idx_array = NULL; col_idx_array = NULL;
        }
        flint_free(model_H);
        flint_free(row_order);
        flint_free(col_order);
        flint_free(row_counts);
        flint_free(col_counts);
        if (use_predicted_candidate && !candidate_ok) {
            if (candidate_constructed) {
                dixon_debug_log("  Predicted candidate failed: size=%ld x %ld, "
                                "evaluated rank=%ld, deficiency=%ld; "
                                "using degree-aware selector\n",
                                predicted, predicted, rank, predicted - rank);
            } else {
                dixon_debug_log("  Predicted candidate failed: could not construct "
                                "%ld x %ld candidate; using degree-aware selector\n",
                                predicted, predicted);
            }
            if (eval_mat_ready) fq_nmod_mat_clear(eval_mat, dixon_poly->ctx);
            if (eval_params) {
                clear_evaluation_parameters(eval_params, 1, dixon_poly->ctx);
            }
            find_fq_optimal_maximal_rank_submatrix(full_matrix, nx_monoms, ndual_monoms,
                                                   &row_idx_array, &col_idx_array,
                                                   &num_rows, &num_cols, npars, -1);
            goto coefficient_matrix_selected;
        }
        if (candidate_ok) {
            dixon_debug_log("  Step 3 heuristic candidate accepted: %ld x %ld%s\n",
                            predicted, predicted, schur_repaired ? " (Schur repair)" : "");
            const char *reorder_env = getenv("DRSOLVE_PREDICT_REORDER");
            if (reorder_env == NULL || strcmp(reorder_env, "0") != 0)
                reorder_fq_selected_minor_by_degree(full_matrix,
                                                    row_idx_array,
                                                    col_idx_array,
                                                    predicted, npars);
            if (eval_mat_ready) fq_nmod_mat_clear(eval_mat, dixon_poly->ctx);
            if (eval_params) {
                clear_evaluation_parameters(eval_params, 1, dixon_poly->ctx);
            }
            goto coefficient_matrix_selected;
        }
        rank = fq_nmod_mat_rank(eval_mat, dixon_poly->ctx);
        slong actual_size = FLINT_MIN(rank, min_size);
        
        row_idx_array = (slong*) flint_malloc(actual_size * sizeof(slong));
        col_idx_array = (slong*) flint_malloc(actual_size * sizeof(slong));
        
        for (slong i = 0; i < actual_size; i++) {
            row_idx_array[i] = i;
            col_idx_array[i] = i;
        }
        num_rows = actual_size;
        num_cols = actual_size;
        
        if (eval_mat_ready) fq_nmod_mat_clear(eval_mat, dixon_poly->ctx);
        if (eval_params) {
            clear_evaluation_parameters(eval_params, 1, dixon_poly->ctx);
        }
    } else {
        slong small_size = 1;
        if (nx_monoms < small_size && ndual_monoms < small_size && 
            expected_rows < small_size && expected_cols < small_size) {
            dixon_debug_log("  Matrix is tiny; taking leading principal block directly...\n");
            slong min_size = FLINT_MIN(nx_monoms, ndual_monoms);
            row_idx_array = (slong*) flint_malloc(min_size * sizeof(slong));
            col_idx_array = (slong*) flint_malloc(min_size * sizeof(slong));
            
            for (slong i = 0; i < min_size; i++) {
                row_idx_array[i] = i;
                col_idx_array[i] = i;
            }
            num_rows = min_size;
            num_cols = min_size;
        } else {
            dixon_debug_log("  Selecting maximal-rank submatrix via specialization heuristics...\n");
            find_fq_optimal_maximal_rank_submatrix(full_matrix, nx_monoms, ndual_monoms,
                                                  &row_idx_array, &col_idx_array, 
                                                  &num_rows, &num_cols,
                                                  npars, -1);
        }
    }

coefficient_matrix_selected: ;

    slong submat_rank = FLINT_MIN(num_rows, num_cols);

    if (submat_rank == 0) {
        dixon_info_log("Warning: Matrix has rank 0\n");
        *matrix_size = 0;
        dixon_maybe_print_parallel_step_time("Step 3",
                                             (double) (clock() - step3_cpu_start) / CLOCKS_PER_SEC,
                                             get_wall_time() - step3_wall_start);

        if (full_matrix) {
            for (slong i = 0; i < nx_monoms; i++) {
                if (full_matrix[i]) {
                    for (slong j = 0; j < ndual_monoms; j++) {
                        if (full_matrix[i][j] != NULL) {
                            fq_mvpoly_clear(full_matrix[i][j]);
                            flint_free(full_matrix[i][j]);
                        }
                    }
                    flint_free(full_matrix[i]);
                }
            }
            flint_free(full_matrix);
        }

        if (row_idx_array) flint_free(row_idx_array);
        if (col_idx_array) flint_free(col_idx_array);
        free_monom_index(x_index, x_hash_size);
        free_monom_index(dual_index, dual_hash_size);
        if (x_monoms) flint_free(x_monoms);
        if (dual_monoms) flint_free(dual_monoms);
        if (d0) flint_free(d0);
        if (d1) flint_free(d1);
        if (full_row_x_powers) flint_free(full_row_x_powers);
        if (full_col_x_powers) flint_free(full_col_x_powers);
        return;
    }

    dixon_info_log("  Submatrix size: %ld x %ld\n", submat_rank, submat_rank);
    dixon_debug_log("  Copying %ld x %ld submatrix...\n",
                    submat_rank, submat_rank);

    clock_t copy_cpu_start = clock();
    double copy_wall_start = get_wall_time();

    /* The selected minor is copied into a second dense array below.  Release
       every unselected row before allocating that array; otherwise the full
       Dixon matrix and the (potentially multi-GB) coefficient matrix coexist
       at the peak.  Selected rows are released one-by-one after their
       selected polynomial objects have been moved below. */
    if (full_matrix && submat_rank > 0) {
        unsigned char *row_keep = (unsigned char *) flint_calloc((size_t) nx_monoms, 1);
        if (row_keep) {
            for (slong i = 0; i < submat_rank; i++) {
                if (row_idx_array[i] >= 0 && row_idx_array[i] < nx_monoms)
                    row_keep[row_idx_array[i]] = 1;
            }
            for (slong i = 0; i < nx_monoms; i++) {
                if (!full_matrix[i] || row_keep[i]) continue;
                for (slong j = 0; j < ndual_monoms; j++) {
                    if (full_matrix[i][j]) {
                        fq_mvpoly_clear(full_matrix[i][j]);
                        flint_free(full_matrix[i][j]);
                        full_matrix[i][j] = NULL;
                    }
                }
                flint_free(full_matrix[i]);
                full_matrix[i] = NULL;
            }
        }
        if (row_keep) flint_free(row_keep);
    }

    if (poly_matrix_out != NULL && npars == 1) {
        fq_nmod_poly_mat_init(*poly_matrix_out, submat_rank, submat_rank,
                              dixon_poly->ctx);
        for (slong i = 0; i < submat_rank; i++) {
            slong src_row = row_idx_array[i];
            for (slong j = 0; j < submat_rank; j++) {
                fq_mvpoly_t *source = full_matrix[src_row][col_idx_array[j]];
                if (source != NULL) {
                    fq_nmod_poly_struct *dst = fq_nmod_poly_mat_entry(
                        *poly_matrix_out, i, j);
                    for (slong t = 0; t < source->nterms; t++) {
                        slong degree = source->terms[t].par_exp
                            ? source->terms[t].par_exp[0] : 0;
                        fq_nmod_poly_set_coeff(dst, degree,
                                               source->terms[t].coeff,
                                               dixon_poly->ctx);
                    }
                }
            }
            /* Conversion is complete for this source row; release its
               multivariate entries before processing the next row. */
            if (full_matrix[src_row] != NULL) {
                for (slong j = 0; j < ndual_monoms; j++) {
                    if (full_matrix[src_row][j] != NULL) {
                        fq_mvpoly_clear(full_matrix[src_row][j]);
                        flint_free(full_matrix[src_row][j]);
                    }
                }
                flint_free(full_matrix[src_row]);
                full_matrix[src_row] = NULL;
            }
        }
        *coeff_matrix = NULL;
    } else {
    *coeff_matrix = (fq_mvpoly_t**) flint_malloc(submat_rank * sizeof(fq_mvpoly_t*));
    for (slong i = 0; i < submat_rank; i++) {
        slong src_row = row_idx_array[i];
        (*coeff_matrix)[i] = (fq_mvpoly_t*) flint_malloc(submat_rank * sizeof(fq_mvpoly_t));
        for (slong j = 0; j < submat_rank; j++) {
            slong src_col = col_idx_array[j];
            fq_mvpoly_t *source = full_matrix[src_row][src_col];
            if (source != NULL) {
                (*coeff_matrix)[i][j] = *source;
                flint_free(source);
                full_matrix[src_row][src_col] = NULL;
            } else {
                fq_mvpoly_init(&(*coeff_matrix)[i][j], 0, npars, dixon_poly->ctx);
            }
        }
        /* All selected entries in this source row have been moved. */
        for (slong j = 0; j < ndual_monoms; j++) {
            if (full_matrix[src_row][j]) {
                fq_mvpoly_clear(full_matrix[src_row][j]);
                flint_free(full_matrix[src_row][j]);
            }
        }
        flint_free(full_matrix[src_row]);
        full_matrix[src_row] = NULL;
    }
    }
    dixon_maybe_print_step_detail_time("Step 3 submatrix copy",
                                       copy_cpu_start,
                                       copy_wall_start);

    for (slong i = 0; i < submat_rank; i++) {
        row_indices[i] = row_idx_array[i];
        col_indices[i] = col_idx_array[i];
        if (extracted_x_power) {
            if (full_row_x_powers)
                *extracted_x_power += full_row_x_powers[row_idx_array[i]];
            if (full_col_x_powers)
                *extracted_x_power += full_col_x_powers[col_idx_array[i]];
        }
    }
    *matrix_size = submat_rank;
    if (mq_profile && poly_matrix_out)
        dixon_mq_step4_prepare(mq_profile, x_monoms, dual_monoms,
                               row_idx_array, col_idx_array, submat_rank, nvars, degrees);
    if (poly_matrix_out == NULL)
        dixon_print_small_dense_submatrix("Maximal Rank Submatrix", *coeff_matrix,
                                          submat_rank, submat_rank,
                                          row_idx_array, col_idx_array,
                                          x_monoms, dual_monoms, nvars,
                                          var_names, par_names, gen_name);

    if (full_matrix) {
        for (slong i = 0; i < nx_monoms; i++) {
            if (full_matrix[i]) {
                for (slong j = 0; j < ndual_monoms; j++) {
                    if (full_matrix[i][j]) {
                        fq_mvpoly_clear(full_matrix[i][j]);
                        flint_free(full_matrix[i][j]);
                    }
                }
                flint_free(full_matrix[i]);
            }
        }
        flint_free(full_matrix);
    }

    if (row_idx_array) flint_free(row_idx_array);
    if (col_idx_array) flint_free(col_idx_array);
    free_monom_index(x_index, x_hash_size);
    free_monom_index(dual_index, dual_hash_size);
    if (x_monoms) flint_free(x_monoms);
    if (dual_monoms) flint_free(dual_monoms);
    if (d0) flint_free(d0);
    if (d1) flint_free(d1);
    if (full_row_x_powers) flint_free(full_row_x_powers);
    if (full_col_x_powers) flint_free(full_col_x_powers);
    dixon_debug_log("  Completed in %.3f seconds\n",
                    get_wall_time() - step3_wall_start);
    dixon_maybe_print_parallel_step_time("Step 3",
                                         (double) (clock() - step3_cpu_start) / CLOCKS_PER_SEC,
                                         get_wall_time() - step3_wall_start);
}

void extract_fq_coefficient_matrix_from_dixon(fq_mvpoly_t ***coeff_matrix,
                                              fq_nmod_poly_mat_t *poly_matrix_out,
                                              slong *row_indices, slong *col_indices,
                                              slong *matrix_size,
                                              slong *extracted_x_power,
                                              const fq_mvpoly_t *dixon_poly,
                                              slong nvars, slong npars,
                                              char **var_names, char **par_names,
                                              const char *gen_name,
                                              const long *degrees, slong num_polys) {
    extract_fq_coefficient_matrix_from_dixon_impl(coeff_matrix, poly_matrix_out,
        row_indices, col_indices, matrix_size, extracted_x_power, dixon_poly,
        nvars, npars, var_names, par_names, gen_name, degrees, num_polys, 0, NULL, NULL, NULL);
}

// Compute determinant of cancellation matrix
void compute_fq_cancel_matrix_det(fq_mvpoly_t *result, fq_mvpoly_t **modified_M_mvpoly,
                                  slong nvars, slong npars, det_method_t method) {
    clock_t start = clock();
    method = dixon_resolve_step1_det_method(modified_M_mvpoly, nvars, npars, method, 0);
    switch (method) {
        case DET_METHOD_INTERPOLATION:
            fq_compute_det_by_interpolation_optimized(result, modified_M_mvpoly,
                                                      nvars + 1, nvars, npars,
                                                      modified_M_mvpoly[0][0].ctx, NULL);
            break;
        case DET_METHOD_KRONECKER:
            compute_fq_det_kronecker(result, modified_M_mvpoly, nvars + 1);
            break;
        case DET_METHOD_KRONECKER_NMOD:
            compute_fq_det_bareiss(result, modified_M_mvpoly, nvars + 1);
            break;
        case DET_METHOD_BALANCED_SPLIT:
            compute_fq_det_balanced_split_experimental(result, modified_M_mvpoly, nvars + 1);
            break;
        case DET_METHOD_HUANG:
            compute_fq_det_huang_interpolation(result, modified_M_mvpoly, nvars + 1);
            break;
        case DET_METHOD_RECURSIVE:
        default:
            compute_fq_det_recursive(result, modified_M_mvpoly, nvars + 1);
            break;
    }
    clock_t end = clock();
    (void) start;
    (void) end;
}

// ============ Get size of Dixon matrix ============

static slong dixon_binomial(slong n, slong k) {
    if (k > n || k < 0) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n - k) k = n - k;
    
    slong result = 1;
    for (slong i = 0; i < k; i++) {
        result = result * (n - i) / (i + 1);
    }
    return result;
}

static void dixon_add_generic_monomials(fq_mvpoly_t *poly,
                                        slong nvars,
                                        slong *exp,
                                        slong pos,
                                        slong remaining,
                                        flint_rand_t state,
                                        const fq_nmod_ctx_t ctx) {
    if (pos == nvars) {
        fq_nmod_t coeff;
        fq_nmod_init(coeff, ctx);
        do {
            fq_nmod_randtest(coeff, state, ctx);
        } while (fq_nmod_is_zero(coeff, ctx));

        slong *var_exp = (slong*) malloc(nvars * sizeof(slong));
        memcpy(var_exp, exp, nvars * sizeof(slong));
        fq_mvpoly_add_term(poly, var_exp, NULL, coeff);

        fq_nmod_clear(coeff, ctx);
        free(var_exp);
        return;
    }

    for (slong d = 0; d <= remaining; d++) {
        exp[pos] = d;
        dixon_add_generic_monomials(poly, nvars, exp, pos + 1, remaining - d, state, ctx);
    }
}

// Calculate actual Dixon matrix size by building the system and extracting matrix
slong dixon_matrix_size(slong nvars, slong degree, ulong prime, slong field_degree) {
    dixon_info_log("\nCalculating Dixon matrix size for n=%ld, d=%ld\n", nvars, degree);
    
    // Initialize field context
    fq_nmod_ctx_t ctx;
    fmpz_t p;
    fmpz_init_set_ui(p, prime);
    fq_nmod_ctx_init(ctx, p, field_degree, "t");
    fmpz_clear(p);
    
    // Generate generic polynomial system
    nvars--;
    slong npolys = nvars + 1;
    fq_mvpoly_t *polys = (fq_mvpoly_t*) malloc(npolys * sizeof(fq_mvpoly_t));
    
    // Generate dense generic polynomials of degree d
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 12345, 67890);
    
    for (slong i = 0; i < npolys; i++) {
        fq_mvpoly_init(&polys[i], nvars, 0, ctx);  // No parameters
                
        // Generate all monomials of degree <= d
        slong total_monomials = dixon_binomial(nvars + degree, degree);
        
        // Enumerate all monomials recursively
        slong *temp_exp = (slong*) calloc(nvars, sizeof(slong));
        for (slong d = 0; d <= degree; d++) {
            dixon_add_generic_monomials(&polys[i], nvars, temp_exp, 0, d, state, ctx);
        }
        free(temp_exp);
    }
    
    // Build cancellation matrix
    fq_mvpoly_t **M_mvpoly;
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, nvars, 0);
    
    // Perform row operations
    fq_mvpoly_t **modified_M_mvpoly;
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, 0);
    
    // Compute determinant to get Dixon polynomial
    fq_mvpoly_t d_poly;
    compute_fq_cancel_matrix_det(&d_poly, modified_M_mvpoly, nvars, 0, DET_METHOD_RECURSIVE);
    
    // Calculate degree bounds
    slong *d0 = (slong*) calloc(nvars, sizeof(slong));
    slong *d1 = (slong*) calloc(nvars, sizeof(slong));
    
    for (slong i = 0; i < d_poly.nterms; i++) {
        if (d_poly.terms[i].var_exp) {
            for (slong j = 0; j < nvars; j++) {
                if (d_poly.terms[i].var_exp[j] > d0[j]) {
                    d0[j] = d_poly.terms[i].var_exp[j];
                }
            }
            for (slong j = 0; j < nvars; j++) {
                if (d_poly.terms[i].var_exp[nvars + j] > d1[j]) {
                    d1[j] = d_poly.terms[i].var_exp[nvars + j];
                }
            }
        }
    }
    
    for (slong i = 0; i < nvars; i++) {
        d0[i]++;
        d1[i]++;
    }
    
    // Collect unique monomials
    monom_t *x_monoms = NULL;
    monom_t *dual_monoms = NULL;
    slong nx_monoms = 0, ndual_monoms = 0;
    
    hash_entry_t **x_index, **dual_index;
    slong x_hash_size, dual_hash_size;
    slong *term_rows = flint_malloc((size_t) d_poly.nterms * sizeof(slong));
    slong *term_cols = flint_malloc((size_t) d_poly.nterms * sizeof(slong));
    collect_unique_monomials(&x_monoms, &nx_monoms,
                            &dual_monoms, &ndual_monoms,
                            &x_index, &x_hash_size, &dual_index, &dual_hash_size,
                            term_rows, term_cols, &d_poly, d0, d1, nvars);
    flint_free(term_rows);
    flint_free(term_cols);
    free_monom_index(x_index, x_hash_size);
    free_monom_index(dual_index, dual_hash_size);
    
    // The actual matrix size (should be square, so take minimum)
    slong matrix_size = FLINT_MAX(nx_monoms, ndual_monoms);
    
    dixon_info_log("Dixon matrix actual size: %ld x %ld\n", nx_monoms, ndual_monoms);
    
    // Cleanup
    if (x_monoms) flint_free(x_monoms);
    if (dual_monoms) flint_free(dual_monoms);
    flint_free(d0);
    flint_free(d1);
    
    fq_mvpoly_clear(&d_poly);
    
    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);
    
    for (slong i = 0; i < npolys; i++) {
        fq_mvpoly_clear(&polys[i]);
    }
    free(polys);
    
    flint_rand_clear(state);
    fq_nmod_ctx_clear(ctx);
    
    return matrix_size;
}

// ============ Main Dixon resultant function ============
static long *dixon_polynomial_degrees(const fq_mvpoly_t *polys, slong npolys,
                                      slong nvars) {
    long *degrees = (long *) flint_calloc((size_t) npolys, sizeof(long));
    for (slong i = 0; i < npolys; i++) {
        for (slong t = 0; t < polys[i].nterms; t++) {
            long degree = 0;
            if (polys[i].terms[t].var_exp)
                for (slong v = 0; v < nvars; v++) degree += polys[i].terms[t].var_exp[v];
            if (degree > degrees[i]) degrees[i] = degree;
        }
    }
    return degrees;
}

void fq_dixon_resultant(fq_mvpoly_t *result, fq_mvpoly_t *polys,
                       slong nvars, slong npars) {
    cleanup_unified_workspace();
    dixon_info_log("\nStep 1: Build Dixon polynomial\n");
    clock_t step1_cpu_start = clock();
    double step1_wall_start = get_wall_time();
    fq_mvpoly_t **M_mvpoly;
    if (dixon_show_step_details()) dixon_debug_log("  Build Cancellation Matrix\n");
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, nvars, npars);
    dixon_print_small_named_dense_matrix("Cancellation Matrix", M_mvpoly,
                                         nvars + 1, nvars + 1,
                                         NULL, NULL, NULL, 1);
    
    // Display analysis of original matrix
    //analyze_fq_matrix_mvpoly(M_mvpoly, nvars + 1, nvars + 1, "Original Cancellation");
    
    fq_mvpoly_t **modified_M_mvpoly;
    if (dixon_show_step_details()) dixon_debug_log("  Perform Matrix Row Operations\n");
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, npars);
    
    fq_mvpoly_t d_poly;
    det_method_t step1_method = DET_METHOD_RECURSIVE;
    if (dixon_global_method_step1 != -1) {
        step1_method = dixon_global_method_step1;
    }
    step1_method = dixon_resolve_step1_det_method(modified_M_mvpoly, nvars, npars,
                                                  step1_method, 1);
    dixon_info_log("  Determinant method: %s\n", dixon_det_method_name(step1_method));
    if (dixon_show_step_details()) {
        dixon_debug_log("  Computing cancellation matrix determinant using %s...\n",
                        dixon_det_method_name(step1_method));
    }
    int projected_verified = dixon_compute_step1(&d_poly, modified_M_mvpoly,
                                polys, nvars, npars, step1_method);
    
    if (g_dixon_verbose_level >= 1 && d_poly.nterms <= 100) {
        dixon_info_log("  Dixon polynomial: %ld terms\n", d_poly.nterms);
        fq_mvpoly_print_expanded(&d_poly, "  DixonPoly", 1);
    } else {
        dixon_info_log("  Dixon polynomial: %ld terms (not shown)\n", d_poly.nterms);
    }
    if (g_dixon_debug_mode) {
        print_dixon_poly_actual_degrees(&d_poly, nvars, npars, NULL, NULL);
    }
    dixon_maybe_print_step_method_time("Step 1",
                                       step1_method,
                                       ((double)(clock() - step1_cpu_start) / CLOCKS_PER_SEC),
                                       get_wall_time() - step1_wall_start);

    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);
    
    fq_mvpoly_t **coeff_matrix = NULL;
    fq_nmod_poly_mat_t poly_matrix;
    int use_poly_matrix = (npars == 1 && (dixon_global_method_step4 == -1 ||
                           (g_dixon_mq_step4_schur && dixon_global_method_step4 == DET_METHOD_KRONECKER)));
    int use_prime_matrix = use_poly_matrix && projected_verified && fq_nmod_ctx_degree(polys[0].ctx)==1;
    nmod_poly_mat_t prime_matrix;
    slong *row_indices = use_prime_matrix ? NULL : (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong *col_indices = use_prime_matrix ? NULL : (slong*) flint_malloc(d_poly.nterms * sizeof(slong));
    slong matrix_size;
    slong extracted_x_power = 0;
    
    dixon_mq_step4_profile mq_profile = {0};
    int try_mq_step4 = g_dixon_mq_step4_schur && use_poly_matrix &&
                      dixon_mq_step4_eligible(polys, nvars, npars);
    long *rank_degrees = dixon_polynomial_degrees(polys, nvars + 1, nvars);
    extract_fq_coefficient_matrix_from_dixon_impl(&coeff_matrix,
                                            use_poly_matrix && !use_prime_matrix ? &poly_matrix : NULL,
                                            row_indices, col_indices,
                                            &matrix_size, &extracted_x_power, &d_poly, nvars, npars,
                                            NULL, NULL, NULL, rank_degrees, nvars + 1, projected_verified,
                                            try_mq_step4 ? &mq_profile : NULL,
                                            use_prime_matrix ? &prime_matrix : NULL,
                                            use_prime_matrix ? &d_poly : NULL);
    flint_free(rank_degrees);

    if (matrix_size > 0 && use_poly_matrix) {
        dixon_info_log("\nStep 4: Compute resultant\n");
        det_method_t coeff_method = dixon_global_method_step4 != -1
            ? dixon_global_method_step4 : DET_METHOD_KRONECKER;
        clock_t step4_cpu_start = clock();
        double step4_wall_start = get_wall_time();
        dixon_info_log("  Determinant method: %s\n",
                       dixon_det_method_name(coeff_method));
        fq_nmod_poly_t det_poly;
        fq_nmod_poly_init(det_poly, polys[0].ctx);
        if (g_dixon_mq_step4_schur && !mq_profile.size)
            dixon_info_log("  MQ Step 4 Schur: no eligible complement profile; using original determinant backend\n");
        if (use_prime_matrix)
            dixon_mq_native_det(det_poly, prime_matrix, &mq_profile, polys[0].ctx);
        else if (!dixon_mq_step4_try(det_poly, poly_matrix, &mq_profile, polys[0].ctx))
            fq_nmod_poly_mat_det_iter(det_poly, poly_matrix, polys[0].ctx);
        dixon_mq_step4_profile_clear(&mq_profile);
        fq_mvpoly_init(result, 0, 1, polys[0].ctx);
        for (slong i = 0; i <= fq_nmod_poly_degree(det_poly, polys[0].ctx); i++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, polys[0].ctx);
            fq_nmod_poly_get_coeff(coeff, det_poly, i, polys[0].ctx);
            if (!fq_nmod_is_zero(coeff, polys[0].ctx)) {
                slong par_exp[1] = {i + extracted_x_power};
                fq_mvpoly_add_term_fast(result, NULL, par_exp, coeff);
            }
            fq_nmod_clear(coeff, polys[0].ctx);
        }
        fq_nmod_poly_clear(det_poly, polys[0].ctx);
        if (use_prime_matrix) nmod_poly_mat_clear(prime_matrix);
        else fq_nmod_poly_mat_clear(poly_matrix, polys[0].ctx);
        dixon_maybe_print_step_method_time("Step 4", coeff_method,
                                           (double) (clock() - step4_cpu_start) / CLOCKS_PER_SEC,
                                           get_wall_time() - step4_wall_start);
        if (g_dixon_verbose_level >= 1 && result->nterms < 100) {
            fq_mvpoly_print(result, "  Final Resultant");
        } else {
            dixon_info_log("  Final resultant too large to display (%ld terms)\n",
                           result->nterms);
        }
        fq_mvpoly_make_monic(result);
        if (g_dixon_verbose_level >= 1)
            print_resultant_summary(result, NULL, 0);
    } else if (matrix_size > 0) {
        slong postselection_x_power =
            extract_fq_matrix_x_content(coeff_matrix, matrix_size, npars);
        extracted_x_power += postselection_x_power;
        if (postselection_x_power > 0) {
            dixon_info_log("  Pre-determinant row/column x-content: x^%ld\n",
                           postselection_x_power);
        }
        dixon_info_log("\nStep 4: Compute resultant\n");
        clock_t step4_cpu_start = clock();
        double step4_wall_start = get_wall_time();
        
        slong res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars+1, nvars, npars);
        dixon_debug_log("  Degree bound: %ld\n", res_deg_bound);
        
        ulong field_size = 1;
        for (slong i = 0; i < fq_nmod_ctx_degree(polys[0].ctx); i++) {
            field_size *= fq_nmod_ctx_prime(polys[0].ctx);
        }
        
        det_method_t coeff_method; // DET_METHOD_RECURSIVE DET_METHOD_KRONECKER DET_METHOD_INTERPOLATION DET_METHOD_HUANG
        #ifdef _OPENMP
        if (npars > 1) {
            coeff_method = DET_METHOD_INTERPOLATION;
        } else 
        #endif
        if (matrix_size < 9) {
            coeff_method = DET_METHOD_RECURSIVE;
        } else {
            coeff_method = DET_METHOD_KRONECKER;
        }
        //coeff_method = DET_METHOD_INTERPOLATION;
        if (dixon_global_method_step4 != -1) {
            coeff_method = dixon_global_method_step4;
        }
        dixon_info_log("  Determinant method: %s\n", dixon_det_method_name(coeff_method));
        
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                         npars, polys[0].ctx, coeff_method, res_deg_bound);
        fq_mvpoly_multiply_by_x_power_inplace(result, extracted_x_power);
        if (extracted_x_power > 0)
            dixon_info_log("  Restored extracted x-content: x^%ld\n", extracted_x_power);
        dixon_maybe_print_step_method_time("Step 4",
                                           coeff_method,
                                           ((double)(clock() - step4_cpu_start) / CLOCKS_PER_SEC),
                                           get_wall_time() - step4_wall_start);
        
        if (g_dixon_verbose_level >= 1 && result->nterms < 100) {
            fq_mvpoly_print(result, "  Final Resultant");
        } else {
            dixon_info_log("  Final resultant too large to display (%ld terms)\n", result->nterms);
        }

        fq_mvpoly_make_monic(result);
        if (g_dixon_verbose_level >= 1) print_resultant_summary(result, NULL, 0);
        // Cleanup coefficient matrix
        for (slong i = 0; i < matrix_size; i++) {
            for (slong j = 0; j < matrix_size; j++) {
                fq_mvpoly_clear(&coeff_matrix[i][j]);
            }
            flint_free(coeff_matrix[i]);
        }
        flint_free(coeff_matrix);
    } else {
        fq_mvpoly_init(result, 0, npars, polys[0].ctx);
        dixon_info_log("Warning: Empty coefficient matrix, resultant is 0\n");
    }
    
    // Cleanup
    flint_free(row_indices);
    flint_free(col_indices);
    
    fq_mvpoly_clear(&d_poly);

    dixon_info_log("\n=== Dixon Resultant Computation Complete ===\n");
}

void fq_dixon_resultant_with_names(fq_mvpoly_t *result, fq_mvpoly_t *polys, 
                                  slong nvars, slong npars,
                                  char **var_names, char **par_names, 
                                  const char *gen_name) {
    cleanup_unified_workspace();
    
    dixon_info_log("\nStep 1: Build Dixon polynomial\n");
    clock_t step1_cpu_start = clock();
    double step1_wall_start = get_wall_time();
    fq_mvpoly_t **M_mvpoly;
    dixon_debug_log("  Build Cancellation Matrix\n");
    build_fq_cancellation_matrix_mvpoly(&M_mvpoly, polys, nvars, npars);
    dixon_print_small_named_dense_matrix("Cancellation Matrix", M_mvpoly,
                                         nvars + 1, nvars + 1,
                                         var_names, par_names, gen_name, 1);
    
    fq_mvpoly_t **modified_M_mvpoly;
    dixon_debug_log("  Perform Matrix Row Operations\n");
    perform_fq_matrix_row_operations_mvpoly(&modified_M_mvpoly, &M_mvpoly, nvars, npars);

    fq_mvpoly_t d_poly;
    det_method_t step1_method = DET_METHOD_RECURSIVE;
    if (dixon_global_method_step1 != -1) {
        step1_method = dixon_global_method_step1;
    }
    step1_method = dixon_resolve_step1_det_method(modified_M_mvpoly, nvars, npars,
                                                  step1_method, 1);
    dixon_info_log("  Determinant method: %s\n", dixon_det_method_name(step1_method));
    dixon_debug_log("  Computing cancellation matrix determinant using %s...\n",
                    dixon_det_method_name(step1_method));
    int projected_verified = dixon_compute_step1(&d_poly, modified_M_mvpoly,
                                polys, nvars, npars, step1_method);
    
    if (g_dixon_verbose_level >= 1 && d_poly.nterms <= 100) {
        dixon_info_log("  Dixon polynomial: %ld terms\n", d_poly.nterms);
        fq_mvpoly_print_with_names(&d_poly, "  DixonPoly", var_names, par_names, gen_name, 1);
    } else {
        dixon_info_log("  Dixon polynomial: %ld terms (not shown)\n", d_poly.nterms);
    }
    if (g_dixon_debug_mode) {
        print_dixon_poly_actual_degrees(&d_poly, nvars, npars, var_names, par_names);
    }
    dixon_maybe_print_step_method_time("Step 1",
                                       step1_method,
                                       ((double)(clock() - step1_cpu_start) / CLOCKS_PER_SEC),
                                       get_wall_time() - step1_wall_start);
    
    fq_mvpoly_t **coeff_matrix = NULL;
    fq_nmod_poly_mat_t poly_matrix;
    int use_poly_matrix = (npars == 1 && (dixon_global_method_step4 == -1 ||
                           (g_dixon_mq_step4_schur && dixon_global_method_step4 == DET_METHOD_KRONECKER)));
    int use_prime_matrix = use_poly_matrix && projected_verified && fq_nmod_ctx_degree(polys[0].ctx)==1;
    nmod_poly_mat_t prime_matrix;
    slong max_indices = d_poly.nterms > 0 ? d_poly.nterms : 1;
    slong *row_indices = use_prime_matrix ? NULL : (slong*) flint_malloc(max_indices * sizeof(slong));
    slong *col_indices = use_prime_matrix ? NULL : (slong*) flint_malloc(max_indices * sizeof(slong));
    slong matrix_size;
    slong extracted_x_power = 0;
    
    dixon_mq_step4_profile mq_profile = {0};
    int try_mq_step4 = g_dixon_mq_step4_schur && use_poly_matrix &&
                      dixon_mq_step4_eligible(polys, nvars, npars);
    long *rank_degrees = dixon_polynomial_degrees(polys, nvars + 1, nvars);
    extract_fq_coefficient_matrix_from_dixon_impl(&coeff_matrix,
                                            use_poly_matrix && !use_prime_matrix ? &poly_matrix : NULL,
                                            row_indices, col_indices,
                                            &matrix_size, &extracted_x_power, &d_poly, nvars, npars,
                                            var_names, par_names, gen_name,
                                            rank_degrees, nvars + 1, projected_verified,
                                            try_mq_step4 ? &mq_profile : NULL,
                                            use_prime_matrix ? &prime_matrix : NULL,
                                            use_prime_matrix ? &d_poly : NULL);
    flint_free(rank_degrees);

    if (matrix_size > 0 && use_poly_matrix) {
        dixon_info_log("\nStep 4: Compute resultant\n");
        det_method_t coeff_method = dixon_global_method_step4 != -1
            ? dixon_global_method_step4 : DET_METHOD_KRONECKER;
        clock_t step4_cpu_start = clock();
        double step4_wall_start = get_wall_time();
        dixon_info_log("  Determinant method: %s\n",
                       dixon_det_method_name(coeff_method));
        fq_nmod_poly_t det_poly;
        fq_nmod_poly_init(det_poly, polys[0].ctx);
        if (g_dixon_mq_step4_schur && !mq_profile.size)
            dixon_info_log("  MQ Step 4 Schur: no eligible complement profile; using original determinant backend\n");
        if (use_prime_matrix)
            dixon_mq_native_det(det_poly, prime_matrix, &mq_profile, polys[0].ctx);
        else if (!dixon_mq_step4_try(det_poly, poly_matrix, &mq_profile, polys[0].ctx))
            fq_nmod_poly_mat_det_iter(det_poly, poly_matrix, polys[0].ctx);
        dixon_mq_step4_profile_clear(&mq_profile);
        fq_mvpoly_init(result, 0, 1, polys[0].ctx);
        for (slong i = 0; i <= fq_nmod_poly_degree(det_poly, polys[0].ctx); i++) {
            fq_nmod_t coeff;
            fq_nmod_init(coeff, polys[0].ctx);
            fq_nmod_poly_get_coeff(coeff, det_poly, i, polys[0].ctx);
            if (!fq_nmod_is_zero(coeff, polys[0].ctx)) {
                slong par_exp[1] = {i + extracted_x_power};
                fq_mvpoly_add_term_fast(result, NULL, par_exp, coeff);
            }
            fq_nmod_clear(coeff, polys[0].ctx);
        }
        fq_nmod_poly_clear(det_poly, polys[0].ctx);
        if (use_prime_matrix) nmod_poly_mat_clear(prime_matrix);
        else fq_nmod_poly_mat_clear(poly_matrix, polys[0].ctx);
        dixon_maybe_print_step_method_time("Step 4", coeff_method,
                                           (double) (clock() - step4_cpu_start) / CLOCKS_PER_SEC,
                                           get_wall_time() - step4_wall_start);
        if (g_dixon_verbose_level >= 1 && result->nterms < 100) {
            fq_mvpoly_print_with_names(result, "  Final Resultant", NULL,
                                       par_names, gen_name, 0);
        } else {
            dixon_info_log("  Final resultant too large to display (%ld terms)\n",
                           result->nterms);
        }
        fq_mvpoly_make_monic(result);
        if (g_dixon_verbose_level >= 1)
            print_resultant_summary(result, par_names, npars);
    } else if (matrix_size > 0) {
        slong postselection_x_power =
            extract_fq_matrix_x_content(coeff_matrix, matrix_size, npars);
        extracted_x_power += postselection_x_power;
        if (postselection_x_power > 0) {
            const char *content_var = (par_names && par_names[0]) ? par_names[0] : "x";
            dixon_info_log("  Pre-determinant row/column %s-content: %s^%ld\n",
                           content_var, content_var, postselection_x_power);
        }
        dixon_info_log("\nStep 4: Compute resultant\n");
        clock_t step4_cpu_start = clock();
        double step4_wall_start = get_wall_time();
        
        slong res_deg_bound = compute_fq_dixon_resultant_degree_bound(polys, nvars+1, nvars, npars);
        dixon_debug_log("  Degree bound: %ld\n", res_deg_bound);
        
        det_method_t coeff_method;
        #ifdef _OPENMP
        if (npars > 1) {
            coeff_method = DET_METHOD_INTERPOLATION;
            if (matrix_size < 10) {
                coeff_method = DET_METHOD_RECURSIVE;
            }
        } else 
        #endif
        if (matrix_size < 10) {
            coeff_method = DET_METHOD_RECURSIVE;
        } else {
            coeff_method = DET_METHOD_KRONECKER;
        }
        if (dixon_global_method_step4 != -1) {
            coeff_method = dixon_global_method_step4;
        }
        dixon_info_log("  Determinant method: %s\n", dixon_det_method_name(coeff_method));
        
        compute_fq_coefficient_matrix_det(result, coeff_matrix, matrix_size,
                                         npars, polys[0].ctx, coeff_method, res_deg_bound);
        fq_mvpoly_multiply_by_x_power_inplace(result, extracted_x_power);
        if (extracted_x_power > 0) {
            const char *content_var = (par_names && par_names[0]) ? par_names[0] : "x";
            dixon_info_log("  Restored extracted %s-content: %s^%ld\n",
                           content_var, content_var, extracted_x_power);
        }
        dixon_maybe_print_step_method_time("Step 4",
                                           coeff_method,
                                           ((double)(clock() - step4_cpu_start) / CLOCKS_PER_SEC),
                                           get_wall_time() - step4_wall_start);
        
        if (g_dixon_verbose_level >= 1 && result->nterms < 100) {
            fq_mvpoly_print_with_names(result, "  Final Resultant", NULL, par_names, gen_name, 0);
        } else {
            dixon_info_log("  Final resultant too large to display (%ld terms)\n", result->nterms);
        }
        fq_mvpoly_make_monic(result);
        if (g_dixon_verbose_level >= 1) print_resultant_summary(result, par_names, npars);
        
        for (slong i = 0; i < matrix_size; i++) {
            for (slong j = 0; j < matrix_size; j++) {
                fq_mvpoly_clear(&coeff_matrix[i][j]);
            }
            flint_free(coeff_matrix[i]);
        }
        flint_free(coeff_matrix);
    } else {
        fq_mvpoly_init(result, 0, npars, polys[0].ctx);
        dixon_info_log("Warning: Empty coefficient matrix, resultant is 0\n");
    }
    
    flint_free(row_indices);
    flint_free(col_indices);
    
    fq_mvpoly_clear(&d_poly);
    
    for (slong i = 0; i <= nvars; i++) {
        for (slong j = 0; j <= nvars; j++) {
            fq_mvpoly_clear(&M_mvpoly[i][j]);
            fq_mvpoly_clear(&modified_M_mvpoly[i][j]);
        }
        flint_free(M_mvpoly[i]);
        flint_free(modified_M_mvpoly[i]);
    }
    flint_free(M_mvpoly);
    flint_free(modified_M_mvpoly);

    dixon_info_log("\n=== Dixon Resultant Computation Complete ===\n");
}



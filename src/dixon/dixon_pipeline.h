/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DIXON_PIPELINE_H
#define DIXON_PIPELINE_H
#include "dixon_flint.h"
#include "dixon_scalar_screen.h"

/* Scoped per-thread screen destination, installed only by the screening API. */
dixon_scalar_screen_report_t *dixon_scalar_screen_set_report(dixon_scalar_screen_report_t *report);
int dixon_scalar_screen_active(void);
int dixon_scalar_screen_value(void);

/* Borrowed matrix/labels; returned indices are owned by the caller. */
void dixon_select_submatrix(unified_mpoly_struct ***matrix, slong nrows, slong ncols,
    const monom_t *rows, const monom_t *cols, slong nvars, slong npars,
    const long *degrees, slong num_polys, slong ksy_constant_col,
    const fq_nmod_ctx_t ctx, slong **selected_rows, slong **selected_cols, slong *size);
long *dixon_polynomial_degrees(const unified_mpoly_struct *polys, slong npolys, slong nvars);

/* Optional recursive native backend metadata; raw indices survive support trimming. */
typedef struct {
    const slong *rows, *cols;
    slong y_degree;
    int enabled;
} dixon_recursive_det_options;
void dixon_compute_dense_resultant(unified_mpoly_struct *result,
    unified_mpoly_struct **matrix, slong size, unified_mpoly_struct *polys,
    slong nvars, slong npars, slong extracted_power, char **par_names,
    const dixon_recursive_det_options *recursive);
#endif

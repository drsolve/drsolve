/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Private bridge between recursive construction and checked Step 4. */
#ifndef DIXON_NATIVE_H
#define DIXON_NATIVE_H
#include "dixon_flint.h"

int dixon_bivariate_native_eligible(const unified_mpoly_struct *polys, slong nvars, slong npars);
void dixon_bivariate_native_selected_det(unified_mpoly_struct *result,
                                         unified_mpoly_struct **matrix, slong size,
                                         const slong *rows, const slong *cols, slong y_degree,
                                         const unified_mpoly_struct *polys);
#endif

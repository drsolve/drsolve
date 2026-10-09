/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DIXON_SCALAR_SCREEN_H
#define DIXON_SCALAR_SCREEN_H

#include <flint/flint.h>

/* A one-sided affine consistency test, not a scalar resultant or root count. */
typedef enum {
    DIXON_SCREEN_NOT_RUN = 0,
    DIXON_SCREEN_NO_COMMON_ZERO,
    DIXON_SCREEN_INCONCLUSIVE
} dixon_scalar_screen_status_t;

typedef struct {
    dixon_scalar_screen_status_t status;
    slong nrows, ncols;
    slong rank, nonconstant_rank;
    slong constant_col; /* Identified from monomial exponents; -1 if absent. */
} dixon_scalar_screen_report_t;

#endif

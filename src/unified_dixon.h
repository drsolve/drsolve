/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef UNIFIED_DIXON_H
#define UNIFIED_DIXON_H
#include "unified_mpoly_interface.h"

/* Map source variables to destination variables; rings share a field. */
void unified_mpoly_reindex(unified_mpoly_t out, const unified_mpoly_t in, const slong *map);
/* Prime fields: n+1 input polynomials, first n variables eliminated;
 * the result ring has only parameters. Returns 0 on unsupported input/failure. */
int unified_dixon_resultant(unified_mpoly_t out, unified_mpoly_t *polys, slong n);
/* Matrix helpers: specialization-based degree-aware selection, and exact
 * weighted-degree Schur reduction in place. The latter returns core size h;
 * det(original) = factor * det(matrix[0:h,0:h]), including permutation signs. */
slong unified_dixon_select_minor(slong *rows, slong *cols, unified_mpoly_t **matrix, slong nr,
                                 slong nc, const slong *rd, const slong *cd, slong sigma);
slong unified_dixon_schur(unified_mpoly_t **matrix, slong n, unified_mpoly_t factor, slong *rd,
                          slong *cd, slong sigma);
#endif

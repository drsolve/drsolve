/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DIXON_RANK_PREDICTOR_H
#define DIXON_RANK_PREDICTOR_H

#include <stdio.h>
#include <flint/fmpz.h>

/* Generic dense Dixon rank for n positive equation degrees (n-1 eliminands).
 * Pass NULL for details or a stream for the boundary deficits and totals.
 * Return 1 on success, 0 for unsupported inputs. Caller initializes rank. */
int dixon_generic_rank_from_degrees(fmpz_t rank, const long *degrees, slong n, FILE *details);

#endif

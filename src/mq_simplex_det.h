/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef MQ_SIMPLEX_DET_H
#define MQ_SIMPLEX_DET_H
#include "fq_mvpoly.h"

typedef struct {
    slong points, threads;
    double setup, entry_eval, determinants, interpolation, packing, total;
    const char *reason;
} mq_simplex_stats;

/* Full determinant of a prime-field MQ divided-difference matrix. Exact nodes
 * 0..n+1 require p>n+1. Rejects unsupported shapes/workspaces without touching
 * result. On success initializes result; caller owns it. Thread count follows
 * OpenMP's configured maximum. No global polynomial cache is used. */
int compute_fq_det_mq_simplex(fq_mvpoly_t *result, fq_mvpoly_t **matrix,
                            slong n, mq_simplex_stats *stats);
#endif

/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef MQ_PENCIL_DET_H
#define MQ_PENCIL_DET_H
#include "fq_mvpoly.h"

typedef struct {
    slong size, threads, peak_terms;
    double normalization, recurrence, assembly, total;
    const char *reason;
} mq_pencil_stats;

/* Exact Faddeev-LeVerrier degree recurrence for MQ divided differences.
 * Prime characteristic must exceed n-1, and the parameter coefficient matrix
 * of the linear rows must have full row rank. On rejection result is untouched;
 * on success result is initialized and owned by the caller. */
int compute_fq_det_mq_pencil(fq_mvpoly_t *result, fq_mvpoly_t **matrix,
                           slong n, mq_pencil_stats *stats);
/* Project throughout the recurrence to the downward closure, then return
 * exact target coefficients. Same ownership/rejection contract as above. */
int compute_fq_det_mq_pencil_projected(fq_mvpoly_t *result,fq_mvpoly_t **matrix,
    slong n,const slong *rows,slong nr,const slong *cols,slong nc,mq_pencil_stats *stats);
#endif

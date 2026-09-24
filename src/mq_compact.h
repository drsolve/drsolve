/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_MQ_COMPACT_H
#define DRSOLVE_MQ_COMPACT_H
#include "fq_mvpoly.h"

/* Prime-field, one-parameter projected polynomial. Row supports and column
 * supports follow first occurrence in the canonical lex polynomial. Terms
 * are contiguous per row; all indices and degrees retain full slong width.
 * Zero-initialize before use; clear owns every array, but no field context. */
typedef struct { slong column, degree; ulong coefficient; } fq_mq_compact_term;
typedef struct {
    slong nvars, nrows, ncols, nterms, degree;
    ulong prime;
    slong *rows, *cols, *offset;
    fq_mq_compact_term *terms;
} fq_mq_compact;

void fq_mq_compact_clear(fq_mq_compact *p);
/* Initialize a conventional polynomial, for repair/fallback or inspection. */
void fq_mq_compact_materialize(fq_mvpoly_t *out, const fq_mq_compact *p,
                              const fq_nmod_ctx_t ctx);
/* Same eligibility and exact projection as compute_fq_det_mq_projected_rect.
 * On rejection leave the zero-initialized output untouched. */
int compute_fq_det_mq_compact(fq_mq_compact *out, fq_mvpoly_t **matrix,
    slong size, const slong *rows, slong nr, const slong *cols, slong nc);
#endif

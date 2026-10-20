/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_MQ_QUOTIENT_H
#define DRSOLVE_MQ_QUOTIENT_H
#include "dr_mpoly.h"
typedef struct {
    double leading, evaluation, interpolation, determinant, total;
    slong core_size, samples;
    const char *reason;
} mq_quotient_stats;
/* Experimental prime-field quadratic resultant. Success initializes result;
 * rejection leaves it untouched. Calls are serialized; samples run in parallel. */
int dr_mq_quotient_resultant(unified_mpoly_struct *result,
    const unified_mpoly_struct *polys, slong nvars, slong npars,
    slong workers, mq_quotient_stats *stats);
/* General positive degrees, prime field only. Rejections leave result untouched.
 * Grouped quotient recurrence; parallel samples, generic ranks and fixed bases.
 * Calls are serialized. workers controls sample concurrency, with no memory budget. */
int dr_degree_quotient_resultant(unified_mpoly_struct *result,
    const unified_mpoly_struct *polys, slong nvars, slong npars,
    slong workers, mq_quotient_stats *stats);
#endif

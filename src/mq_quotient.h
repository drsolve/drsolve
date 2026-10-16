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
#endif

/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef MQ_COEFFICIENT_FILTER_H
#define MQ_COEFFICIENT_FILTER_H
#include <flint/nmod_mpoly.h>
/* Immutable after construction; safe to apply to independent polynomials in
 * parallel. The closure ignores parameter degree. NULL means no projection. */
void *mq_coefficient_filter_create(slong nvars,const slong *rows,slong nr,
    const slong *cols,slong nc,const nmod_mpoly_ctx_t ctx,slong degree_bound);
void mq_coefficient_filter_apply(nmod_mpoly_t poly,const nmod_mpoly_ctx_t ctx,
    const void *filter,int target);
void mq_coefficient_filter_destroy(void *filter);
/* Fix exponent packing before parallel reads. The degree bound was supplied
 * at construction; the filter uses the same packing for fast axis lookups. */
void mq_coefficient_filter_repack(nmod_mpoly_t poly,const nmod_mpoly_ctx_t ctx,
    const void *filter);
#endif

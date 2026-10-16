/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Private ABI used only by the standalone DAG experiment object. */
#ifndef MQ_DAG_HOOK_H
#define MQ_DAG_HOOK_H
#include <stdint.h>
#include <flint/nmod_mpoly.h>
/* On success adopts children and factors; caller frees map and owns output.
 * This is allowed only when the pending root will also be consumed. */
int mq_dag_checkpoint(void *children,int narrow,size_t child_count,
    ulong *factors,size_t shifts,const uint32_t *map,size_t count,
    slong n,void **output);
int mq_dag_root(const void *children,int narrow,size_t child_count,
    const ulong *factors,size_t shifts,const uint32_t *map,
    const ulong *keys,size_t count,slong words,flint_bitcnt_t bits,
    slong n,const nmod_mpoly_ctx_t ctx);
#endif

/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_MQ_FQ_COEFFS_H
#define DRSOLVE_MQ_FQ_COEFFS_H
#include <flint/fq_nmod.h>
#include <string.h>

/* Fixed-width polynomial-basis vectors: d limbs per extension element, no
 * per-element allocation. Borrowed fq views are read-only and never cleared.
 * Small-degree scalar multiplication uses its linear operator, reduced by
 * FLINT with the actual field modulus. Larger degrees keep fq_nmod_mul. */
#define MQ_FQ_OPERATOR_DEGREE 8
static slong mq_fq_length(const ulong *a, slong d)
{
    while (d && !a[d-1]) d--;
    return d;
}
static fq_nmod_struct mq_fq_view(const ulong *a, const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    fq_nmod_struct view;
    view.coeffs = (ulong *)a; view.alloc = d;
    view.length = mq_fq_length(a,d); view.mod = fq_nmod_ctx_modulus(ctx)->mod;
    return view;
}
static void mq_fq_store(ulong *dest, const fq_nmod_t source, slong d)
{
    FLINT_ASSERT(source->length <= d);
    if (source->length) memcpy(dest,source->coeffs,source->length*sizeof(ulong));
    memset(dest+source->length,0,(d-source->length)*sizeof(ulong));
}
typedef struct {
    fq_nmod_struct scalar;
    ulong matrix[MQ_FQ_OPERATOR_DEGREE*MQ_FQ_OPERATOR_DEGREE];
    int lazy;
} mq_fq_operator;

static void mq_fq_operator_prepare(mq_fq_operator *op, const ulong *scalar,
    fq_nmod_t scratch, const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    op->scalar = mq_fq_view(scalar,ctx);
    if (op->scalar.length <= 1 || d > MQ_FQ_OPERATOR_DEGREE) return;
    ulong largest = fq_nmod_ctx_prime(ctx)-1;
    op->lazy = largest <= (UWORD_MAX/(ulong)d)/largest;
    fq_nmod_set(scratch,&op->scalar,ctx);
    for (slong j = 0; j < d; j++) {
        for (slong i = 0; i < d; i++)
            op->matrix[i*d+j] = i < scratch->length ? scratch->coeffs[i] : 0;
        if (j+1 < d) {
            nmod_poly_shift_left(scratch,scratch,1);
            fq_nmod_reduce(scratch,ctx);
        }
    }
}
/* dst and src must be distinct. scratch is private to this worker. */
static void mq_fq_addmul(ulong *dst, const mq_fq_operator *op, const ulong *src,
    int subtract, fq_nmod_t scratch, const fq_nmod_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx);
    nmod_t mod = fq_nmod_ctx_modulus(ctx)->mod;
    if (!op->scalar.length) return;
    if (op->scalar.length == 1) {
        ulong scalar = op->scalar.coeffs[0];
        if (subtract) scalar = nmod_neg(scalar,mod);
        for (slong i = 0; i < d; i++)
            dst[i] = nmod_add(dst[i],nmod_mul(scalar,src[i],mod),mod);
    } else if (d <= MQ_FQ_OPERATOR_DEGREE) {
        for (slong i = 0; i < d; i++) {
            ulong product = 0;
            if (op->lazy) {
                for (slong j = 0; j < d; j++) product += op->matrix[i*d+j]*src[j];
                NMOD_RED(product,product,mod);
            } else {
                for (slong j = 0; j < d; j++)
                    product = nmod_add(product,nmod_mul(op->matrix[i*d+j],src[j],mod),mod);
            }
            dst[i] = subtract ? nmod_sub(dst[i],product,mod) : nmod_add(dst[i],product,mod);
        }
    } else {
        fq_nmod_struct source = mq_fq_view(src,ctx);
        fq_nmod_mul(scratch,&op->scalar,&source,ctx);
        for (slong i = 0; i < scratch->length; i++)
            dst[i] = subtract ? nmod_sub(dst[i],scratch->coeffs[i],mod)
                              : nmod_add(dst[i],scratch->coeffs[i],mod);
    }
}
#endif

/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_ZECH_MPOLY_ACCESS_H
#define DRSOLVE_ZECH_MPOLY_ACCESS_H
#include <flint/fq_zech_mpoly.h>

/* Some FLINT builds declare but do not export the Zech term accessors.
 * Use FLINT's common packed-monomial operations on its native storage. */
static inline void dr_zech_get_exp(slong *exp, const fq_zech_mpoly_t p, slong i,
                                   const fq_zech_mpoly_ctx_t ctx)
{
    slong words = mpoly_words_per_exp(p->bits, ctx->minfo);
    mpoly_get_monomial_si(exp, p->exps + i * words, p->bits, ctx->minfo);
}

static inline void dr_zech_set_exp(fq_zech_mpoly_t p, slong i, const ulong *exp,
                                   const fq_zech_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(mpoly_exp_bits_required_ui(exp, ctx->minfo), ctx->minfo);
    fq_zech_mpoly_fit_bits(p, bits, ctx);
    slong words = mpoly_words_per_exp(p->bits, ctx->minfo);
    mpoly_set_monomial_ui(p->exps + i * words, exp, p->bits, ctx->minfo);
}

static inline void dr_zech_push(fq_zech_mpoly_t p, const fq_zech_t c, const ulong *exp,
                                const fq_zech_mpoly_ctx_t ctx)
{
    fq_zech_mpoly_fit_length(p, p->length + 1, ctx);
    dr_zech_set_exp(p, p->length, exp, ctx);
    p->coeffs[p->length++] = *c;
}

static inline void dr_zech_get_coeff(fq_zech_t c, const fq_zech_mpoly_t p, const ulong *exp,
                                     const fq_zech_mpoly_ctx_t ctx)
{
    slong i = mpoly_monomial_index_ui(p->exps, p->bits, p->length, exp, ctx->minfo);
    if (i < 0)
        fq_zech_zero(c, ctx->fqctx);
    else
        *c = p->coeffs[i];
}

static inline void dr_zech_set_coeff(fq_zech_mpoly_t p, const fq_zech_t c, const ulong *exp,
                                     const fq_zech_mpoly_ctx_t ctx)
{
    slong i = mpoly_monomial_index_ui(p->exps, p->bits, p->length, exp, ctx->minfo);
    if (i < 0) {
        if (!fq_zech_is_zero(c, ctx->fqctx)) {
            dr_zech_push(p, c, exp, ctx);
            fq_zech_mpoly_sort_terms(p, ctx);
        }
    } else if (!fq_zech_is_zero(c, ctx->fqctx)) {
        p->coeffs[i] = *c;
    } else {
        slong tail = p->length - i - 1;
        slong words = mpoly_words_per_exp(p->bits, ctx->minfo);
        memmove(p->coeffs + i, p->coeffs + i + 1, tail * sizeof(*p->coeffs));
        memmove(p->exps + i * words, p->exps + (i + 1) * words, tail * words * sizeof(*p->exps));
        p->length--;
    }
}
#endif

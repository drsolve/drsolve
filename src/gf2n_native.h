/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_GF2N_NATIVE_H
#define DRSOLVE_GF2N_NATIVE_H
#include "gf2n_field.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The first member is a borrowed, read-only FLINT view. Existing algorithms
 * may pass it to FLINT; only FIELD_ID_GF2N_* may cast it back to this object.
 * Clear the wrapper, never the borrowed view. The source context must outlive it.
 * All arithmetic metadata is immutable after construction, including tables. */
typedef struct {
    fq_nmod_ctx_struct fq;
    unsigned n;
    int clmul, fold;
    gf2128_t modulus, mask;
    uint64_t reciprocal; /* low n bits of floor(x^(2n) / modulus), n<=64 */
    gf2128_t (*reduction)[256];
    uint8_t *byte_mul;
    uint8_t (*byte_simd)[32];
    uint16_t *log, *exp, *inverse;
    void *shared; /* Reference-counted immutable precomputation, keyed by modulus. */
} gf2n_native_ctx;

gf2n_native_ctx *gf2n_native_create(const fq_nmod_ctx_t fq);
void gf2n_native_destroy(gf2n_native_ctx *ctx);
gf2128_t gf2n_native_mul(gf2128_t a, gf2128_t b, const gf2n_native_ctx *ctx);
gf2128_t gf2n_native_sqr(gf2128_t a, const gf2n_native_ctx *ctx);
gf2128_t gf2n_native_inv(gf2128_t a, const gf2n_native_ctx *ctx);
gf2128_t gf2n_native_from_fq(const fq_nmod_t a);
void gf2n_native_to_fq(fq_nmod_t out, gf2128_t a, const fq_nmod_ctx_t fq);

#if defined(__PCLMUL__)
/* Polynomial Barrett reduction. Retain the monic terms implicitly, including
 * n=64, so arbitrary dense moduli need no coefficient-by-coefficient lookup.
 * Inputs have degree < n; the approximate quotient leaves at most bit n. */
static inline uint64_t gf2n_native_reduce_word_hw(__m128i p,
                                                const gf2n_native_ctx *c)
{
    unsigned n=c->n;
    uint64_t lo=(uint64_t)_mm_cvtsi128_si64(p);
    uint64_t hi=(uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(p,8));
    uint64_t q=(lo>>(n-1))|(hi<<(65-n));
    p=_mm_clmulepi64_si128(_mm_cvtsi64_si128(q),_mm_cvtsi64_si128(c->reciprocal),0);
    uint64_t qlo=(uint64_t)_mm_cvtsi128_si64(p);
    uint64_t qhi=(uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(p,8));
    q=(q>>1)^(n<63 ? (qlo>>(n+1))|(qhi<<(63-n)) : qhi>>(n-63));
    p=_mm_clmulepi64_si128(_mm_cvtsi64_si128(q),_mm_cvtsi64_si128(c->modulus.low),0);
    uint64_t r=lo^(uint64_t)_mm_cvtsi128_si64(p),carry;
    if(n<64) { r^=q<<n; carry=(r>>n)&1; }
    else carry=(hi^(uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(p,8))^q)&1;
    return (r^(-carry&c->modulus.low))&c->mask.low;
}
static inline uint64_t gf2n_native_mul_word_hw(uint64_t a,uint64_t b,
                                             const gf2n_native_ctx *c)
{
    return gf2n_native_reduce_word_hw(_mm_clmulepi64_si128(
        _mm_cvtsi64_si128(a),_mm_cvtsi64_si128(b),0),c);
}

#endif

static inline uint64_t gf2n_native_mul_word(uint64_t a, uint64_t b,
                                           const gf2n_native_ctx *ctx)
{
    if (ctx->byte_mul) return ctx->byte_mul[(a << ctx->n) | b];
    if (ctx->log) return a && b ? ctx->exp[ctx->log[a]+ctx->log[b]] : 0;
#if defined(__PCLMUL__)
    if(ctx->clmul) return gf2n_native_mul_word_hw(a,b,ctx);
#endif
    return gf2n_native_mul((gf2128_t){a,0}, (gf2128_t){b,0}, ctx).low;
}
static inline gf2128_t gf2n_native_mul_fast(gf2128_t a,gf2128_t b,
                                            const gf2n_native_ctx *ctx)
{
    if(ctx->n<=64) return (gf2128_t){gf2n_native_mul_word(a.low,b.low,ctx),0};
#if defined(__PCLMUL__)
    if (ctx->clmul && ctx->n==128 && !ctx->modulus.high) {
        __m128i x=_mm_set_epi64x(a.high,a.low), y=_mm_set_epi64x(b.high,b.low);
        __m128i lo=_mm_clmulepi64_si128(x,y,0), hi=_mm_clmulepi64_si128(x,y,0x11);
        __m128i cross=_mm_clmulepi64_si128(_mm_cvtsi64_si128(a.low^a.high),
                                         _mm_cvtsi64_si128(b.low^b.high),0);
        cross=_mm_xor_si128(cross,_mm_xor_si128(lo,hi));
        lo=_mm_xor_si128(lo,_mm_slli_si128(cross,8));
        hi=_mm_xor_si128(hi,_mm_srli_si128(cross,8));
        __m128i modulus=_mm_cvtsi64_si128(ctx->modulus.low);
        __m128i t=_mm_clmulepi64_si128(hi,modulus,1);
        hi=_mm_xor_si128(hi,_mm_srli_si128(t,8));
        lo=_mm_xor_si128(lo,_mm_slli_si128(t,8));
        lo=_mm_xor_si128(lo,_mm_clmulepi64_si128(hi,modulus,0));
        return (gf2128_t){(uint64_t)_mm_cvtsi128_si64(lo),
                          (uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(lo,8))};
    }
#endif
    return gf2n_native_mul(a,b,ctx);
}
static inline uint64_t gf2n_native_inv_word(uint64_t a, const gf2n_native_ctx *ctx)
{
    return ctx->inverse ? ctx->inverse[a] : gf2n_native_inv((gf2128_t){a,0},ctx).low;
}
#ifdef __cplusplus
}
#endif
#endif

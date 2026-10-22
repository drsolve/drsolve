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

static inline uint64_t gf2n_native_mul_word(uint64_t a, uint64_t b,
                                           const gf2n_native_ctx *ctx)
{
    if (ctx->byte_mul) return ctx->byte_mul[(a << ctx->n) | b];
    if (ctx->log) return a && b ? ctx->exp[ctx->log[a]+ctx->log[b]] : 0;
    return gf2n_native_mul((gf2128_t){a,0}, (gf2128_t){b,0}, ctx).low;
}
static inline uint64_t gf2n_native_inv_word(uint64_t a, const gf2n_native_ctx *ctx)
{
    return ctx->inverse ? ctx->inverse[a] : gf2n_native_inv((gf2128_t){a,0},ctx).low;
}
#ifdef __cplusplus
}
#endif
#endif

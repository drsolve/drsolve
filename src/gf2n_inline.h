/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_GF2N_INLINE_H
#define DRSOLVE_GF2N_INLINE_H
#include "gf2n_field.h"

/* Keep the arithmetic visible in bulk loops. The public dispatch pointers
 * still select the portable implementation on builds/CPUs without CLMUL,
 * and explicit software overrides remain effective. */
#if defined(__PCLMUL__)
static inline gf232_t gf232_reduce_hw(__m128i p)
{
    uint64_t v=(uint64_t)_mm_cvtsi128_si64(p);
    __m128i modulus=_mm_cvtsi64_si128(0x8299);
    p=_mm_clmulepi64_si128(_mm_cvtsi64_si128(v>>32),modulus,0);
    v=(uint32_t)v^(uint64_t)_mm_cvtsi128_si64(p);
    p=_mm_clmulepi64_si128(_mm_cvtsi64_si128(v>>32),modulus,0);
    return (gf232_t){(uint32_t)v^(uint32_t)_mm_cvtsi128_si64(p)};
}
static inline gf264_t gf264_reduce_hw(__m128i p)
{
    uint64_t lo=(uint64_t)_mm_cvtsi128_si64(p);
    __m128i modulus=_mm_cvtsi64_si128(UINT64_C(0x247f43cb7));
    p=_mm_clmulepi64_si128(_mm_srli_si128(p,8),modulus,0);
    lo^=(uint64_t)_mm_cvtsi128_si64(p);
    p=_mm_clmulepi64_si128(_mm_srli_si128(p,8),modulus,0);
    lo^=(uint64_t)_mm_cvtsi128_si64(p);
    uint64_t hi=(uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(p,8));
    lo^=(-(hi&1)&UINT64_C(0x247f43cb7))^(-((hi>>1)&1)&(UINT64_C(0x247f43cb7)<<1));
    return (gf264_t){lo};
}
static inline gf232_t gf232_product_hw(const gf232_t *a,const gf232_t *b)
{
    return gf232_reduce_hw(_mm_clmulepi64_si128(_mm_cvtsi64_si128(a->value),_mm_cvtsi64_si128(b->value),0));
}
static inline gf264_t gf264_product_hw(const gf264_t *a,const gf264_t *b)
{
    return gf264_reduce_hw(_mm_clmulepi64_si128(_mm_cvtsi64_si128(a->value),_mm_cvtsi64_si128(b->value),0));
}
static inline gf2128_t gf2128_product_hw(const gf2128_t *a,const gf2128_t *b)
{
    __m128i x=_mm_set_epi64x(a->high,a->low),y=_mm_set_epi64x(b->high,b->low);
    __m128i lo=_mm_clmulepi64_si128(x,y,0),hi=_mm_clmulepi64_si128(x,y,0x11);
    __m128i cross=_mm_clmulepi64_si128(_mm_cvtsi64_si128(a->low^a->high),
                                      _mm_cvtsi64_si128(b->low^b->high),0);
    cross=_mm_xor_si128(cross,_mm_xor_si128(lo,hi));
    lo=_mm_xor_si128(lo,_mm_slli_si128(cross,8));
    hi=_mm_xor_si128(hi,_mm_srli_si128(cross,8));
    __m128i modulus=_mm_cvtsi64_si128(0x87);
    __m128i t=_mm_clmulepi64_si128(hi,modulus,1);
    hi=_mm_xor_si128(hi,_mm_srli_si128(t,8));
    lo=_mm_xor_si128(lo,_mm_slli_si128(t,8));
    lo=_mm_xor_si128(lo,_mm_clmulepi64_si128(hi,modulus,0));
    return (gf2128_t){(uint64_t)_mm_cvtsi128_si64(lo),(uint64_t)_mm_cvtsi128_si64(_mm_srli_si128(lo,8))};
}
#endif
static inline gf232_t gf232_mul_inline(const gf232_t *a,const gf232_t *b)
{
#if defined(__PCLMUL__)
    if(gf232_mul==gf232_mul_pclmul) return gf232_product_hw(a,b);
#endif
    return gf232_mul(a,b);
}
static inline gf264_t gf264_mul_inline(const gf264_t *a,const gf264_t *b)
{
#if defined(__PCLMUL__)
    if(gf264_mul==gf264_mul_pclmul) return gf264_product_hw(a,b);
#endif
    return gf264_mul(a,b);
}
static inline gf2128_t gf2128_mul_inline(const gf2128_t *a,const gf2128_t *b)
{
#if defined(__PCLMUL__)
    if(gf2128_mul==gf2128_mul_clmul) return gf2128_product_hw(a,b);
#endif
    return gf2128_mul(a,b);
}
#endif

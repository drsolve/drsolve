/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_QUOTIENT_BYTE_H
#define DRSOLVE_QUOTIENT_BYTE_H
#include <stdint.h>
#include <flint/flint.h>
#if defined(__AVX2__)
#include <immintrin.h>
#elif defined(__SSSE3__)
#include <tmmintrin.h>
#endif

/* Multiplication is F_2-linear: a*b = a*(b&15) XOR a*(b&240).
 * Two 16-entry shuffle tables replace one dependent scalar lookup per byte.
 * GF(16) uses the same kernel with zero high nibbles. */
static void quotient_byte_tables(uint8_t tables[256][32], int order,
                                 const uint8_t *(*scalar_row)(uint8_t))
{
    for (int a=0;a<order;a++) {
        const uint8_t *row=scalar_row(a);
        for (int b=0;b<16;b++) {
            tables[a][b]=row[b];
            tables[a][16+b]=order==256 ? row[b<<4] : 0;
        }
    }
}

static void quotient_byte_addmul(uint8_t *v, const uint8_t *w, slong n,
                                 const uint8_t table[32])
{
    slong j=0;
#if defined(__AVX2__)
    __m256i lo=_mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)table));
    __m256i hi=_mm256_broadcastsi128_si256(_mm_loadu_si128((const __m128i *)(table+16)));
    __m256i mask=_mm256_set1_epi8(15);
    for (;j+32<=n;j+=32) {
        __m256i b=_mm256_loadu_si256((const __m256i *)(w+j));
        __m256i p=_mm256_xor_si256(_mm256_shuffle_epi8(lo,_mm256_and_si256(b,mask)),
            _mm256_shuffle_epi8(hi,_mm256_and_si256(_mm256_srli_epi16(b,4),mask)));
        _mm256_storeu_si256((__m256i *)(v+j),
            _mm256_xor_si256(_mm256_loadu_si256((const __m256i *)(v+j)),p));
    }
#elif defined(__SSSE3__)
    __m128i lo=_mm_loadu_si128((const __m128i *)table);
    __m128i hi=_mm_loadu_si128((const __m128i *)(table+16));
    __m128i mask=_mm_set1_epi8(15);
    for (;j+16<=n;j+=16) {
        __m128i b=_mm_loadu_si128((const __m128i *)(w+j));
        __m128i p=_mm_xor_si128(_mm_shuffle_epi8(lo,_mm_and_si128(b,mask)),
            _mm_shuffle_epi8(hi,_mm_and_si128(_mm_srli_epi16(b,4),mask)));
        _mm_storeu_si128((__m128i *)(v+j),
            _mm_xor_si128(_mm_loadu_si128((const __m128i *)(v+j)),p));
    }
#endif
    for (;j<n;j++) v[j]^=table[w[j]&15]^table[16+(w[j]>>4)];
}
#endif

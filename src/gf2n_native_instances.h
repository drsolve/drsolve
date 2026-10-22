/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Define MQGF_PREFIX, MQGF_CONTEXT and MQGF_TEMPLATE before including. */
#define GFN_JOIN_(a,b) a##b
#define GFN_JOIN(a,b) GFN_JOIN_(a,b)

#define MQGF_TYPE uint8_t
#define MQGF_NAME GFN_JOIN(MQGF_PREFIX,_gf2n8)
#define MQGF_MEMBER binary8
#define MQGF_ZERO(a) ((a)==0)
#define MQGF_ADD(a,b) ((a)^(b))
#ifdef MQGF_USE_BYTE
#define MQGF_BYTE
#endif
#define MQGF_PREPARE(a) const gf2n_native_ctx *nc=(const gf2n_native_ctx *)(MQGF_CONTEXT); const uint8_t *scalar_row=nc->byte_mul+((size_t)(a)<<nc->n)
#define MQGF_MUL(a,b) scalar_row[b]
#include MQGF_TEMPLATE

#define MQGF_TYPE uint16_t
#define MQGF_NAME GFN_JOIN(MQGF_PREFIX,_gf2n16)
#define MQGF_MEMBER binary16
#define MQGF_ZERO(a) ((a)==0)
#define MQGF_ADD(a,b) ((a)^(b))
#define MQGF_PREPARE(a) const gf2n_native_ctx *nc=(const gf2n_native_ctx *)(MQGF_CONTEXT)
#define MQGF_MUL(a,b) gf2n_native_mul_word((a),(b),nc)
#include MQGF_TEMPLATE

#define MQGF_TYPE uint64_t
#define MQGF_NAME GFN_JOIN(MQGF_PREFIX,_gf2n64)
#define MQGF_MEMBER binary64
#define MQGF_ZERO(a) ((a)==0)
#define MQGF_ADD(a,b) ((a)^(b))
#define MQGF_PREPARE(a) const gf2n_native_ctx *nc=(const gf2n_native_ctx *)(MQGF_CONTEXT)
#define MQGF_MUL(a,b) gf2n_native_mul_word((a),(b),nc)
#include MQGF_TEMPLATE

#define MQGF_TYPE gf2128_t
#define MQGF_NAME GFN_JOIN(MQGF_PREFIX,_gf2n128)
#define MQGF_MEMBER binary128
#define MQGF_ZERO(a) (!((a).low|(a).high))
#define MQGF_ADD(a,b) (gf2128_t){(a).low^(b).low,(a).high^(b).high}
#define MQGF_PREPARE(a) const gf2n_native_ctx *nc=(const gf2n_native_ctx *)(MQGF_CONTEXT)
#define MQGF_MUL(a,b) gf2n_native_mul((a),(b),nc)
#include MQGF_TEMPLATE

#undef GFN_JOIN
#undef GFN_JOIN_
#undef MQGF_PREFIX
#undef MQGF_CONTEXT
#undef MQGF_TEMPLATE
#undef MQGF_USE_BYTE

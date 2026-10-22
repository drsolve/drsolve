/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Compact coefficient instances of the same matrix algorithms. */

#define UM_NAME(x) b8_##x
#define UM_E uint8_t
#define UM_BYTE 1
#define UM_INIT(x) (*(x)=0)
#define UM_CLEAR(x) ((void)0)
#define UM_ZERO(x) (*(x)=0)
#define UM_ONE(x) (*(x)=1)
#define UM_SET(x,y) (*(x)=*(y))
#define UM_ISZERO(x) (*(x)==0)
#define UM_ADD(x,y,z) (*(x)=(*(y)^*(z)))
#define UM_SUB(x,y,z) UM_ADD(x,y,z)
#define UM_NEG(x,y) UM_SET(x,y)
#define UM_MUL(x,y,z) (*(x)=((const gf2n_native_ctx *)ctx->ctx.fq_ctx)->byte_mul[((size_t)*(y)<<((const gf2n_native_ctx *)ctx->ctx.fq_ctx)->n)|*(z)])
#define UM_INV(x,y) (*(x)=((const gf2n_native_ctx *)ctx->ctx.fq_ctx)->inverse[*(y)])
#include "unified_mat_poly_template.h"
#undef UM_NAME
#undef UM_E
#undef UM_BYTE
#undef UM_INIT
#undef UM_CLEAR
#undef UM_ZERO
#undef UM_ONE
#undef UM_SET
#undef UM_ISZERO
#undef UM_ADD
#undef UM_SUB
#undef UM_NEG
#undef UM_MUL
#undef UM_INV

#define UM_NAME(x) u16_##x
#define UM_E uint16_t
#define UM_BYTE 0
#define UM_INIT(x) (*(x)=0)
#define UM_CLEAR(x) ((void)0)
#define UM_ZERO(x) (*(x)=0)
#define UM_ONE(x) (*(x)=1)
#define UM_SET(x,y) (*(x)=*(y))
#define UM_ISZERO(x) (*(x)==0)
#define UM_ADD(x,y,z) (*(x)=(*(y)^*(z)))
#define UM_SUB(x,y,z) UM_ADD(x,y,z)
#define UM_NEG(x,y) UM_SET(x,y)
#define UM_MUL(x,y,z) (*(x)=gf216_mul(*(y),*(z)))
#define UM_INV(x,y) (*(x)=gf216_inv(*(y)))
#include "unified_mat_poly_template.h"
#undef UM_NAME
#undef UM_E
#undef UM_BYTE
#undef UM_INIT
#undef UM_CLEAR
#undef UM_ZERO
#undef UM_ONE
#undef UM_SET
#undef UM_ISZERO
#undef UM_ADD
#undef UM_SUB
#undef UM_NEG
#undef UM_MUL
#undef UM_INV

#define UM_NAME(x) u32_##x
#define UM_E gf232_t
#define UM_BYTE 0
#define UM_INIT(x) (*(x)=(gf232_t){0})
#define UM_CLEAR(x) ((void)0)
#define UM_ZERO(x) (*(x)=(gf232_t){0})
#define UM_ONE(x) (*(x)=(gf232_t){1})
#define UM_SET(x,y) (*(x)=*(y))
#define UM_ISZERO(x) ((x)->value==0)
#define UM_ADD(x,y,z) (*(x)=(gf232_t){(y)->value^(z)->value})
#define UM_SUB(x,y,z) UM_ADD(x,y,z)
#define UM_NEG(x,y) UM_SET(x,y)
#define UM_MUL(x,y,z) (*(x)=gf232_mul(y,z))
#define UM_INV(x,y) (*(x)=gf232_inv(y))
#include "unified_mat_poly_template.h"
#undef UM_NAME
#undef UM_E
#undef UM_BYTE
#undef UM_INIT
#undef UM_CLEAR
#undef UM_ZERO
#undef UM_ONE
#undef UM_SET
#undef UM_ISZERO
#undef UM_ADD
#undef UM_SUB
#undef UM_NEG
#undef UM_MUL
#undef UM_INV

#define UM_NAME(x) u64_##x
#define UM_E gf264_t
#define UM_BYTE 0
#define UM_INIT(x) (*(x)=(gf264_t){0})
#define UM_CLEAR(x) ((void)0)
#define UM_ZERO(x) (*(x)=(gf264_t){0})
#define UM_ONE(x) (*(x)=(gf264_t){1})
#define UM_SET(x,y) (*(x)=*(y))
#define UM_ISZERO(x) ((x)->value==0)
#define UM_ADD(x,y,z) (*(x)=(gf264_t){(y)->value^(z)->value})
#define UM_SUB(x,y,z) UM_ADD(x,y,z)
#define UM_NEG(x,y) UM_SET(x,y)
#define UM_MUL(x,y,z) (*(x)=gf264_mul(y,z))
#define UM_INV(x,y) (*(x)=gf264_inv(y))
#include "unified_mat_poly_template.h"
#undef UM_NAME
#undef UM_E
#undef UM_BYTE
#undef UM_INIT
#undef UM_CLEAR
#undef UM_ZERO
#undef UM_ONE
#undef UM_SET
#undef UM_ISZERO
#undef UM_ADD
#undef UM_SUB
#undef UM_NEG
#undef UM_MUL
#undef UM_INV

#define UM_NAME(x) u128_##x
#define UM_E gf2128_t
#define UM_BYTE 0
#define UM_INIT(x) (*(x)=(gf2128_t){0,0})
#define UM_CLEAR(x) ((void)0)
#define UM_ZERO(x) (*(x)=(gf2128_t){0,0})
#define UM_ONE(x) (*(x)=(gf2128_t){1,0})
#define UM_SET(x,y) (*(x)=*(y))
#define UM_ISZERO(x) !((x)->low|(x)->high)
#define UM_ADD(x,y,z) (*(x)=(gf2128_t){(y)->low^(z)->low,(y)->high^(z)->high})
#define UM_SUB(x,y,z) UM_ADD(x,y,z)
#define UM_NEG(x,y) UM_SET(x,y)
#define UM_MUL(x,y,z) (*(x)=gf2128_mul(y,z))
#define UM_INV(x,y) (*(x)=gf2128_inv(y))
#include "unified_mat_poly_template.h"
#undef UM_NAME
#undef UM_E
#undef UM_BYTE
#undef UM_INIT
#undef UM_CLEAR
#undef UM_ZERO
#undef UM_ONE
#undef UM_SET
#undef UM_ISZERO
#undef UM_ADD
#undef UM_SUB
#undef UM_NEG
#undef UM_MUL
#undef UM_INV

#define UM_NAME(x) gf2n16_##x
#define UM_E uint16_t
#define UM_BYTE 0
#define UM_INIT(x) (*(x)=0)
#define UM_CLEAR(x) ((void)0)
#define UM_ZERO(x) (*(x)=0)
#define UM_ONE(x) (*(x)=1)
#define UM_SET(x,y) (*(x)=*(y))
#define UM_ISZERO(x) (*(x)==0)
#define UM_ADD(x,y,z) (*(x)=(*(y)^*(z)))
#define UM_SUB(x,y,z) UM_ADD(x,y,z)
#define UM_NEG(x,y) UM_SET(x,y)
#define UM_MUL(x,y,z) (*(x)=gf2n_native_mul_word(*(y),*(z),(const gf2n_native_ctx *)ctx->ctx.fq_ctx))
#define UM_INV(x,y) (*(x)=gf2n_native_inv_word(*(y),(const gf2n_native_ctx *)ctx->ctx.fq_ctx))
#include "unified_mat_poly_template.h"
#undef UM_NAME
#undef UM_E
#undef UM_BYTE
#undef UM_INIT
#undef UM_CLEAR
#undef UM_ZERO
#undef UM_ONE
#undef UM_SET
#undef UM_ISZERO
#undef UM_ADD
#undef UM_SUB
#undef UM_NEG
#undef UM_MUL
#undef UM_INV

#define UM_NAME(x) gf2n64_##x
#define UM_E uint64_t
#define UM_BYTE 0
#define UM_INIT(x) (*(x)=0)
#define UM_CLEAR(x) ((void)0)
#define UM_ZERO(x) (*(x)=0)
#define UM_ONE(x) (*(x)=1)
#define UM_SET(x,y) (*(x)=*(y))
#define UM_ISZERO(x) (*(x)==0)
#define UM_ADD(x,y,z) (*(x)=(*(y)^*(z)))
#define UM_SUB(x,y,z) UM_ADD(x,y,z)
#define UM_NEG(x,y) UM_SET(x,y)
#define UM_MUL(x,y,z) (*(x)=gf2n_native_mul_word(*(y),*(z),(const gf2n_native_ctx *)ctx->ctx.fq_ctx))
#define UM_INV(x,y) (*(x)=gf2n_native_inv_word(*(y),(const gf2n_native_ctx *)ctx->ctx.fq_ctx))
#include "unified_mat_poly_template.h"
#undef UM_NAME
#undef UM_E
#undef UM_BYTE
#undef UM_INIT
#undef UM_CLEAR
#undef UM_ZERO
#undef UM_ONE
#undef UM_SET
#undef UM_ISZERO
#undef UM_ADD
#undef UM_SUB
#undef UM_NEG
#undef UM_MUL
#undef UM_INV

#define UM_NAME(x) gf2n128_##x
#define UM_E gf2128_t
#define UM_BYTE 0
#define UM_INIT(x) (*(x)=(gf2128_t){0,0})
#define UM_CLEAR(x) ((void)0)
#define UM_ZERO(x) (*(x)=(gf2128_t){0,0})
#define UM_ONE(x) (*(x)=(gf2128_t){1,0})
#define UM_SET(x,y) (*(x)=*(y))
#define UM_ISZERO(x) !((x)->low|(x)->high)
#define UM_ADD(x,y,z) (*(x)=(gf2128_t){(y)->low^(z)->low,(y)->high^(z)->high})
#define UM_SUB(x,y,z) UM_ADD(x,y,z)
#define UM_NEG(x,y) UM_SET(x,y)
#define UM_MUL(x,y,z) (*(x)=gf2n_native_mul(*(y),*(z),(const gf2n_native_ctx *)ctx->ctx.fq_ctx))
#define UM_INV(x,y) (*(x)=gf2n_native_inv(*(y),(const gf2n_native_ctx *)ctx->ctx.fq_ctx))
#include "unified_mat_poly_template.h"
#undef UM_NAME
#undef UM_E
#undef UM_BYTE
#undef UM_INIT
#undef UM_CLEAR
#undef UM_ZERO
#undef UM_ONE
#undef UM_SET
#undef UM_ISZERO
#undef UM_ADD
#undef UM_SUB
#undef UM_NEG
#undef UM_MUL
#undef UM_INV

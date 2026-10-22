/* SPDX-License-Identifier: GPL-2.0-or-later */
#define DP_NAME(x) minor_gf28_##x
#define DP_POLY gf28_mpoly_struct
#define DP_POLY_T gf28_mpoly_t
#define DP_CTX_T gf28_mpoly_ctx_t
#define DP_E uint8_t
#define DP_SET gf28_mpoly_set
#define DP_INIT gf28_mpoly_init
#define DP_CLEAR gf28_mpoly_clear
#define DP_CTX_INIT gf28_mpoly_ctx_init
#define DP_CTX_CLEAR gf28_mpoly_ctx_clear
#define DP_MUL(o,a,b,c) (result->field_id==FIELD_ID_GF24 ? gf24_mpoly_mul(o,a,b,c) : gf28_mpoly_mul(o,a,b,c))
#define DP_IMPORT(o,a,c,nc) (result->field_id==FIELD_ID_GF24 ? fq_nmod_mpoly_to_gf24_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c)) : fq_nmod_mpoly_to_gf28_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c)))
#define DP_EXPORT(o,a,c,nc) (result->field_id==FIELD_ID_GF24 ? gf24_mpoly_to_fq_nmod_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c)) : gf28_mpoly_to_fq_nmod_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c)))
#define DP_ADD(a,b) ((a)^(b))
#define DP_ZERO(a) (!(a))
#include "gf2n_minor_dp_template.h"

#define DP_NAME(x) minor_gf216_##x
#define DP_POLY gf216_mpoly_struct
#define DP_POLY_T gf216_mpoly_t
#define DP_CTX_T gf216_mpoly_ctx_t
#define DP_E uint16_t
#define DP_SET gf216_mpoly_set
#define DP_INIT gf216_mpoly_init
#define DP_CLEAR gf216_mpoly_clear
#define DP_CTX_INIT gf216_mpoly_ctx_init
#define DP_CTX_CLEAR gf216_mpoly_ctx_clear
#define DP_MUL gf216_mpoly_mul
#define DP_IMPORT(o,a,c,nc) fq_nmod_mpoly_to_gf216_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c))
#define DP_EXPORT(o,a,c,nc) gf216_mpoly_to_fq_nmod_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c))
#define DP_ADD(a,b) ((a)^(b))
#define DP_ZERO(a) (!(a))
#include "gf2n_minor_dp_template.h"

#define DP_NAME(x) minor_gf232_##x
#define DP_POLY gf232_mpoly_struct
#define DP_POLY_T gf232_mpoly_t
#define DP_CTX_T gf232_mpoly_ctx_t
#define DP_E gf232_t
#define DP_SET gf232_mpoly_set
#define DP_INIT gf232_mpoly_init
#define DP_CLEAR gf232_mpoly_clear
#define DP_CTX_INIT gf232_mpoly_ctx_init
#define DP_CTX_CLEAR gf232_mpoly_ctx_clear
#define DP_MUL gf232_mpoly_mul
#define DP_IMPORT(o,a,c,nc) fq_nmod_mpoly_to_gf232_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c))
#define DP_EXPORT(o,a,c,nc) gf232_mpoly_to_fq_nmod_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c))
#define DP_ADD(a,b) gf232_add(&(a),&(b))
#define DP_ZERO(a) gf232_is_zero(&(a))
#include "gf2n_minor_dp_template.h"

#define DP_NAME(x) minor_gf264_##x
#define DP_POLY gf264_mpoly_struct
#define DP_POLY_T gf264_mpoly_t
#define DP_CTX_T gf264_mpoly_ctx_t
#define DP_E gf264_t
#define DP_SET gf264_mpoly_set
#define DP_INIT gf264_mpoly_init
#define DP_CLEAR gf264_mpoly_clear
#define DP_CTX_INIT gf264_mpoly_ctx_init
#define DP_CTX_CLEAR gf264_mpoly_ctx_clear
#define DP_MUL gf264_mpoly_mul
#define DP_IMPORT(o,a,c,nc) fq_nmod_mpoly_to_gf264_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c))
#define DP_EXPORT(o,a,c,nc) gf264_mpoly_to_fq_nmod_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c))
#define DP_ADD(a,b) gf264_add(&(a),&(b))
#define DP_ZERO(a) gf264_is_zero(&(a))
#include "gf2n_minor_dp_template.h"

#define DP_NAME(x) minor_gf2128_##x
#define DP_POLY gf2128_mpoly_struct
#define DP_POLY_T gf2128_mpoly_t
#define DP_CTX_T gf2128_mpoly_ctx_t
#define DP_E gf2128_t
#define DP_SET gf2128_mpoly_set
#define DP_INIT gf2128_mpoly_init
#define DP_CLEAR gf2128_mpoly_clear
#define DP_CTX_INIT gf2128_mpoly_ctx_init
#define DP_CTX_CLEAR gf2128_mpoly_ctx_clear
#define DP_MUL gf2128_mpoly_mul
#define DP_IMPORT(o,a,c,nc) fq_nmod_mpoly_to_gf2128_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c))
#define DP_EXPORT(o,a,c,nc) gf2128_mpoly_to_fq_nmod_mpoly(o,a,(c)->field_ctx->ctx.fq_ctx,GET_FQ_CTX(c))
#define DP_ADD(a,b) gf2128_add(&(a),&(b))
#define DP_ZERO(a) gf2128_is_zero(&(a))
#include "gf2n_minor_dp_template.h"

#include "gf2n_minor_dynamic.h"

static int compute_native_mpoly_det_dp(unified_mpoly_t result,
    unified_mpoly_t **matrix, slong size, unified_mpoly_ctx_t ctx,
    int use_parallel, const ulong choose[FLINT_BITS][FLINT_BITS])
{
#define NATIVE_DP(n) minor_gf2##n##_det(result,matrix,size,ctx,use_parallel,choose)
    switch(result->field_id) {
    case FIELD_ID_GF24: case FIELD_ID_GF28: return NATIVE_DP(8);
    case FIELD_ID_GF216: return NATIVE_DP(16);
    case FIELD_ID_GF232: return NATIVE_DP(32);
    case FIELD_ID_GF264: return NATIVE_DP(64);
    case FIELD_ID_GF2128: return NATIVE_DP(128);
    case FIELD_ID_GF2N8: case FIELD_ID_GF2N16:
    case FIELD_ID_GF2N64: case FIELD_ID_GF2N128:
        return minor_dynamic_det(result,matrix,size,ctx,use_parallel,choose);
    default: return 0;
    }
#undef NATIVE_DP
}

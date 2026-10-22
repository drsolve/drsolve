/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_mat_poly.h"
#include "unified_byte.h"
#include <assert.h>
#include <limits.h>
#ifdef _OPENMP
#include <omp.h>
#endif

extern double get_wall_time(void);
extern int g_dixon_verbose_level;

static const void *unified_field_data(const field_ctx_t *ctx)
{
    if (ctx->field_id == FIELD_ID_NMOD) return &ctx->ctx.nmod_ctx;
    if (ctx->field_id == FIELD_ID_FQ_ZECH) return ctx->ctx.zech_ctx;
    return ctx->ctx.fq_ctx;
}

static const uint8_t (*unified_byte_multipliers(const field_ctx_t *ctx))[32]
{
    if(ctx->field_id==FIELD_ID_GF2N8)
        return ((const gf2n_native_ctx *)ctx->ctx.fq_ctx)->byte_simd;
    /* Native byte arithmetic uses fixed standard fields; fq modulus
     * isomorphisms are applied only at the conversion boundary. */
    static _Thread_local uint8_t tables[2][256][32];
    static _Thread_local int ready[2];
    int which=ctx->field_id==FIELD_ID_GF28;
    if(!ready[which]) {
        unified_byte_tables(tables[which],which ? 256 : 16,
                            which ? gf28_get_scalar_row : gf24_get_scalar_row);
        ready[which]=1;
    }
    return tables[which];
}

int unified_mat_poly_is_byte(const unified_mat_poly_t a)
{
    return a->ctx->field_id == FIELD_ID_GF24 || a->ctx->field_id == FIELD_ID_GF28 || a->ctx->field_id == FIELD_ID_GF2N8;
}

void unified_mat_poly_init(unified_mat_poly_t a, slong r, slong c, field_ctx_t *ctx)
{
    assert(r >= 0 && c >= 0);
    *a = (unified_mat_poly_struct){.r=r, .c=c, .ctx=ctx};
}

/* Instantiate complete operations, not individual coefficient callbacks.
 * The byte convolution selects a multiplication table once per scalar. */
#define UM_NAME(x) byte_##x
#define UM_E uint8_t
#define UM_BYTE 1
#define UM_INIT(x) (*(x)=0)
#define UM_CLEAR(x) ((void)0)
#define UM_ZERO(x) (*(x)=0)
#define UM_ONE(x) (*(x)=1)
#define UM_SET(x,y) (*(x)=*(y))
#define UM_ISZERO(x) (*(x)==0)
#define UM_ADD(x,y,z) (*(x)=*(y)^*(z))
#define UM_SUB(x,y,z) UM_ADD(x,y,z)
#define UM_NEG(x,y) UM_SET(x,y)
#define UM_MUL(x,y,z) (*(x)=ctx->field_id==FIELD_ID_GF24 ? gf24_mul(*(y),*(z)) : gf28_mul(*(y),*(z)))
#define UM_INV(x,y) (*(x)=ctx->field_id==FIELD_ID_GF24 ? gf24_inv(*(y)) : gf28_inv(*(y)))
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

#define UM_NAME(x) generic_##x
#define UM_E field_elem_u
#define UM_BYTE 0
#define UM_INIT(x) field_init_elem(x,ctx->field_id,unified_field_data(ctx))
#define UM_CLEAR(x) field_clear_elem(x,ctx->field_id,unified_field_data(ctx))
#define UM_ZERO(x) field_set_zero(x,ctx->field_id,unified_field_data(ctx))
#define UM_ONE(x) field_set_one(x,ctx->field_id,unified_field_data(ctx))
#define UM_SET(x,y) field_set_elem(x,y,ctx->field_id,unified_field_data(ctx))
#define UM_ISZERO(x) field_is_zero(x,ctx->field_id,unified_field_data(ctx))
#define UM_ADD(x,y,z) field_add(x,y,z,ctx->field_id,unified_field_data(ctx))
#define UM_SUB(x,y,z) field_sub(x,y,z,ctx->field_id,unified_field_data(ctx))
#define UM_NEG(x,y) field_neg(x,y,ctx->field_id,unified_field_data(ctx))
#define UM_MUL(x,y,z) field_mul(x,y,z,ctx->field_id,unified_field_data(ctx))
#define UM_INV(x,y) field_inv(x,y,ctx->field_id,unified_field_data(ctx))
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

#include "unified_mat_poly_native.h"

/* Dispatch once per complete operation, outside coefficient loops. */
#define UM_DISPATCH(a,op,...) do { \
    switch((a)->ctx->field_id) { \
    case FIELD_ID_GF24: case FIELD_ID_GF28: byte_##op(__VA_ARGS__); break; \
    case FIELD_ID_GF2N8: b8_##op(__VA_ARGS__); break; \
    case FIELD_ID_GF216: u16_##op(__VA_ARGS__); break; \
    case FIELD_ID_GF232: u32_##op(__VA_ARGS__); break; \
    case FIELD_ID_GF264: u64_##op(__VA_ARGS__); break; \
    case FIELD_ID_GF2128: u128_##op(__VA_ARGS__); break; \
    case FIELD_ID_GF2N16: gf2n16_##op(__VA_ARGS__); break; \
    case FIELD_ID_GF2N64: gf2n64_##op(__VA_ARGS__); break; \
    case FIELD_ID_GF2N128: gf2n128_##op(__VA_ARGS__); break; \
    default: generic_##op(__VA_ARGS__); break; \
    } \
} while(0)

void unified_mat_poly_clear(unified_mat_poly_t a)
{
    UM_DISPATCH(a,clear,a);
}

static void unified_mat_poly_resize(unified_mat_poly_t a, slong length)
{
    UM_DISPATCH(a,resize,a,length);
}

static void unified_mat_poly_put(unified_mat_poly_t a, slong index, const field_elem_u *v)
{
    switch(a->ctx->field_id) {
    case FIELD_ID_GF216: ((uint16_t *)a->data)[index]=v->gf216; break;
    case FIELD_ID_GF232: ((gf232_t *)a->data)[index]=v->gf232; break;
    case FIELD_ID_GF264: ((gf264_t *)a->data)[index]=v->gf264; break;
    case FIELD_ID_GF2128: ((gf2128_t *)a->data)[index]=v->gf2128; break;
    case FIELD_ID_GF2N16: ((uint16_t *)a->data)[index]=v->binary16; break;
    case FIELD_ID_GF2N64: ((uint64_t *)a->data)[index]=v->binary64; break;
    case FIELD_ID_GF2N128: ((gf2128_t *)a->data)[index]=v->binary128; break;
    case FIELD_ID_GF24: ((uint8_t *)a->data)[index]=v->gf24; break;
    case FIELD_ID_GF28: ((uint8_t *)a->data)[index]=v->gf28; break;
    case FIELD_ID_GF2N8: ((uint8_t *)a->data)[index]=v->binary8; break;
    default: field_set_elem((field_elem_u *)a->data+index,v,a->ctx->field_id,unified_field_data(a->ctx)); break;
    }
}
static void unified_mat_poly_get(field_elem_u *v, const unified_mat_poly_t a, slong index)
{
    switch(a->ctx->field_id) {
    case FIELD_ID_GF216: v->gf216=((uint16_t *)a->data)[index]; break;
    case FIELD_ID_GF232: v->gf232=((gf232_t *)a->data)[index]; break;
    case FIELD_ID_GF264: v->gf264=((gf264_t *)a->data)[index]; break;
    case FIELD_ID_GF2128: v->gf2128=((gf2128_t *)a->data)[index]; break;
    case FIELD_ID_GF2N16: v->binary16=((uint16_t *)a->data)[index]; break;
    case FIELD_ID_GF2N64: v->binary64=((uint64_t *)a->data)[index]; break;
    case FIELD_ID_GF2N128: v->binary128=((gf2128_t *)a->data)[index]; break;
    case FIELD_ID_GF24: v->gf24=((uint8_t *)a->data)[index]; break;
    case FIELD_ID_GF28: v->gf28=((uint8_t *)a->data)[index]; break;
    case FIELD_ID_GF2N8: v->binary8=((uint8_t *)a->data)[index]; break;
    default: field_set_elem(v,(field_elem_u *)a->data+index,a->ctx->field_id,unified_field_data(a->ctx)); break;
    }
}

void unified_mat_poly_from_fq(unified_mat_poly_t a, const fq_nmod_poly_mat_t b,
                              const fq_nmod_ctx_t ctx)
{
    assert(a->r==b->r && a->c==b->c);
    (void)ctx; /* Conversions use the borrowed context in a->ctx. */
    slong len=0;
    for(slong i=0;i<a->r;i++) for(slong j=0;j<a->c;j++)
        len=FLINT_MAX(len,fq_nmod_poly_mat_entry(b,i,j)->length);
    unified_mat_poly_resize(a,0); unified_mat_poly_resize(a,len);
    field_elem_u v; field_init_elem(&v,a->ctx->field_id,unified_field_data(a->ctx));
    for(slong i=0;i<a->r;i++) for(slong j=0;j<a->c;j++) {
        const fq_nmod_poly_struct *p=fq_nmod_poly_mat_entry(b,i,j);
        for(slong k=0;k<p->length;k++) {
            fq_nmod_to_field_elem(&v,p->coeffs+k,a->ctx);
            unified_mat_poly_put(a,(k*a->r+i)*a->c+j,&v);
        }
    }
    field_clear_elem(&v,a->ctx->field_id,unified_field_data(a->ctx));
}

void unified_mat_poly_to_fq(fq_nmod_poly_mat_t a, const unified_mat_poly_t b,
                            const fq_nmod_ctx_t ctx)
{
    assert(a->r==b->r && a->c==b->c);
    field_elem_u v; field_init_elem(&v,b->ctx->field_id,unified_field_data(b->ctx));
    fq_nmod_t x; fq_nmod_init(x,ctx);
    for(slong i=0;i<b->r;i++) for(slong j=0;j<b->c;j++) {
        fq_nmod_poly_struct *p=fq_nmod_poly_mat_entry(a,i,j);
        fq_nmod_poly_zero(p,ctx);
        for(slong k=0;k<b->length;k++) {
            unified_mat_poly_get(&v,b,(k*b->r+i)*b->c+j);
            field_elem_to_fq_nmod(x,&v,b->ctx);
            fq_nmod_poly_set_coeff(p,k,x,ctx);
        }
    }
    fq_nmod_clear(x,ctx);
    field_clear_elem(&v,b->ctx->field_id,unified_field_data(b->ctx));
}

void unified_mat_poly_from_unified(unified_mat_poly_t a, const unified_poly_mat_t b)
{
    assert(a->r==b->r && a->c==b->c && a->ctx==b->ctx);
    slong len=0;
    for(slong i=0;i<a->r;i++) for(slong j=0;j<a->c;j++) len=FLINT_MAX(len,b->rows[i][j].length);
    unified_mat_poly_resize(a,0); unified_mat_poly_resize(a,len);
    for(slong i=0;i<a->r;i++) for(slong j=0;j<a->c;j++)
        for(slong k=0;k<b->rows[i][j].length;k++)
            unified_mat_poly_put(a,(k*a->r+i)*a->c+j,b->rows[i][j].coeffs+k);
}

void unified_mat_poly_to_unified(unified_poly_mat_t a, const unified_mat_poly_t b)
{
    assert(a->r==b->r && a->c==b->c && a->ctx==b->ctx);
    for(slong i=0;i<a->r;i++) for(slong j=0;j<a->c;j++) {
        unified_poly_struct *p=a->rows[i]+j;
        unified_poly_fit_length(p,b->length);
        for(slong k=0;k<b->length;k++) unified_mat_poly_get(p->coeffs+k,b,(k*b->r+i)*b->c+j);
        p->length=b->length; unified_poly_normalise(p);
    }
}

void unified_mat_poly_mul_window(unified_mat_poly_t out, const unified_mat_poly_t a,
                                 const unified_mat_poly_t b, slong start, slong length)
{
    assert(start>=0 && length>=0 && a->c==b->r && out->r==a->r && out->c==b->c);
    assert(out->ctx==a->ctx && a->ctx==b->ctx);
    UM_DISPATCH(a,mul_window,out,a,b,start,length);
}

void unified_mat_poly_mul(unified_mat_poly_t out, const unified_mat_poly_t a,
                          const unified_mat_poly_t b)
{
    unified_mat_poly_mul_window(out,a,b,0,a->length && b->length ? a->length+b->length-1 : 0);
}

void unified_mat_poly_pmbasis(unified_mat_poly_t out, slong *shift,
                              const unified_mat_poly_t a, slong order)
{
    assert(order>=0 && out->r==a->r && out->c==a->r && out->ctx==a->ctx);
    UM_DISPATCH(a,pmbasis,out,shift,a,order);
}

int unified_mat_poly_kernel(unified_mat_poly_t out, const unified_mat_poly_t a)
{
    assert(a->r>=a->c && out->r==a->r-a->c && out->c==a->r && out->ctx==a->ctx);
    unified_mat_poly_t transpose;unified_mat_poly_init(transpose,a->c,a->r,a->ctx);
    int ok;
    switch(a->ctx->field_id) {
    case FIELD_ID_GF2N8:
        b8_select(transpose,a,NULL,NULL,1);
        ok=b8_weak_popov(transpose)>=0 && b8_kernel(out,a); break;
    case FIELD_ID_GF24: case FIELD_ID_GF28:
        byte_select(transpose,a,NULL,NULL,1);
        ok=byte_weak_popov(transpose)>=0 && byte_kernel(out,a); break;
    case FIELD_ID_GF216:
        u16_select(transpose,a,NULL,NULL,1);
        ok=u16_weak_popov(transpose)>=0 && u16_kernel(out,a); break;
    case FIELD_ID_GF232:
        u32_select(transpose,a,NULL,NULL,1);
        ok=u32_weak_popov(transpose)>=0 && u32_kernel(out,a); break;
    case FIELD_ID_GF264:
        u64_select(transpose,a,NULL,NULL,1);
        ok=u64_weak_popov(transpose)>=0 && u64_kernel(out,a); break;
    case FIELD_ID_GF2128:
        u128_select(transpose,a,NULL,NULL,1);
        ok=u128_weak_popov(transpose)>=0 && u128_kernel(out,a); break;
    case FIELD_ID_GF2N16:
        gf2n16_select(transpose,a,NULL,NULL,1);
        ok=gf2n16_weak_popov(transpose)>=0 && gf2n16_kernel(out,a); break;
    case FIELD_ID_GF2N64:
        gf2n64_select(transpose,a,NULL,NULL,1);
        ok=gf2n64_weak_popov(transpose)>=0 && gf2n64_kernel(out,a); break;
    case FIELD_ID_GF2N128:
        gf2n128_select(transpose,a,NULL,NULL,1);
        ok=gf2n128_weak_popov(transpose)>=0 && gf2n128_kernel(out,a); break;
    default:
        generic_select(transpose,a,NULL,NULL,1);
        ok=generic_weak_popov(transpose)>=0 && generic_kernel(out,a); break;
    }
    unified_mat_poly_clear(transpose);return ok;
}

int unified_mat_poly_det_hnf(unified_poly_t det, const unified_mat_poly_t a)
{
    if(a->r!=a->c || det->ctx!=a->ctx) return 0;
    unified_mat_poly_t result; unified_mat_poly_init(result,1,1,a->ctx);
    int ok;
    switch(a->ctx->field_id) {
    case FIELD_ID_GF2N8: ok=b8_det(result,a); break;
    case FIELD_ID_GF24: case FIELD_ID_GF28: ok=byte_det(result,a); break;
    case FIELD_ID_GF216: ok=u16_det(result,a); break;
    case FIELD_ID_GF232: ok=u32_det(result,a); break;
    case FIELD_ID_GF264: ok=u64_det(result,a); break;
    case FIELD_ID_GF2128: ok=u128_det(result,a); break;
    case FIELD_ID_GF2N16: ok=gf2n16_det(result,a); break;
    case FIELD_ID_GF2N64: ok=gf2n64_det(result,a); break;
    case FIELD_ID_GF2N128: ok=gf2n128_det(result,a); break;
    default: ok=generic_det(result,a); break;
    }
    if(ok) {
        unified_poly_fit_length(det,result->length);
        for(slong k=0;k<result->length;k++) unified_mat_poly_get(det->coeffs+k,result,k);
        det->length=result->length;
    }
    unified_mat_poly_clear(result);
    return ok;
}

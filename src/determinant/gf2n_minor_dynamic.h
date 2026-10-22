/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Arbitrary-modulus minors use the same packed exponent layout as GF(2^128).
 * No basis conversion and no expanded FLINT coefficients inside the DP. */
#include "../field/gf2n_mpoly_sparse.h"

static void minor_dynamic_fit(gf2128_mpoly_t p, slong len,
    flint_bitcnt_t bits, const gf2128_mpoly_ctx_t ctx)
{
    slong N=mpoly_words_per_exp(bits,ctx->minfo);
    if(p->coeffs_alloc<len) {
        p->coeffs_alloc=FLINT_MAX(len,2*p->coeffs_alloc);
        p->coeffs=flint_realloc(p->coeffs,p->coeffs_alloc*sizeof(*p->coeffs));
    }
    if(p->exps_alloc<N*len) {
        p->exps_alloc=FLINT_MAX(N*len,2*p->exps_alloc);
        p->exps=flint_realloc(p->exps,p->exps_alloc*sizeof(ulong));
    }
    p->bits=bits;
}
static void minor_dynamic_import(gf2128_mpoly_t out, const fq_nmod_mpoly_t a,
    const unified_mpoly_ctx_struct *ctx, const gf2128_mpoly_ctx_t nc)
{
    slong n=fq_nmod_ctx_degree(ctx->field_ctx->ctx.fq_ctx);
    minor_dynamic_fit(out,a->length,a->bits,nc);
    slong N=mpoly_words_per_exp(a->bits,nc->minfo);
    if(a->length) memcpy(out->exps,a->exps,N*a->length*sizeof(ulong));
    for(slong i=0;i<a->length;i++) {
        gf2128_t v={0,0};
        for(slong j=0;j<n;j++) if(a->coeffs[i*n+j]) {
            if(j<64) v.low|=UINT64_C(1)<<j;
            else v.high|=UINT64_C(1)<<(j-64);
        }
        out->coeffs[i]=v;
    }
    out->length=a->length;
}
static void minor_dynamic_export(fq_nmod_mpoly_t out, const gf2128_mpoly_t a,
    const unified_mpoly_ctx_struct *ctx, const gf2128_mpoly_ctx_t nc)
{
    slong n=fq_nmod_ctx_degree(ctx->field_ctx->ctx.fq_ctx);
    fq_nmod_mpoly_fit_length_reset_bits(out,a->length,a->bits,GET_FQ_CTX(ctx));
    slong N=mpoly_words_per_exp(a->bits,nc->minfo);
    if(a->length) memcpy(out->exps,a->exps,N*a->length*sizeof(ulong));
    for(slong i=0;i<a->length;i++) for(slong j=0;j<n;j++)
        out->coeffs[i*n+j]=j<64 ? (a->coeffs[i].low>>j)&1 : (a->coeffs[i].high>>(j-64))&1;
    _fq_nmod_mpoly_set_length(out,a->length,GET_FQ_CTX(ctx));
}
static int minor_dynamic_mul(gf2128_mpoly_t out, const gf2128_mpoly_t a,
    const gf2128_mpoly_t b, const gf2128_mpoly_ctx_t ctx, const gf2n_native_ctx *field)
{
    /* DP outputs are separate from both inputs, so retain reusable capacity. */
    FLINT_ASSERT(out!=a && out!=b);
    out->length=0;
    if(!a->length || !b->length) return 1;
    const gf2128_mpoly_struct *x=a,*y=b;
    if(x->length>y->length) {x=b; y=a;}
    ulong unused,main;
    flint_bitcnt_t bits=gf2n_product_bits(x->exps,x->length,x->bits,
        y->exps,y->length,y->bits,ctx->minfo,&unused,&main);
    slong N=mpoly_words_per_exp(bits,ctx->minfo);
    ulong *xe=flint_malloc(N*x->length*sizeof(ulong));
    ulong *ye=flint_malloc(N*y->length*sizeof(ulong));
    int ok=mpoly_repack_monomials(xe,bits,x->exps,x->bits,x->length,ctx->minfo) &&
           mpoly_repack_monomials(ye,bits,y->exps,y->bits,y->length,ctx->minfo);
    if(!ok) {flint_free(xe); flint_free(ye); return 0;}
    ulong *keys=flint_malloc(N*x->length*sizeof(ulong));
    ulong *mask=flint_malloc(N*sizeof(ulong)),*exp=flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(mask,N,bits,ctx->minfo);
    gf2n_product_heap_entry *heap=flint_malloc(x->length*sizeof(*heap));
    slong len=x->length;
    for(slong i=0;i<len;i++) {
        heap[i]=(gf2n_product_heap_entry){i,0,keys+N*i};
        gf2n_monomial_add(heap[i].exp,xe+N*i,ye,N,bits);
    }
    for(slong i=len/2;i-- >0;) gf2n_heap_sift_down(heap,len,i,N,mask);
    while(len) {
        memcpy(exp,heap[0].exp,N*sizeof(ulong));
        gf2128_t c={0,0};
        do {
            slong i=heap[0].row,j=heap[0].col;
            gf2128_t p=gf2n_native_mul(x->coeffs[i],y->coeffs[j],field);
            c.low^=p.low; c.high^=p.high;
            if(++j<y->length) {
                heap[0].col=j;
                gf2n_monomial_add(heap[0].exp,xe+N*i,ye+N*j,N,bits);
            } else heap[0]=heap[--len];
            if(len) gf2n_heap_sift_down(heap,len,0,N,mask);
        } while(len && mpoly_monomial_equal(exp,heap[0].exp,N));
        if(c.low|c.high) {
            minor_dynamic_fit(out,out->length+1,bits,ctx);
            out->coeffs[out->length]=c;
            memcpy(out->exps+N*out->length++,exp,N*sizeof(ulong));
        }
    }
    flint_free(heap); flint_free(mask); flint_free(exp); flint_free(keys);
    flint_free(xe); flint_free(ye); return 1;
}
#define DP_NAME(x) minor_dynamic_##x
#define DP_POLY gf2128_mpoly_struct
#define DP_POLY_T gf2128_mpoly_t
#define DP_CTX_T gf2128_mpoly_ctx_t
#define DP_E gf2128_t
#define DP_SET gf2128_mpoly_set
#define DP_INIT gf2128_mpoly_init
#define DP_CLEAR gf2128_mpoly_clear
#define DP_CTX_INIT gf2128_mpoly_ctx_init
#define DP_CTX_CLEAR gf2128_mpoly_ctx_clear
#define DP_MUL(o,a,b,nc) minor_dynamic_mul(o,a,b,nc,(const gf2n_native_ctx *)ctx->field_ctx->ctx.fq_ctx)
#define DP_IMPORT minor_dynamic_import
#define DP_EXPORT minor_dynamic_export
#define DP_ADD(a,b) gf2128_add(&(a),&(b))
#define DP_ZERO(a) gf2128_is_zero(&(a))
#include "gf2n_minor_dp_template.h"

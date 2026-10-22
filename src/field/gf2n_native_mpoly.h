/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Reuse the exponent planner and heap from the fixed-width backends. */
#include "gf2n_mpoly_sparse.h"

static int gf2n_native_mpoly_mul(fq_nmod_mpoly_t res,
    const fq_nmod_mpoly_t a, const fq_nmod_mpoly_t b,
    const fq_nmod_mpoly_ctx_t ctx, const gf2n_native_ctx *field)
{
    if(!a->length || !b->length) { fq_nmod_mpoly_zero(res,ctx); return 1; }
    const fq_nmod_mpoly_struct *x=a,*y=b;
    if(x->length>y->length) { x=b; y=a; }
    ulong unused,main;
    flint_bitcnt_t bits=gf2n_product_bits(x->exps,x->length,x->bits,
        y->exps,y->length,y->bits,ctx->minfo,&unused,&main);
    slong N=mpoly_words_per_exp(bits,ctx->minfo);
    ulong *xe=flint_malloc(N*x->length*sizeof(ulong));
    ulong *ye=flint_malloc(N*y->length*sizeof(ulong));
    int ok=mpoly_repack_monomials(xe,bits,x->exps,x->bits,x->length,ctx->minfo) &&
           mpoly_repack_monomials(ye,bits,y->exps,y->bits,y->length,ctx->minfo);
    if(!ok) { flint_free(xe); flint_free(ye); return 0; }
    gf2128_t *xc=flint_malloc((x->length+y->length)*sizeof(*xc)),*yc=xc+x->length;
    unsigned n=field->n;
    for(slong i=0;i<x->length+y->length;i++) {
        const ulong *src=i<x->length ? x->coeffs+i*n : y->coeffs+(i-x->length)*n;
        xc[i]=(gf2128_t){0,0};
        for(unsigned j=0;j<n;j++) if(src[j]) {
            if(j<64) xc[i].low|=UINT64_C(1)<<j;
            else xc[i].high|=UINT64_C(1)<<(j-64);
        }
    }
    ulong *keys=flint_malloc(N*x->length*sizeof(ulong));
    ulong *mask=flint_malloc(N*sizeof(ulong)),*exp=flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(mask,N,bits,ctx->minfo);
    gf2n_product_heap_entry *heap=flint_malloc(x->length*sizeof(*heap));
    slong len=x->length;
    for(slong i=0;i<len;i++) {
        heap[i]=(gf2n_product_heap_entry){i,0,keys+N*i};
        gf2n_monomial_add(heap[i].exp,xe+N*i,ye,N,bits);
    }
    for(slong i=len/2;i-- > 0;) gf2n_heap_sift_down(heap,len,i,N,mask);
    fq_nmod_mpoly_t out; fq_nmod_mpoly_init(out,ctx);
    fq_nmod_mpoly_fit_length_reset_bits(out,16,bits,ctx);
    while(len) {
        memcpy(exp,heap[0].exp,N*sizeof(ulong));
        gf2128_t coeff={0,0};
        do {
            slong i=heap[0].row,j=heap[0].col;
            gf2128_t p=gf2n_native_mul(xc[i],yc[j],field);
            coeff.low^=p.low; coeff.high^=p.high;
            if(++j<y->length) {
                heap[0].col=j;
                gf2n_monomial_add(heap[0].exp,xe+N*i,ye+N*j,N,bits);
            } else heap[0]=heap[--len];
            if(len) gf2n_heap_sift_down(heap,len,0,N,mask);
        } while(len && mpoly_monomial_equal(exp,heap[0].exp,N));
        if(coeff.low|coeff.high) {
            slong k=out->length;
            fq_nmod_mpoly_fit_length(out,k+1,ctx);
            for(unsigned j=0;j<n;j++)
                out->coeffs[k*n+j]=j<64 ? (coeff.low>>j)&1 : (coeff.high>>(j-64))&1;
            memcpy(out->exps+N*k,exp,N*sizeof(ulong));
            _fq_nmod_mpoly_set_length(out,k+1,ctx);
        }
    }
    fq_nmod_mpoly_swap(res,out,ctx); fq_nmod_mpoly_clear(out,ctx);
    flint_free(heap); flint_free(exp); flint_free(mask); flint_free(keys);
    flint_free(xc); flint_free(xe); flint_free(ye);
    return 1;
}

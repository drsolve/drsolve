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
    ulong *mask=flint_malloc(N*sizeof(ulong)),*exp=flint_malloc(N*sizeof(ulong));
    mpoly_get_cmpmask(mask,N,bits,ctx->minfo);
    gf2n_product_iterator iterator;
    gf2n_product_iterator_init(&iterator,xe,x->length,ye,y->length,N,bits,mask);
    fq_nmod_mpoly_t out; fq_nmod_mpoly_init(out,ctx);
    fq_nmod_mpoly_fit_length_reset_bits(out,16,bits,ctx);
    while(gf2n_product_iterator_active(&iterator)) {
        memcpy(exp,gf2n_product_iterator_key(&iterator),N*sizeof(ulong));
        gf2128_t coeff={0,0};
        do {
            slong i,j;gf2n_product_iterator_pop(&iterator,&i,&j);
            gf2128_t p=gf2n_native_mul_fast(xc[i],yc[j],field);
            coeff.low^=p.low; coeff.high^=p.high;
        } while(gf2n_product_iterator_active(&iterator) &&
                mpoly_monomial_equal(exp,gf2n_product_iterator_key(&iterator),N));
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
    gf2n_product_iterator_clear(&iterator); flint_free(exp); flint_free(mask);
    flint_free(xc); flint_free(xe); flint_free(ye);
    return 1;
}

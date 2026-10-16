/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include "quotient_solver.h"
#include "fq_unified_interface.h"

static void check(const fq_nmod_ctx_t field, const char **names, slong n,
                  const char **input, slong m, const char *raw,
                  const char *expected, slong limit)
{
    fq_nmod_mpoly_ctx_t ctx; fq_nmod_mpoly_ctx_init(ctx,n,ORD_DEGLEX,field);
    fq_nmod_mpoly_struct fs[m];
    for (slong i=0;i<m;i++) {
        fq_nmod_mpoly_init(fs+i,ctx);
        assert(!fq_nmod_mpoly_set_str_pretty(fs+i,input[i],names,ctx));
    }
    fq_nmod_mpoly_t d,s,r,w,product;
    fq_nmod_mpoly_init(d,ctx); fq_nmod_mpoly_init(s,ctx); fq_nmod_mpoly_init(r,ctx);
    fq_nmod_mpoly_init(w,ctx); fq_nmod_mpoly_init(product,ctx);
    assert(!fq_nmod_mpoly_set_str_pretty(d,raw,names,ctx));
    assert(!fq_nmod_mpoly_set_str_pretty(w,expected,names,ctx));
    slong used;
    assert(quotient_sift_resultant_extension(s,r,d,fs,m,ctx,limit,&used));
    if (!fq_nmod_mpoly_equal(s,w,ctx)) {
        char *got=fq_nmod_mpoly_get_str_pretty(s,names,ctx);
        fprintf(stderr,"Expected %s, got %s\n",expected,got); flint_free(got); abort();
    }
    assert(used<=limit);
    fq_nmod_mpoly_mul(product,s,r,ctx);
    assert(fq_nmod_mpoly_equal(product,d,ctx));
    for (slong i=0;i<m;i++) fq_nmod_mpoly_clear(fs+i,ctx);
    fq_nmod_mpoly_clear(d,ctx); fq_nmod_mpoly_clear(s,ctx); fq_nmod_mpoly_clear(r,ctx);
    fq_nmod_mpoly_clear(w,ctx); fq_nmod_mpoly_clear(product,ctx); fq_nmod_mpoly_ctx_clear(ctx);
}

static void cases(const fq_nmod_ctx_t field)
{
    const char *names[]={"x","a","b"};
    const char *nonradical[]={"x-t","(a-t)^2"};
    const char *raw="(a-t)^4*(a-t-1)";
    check(field,names,2,nonradical,2,raw,"(a-t)^2",4);
    check(field,names,2,nonradical,2,raw,raw,1);
    /* General parameter support and a positive-dimensional ideal. */
    const char *positive[]={"(a-t)*x","(a-t)*(x-1)"};
    check(field,names,3,positive,2,"t*(a-t)^2*(b-t)","t*(a-t)",3);
    /* A bounded closure relation can also certify a higher-degree multiple. */
    const char *high[]={"x","(a-t)^8"};
    check(field,names,2,high,2,"(a-t)^12*(a-t-1)","(a-t)^8",8);
    const char *empty[]={"x","x-t"};
    check(field,names,2,empty,2,"a^2","1",3);
    check(field,names,2,empty,2,"0","0",3);
    check(field,names,2,nonradical,2,"t+1","t+1",3);
    /* A zero generator must not become a constant contradiction. */
    const char *zero[]={"0","(a-t)^2"};
    check(field,names,2,zero,2,"(a-t)^3","(a-t)^2",3);
}

static void coupled(const fq_nmod_ctx_t field)
{
    const char *names[]={"x","y","b","a"};
    const char *fs[]={"x^2+y+b+(a-t)","x*y*b+x*(a-t)-1","x*(a-t)*b+(a-t)"};
    const char *raw="(a-t)^3*b^3*(b^3+(a-t)*b^2-b^2-(a-t)*b+1)";
    check(field,names,4,fs,3,raw,"(a-t)*(b^3+(a-t)*b^2-b^2-(a-t)*b+1)",6);
}

int main(void)
{
    const slong ds[]={4,8,16,32,64,128};
    const uint64_t mods[]={0x13,0x11d,0x1002d,0x100008299,UINT64_C(0x247f43cb7),0x87};
    for (int i=0;i<6;i++) {
        nmod_poly_t mod; nmod_poly_init(mod,2); nmod_poly_set_coeff_ui(mod,ds[i],1);
        for (int j=0;j<64;j++) if ((mods[i]>>j)&1) nmod_poly_set_coeff_ui(mod,j,1);
        fq_nmod_ctx_t field; fq_nmod_ctx_init_modulus(field,mod,"t");
        field_ctx_t backend; field_ctx_init(&backend,field);
        assert(backend.field_id==FIELD_ID_GF24+i); field_ctx_clear(&backend);
        cases(field); if (i==1) coupled(field);
        fq_nmod_ctx_clear(field); nmod_poly_clear(mod);
    }
    const ulong ps[]={3,5,1031};
    for (int i=0;i<3;i++) {
        fq_nmod_ctx_t field; fq_nmod_ctx_init_ui(field,ps[i],2,"t");
        cases(field);
        if (!i) {
            coupled(field);
            /* t is nonsquare in this F_9 context: sifting must not use the
             * absence of F_9-points to cancel the genuine a-t factor. */
            const char *names[]={"x","a"},*fs[]={"x^2-t","a-t"};
            check(field,names,2,fs,2,"(a-t)*(a-t-1)","a-t",4);
        }
        fq_nmod_ctx_clear(field);
    }
    flint_cleanup(); puts("Quotient extension sift: membership, multiplicity, extension points and limits passed");
    return 0;
}

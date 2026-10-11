/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include "quotient_solver.h"

static void check(const char **names, slong n, const char **input, slong m,
                  const char *raw, const char *expected, slong limit)
{
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx,n,ORD_DEGLEX,65537);
    nmod_mpoly_struct *fs=flint_malloc(m*sizeof(*fs));
    for (slong i=0;i<m;i++) {
        nmod_mpoly_init(fs+i,ctx);
        assert(!nmod_mpoly_set_str_pretty(fs+i,input[i],names,ctx));
    }
    nmod_mpoly_t d,s,r,w,product;
    nmod_mpoly_init(d,ctx); nmod_mpoly_init(s,ctx); nmod_mpoly_init(r,ctx);
    nmod_mpoly_init(w,ctx); nmod_mpoly_init(product,ctx);
    assert(!nmod_mpoly_set_str_pretty(d,raw,names,ctx));
    assert(!nmod_mpoly_set_str_pretty(w,expected,names,ctx));
    slong used;
    assert(quotient_sift_resultant(s,r,d,fs,m,ctx,limit,&used));
    if (!nmod_mpoly_equal(s,w,ctx)) {
        char *got=nmod_mpoly_get_str_pretty(s,names,ctx);
        fprintf(stderr,"Expected %s, got %s\n",expected,got); flint_free(got);
        abort();
    }
    nmod_mpoly_mul(product,s,r,ctx);
    assert(nmod_mpoly_equal(product,d,ctx));
    for (slong i=0;i<m;i++) nmod_mpoly_clear(fs+i,ctx);
    flint_free(fs);
    nmod_mpoly_clear(d,ctx); nmod_mpoly_clear(s,ctx); nmod_mpoly_clear(r,ctx);
    nmod_mpoly_clear(w,ctx); nmod_mpoly_clear(product,ctx); nmod_mpoly_ctx_clear(ctx);
}

int main(void)
{
    const char *names[]={"x","y","b","a"};
    const char *fs[]={"x^2+y+b+a","x*y*b+x*a-1","x*a*b+a"};
    const char *d="a^3*b^3*(b^3+a*b^2-b^2-a*b+1)";
    check(names,4,fs,3,d,"a*(b^3+a*b^2-b^2-a*b+1)",8);
    check(names,4,fs,3,d,d,1); /* Budget exhaustion preserves input. */
    const char *small_names[]={"x","t"};
    const char *nonradical[]={"x","t^2"};
    check(small_names,2,nonradical,2,"t^4*(t-1)","t^2",6);
    /* 3 is a nonsquare modulo 65537: lack of rational points must not
     * erase the genuine factor t. */
    const char *extension_points[]={"x^2-3","t"};
    check(small_names,2,extension_points,2,"t*(t-1)","t",4);
    const char *inconsistent[]={"x","x-1"};
    check(small_names,2,inconsistent,2,"t^2","1",3);
    check(small_names,2,inconsistent,2,"0","0",3);
    check(small_names,2,nonradical,2,"7","7",3);
    puts("Quotient sift: exact cancellation, multiplicity, extension points and budget tests passed");
    return 0;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "../../resultant/bezout_prefix.h"
#include <flint/ulong_extras.h>
#include <stdlib.h>
#include <stdio.h>
int main(void)
{
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,123,456);
    ulong primes[]={2,17,65537,998244353,1000000007,
#if FLINT_BITS == 64
    UWORD(28407454060060787),UWORD(18446744069414584321),
#endif
};
    slong degrees[]={1,2,3,8,16,31,64,128,257,512};long checks=0,blocks=0,fallbacks=0;
    for(size_t p=0;p<sizeof(primes)/sizeof(*primes);p++) {
        nmod_poly_t f,g,expected,out,tmp,old,alias,rev,moments;
        nmod_poly_struct *polys[]={f,g,expected,out,tmp,old,alias,rev,moments};
        for(size_t j=0;j<9;j++)nmod_poly_init(polys[j],primes[p]);
        for(size_t nd=0;nd<10;nd++)for(int trial=0;trial<8;trial++) {
            slong n=degrees[nd];
            nmod_poly_randtest(f,rng,n+1);nmod_poly_set_coeff_ui(f,n,1+n_randint(rng,primes[p]-1));
            nmod_poly_randtest(g,rng,1+n_randint(rng,n));
            if(trial==0)nmod_poly_zero(g);
            int invertible=nmod_poly_invmod(tmp,g,f);
            slong limits[]={1,2,7,2*n-1};
            for(size_t l=0;l<4;l++)for(slong cutoff=0;cutoff<=64;cutoff+=16) {
                slong limit=limits[l];dr_bezout_prefix_stats stats;
                nmod_poly_zero(out);nmod_poly_set_coeff_ui(out,0,7);nmod_poly_set(old,out);
                int ok=dr_bezout_inverse_prefix(out,f,g,limit,cutoff,&stats);
                blocks+=stats.hgcd_blocks;fallbacks+=stats.full_inverse_fallbacks;
                if(invertible) {
                    nmod_poly_reverse(expected,tmp,n);nmod_poly_truncate(expected,limit);
                    if(ok!=1 || !nmod_poly_equal(expected,out))goto fail;
                    nmod_poly_set(alias,f);
                    if(dr_bezout_inverse_prefix(alias,alias,g,limit,cutoff,NULL)!=1 || !nmod_poly_equal(alias,out))goto fail;
                    nmod_poly_set(alias,g);
                    if(dr_bezout_inverse_prefix(alias,f,alias,limit,cutoff,NULL)!=1 || !nmod_poly_equal(alias,out))goto fail;
                    nmod_poly_reverse(rev,f,n+1);nmod_poly_div_series(expected,out,rev,limit);
                    if(dr_bezout_left_moments(moments,f,g,limit,cutoff,NULL)!=1 || !nmod_poly_equal(moments,expected))goto fail;
                    if(dr_bezout_left_moments_profile(moments,f,g,limit,cutoff,&stats,1)!=1 || !nmod_poly_equal(moments,expected))goto fail;
                    nmod_poly_set(alias,f);
                    if(dr_bezout_left_moments(alias,alias,g,limit,cutoff,NULL)!=1 || !nmod_poly_equal(alias,expected))goto fail;
                    nmod_poly_set(alias,g);
                    if(dr_bezout_left_moments(alias,f,alias,limit,cutoff,NULL)!=1 || !nmod_poly_equal(alias,expected))goto fail;
                } else if(ok!=0 || !nmod_poly_equal(out,old))goto fail;
                checks++;continue;
fail:
                fprintf(stderr,"FAIL p=%lu n=%ld trial=%d limit=%ld cutoff=%ld inverse=%d status=%d\n",primes[p],n,trial,limit,cutoff,invertible,ok);return 1;
            }
        }
        nmod_poly_set(old,out);
        if(dr_bezout_inverse_prefix(out,f,g,0,64,NULL)!=-1 || !nmod_poly_equal(old,out) ||
           dr_bezout_inverse_prefix(out,f,f,1,64,NULL)!=-1 || !nmod_poly_equal(old,out) ||
           dr_bezout_inverse_prefix(out,f,g,1,-1,NULL)!=-1 || !nmod_poly_equal(old,out))return 1;
        for(size_t j=0;j<9;j++)nmod_poly_clear(polys[j]);
    }
    printf("C_PREFIX_CHECKS=%ld HGCD_BLOCKS=%ld FALLBACKS=%ld PASS\n",checks,blocks,fallbacks);
    flint_rand_clear(rng);flint_cleanup();return 0;
}

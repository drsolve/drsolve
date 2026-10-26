/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include <stdio.h>
#include "bezout_block.h"
#include "unified_mpoly_resultant.h"
#include "dixon_flint.h"

int main(void)
{
    g_dixon_verbose_level=0;g_dixon_bezout_block=0;
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,711,98);
    ulong primes[]={17,65537,998244353,1000000007,
#if FLINT_BITS == 64
        UWORD(18446744069414584321),
#endif
    };
    slong passed=0,fallbacks=0,extended=0,unbalanced=0;
    for(size_t pi=0;pi<sizeof(primes)/sizeof(*primes);pi++) {
        fmpz_t prime;fmpz_init_set_ui(prime,primes[pi]);
        field_ctx_t field;field_ctx_init_prime(&field,prime);
        unified_mpoly_ctx_t ctx=unified_mpoly_ctx_init(2,ORD_LEX,&field);
        unified_mpoly_t a=unified_mpoly_init(ctx),b=unified_mpoly_init(ctx),ref=unified_mpoly_init(ctx);
        unified_mpoly_t got=unified_mpoly_init(ctx),saved=unified_mpoly_init(ctx);
        nmod_poly_t expected,actual,sentinel;
        nmod_poly_init(expected,primes[pi]);nmod_poly_init(actual,primes[pi]);nmod_poly_init(sentinel,primes[pi]);
        nmod_poly_set_coeff_ui(sentinel,2,3);
        for(int shape=0;shape<4;shape++)
        for(slong var=0;var<2;var++)for(slong trial=0;trial<9;trial++) {
            slong d=pi==0?2:2+3*(trial%4);
            if(shape && trial==0)d=1;
            unified_mpoly_zero(a);unified_mpoly_zero(b);
            int large=trial%2;
            for(int k=0;k<2;k++) {
              slong degree=(shape && k==large)?2*d+1:d;
              slong xdegree=(shape>=2 && k==large)?(shape==2?d-1:d+1):degree;
              for(ulong x=0;x<=(ulong)xdegree;x++)for(ulong y=0;y<=(ulong)degree-x;y++) {
                ulong e[2];e[var]=x;e[1-var]=y;
                ulong c=n_randint(rng,primes[pi]);
                if(x==(ulong)xdegree || y==(ulong)degree)c=1+n_randint(rng,primes[pi]-1);
                nmod_mpoly_set_coeff_ui_ui(GET_NMOD_POLY(k?b:a),c,e,GET_NMOD_CTX(ctx));
              }
            }
            assert(nmod_mpoly_resultant(GET_NMOD_POLY(ref),GET_NMOD_POLY(a),GET_NMOD_POLY(b),var,GET_NMOD_CTX(ctx)));
            assert(nmod_mpoly_get_nmod_poly(expected,GET_NMOD_POLY(ref),1-var,GET_NMOD_CTX(ctx)));
            nmod_poly_set(actual,sentinel);dr_bezout_block_stats stats;
            int ok=dr_nmod_mpoly_resultant_bezout_block_try(actual,GET_NMOD_POLY(a),GET_NMOD_POLY(b),var,GET_NMOD_CTX(ctx),0,rng,&stats);
            if(ok){assert(nmod_poly_equal(actual,expected));passed++;unbalanced+=shape!=0;extended+=stats.verification_extension_degree>1;}
            else {assert(nmod_poly_equal(actual,sentinel));fallbacks++;}
            if(pi>0 && trial==3 && var==0) {
                slong blocks[]={1,d,d-1,d+1,-1};
                for(int j=0;j<5;j++) {
                    nmod_poly_set(actual,sentinel);
                    int accepted=dr_nmod_mpoly_resultant_bezout_block_try(actual,GET_NMOD_POLY(a),GET_NMOD_POLY(b),var,GET_NMOD_CTX(ctx),blocks[j],rng,NULL);
                    if(accepted) {assert(j<3);assert(nmod_poly_equal(actual,expected));}
                    else assert(nmod_poly_equal(actual,sentinel));
                }
            }
            /* Generic dispatch must also work when the opt-in declines. */
            g_dixon_bezout_block=1;
            assert(unified_mpoly_resultant(got,a,b,var,ctx));assert(unified_mpoly_equal(got,ref));
            g_dixon_bezout_block=0;
            unified_mpoly_set(saved,a);
            ok=unified_mpoly_resultant_bezout_block_try(a,a,b,var,ctx);
            if(ok)assert(unified_mpoly_equal(a,ref));else assert(unified_mpoly_equal(a,saved));
            unified_mpoly_set(a,saved);unified_mpoly_set(saved,b);
            ok=unified_mpoly_resultant_bezout_block_try(b,a,b,var,ctx);
            if(ok)assert(unified_mpoly_equal(b,ref));else assert(unified_mpoly_equal(b,saved));
            unified_mpoly_set(b,saved);
            nmod_poly_set(actual,sentinel);
            assert(!dr_nmod_mpoly_resultant_bezout_block_try(actual,GET_NMOD_POLY(a),GET_NMOD_POLY(a),var,GET_NMOD_CTX(ctx),0,rng,NULL));
            assert(nmod_poly_equal(actual,sentinel));
        }
        nmod_poly_set(actual,sentinel);
        assert(!dr_nmod_mpoly_resultant_bezout_block_try(actual,GET_NMOD_POLY(a),GET_NMOD_POLY(b),2,GET_NMOD_CTX(ctx),0,rng,NULL));
        assert(nmod_poly_equal(actual,sentinel));
        /* Analytic witness with resultant degree above 65536. No large-degree
         * reference resultant is needed: Res(x^257+x-y^257,x^256)=y^65792. */
        if(pi==1) {
            unified_mpoly_zero(a);unified_mpoly_zero(b);
            ulong e[2]={257,0};
            nmod_mpoly_set_coeff_ui_ui(GET_NMOD_POLY(a),1,e,GET_NMOD_CTX(ctx));
            e[0]=1;nmod_mpoly_set_coeff_ui_ui(GET_NMOD_POLY(a),1,e,GET_NMOD_CTX(ctx));
            e[0]=0;e[1]=257;
            nmod_mpoly_set_coeff_ui_ui(GET_NMOD_POLY(a),primes[pi]-1,e,GET_NMOD_CTX(ctx));
            e[0]=256;e[1]=0;
            nmod_mpoly_set_coeff_ui_ui(GET_NMOD_POLY(b),1,e,GET_NMOD_CTX(ctx));
            nmod_poly_zero(expected);nmod_poly_set_coeff_ui(expected,65792,1);
            assert(dr_nmod_mpoly_resultant_bezout_block_try(actual,GET_NMOD_POLY(a),GET_NMOD_POLY(b),0,GET_NMOD_CTX(ctx),0,rng,NULL));
            assert(nmod_poly_equal(actual,expected));
        }
        nmod_poly_clear(expected);nmod_poly_clear(actual);nmod_poly_clear(sentinel);
        unified_mpoly_clear(a);unified_mpoly_clear(b);unified_mpoly_clear(ref);unified_mpoly_clear(got);unified_mpoly_clear(saved);
        unified_mpoly_ctx_clear(ctx);field_ctx_clear(&field);fmpz_clear(prime);
    }
    assert(passed>=40 && unbalanced>=150 && extended>0);
    flint_rand_clear(rng);flint_cleanup();
    printf("Bezout block: %ld successful exact coefficient comparisons (%ld unequal degrees), %ld safe fallbacks, %ld extension checks; aliases and dispatch PASS\n",passed,unbalanced,fallbacks,extended);
    return 0;
}

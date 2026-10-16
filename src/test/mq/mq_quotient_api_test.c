/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_quotient.h"
#include <assert.h>
#include <flint/nmod_mpoly.h>

int main(void)
{
    fmpz_t prime;fmpz_init_set_ui(prime,65537);
    fq_nmod_ctx_t ctx;fq_nmod_ctx_init(ctx,prime,1,"a");
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,12345,54321);
    fq_nmod_t coefficient;fq_nmod_init(coefficient,ctx);
    flint_set_num_threads(2);
    const slong counts[]={1,2,3,4,5,3};
    for(size_t trial=0;trial<sizeof(counts)/sizeof(counts[0]);trial++) {
        slong nv=counts[trial],n=nv+1,*exp=flint_calloc(n,sizeof(slong));
        unified_mpoly_struct *p=flint_calloc(n,sizeof(*p)),*zero=flint_calloc(n,sizeof(*zero)),a={0},b={0};
        for(slong f=0;f<n;f++) {
            dr_mpoly_init(p+f,nv,1,ctx);dr_mpoly_init(zero+f,nv,1,ctx);
            for(slong i=0;i<=n;i++)for(slong j=i;j<=n;j++) {
                if(i<n)exp[i]++;if(j<n)exp[j]++;
                fq_nmod_set_ui(coefficient,n_randint(rng,65537),ctx);
                dr_mpoly_add_term_fast(p+f,exp,exp+nv,coefficient);
                if(i<n)exp[i]--;if(j<n)exp[j]--;
            }
            dr_mpoly_normalize(p+f);
        }
        mq_quotient_stats stats;
        assert(dr_mq_quotient_resultant(&a,p,nv,1,2,&stats));
        assert(flint_get_num_threads()==2);
        assert(dr_mq_quotient_resultant(&b,p,nv,1,1,&stats));
        assert(nmod_mpoly_is_canonical(GET_NMOD_POLY(&a),GET_NMOD_CTX(a.ctx_ptr)));
        assert(nmod_mpoly_equal(GET_NMOD_POLY(&a),GET_NMOD_POLY(&b),GET_NMOD_CTX(a.ctx_ptr)));
        /* Failure must leave a pre-existing output intact and release frames. */
        assert(!dr_mq_quotient_resultant(&a,zero,nv,1,2,&stats));
        assert(stats.reason && flint_get_num_threads()==2);
        assert(nmod_mpoly_equal(GET_NMOD_POLY(&a),GET_NMOD_POLY(&b),GET_NMOD_CTX(a.ctx_ptr)));
        dr_mpoly_clear(&b);
        assert(dr_mq_quotient_resultant(&b,p,nv,1,2,&stats));
        assert(nmod_mpoly_equal(GET_NMOD_POLY(&a),GET_NMOD_POLY(&b),GET_NMOD_CTX(a.ctx_ptr)));
        assert(!dr_mq_quotient_resultant(&a,p,16,1,2,&stats));
        assert(nmod_mpoly_equal(GET_NMOD_POLY(&a),GET_NMOD_POLY(&b),GET_NMOD_CTX(a.ctx_ptr)));
        dr_mpoly_clear(&a);dr_mpoly_clear(&b);
        for(slong f=0;f<n;f++){dr_mpoly_clear(p+f);dr_mpoly_clear(zero+f);}
        flint_free(p);flint_free(zero);flint_free(exp);
    }
    fq_nmod_clear(coefficient,ctx);flint_rand_clear(rng);fq_nmod_ctx_clear(ctx);
    fmpz_clear(prime);flint_cleanup_master();puts("MQ quotient API checks passed");return 0;
}

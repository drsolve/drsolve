/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Exercise ownership, repeated calls, parallel sampling, and rejection after
 * admission. CLI tests independently compare the polynomials with Dixon. */
#include "mq_quotient.h"
#include "dixon_flint.h"
#include <assert.h>
#include <string.h>

static void check_field(ulong characteristic,slong degree)
{
    fmpz_t prime;fmpz_init_set_ui(prime,characteristic);
    fq_nmod_ctx_t ctx;fq_nmod_ctx_init(ctx,prime,degree,"t");
    fq_nmod_t c;fq_nmod_init(c,ctx);
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,12345,54321);
    unified_mpoly_struct p[3]={{0}},a={0},b={0},singular[3],bad[3]={{0}};
    mq_quotient_stats stats;
    int admitted=0;
    for(int attempt=0;attempt<10 && !admitted;attempt++) {
        for(slong f=0;f<3;f++) {
            dr_mpoly_init(p+f,2,1,ctx);
            slong exp[3]={0};
            for(slong i=0;i<=3;i++)for(slong j=i;j<=3;j++) {
                if(i<3)exp[i]++;if(j<3)exp[j]++;
                fq_nmod_randtest(c,rng,ctx);
                dr_mpoly_add_term_fast(p+f,exp,exp+2,c);
                if(i<3)exp[i]--;if(j<3)exp[j]--;
            }
            dr_mpoly_normalize(p+f);
        }
        admitted=dr_mq_quotient_resultant(&a,p,2,1,1,&stats);
        if(!admitted)for(slong f=0;f<3;f++)dr_mpoly_clear(p+f);
    }
    assert(admitted && stats.reason==NULL && stats.core_size==3);
    assert(flint_get_num_threads()==2);
    assert(dr_degree_quotient_resultant(&b,p,2,1,4,&stats));
    assert(unified_mpoly_equal(&a,&b));
    for(slong f=0;f<3;f++)singular[f]=p[0]; /* borrowed */
    for(int retry=0;retry<2;retry++) {
        assert(!dr_extension_quotient_resultant(&b,singular,2,1,4,&stats));
        assert(strstr(stats.reason,"rank") && unified_mpoly_equal(&a,&b));
        assert(flint_get_num_threads()==2);
    }
    /* Multiplication by (parameter-1) preserves the leading basis at zero,
     * then destroys rank at a worker's interpolation point. */
    for(slong f=0;f<3;f++) {
        dr_mpoly_init(bad+f,2,1,ctx);
        for(slong j=0;j<dr_mpoly_length(p+f);j++) {
            DR_MPOLY_TERM(term,p+f,j);
            slong exp=term.par_exp[0];
            fq_nmod_neg(c,term.coeff,ctx);
            dr_mpoly_add_term_fast(bad+f,term.var_exp,&exp,c);
            exp++;dr_mpoly_add_term_fast(bad+f,term.var_exp,&exp,term.coeff);
        }
        dr_mpoly_normalize(bad+f);
    }
    assert(!dr_degree_quotient_resultant(&b,bad,2,1,4,&stats));
    assert(strstr(stats.reason,"rank") && unified_mpoly_equal(&a,&b));
    assert(flint_get_num_threads()==2);
    assert(!dr_extension_quotient_resultant(&b,p,2,2,4,&stats));
    assert(unified_mpoly_equal(&a,&b));
    /* Reuse after failures must release all worker and field contexts. */
    dr_mpoly_clear(&b);
    assert(dr_extension_quotient_resultant(&b,p,2,1,4,NULL));
    assert(unified_mpoly_equal(&a,&b));
    dr_mpoly_clear(&a);dr_mpoly_clear(&b);
    for(slong f=0;f<3;f++){dr_mpoly_clear(p+f);dr_mpoly_clear(bad+f);}
    fq_nmod_clear(c,ctx);flint_rand_clear(rng);fq_nmod_ctx_clear(ctx);fmpz_clear(prime);
}
/* Beyond both the old eight-exponent packing and the MQ packed-key range.
 * A triangular linear system has the independently known resultant y - t. */
static void check_dynamic_keys(ulong characteristic,slong degree)
{
    enum { n=17, nv=n-1 };
    fmpz_t prime;fmpz_init_set_ui(prime,characteristic);
    fq_nmod_ctx_t ctx;fq_nmod_ctx_init(ctx,prime,degree,"t");
    fq_nmod_t c;fq_nmod_init(c,ctx);
    unified_mpoly_struct *p=flint_calloc(n,sizeof(*p)),a={0},expected={0};
    slong e[nv];memset(e,0,sizeof(e));
    for(slong f=0;f<nv;f++) {
        dr_mpoly_init(p+f,nv,1,ctx);
        slong pe=0;e[f]=1;fq_nmod_one(c,ctx);
        dr_mpoly_add_term_fast(p+f,e,&pe,c);e[f]=0;
        pe=1;fq_nmod_neg(c,c,ctx);dr_mpoly_add_term_fast(p+f,e,&pe,c);
        dr_mpoly_normalize(p+f);
    }
    dr_mpoly_init(p+nv,nv,1,ctx);dr_mpoly_init(&expected,0,1,ctx);
    slong pe=1;fq_nmod_one(c,ctx);
    dr_mpoly_add_term_fast(p+nv,e,&pe,c);dr_mpoly_add_term_fast(&expected,NULL,&pe,c);
    pe=0;fq_nmod_gen(c,ctx);fq_nmod_neg(c,c,ctx);
    dr_mpoly_add_term_fast(p+nv,e,&pe,c);dr_mpoly_add_term_fast(&expected,NULL,&pe,c);
    dr_mpoly_normalize(p+nv);dr_mpoly_normalize(&expected);
    mq_quotient_stats stats;
    assert(dr_degree_quotient_resultant(&a,p,nv,1,4,&stats));
    assert(unified_mpoly_equal(&a,&expected));
    assert(!dr_degree_quotient_resultant(&a,p,FLINT_BITS,1,1,&stats));
    assert(strstr(stats.reason,"address capacity") && unified_mpoly_equal(&a,&expected));
    dr_mpoly_clear(&a);dr_mpoly_clear(&expected);
    for(slong f=0;f<n;f++)dr_mpoly_clear(p+f);flint_free(p);
    fq_nmod_clear(c,ctx);fq_nmod_ctx_clear(ctx);fmpz_clear(prime);
}
int main(void)
{
    g_dixon_verbose_level=0;flint_set_num_threads(2);
    const slong binary_degrees[]={4,8,16,32,64,128,5};
    for(size_t i=0;i<sizeof(binary_degrees)/sizeof(*binary_degrees);i++)
        check_field(2,binary_degrees[i]);
    check_field(3,3);check_field(101,3);
    check_dynamic_keys(2,8);check_dynamic_keys(101,3);check_dynamic_keys(65537,1);
    flint_cleanup_master();puts("Extension quotient API tests passed");return 0;
}

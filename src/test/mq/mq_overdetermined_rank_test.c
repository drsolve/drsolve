/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Measure truncated Macaulay ranks, not a complete overdetermined solver.
 * Eight affine variables, eight/nine quadratics sharing the origin. */
#include "../../determinant/mq_quotient_impl.h"

static slong relation_rank(const grade *g,const grade *lower,const ulong *q,slong equations)
{
    nmod_mat_t a;nmod_mat_init(a,equations*lower->s,g->s,prime);
    for(slong f=0;f<equations;f++)for(slong r=0;r<lower->s;r++)
        for(slong i=0;i<N;i++)for(slong j=i;j<N;j++) {
            uint64_t key=lower->keys[r]+(UINT64_C(1)<<(key_bits*i))+(UINT64_C(1)<<(key_bits*j));
            nmod_mat_entry(a,f*lower->s+r,index_key(g,key))=q[(f*N+i)*N+j];
        }
    slong rank=nmod_mat_rank(a);nmod_mat_clear(a);return rank;
}
int main(void)
{
    assert(quotient_dimensions(9));prime=65537;nmod_init(&mod,prime);
    flint_set_num_threads(1);
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,12345,54321);
    ulong *independent=flint_calloc(N*N*N,sizeof(ulong));
    ulong *redundant=flint_calloc(N*N*N,sizeof(ulong));
    for(slong f=0;f<N;f++)for(slong i=0;i<N;i++)for(slong j=i;j<N;j++)
        if(i!=N-1 || j!=N-1)independent[(f*N+i)*N+j]=n_randint(rng,prime);
    memcpy(redundant,independent,N*N*N*sizeof(ulong));
    for(slong i=0;i<N;i++)for(slong j=i;j<N;j++) {
        ulong c=0;for(slong f=0;f<N-1;f++)c=nmod_add(c,independent[(f*N+i)*N+j],mod);
        redundant[((N-1)*N+i)*N+j]=c;
    }
    grade *g=flint_calloc(L+1,sizeof(grade));
    for(slong d=0;d<=L;d++)grade_support(g+d,d);
    for(slong d=2;d<=L;d++) {
        slong a=relation_rank(g+d,g+d-2,independent,8);
        slong b=relation_rank(g+d,g+d-2,redundant,9);
        slong c=relation_rank(g+d,g+d-2,independent,9);
        assert(a==b && c>a);
        printf("{\"degree\":%ld,\"columns\":%ld,\"rank_8\":%ld,\"rank_9_redundant\":%ld,\"rank_9_independent\":%ld,\"remaining_8\":%ld,\"remaining_9\":%ld}\n",
            d,g[d].s,a,b,c,g[d].s-a,g[d].s-c);
    }
    clear_grades(g);flint_free(g);flint_free(independent);flint_free(redundant);
    flint_rand_clear(rng);flint_cleanup_master();return 0;
}

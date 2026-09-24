/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "../determinant/fq_mpoly_mat_det.c"
#include <assert.h>

static void verify(const mq_shared_support *previous, const mq_shared_support *shifts,
                   const mq_det_filter *filter, const nmod_mpoly_ctx_t ctx, slong n,
                   mq_shared_support *actual)
{
    mq_shared_support reference;
    uint32_t *expected=mq_shared_next(&reference,previous,shifts,filter,ctx,1,0), *map=NULL;
    assert(mq_rank_next(actual,&map,previous,shifts,filter,ctx,n,2,1,1,1));
    for(size_t i=0;i<previous->count*shifts->count;i++) {
        assert((expected[i]==UINT32_MAX)==(map[i]==UINT32_MAX));
        if(map[i]!=UINT32_MAX) {
            assert(map[i]<actual->count);
            assert(!memcmp(actual->keys+(size_t)map[i]*actual->words,
                           reference.keys+(size_t)expected[i]*reference.words,
                           actual->words*sizeof(ulong)));
        }
    }
    for(size_t i=0;i<actual->count;i++)
        assert(mq_shared_keep(actual->keys+i*actual->words,filter,ctx,actual->bits));
    flint_free(map); flint_free(expected); mq_shared_support_clear(&reference);
}
static void check(slong m)
{
    slong n=m+1,nv=2*m+1;
    slong *targets=flint_calloc(2*m,sizeof(slong));
    targets[0]=targets[m-1]=1;
    targets[m+1]=1;
    mq_det_filter filter; assert(mq_filter_init(&filter,m,targets,2,targets,2));
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx,nv,ORD_LEX,101);
    mq_filter_prepare_packed(&filter,ctx,n+1);
    mq_shared_support prev,shifts;
    mq_shared_support_init(&prev,ctx,n+1);
    mq_shared_support_init(&shifts,ctx,n+1);
    ulong exp[FLINT_BITS]={0},key[3]={0};
    for(slong v=0;v<nv;v++) {
        exp[v]=1; mpoly_set_monomial_ui(key,exp,prev.bits,ctx->minfo);
        mq_shared_insert(&prev,key); mq_shared_insert(&shifts,key); exp[v]=0;
    }
    /* Sources outside the ideal must be safely rejected, including those
     * arriving from an unfiltered hash fallback. */
    exp[0]=3; mpoly_set_monomial_ui(key,exp,prev.bits,ctx->minfo);
    mq_shared_insert(&prev,key);
    mq_shared_support ranked, hashed, again;
    verify(&prev,&shifts,&filter,ctx,n,&ranked);
    uint32_t *map=mq_shared_next(&hashed,&ranked,&shifts,NULL,ctx,1,0);
    flint_free(map);
    verify(&hashed,&shifts,&filter,ctx,n,&again);
    mq_shared_support sentinel; memset(&sentinel,0,sizeof(sentinel));
    sentinel.count=123; map=NULL;
    assert(!mq_rank_next(&sentinel,&map,&prev,&shifts,&filter,ctx,n,2,1,SIZE_MAX,0));
    assert(sentinel.count==123 && !map);
    assert(!mq_rank_next(&sentinel,&map,&prev,&shifts,NULL,ctx,n,2,1,1,0));
    mq_shared_support_clear(&again); mq_shared_support_clear(&hashed);
    mq_shared_support_clear(&ranked); mq_shared_support_clear(&prev);
    mq_shared_support_clear(&shifts); mq_filter_clear(&filter);
    nmod_mpoly_ctx_clear(ctx); flint_free(targets);
}
int main(void)
{
    omp_set_num_threads(4);
    check(3); check(6); check(9); /* Native one-, two-, three-word layouts. */
    flint_cleanup_master(); puts("MQ direct-index transitions and fallback PASS"); return 0;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Direct x/y/parameter indices for projected shared MQ DP. Axis keys retain
 * native multiword FLINT packing; no new single-word exponent restriction. */
typedef struct { ulong key[3]; slong degree; } mq_rank_axis;
static int mq_rank_compare(const void *a, const void *b)
{
    const mq_rank_axis *x=a, *y=b;
    for (slong w=3; w-- >0;)
        if(x->key[w]!=y->key[w]) return x->key[w]<y->key[w] ? -1 : 1;
    return 0;
}
static uint32_t mq_rank_find(const mq_rank_axis *axis, size_t count, const ulong *key)
{
    mq_rank_axis query={{key[0],key[1],key[2]},0};
    size_t lo=0, hi=count;
    while(lo<hi) {
        size_t mid=lo+(hi-lo)/2;
        int cmp=mq_rank_compare(axis+mid,&query);
        if(cmp<0) lo=mid+1; else hi=mid;
    }
    return lo<count && !mq_rank_compare(axis+lo,&query) ? (uint32_t)lo : UINT32_MAX;
}
static int mq_rank_enabled(void)
{
#ifdef DRSOLVE_MQ_LAYOUT_TEST
    /* Only mode 12 follows the default; preserve explicit hash ablations. */
    if(mq_shared_test_variant!=8) return mq_shared_test_variant==9;
#endif
    return g_dixon_mq_step1_rank;
}
/* Failure leaves next/map untouched; the caller uses ordinary hash building.
 * Extra envelope slots contain zero coefficients until reached by recurrence. */
static int mq_rank_next(mq_shared_support *next, uint32_t **map_out,
                        const mq_shared_support *prev, const mq_shared_support *shifts,
                        const mq_det_filter *f, const nmod_mpoly_ctx_t ctx,
                        slong n, slong k, size_t minors, size_t previous_minors, int parallel)
{
    if(!f) return 0;
    const size_t budget=DRSOLVE_MQ_SHARED_WORKSPACE_BYTES;
    slong nv=ctx->minfo->nvars, m=f->nvars, degree=0;
    ulong active[FLINT_BITS]={0}, exp[FLINT_BITS];
    /* The admitted simplex estimate also bounds this envelope. Its variables
     * are a subset of the processed row union, with the same degree budget. */
    const mq_shared_support *sets[2]={prev,shifts};
    for(int set=0;set<2;set++) {
        slong maximum=0;
        for(size_t i=0;i<sets[set]->count;i++) {
            slong d=0;
            mpoly_get_monomial_ui(exp,sets[set]->keys+i*prev->words,prev->bits,ctx->minfo);
            for(slong v=0;v<nv;v++) { active[v]|=exp[v]; d+=exp[v]; }
            maximum=FLINT_MAX(maximum,d);
        }
        degree+=maximum;
    }
    if(!active[nv-1]) return 0;
    size_t allocated_axes=(size_t)f->rows.count+f->cols.count;
    if(allocated_axes>budget/(8*sizeof(mq_rank_axis))) return 0;
    mq_rank_axis *axis[2]; size_t length[2]={0,0};
    const mq_monom_set *ideals[2]={&f->rows,&f->cols};
    for(int side=0;side<2;side++) {
        axis[side]=flint_malloc(ideals[side]->count*sizeof(mq_rank_axis));
        for(slong i=0;i<ideals[side]->alloc;i++) if(ideals[side]->keys[i]) {
            ulong code=ideals[side]->keys[i]-1;
            memset(exp,0,nv*sizeof(ulong));
            slong d=0; int allowed=1;
            for(slong v=0;v<m;v++) {
                ulong e=(code>>(v*f->bits))&f->digit_mask;
                exp[side*m+v]=e; d+=e;
                if(e && !active[side*m+v]) allowed=0;
            }
            if(!allowed || d>degree) continue;
            mq_rank_axis *entry=axis[side]+length[side]++;
            memset(entry,0,sizeof(*entry)); entry->degree=d;
            mpoly_set_monomial_ui(entry->key,exp,prev->bits,ctx->minfo);
        }
        qsort(axis[side],length[side],sizeof(mq_rank_axis),mq_rank_compare);
    }
    size_t nx=length[0], ny=length[1], count=0;
    if(nx && ny>budget/sizeof(uint32_t)/nx) goto reject;
    for(size_t x=0;x<nx;x++) for(size_t y=0;y<ny;y++) {
        slong d=degree-axis[0][x].degree-axis[1][y].degree;
        if(d>=0) count+=(size_t)d+1;
    }
    if(count>=UINT32_MAX || count>budget/sizeof(ulong)) goto reject;
    /* Account conservatively for simultaneous coefficients, packing allowance,
     * retained keys, full transition map, and all direct-index scratch. All
     * multiplicands below are bounded by admission or the checks above. */
    long double bytes=65536.0L+mq_shared_support_bytes(prev)+mq_shared_support_bytes(shifts)
        +(long double)previous_minors*prev->count*sizeof(ulong)
        +(long double)minors*count*sizeof(ulong)
        +(long double)n*shifts->count*sizeof(ulong)
        +(long double)shifts->count*prev->count*sizeof(uint32_t)
        +(long double)count*prev->words*sizeof(ulong)
        +(long double)allocated_axes*sizeof(mq_rank_axis)
        +(long double)nx*ny*sizeof(uint32_t)
        +(long double)3*prev->count*sizeof(uint32_t)
        +(long double)shifts->count*(nx+ny+1)*sizeof(uint32_t);
    if(k==n) bytes+=2.0L*count*(prev->words+1)*sizeof(ulong);
    if(bytes>budget) goto reject;
    uint32_t *base=flint_malloc(nx*ny*sizeof(uint32_t));
    size_t pos=0;
    for(size_t x=0;x<nx;x++) for(size_t y=0;y<ny;y++) {
        slong d=degree-axis[0][x].degree-axis[1][y].degree;
        base[x*ny+y]=d<0 ? UINT32_MAX : (uint32_t)pos;
        if(d>=0) pos+=d+1;
    }
    memset(next,0,sizeof(*next));
    next->bits=prev->bits; next->words=prev->words;
    next->count=next->capacity=count;
    next->keys=flint_malloc(count*prev->words*sizeof(ulong));
    slong tw,ts;
    mpoly_gen_offset_shift_sp(&tw,&ts,nv-1,prev->bits,ctx->minfo);
    #pragma omp parallel for if(parallel) schedule(static)
    for(size_t x=0;x<nx;x++) for(size_t y=0;y<ny;y++) {
        slong d=degree-axis[0][x].degree-axis[1][y].degree;
        for(slong t=0;t<=d;t++) {
            ulong *key=next->keys+((size_t)base[x*ny+y]+t)*prev->words;
            for(slong w=0;w<prev->words;w++) key[w]=axis[0][x].key[w]+axis[1][y].key[w];
            key[tw]+=(ulong)t<<ts;
        }
    }
    uint32_t *source=flint_malloc(3*prev->count*sizeof(uint32_t));
    for(size_t b=0;b<prev->count;b++) {
        ulong ex[FLINT_BITS], key[3]={0};
        mpoly_get_monomial_ui(exp,prev->keys+b*prev->words,prev->bits,ctx->minfo);
        for(int side=0;side<2;side++) {
            memset(ex,0,nv*sizeof(ulong)); memcpy(ex+side*m,exp+side*m,m*sizeof(ulong));
            mpoly_set_monomial_ui(key,ex,prev->bits,ctx->minfo);
            source[3*b+side]=mq_rank_find(axis[side],length[side],key);
            /* A preceding hash fallback may retain zero padding outside the
             * ideal. Such sources cannot reach a target and are dropped. */
        }
        source[3*b+2]=exp[nv-1];
    }
    uint32_t *xm=flint_malloc(shifts->count*nx*sizeof(uint32_t));
    uint32_t *ym=flint_malloc(shifts->count*ny*sizeof(uint32_t));
    uint32_t *tm=flint_malloc(shifts->count*sizeof(uint32_t));
    for(size_t a=0;a<shifts->count;a++) {
        mpoly_get_monomial_ui(exp,shifts->keys+a*prev->words,prev->bits,ctx->minfo);
        tm[a]=exp[nv-1];
        for(int side=0;side<2;side++) {
            ulong ex[FLINT_BITS]={0}, key[3]={0};
            memcpy(ex+side*m,exp+side*m,m*sizeof(ulong));
            mpoly_set_monomial_ui(key,ex,prev->bits,ctx->minfo);
            uint32_t *dest=(side?ym:xm)+a*length[side];
            for(size_t i=0;i<length[side];i++) {
                ulong sum[3];
                for(slong w=0;w<3;w++) sum[w]=axis[side][i].key[w]+key[w];
                dest[i]=mq_rank_find(axis[side],length[side],sum);
            }
        }
    }
    uint32_t *map=flint_malloc(shifts->count*prev->count*sizeof(uint32_t));
    #pragma omp parallel for if(parallel) schedule(static)
    for(size_t a=0;a<shifts->count;a++) for(size_t b=0;b<prev->count;b++) {
        if(source[3*b]==UINT32_MAX || source[3*b+1]==UINT32_MAX) {
            map[a*prev->count+b]=UINT32_MAX; continue;
        }
        uint32_t x=xm[a*nx+source[3*b]], y=ym[a*ny+source[3*b+1]], t=tm[a]+source[3*b+2];
        uint32_t id=UINT32_MAX;
        if(x!=UINT32_MAX && y!=UINT32_MAX &&
           axis[0][x].degree+axis[1][y].degree+t<=degree)
            id=base[(size_t)x*ny+y]+t;
        map[a*prev->count+b]=id;
    }
    *map_out=map;
    flint_free(base); flint_free(source); flint_free(xm); flint_free(ym); flint_free(tm);
    flint_free(axis[0]); flint_free(axis[1]);
    return 1;
reject:
    flint_free(axis[0]); flint_free(axis[1]); return 0;
}

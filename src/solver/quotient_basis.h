/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Shared degree-ordered monomial indexing; private to the quotient backends. */
#ifndef DRSOLVE_QUOTIENT_BASIS_H
#define DRSOLVE_QUOTIENT_BASIS_H
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <flint/flint.h>

typedef struct {
    slong n, degree, count, capacity, slots;
    ulong *exps;
    slong *degrees, *lookup;
} basis_t;

static uint64_t code(const ulong *e, slong n)
{
    uint64_t h = UINT64_C(1469598103934665603);
    for (slong i = 0; i < n; i++) { h ^= e[i]; h *= UINT64_C(1099511628211); }
    return h;
}

static slong hash_slot(uint64_t key, slong slots)
{
    key ^= key>>30; key *= UINT64_C(0xbf58476d1ce4e5b9);
    key ^= key>>27; key *= UINT64_C(0x94d049bb133111eb);
    return (slong)((key^(key>>31)) & (uint64_t)(slots-1));
}

static void basis_enumerate(basis_t *b, ulong *e, slong pos, slong left, slong degree)
{
    if (pos==b->n-1) {
        assert(b->count<b->capacity); e[pos]=left;
        memcpy(b->exps+b->count*b->n,e,b->n*sizeof(ulong));
        b->degrees[b->count++]=degree; return;
    }
    for (slong i=left;i>=0;i--) { e[pos]=i; basis_enumerate(b,e,pos+1,left-i,degree); }
}

static void basis_init(basis_t *b, slong n, slong degree)
{
    memset(b,0,sizeof(*b)); b->n=n; b->degree=degree;
    slong count=1;
    for (slong i=1;i<=degree;i++) count=count*(n+i)/i;
    b->capacity=count;
    b->exps=flint_malloc(count*n*sizeof(ulong)); b->degrees=flint_malloc(count*sizeof(slong));
    ulong e[n]; memset(e,0,sizeof(e));
    for (slong d=degree;d>=0;d--) basis_enumerate(b,e,0,d,d);
    assert(b->count==count); b->slots=1;
    while (b->slots<2*count) b->slots*=2;
    b->lookup=flint_malloc(b->slots*sizeof(slong));
    for (slong j=0;j<b->slots;j++) b->lookup[j]=-1;
    for (slong j=0;j<count;j++) {
        uint64_t k=code(b->exps+j*n,n); slong slot=hash_slot(k,b->slots);
        while (b->lookup[slot]>=0) slot=(slot+1)&(b->slots-1);
        b->lookup[slot]=j;
    }
}

static slong basis_find(const basis_t *b, const ulong *e)
{
    uint64_t k=code(e,b->n); slong slot=hash_slot(k,b->slots);
    while (b->lookup[slot]>=0) {
        if (!memcmp(b->exps+b->lookup[slot]*b->n,e,b->n*sizeof(ulong))) return b->lookup[slot];
        slot=(slot+1)&(b->slots-1);
    }
    return -1;
}

static void basis_clear(basis_t *b)
{
    flint_free(b->exps); flint_free(b->degrees); flint_free(b->lookup);
}

static slong monomial_count_bounded(slong n, slong d, slong limit)
{
    slong count=1;
    for(slong i=1;i<=d;i++) {
        if(n>WORD_MAX-i || count>WORD_MAX/(n+i)) return 0;
        count=count*(n+i)/i;
        if(count>limit) return 0;
    }
    return count;
}

#endif

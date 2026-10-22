/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "gf2n_native.h"
#include <pthread.h>

/* Wrappers retain their own FLINT view. Only immutable arithmetic data is
 * shared, and it is released with the last live context for that modulus. */
typedef struct gf2n_precomputed {
    gf2n_native_ctx data;
    size_t references;
    struct gf2n_precomputed *next;
} gf2n_precomputed;
static gf2n_precomputed *precomputed;
static pthread_mutex_t precomputed_lock=PTHREAD_MUTEX_INITIALIZER;

static gf2128_t add(gf2128_t a, gf2128_t b)
{ return (gf2128_t){a.low^b.low,a.high^b.high}; }

static gf2128_t times_x(gf2128_t a, const gf2n_native_ctx *c)
{
    unsigned k=c->n-1;
    uint64_t carry=k<64 ? (a.low>>k)&1 : (a.high>>(k-64))&1;
    a=(gf2128_t){a.low<<1,(a.high<<1)|(a.low>>63)};
    a.low &= c->mask.low; a.high &= c->mask.high;
    return carry ? add(a,c->modulus) : a;
}

static gf2128_t clmul_soft(uint64_t a, uint64_t b)
{
    gf2128_t r={0,0};
    while (b) {
        unsigned s=__builtin_ctzll(b); b &= b-1;
        r.low ^= a<<s;
        if (s) r.high ^= a>>(64-s);
    }
    return r;
}
#if DIXON_X86_SIMD
__attribute__((target("pclmul,sse2")))
static gf2128_t clmul_hw(uint64_t a, uint64_t b)
{
    __m128i p=_mm_clmulepi64_si128(_mm_set_epi64x(0,a),_mm_set_epi64x(0,b),0);
    gf2128_t r; _mm_storeu_si128((__m128i *)&r,p); return r;
}
#endif
static gf2128_t clmul(uint64_t a,uint64_t b,const gf2n_native_ctx *c)
{
#if DIXON_X86_SIMD
    if(c->clmul) return clmul_hw(a,b);
#endif
    return clmul_soft(a,b);
}
static void product(uint64_t p[4],gf2128_t a,gf2128_t b,const gf2n_native_ctx *c)
{
    gf2128_t lo=clmul(a.low,b.low,c);
    p[0]=lo.low; p[1]=lo.high; p[2]=p[3]=0;
    if(c->n<=64) return;
    gf2128_t hi=clmul(a.high,b.high,c);
    gf2128_t cross=clmul(a.low^a.high,b.low^b.high,c);
    cross=add(cross,add(lo,hi));
    p[1]^=cross.low; p[2]=hi.low^cross.high; p[3]=hi.high;
}
static gf2128_t high_part(const uint64_t p[4],unsigned n)
{
    unsigned w=n/64,s=n%64;
    if(!s) return (gf2128_t){p[w],p[w+1]};
    return (gf2128_t){(p[w]>>s)|(p[w+1]<<(64-s)),
                      (p[w+1]>>s)|(p[w+2]<<(64-s))};
}
static gf2128_t reduce(uint64_t p[4],const gf2n_native_ctx *c)
{
    gf2128_t low={p[0]&c->mask.low,p[1]&c->mask.high};
    gf2128_t high=high_part(p,c->n);
    if(c->fold) {
        while(high.low || high.high) {
            product(p,high,c->modulus,c);
            low.low ^= p[0]&c->mask.low; low.high ^= p[1]&c->mask.high;
            high=high_part(p,c->n);
        }
    } else {
        for(unsigned i=0;i<(c->n+7)/8;i++) {
            unsigned b=i<8 ? (high.low>>(8*i))&255 : (high.high>>(8*(i-8)))&255;
            low=add(low,c->reduction[i][b]);
        }
    }
    return low;
}
gf2128_t gf2n_native_mul(gf2128_t a,gf2128_t b,const gf2n_native_ctx *c)
{
    if(c->n<=16) return (gf2128_t){gf2n_native_mul_word(a.low,b.low,c),0};
    uint64_t p[4]; product(p,a,b,c); return reduce(p,c);
}
gf2128_t gf2n_native_sqr(gf2128_t a,const gf2n_native_ctx *c)
{
    if(c->n<=16) return gf2n_native_mul(a,a,c);
    gf2128_t lo=clmul(a.low,a.low,c),hi=clmul(a.high,a.high,c);
    uint64_t p[4]={lo.low,lo.high,hi.low,hi.high}; return reduce(p,c);
}
gf2128_t gf2n_native_inv(gf2128_t a,const gf2n_native_ctx *c)
{
    if(!(a.low|a.high)) return a;
    if(c->inverse) return (gf2128_t){c->inverse[a.low],0};
    /* Binary addition chain for a^(2^(n-1)-1), then one square.
     * At most floor(log2(n-1))+popcount(n-1)-1 multiplications. */
    unsigned target=c->n-1,k=1,bit=1;
    if(!target) return a;
    while(bit<=target/2) bit<<=1;
    gf2128_t r=a;
    while((bit>>=1)) {
        gf2128_t t=r;
        for(unsigned i=0;i<k;i++) t=gf2n_native_sqr(t,c);
        r=gf2n_native_mul(r,t,c); k*=2;
        if(target&bit) { r=gf2n_native_mul(gf2n_native_sqr(r,c),a,c); k++; }
    }
    return gf2n_native_sqr(r,c);
}

gf2128_t gf2n_native_from_fq(const fq_nmod_t a)
{
    gf2128_t r={0,0};
    for(slong i=0;i<a->length;i++) if(a->coeffs[i]&1) {
        if(i<64) r.low|=UINT64_C(1)<<i;
        else r.high|=UINT64_C(1)<<(i-64);
    }
    return r;
}
void gf2n_native_to_fq(fq_nmod_t out,gf2128_t a,const fq_nmod_ctx_t fq)
{
    fq_nmod_zero(out,fq);
    while(a.low) { unsigned i=__builtin_ctzll(a.low); a.low&=a.low-1; nmod_poly_set_coeff_ui(out,i,1); }
    while(a.high) { unsigned i=__builtin_ctzll(a.high); a.high&=a.high-1; nmod_poly_set_coeff_ui(out,64+i,1); }
}
static uint16_t small_mul(uint16_t a,uint16_t b,const gf2n_native_ctx *c)
{
    gf2128_t x={a,0},r={0,0};
    while(b) { if(b&1) r=add(r,x); b>>=1; x=times_x(x,c); }
    return r.low;
}
gf2n_native_ctx *gf2n_native_create(const fq_nmod_ctx_t fq)
{
    unsigned n=fq_nmod_ctx_degree(fq);
    FLINT_ASSERT(fq_nmod_ctx_prime(fq)==2 && n>=1 && n<=128);
    gf2n_native_ctx *c=flint_calloc(1,sizeof(*c));
    c->fq=*fq; c->n=n; c->clmul=has_pclmulqdq();
    c->mask.low=n>=64 ? UINT64_MAX : (UINT64_C(1)<<n)-1;
    c->mask.high=n==128 ? UINT64_MAX : n>64 ? (UINT64_C(1)<<(n-64))-1 : 0;
    const nmod_poly_struct *f=fq_nmod_ctx_modulus(fq);
    unsigned degree=0;
    for(unsigned i=0;i<n;i++) if(nmod_poly_get_coeff_ui(f,i)) {
        degree=i;
        if(i<64) c->modulus.low|=UINT64_C(1)<<i;
        else c->modulus.high|=UINT64_C(1)<<(i-64);
    }
    pthread_mutex_lock(&precomputed_lock);
    for(gf2n_precomputed *p=precomputed;p;p=p->next)
        if(p->data.n==n && p->data.modulus.low==c->modulus.low &&
           p->data.modulus.high==c->modulus.high) {
            *c=p->data; c->fq=*fq; c->shared=p; p->references++;
            pthread_mutex_unlock(&precomputed_lock);
            return c;
        }
    if(n<=16) {
        unsigned q=1U<<n,order=q-1;
        c->log=flint_calloc(q,sizeof(uint16_t));
        c->exp=flint_calloc(2*q,sizeof(uint16_t));
        c->inverse=flint_calloc(q,sizeof(uint16_t));
        /* Find a generator in this exact polynomial basis. */
        for(unsigned g=1;g<q;g++) {
            unsigned x=1,i=0;
            do { c->exp[i]=x; c->log[x]=i; x=small_mul(x,g,c); i++; } while(x!=1 && i<order);
            if(i==order && x==1) break;
        }
        for(unsigned i=order;i<2*q;i++) c->exp[i]=c->exp[i-order];
        for(unsigned a=1;a<q;a++) c->inverse[a]=c->exp[order-c->log[a]];
        if(n<=8) {
            c->byte_mul=flint_calloc(q*q,1);
            c->byte_simd=flint_calloc(256,sizeof(*c->byte_simd));
            for(unsigned a=1;a<q;a++) for(unsigned b=1;b<q;b++)
                c->byte_mul[a*q+b]=c->exp[c->log[a]+c->log[b]];
            for(unsigned a=0;a<q;a++) for(unsigned b=0;b<16;b++) {
                c->byte_simd[a][b]=b<q ? c->byte_mul[a*q+b] : 0;
                c->byte_simd[a][16+b]=(b<<4)<q ? c->byte_mul[a*q+(b<<4)] : 0;
            }
        }
    } else {
        c->fold=degree<n/2;
        if(!c->fold) {
            c->reduction=flint_calloc((n+7)/8,sizeof(*c->reduction));
            gf2128_t power=c->modulus;
            for(unsigned i=0;i<(n+7)/8;i++) for(unsigned bit=0;bit<8;bit++) {
                for(unsigned j=0;j<(1U<<bit);j++)
                    c->reduction[i][j|(1U<<bit)]=add(c->reduction[i][j],power);
                power=times_x(power,c);
            }
        }
    }
    gf2n_precomputed *p=flint_malloc(sizeof(*p));
    p->data=*c;
    memset(&p->data.fq,0,sizeof(p->data.fq));
    p->references=1; p->next=precomputed; precomputed=p; c->shared=p;
    pthread_mutex_unlock(&precomputed_lock);
    return c;
}
void gf2n_native_destroy(gf2n_native_ctx *c)
{
    if(!c) return;
    gf2n_precomputed *p=c->shared;
    pthread_mutex_lock(&precomputed_lock);
    if(--p->references==0) {
        gf2n_precomputed **link=&precomputed;
        while(*link!=p) link=&(*link)->next;
        *link=p->next;
        flint_free(p->data.reduction); flint_free(p->data.byte_mul); flint_free(p->data.byte_simd);
        flint_free(p->data.log); flint_free(p->data.exp); flint_free(p->data.inverse); flint_free(p);
    }
    pthread_mutex_unlock(&precomputed_lock);
    flint_free(c);
}

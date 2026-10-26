/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "bezout_prefix.h"
#include <flint/nmod.h>

/* Leading coefficients only, in reverse order. */
static void reverse_high(nmod_poly_t out, const nmod_poly_t p, slong limit)
{
    slong len=FLINT_MIN(p->length,limit);
    nmod_poly_fit_length(out,len);
    for(slong i=0;i<len;i++)out->coeffs[i]=p->coeffs[p->length-1-i];
    out->length=len;_nmod_poly_normalise(out);
}

/* Given exact degrees and leading prefixes of t0,t1, form a*t0-b*t1.
 * Reject possible leading cancellation: knowing a short prefix alone would
 * not justify guessing the resulting degree or missing trailing coefficients. */
static int combine(nmod_poly_t out, slong *degree,
    const nmod_poly_t a, const nmod_poly_t t0, slong d0,
    const nmod_poly_t b, const nmod_poly_t t1, slong d1,
    slong limit, nmod_poly_t reverse, nmod_poly_t product)
{
    slong da=a->length && d0>=0 ? nmod_poly_degree(a)+d0 : -1;
    slong db=b->length && d1>=0 ? nmod_poly_degree(b)+d1 : -1;
    if(da>=0 && da==db)return 0;
    *degree=FLINT_MAX(da,db);nmod_poly_zero(out);
    if(da>=0 && *degree-da<limit) {
        slong shift=*degree-da;
        reverse_high(reverse,a,limit-shift);
        nmod_poly_mullow(product,reverse,t0,limit-shift);
        nmod_poly_shift_left(out,product,shift);
    }
    if(db>=0 && *degree-db<limit) {
        slong shift=*degree-db;
        reverse_high(reverse,b,limit-shift);
        nmod_poly_mullow(product,reverse,t1,limit-shift);
        nmod_poly_shift_left(product,product,shift);
        nmod_poly_sub(out,out,product);
    }
    return *degree<0 || (out->length && out->coeffs[0]);
}

int dr_bezout_inverse_prefix(nmod_poly_t out, const nmod_poly_t f,
    const nmod_poly_t g, slong precision, slong cutoff, dr_bezout_prefix_stats *stats)
{
    dr_bezout_prefix_stats local={0,0,0};if(!stats)stats=&local;*stats=local;
    slong n=nmod_poly_degree(f);
    if(n<1 || nmod_poly_degree(g)>=n || precision<1 || cutoff<0 ||
       out->mod.n!=f->mod.n || g->mod.n!=f->mod.n)return -1;
    nmod_poly_t u,v,q,r,t0,t1,new0,new1,rev,prod,answer;
    nmod_poly_t m11,m12,m21,m22,A,B,one;
    nmod_poly_struct *all[]={u,v,q,r,t0,t1,new0,new1,rev,prod,answer,
                            m11,m12,m21,m22,A,B,one};
    for(size_t i=0;i<sizeof(all)/sizeof(*all);i++)nmod_poly_init_mod(all[i],f->mod);
    nmod_poly_set(u,f);nmod_poly_set(v,g);nmod_poly_one(t1);nmod_poly_one(one);
    slong d0=-1,d1=0;int ok=0;
    while(v->length>1) {
        if(cutoff && u->length>=cutoff && nmod_poly_degree(v)>=((nmod_poly_degree(u)+1)/2)) {
            /* FLINT returns [u;v]=M*[A;B], det(M)=sign.
             * Propagate just the leading prefixes via M^{-1}. The internal
             * HGCD matrix is full; the cumulative cofactors remain truncated. */
            slong sign=nmod_poly_hgcd(m11,m12,m21,m22,A,B,u,v);
            slong nd0,nd1;stats->hgcd_blocks++;
            if(!combine(new0,&nd0,m22,t0,d0,m12,t1,d1,precision,rev,prod) ||
               !combine(new1,&nd1,m21,t0,d0,m11,t1,d1,precision,rev,prod))goto fallback;
            if(sign<0)nmod_poly_neg(new0,new0);
            if(sign>0)nmod_poly_neg(new1,new1);
            nmod_poly_swap(t0,new0);nmod_poly_swap(t1,new1);d0=nd0;d1=nd1;
            nmod_poly_swap(u,A);nmod_poly_swap(v,B);
            if(!v->length) {
                if(u->length!=1)goto done;
                nmod_poly_swap(u,v);nmod_poly_swap(t0,t1);d1=d0;
            }
        } else {
            nmod_poly_divrem(q,r,u,v);stats->euclidean_steps++;
            if(!r->length)goto done;
            slong next;
            if(!combine(new1,&next,one,t0,d0,q,t1,d1,precision,rev,prod))goto fallback;
            nmod_poly_swap(t0,t1);nmod_poly_swap(t1,new1);d0=d1;d1=next;
            nmod_poly_swap(u,v);nmod_poly_swap(v,r);
        }
    }
    if(!v->length || d1<0 || d1>n-1)goto done;
    if(n-1-d1<precision) {
        nmod_poly_truncate(t1,precision-(n-1-d1));
        nmod_poly_scalar_mul_nmod(answer,t1,nmod_inv(v->coeffs[0],f->mod));
        nmod_poly_shift_left(answer,answer,n-1-d1);
    }
    ok=1;goto publish;
fallback:
    stats->full_inverse_fallbacks++;
    ok=nmod_poly_invmod(answer,g,f);
    if(ok){nmod_poly_reverse(answer,answer,n);nmod_poly_truncate(answer,precision);}
publish:
    if(ok)nmod_poly_swap(out,answer);
done:
    for(size_t i=0;i<sizeof(all)/sizeof(*all);i++)nmod_poly_clear(all[i]);
    return ok;
}

int dr_bezout_left_moments(nmod_poly_t out, const nmod_poly_t f,
    const nmod_poly_t g, slong precision, slong cutoff, dr_bezout_prefix_stats *stats)
{
    nmod_poly_t prefix,rev,result;
    nmod_poly_init_mod(prefix,f->mod);nmod_poly_init_mod(rev,f->mod);nmod_poly_init_mod(result,f->mod);
    int ok=dr_bezout_inverse_prefix(prefix,f,g,precision,cutoff,stats);
    if(out->mod.n!=f->mod.n)ok=-1;
    if(ok==1) {
        reverse_high(rev,f,precision);
        nmod_poly_div_series(result,prefix,rev,precision);
        nmod_poly_swap(out,result);
    }
    nmod_poly_clear(prefix);nmod_poly_clear(rev);nmod_poly_clear(result);
    return ok;
}

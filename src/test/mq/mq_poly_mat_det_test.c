/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_poly_mat_det.h"
#include <assert.h>
#include <stdio.h>
#include <flint/ulong_extras.h>

static void check_coefficients(const nmod_poly_mat_t B,const slong *rd,const slong *cd,slong h,slong sigma)
{
    mq_coeff_mat a={0};assert(mq_coeff_mat_init(&a,B->r,rd,cd,sigma,B->modulus));
    for(slong i=0;i<B->r;i++)for(slong j=0;j<B->c;j++) {
        const nmod_poly_struct *p=nmod_poly_mat_entry(B,i,j);
        assert(p->length<=FLINT_MAX(0,sigma-rd[i]-cd[j]+1));
        for(slong d=0;d<p->length;d++)a.values[a.offset[i*B->r+j]+d]=p->coeffs[d];
    }
    nmod_poly_mat_t expected,actual;nmod_poly_mat_init(expected,h,h,B->modulus);
    nmod_poly_mat_init(actual,h,h,B->modulus);
    nmod_poly_one(nmod_poly_mat_entry(expected,0,0));nmod_poly_mat_set(actual,expected);
    ulong ef=42,af=42;
    int want=nmod_poly_mat_mq_schur(expected,&ef,B,rd,cd,h,sigma);
    int got=nmod_poly_mat_mq_schur_coeffs(actual,&af,&a,h);
    assert(want==got && ef==af && nmod_poly_mat_equal(expected,actual));
    assert(!nmod_poly_mat_mq_schur_coeffs(actual,&af,&a,0));
    nmod_poly_mat_clear(expected);nmod_poly_mat_clear(actual);mq_coeff_mat_clear(&a);
}

/* Complement [[a(t),b],[c,0]]: nonconstant E with constant nonzero det.
 * This checks exact scaling/sign, including characteristic two, as well as
 * safe rejection without corrupting the caller's fallback input/output. */
static void check_degree_blocks(ulong prime)
{
    const slong n=15,h=3;
    slong rd[15]={0,1,2},cd[15]={0,1,2};
    for(slong i=h;i<n;i++) {
        rd[i]=i<h+3?1:i<h+9?2:3;
        cd[i]=4-rd[i];
    }
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,prime,971);
    nmod_poly_mat_t B,saved,core,E;
    nmod_poly_mat_init(B,n,n,prime);nmod_poly_mat_init(saved,n,n,prime);
    nmod_poly_mat_init(core,h,h,prime);nmod_poly_mat_init(E,n-h,n-h,prime);
    for(slong i=0;i<n;i++)for(slong j=0;j<n;j++) {
        slong bound=4-rd[i]-cd[j];
        for(slong k=0;k<=bound;k++)nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B,i,j),k,n_randint(rng,prime));
        if(i>=h && j>=h && bound==0) {
            if(i>j)nmod_poly_zero(nmod_poly_mat_entry(B,i,j));
            if(i==j)nmod_poly_one(nmod_poly_mat_entry(B,i,j));
        }
    }
    /* Independently permute both axes, including the degree groups. */
    for(slong i=n-1;i>h;i--) {
        slong r=h+n_randint(rng,i-h+1),c=h+n_randint(rng,i-h+1),tmp;
        for(slong j=0;j<n;j++)nmod_poly_swap(nmod_poly_mat_entry(B,i,j),nmod_poly_mat_entry(B,r,j));
        for(slong j=0;j<n;j++)nmod_poly_swap(nmod_poly_mat_entry(B,j,i),nmod_poly_mat_entry(B,j,c));
        tmp=rd[i];rd[i]=rd[r];rd[r]=tmp;tmp=cd[i];cd[i]=cd[c];cd[c]=tmp;
    }
    check_coefficients(B,rd,cd,h,4);
    nmod_poly_mat_set(saved,B);
    nmod_poly_t expected,got,factor_poly;
    nmod_poly_init(expected,prime);nmod_poly_init(got,prime);nmod_poly_init(factor_poly,prime);
    nmod_poly_mat_det(expected,B);
    for(slong i=h;i<n;i++)for(slong j=h;j<n;j++)nmod_poly_set(nmod_poly_mat_entry(E,i-h,j-h),nmod_poly_mat_entry(B,i,j));
    nmod_poly_mat_det(factor_poly,E);
    assert(nmod_poly_degree(factor_poly)==0);
    ulong factor=0;
    assert(nmod_poly_mat_mq_schur(core,&factor,B,rd,cd,h,4));
    assert(factor==nmod_poly_get_coeff_ui(factor_poly,0));
    nmod_poly_mat_det(got,core);nmod_poly_scalar_mul_nmod(got,got,factor);
    assert(nmod_poly_equal(got,expected) && nmod_poly_mat_equal(B,saved));
    /* Empty complement returns the original core and unit factor. */
    assert(nmod_poly_mat_mq_schur(saved,&factor,B,rd,cd,n,4));
    assert(factor==1 && nmod_poly_mat_equal(saved,B));
    nmod_poly_clear(expected);nmod_poly_clear(got);nmod_poly_clear(factor_poly);
    nmod_poly_mat_clear(E);nmod_poly_mat_clear(core);nmod_poly_mat_clear(saved);nmod_poly_mat_clear(B);
    flint_rand_clear(rng);
}

static void check_adaptive(void)
{
    flint_rand_t state;flint_rand_init(state);
    const ulong primes[]={2,101,1073741831};
    for(int trial=0;trial<90;trial++) {
        slong n=3+trial%8, sigma=trial%9==0 ? 0 : 6;
        slong rd[10],cd[10];
        nmod_poly_mat_t B,saved,core,old;
        nmod_poly_mat_init(B,n,n,primes[trial%3]);
        nmod_poly_mat_init(saved,n,n,primes[trial%3]);
        nmod_poly_mat_init(core,1,1,primes[trial%3]);
        nmod_poly_mat_init(old,1,1,primes[trial%3]);
        nmod_poly_one(nmod_poly_mat_entry(core,0,0));nmod_poly_mat_set(old,core);
        for(slong i=0;i<n;i++) {
            rd[i]=sigma ? n_randint(state,7) : 0;
            cd[i]=sigma ? n_randint(state,7) : 0;
        }
        for(slong i=0;i<n;i++)for(slong j=0;j<n;j++) {
            slong degree=sigma-rd[i]-cd[j];
            for(slong k=0;k<=degree;k++) {
                /* Include completely singular degree-diagonal layers. */
                ulong v=(trial%4==0 && degree==0) ? 0 : n_randint(state,primes[trial%3]);
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B,i,j),k,v);
            }
        }
        nmod_poly_mat_set(saved,B);
        ulong factor=42;
        int ok=nmod_poly_mat_mq_schur_adaptive(core,&factor,B,rd,cd,sigma);
        assert(nmod_poly_mat_equal(B,saved));
        if(ok) {
            assert(core->r<n);
            nmod_poly_t expected,got;nmod_poly_init(expected,B->modulus);nmod_poly_init(got,B->modulus);
            nmod_poly_mat_det(expected,B);nmod_poly_mat_det(got,core);
            nmod_poly_scalar_mul_nmod(got,got,factor);
            assert(nmod_poly_equal(expected,got));
            nmod_poly_clear(expected);nmod_poly_clear(got);
        } else assert(factor==42 && nmod_poly_mat_equal(core,old));
        /* Invalid degree metadata must fail without touching output. */
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B,0,0),sigma+1,1);
        nmod_poly_mat_clear(old);nmod_poly_mat_init(old,core->r,core->c,B->modulus);nmod_poly_mat_set(old,core);
        ulong before=factor;
        assert(!nmod_poly_mat_mq_schur_adaptive(core,&factor,B,rd,cd,sigma));
        assert(factor==before && nmod_poly_mat_equal(core,old));
        nmod_poly_mat_clear(B);nmod_poly_mat_clear(saved);nmod_poly_mat_clear(core);nmod_poly_mat_clear(old);
    }
    flint_rand_clear(state);
    puts("Adaptive Schur: randomized exact determinants and rejection preservation PASS");
}

int main(void)
{
    check_adaptive();
    ulong primes[] = {2, 3, 257, 65537};
    for (slong z = 0; z < 4; z++)
    {
        check_degree_blocks(primes[z]);
        ulong p = primes[z], factor = 0;
        slong d[] = {1, 1, 2};
        nmod_poly_mat_t B, saved, core;
        nmod_poly_mat_init(B, 3, 3, p);
        nmod_poly_mat_init(saved, 3, 3, p);
        nmod_poly_mat_init(core, 1, 1, p);
        for (slong i = 0; i < 3; i++)
            for (slong j = 0; j < 3; j++)
            {
                slong degree = 3 - d[i] - d[j];
                for (slong k = 0; k <= degree; k++)
                    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, i, j), k,
                                           (i + 2 * j + k + 1) % p);
            }
        nmod_poly_one(nmod_poly_mat_entry(B, 1, 2));
        nmod_poly_one(nmod_poly_mat_entry(B, 2, 1));
        check_coefficients(B,d,d,1,3);
        nmod_poly_mat_set(saved, B);
        nmod_poly_t expected, got;
        nmod_poly_init(expected, p);
        nmod_poly_init(got, p);
        nmod_poly_mat_det(expected, B);
        assert(nmod_poly_mat_mq_schur(core, &factor, B, d, d, 1, 3));
        nmod_poly_mat_det(got, core);
        nmod_poly_scalar_mul_nmod(got, got, factor);
        assert(nmod_poly_equal(expected, got));
        assert(nmod_poly_mat_equal(B, saved));
        assert(factor == p - 1);
        nmod_poly_zero(nmod_poly_mat_entry(B, 2, 1));
        check_coefficients(B,d,d,1,3);
        nmod_poly_mat_set(saved, B);
        nmod_poly_one(nmod_poly_mat_entry(core, 0, 0));
        factor = 42;
        assert(!nmod_poly_mat_mq_schur(core, &factor, B, d, d, 1, 3));
        assert(factor == 42 && nmod_poly_is_one(nmod_poly_mat_entry(core, 0, 0)) &&
               nmod_poly_mat_equal(B, saved));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, 2, 2), 1, 1);
        assert(!nmod_poly_mat_mq_schur(core, &factor, B, d, d, 1, 3));
        assert(factor == 42 && nmod_poly_is_one(nmod_poly_mat_entry(core, 0, 0)));
        assert(!nmod_poly_mat_mq_schur(B, &factor, B, d, d, 3, 3));
        /* B can be nonsingular even when this specified E is singular. */
        nmod_poly_mat_zero(B);
        nmod_poly_one(nmod_poly_mat_entry(B, 0, 2));
        nmod_poly_one(nmod_poly_mat_entry(B, 2, 0));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, 1, 1), 1, 1);
        nmod_poly_mat_det(expected, B);
        assert(nmod_poly_degree(expected) == 1);
        check_coefficients(B,d,d,1,3);
        nmod_poly_mat_set(saved, B);
        assert(!nmod_poly_mat_mq_schur(core, &factor, B, d, d, 1, 3));
        assert(nmod_poly_mat_equal(B, saved) && factor == 42);
        /* Fixed monomial ordering fails, but same-degree pivots repair both
         * axes. Check the exact determinant, not just nonsingularity. */
        assert(nmod_poly_mat_mq_schur_repair(core, &factor, B, d, d, 1, 3));
        nmod_poly_mat_det(got, core);
        nmod_poly_scalar_mul_nmod(got, got, factor);
        assert(nmod_poly_equal(expected, got) && nmod_poly_mat_equal(B, saved));
        /* Only a row swap is needed here: its odd sign must be included. */
        nmod_poly_mat_zero(B);
        nmod_poly_one(nmod_poly_mat_entry(B, 0, 2));
        nmod_poly_one(nmod_poly_mat_entry(B, 2, 1));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, 1, 0), 1, 1);
        check_coefficients(B,d,d,1,3);
        nmod_poly_mat_set(saved, B);
        nmod_poly_mat_det(expected, B);
        assert(!nmod_poly_mat_mq_schur(core, &factor, B, d, d, 1, 3));
        assert(nmod_poly_mat_mq_schur_repair(core, &factor, B, d, d, 1, 3));
        nmod_poly_mat_det(got, core);
        nmod_poly_scalar_mul_nmod(got, got, factor);
        assert(nmod_poly_equal(expected, got) && nmod_poly_mat_equal(B, saved));
        /* A genuinely deficient full layer cannot be repaired. Preserve all
         * caller-owned outputs on this failure and on degree violations. */
        nmod_poly_zero(nmod_poly_mat_entry(B, 2, 1));
        check_coefficients(B,d,d,1,3);
        nmod_poly_mat_set(saved, B);
        nmod_poly_one(nmod_poly_mat_entry(core, 0, 0)); factor = 42;
        assert(!nmod_poly_mat_mq_schur_repair(core, &factor, B, d, d, 1, 3));
        assert(factor == 42 && nmod_poly_is_one(nmod_poly_mat_entry(core, 0, 0)) &&
               nmod_poly_mat_equal(B, saved));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, 2, 2), 1, 1);
        assert(!nmod_poly_mat_mq_schur_repair(core, &factor, B, d, d, 1, 3));
        assert(factor == 42 && nmod_poly_is_one(nmod_poly_mat_entry(core, 0, 0)));
        nmod_poly_clear(expected);
        nmod_poly_clear(got);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(saved);
        nmod_poly_mat_clear(core);
    }
    puts("MQ Schur: exact determinant/scaling, characteristic two, rejection and input "
         "preservation PASS");
    flint_cleanup_master();
    return 0;
}

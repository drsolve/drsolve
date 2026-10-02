/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include <stdlib.h>
#include <flint/nmod_poly_mat.h>
#include "../../pml_det/pml_det.h"
#include "../../pml_det/src/nmod_poly_mat_multiply.h"

static void window_reference(nmod_poly_mat_t out, const nmod_poly_mat_t a,
                             const nmod_poly_mat_t b, slong start, slong count)
{
    nmod_poly_mat_mul_classical(out,a,b);
    nmod_poly_mat_shift_right(out,out,start);
    nmod_poly_mat_truncate(out,count);
}

static void windows(flint_rand_t state)
{
    const ulong primes[] = {257,65537,998244353,
#if FLINT_BITS == 64
        UWORD(18446744069414584321),
#endif
    };
    for (slong p = 0; p < (slong)(sizeof(primes)/sizeof(*primes)); p++)
        for (slong run = 0; run < 90; run++) {
            slong m=1+run%5, k=1+run%4, n=1+run%3;
            slong start=run%33, count=run%39;
            nmod_poly_mat_t a,b,c,ref;
            nmod_poly_mat_init(a,m,k,primes[p]); nmod_poly_mat_init(b,k,n,primes[p]);
            nmod_poly_mat_init(c,m,n,primes[p]); nmod_poly_mat_init(ref,m,n,primes[p]);
            nmod_poly_mat_randtest(a,state,run%53);
            nmod_poly_mat_randtest(b,state,run%37);
            if (run%11==0) {
                for (slong i=0;i<m;i++) for (slong j=0;j<k;j++)
                    for (slong d=0;d<53;d++)
                        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(a,i,j),d,primes[p]-1);
            }
            window_reference(ref,a,b,start,count);
            assert(drsolve_nmod_poly_mat_mul_window_ntt(c,a,b,start,count));
            assert(nmod_poly_mat_equal(c,ref));
            if (n==k) {
                assert(drsolve_nmod_poly_mat_mul_window_ntt(a,a,b,start,count));
                assert(nmod_poly_mat_equal(a,ref));
            } else if (m==k) {
                assert(drsolve_nmod_poly_mat_mul_window_ntt(b,a,b,start,count));
                assert(nmod_poly_mat_equal(b,ref));
            }
            nmod_poly_mat_clear(ref); nmod_poly_mat_clear(c);
            nmod_poly_mat_clear(b); nmod_poly_mat_clear(a);
        }
    /* In-place squaring, length-one transforms, and unsupported transforms. */
    for (int run=0;run<4;run++) {
        ulong prime=run==3?19:65537;
        nmod_poly_mat_t a,ref,copy;
        nmod_poly_mat_init(a,3,3,prime); nmod_poly_mat_init(ref,3,3,prime);
        nmod_poly_mat_init(copy,3,3,prime);
        nmod_poly_mat_randtest(a,state,run==0?1:31);
        nmod_poly_mat_set(copy,a);
        window_reference(ref,a,a,0,61);
        int ok=drsolve_nmod_poly_mat_mul_window_ntt(a,a,a,0,61);
        assert(ok==(run!=3));
        assert(nmod_poly_mat_equal(a,ok?ref:copy));
        nmod_poly_mat_clear(copy); nmod_poly_mat_clear(ref); nmod_poly_mat_clear(a);
    }
    /* Windows sharing storage with inputs are safe: all input transforms
     * finish before writing a single output coefficient. */
    nmod_poly_mat_t a,b,out,ref;
    nmod_poly_mat_init(a,5,5,65537); nmod_poly_mat_init(b,5,3,65537);
    nmod_poly_mat_init(ref,5,3,65537); nmod_poly_mat_randtest(a,state,45);
    nmod_poly_mat_randtest(b,state,51); nmod_poly_mat_window_init(out,a,0,0,5,3);
    window_reference(ref,a,b,17,29);
    assert(drsolve_nmod_poly_mat_mul_window_ntt(out,a,b,17,29));
    assert(nmod_poly_mat_equal(out,ref));
    nmod_poly_mat_window_clear(out); nmod_poly_mat_clear(ref);
    nmod_poly_mat_clear(b); nmod_poly_mat_clear(a);
}

static void geometric_bounds(flint_rand_t state)
{
    nmod_poly_mat_t a,b,c,ref;
    nmod_poly_mat_init(a,3,3,65537); nmod_poly_mat_init(b,3,2,65537);
    nmod_poly_mat_init(c,3,2,65537); nmod_poly_mat_init(ref,3,2,65537);
    /* PMBasis can retain input coefficients beyond its approximation order,
     * and a shifted approximant can have degree greater than start. */
    for (slong run=0;run<3;run++) {
        nmod_poly_mat_randtest(a,state,run==0?9:25);
        nmod_poly_mat_randtest(b,state,run==2?80:17);
        window_reference(ref,a,b,8,9);
        nmod_poly_mat_middle_product_geometric(c,a,b,8,8);
        assert(nmod_poly_mat_equal(c,ref));
    }
    nmod_poly_mat_clear(ref); nmod_poly_mat_clear(c);
    nmod_poly_mat_clear(b); nmod_poly_mat_clear(a);
}

static void determinants(flint_rand_t state)
{
    const ulong primes[]={2,3,17,257,65537,1000003};
    for (slong p=0;p<6;p++) for (slong run=0;run<5;run++) {
        slong n=9+run;
        nmod_poly_mat_t a,copy;
        nmod_poly_t got,ref;
        nmod_poly_mat_init(a,n,n,primes[p]); nmod_poly_mat_init(copy,n,n,primes[p]);
        nmod_poly_init(got,primes[p]); nmod_poly_init(ref,primes[p]);
        for (slong i=0;i<n;i++) for (slong j=0;j<n;j++)
            nmod_poly_randtest(nmod_poly_mat_entry(a,i,j),state,2+(run==4?4*(n-i):(n-i)/2));
        if (run==1) {
            /* An odd permutation and a nonconstant common row factor. */
            for(slong j=0;j<n;j++) nmod_poly_swap(nmod_poly_mat_entry(a,0,j),nmod_poly_mat_entry(a,1,j));
            for(slong j=0;j<n;j++) nmod_poly_shift_left(nmod_poly_mat_entry(a,0,j),nmod_poly_mat_entry(a,0,j),2);
        }
        if (run==2) for(slong j=0;j<n;j++)
            nmod_poly_set(nmod_poly_mat_entry(a,0,j),nmod_poly_mat_entry(a,1,j));
        if (run==3) for(slong j=0;j<n;j++) nmod_poly_zero(nmod_poly_mat_entry(a,0,j));
        nmod_poly_mat_set(copy,a);
        nmod_poly_mat_det_fflu(ref,a);
        int ok=nmod_poly_mat_det_hnf(got,a);
        /* The HNF backend may reject rank-deficient splits. The caller's
         * existing fallback handles rejection; an accepted answer must match. */
        if(ok) assert(nmod_poly_equal(got,ref));
        if(run==0 || run==1 || run==4) assert(ok);
        assert(nmod_poly_mat_equal(a,copy));
        nmod_poly_clear(ref); nmod_poly_clear(got);
        nmod_poly_mat_clear(copy); nmod_poly_mat_clear(a);
    }
}

int main(void)
{
    flint_rand_t state; flint_rand_init(state);
    windows(state); geometric_bounds(state); determinants(state);
    flint_rand_clear(state);
    puts("Prime Step 4: NTT windows/aliasing/fallback, geometric bounds, degree-order determinants PASS");
    return 0;
}

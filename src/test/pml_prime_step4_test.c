/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <flint/nmod_poly_mat.h>
#include "../../pml_det/pml_det.h"
#include "../../pml_det/src/nmod_poly_mat_multiply.h"
#include "../../pml_det/src/nmod_mat_poly.h"
#include "../../pml_det/src/nmod_poly_mat_utils.h"

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

/* Large enough to enter the parallel dispatch, with coefficient-level
 * references and both direct and offset-window input/output aliasing. */
static void parallel_windows(flint_rand_t state)
{
    const ulong primes[] = {257,65537,998244353,
#if FLINT_BITS == 64
        UWORD(18446744069414584321),
#endif
    };
#ifdef _OPENMP
    int saved = omp_get_max_threads();
#endif
    for (size_t p=0;p<sizeof(primes)/sizeof(*primes);p++) {
        nmod_poly_mat_t a,b,c,ref,copy,window;
        nmod_poly_mat_init(a,24,24,primes[p]);
        nmod_poly_mat_init(b,24,24,primes[p]);
        nmod_poly_mat_init(c,24,24,primes[p]);
        nmod_poly_mat_init(ref,24,24,primes[p]);
        nmod_poly_mat_init(copy,25,25,primes[p]);
        nmod_poly_mat_window_init(window,copy,1,1,25,25);
        nmod_poly_mat_randtest(a,state,65);
        nmod_poly_mat_randtest(b,state,81);
        window_reference(ref,a,b,17,73);
        for (int threads=1;threads<=4;threads*=2) {
#ifdef _OPENMP
            omp_set_num_threads(threads);
#endif
            assert(drsolve_nmod_poly_mat_mul_window_ntt(c,a,b,17,73));
            assert(nmod_poly_mat_equal(c,ref));
            nmod_poly_mat_set(c,a);
            assert(drsolve_nmod_poly_mat_mul_window_ntt(c,c,b,17,73));
            assert(nmod_poly_mat_equal(c,ref));
            nmod_poly_mat_set(window,b);
            assert(drsolve_nmod_poly_mat_mul_window_ntt(window,a,window,17,73));
            assert(nmod_poly_mat_equal(window,ref));
        }
#ifdef _OPENMP
        /* Existing teams must not recursively create an NTT worker team. */
        #pragma omp parallel num_threads(2)
        {
            #pragma omp single
            {
                assert(drsolve_nmod_poly_mat_mul_window_ntt(c,a,b,17,73));
                assert(nmod_poly_mat_equal(c,ref));
            }
        }
#endif
        window_reference(ref,a,a,17,73);
        assert(drsolve_nmod_poly_mat_mul_window_ntt(a,a,a,17,73));
        assert(nmod_poly_mat_equal(a,ref));
        nmod_poly_mat_window_clear(window); nmod_poly_mat_clear(copy);
        nmod_poly_mat_clear(ref); nmod_poly_mat_clear(c);
        nmod_poly_mat_clear(b); nmod_poly_mat_clear(a);
    }
#ifdef _OPENMP
    omp_set_num_threads(saved);
#endif
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

/* Compare parallel coefficient-block updates with the serial basis, and
 * independently check that the basis annihilates the input to the order. */
static void parallel_mbasis(flint_rand_t state)
{
#ifdef _OPENMP
    int saved = omp_get_max_threads();
#endif
    const ulong primes[] = {19,257,65537};
    for (int run=0;run<6;run++) {
        const slong m=80,n=40,order=16;
        ulong prime=primes[run%3];
        nmod_poly_mat_t input,ref,got,product;
        nmod_mat_poly_t matp,app;
        slong shifts[80],expected[80];
        nmod_poly_mat_init(input,m,n,prime);
        nmod_poly_mat_init(ref,m,m,prime);
        nmod_poly_mat_init(got,m,m,prime);
        nmod_poly_mat_init(product,m,n,prime);
        nmod_poly_mat_randtest(input,state,order);
        if (run>=3) {
            for (slong j=0;j<n;j++) {
                nmod_poly_zero(nmod_poly_mat_entry(input,0,j));
                nmod_poly_set(nmod_poly_mat_entry(input,1,j),nmod_poly_mat_entry(input,2,j));
            }
        }
        nmod_mat_poly_init(matp,m,n,prime);
        nmod_mat_poly_set_trunc_from_poly_mat(matp,input,order);
        nmod_mat_poly_init(app,m,m,prime);
        for (int threads=1;threads<=4;threads*=2) {
#ifdef _OPENMP
            omp_set_num_threads(threads);
#endif
            for (slong i=0;i<m;i++) shifts[i]=run>=3 ? i%7-3 : 0;
            nmod_mat_poly_mbasis_resupdate(app,shifts,matp,order);
            nmod_poly_mat_set_from_mat_poly(got,app);
            if (threads==1) {
                nmod_poly_mat_set(ref,got);
                for (slong i=0;i<m;i++) expected[i]=shifts[i];
            } else {
                assert(nmod_poly_mat_equal(got,ref));
                for (slong i=0;i<m;i++) assert(expected[i]==shifts[i]);
            }
        }
        nmod_poly_mat_mul_classical(product,got,input);
        nmod_poly_mat_truncate(product,order);
        assert(nmod_poly_mat_is_zero(product));
        nmod_mat_poly_clear(app); nmod_mat_poly_clear(matp);
        nmod_poly_mat_clear(product); nmod_poly_mat_clear(got);
        nmod_poly_mat_clear(ref); nmod_poly_mat_clear(input);
    }
#ifdef _OPENMP
    omp_set_num_threads(saved);
#endif
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
    windows(state); parallel_windows(state); geometric_bounds(state); parallel_mbasis(state); determinants(state);
    flint_rand_clear(state);
    puts("Prime Step 4: NTT serial/parallel windows/aliasing/fallback, geometric bounds, parallel M-Basis, degree-order determinants PASS");
    return 0;
}

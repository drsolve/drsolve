/* Differential tests for full and arbitrary-window matrix products. */
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "../test_env.h"
#include <flint/nmod_poly_mat.h>
#include "nmod_poly_mat_multiply.h"
#ifdef _OPENMP
#include <omp.h>
#endif
static void trial(ulong p, slong m, slong k, slong n, slong alen, slong blen,
                  slong start, slong count, flint_rand_t rng)
{
    nmod_poly_mat_t a,b,c,ref,full,save;
    nmod_poly_mat_init(a,m,k,p); nmod_poly_mat_init(b,k,n,p);
    nmod_poly_mat_init(c,m,n,p); nmod_poly_mat_init(ref,m,n,p);
    nmod_poly_mat_init(full,m,n,p); nmod_poly_mat_init(save,m,n,p);
    nmod_poly_mat_randtest(a,rng,alen); nmod_poly_mat_randtest(b,rng,blen);
    nmod_poly_mat_mul(full,a,b);
    nmod_poly_mat_shift_right(ref,full,start); nmod_poly_mat_truncate(ref,count);
    for (int enabled=0;enabled<2;enabled++) {
        drsolve_pml_set_fft_enabled(enabled);
        drsolve_nmod_poly_mat_mul(c,a,b); assert(nmod_poly_mat_equal(c,full));
        drsolve_nmod_poly_mat_middle_product(c,a,b,start,count);
        assert(nmod_poly_mat_equal(c,ref));
        if(m==k && k==n) {
            nmod_poly_mat_set(save,a);
            drsolve_nmod_poly_mat_middle_product(save,save,b,start,count);
            assert(nmod_poly_mat_equal(save,ref));
            nmod_poly_mat_set(save,b);
            drsolve_nmod_poly_mat_middle_product(save,a,save,start,count);
            assert(nmod_poly_mat_equal(save,ref));
            nmod_poly_mat_set(save,a);
            drsolve_nmod_poly_mat_mul(save,save,b);
            assert(nmod_poly_mat_equal(save,full));
        }
    }
    nmod_poly_mat_clear(a);nmod_poly_mat_clear(b);nmod_poly_mat_clear(c);
    nmod_poly_mat_clear(ref);nmod_poly_mat_clear(full);nmod_poly_mat_clear(save);
}
int main(void)
{
    flint_rand_t rng; flint_rand_init(rng);flint_rand_set_seed(rng,345,789);
#ifdef _OPENMP
    flint_set_num_threads(omp_get_max_threads());
#endif
    /* Exercise even tiny FFT inputs, including with the Windows CRT. */
    if (dr_test_setenv("DRSOLVE_PML_MUL", "fft", 1) != 0) {
        perror("DRSOLVE_PML_MUL");
        return 1;
    }
    ulong primes[]={2,17,65537,998244353,1000000007,
#if FLINT_BITS == 64
        UWORD(18446744069414584321),
#endif
    };
    slong lens[]={0,1,7,63,128,257,513};
    int cases=0;
    for(unsigned pi=0;pi<sizeof(primes)/sizeof(*primes);pi++)
      for(unsigned j=0;j<sizeof(lens)/sizeof(*lens);j++)
        for(int shape=0;shape<2;shape++) {
          slong a=lens[j], b=2*a+1;
          slong starts[]={0,a/2,a,a+b+3};
          slong counts[]={0,3,a+1,2*a+7};
          for(int w=0;w<4;w++) {
            trial(primes[pi],shape?3:2,2,shape?4:2,a,b,starts[w],counts[w],rng);
            ++cases;
          }
        }
    trial(998244353,6,6,6,1025,2049,1024,1025,rng);++cases;
    trial(1000000007,6,6,6,1025,2049,1024,1025,rng);++cases;
    printf("PML FFT available=%d: %d cases, full/window/aliases/on-off PASS\n",drsolve_pml_fft_available(),cases);
    flint_rand_clear(rng);flint_cleanup_master();return 0;
}

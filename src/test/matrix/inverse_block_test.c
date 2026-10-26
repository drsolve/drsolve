/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include <stdio.h>
#include "inverse_block_det.h"
#include "dixon_flint.h"
#include <flint/ulong_extras.h>

/* Exhaustive degree assignment, independent of the Hungarian implementation. */
static slong brute(const nmod_poly_mat_t a, slong row, ulong used)
{
    if (row==a->r) return 0;
    slong best=0;
    for (slong j=0;j<a->c;j++) if (!(used & (UWORD(1)<<j))) {
        slong d=FLINT_MAX(0,nmod_poly_degree(nmod_poly_mat_entry(a,row,j)));
        best=FLINT_MAX(best,d+brute(a,row+1,used|(UWORD(1)<<j)));
    }
    return best;
}

int main(void)
{
#ifndef HAVE_PML
    puts("inverse-block: PML unavailable; backend disabled");return 0;
#else
    /* Pull the logging state from libdrsolve before the static PML archive,
       whose determinant routines reference it, and keep this test quiet. */
    g_dixon_verbose_level = 0;
    g_dixon_debug_mode = 0;
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,4321,9123);
    ulong primes[]={257,65537,1073741827,
#if FLINT_BITS == 64
        UWORD(18446744069414584321),
#endif
    };
    slong successes=0;
    for (size_t pi=0;pi<sizeof(primes)/sizeof(*primes);pi++) {
        ulong prime=primes[pi];
        for (slong n=2;n<=6;n++) {
            nmod_poly_mat_t a,saved;nmod_poly_mat_init(a,n,n,prime);nmod_poly_mat_init(saved,n,n,prime);
            nmod_poly_t ref,got,sentinel;nmod_poly_init(ref,prime);nmod_poly_init(got,prime);nmod_poly_init(sentinel,prime);
            nmod_poly_set_coeff_ui(sentinel,9,17);
            for (int shape=0;shape<3;shape++) {
                for (slong i=0;i<n;i++) for (slong j=0;j<n;j++) {
                    nmod_poly_struct *entry=nmod_poly_mat_entry(a,i,j);nmod_poly_zero(entry);
                    slong degree=shape==0 ? 3 : 1+n_randint(rng,7);
                    if (shape==2 && i>j) continue;
                    for (slong k=0;k<=degree;k++) nmod_poly_set_coeff_ui(entry,k,n_randint(rng,prime));
                }
                slong bound=dr_nmod_poly_mat_det_degree_bound(a);
                assert(bound==brute(a,0,0));
                nmod_poly_mat_set(saved,a);nmod_poly_mat_det(ref,a);
                assert(nmod_poly_degree(ref)<=bound);
                for (slong m=0;m<=n;m++) {
                    dr_inverse_block_stats stats;nmod_poly_set(got,sentinel);
                    int ok=dr_nmod_poly_mat_det_inverse_block_try(got,a,m,(m&1)?4:1,rng,&stats);
                    assert(nmod_poly_mat_equal(a,saved));
                    if (ok) {assert(nmod_poly_equal(got,ref));assert(stats.checks>=8);successes++;}
                    else {assert(nmod_poly_equal(got,sentinel));assert(shape!=0);}
                }
            }
            nmod_poly_mat_zero(a);
            for (slong i=0;i<n;i++) nmod_poly_set_coeff_ui(nmod_poly_mat_entry(a,i,i),1,1);
            /* A tiny inverse projection loses multiplicity for t I_n. */
            nmod_poly_set(got,sentinel);
            assert(!dr_nmod_poly_mat_det_inverse_block_try(got,a,1,2,rng,NULL));
            assert(nmod_poly_equal(got,sentinel));
            assert(dr_nmod_poly_mat_det_inverse_block_try(got,a,n,2,rng,NULL));
            nmod_poly_zero(ref);nmod_poly_set_coeff_ui(ref,n,1);assert(nmod_poly_equal(got,ref));
            nmod_poly_mat_zero(a);nmod_poly_set(got,sentinel);
            assert(!dr_nmod_poly_mat_det_inverse_block_try(got,a,0,1,rng,NULL));
            assert(nmod_poly_equal(got,sentinel));
            nmod_poly_one(nmod_poly_mat_entry(a,0,0));
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(a,0,1),2,1);
            assert(!dr_nmod_poly_mat_det_inverse_block_try(got,a,0,1,rng,NULL));
            assert(nmod_poly_equal(got,sentinel));
            assert(!dr_nmod_poly_mat_det_inverse_block_try(got,a,n+1,1,rng,NULL));
            assert(!dr_nmod_poly_mat_det_inverse_block_try(got,a,-1,1,rng,NULL));
            nmod_poly_clear(ref);nmod_poly_clear(got);nmod_poly_clear(sentinel);
            nmod_poly_mat_clear(a);nmod_poly_mat_clear(saved);
        }
    }
    for (ulong p=15;p<=17;p+=2) {
        nmod_poly_mat_t a;nmod_poly_mat_init(a,2,2,p);
        nmod_poly_t got;nmod_poly_init(got,p);nmod_poly_one(got);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(a,0,0),5,1);
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(a,1,1),5,1);
        assert(!dr_nmod_poly_mat_det_inverse_block_try(got,a,1,1,rng,NULL));
        assert(nmod_poly_is_one(got));nmod_poly_clear(got);nmod_poly_mat_clear(a);
    }
    assert(dr_inverse_block_choose(381,3125,8)<dr_inverse_block_choose(381,3125,1));
    assert(successes>100);printf("inverse-block: %ld exact matrix comparisons passed\n",successes);
    flint_rand_clear(rng);flint_cleanup();return 0;
#endif
}

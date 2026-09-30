/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "fq_poly_mat_det.h"
#include "dixon_flint.h"
#include "fq_multivariate_interpolation.h"
#include <assert.h>
#include <flint/ulong_extras.h>

static void check(ulong prime, slong n, slong degree, flint_rand_t rng)
{
    nmod_poly_mat_t a,saved;
    nmod_poly_mat_init(a,n,n,prime); nmod_poly_mat_init(saved,n,n,prime);
    for (slong i=0;i<n;i++) for (slong j=0;j<n;j++)
        nmod_poly_randtest(nmod_poly_mat_entry(a,i,j),rng,degree+1);
    nmod_poly_mat_set(saved,a);
    nmod_poly_t expected,got,sentinel;
    nmod_poly_init(expected,prime); nmod_poly_init(got,prime); nmod_poly_init(sentinel,prime);
    nmod_poly_set_coeff_ui(got,0,1); nmod_poly_set(sentinel,got);
    fq_nmod_poly_mat_det_set_method(FQ_NMOD_POLY_DET_METHOD_AUTO);
    nmod_poly_mat_det_fflu(expected,a);
    int admitted = !n || (ulong)(n*degree) < prime;
    int ok = dixon_nmod_poly_mat_det_interpolate(got,a,n*degree);
    assert(ok == admitted);
    assert(nmod_poly_equal(got,admitted ? expected : sentinel));
    if (admitted) {
        assert(dixon_nmod_poly_mat_det_interpolate(got,a,-1));
        assert(nmod_poly_equal(got,expected));
    }
    /* Explicit interpolation also handles small-field rejection by fallback. */
    assert(nmod_poly_mat_equal(a,saved));
    fq_nmod_poly_mat_det_set_method(FQ_NMOD_POLY_DET_METHOD_INTERP);
    dixon_nmod_poly_mat_det(got,a);
    if ((ulong)(n*degree) >= prime && n) {
        for (ulong x=0;x<prime;x++)
            assert(nmod_poly_evaluate_nmod(got,x)==nmod_poly_evaluate_nmod(expected,x));
        assert(nmod_poly_degree(got)<(slong)prime);
    } else assert(nmod_poly_equal(got,expected));
    nmod_poly_clear(sentinel); nmod_poly_clear(got); nmod_poly_clear(expected);
    nmod_poly_mat_clear(saved); nmod_poly_mat_clear(a);
}

static void check_generic_prime_interpolation(ulong prime)
{
    fq_nmod_ctx_t ctx; fq_nmod_ctx_init_ui(ctx,prime,1,"a");
    slong n = FLINT_MIN(prime,20);
    fq_nmod_t *xs = flint_malloc(n*sizeof(fq_nmod_t));
    fq_nmod_t *ys = flint_malloc(n*sizeof(fq_nmod_t));
    ulong *nx = flint_malloc(n*sizeof(ulong)), *ny = flint_malloc(n*sizeof(ulong));
    for (slong i=0;i<n;i++) {
        nx[i]=i; ny[i]=(i*i*i+2*i+1)%prime;
        fq_nmod_init(xs[i],ctx); fq_nmod_init(ys[i],ctx);
        fq_nmod_set_ui(xs[i],nx[i],ctx); fq_nmod_set_ui(ys[i],ny[i],ctx);
    }
    nmod_poly_t expected; nmod_poly_init(expected,prime);
    nmod_poly_interpolate_nmod_vec_fast(expected,nx,ny,n);
    fq_nmod_poly_t got; fq_nmod_poly_init(got,ctx);
    fq_lagrange_interpolation_optimized(got,xs,ys,n,ctx);
    assert(fq_nmod_poly_length(got,ctx)==expected->length);
    for (slong i=0;i<expected->length;i++)
        assert(nmod_poly_get_coeff_ui(got->coeffs+i,0)==nmod_poly_get_coeff_ui(expected,i));
    fq_nmod_poly_clear(got,ctx); nmod_poly_clear(expected);
    for (slong i=0;i<n;i++) { fq_nmod_clear(xs[i],ctx); fq_nmod_clear(ys[i],ctx); }
    flint_free(xs); flint_free(ys); flint_free(nx); flint_free(ny); fq_nmod_ctx_clear(ctx);
}

int main(void)
{
    g_dixon_verbose_level=0;
    check_generic_prime_interpolation(2);
    check_generic_prime_interpolation(13);
    check_generic_prime_interpolation(65537);
    flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng,772,913);
    ulong primes[]={2,13,101,65537};
    for (slong p=0;p<4;p++) for (slong n=0;n<=8;n++)
        for (slong d=0;d<=4;d++) check(primes[p],n,d,rng);
    nmod_poly_mat_t a,saved;
    nmod_poly_mat_init(a,2,2,65537); nmod_poly_mat_init(saved,2,2,65537);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(a,0,0),2,1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(a,1,1),1,1);
    nmod_poly_mat_set(saved,a);
    nmod_poly_t got; nmod_poly_init(got,65537); nmod_poly_one(got);
    fq_nmod_poly_mat_det_set_method(FQ_NMOD_POLY_DET_METHOD_AUTO);
    /* A bad caller bound fails verification without changing output/input. */
    assert(!dixon_nmod_poly_mat_det_interpolate(got,a,1));
    assert(nmod_poly_is_one(got) && nmod_poly_mat_equal(a,saved));
    assert(dixon_nmod_poly_mat_det_interpolate(got,a,3));
    assert(nmod_poly_degree(got)==3 && nmod_poly_get_coeff_ui(got,3)==1);
    /* Identically singular matrices are successful zero interpolants. */
    nmod_poly_zero(nmod_poly_mat_entry(a,1,1));
    assert(dixon_nmod_poly_mat_det_interpolate(got,a,3) && nmod_poly_is_zero(got));
    nmod_poly_clear(got); nmod_poly_mat_clear(saved); nmod_poly_mat_clear(a);
    fq_nmod_poly_mat_det_set_method(FQ_NMOD_POLY_DET_METHOD_AUTO);
    flint_rand_clear(rng);
    puts("Polynomial-matrix interpolation: FFLU equality, singular points, zero determinants, bounds, input preservation and small-field fallback PASS");
    flint_cleanup_master(); return 0;
}

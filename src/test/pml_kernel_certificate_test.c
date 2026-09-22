/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Keep private certificate counters and helper visible only in this executable. */
#include "../../pml_det/src/nmod_poly_mat_extra/kernel.c"
#include <assert.h>

static int certificate(const nmod_poly_mat_t A, slong *input_shift, slong order,
                        slong expected, int should_pass)
{
    nmod_poly_mat_t AT, PT, N;
    slong n = A->c, q = -1;
    slong *shift = flint_malloc((size_t) n * sizeof(slong));
    slong *degrees = flint_malloc((size_t) n * sizeof(slong));
    memcpy(shift, input_shift, (size_t) n * sizeof(slong));
    nmod_poly_mat_init(AT, n, A->r, A->modulus);
    nmod_poly_mat_init(PT, n, n, A->modulus);
    nmod_poly_mat_transpose(AT, A);
    nmod_poly_mat_pmbasis(PT, shift, AT, order);
    int hit = _nmod_zls_degree_rank_certificate(N, degrees, A, PT, input_shift,
                                               shift, order, &q, 1);
    assert(hit == should_pass);
    if (hit) {
        assert(q == expected);
        if (q) {
            nmod_poly_mat_t product;
            nmod_poly_mat_init(product, A->r, q, A->modulus);
            nmod_poly_mat_mul(product, A, N);
            assert(nmod_poly_mat_is_zero(product));
            nmod_poly_mat_clear(product);
            nmod_poly_mat_clear(N);
        }
    } else assert(q == -1);
    nmod_poly_mat_clear(PT); nmod_poly_mat_clear(AT);
    flint_free(degrees); flint_free(shift);
    return hit;
}

static void edge_cases(void)
{
    nmod_poly_mat_t A;
    slong shift[3] = {0,0,0};
    nmod_poly_mat_init(A, 2, 3, 257);
    certificate(A, shift, 1, 3, 1); /* Zero matrix: no rank evaluation needed. */
    nmod_poly_mat_clear(A);
    nmod_poly_mat_init(A, 2, 2, 257);
    nmod_poly_mat_one(A);
    certificate(A, shift, 1, 0, 1); /* Zero nullity: do not initialize N. */
    nmod_poly_mat_clear(A);
    nmod_poly_mat_init(A, 1, 2, 257);
    nmod_poly_one(nmod_poly_mat_entry(A,0,0));
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A,0,1),1,1);
    shift[1] = 1;
    slong insufficient = g_nmod_kernel_zls_profile.certificate_insufficient_rows;
    certificate(A, shift, 1, 0, 0); /* shift == order is NOT a zero-product proof. */
    assert(g_nmod_kernel_zls_profile.certificate_insufficient_rows == insufficient + 1);
    nmod_poly_zero(nmod_poly_mat_entry(A,0,0));
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(A,0,0),1,1);
    shift[1] = 0;
    slong invalid = g_nmod_kernel_zls_profile.certificate_invalid_bounds;
    certificate(A, shift, 1, 0, 0); /* [x,x], I is an order-1 approximant. */
    assert(g_nmod_kernel_zls_profile.certificate_invalid_bounds == invalid + 1);
    nmod_poly_mat_clear(A);
}

/* A=h*[I B].  Every kernel vector is [-B;I]v, so the bottom q x q block of
   any genuine polynomial kernel basis must be UNIMODULAR.  Merely checking
   A*N=0 and rational rank would miss spurious polynomial factors. */
static void structured_case(flint_rand_t state, ulong prime, slong m, slong q,
                            int bad_points, int shifted)
{
    slong n = m + q;
    nmod_poly_mat_t A, N, product, bottom;
    nmod_poly_t h, determinant;
    nmod_poly_init(h, prime); nmod_poly_init(determinant, prime);
    nmod_poly_one(h);
    if (bad_points) {
        /* (x-1)(x-2); over F_2 this vanishes at both base-field points. */
        nmod_poly_set_coeff_ui(h,0,2 % prime);
        nmod_poly_set_coeff_ui(h,1,(prime - (3 % prime)) % prime);
        nmod_poly_set_coeff_ui(h,2,1);
    }
    nmod_poly_mat_init(A,m,n,prime);
    nmod_poly_mat_init(N,n,n,prime);
    for (slong i=0;i<m;i++) {
        nmod_poly_set(nmod_poly_mat_entry(A,i,i),h);
        for(slong j=m;j<n;j++) {
            nmod_poly_randtest(nmod_poly_mat_entry(A,i,j),state,4);
            nmod_poly_mul(nmod_poly_mat_entry(A,i,j),nmod_poly_mat_entry(A,i,j),h);
        }
    }
    slong *shift = flint_malloc((size_t)n*sizeof(slong));
    slong *deg = flint_malloc((size_t)n*sizeof(slong));
    slong *actual = flint_malloc((size_t)n*sizeof(slong));
    nmod_poly_mat_column_degree(shift,A,NULL);
    for(slong i=0;i<n;i++) {
        shift[i]=FLINT_MAX(0,shift[i]) + (shifted ? (i*7)%5 : 0);
        deg[i] = -999;
    }
    slong dim=nmod_poly_mat_kernel_zls(N,deg,A,shift,3.0);
    assert(dim==q);
    nmod_poly_mat_init(product,m,n,prime);
    nmod_poly_mat_mul(product,A,N);
    assert(nmod_poly_mat_is_zero(product));
    nmod_poly_mat_column_degree(actual,N,shift);
    for(slong j=0;j<q;j++) assert(actual[j]==deg[j]);
    nmod_poly_mat_init(bottom,q,q,prime);
    for(slong i=0;i<q;i++) for(slong j=0;j<q;j++)
        nmod_poly_set(nmod_poly_mat_entry(bottom,i,j),nmod_poly_mat_entry(N,m+i,j));
    nmod_poly_mat_det(determinant,bottom);
    assert(nmod_poly_degree(determinant)==0 && !nmod_poly_is_zero(determinant));
    nmod_poly_mat_clear(bottom); nmod_poly_mat_clear(product);
    nmod_poly_mat_clear(N); nmod_poly_mat_clear(A);
    nmod_poly_clear(determinant); nmod_poly_clear(h);
    flint_free(actual); flint_free(deg); flint_free(shift);
}

int main(void)
{
    g_dixon_verbose_level=3; g_dixon_debug_mode=1;
    nmod_poly_mat_kernel_zls_profile_reset();
    edge_cases();
    nmod_poly_mat_kernel_zls_profile_reset();
    flint_rand_t state; flint_rand_init(state); flint_rand_set_seed(state,19371,87129);
    const ulong primes[4]={2,3,257,65537};
    for(slong t=0;t<160;t++)
        structured_case(state,primes[t%4],1+t%4,1+(t/4)%4,t%3==0,t%2);
    const nmod_poly_mat_kernel_zls_profile_t *p=&g_nmod_kernel_zls_profile;
    assert(p->certificate_hits>0 && p->certificate_rank_misses>0);
    assert(p->residual_matrices_skipped==p->certificate_hits);
    assert(p->residual_rows_skipped>=p->certified_kernel_rows);
    double self=0,exclusive=0; slong calls=0,hits=0;
    for(slong d=0;d<NMOD_ZLS_PROFILE_DEPTHS;d++) {
        self+=p->by_depth[d].self; exclusive+=p->by_depth[d].exclusive;
        calls+=p->by_depth[d].calls; hits+=p->by_depth[d].certificate_hits;
        assert(p->by_depth[d].self>=-1e-9 && p->by_depth[d].exclusive>=-1e-9);
    }
    assert(calls==p->zls_calls && hits==p->certificate_hits);
    assert(fabs(self-p->by_depth[0].inclusive)<1e-7);
    assert(fabs(self-p->zls_self_time)<1e-7);
    assert(fabs(exclusive-p->zls_exclusive_time)<1e-7);
    assert(g_nmod_kernel_zls_depth==0 && !g_nmod_zls_profile_frame);
    nmod_poly_mat_kernel_zls_profile_print();
    printf("Kernel certificate: 164 cases passed; hits=%ld, misses=%ld, residual rows skipped=%ld\n",
           p->certificate_hits,p->certificate_rank_misses,p->residual_rows_skipped);
    /* Profiling must not control the mathematical shortcut. */
    g_dixon_verbose_level=0; g_dixon_debug_mode=0;
    structured_case(state,257,3,2,0,1);
    structured_case(state,257,3,2,1,1);
    flint_rand_clear(state); flint_cleanup(); return 0;
}

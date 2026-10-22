/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "fq_nmod_roots.h"
#include "gf2n_native.h"
#include "../../solver/gf2n_root_poly.h"
#include <assert.h>

int main(void)
{
    flint_rand_t state;
    flint_rand_init(state);
    const slong ns[] = {1,4,8,13,16,32,63,64,65,127,128};
    for (size_t k=0; k<sizeof(ns)/sizeof(*ns); ++k) {
        printf("Testing n=%ld\n", ns[k]); fflush(stdout);
        nmod_poly_t modulus;
        nmod_poly_init(modulus, 2);
        nmod_poly_randtest_monic_irreducible(modulus, state, ns[k]+1);
        fq_nmod_ctx_t ctx;
        fq_nmod_ctx_init_modulus(ctx, modulus, "t");
        fq_nmod_poly_t f, linear;
        fq_nmod_poly_init(f, ctx);
        fq_nmod_poly_init(linear, ctx);
        fq_nmod_t r, v;
        fq_nmod_init(r, ctx);
        fq_nmod_init(v, ctx);
        /* Check the multiplication itself, including zero padding and both
         * hardware and software arithmetic, independently of root extraction. */
        gf2n_native_ctx *native=gf2n_native_create(ctx);
        gf2128_t *buf=flint_calloc(8*256,sizeof(*buf));
        gf2128_t *a=buf, *b=a+256, *product=b+256, *scratch=product+512;
        fq_nmod_poly_randtest(f,state,130,ctx);
        fq_nmod_poly_randtest(linear,state,194,ctx);
        for (slong i=0; i<f->length; ++i) a[i]=gf2n_native_from_fq(f->coeffs+i);
        for (slong i=0; i<linear->length; ++i) b[i]=gf2n_native_from_fq(linear->coeffs+i);
        fq_nmod_poly_mul(f,f,linear,ctx);
        for (int software=0; software<2; ++software) {
            native->clmul=software ? 0 : has_pclmulqdq();
            gf2n_root_mul(product,a,b,256,scratch,native);
            for (slong i=0; i<512; ++i) {
                gf2n_native_to_fq(v,product[i],ctx);
                fq_nmod_poly_get_coeff(r,f,i,ctx);
                assert(fq_nmod_equal(r,v,ctx));
            }
        }
        flint_free(buf); gf2n_native_destroy(native);
        for (int trial=0; trial<6; ++trial) {
            fq_nmod_poly_randtest(f, state, trial == 3 ? (ns[k] == 32 ? 1002 : ns[k] == 128 ? 1026 : 258) : 36, ctx);
            fq_nmod_randtest_not_zero(r, state, ctx);
            if (trial == 2) { /* Purely inseparable, including zero root. */
                fq_nmod_poly_zero(f, ctx);
                fq_nmod_poly_set_coeff(f, 64, r, ctx);
            } else if (trial == 5) {
                fq_nmod_poly_zero(f, ctx);
                fq_nmod_poly_set_coeff(f, ns[k] == 1 ? 8193 : 31, r, ctx);
                fq_nmod_poly_set_coeff(f, 0, r, ctx);
            } else if (trial == 4) {
                fq_nmod_poly_one(f, ctx);
                for (slong j=0; j<40; ++j) {
                    fq_nmod_randtest(r, state, ctx);
                    fq_nmod_poly_gen(linear, ctx);
                    fq_nmod_poly_set_coeff(linear, 0, r, ctx);
                    fq_nmod_poly_mul(f, f, linear, ctx);
                }
            } else {
                fq_nmod_poly_gen(linear, ctx);
                fq_nmod_poly_set_coeff(linear, 0, r, ctx);
                for (int j=0; j<3; ++j) fq_nmod_poly_mul(f, f, linear, ctx);
            }
            fq_nmod_poly_factor_t expected;
            fq_nmod_poly_factor_init(expected, ctx);
            fq_nmod_poly_roots(expected, f, 1, ctx);
            for (int multiplicity=0; multiplicity<=1; ++multiplicity) {
                fq_nmod_roots_t got;
                fq_nmod_roots_init(got, ctx);
                our_fq_nmod_poly_roots(got, f, multiplicity, ctx);
                assert(got->num == expected->num);
                for (slong i=0; i<got->num; ++i) {
                    slong match=-1;
                    for (slong j=0; j<expected->num; ++j) {
                        fq_nmod_poly_evaluate_fq_nmod(v, expected->poly+j, got->roots+i, ctx);
                        if (fq_nmod_is_zero(v,ctx)) match=j;
                    }
                    assert(match>=0);
                    assert(got->mult[i] == (multiplicity ? expected->exp[match] : 1));
                    for (slong j=0; j<i; ++j) assert(!fq_nmod_equal(got->roots+i,got->roots+j,ctx));
                }
                fq_nmod_roots_clear(got, ctx);
            }
            fq_nmod_poly_factor_clear(expected, ctx);
        }
        fq_nmod_clear(r,ctx); fq_nmod_clear(v,ctx);
        fq_nmod_poly_clear(f,ctx); fq_nmod_poly_clear(linear,ctx);
        fq_nmod_ctx_clear(ctx); nmod_poly_clear(modulus);
    }
    flint_rand_clear(state);
    puts("GF(2^n) roots and multiplicities match FLINT.");
    return 0;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_mpoly_interface.h"
#include <assert.h>

/* Exercise the native entry points as well as unified dispatch: a silent
   fallback to FLINT must not be enough to pass these tests. */
#define CHECK_NATIVE(FIELD) do { \
    FIELD##_mpoly_ctx_t nc; \
    FIELD##_mpoly_ctx_init(nc, nvars, order); \
    FIELD##_mpoly_t x, y, z; \
    FIELD##_mpoly_init(x, nc); FIELD##_mpoly_init(y, nc); FIELD##_mpoly_init(z, nc); \
    fq_nmod_mpoly_to_##FIELD##_mpoly(x, a, field, fc); \
    fq_nmod_mpoly_to_##FIELD##_mpoly(y, b, field, fc); \
    assert(FIELD##_mpoly_mul_sparse(z, x, y, nc)); \
    slong N = mpoly_words_per_exp(z->bits, nc->minfo); \
    ulong *mask = flint_malloc(N * sizeof(ulong)); \
    mpoly_get_cmpmask(mask, N, z->bits, nc->minfo); \
    for (slong k = 1; k < z->length; ++k) \
        assert(mpoly_monomial_cmp(z->exps + N*(k-1), z->exps + N*k, N, mask) > 0); \
    flint_free(mask); \
    FIELD##_mpoly_to_fq_nmod_mpoly(got, z, field, fc); \
    assert(fq_nmod_mpoly_equal(got, expected, fc)); \
    assert(FIELD##_mpoly_mul(z, x, y, nc)); \
    FIELD##_mpoly_to_fq_nmod_mpoly(got, z, field, fc); \
    assert(fq_nmod_mpoly_equal(got, expected, fc)); \
    assert(FIELD##_mpoly_mul_sparse(x, x, y, nc)); \
    FIELD##_mpoly_to_fq_nmod_mpoly(got, x, field, fc); \
    assert(fq_nmod_mpoly_equal(got, expected, fc)); \
    fq_nmod_mpoly_to_##FIELD##_mpoly(x, a, field, fc); \
    assert(FIELD##_mpoly_mul(y, x, y, nc)); \
    FIELD##_mpoly_to_fq_nmod_mpoly(got, y, field, fc); \
    assert(fq_nmod_mpoly_equal(got, expected, fc)); \
    fq_nmod_mpoly_to_##FIELD##_mpoly(x, a, field, fc); \
    fq_nmod_mpoly_mul_johnson(got, a, a, fc); \
    assert(FIELD##_mpoly_mul_sparse(x, x, x, nc)); \
    FIELD##_mpoly_to_fq_nmod_mpoly(squared, x, field, fc); \
    assert(fq_nmod_mpoly_equal(got, squared, fc)); \
    FIELD##_mpoly_clear(x, nc); FIELD##_mpoly_clear(y, nc); FIELD##_mpoly_clear(z, nc); \
    FIELD##_mpoly_ctx_clear(nc); \
} while (0)

static void check_case(slong degree, uint64_t modulus_low, slong nvars,
                       ordering_t order, int limit_k)
{
    nmod_poly_t modulus;
    nmod_poly_init(modulus, 2);
    nmod_poly_set_coeff_ui(modulus, degree, 1);
    for (slong i = 0; i < 64; ++i)
        if ((modulus_low >> i) & 1) nmod_poly_set_coeff_ui(modulus, i, 1);
    fq_nmod_ctx_t field;
    fq_nmod_ctx_init_modulus(field, modulus, "a");
    nmod_poly_clear(modulus);
    field_ctx_t native_field;
    field_ctx_init(&native_field, field);
    unified_mpoly_ctx_t uc = unified_mpoly_ctx_init(nvars, order, &native_field);
    fq_nmod_mpoly_ctx_struct *fc = GET_FQ_CTX(uc);
    unified_mpoly_t u = unified_mpoly_init(uc), v = unified_mpoly_init(uc);
    unified_mpoly_t w = unified_mpoly_init(uc);
    fq_nmod_mpoly_t a, b, expected, got, squared;
    fq_nmod_mpoly_init(a, fc); fq_nmod_mpoly_init(b, fc);
    fq_nmod_mpoly_init(expected, fc); fq_nmod_mpoly_init(got, fc);
    fq_nmod_mpoly_init(squared, fc);
    flint_rand_t random;
    flint_rand_init(random);
    flint_rand_set_seed(random, degree * 111 + nvars, 1234 + order);
    gf2n_mpoly_set_array_limit_k(limit_k);
    for (int trial = 0; trial < 18; ++trial) {
        fq_nmod_mpoly_randtest_bound(a, random, 30, 4, fc);
        fq_nmod_mpoly_randtest_bound(b, random, 23, 4, fc);
        if (trial == 0) fq_nmod_mpoly_zero(a, fc);
        if (trial == 1) fq_nmod_mpoly_one(b, fc);
        if (trial == 2) fq_nmod_mpoly_set(b, a, fc);
        if (trial == 3) fq_nmod_mpoly_zero(b, fc);
        if (trial == 4) { fq_nmod_mpoly_one(a, fc); fq_nmod_mpoly_one(b, fc); }
        /* Large gaps; pack growth; >2 words; and >64-bit exponents. */
        if (trial >= 9) {
            flint_bitcnt_t bits[] = {7, 15, 31, 63, 64, 65, 80, 127, 130};
            fq_nmod_mpoly_randtest_bits(a, random, 13, bits[trial - 9], fc);
            fq_nmod_mpoly_randtest_bits(b, random, 17,
                                      trial % 2 ? 5 : bits[trial - 9], fc);
        }
        fq_nmod_mpoly_mul_johnson(expected, a, b, fc);
        switch (degree) {
            case 4: assert(native_field.field_id == FIELD_ID_GF24); CHECK_NATIVE(gf24); break;
            case 8: assert(native_field.field_id == FIELD_ID_GF28); CHECK_NATIVE(gf28); break;
            case 16: assert(native_field.field_id == FIELD_ID_GF216); CHECK_NATIVE(gf216); break;
            case 32: assert(native_field.field_id == FIELD_ID_GF232); CHECK_NATIVE(gf232); break;
            case 64: assert(native_field.field_id == FIELD_ID_GF264); CHECK_NATIVE(gf264); break;
            case 128: assert(native_field.field_id == FIELD_ID_GF2128); CHECK_NATIVE(gf2128); break;
            default: abort();
        }
        fq_nmod_mpoly_set(GET_FQ_POLY(u), a, fc);
        fq_nmod_mpoly_set(GET_FQ_POLY(v), b, fc);
        assert(unified_mpoly_mul(w, u, v));
        assert(fq_nmod_mpoly_equal(GET_FQ_POLY(w), expected, fc));
        assert(unified_mpoly_mul(u, u, v));
        assert(fq_nmod_mpoly_equal(GET_FQ_POLY(u), expected, fc));
    }
    flint_rand_clear(random);
    fq_nmod_mpoly_clear(a, fc); fq_nmod_mpoly_clear(b, fc);
    fq_nmod_mpoly_clear(expected, fc); fq_nmod_mpoly_clear(got, fc);
    fq_nmod_mpoly_clear(squared, fc);
    unified_mpoly_clear(u); unified_mpoly_clear(v); unified_mpoly_clear(w);
    unified_mpoly_ctx_clear(uc);
    field_ctx_clear(&native_field);
    fq_nmod_ctx_clear(field);
}

int main(void)
{
    const slong degrees[] = {4, 8, 16, 32, 64, 128};
    const uint64_t moduli[] = {0x13, 0x11d, 0x1002d, 0x100008299,
                              UINT64_C(0x247f43cb7), 0x87};
    const ordering_t orders[] = {ORD_LEX, ORD_DEGLEX, ORD_DEGREVLEX};
    const slong variables[] = {1, 3, 13, 21};
    for (int f = 0; f < 6; ++f) {
        for (int o = 0; o < 3; ++o)
            for (int n = 0; n < 4; ++n)
                for (int limit = -1; limit <= 0; ++limit)
                    check_case(degrees[f], moduli[f], variables[n], orders[o], limit);
        printf("GF(2^%ld) sparse/automatic multiplication: PASS\n", degrees[f]);
    }
    /* GF(256) also supports a non-native defining polynomial via conversion. */
    check_case(8, 0x11b, 3, ORD_DEGREVLEX, 0);
    check_case(8, 0x11d, 3, ORD_LEX, 0);
    gf2n_mpoly_set_array_limit_k(-1);
    flint_cleanup_master();
    return 0;
}

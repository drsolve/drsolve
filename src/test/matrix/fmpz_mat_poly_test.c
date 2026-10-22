/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_mat_poly.h"
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_mat.h>
#include <assert.h>

static void fill(unified_poly_mat_t a, slong len, flint_rand_t rng)
{
    field_elem_u v;
    field_init_elem(&v, a->ctx->field_id, field_ctx_data(a->ctx));
    for (slong i = 0; i < a->r; i++)
        for (slong j = 0; j < a->c; j++)
            for (slong k = 0; k < len; k++) {
                fmpz_randm(&v.big, rng, fmpz_mod_ctx_modulus(a->ctx->ctx.big_ctx));
                unified_poly_set_coeff(a->rows[i] + j, k, &v);
            }
    field_clear_elem(&v, a->ctx->field_id, field_ctx_data(a->ctx));
}
static void convert(fmpz_mod_poly_t out, const unified_poly_struct *in, const fmpz_mod_ctx_t ctx)
{
    fmpz_mod_poly_zero(out, ctx);
    for (slong k = 0; k < in->length; k++)
        fmpz_mod_poly_set_coeff_fmpz(out, k, &in->coeffs[k].big, ctx);
}
static void products(field_ctx_t *field, slong len, flint_rand_t rng)
{
    const fmpz_mod_ctx_struct *ctx = field->ctx.big_ctx;
    unified_poly_mat_t a, b, got;
    unified_poly_mat_init(a, 3, 2, field);
    unified_poly_mat_init(b, 2, 3, field);
    unified_poly_mat_init(got, 3, 3, field);
    fill(a, len, rng);
    fill(b, len + 3, rng);
    fmpz_mod_poly_t ref[9], x, y, t, z;
    fmpz_mod_poly_init(x, ctx);
    fmpz_mod_poly_init(y, ctx);
    fmpz_mod_poly_init(t, ctx);
    fmpz_mod_poly_init(z, ctx);
    for (slong i = 0; i < 3; i++)
        for (slong j = 0; j < 3; j++) {
            slong q = 3 * i + j;
            fmpz_mod_poly_init(ref[q], ctx);
            for (slong k = 0; k < 2; k++) {
                convert(x, a->rows[i] + k, ctx);
                convert(y, b->rows[k] + j, ctx);
                fmpz_mod_poly_mul(t, x, y, ctx);
                fmpz_mod_poly_add(ref[q], ref[q], t, ctx);
            }
        }
    unified_mat_poly_t ua, ub, uc;
    unified_mat_poly_init(ua, 3, 2, field);
    unified_mat_poly_init(ub, 2, 3, field);
    unified_mat_poly_init(uc, 3, 3, field);
    unified_mat_poly_from_unified(ua, a);
    unified_mat_poly_from_unified(ub, b);
    slong starts[] = {0, 1, len, 2 * len + 8};
    for (int test = 0; test < 5; test++) {
        slong start = test == 4 ? 0 : starts[test], count = test == 4 ? 0 : len + 1;
        unified_mat_poly_mul_window(uc, ua, ub, start, count);
        unified_mat_poly_to_unified(got, uc);
        for (slong i = 0; i < 3; i++)
            for (slong j = 0; j < 3; j++) {
                fmpz_mod_poly_shift_right(t, ref[3 * i + j], start, ctx);
                fmpz_mod_poly_truncate(t, count, ctx);
                convert(z, got->rows[i] + j, ctx);
                assert(fmpz_mod_poly_equal(t, z, ctx));
            }
    }
    unified_mat_poly_clear(uc);
    unified_mat_poly_clear(ub);
    unified_mat_poly_clear(ua);
    for (int i = 0; i < 9; i++)
        fmpz_mod_poly_clear(ref[i], ctx);
    fmpz_mod_poly_clear(x, ctx);
    fmpz_mod_poly_clear(y, ctx);
    fmpz_mod_poly_clear(t, ctx);
    fmpz_mod_poly_clear(z, ctx);
    unified_poly_mat_clear(got);
    unified_poly_mat_clear(b);
    unified_poly_mat_clear(a);
}
static void determinant(field_ctx_t *field, flint_rand_t rng)
{
    const fmpz_mod_ctx_struct *ctx = field->ctx.big_ctx;
    unified_poly_mat_t a, got;
    unified_poly_mat_init(a, 4, 4, field);
    unified_poly_mat_init(got, 4, 4, field);
    fill(a, 5, rng);
    unified_mat_poly_t ua;
    unified_mat_poly_init(ua, 4, 4, field);
    unified_mat_poly_from_unified(ua, a);
    unified_poly_t det;
    unified_poly_init(det, field);
    assert(unified_mat_poly_det_hnf(det, ua));
    fmpz_mod_poly_t f, d;
    fmpz_mod_poly_init(f, ctx);
    fmpz_mod_poly_init(d, ctx);
    convert(d, det, ctx);
    fmpz_mod_mat_t scalar;
    fmpz_mod_mat_init(scalar, 4, 4, ctx);
    fmpz_t point, expected, actual;
    fmpz_init(point);
    fmpz_init(expected);
    fmpz_init(actual);
    /* Degree <= 16: 17 distinct values certify the complete determinant. */
    assert(d->length <= 17);
    for (slong k = 0; k < 17; k++) {
        fmpz_set_si(point, k);
        for (slong i = 0; i < 4; i++)
            for (slong j = 0; j < 4; j++) {
                convert(f, a->rows[i] + j, ctx);
                fmpz_mod_poly_evaluate_fmpz(fmpz_mod_mat_entry(scalar, i, j), f, point, ctx);
            }
        fmpz_mod_mat_det(expected, scalar, ctx);
        fmpz_mod_poly_evaluate_fmpz(actual, d, point, ctx);
        assert(fmpz_equal(actual, expected));
    }
    /* In-place square uses the same ownership path. */
    unified_mat_poly_mul(ua, ua, ua);
    unified_mat_poly_to_unified(got, ua);
    assert(unified_mat_poly_det_hnf(det, ua));
    convert(f, det, ctx);
    fmpz_mod_poly_mul(d, d, d, ctx);
    assert(fmpz_mod_poly_equal(f, d, ctx));
    fmpz_clear(actual);
    fmpz_clear(expected);
    fmpz_clear(point);
    fmpz_mod_mat_clear(scalar, ctx);
    fmpz_mod_poly_clear(f, ctx);
    fmpz_mod_poly_clear(d, ctx);
    unified_poly_clear(det);
    unified_mat_poly_clear(ua);
    unified_poly_mat_clear(got);
    unified_poly_mat_clear(a);
}
int main(void)
{
    const char *primes[] = {"18446744073709551629", "170141183460469231731687303715884105727"};
    flint_rand_t rng;
    flint_rand_init(rng);
    flint_rand_set_seed(rng, 123, 456);
    for (int i = 0; i < 2; i++) {
        fmpz_t p;
        fmpz_init(p);
        fmpz_set_str(p, primes[i], 10);
        field_ctx_t field;
        field_ctx_init_prime(&field, p);
        products(&field, 0, rng);
        products(&field, 1, rng);
        products(&field, 17, rng);
        products(&field, 65, rng);
        determinant(&field, rng);
        field_ctx_clear(&field);
        fmpz_clear(p);
    }
    flint_rand_clear(rng);
    puts("fmpz matrix polynomial windows, aliasing and HNF passed");
}

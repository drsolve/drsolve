/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "dr_mpoly.h"
#include <assert.h>

static void equal_native(const unified_mpoly_struct *actual, const fq_nmod_mpoly_t expected,
                         fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpoly_t got;
    fq_nmod_mpoly_init(got, ctx);
    dr_mpoly_to_fq_nmod_mpoly(got, actual, ctx);
    assert(fq_nmod_mpoly_equal(got, expected, ctx));
    fq_nmod_mpoly_clear(got, ctx);
}

static void check_field(ulong prime, slong degree)
{
    fq_nmod_ctx_t field;
    fq_nmod_ctx_init_ui(field, prime, degree, "a");
    fq_nmod_mpoly_ctx_t ctx;
    fq_nmod_mpoly_ctx_init(ctx, 3, ORD_LEX, field);
    flint_rand_t random;
    flint_rand_init(random);
    flint_rand_set_seed(random, prime, degree + 417);
    fq_nmod_mpoly_t a, b, reference, factor;
    fq_nmod_mpoly_init(a, ctx);
    fq_nmod_mpoly_init(b, ctx);
    fq_nmod_mpoly_init(reference, ctx);
    fq_nmod_mpoly_init(factor, ctx);
    fq_nmod_t coefficient;
    fq_nmod_init(coefficient, field);
    unified_mpoly_struct x = {0}, y = {0}, result = {0}, moved = {0};
    for (int trial = 0; trial < 12; trial++) {
        fq_nmod_mpoly_randtest_bound(a, random, 35, 4, ctx);
        fq_nmod_mpoly_randtest_bound(b, random, 27, 3, ctx);
        fq_nmod_mpoly_to_dr_mpoly(&x, a, 2, 1, ctx, field);
        fq_nmod_mpoly_to_dr_mpoly(&y, b, 2, 1, ctx, field);
        if (degree == 1)
            assert(x.field_id == FIELD_ID_NMOD);

        fq_nmod_mpoly_add(reference, a, b, ctx);
        dr_mpoly_add(&result, &x, &y);
        equal_native(&result, reference, ctx);
        dr_mpoly_copy(&moved, &x);
        dr_mpoly_add(&moved, &moved, &y);
        equal_native(&moved, reference, ctx);
        fq_nmod_mpoly_sub(reference, a, b, ctx);
        dr_mpoly_sub(&result, &x, &y);
        equal_native(&result, reference, ctx);
        fq_nmod_mpoly_mul(reference, a, b, ctx);
        dr_mpoly_mul(&result, &x, &y);
        equal_native(&result, reference, ctx);
        dr_mpoly_mul(&moved, &x, &y);
        equal_native(&moved, reference, ctx);
        dr_mpoly_mul(&x, &x, &y);
        equal_native(&x, reference, ctx);
        dr_mpoly_move(&result, &moved);
        assert(!moved.ring);
        equal_native(&result, reference, ctx);

        /* In-place scalar operations and parameter evaluation against FLINT. */
        fq_nmod_mpoly_to_dr_mpoly(&x, a, 2, 1, ctx, field);
        fq_nmod_randtest(coefficient, random, field);
        fq_nmod_mpoly_scalar_mul_fq_nmod(reference, a, coefficient, ctx);
        dr_mpoly_scalar_mul(&x, &x, coefficient);
        equal_native(&x, reference, ctx);
        fq_nmod_mpoly_neg(reference, reference, ctx);
        dr_mpoly_neg(&x, &x);
        equal_native(&x, reference, ctx);
        fq_nmod_t one, expected_value, actual_value;
        fq_nmod_init(one, field);
        fq_nmod_init(expected_value, field);
        fq_nmod_init(actual_value, field);
        fq_nmod_one(one, field);
        fq_nmod_struct *values[] = {one, one, coefficient};
        fq_nmod_mpoly_evaluate_all_fq_nmod(expected_value, reference, values, ctx);
        evaluate_dr_mpoly_at_params(actual_value, &x, (const fq_nmod_t *)&coefficient);
        assert(fq_nmod_equal(expected_value, actual_value, field));
        fq_nmod_clear(actual_value, field);
        fq_nmod_clear(expected_value, field);
        fq_nmod_clear(one, field);

        /* Difference division must be exact in every characteristic/backend. */
        fq_nmod_mpoly_gen(factor, 0, ctx);
        fq_nmod_mpoly_gen(reference, 1, ctx);
        fq_nmod_mpoly_sub(factor, factor, reference, ctx);
        fq_nmod_mpoly_mul(reference, a, factor, ctx);
        fq_nmod_mpoly_to_dr_mpoly(&x, reference, 2, 1, ctx, field);
        dr_mpoly_divide_difference(&x, &x, 0, 2, 1);
        equal_native(&x, a, ctx);

        /* Unsorted duplicate insertion, including exact cancellation. */
        dr_mpoly_init(&result, 2, 1, field);
        for (slong i = a->length; i-- > 0;) {
            slong exp[3];
            fq_nmod_mpoly_get_term_exp_si(exp, a, i, ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(coefficient, a, i, ctx);
            dr_mpoly_add_term_fast(&result, exp, exp + 2, coefficient);
            dr_mpoly_add_term_fast(&result, exp, exp + 2, coefficient);
            fq_nmod_neg(coefficient, coefficient, field);
            dr_mpoly_add_term_fast(&result, exp, exp + 2, coefficient);
        }
        dr_mpoly_normalize(&result);
        equal_native(&result, a, ctx);
    }
    /* Explicit field-equation reduction, including the parameter coordinate. */
    if (prime < 300 && degree <= 2) {
        ulong q = 1;
        for (slong i = 0; i < degree; i++)
            q *= prime;
        slong exp[3] = {q + 1, 0, q};
        dr_mpoly_init(&x, 2, 1, field);
        fq_nmod_one(coefficient, field);
        dr_mpoly_add_term_fast(&x, exp, exp + 2, coefficient);
        dr_mpoly_reduce_field_equation(&x);
        fq_nmod_mpoly_zero(reference, ctx);
        ulong reduced[3] = {1 + q % (q - 1), 0, 1};
        fq_nmod_mpoly_set_coeff_fq_nmod_ui(reference, coefficient, reduced, ctx);
        equal_native(&x, reference, ctx);
        dr_mpoly_set_field_equation_reduction(1);
        dr_mpoly_mul(&result, &x, &x);
        dr_mpoly_set_field_equation_reduction(0);
        for (slong i = 0; i < 3; i++)
            reduced[i] = reduced[i] ? 1 + (2 * reduced[i] - 1) % (q - 1) : 0;
        fq_nmod_mpoly_zero(reference, ctx);
        fq_nmod_mpoly_set_coeff_fq_nmod_ui(reference, coefficient, reduced, ctx);
        equal_native(&result, reference, ctx);
    }
    dr_mpoly_clear(&moved);
    dr_mpoly_clear(&result);
    dr_mpoly_clear(&y);
    dr_mpoly_clear(&x);
    fq_nmod_clear(coefficient, field);
    fq_nmod_mpoly_clear(factor, ctx);
    fq_nmod_mpoly_clear(reference, ctx);
    fq_nmod_mpoly_clear(b, ctx);
    fq_nmod_mpoly_clear(a, ctx);
    flint_rand_clear(random);
    fq_nmod_mpoly_ctx_clear(ctx);
    fq_nmod_ctx_clear(field);
}

static void check_zech(void)
{
    fq_nmod_ctx_t fq;
    fq_nmod_ctx_init_ui(fq, 7, 2, "a");
    field_ctx_t field = {0};
    field.field_id = FIELD_ID_FQ_ZECH;
    field.ctx.zech_ctx = flint_malloc(sizeof(fq_zech_ctx_struct));
    fq_zech_ctx_init_modulus(field.ctx.zech_ctx, fq_nmod_ctx_modulus(fq), "a");
    unified_mpoly_ctx_t ctx = unified_mpoly_ctx_init(3, ORD_LEX, &field);
    unified_mpoly_t p = unified_mpoly_init(ctx);
    p->ctx = fq;
    p->nvars = 2;
    p->npars = 1;
    field_elem_u c = {0}, got = {0};
    fq_nmod_t coefficient;
    fq_nmod_init(coefficient, fq);
    fq_nmod_gen(coefficient, fq);
    fq_zech_set_fq_nmod(&c.fq_zech, coefficient, field.ctx.zech_ctx);
    ulong exp[] = {2, 1, 3};
    unified_mpoly_set_coeff_ui(p, &c, exp);
    unified_mpoly_get_coeff_ui(&got, p, exp);
    assert(fq_zech_equal(&got.fq_zech, &c.fq_zech, field.ctx.zech_ctx));
    /* Setting a coefficient replaces it, including replacement by zero. */
    fq_zech_one(&c.fq_zech, field.ctx.zech_ctx);
    unified_mpoly_set_coeff_ui(p, &c, exp);
    unified_mpoly_get_coeff_ui(&got, p, exp);
    assert(fq_zech_is_one(&got.fq_zech, field.ctx.zech_ctx));
    DR_MPOLY_TERM(term, p, 0);
    assert(term.var_exp[0] == 2 && term.par_exp[0] == 3);
    assert(fq_nmod_is_one(term.coeff, fq));
    assert(unified_mpoly_mul(p, p, p));
    for (slong i = 0; i < 3; i++)
        exp[i] *= 2;
    unified_mpoly_get_coeff_ui(&got, p, exp);
    assert(fq_zech_is_one(&got.fq_zech, field.ctx.zech_ctx));
    fq_zech_zero(&c.fq_zech, field.ctx.zech_ctx);
    unified_mpoly_set_coeff_ui(p, &c, exp);
    assert(unified_mpoly_is_zero(p));
    unified_mpoly_clear(p);
    unified_mpoly_ctx_clear(ctx);
    field_ctx_clear(&field);
    fq_nmod_clear(coefficient, fq);
    fq_nmod_ctx_clear(fq);
}

int main(void)
{
    check_field(2, 1);
    check_field(257, 1);
    check_field(UWORD(18446744073709551557), 1);
    check_field(7, 2);
    check_field(2, 4);
    check_field(2, 8);
    check_field(65537, 2);
    check_zech();
    puts("Native polynomial storage: arithmetic, aliasing, ownership, division and reduction PASS");
    flint_cleanup_master();
    return 0;
}

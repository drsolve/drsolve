/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include <flint/fmpz_mod_mat.h>
#include "unified_dixon.h"
#include "unified_mpoly_resultant.h"
#include "unified_mat_poly.h"
#include "dixon_interface_flint.h"

static void check_limb_prime(void)
{
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_str(p, "18446744073709551557", 10);
    field_ctx_t f;
    field_ctx_init_prime(&f, p);
    assert(f.field_id == FIELD_ID_NMOD);
    assert(f.ctx.nmod_ctx.n == fmpz_get_ui(p));
    field_elem_u a, b;
    a.nmod = f.ctx.nmod_ctx.n - 2;
    field_mul(&b, &a, &a, f.field_id, field_ctx_data(&f));
    assert(b.nmod == 4);
    field_ctx_clear(&f);
    fmpz_clear(p);
}

static void check_field(const char *prime)
{
    fmpz_t p, x;
    fmpz_init(p);
    fmpz_init(x);
    assert(!fmpz_set_str(p, prime, 10));
    field_ctx_t f;
    field_ctx_init_prime(&f, p);
    assert(f.field_id == FIELD_ID_FMPZ_MOD);
    const void *ctx = field_ctx_data(&f);
    field_elem_u a, b, c;
    field_init_elem(&a, f.field_id, ctx);
    field_init_elem(&b, f.field_id, ctx);
    field_init_elem(&c, f.field_id, ctx);
    fmpz_sub_ui(&a.big, p, 2);
    field_set_elem(&b, &a, f.field_id, ctx);
    field_mul(&c, &a, &b, f.field_id, ctx);
    assert(fmpz_equal_ui(&c.big, 4));
    field_inv(&c, &a, f.field_id, ctx);
    field_mul(&c, &c, &a, f.field_id, ctx);
    assert(field_is_one(&c, f.field_id, ctx));
    field_neg(&c, &a, f.field_id, ctx);
    assert(fmpz_equal_ui(&c.big, 2));
    field_add(&c, &a, &c, f.field_id, ctx);
    assert(field_is_zero(&c, f.field_id, ctx));
    unified_poly_t u, v;
    unified_poly_init(u, &f);
    unified_poly_init(v, &f);
    unified_poly_set_coeff(u, 0, &a);
    unified_poly_set_coeff(u, 40, &b);
    unified_poly_mul(v, u, u);
    unified_poly_mul(u, u, u);
    for (slong i = 0; i <= 80; i++) {
        unified_poly_get_coeff(&b, u, i);
        unified_poly_get_coeff(&c, v, i);
        assert(field_equal(&b, &c, f.field_id, ctx));
    }
    assert(fmpz_equal_ui(&u->coeffs[40].big, 8));
    unified_poly_clear(u);
    unified_poly_clear(v);
    unified_mpoly_ctx_t ring = unified_mpoly_ctx_init(2, ORD_LEX, &f);
    unified_mpoly_t q = unified_mpoly_init(ring), r = unified_mpoly_init(ring),
                    product = unified_mpoly_init(ring);
    const char *names[] = {"x", "y"};
    assert(!fmpz_mod_mpoly_set_str_pretty(GET_BIG_POLY(q), "x^2+y+1", names, GET_BIG_CTX(ring)));
    assert(!fmpz_mod_mpoly_set_str_pretty(GET_BIG_POLY(r), "x-y", names, GET_BIG_CTX(ring)));
    unified_mpoly_mul(product, q, r);
    assert(unified_mpoly_divides(product, product, r));
    assert(unified_mpoly_equal(product, q));
    assert(unified_mpoly_resultant(product, q, r, 0, ring));
    assert(!fmpz_mod_mpoly_set_str_pretty(GET_BIG_POLY(q), "y^2+y+1", names, GET_BIG_CTX(ring)));
    assert(unified_mpoly_equal(product, q));
    flint_rand_t state;
    flint_rand_init(state);
    unified_mpoly_randtest(product, state, 20, 5);
    flint_rand_clear(state);
    unified_mpoly_clear(q);
    unified_mpoly_clear(r);
    unified_mpoly_clear(product);
    unified_mpoly_ctx_clear(ring);
    field_clear_elem(&a, f.field_id, ctx);
    field_clear_elem(&b, f.field_id, ctx);
    field_clear_elem(&c, f.field_id, ctx);
    ensure_workspace_initialized(&f);
    field_ctx_clear(&f);
    fmpz_clear(p);
    fmpz_clear(x);
}

/* A dense pencil has degree at most four: five independent scalar determinants
 * uniquely check its HNF determinant, including non-small fmpz coefficients. */
static void check_big_hnf(void)
{
    fmpz_t p, expected, actual, x;
    fmpz_init(p);
    fmpz_init(expected);
    fmpz_init(actual);
    fmpz_init(x);
    fmpz_set_str(p, "18446744073709551629", 10);
    field_ctx_t f;
    field_ctx_init_prime(&f, p);
    unified_poly_mat_t polys;
    unified_poly_mat_init(polys, 4, 4, &f);
    field_elem_u c;
    field_init_elem(&c, f.field_id, field_ctx_data(&f));
    for (slong i = 0; i < 4; i++)
        for (slong j = 0; j < 4; j++) {
            fmpz_sub_ui(&c.big, p, (i + 1) * (j + 1) + (i == j ? i + 1 : 0));
            unified_poly_set_coeff(unified_poly_mat_entry(polys, i, j), 0, &c);
            if (i == j) {
                fmpz_one(&c.big);
                unified_poly_set_coeff(unified_poly_mat_entry(polys, i, j), 1, &c);
            }
        }
    unified_mat_poly_t mat;
    unified_mat_poly_init(mat, 4, 4, &f);
    unified_mat_poly_from_unified(mat, polys);
    unified_poly_t det;
    unified_poly_init(det, &f);
    assert(unified_mat_poly_det_hnf(det, mat));
    assert(det->length <= 5);
    fmpz_mod_mat_t numeric;
    fmpz_mod_mat_init(numeric, 4, 4, f.ctx.big_ctx);
    for (slong k = 0; k < 5; k++) {
        for (slong i = 0; i < 4; i++)
            for (slong j = 0; j < 4; j++) {
                fmpz_set_si(x, -((i + 1) * (j + 1) + (i == j ? i + 1 : 0)) + (i == j ? k : 0));
                fmpz_mod_set_fmpz(fmpz_mod_mat_entry(numeric, i, j), x, f.ctx.big_ctx);
            }
        fmpz_mod_mat_det(expected, numeric, f.ctx.big_ctx);
        fmpz_zero(actual);
        fmpz_set_si(x, k);
        for (slong i = det->length; i-- > 0;) {
            fmpz_mod_mul(actual, actual, x, f.ctx.big_ctx);
            fmpz_mod_add(actual, actual, &det->coeffs[i].big, f.ctx.big_ctx);
        }
        assert(fmpz_equal(actual, expected));
    }
    fmpz_mod_mat_clear(numeric, f.ctx.big_ctx);
    unified_poly_clear(det);
    unified_mat_poly_clear(mat);
    unified_poly_mat_clear(polys);
    field_clear_elem(&c, f.field_id, field_ctx_data(&f));
    field_ctx_clear(&f);
    fmpz_clear(p);
    fmpz_clear(expected);
    fmpz_clear(actual);
    fmpz_clear(x);
}

/* Dense total-degree-five systems should not acquire an extra parameter
 * factor merely from the choice of the maximal-rank Dixon minor. */
static void check_quintics(void)
{
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_str(p, "18446744073709551629", 10);
    field_ctx_t f;
    field_ctx_init_prime(&f, p);
    unified_mpoly_ctx_t in = unified_mpoly_ctx_init(3, ORD_LEX, &f);
    unified_mpoly_ctx_t out = unified_mpoly_ctx_init(1, ORD_LEX, &f);
    unified_mpoly_t polys[3], result = unified_mpoly_init(out);
    flint_rand_t rng;
    flint_rand_init(rng);
    flint_rand_set_seed(rng, 12345, 67890);
    field_elem_u c;
    field_init_elem(&c, f.field_id, field_ctx_data(&f));
    for (slong i = 0; i < 3; i++) {
        polys[i] = unified_mpoly_init(in);
        for (ulong x = 0; x <= 5; x++)
            for (ulong y = 0; y <= 5 - x; y++)
                for (ulong z = 0; z <= 5 - x - y; z++) {
                    ulong e[3] = {x, y, z};
                    fmpz_set_ui(&c.big, 1 + n_randint(rng, 100000));
                    unified_mpoly_set_coeff_ui(polys[i], &c, e);
                }
    }
    assert(unified_dixon_resultant(result, polys, 2));
    assert(unified_mpoly_total_degree_si(result) == 125);
    unified_mpoly_t reference = unified_mpoly_init(out);
    int saved_method = dixon_global_method_step4;
    dixon_global_method_step4 = DET_METHOD_INTERPOLATION;
    assert(unified_dixon_resultant(reference, polys, 2));
    dixon_global_method_step4 = saved_method;
    assert(unified_mpoly_equal(result, reference));
    unified_mpoly_clear(reference);
    for (slong i = 0; i < 3; i++)
        unified_mpoly_clear(polys[i]);
    unified_mpoly_clear(result);
    field_clear_elem(&c, f.field_id, field_ctx_data(&f));
    unified_mpoly_ctx_clear(in);
    unified_mpoly_ctx_clear(out);
    field_ctx_clear(&f);
    flint_rand_clear(rng);
    fmpz_clear(p);
}

static void check_result(const char *input, const char *elim, const char **vars, slong nvars,
                         const char *expected)
{
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_str(p, "18446744073709551629", 10);
    char *answer = dixon_str_large_prime(input, elim, p);
    assert(answer);
    fmpz_mod_mpoly_ctx_t ctx;
    fmpz_mod_mpoly_ctx_init(ctx, nvars, ORD_LEX, p);
    fmpz_mod_mpoly_t a, b;
    fmpz_mod_mpoly_init(a, ctx);
    fmpz_mod_mpoly_init(b, ctx);
    assert(!fmpz_mod_mpoly_set_str_pretty(a, answer, vars, ctx));
    assert(!fmpz_mod_mpoly_set_str_pretty(b, expected, vars, ctx));
    if (!fmpz_mod_mpoly_is_zero(a, ctx))
        fmpz_mod_mpoly_make_monic(a, a, ctx);
    if (!fmpz_mod_mpoly_is_zero(b, ctx))
        fmpz_mod_mpoly_make_monic(b, b, ctx);
    if (!fmpz_mod_mpoly_equal(a, b, ctx)) {
        printf("Input: %s\nActual: %s\nExpected: %s\n", input, answer, expected);
        abort();
    }
    fmpz_mod_mpoly_clear(a, ctx);
    fmpz_mod_mpoly_clear(b, ctx);
    fmpz_mod_mpoly_ctx_clear(ctx);
    fmpz_clear(p);
    free(answer);
}
int main(void)
{
    g_dixon_verbose_level = 0;
    check_limb_prime();
    check_field("18446744073709551629");
    check_field("170141183460469231731687303715884105727");
    check_field("57896044618658097711785492504343953926634992332820282019728792003956564819949");
    check_big_hnf();
    check_quintics();
    const char *y[] = {"y"}, *z[] = {"z"}, *yz[] = {"y", "z"};
    check_result("x+y,x-y+1", "x", y, 1, "2*y-1");
    check_result("x^2+y+1,x-y", "x", y, 1, "y^2+y+1");
    check_result("x+y+z,x*y+y*z+z*x,x*y*z+1", "x,y", z, 1, "(z^3+1)^2");
    check_result("x+y*z,x-y+1", "x", yz, 2, "y*z+y-1");
    check_result("x+y,2*x+2*y", "x", y, 1, "0");
    check_result("18446744073709551629*x^3+x+y,x-y", "x", y, 1, "y");
    check_result("x/2+y,x-y", "x", y, 1, "y");
    check_result("x+1,x+2", "x", NULL, 0, "1");
    check_result("x+y,x+z,y+z-1", "x,y", z, 1, "2*z-1");
    puts("fmpz unified backend and direct Dixon tests passed");
    flint_cleanup_master();
    return 0;
}

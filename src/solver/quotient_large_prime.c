/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "quotient_solver.h"
#include "quotient_closure.h"
#include <flint/fmpz_mod_mpoly.h>
#include <flint/fmpz_mod_mat.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_mod_poly_factor.h>

typedef struct {
    slong num_variables, num_solution_sets;
    const fmpz_mod_ctx_struct *ctx;
    fmpz_mat_t points;
} prime_points_t;

static slong prime_input_length(const quotient_input_t *input, slong i)
{
    return ((const fmpz_mod_mpoly_struct *)input->polys)[i].length;
}
static void prime_input_term(field_elem_u *c, ulong *exp, const quotient_input_t *input, slong i,
                             slong t, const quotient_field_t *e)
{
    (void)e;
    const fmpz_mod_mpoly_struct *f = (const fmpz_mod_mpoly_struct *)input->polys + i;
    fmpz_mod_mpoly_get_term_coeff_fmpz(&c->big, f, t, input->ctx);
    fmpz_mod_mpoly_get_term_exp_ui(exp, f, t, input->ctx);
}

static int prime_points(prime_points_t *sols, fmpz_mod_mat_struct *T, const fmpz_mod_mat_t S,
                        slong variable, fmpz *point, fmpz_mod_mpoly_struct *fs, slong m,
                        const fmpz_mod_mpoly_ctx_t ctx)
{
    const fmpz_mod_ctx_struct *fq = sols->ctx;
    slong n = sols->num_variables, dim = S->r, k = S->c;
    if (variable == n) {
        if (k != 1 || sols->num_solution_sets >= T[0].r)
            return 0;
        fmpz *values[n];
        for (slong j = 0; j < n; j++)
            values[j] = point + j;
        fmpz_t ev;
        fmpz_init(ev);
        int ok = 1;
        for (slong i = 0; i < m && ok; i++) {
            fmpz_mod_mpoly_evaluate_all_fmpz(ev, fs + i, values, ctx);
            ok = fmpz_is_zero(ev);
        }
        fmpz_clear(ev);
        if (!ok)
            return 0;
        slong row = sols->num_solution_sets++;
        for (slong j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(sols->points, row, j), point + j);
        return 1;
    }
    fmpz_mod_mat_t ST, SB, inv, TS, rows, A;
    fmpz_mod_mat_init(ST, k, dim, fq);
    fmpz_mod_mat_transpose(ST, S, fq);
    slong rank = fmpz_mod_mat_rref(ST, ST, fq);
    if (rank != k) {
        fmpz_mod_mat_clear(ST, fq);
        return 0;
    }
    slong *rr = flint_malloc(k * sizeof(slong));
    for (slong i = 0; i < k; i++) {
        slong j = 0;
        while (j < dim && fmpz_is_zero(fmpz_mod_mat_entry(ST, i, j)))
            j++;
        rr[i] = j;
    }
    fmpz_mod_mat_clear(ST, fq);
    fmpz_mod_mat_init(SB, k, k, fq);
    fmpz_mod_mat_init(inv, k, k, fq);
    for (slong i = 0; i < k; i++)
        for (slong j = 0; j < k; j++)
            fmpz_set(fmpz_mod_mat_entry(SB, i, j), fmpz_mod_mat_entry(S, rr[i], j));
    int ok = fmpz_mod_mat_inv(inv, SB, fq);
    fmpz_mod_mat_clear(SB, fq);
    if (!ok) {
        flint_free(rr);
        fmpz_mod_mat_clear(inv, fq);
        return 0;
    }
    fmpz_mod_mat_init(TS, dim, k, fq);
    fmpz_mod_mat_mul(TS, T + variable, S, fq);
    fmpz_mod_mat_init(rows, k, k, fq);
    for (slong i = 0; i < k; i++)
        for (slong j = 0; j < k; j++)
            fmpz_set(fmpz_mod_mat_entry(rows, i, j), fmpz_mod_mat_entry(TS, rr[i], j));
    flint_free(rr);
    fmpz_mod_mat_clear(TS, fq);
    fmpz_mod_mat_init(A, k, k, fq);
    fmpz_mod_mat_mul(A, inv, rows, fq);
    fmpz_mod_mat_clear(inv, fq);
    fmpz_mod_mat_clear(rows, fq);
    fmpz_mod_poly_t mu;
    fmpz_mod_poly_init(mu, fq);
    fmpz_mod_mat_minpoly(mu, A, fq);
    fmpz_mod_poly_factor_t roots;
    fmpz_mod_poly_factor_init(roots, fq);
    fmpz_mod_poly_roots(roots, mu, 0, fq);
    fmpz_mod_poly_clear(mu, fq);
    fmpz_t a;
    fmpz_init(a);
    for (slong r = 0; r < roots->num && ok; r++) {
        fmpz_mod_poly_get_coeff_fmpz(a, roots->poly + r, 0, fq);
        fmpz_mod_neg(a, a, fq);
        fmpz_set(point + variable, a);
        fmpz_mod_mat_t shifted, K, basis, next;
        fmpz_mod_mat_init_set(shifted, A, fq);
        for (slong j = 0; j < k; j++)
            fmpz_mod_sub(fmpz_mod_mat_entry(shifted, j, j), fmpz_mod_mat_entry(shifted, j, j), a,
                         fq);
        fmpz_mod_mat_init(K, k, k, fq);
        slong nullity = fmpz_mod_mat_nullspace(K, shifted, fq);
        fmpz_mod_mat_clear(shifted, fq);
        if (!nullity) {
            fmpz_mod_mat_clear(K, fq);
            ok = 0;
            break;
        }
        fmpz_mod_mat_window_init(basis, K, 0, 0, k, nullity, fq);
        fmpz_mod_mat_init(next, dim, nullity, fq);
        fmpz_mod_mat_mul(next, S, basis, fq);
        fmpz_mod_mat_window_clear(basis, fq);
        fmpz_mod_mat_clear(K, fq);
        ok = prime_points(sols, T, next, variable + 1, point, fs, m, ctx);
        fmpz_mod_mat_clear(next, fq);
    }
    fmpz_clear(a);
    fmpz_mod_poly_factor_clear(roots, fq);
    fmpz_mod_mat_clear(A, fq);
    return ok;
}

static int prime_degree_fits(slong n, slong d, slong memory_mb, const quotient_field_t *e,
                             slong *columns)
{
    slong N = monomial_count_bounded(n, d, WORD_MAX / (slong)sizeof(field_elem_u));
    if (!N)
        return -1;
    /* Echelon relations, normal forms, operators, certification products,
     * and the recursive FLINT eigenspaces. Account for multiprecision limb allocations. */
    long double scalar =
        sizeof(field_elem_u) +
        3.0L * ((fmpz_bits(fmpz_mod_ctx_modulus(e->arithmetic)) + FLINT_BITS - 1) / FLINT_BITS) *
            sizeof(ulong);
    long double bytes = (16.0L + 4.0L * n) * N * N * scalar;
    if (bytes > WORD_MAX)
        return -1;
    if (memory_mb && bytes > (long double)memory_mb * 1024 * 1024)
        return 0;
    *columns = N;
    return 1;
}

int solve_by_quotient_closure_large_prime(char **polys, slong count, const char **names, slong n,
                                          const fmpz_t prime, slong max_degree, slong memory_mb,
                                          fmpz_mat_t points, char **summary, char **error)
{
    if (max_degree < 0 || memory_mb < 0 || n <= 0 || count < n || fmpz_bits(prime) <= FLINT_BITS) {
        *error = strdup("Invalid large-prime quotient closure arguments");
        return 0;
    }
    quotient_field_t e;
    memset(&e, 0, sizeof(e));
    field_ctx_init_prime(&e.field, prime);
    e.arithmetic = field_ctx_data(&e.field);
    const fmpz_mod_ctx_struct *fq = e.arithmetic;
    fmpz_mod_mpoly_ctx_t ctx;
    fmpz_mod_mpoly_ctx_init(ctx, n, ORD_DEGLEX, prime);
    fmpz_mod_mpoly_struct *fs = flint_malloc(count * sizeof(*fs));
    quotient_input_t input = {fs, ctx, prime_input_length, prime_input_term};
    prime_points_t storage = {0};
    prime_points_t *sols = &storage;
    sols->ctx = fq;
    sols->num_variables = n;
    fmpz_mat_init(sols->points, 0, n);
    slong *degrees = flint_malloc(count * sizeof(slong)), start = 1, initialized = 0;
    slong used = 0, columns = 0, dim = 0;
    int success = 0, certified = 0;
    quotient_closure_t previous;
    int have_previous = 0;
    unified_mat_struct *operators = NULL;
    const char *reason = "degree limit reached";
    for (slong i = 0; i < count; i++) {
        fmpz_mod_mpoly_init(fs + i, ctx);
        initialized++;
        if (fmpz_mod_mpoly_set_str_pretty(fs + i, polys[i], names, ctx) ||
            !fmpz_mod_mpoly_total_degree_fits_si(fs + i, ctx)) {
            *error = strdup("Could not parse polynomial for quotient closure (or degree overflow)");
            goto done;
        }
        degrees[i] = fs[i].length ? fmpz_mod_mpoly_total_degree_si(fs + i, ctx) : -1;
        start = FLINT_MAX(start, degrees[i]);
        if (degrees[i] == 0)
            certified = 1;
    }
    for (slong d = start; !certified && (!max_degree || d <= max_degree); d++) {
        int fits = prime_degree_fits(n, d, memory_mb, &e, &columns);
        if (fits != 1) {
            reason =
                fits < 0 ? "matrix dimensions exceed addressable range" : "memory budget reached";
            break;
        }
        used = d;
        quotient_closure_t cl;
        quotient_closure_init(&cl, n, d);
        if (have_previous) {
            quotient_closure_seed(&cl, &previous, &e);
            quotient_closure_clear(&previous, &e);
            have_previous = 0;
        } else
            insert_generators(&cl, &input, count, degrees, &e);
        extend_relations(&cl, &e);
        dim = cl.basis.count - cl.rank;
        if (!dim)
            certified = 1;
        else {
            operators = build_operators(&cl, &input, count, &e);
            certified = operators != NULL;
        }
        if (certified)
            quotient_closure_clear(&cl, &e);
        else {
            previous = cl;
            have_previous = 1;
        }
        if (d == WORD_MAX) {
            reason = "degree exceeds integer range";
            break;
        }
    }
    if (!certified) {
        char message[384];
        snprintf(message, sizeof(message),
                 "Quotient closure incomplete: %s (degree limit %ld, memory budget %ld MiB; 0 = "
                 "unlimited; last attempted degree %ld). No conclusion about existence or "
                 "dimension of solutions.",
                 reason, max_degree, memory_mb, used);
        *error = strdup(message);
        goto done;
    }
    if (dim) {
        fmpz_mod_mat_struct *T = flint_malloc(n * sizeof(*T));
        for (slong v = 0; v < n; v++) {
            fmpz_mod_mat_init(T + v, dim, dim, fq);
            for (slong i = 0; i < dim; i++)
                for (slong j = 0; j < dim; j++)
                    fmpz_set(fmpz_mod_mat_entry(T + v, i, j), &ENTRY(operators + v, i, j)->big);
        }
        fmpz_mat_clear(sols->points);
        fmpz_mat_init(sols->points, dim, n);
        fmpz *point = flint_malloc(n * sizeof(*point));
        for (slong i = 0; i < n; i++)
            fmpz_init(point + i);
        fmpz_mod_mat_t S;
        fmpz_mod_mat_init(S, dim, dim, fq);
        fmpz_mod_mat_one(S, fq);
        success = prime_points(sols, T, S, 0, point, fs, count, ctx);
        fmpz_mod_mat_clear(S, fq);
        for (slong i = 0; i < n; i++) {
            fmpz_clear(point + i);
            fmpz_mod_mat_clear(T + i, fq);
        }
        flint_free(point);
        flint_free(T);
        if (!success) {
            *error = strdup(
                "Quotient multiplication representation failed during rational-point extraction");
            goto done;
        }
    }
    success = 1;
    fmpz_mat_clear(points);
    fmpz_mat_init(points, sols->num_solution_sets, n);
    for (slong i = 0; i < sols->num_solution_sets; i++)
        for (slong j = 0; j < n; j++)
            fmpz_set(fmpz_mat_entry(points, i, j), fmpz_mat_entry(sols->points, i, j));
    char report[384];
    snprintf(report, sizeof(report),
             "Quotient algebra over large prime: certified dimension %ld, degree %ld, %ld "
             "monomials; %ld rational point(s)",
             dim, used, columns, sols->num_solution_sets);
    *summary = strdup(report);
done:
    if (have_previous)
        quotient_closure_clear(&previous, &e);
    if (operators) {
        for (slong i = 0; i < n; i++)
            unified_mat_clear(operators + i);
        flint_free(operators);
    }
    for (slong i = 0; i < initialized; i++)
        fmpz_mod_mpoly_clear(fs + i, ctx);
    flint_free(fs);
    flint_free(degrees);
    fmpz_mat_clear(sols->points);
    fmpz_mod_mpoly_ctx_clear(ctx);
    field_ctx_clear(&e.field);
    return success;
}

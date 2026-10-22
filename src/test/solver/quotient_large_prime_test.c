/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "quotient_solver.h"
#include <flint/fmpz_mod_mpoly.h>
#include <assert.h>

static const char *names[] = {"x", "y"};
static void check(char **polys, slong m, const fmpz_t prime, slong expected, const char *dimension,
                  slong degree, slong memory)
{
    fmpz_mat_t points;
    fmpz_mat_init(points, 1, 1);
    fmpz_set_ui(fmpz_mat_entry(points, 0, 0), 42);
    char *summary = NULL, *error = NULL;
    int ok = solve_by_quotient_closure_large_prime(polys, m, names, 2, prime, degree, memory,
                                                   points, &summary, &error);
    if (expected < 0) {
        assert(!ok && error && strstr(error, "incomplete"));
        assert(points->r == 1 && points->c == 1);
        assert(fmpz_equal_ui(fmpz_mat_entry(points, 0, 0), 42));
    } else {
        if (!ok)
            fprintf(stderr, "%s\n", error);
        assert(ok && !error && points->r == expected && points->c == 2);
        assert(strstr(summary, dimension));
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_ctx_init(ctx, 2, ORD_DEGLEX, prime);
        fmpz_mod_mpoly_t f;
        fmpz_mod_mpoly_init(f, ctx);
        fmpz_t value;
        fmpz_init(value);
        for (slong i = 0; i < expected; i++) {
            fmpz *values[] = {fmpz_mat_entry(points, i, 0), fmpz_mat_entry(points, i, 1)};
            for (slong j = 0; j < m; j++) {
                assert(!fmpz_mod_mpoly_set_str_pretty(f, polys[j], names, ctx));
                fmpz_mod_mpoly_evaluate_all_fmpz(value, f, values, ctx);
                assert(fmpz_is_zero(value));
            }
            for (slong j = 0; j < i; j++)
                assert(!fmpz_equal(values[0], fmpz_mat_entry(points, j, 0)) ||
                       !fmpz_equal(values[1], fmpz_mat_entry(points, j, 1)));
        }
        fmpz_clear(value);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }
    free(summary);
    free(error);
    fmpz_mat_clear(points);
}
int main(void)
{
    fmpz_t prime;
    fmpz_init(prime);
    const char *primes[] = {"18446744073709551629", "170141183460469231731687303715884105727"};
    for (int i = 0; i < 2; i++) {
        fmpz_set_str(prime, primes[i], 10);
        char *grid[] = {"x^2-x", "y^2-y", "(x^2-x)*y"};
        check(grid, 3, prime, 4, "dimension 4,", 8, 0);
        char *fat[] = {"(x-1)^2", "(x-1)*(y-2)", "(y-2)^2"};
        check(fat, 3, prime, 1, "dimension 3,", 8, 0);
        char *none[] = {"x", "y", "x+y-1"};
        check(none, 3, prime, 0, "dimension 0,", 8, 0);
        char *positive[] = {"x*y", "x^2*y"};
        check(positive, 2, prime, -1, NULL, 4, 0);
        check(grid, 3, prime, -1, NULL, 1, 0);
        char *large[] = {"x^20-1", "y^20-1"};
        check(large, 2, prime, -1, NULL, 0, 1);
    }
    /* -1 is a nonsquare for this prime: nonzero quotient, no rational points. */
    char *extension_only[] = {"x^2+1", "y"};
    check(extension_only, 2, prime, 0, "dimension 2,", 8, 0);
    fmpz_clear(prime);
    puts("Large-prime quotient closure checks passed.");
    return 0;
}

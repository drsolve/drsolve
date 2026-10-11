/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include "polynomial_system_solver.h"

/* x+y=0, y^2=1: x=0 is an impossible branch; x=1,-1 both lift.
 * Supplying base roots directly tests certified empty fibers as well as points. */
static void check(polynomial_solver_method_t method, slong limit,
                  int empty_only, const fq_nmod_ctx_t ctx)
{
    char *polys[] = {"x+y", "y^2-1"};
    variable_info_t vars[] = {{"x", 1, 0}, {"y", 2, 1}};
    fq_nmod_t roots[3];
    for (slong i = 0; i < 3; i++) fq_nmod_init(roots[i], ctx);
    fq_nmod_zero(roots[0], ctx);
    fq_nmod_one(roots[1], ctx);
    fq_nmod_neg(roots[2], roots[1], ctx);
    polynomial_solutions_t sols;
    polynomial_solutions_init(&sols, 2, ctx);
    sols.variable_names[0] = strdup("x");
    sols.variable_names[1] = strdup("y");
    polynomial_solver_set_method(method);
    polynomial_solver_set_quotient_limits(limit, 0);
    assert(solve_by_back_substitution_recursive_enhanced(
        polys, 2, vars, 2, roots, empty_only ? 1 : 3, &sols));
    assert(!sols.error_message);
    assert(sols.num_solution_sets == (empty_only ? 0 : 2));
    assert(sols.has_no_solutions == empty_only);
    int quotient_steps = 0, fallback_steps = 0, dixon_steps = 0;
    for (slong i = 0; i < sols.num_resultant_steps; i++) {
        quotient_steps += strstr(sols.resultant_steps[i], "-> quotient algebra") != NULL;
        fallback_steps += strstr(sols.resultant_steps[i], "Quotient failed;") != NULL;
        dixon_steps += strstr(sols.resultant_steps[i], "-> solve reduced eq") != NULL;
    }
    int hybrid = method == POLYNOMIAL_SOLVER_AUTO;
    assert(quotient_steps == (hybrid ? (empty_only ? 1 : 3) : 0));
    if (hybrid && !limit) {
        assert(fallback_steps == 0);
        assert(dixon_steps == 0); /* Includes the certified empty branch. */
    } else {
        assert(dixon_steps > 0);
        assert((fallback_steps > 0) == hybrid);
    }
    fq_nmod_t sum, square;
    fq_nmod_init(sum, ctx); fq_nmod_init(square, ctx);
    for (slong i = 0; i < sols.num_solution_sets; i++) {
        assert(sols.solutions_per_var[2*i] == 1);
        assert(sols.solutions_per_var[2*i+1] == 1);
        fq_nmod_add(sum, sols.solution_sets[i][0][0], sols.solution_sets[i][1][0], ctx);
        fq_nmod_sqr(square, sols.solution_sets[i][1][0], ctx);
        assert(fq_nmod_is_zero(sum, ctx));
        assert(fq_nmod_is_one(square, ctx));
    }
    if (!empty_only)
        assert(!fq_nmod_equal(sols.solution_sets[0][0][0], sols.solution_sets[1][0][0], ctx));
    fq_nmod_clear(sum, ctx); fq_nmod_clear(square, ctx);
    polynomial_solutions_clear(&sols);
    for (slong i = 0; i < 3; i++) fq_nmod_clear(roots[i], ctx);
}

int main(void)
{
    fmpz_t prime;
    fmpz_init_set_ui(prime, 101);
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, prime, 1, "t");
    for (int empty = 0; empty <= 1; empty++) {
        check(POLYNOMIAL_SOLVER_AUTO, 0, empty, ctx);
        check(POLYNOMIAL_SOLVER_AUTO, 1, empty, ctx);
        check(POLYNOMIAL_SOLVER_DIXON, 0, empty, ctx);
    }
    polynomial_solver_set_method(POLYNOMIAL_SOLVER_AUTO);
    polynomial_solver_set_quotient_limits(0, 0);
    fq_nmod_ctx_clear(ctx); fmpz_clear(prime);
    puts("Quotient back-substitution tests passed");
    return 0;
}

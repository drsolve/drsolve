/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "polynomial_system_solver.h"

void test_polynomial_solver(void) {
    printf("\n=== Testing Enhanced Polynomial System Solver ===\n");
    
    fq_nmod_ctx_t ctx;
    mp_limb_t prime = 7;
    fmpz_t p;
    fmpz_init_set_ui(p, prime);
    fq_nmod_ctx_init(ctx, p, 1, "t");  
    
    // Test 1: Simple linear system
    printf("\n--- Test 1: Simple linear system ---\n");
    char *linear_polys[] = {
        "x + y - 3",
        "2*x - y - 0"
    };
    
    polynomial_solutions_t *sols1 = solve_polynomial_system_array(linear_polys, 2, ctx);
    print_polynomial_solutions(sols1);
    polynomial_solutions_clear(sols1);
    free(sols1);
    
    // Test 2: Example system from user (should demonstrate robust equation selection)
    printf("\n--- Test 2: User example system (x^2*y^9*z+x*y, x*y+z, x^4+z+1) ---\n");
    const char *user_system = "x^2*y^9*z+x*y, x*y+z, x^4+z+1";
    
    polynomial_solutions_t *sols2 = solve_polynomial_system_string(user_system, ctx);
    print_polynomial_solutions(sols2);
    polynomial_solutions_clear(sols2);
    free(sols2);
    
    // Test 3: System that should trigger "dimension > 0" error
    printf("\n--- Test 3: System with dimension > 0 ---\n");
    const char *high_dim_system = "x*y, x*y, x^2 + 1";
    
    polynomial_solutions_t *sols3 = solve_polynomial_system_string(high_dim_system, ctx);
    print_polynomial_solutions(sols3);
    polynomial_solutions_clear(sols3);
    free(sols3);
    
    printf("=== Enhanced Testing Complete ===\n");
}

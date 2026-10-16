/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "rational_system_solver.h"

void test_rational_polynomial_solver(void) {
    printf("\n=== Testing Rational Polynomial System Solver ===\n");
    
    char *test1[] = {
        "x + y - 3",
        "x - y + 1"
    };
    
    printf("\nTest 1: Linear system\n");
    rational_solutions_t *sols1 = solve_rational_polynomial_system_array(test1, 2);
    if (sols1) {
        print_rational_solutions(sols1);
        rational_solutions_clear(sols1);
        free(sols1);
    }
    
    char *test2[] = {
        "x^2 - 2",
        "x - 1"
    };
    
    printf("\nTest 2: Quadratic with no rational solution\n");
    rational_solutions_t *sols2 = solve_rational_polynomial_system_array(test2, 2);
    if (sols2) {
        print_rational_solutions(sols2);
        rational_solutions_clear(sols2);
        free(sols2);
    }

    char *test3[] = {
        "x^2 - 2"
    };

    printf("\nTest 3: Univariate irrational real roots\n");
    rational_solutions_t *sols3 = solve_rational_polynomial_system_array(test3, 1);
    if (sols3) {
        print_rational_solutions(sols3);
        rational_solutions_clear(sols3);
        free(sols3);
    }
    
    printf("\n=== Test Complete ===\n");
}


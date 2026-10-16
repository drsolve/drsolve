/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "fq_nmod_roots.h"

// Benchmark functions

static double benchmark_nmod_roots(slong degree, int num_tests) {
    printf("\n=== nmod_poly CZ Root Finding Test ===\n");
    
    nmod_poly_t poly;
    nmod_roots_t roots;
    flint_rand_t state;
    
    nmod_poly_init(poly, PRIME);
    nmod_roots_init(roots);
    flint_rand_init(state);
    
    double total_time = 0.0;
    slong total_roots = 0;
    
    for (int i = 0; i < num_tests; i++) {
        generate_nmod_poly(poly, state, degree, PRIME);
        
        roots->num = 0;
        double start = get_time_roots();
        slong num_roots = our_nmod_poly_roots(roots, poly, 1);
        double end = get_time_roots();
        
        printf("nmod test %d: %.6f seconds, found %ld roots\n", i + 1, end - start, num_roots);
        total_time += (end - start);
        total_roots += num_roots;
    }
    
    double avg = total_time / num_tests;
    printf("nmod average: %.6f seconds, average roots: %.1f\n", avg, (double)total_roots / num_tests);
    
    nmod_poly_clear(poly);
    nmod_roots_clear(roots);
    flint_rand_clear(state);
    
    return avg;
}

static double benchmark_fq_nmod_roots(slong degree, slong extension, int num_tests) {
    printf("\n=== fq_nmod_poly CZ Root Finding Test (F_{%d^%ld}) ===\n", 2, extension);
    
    fq_nmod_ctx_t ctx;
    fq_nmod_poly_t poly;
    fq_nmod_roots_t roots;
    flint_rand_t state;
    
    // Initialize finite field context
    fq_nmod_ctx_init_ui(ctx, 2, extension, "a");
    fq_nmod_poly_init(poly, ctx);
    fq_nmod_roots_init(roots, ctx);
    flint_rand_init(state);
    
    double total_time = 0.0;
    slong total_roots = 0;
    
    for (int i = 0; i < num_tests; i++) {
        generate_fq_nmod_poly(poly, state, degree, ctx);
        
        roots->num = 0;
        double start = get_time_roots();
        slong num_roots = our_fq_nmod_poly_roots(roots, poly, 1, ctx);
        double end = get_time_roots();
        
        printf("fq_nmod test %d: %.6f seconds, found %ld roots\n", i + 1, end - start, num_roots);
        total_time += (end - start);
        total_roots += num_roots;
    }
    
    double avg = total_time / num_tests;
    printf("fq_nmod average: %.6f seconds, average roots: %.1f\n", avg, (double)total_roots / num_tests);
    
    fq_nmod_poly_clear(poly, ctx);
    fq_nmod_roots_clear(roots, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_rand_clear(state);
    
    return avg;
}

static double benchmark_flint_fq_nmod_factor(slong degree, slong extension, int num_tests) {
    printf("\n=== FLINT fq_nmod_poly Factorization Test (F_{%d^%ld}) ===\n", 2, extension);
    
    fq_nmod_ctx_t ctx;
    fq_nmod_poly_t poly;
    fq_nmod_poly_factor_t factors;
    flint_rand_t state;
    
    fq_nmod_ctx_init_ui(ctx, 2, extension, "a");
    fq_nmod_poly_init(poly, ctx);
    fq_nmod_poly_factor_init(factors, ctx);
    flint_rand_init(state);
    
    double total_time = 0.0;
    slong total_roots = 0;
    
    for (int i = 0; i < num_tests; i++) {
        generate_fq_nmod_poly(poly, state, degree, ctx);
        
        fq_nmod_poly_factor_clear(factors, ctx);
        fq_nmod_poly_factor_init(factors, ctx);
        
        fq_nmod_t lead_coeff;
        fq_nmod_init(lead_coeff, ctx);
        
        double start = get_time_roots();
        fq_nmod_poly_factor(factors, lead_coeff, poly, ctx);
        double end = get_time_roots();
        
        slong linear_factors = 0;
        for (slong j = 0; j < factors->num; j++) {
            if (fq_nmod_poly_degree(factors->poly + j, ctx) == 1) {
                linear_factors += factors->exp[j];
            }
        }
        
        printf("FLINT fq_nmod test %d: %.6f seconds, found %ld roots\n", i + 1, end - start, linear_factors);
        total_time += (end - start);
        total_roots += linear_factors;
        
        fq_nmod_clear(lead_coeff, ctx);
    }
    
    double avg = total_time / num_tests;
    printf("FLINT fq_nmod average: %.6f seconds, average roots: %.1f\n", avg, (double)total_roots / num_tests);
    
    fq_nmod_poly_clear(poly, ctx);
    fq_nmod_poly_factor_clear(factors, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_rand_clear(state);
    
    return avg;
}

// Test and verification functions

void test_fq_nmod_correctness(void) {
    printf("\n=== fq_nmod_poly Correctness Verification Test ===\n");
    
    fq_nmod_ctx_t ctx;
    fq_nmod_poly_t poly;
    fq_nmod_roots_t our_roots;
    flint_rand_t state;
    
    // Use F_5^2
    fq_nmod_ctx_init_ui(ctx, 5, 2, "a");
    fq_nmod_poly_init(poly, ctx);
    fq_nmod_roots_init(our_roots, ctx);
    flint_rand_init(state);
    
    printf("Test field: F_{5^2} = F_25\n");
    
    // Construct a simple polynomial: x^2 - 1 = (x-1)(x+1)
    fq_nmod_poly_zero(poly, ctx);
    fq_nmod_t one, neg_one;
    fq_nmod_init(one, ctx);
    fq_nmod_init(neg_one, ctx);
    fq_nmod_one(one, ctx);
    fq_nmod_set_ui(neg_one, 4, ctx);  // -1 in F_5 is 4
    
    fq_nmod_poly_set_coeff(poly, 2, one, ctx);      // x^2
    fq_nmod_poly_set_coeff(poly, 0, neg_one, ctx);  // -1 (in F_5, -1 = 4)
    
    printf("Test polynomial: x^2 - 1\n");
    printf("Expected roots: 1, 4 (in F_5)\n");
    
    our_roots->num = 0;
    slong num_roots = our_fq_nmod_poly_roots(our_roots, poly, 1, ctx);
    
    printf("Our found roots: ");
    for (slong i = 0; i < our_roots->num; i++) {
        fq_nmod_print_pretty(our_roots->roots + i, ctx);
        printf(" ");
    }
    printf("(total %ld)\n", num_roots);
    
    // Verify root correctness
    for (slong i = 0; i < our_roots->num; i++) {
        fq_nmod_t value;
        fq_nmod_init(value, ctx);
        fq_nmod_poly_evaluate_fq_nmod(value, poly, our_roots->roots + i, ctx);
        printf("Verify f(");
        fq_nmod_print_pretty(our_roots->roots + i, ctx);
        printf(") = ");
        fq_nmod_print_pretty(value, ctx);
        printf("\n");
        fq_nmod_clear(value, ctx);
    }
    
    fq_nmod_clear(one, ctx);
    fq_nmod_clear(neg_one, ctx);
    fq_nmod_poly_clear(poly, ctx);
    fq_nmod_roots_clear(our_roots, ctx);
    fq_nmod_ctx_clear(ctx);
    flint_rand_clear(state);
}

void test_fq_nmod_benchmarks(void) {
    printf("=== Fixed Version Unified CZ Root Finding Algorithm Test ===\n");
    printf("Comparing nmod_poly and fq_nmod_poly versions\n");
    printf("====================================\n");
    
    test_fq_nmod_correctness();
    
    slong degrees[] = {100, 500, 1000};
    int num_degrees = sizeof(degrees) / sizeof(degrees[0]);
    int num_tests = 3;
    
    printf("\n=== nmod_poly Test (F_%llu) ===\n", PRIME);
    for (int i = 0; i < num_degrees; i++) {
        printf("\n--- Degree %ld ---\n", degrees[i]);
        benchmark_nmod_roots(degrees[i], num_tests);
    }
    
    printf("\n=== fq_nmod_poly Test (F_{%llu^2}) ===\n", SMALL_PRIME);
    for (int i = 0; i < num_degrees; i++) {
        printf("\n--- Degree %ld ---\n", degrees[i]);
        double our_time = benchmark_fq_nmod_roots(degrees[i], 8, num_tests);
        double flint_time = benchmark_flint_fq_nmod_factor(degrees[i], 8, num_tests);
        
        double ratio = (flint_time > 0) ? our_time / flint_time : 0;
        printf("Ratio (ours/FLINT): %.2fx\n", ratio);
    }
}

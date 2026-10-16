/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "resultant_with_ideal_reduction.h"

void test_iterative_elimination_str2(void) {
    printf("\n================================================\n");
    printf("Test: Simplified Iterative Elimination\n");
    printf("================================================\n");
    
    fq_nmod_ctx_t ctx;
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, 2147483489);
    fq_nmod_ctx_init(ctx, p, 1, "t");
    fmpz_clear(p);
    
    // Define ideal generators in equation format
    const char *ideal_gens_str = "a2^3 = 2*a1 + 1; a3^3 = a1*a2 + 3; a4^3 = a1 + a2*a3 + 5";
    const char *all_vars_str = "a1, a2, a3, a4";
    
    // Step 1: Eliminate a4
    const char *step1_polys_str = "a1^2 + a2^2 + a3^2 + a4^2 - 100, a4^3 - a1 - a2*a3 - 5";
    const char *step1_elim_str = "a4";
    
    char *result1 = elimination_with_ideal_reduction_str(step1_polys_str, step1_elim_str, 
                                                  ideal_gens_str, ctx);
    printf("After eliminating a4: %s\n\n", result1);
    
    // Step 2: Eliminate a3
    char step2_polys_str[10000];
    snprintf(step2_polys_str, sizeof(step2_polys_str), "%s, a3^3 - a1*a2 - 3", result1);
    const char *step2_elim_str = "a3";
    
    char *result2 = elimination_with_ideal_reduction_str(step2_polys_str, step2_elim_str,
                                                  ideal_gens_str, ctx);
    printf("After eliminating a3: %s\n\n", result2);
    
    // Step 3: Eliminate a2
    char step3_polys_str[10000];
    snprintf(step3_polys_str, sizeof(step3_polys_str), "%s, a2^3 - 2*a1 - 1", result2);
    const char *step3_elim_str = "a2";
    
    char *final_result = elimination_with_ideal_reduction_str(step3_polys_str, step3_elim_str,
                                                       ideal_gens_str, ctx);
    printf("Final univariate polynomial in a1: %s\n", final_result);
    
    // Cleanup
    free(result1);
    free(result2);
    free(final_result);
    fq_nmod_ctx_clear(ctx);
}

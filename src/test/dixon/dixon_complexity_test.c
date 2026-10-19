/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "dixon_complexity.h"
#include <assert.h>
#include <math.h>

/* Count monomial vectors independently by multiplying truncated geometric
 * series. Then enumerate column subsets instead of using binomial formulas. */
static void check_mq_layered_bounds(void)
{
    fmpz_t q;
    fmpz_init_set_ui(q, 257);
    for (int n = 2; n <= 10; n++) {
        long degrees[10];
        for (int i = 0; i < n; i++) degrees[i] = 2;
        dixon_complexity_report_t report;
        dixon_complexity_report_from_degrees(&report, degrees, n, n, n-1, 1, q, q, DIXON_OMEGA, 0);
        assert(report.step1_mq_bounds_applicable);
        unsigned long long terms[12] = {1}, cumulative[12];
        for (int v = 0; v < 2*n-1; v++)
            for (int d = 1; d <= n+1; d++) terms[d] += terms[d-1];
        cumulative[0] = terms[0];
        for (int d = 1; d <= n+1; d++) cumulative[d] = cumulative[d-1] + terms[d];
        long double work = 0;
        for (unsigned mask = 1; mask < (1U << n); mask++) {
            int k = __builtin_popcount(mask);
            int entry_terms = k == n ? (n+1)*(n+2)/2 : 2*n;
            work += (long double) k * entry_terms * cumulative[k-1];
        }
        assert(fabs(report.step1_mq_total_support_log2 - log2((double) cumulative[n+1])) < 1e-9);
        assert(fabs(report.step1_mq_layered_log2 - (double) log2l(work)) < 1e-9);
        assert(report.step1_mq_layered_log2 <= report.step1_mq_uniform_log2);
        double probe=pow(n,DIXON_OMEGA)+3.0*n*(n+1)*(n+2)/2+4.0*n*n*(n-1);
        double transforms=4.0*(2*n-1)*(n+2)*(n+2);
        assert(report.step1_mq_simplex_extension==1);
        assert(fabs(report.step1_mq_simplex_log2-log2((double)cumulative[n+1]*(probe+transforms)))<1e-9);
        assert(isfinite(report.step4_mq_schur_log2));
        assert(report.step4_mq_schur_log2>=report.step4_mq_core_log2);
        assert(report.step4_mq_schur_log2>=report.step4_mq_schur_formation_log2);
        assert(report.step4_mq_schur_log2>=report.step4_mq_verification_log2);

        /* New diagnostics must not change the original method selection. */
        double legacy[] = {report.step1_direct_log2, report.step1_direct_mpoly_log2,
            report.step1_direct_mpoly_split_log2, report.step1_bareiss_log2,
            report.step1_ordinary_interp_log2, report.step1_hnf_log2,
            report.step1_sparse_log2, report.step12_recursive_log2};
        double best = INFINITY;
        for (unsigned i = 0; i < sizeof(legacy)/sizeof(*legacy); i++)
            if (legacy[i] < best) best = legacy[i];
        assert(fabs(report.step1_best_log2 - best) < 1e-9);
        assert(report.step1_mq_peak_layer >= 1 && report.step1_mq_peak_layer <= n);
        if (n == 2) assert(work == 56);
        if (n == 3) assert(work == 864);
        degrees[0] = 3;
        dixon_complexity_report_from_degrees(&report, degrees, n, n, n-1, 1, q, q, DIXON_OMEGA, 0);
        assert(!report.step1_mq_bounds_applicable && isinf(report.step1_mq_layered_log2));
        degrees[0] = 2;
        dixon_complexity_report_from_degrees(&report, degrees, n, n+1, n-1, 2, q, q, DIXON_OMEGA, 0);
        assert(!report.step1_mq_bounds_applicable && isinf(report.step1_mq_layered_log2));
    }
    long small_degrees[]={2,2,2,2,2,2};
    dixon_complexity_report_t large,small;
    dixon_complexity_report_from_degrees(&large,small_degrees,6,6,5,1,q,q,DIXON_OMEGA,0);
    fmpz_set_ui(q,2);
    dixon_complexity_report_from_degrees(&small,small_degrees,6,6,5,1,q,q,DIXON_OMEGA,0);
    assert(small.step1_mq_simplex_extension==3);
    assert(fabs(small.step1_mq_simplex_log2-large.step1_mq_simplex_log2-2*log2(3))<1e-9);
    fmpz_clear(q);
    puts("MQ layered estimates: monomial/subset enumeration and applicability checks passed");
}

static void check_mixed_schur_bounds(void)
{
    fmpz_t q; fmpz_init_set_ui(q,1073741831);
    const long cases[2][4]={{3,9,9,9},{9,3,9,9}};
    for(int i=0;i<2;i++) {
        dixon_complexity_report_t r;
        dixon_complexity_report_from_degrees(&r,cases[i],4,4,3,1,q,q,DIXON_OMEGA,1);
        /* Independently expanded positive Hilbert series for (3,9,9,9). */
        const long H[]={1,3,6,9,12,15,18,21,24,24,21,15,9,3};
        long h=0,delta=0;
        for(int d=0;d<14;d++){h+=H[d];delta+=H[d]*(27-2*d);}
        assert(h==181 && delta==2187);
        assert(r.step4_mq_core_size==h && r.step4_schur_degree_bound==delta);
        assert(fabs(r.step4_mq_core_log2-(DIXON_OMEGA*log2(181.0)+log2(13.0)))<1e-9);
        assert(isfinite(r.step4_mq_schur_log2));
        assert(r.step4_mq_schur_log2>=r.step4_mq_schur_formation_log2);
        assert(r.step4_mq_schur_log2>=r.step4_mq_core_log2);
        assert(r.step4_mq_schur_log2>=r.step4_mq_verification_log2);
        assert(r.step4_mq_schur_log2==fmax(r.step4_mq_schur_formation_log2,
            fmax(r.step4_mq_core_log2,r.step4_mq_verification_log2)));
        assert(fabs(r.step3_schur_selection_log2-DIXON_OMEGA*log2(869.0))<1e-9);
        /* Cached degree checks scan matrix entries, not all coefficients. */
        assert(fabs(r.step4_mq_verification_log2-2*log2(869.0))<1e-9);
        assert(!r.step1_mq_bounds_applicable && isinf(r.step1_mq_simplex_log2));
        dixon_complexity_report_from_degrees(&r,cases[i],4,5,3,2,q,q,DIXON_OMEGA,1);
        assert(isinf(r.step4_mq_schur_log2));
    }
    dixon_complexity_report_t r;
    dixon_complexity_report_from_degrees(&r,cases[0],4,4,3,1,q,q,3.0,1);
    /* Classical field matrix multiplication, quasi-linear polynomial degree:
     * direct summation of block products gives this integer operation proxy. */
    assert(fabs(r.step4_mq_schur_formation_log2-log2(651887135.0))<1e-9);
    assert(fabs(r.step3_schur_selection_log2-3.0*log2(869.0))<1e-9);
    fmpz_zero(q);
    dixon_complexity_report_from_degrees(&r,cases[0],4,4,3,1,q,q,DIXON_OMEGA,1);
    assert(isinf(r.step4_mq_schur_log2));
    fmpz_clear(q);
    puts("Mixed-degree Schur: Hilbert core, graded degree budget, costs and applicability PASS");
}

int test_dixon_complexity(void) {
    check_mixed_schur_bounds();
    check_mq_layered_bounds();
    // Test data
    long a1[] = {1000, 1000, 1000, 1001, 1002, 1003};
    int len = sizeof(a1) / sizeof(a1[0]);
    double omega = DIXON_OMEGA;
    
    printf("Dixon Complexity Results (Hessenberg Method):\n");
    for (int n = 5; n < 10; n++) {
        double complexity = dixon_complexity(a1, len, n, omega);
        printf("n=%d: %.6f\n", n, complexity);
    }
    
    // Test dixon_size function
    printf("\nTesting dixon_size with Hessenberg method:\n");
    fmpz_t test_result;
    fmpz_init(test_result);
    dixon_size(test_result, a1, len, 1);
    fmpz_clear(test_result);
    
    return 0;
}

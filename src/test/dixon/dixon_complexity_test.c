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
    /* Newton inversion and products use the paper's sigma*rho^omega budget. */
    assert(fabs(r.step4_mq_schur_formation_log2-(log2(27.0)+3.0*log2(869.0)))<1e-9);
    assert(fabs(r.step3_schur_selection_log2-3.0*log2(869.0))<1e-9);
    fmpz_zero(q);
    dixon_complexity_report_from_degrees(&r,cases[0],4,4,3,1,q,q,DIXON_OMEGA,1);
    assert(isfinite(r.step4_mq_schur_log2));
    assert(isinf(r.step1_sparse_log2));
    fmpz_clear(q);
    puts("Mixed-degree Schur: Hilbert core, graded degree budget, costs and applicability PASS");
}

/* Independent reference: enumerate prefix extensions and expand Hilbert
 * numerator subsets. This retains the old slow counting method for small cases. */
static long small_binomial(long n, long k)
{
    if (n < k || k < 0) return 0;
    if (k > n-k) k=n-k;
    long value=1;
    for (long j=1;j<=k;j++) value=value*(n-k+j)/j;
    return value;
}

static void check_profile_reference(const long *degrees, int n)
{
    long sorted[8], R[128]={1}, next[128], H[128]={0};
    for(int i=0;i<n;i++) sorted[i]=degrees[i];
    for(int i=0;i<n;i++) for(int j=i+1;j<n;j++)
        if(sorted[j]>sorted[i]) {long t=sorted[i];sorted[i]=sorted[j];sorted[j]=t;}
    int cap=0;
    for(int k=0;k<n-1;k++) {
        int old=cap;cap+=sorted[k]-1;
        assert(cap<128);
        memset(next,0,sizeof(next));
        for(int a=0;a<=old;a++) for(int b=0;a+b<=cap;b++) next[a+b]+=R[a];
        memcpy(R,next,sizeof(R));
    }
    int sigma=1;
    for(int i=0;i<n;i++) sigma+=degrees[i]-1;
    int nh=0;
    for(int e=0;e<128;e++) {
        long value=0;
        for(unsigned subset=0;subset<(1U<<n);subset++) {
            int rem=e;
            for(int i=0;i<n;i++) if(subset&(1U<<i)) rem-=degrees[i];
            if(rem>=0) value+=(__builtin_popcount(subset)%2?-1:1)*small_binomial(n-2+rem,n-2);
        }
        if(value<=0) break;
        H[nh++]=value;
    }
    long h=0,delta=0,rank=0;
    for(int e=0;e<nh;e++) {h+=H[e];delta+=H[e]*(sigma-2*e);}
    rank=h;
    for(int e=0;e<=cap;e++) {
        int other=sigma-e,low=FLINT_MIN(e,other);
        if(other>=0 && other<=cap)
            rank+=FLINT_MAX(0,FLINT_MIN(R[e],R[other])-(low<nh?H[low]:0));
    }
    slong *actual_R=NULL,*actual_H=NULL,nr,actual_nh,actual_sigma,actual_rank;
    assert(dixon_rank_profile_from_degrees(&actual_R,&nr,&actual_H,&actual_nh,
        &actual_sigma,&actual_rank,degrees,n,n-1));
    assert(nr==cap+1 && actual_nh==nh && actual_sigma==sigma && actual_rank==rank);
    for(int e=0;e<nr;e++) assert(actual_R[e]==R[e]);
    for(int e=0;e<nh;e++) assert(actual_H[e]==H[e]);
    flint_free(actual_R);flint_free(actual_H);
    fmpz_t q;fmpz_init_set_ui(q,65537);
    dixon_complexity_report_t r;
    dixon_complexity_report_from_degrees(&r,degrees,n,n,n-1,1,q,q,DIXON_OMEGA,1);
    assert(r.step4_rank_size==rank && r.step4_mq_core_size==h);
    assert(r.step4_schur_degree_bound==delta);
    assert(fabs(r.step4_mq_schur_formation_log2-(log2((double)sigma)+DIXON_OMEGA*log2((double)rank)))<1e-9);
    assert(fabs(r.step4_mq_core_log2-(DIXON_OMEGA*log2((double)h)+log2((double)((delta+h-1)/h))))<1e-9);
    assert(fabs(r.step3_schur_selection_log2-fmax(DIXON_OMEGA*log2((double)rank),
        log2(sigma+1.0)+2*log2((double)rank)))<1e-9);
    fmpz_clear(q);
}

static void check_fast_rank_counts(void)
{
    long degrees[8];
    for(int n=2;n<=8;n++) for(int d=1;d<=12;d++) {
        for(int i=0;i<n;i++) degrees[i]=d;
        check_profile_reference(degrees,n);
    }
    for(int n=3;n<=8;n++) for(int trial=0;trial<12;trial++) {
        for(int i=0;i<n;i++) degrees[i]=1+(trial*3+i*5+i*i)%7;
        check_profile_reference(degrees,n);
    }
    /* Finite sums support degrees whose rank and core exceed a machine word. */
    if(sizeof(long)>=8) {
        const long big[]={4000000000L,4000000000L,4000000000L};
        fmpz_t q,rank;fmpz_init_set_ui(q,65537);fmpz_init(rank);
        fmpz_set_ui(rank,4000000000UL);fmpz_mul(rank,rank,rank);
        fmpz_mul_ui(rank,rank,5);fmpz_fdiv_q_ui(rank,rank,4);
        dixon_complexity_report_t r;
        dixon_complexity_report_from_degrees(&r,big,3,3,2,1,q,q,DIXON_OMEGA,1);
        assert(r.step4_rank_model_applicable && r.step4_rank_size==WORD_MAX);
        assert(fabs(r.step4_rank_size_log2-log2(fmpz_get_d(rank)))<1e-9);
        assert(isfinite(r.step4_mq_schur_log2));
        fmpz_clear(q);fmpz_clear(rank);
    }
    {
        const long degrees[]={2,2,2,2,2,2};
        fmpz_t p,q;fmpz_init_set_si(p,-1);fmpz_init(q);
        dixon_complexity_report_t generic,finite;
        dixon_complexity_report_from_degrees(&generic,degrees,6,6,5,1,p,q,DIXON_OMEGA,1);
        fmpz_set_ui(p,65537);fmpz_set_ui(q,65537);
        dixon_complexity_report_from_degrees(&finite,degrees,6,6,5,1,p,q,DIXON_OMEGA,1);
        assert(generic.step1_mq_simplex_extension==1);
        assert(generic.step1_mq_simplex_log2==finite.step1_mq_simplex_log2);
        assert(generic.step4_mq_schur_log2==finite.step4_mq_schur_log2);
        assert(generic.step4_rank_size==finite.step4_rank_size);
        assert(generic.step1_sparse_log2==generic.step1_sparse_term_bound_log2+generic.step1_sparse_slp_length_log2);
        assert(generic.step4_sparse_log2==generic.step4_sparse_term_bound_log2+generic.step4_sparse_slp_length_log2);
        fmpz_clear(p);fmpz_clear(q);
    }
    puts("Fast rank counts: independent uniform/mixed profiles and large-degree budgets PASS");
}

int test_dixon_complexity(void) {
    check_fast_rank_counts();
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

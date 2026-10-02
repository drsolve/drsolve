/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Compare every recursive coefficient against both the legacy arithmetic
 * and the independent cancellation determinant, including sparse inputs. */
#include "../dixon/dixon_recursive.c"
#include <assert.h>
#include <flint/ulong_extras.h>

static void random_terms(unified_mpoly_struct *p, slong *exp, slong pos, slong left,
                         flint_rand_t rng, int sparse)
{
    if (pos < p->nvars+p->npars) {
        for (slong d=0; d<=left; d++) {
            exp[pos]=d;
            random_terms(p,exp,pos+1,left-d,rng,sparse);
        }
        return;
    }
    ulong prime=fq_nmod_ctx_prime(p->ctx);
    ulong value=n_randint(rng,prime);
    if (sparse && n_randint(rng,3)) value=0;
    fq_nmod_t c; fq_nmod_init(c,p->ctx); fq_nmod_set_ui(c,value,p->ctx);
    if (value)
        dr_mpoly_add_term_fast(p, exp, exp + p->nvars, c);
    fq_nmod_clear(c,p->ctx);
}

static int equal_poly(const unified_mpoly_struct *a, const unified_mpoly_struct *b)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx,a->nvars+a->npars,ORD_LEX,fq_nmod_ctx_prime(a->ctx));
    nmod_mpoly_t x,y; nmod_mpoly_init(x,ctx); nmod_mpoly_init(y,ctx);
    dr_mpoly_to_nmod_mpoly(x, a, ctx);
    dr_mpoly_to_nmod_mpoly(y, b, ctx);
    int same=nmod_mpoly_equal(x,y,ctx);
    nmod_mpoly_clear(y,ctx); nmod_mpoly_clear(x,ctx); nmod_mpoly_ctx_clear(ctx);
    return same;
}

static void check(ulong prime, slong degree, int sparse, slong npars)
{
    fq_nmod_ctx_t ctx; fq_nmod_ctx_init_ui(ctx,prime,1,"a");
    flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng,prime+degree,17+sparse);
    unified_mpoly_struct p[3] = {0};
    const unified_mpoly_struct *ptrs[3];
    slong exp[4]={0}, degrees[2];
    for (slong i=0; i<3; i++) {
        dr_mpoly_init(p + i, 2, npars, ctx);
        ptrs[i] = p + i;
        random_terms(p+i,exp,0,degree,rng,sparse);
        /* Keep both recursive degree bounds positive even for sparse F_2. */
        fq_nmod_t one; fq_nmod_init(one,ctx); fq_nmod_one(one,ctx);
        slong x[2]={degree,0}, y[2]={0,degree}, par[2]={0,0};
        dr_mpoly_add_term_fast(p + i, x, par, one);
        dr_mpoly_add_term_fast(p + i, y, par, one);
        fq_nmod_clear(one,ctx);
    }
    fast_dixon_compute_degree_bounds(degrees,ptrs,3,2);
    fast_dixon_matrix_t legacy,native;
    fast_dixon_profile_reset(2); fast_dixon_subproblem_cache_reset_all();
    setenv("DRSOLVE_FAST_NATIVE","0",1);
    fast_dixon_build_matrix(&legacy,ptrs,degrees,2,0,npars,ctx);
    setenv("DRSOLVE_FAST_NATIVE","1",1);
    fast_dixon_build_matrix(&native,ptrs,degrees,2,0,npars,ctx);
    assert(native.rows==legacy.rows && native.cols==legacy.cols);
    for (slong i=0; i<native.rows*native.cols; i++)
        assert(equal_poly(native.entries+i,legacy.entries+i));

    /* Reconstruct the polynomial using the recursive matrix's raw labels. */
    unified_mpoly_struct actual = {0}, expected = {0}, **cancel, **modified;
    dr_mpoly_init(&actual, 4, npars, ctx);
    for (slong r=0; r<native.rows; r++) for (slong c=0; c<native.cols; c++) {
        slong label[4]={r/(2*degrees[1]),r%(2*degrees[1]),c/degrees[1],c%degrees[1]};
        const unified_mpoly_struct *entry = &FAST_DIXON_ENTRY(&native, r, c);
        for (slong t = 0; t < dr_mpoly_length(entry); t++) {
            DR_MPOLY_TERM(term_1, entry, t);
            dr_mpoly_add_term_fast(&actual, label, term_1.par_exp, term_1.coeff);
        }
    }
    build_fq_cancellation_matrix(&cancel, p, 2, npars);
    perform_fq_matrix_row_operations(&modified, &cancel, 2, npars);
    compute_fq_cancel_matrix_det(&expected,modified,2,npars,DET_METHOD_RECURSIVE);
    /* Recursive convention can differ by the global Dixon determinant sign. */
    if (!equal_poly(&actual,&expected)) {
        dr_mpoly_neg(&actual, &actual);
        assert(equal_poly(&actual,&expected));
    }
    for (slong i=0; i<3; i++) {
        for (slong j = 0; j < 3; j++) {
            dr_mpoly_clear(&cancel[i][j]);
            dr_mpoly_clear(&modified[i][j]);
        }
        flint_free(cancel[i]); flint_free(modified[i]);
    }
    flint_free(cancel); flint_free(modified);
    dr_mpoly_clear(&expected);
    dr_mpoly_clear(&actual);
    fast_dixon_matrix_clear(&native); fast_dixon_matrix_clear(&legacy);
    for (slong i = 0; i < 3; i++)
        dr_mpoly_clear(p + i);
    flint_rand_clear(rng); fq_nmod_ctx_clear(ctx);
}

int main(void)
{
    g_dixon_verbose_level=0; omp_set_num_threads(4);
    for (slong d=2; d<=5; d++) {
        check(65537,d,0,1); check(101,d,1,1);
    }
    check(2,3,0,1); check(65537,3,1,2);
    unsetenv("DRSOLVE_FAST_NATIVE");
    puts("Recursive native/legacy coefficients and independent cancellation determinant PASS");
    flint_cleanup_master(); return 0;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
#define main mq_native_test_main
#include "mq_gf2n_test.c"
#undef main

static void check_coefficients(const fq_nmod_ctx_t ctx, flint_rand_t state)
{
    slong d = fq_nmod_ctx_degree(ctx);
    ulong *a = flint_malloc(d*sizeof(ulong)), *b = flint_malloc(d*sizeof(ulong));
    ulong *actual = flint_malloc(d*sizeof(ulong));
    fq_nmod_t x,y,initial,expected,scratch;
    fq_nmod_init(x,ctx); fq_nmod_init(y,ctx); fq_nmod_init(initial,ctx);
    fq_nmod_init(expected,ctx); fq_nmod_init(scratch,ctx);
    for (int trial = 0; trial < 40; trial++) {
        fq_nmod_randtest(x,state,ctx); fq_nmod_randtest(y,state,ctx);
        fq_nmod_randtest(initial,state,ctx);
        if (trial < 3) fq_nmod_set_ui(x,trial,ctx);
        if (trial == 3) { /* Exercises worst-case full-word dot products. */
            for (slong j = 0; j < d; j++) {
                nmod_poly_set_coeff_ui(x,j,fq_nmod_ctx_prime(ctx)-1);
                nmod_poly_set_coeff_ui(y,j,fq_nmod_ctx_prime(ctx)-1);
            }
        }
        mq_fq_store(a,x,d); mq_fq_store(b,y,d);
        mq_fq_operator op; mq_fq_operator_prepare(&op,a,scratch,ctx);
        for (int subtract = 0; subtract < 2; subtract++) {
            mq_fq_store(actual,initial,d);
            fq_nmod_mul(expected,x,y,ctx);
            if (subtract) fq_nmod_sub(expected,initial,expected,ctx);
            else fq_nmod_add(expected,initial,expected,ctx);
            mq_fq_addmul(actual,&op,b,subtract,scratch,ctx);
            fq_nmod_struct view = mq_fq_view(actual,ctx);
            assert(fq_nmod_equal(&view,expected,ctx));
        }
    }
    fq_nmod_clear(x,ctx); fq_nmod_clear(y,ctx); fq_nmod_clear(initial,ctx);
    fq_nmod_clear(expected,ctx); fq_nmod_clear(scratch,ctx);
    flint_free(a); flint_free(b); flint_free(actual);
}

static void check_rank_retry(const fq_nmod_ctx_t ctx)
{
    unified_mpoly_struct diagonal[4] = {0};
    unified_mpoly_struct *cells[16] = {0}, **matrix[4];
    slong indices[] = {0,1,2,3}, power = 1;
    fq_nmod_t point[1], coefficient;
    fq_nmod_init(point[0],ctx); fq_nmod_init(coefficient,ctx);
    fq_nmod_set_ui(point[0],2,ctx);
    for (slong i = 0; i < 4; i++) {
        matrix[i] = cells+4*i; cells[5*i] = diagonal+i;
        dr_mpoly_init(diagonal+i,0,1,ctx);
        fq_nmod_one(coefficient,ctx); power = 1;
        dr_mpoly_add_term_fast(diagonal+i,NULL,&power,coefficient);
        fq_nmod_neg(coefficient,point[0],ctx); power = 0;
        dr_mpoly_add_term_fast(diagonal+i,NULL,&power,coefficient);
    }
    assert(dixon_fq_candidate_rank(matrix,indices,indices,4,1,point,ctx) == 4);
    fq_nmod_set_ui(coefficient,2,ctx);
    assert(!fq_nmod_equal(point[0],coefficient,ctx));
    /* A truly deficient candidate must not manufacture rank or change point. */
    matrix[3][3] = NULL; fq_nmod_gen(point[0],ctx); fq_nmod_set(coefficient,point[0],ctx);
    assert(dixon_fq_candidate_rank(matrix,indices,indices,4,1,point,ctx) == 3);
    assert(fq_nmod_equal(point[0],coefficient,ctx));
    for (slong i = 0; i < 4; i++) dr_mpoly_clear(diagonal+i);
    fq_nmod_clear(point[0],ctx); fq_nmod_clear(coefficient,ctx);
}

int main(void)
{
    const ulong primes[] = {3,5,17,257,65537,2,2,3,UWORD(18446744073709551557)};
    const slong degrees[] = {2,3,2,3,2,7,9,9,2};
    g_dixon_verbose_level = 0; g_dixon_det_cache_limit = 1024;
    g_dixon_mq_step1_filter = 1; g_dixon_mq_step1_shared = 1;
    flint_rand_t state; flint_rand_init(state); flint_rand_set_seed(state,713,912);
    for (size_t f = 0; f < sizeof(primes)/sizeof(*primes); f++) {
        fq_nmod_ctx_t ctx; fq_nmod_ctx_init_ui(ctx,primes[f],degrees[f],"a");
        check_coefficients(ctx,state); check_projection(ctx,state); check_schur(ctx,state); check_rank_retry(ctx);
        printf("MQ F_%lu^%ld: scalar arithmetic, signed projection, rank, Schur and fallback PASS\n",primes[f],degrees[f]);
        fq_nmod_ctx_clear(ctx);
    }
    /* Same characteristic/degree, distinct defining polynomial and basis. */
    nmod_poly_t modulus; nmod_poly_init(modulus,17);
    nmod_poly_set_coeff_ui(modulus,2,1); nmod_poly_set_coeff_ui(modulus,0,3);
    assert(nmod_poly_is_irreducible(modulus));
    fq_nmod_ctx_t alternate; fq_nmod_ctx_init_modulus(alternate,modulus,"b");
    check_coefficients(alternate,state); check_projection(alternate,state); check_schur(alternate,state);
    fq_nmod_ctx_clear(alternate); nmod_poly_clear(modulus);
    flint_rand_clear(state); flint_cleanup_master(); puts("MQ general extension tests PASS"); return 0;
}

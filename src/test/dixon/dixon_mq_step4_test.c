/* SPDX-License-Identifier: GPL-2.0-or-later */
#define main mq_filter_existing_main
#include "dixon_mq_filter_test.c"
#undef main

static void same_result(const unified_mpoly_struct *a, const unified_mpoly_struct *b)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, 1, ORD_LEX, fq_nmod_ctx_prime(a->ctx));
    nmod_mpoly_t x, y;
    nmod_mpoly_init(x, ctx); nmod_mpoly_init(y, ctx);
    dr_mpoly_to_nmod_mpoly(x, a, ctx);
    dr_mpoly_to_nmod_mpoly(y, b, ctx);
    assert(nmod_mpoly_equal(x, y, ctx));
    nmod_mpoly_clear(x, ctx); nmod_mpoly_clear(y, ctx); nmod_mpoly_ctx_clear(ctx);
}

static void check_selected_matrix(unified_mpoly_struct *polys, slong m, const long *degrees)
{
    unified_mpoly_struct **matrix, **a, full = {0};
    build_fq_cancellation_matrix(&matrix, polys, m, 1);
    perform_fq_matrix_row_operations(&a, &matrix, m, 1);
    compute_fq_cancel_matrix_det(&full, a, m, 1, DET_METHOD_RECURSIVE);
    unified_mpoly_struct **unused = NULL;
    fq_nmod_poly_mat_t B;
    slong *ri = flint_malloc(dr_mpoly_length(&(full)) * sizeof(slong));
    slong *ci = flint_malloc(dr_mpoly_length(&(full)) * sizeof(slong));
    slong size, content;
    dixon_mq_step4_profile p = {0};
    extract_fq_coefficient_matrix_from_dixon_impl(&unused, &B, ri, ci, &size,
        &content, &full, m, 1, NULL, NULL, NULL, degrees, m+1, 0, &p,NULL,NULL);
    assert(p.size == size && p.size > p.h);
    fq_nmod_poly_t expected, got;
    fq_nmod_poly_init(expected, polys[0].ctx); fq_nmod_poly_init(got, polys[0].ctx);
    fq_nmod_poly_mat_det_iter(expected, B, polys[0].ctx);
    assert(dixon_mq_step4_try(got, B, &p, polys[0].ctx));
    assert(fq_nmod_poly_equal(expected, got, polys[0].ctx));
    /* Independently permute an odd number of rows only. Update the selected
     * labels, as Step 3 does; this catches missing external permutation signs. */
    for (slong j = 0; j < size; j++)
        fq_nmod_poly_swap(fq_nmod_poly_mat_entry(B,0,j), fq_nmod_poly_mat_entry(B,1,j), polys[0].ctx);
    for (slong i = 0; i < size; i++) {
        if (p.rows[i] == 0) p.rows[i] = 1;
        else if (p.rows[i] == 1) p.rows[i] = 0;
    }
    p.odd ^= 1;
    fq_nmod_poly_neg(expected, expected, polys[0].ctx);
    assert(dixon_mq_step4_try(got, B, &p, polys[0].ctx));
    assert(fq_nmod_poly_equal(expected, got, polys[0].ctx));
    /* Rejection leaves the caller's determinant and source intact. */
    p.sigma = 0;
    assert(!dixon_mq_step4_try(got, B, &p, polys[0].ctx));
    assert(fq_nmod_poly_equal(expected, got, polys[0].ctx));
    fq_nmod_poly_mat_det_iter(got, B, polys[0].ctx);
    assert(fq_nmod_poly_equal(expected, got, polys[0].ctx));
    fq_nmod_poly_clear(got, polys[0].ctx); fq_nmod_poly_clear(expected, polys[0].ctx);
    fq_nmod_poly_mat_clear(B, polys[0].ctx); dixon_mq_step4_profile_clear(&p);
    flint_free(ri);
    flint_free(ci);
    dr_mpoly_clear(&full);
    for(slong i=0;i<=m;i++) {
        for (slong j = 0; j <= m; j++) {
            dr_mpoly_clear(&matrix[i][j]);
            dr_mpoly_clear(&a[i][j]);
        }
        flint_free(matrix[i]);flint_free(a[i]);
    }
    flint_free(matrix);flint_free(a);
}

static void add_dense_terms(unified_mpoly_struct *p, slong *exp, slong axis,
                            slong remaining, flint_rand_t state)
{
    if (axis == p->nvars + 1) {
        fq_nmod_t c;
        fq_nmod_init(c, p->ctx);
        fq_nmod_randtest_not_zero(c, state, p->ctx);
        dr_mpoly_add_term_fast(p, exp, exp + p->nvars, c);
        fq_nmod_clear(c, p->ctx);
        return;
    }
    for (slong d = 0; d <= remaining; d++) {
        exp[axis] = d;
        add_dense_terms(p, exp, axis + 1, remaining - d, state);
    }
}

static void check_general_degrees(const long *degrees, slong m,
                                   const fq_nmod_ctx_t ctx, flint_rand_t state)
{
    unified_mpoly_struct *p = flint_calloc(m + 1, sizeof(*p));
    slong *exp = flint_calloc(m + 1, sizeof(*exp));
    for (slong i = 0; i <= m; i++) {
        dr_mpoly_init(p + i, m, 1, ctx);
        add_dense_terms(p + i, exp, 0, degrees[i], state);
    }
    assert(dixon_mq_step4_eligible(p, m, 1));
    assert(!dixon_mq_step4_eligible(p, m, 2));
    check_selected_matrix(p, m, degrees);
    unified_mpoly_struct baseline = {0}, compressed = {0};
    g_dixon_mq_step4_schur = 0;
    fq_dixon_resultant(&baseline, p, m, 1);
    g_dixon_mq_step4_schur = 1;
    fq_dixon_resultant_with_names(&compressed, p, m, 1, NULL, NULL, NULL);
    fq_nmod_mpoly_ctx_t mctx;
    fq_nmod_mpoly_ctx_init(mctx, 1, ORD_LEX, ctx);
    fq_nmod_mpoly_t a, b;
    fq_nmod_mpoly_init(a, mctx); fq_nmod_mpoly_init(b, mctx);
    dr_mpoly_to_fq_nmod_mpoly(a, &baseline, mctx);
    dr_mpoly_to_fq_nmod_mpoly(b, &compressed, mctx);
    assert(fq_nmod_mpoly_equal(a, b, mctx));
    fq_nmod_mpoly_clear(a, mctx); fq_nmod_mpoly_clear(b, mctx);
    fq_nmod_mpoly_ctx_clear(mctx);
    dr_mpoly_clear(&baseline); dr_mpoly_clear(&compressed);
    for (slong i = 0; i <= m; i++) dr_mpoly_clear(p + i);
    flint_free(p); flint_free(exp);
}

int main(void)
{
    omp_set_num_threads(4);g_dixon_verbose_level=0;g_dixon_det_cache_limit=100000;
    flint_rand_t state;flint_rand_init(state);flint_rand_set_seed(state,24092026,808);
    fq_nmod_ctx_t ctx;fq_nmod_ctx_init_ui(ctx,65537,1,"a");
    for(slong m=3;m<=5;m++) {
        unified_mpoly_struct *p = random_mq(m, ctx, state), baseline = {0}, compressed = {0};
        assert(dixon_mq_step4_eligible(p,m,1));
        long degrees[6] = {2,2,2,2,2,2};
        check_selected_matrix(p,m,degrees);
        g_dixon_mq_step1_simplex=0;
        g_dixon_mq_step4_schur=0;
        fq_dixon_resultant(&baseline,p,m,1);
        g_dixon_mq_step1_simplex=1;
        g_dixon_mq_step4_schur=1;
        fq_dixon_resultant_with_names(&compressed,p,m,1,NULL,NULL,NULL);
        same_result(&baseline,&compressed);
        dr_mpoly_clear(&baseline);
        dr_mpoly_clear(&compressed);
        for (slong i = 0; i <= m; i++)
            dr_mpoly_clear(p + i);
        flint_free(p);
    }
    /* Three equations in x,y,t: compare the full determinant, including sign,
     * with checked compression for several higher total degrees. */
    for (slong d = 3; d <= 6; d++) {
        unified_mpoly_struct p[3] = {0}, baseline = {0}, compressed = {0};
        fq_nmod_t c; fq_nmod_init(c, ctx);
        for (slong i = 0; i < 3; i++) {
            dr_mpoly_init(p + i, 2, 1, ctx);
            for (slong x = 0; x <= d; x++)
                for (slong y = 0; y <= d-x; y++)
                    for (slong t = 0; t <= d-x-y; t++) {
                        slong exp[2] = {x,y}, par[1] = {t};
                        fq_nmod_set_ui(c, 1+n_randint(state,65536), ctx);
                        dr_mpoly_add_term_fast(p + i, exp, par, c);
                    }
        }
        assert(dixon_mq_step4_eligible(p,2,1));
        long degrees[3] = {d,d,d};
        check_selected_matrix(p,2,degrees);
        g_dixon_mq_step1_simplex=0;
        g_dixon_mq_step4_schur=0;
        fq_dixon_resultant(&baseline,p,2,1);
        g_dixon_mq_step4_schur=1;
        fq_dixon_resultant_with_names(&compressed,p,2,1,NULL,NULL,NULL);
        same_result(&baseline,&compressed);
        dr_mpoly_clear(&baseline);
        dr_mpoly_clear(&compressed);
        /* A parameter outside the shared degree budget must be rejected. */
        slong exp[2] = {0,0}, par[1] = {d+1};
        dr_mpoly_add_term_fast(p, exp, par, c);
        assert(!dixon_mq_step4_eligible(p,2,1));
        for (slong i = 0; i < 3; i++)
            dr_mpoly_clear(p + i);
        fq_nmod_clear(c,ctx);
    }
    g_dixon_mq_step1_simplex=0;
    const long quartics[] = {4,4,4,4};
    const long mixed[] = {2,3,4,2};
    const long linear[] = {1,2,3,2};
    const long five[] = {3,2,2,2,2};
    check_general_degrees(quartics,3,ctx,state);
    check_general_degrees(mixed,3,ctx,state);
    check_general_degrees(linear,3,ctx,state);
    check_general_degrees(five,4,ctx,state);
    fq_nmod_ctx_t extension;
    fq_nmod_ctx_init_ui(extension,257,2,"a");
    check_general_degrees(mixed,3,extension,state);
    fq_nmod_ctx_clear(extension);
    fq_nmod_ctx_clear(ctx);flint_rand_clear(state);
    puts("MQ/high/mixed-degree Step 4 integration: selected labels, odd row permutation, unchanged fallback input, both resultant APIs PASS");
    flint_cleanup_master();return 0;
}

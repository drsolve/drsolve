/* SPDX-License-Identifier: GPL-2.0-or-later */
#define main mq_filter_existing_main
#include "dixon_mq_filter_test.c"
#undef main

static void same_result(const fq_mvpoly_t *a, const fq_mvpoly_t *b)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, 1, ORD_LEX, fq_nmod_ctx_prime(a->ctx));
    nmod_mpoly_t x, y;
    nmod_mpoly_init(x, ctx); nmod_mpoly_init(y, ctx);
    fq_mvpoly_to_nmod_mpoly(x, a, ctx); fq_mvpoly_to_nmod_mpoly(y, b, ctx);
    assert(nmod_mpoly_equal(x, y, ctx));
    nmod_mpoly_clear(x, ctx); nmod_mpoly_clear(y, ctx); nmod_mpoly_ctx_clear(ctx);
}

static void check_selected_matrix(fq_mvpoly_t *polys, slong m)
{
    fq_mvpoly_t **matrix, **a, full;
    build_fq_cancellation_matrix_mvpoly(&matrix, polys, m, 1);
    perform_fq_matrix_row_operations_mvpoly(&a, &matrix, m, 1);
    compute_fq_cancel_matrix_det(&full, a, m, 1, DET_METHOD_RECURSIVE);
    fq_mvpoly_t **unused = NULL;
    fq_nmod_poly_mat_t B;
    slong *ri = flint_malloc(full.nterms * sizeof(slong));
    slong *ci = flint_malloc(full.nterms * sizeof(slong));
    long deg[8]; for (slong i = 0; i <= m; i++) deg[i] = 2;
    slong size, content;
    dixon_mq_step4_profile p = {0};
    extract_fq_coefficient_matrix_from_dixon_impl(&unused, &B, ri, ci, &size,
        &content, &full, m, 1, NULL, NULL, NULL, deg, m+1, 0, &p);
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
    flint_free(ri); flint_free(ci); fq_mvpoly_clear(&full);
    for(slong i=0;i<=m;i++) {
        for(slong j=0;j<=m;j++){fq_mvpoly_clear(&matrix[i][j]);fq_mvpoly_clear(&a[i][j]);}
        flint_free(matrix[i]);flint_free(a[i]);
    }
    flint_free(matrix);flint_free(a);
}

int main(void)
{
    omp_set_num_threads(4);g_dixon_verbose_level=0;g_dixon_det_cache_limit=100000;
    flint_rand_t state;flint_rand_init(state);flint_rand_set_seed(state,24092026,808);
    fq_nmod_ctx_t ctx;fq_nmod_ctx_init_ui(ctx,65537,1,"a");
    for(slong m=3;m<=5;m++) {
        fq_mvpoly_t *p=random_mq(m,ctx,state),baseline,compressed;
        assert(dixon_mq_step4_eligible(p,m,1));
        check_selected_matrix(p,m);
        g_dixon_mq_step1_simplex=0;
        g_dixon_mq_step4_schur=0;
        fq_dixon_resultant(&baseline,p,m,1);
        g_dixon_mq_step1_simplex=1;
        g_dixon_mq_step4_schur=1;
        fq_dixon_resultant_with_names(&compressed,p,m,1,NULL,NULL,NULL);
        same_result(&baseline,&compressed);
        fq_mvpoly_clear(&baseline);fq_mvpoly_clear(&compressed);
        for(slong i=0;i<=m;i++)fq_mvpoly_clear(p+i);
        flint_free(p);
    }
    g_dixon_mq_step1_simplex=0;
    fq_nmod_ctx_clear(ctx);flint_rand_clear(state);
    puts("MQ Step 4 integration: selected labels, odd row permutation, unchanged fallback input, both resultant APIs PASS");
    flint_cleanup_master();return 0;
}

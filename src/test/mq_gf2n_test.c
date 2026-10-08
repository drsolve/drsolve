/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "../dixon/dixon_flint.c"
#include <assert.h>

static void equal_poly(const unified_mpoly_struct *a, const unified_mpoly_struct *b)
{
    fq_nmod_mpoly_ctx_t ctx;
    fq_nmod_mpoly_ctx_init(ctx, a->nvars+a->npars, ORD_LEX, a->ctx);
    fq_nmod_mpoly_t x, y;
    fq_nmod_mpoly_init(x, ctx); fq_nmod_mpoly_init(y, ctx);
    dr_mpoly_to_fq_nmod_mpoly(x, a, ctx); dr_mpoly_to_fq_nmod_mpoly(y, b, ctx);
    assert(fq_nmod_mpoly_equal(x, y, ctx));
    fq_nmod_mpoly_clear(x, ctx); fq_nmod_mpoly_clear(y, ctx); fq_nmod_mpoly_ctx_clear(ctx);
}

static void check_projection(const fq_nmod_ctx_t ctx, flint_rand_t state)
{
    const slong n = 3;
    unified_mpoly_struct polys[4] = {0}, **matrix, **modified, full = {0};
    fq_nmod_t c; fq_nmod_init(c, ctx);
    for (slong k = 0; k <= n; k++) {
        dr_mpoly_init(polys+k, n, 1, ctx);
        for (slong i = -1; i <= n; i++) for (slong j = i; j <= n; j++) {
            slong exp[4] = {0};
            if (i >= 0) exp[i]++;
            if (j >= 0) exp[j]++;
            fq_nmod_randtest_not_zero(c, state, ctx);
            dr_mpoly_add_term_fast(polys+k, exp, exp+n, c);
        }
    }
    build_fq_cancellation_matrix(&matrix, polys, n, 1);
    perform_fq_matrix_row_operations(&modified, &matrix, n, 1);
    compute_fq_cancel_matrix_det(&full, modified, n, 1, DET_METHOD_RECURSIVE);
    slong rows[15], cols[21];
    for (slong i = 0; i < 7; i++) {
        DR_MPOLY_TERM(term, &full, i*dr_mpoly_length(&full)/7);
        if (i < 5) memcpy(rows+i*n, term.var_exp, n*sizeof(slong));
        memcpy(cols+i*n, term.var_exp+n, n*sizeof(slong));
    }
    for (int indexed = 0; indexed <= 1; indexed++)
    for (int threads = 1; threads <= 4; threads += 3) {
        g_dixon_mq_step1_rank = indexed;
        omp_set_num_threads(threads);
        unified_mpoly_struct actual = {0}, expected = {0};
        assert(compute_fq_det_mq_projected_rect(&actual, modified, n+1, rows, 5, cols, 7));
        assert(fq_mq_project_full(&expected, &full, rows, 5, cols, 7));
        equal_poly(&actual, &expected);
        dr_mpoly_clear(&actual); dr_mpoly_clear(&expected);
    }
    /* Empty target block and exact cancellation of two equal columns. */
    slong absent[3] = {4,4,4};
    unified_mpoly_struct zero = {0};
    assert(compute_fq_det_mq_projected_rect(&zero, modified, n+1, absent, 1, cols, 7));
    assert(dr_mpoly_length(&zero) == 0); dr_mpoly_clear(&zero);
    unified_mpoly_struct projected = {0};
    int verified = dixon_try_mq_projection(&projected, modified, polys, n, 1, DET_METHOD_RECURSIVE);
    assert(verified);
    dr_mpoly_clear(&projected);
    slong saved_limit = g_dixon_det_cache_limit;
    g_dixon_det_cache_limit = 1;
    assert(!compute_fq_det_mq_projected_rect(&projected, modified, n+1, rows, 5, cols, 7));
    assert(!projected.ring);
    g_dixon_det_cache_limit = saved_limit;
    int saved_shared = g_dixon_mq_step1_shared;
    g_dixon_mq_step1_shared = 0;
    assert(!compute_fq_det_mq_projected_rect(&projected, modified, n+1, rows, 5, cols, 7));
    g_dixon_mq_step1_shared = saved_shared;
    g_field_equation_reduction = 1;
    assert(!compute_fq_det_mq_projected_rect(&projected, modified, n+1, rows, 5, cols, 7));
    g_field_equation_reduction = 0;
    for (slong i = 0; i <= n; i++) dr_mpoly_copy(&modified[i][1], &modified[i][0]);
    assert(compute_fq_det_mq_projected_rect(&zero, modified, n+1, rows, 5, cols, 7));
    assert(dr_mpoly_length(&zero) == 0); dr_mpoly_clear(&zero);
    assert(!dixon_try_mq_projection(&projected, modified, polys, n, 1, DET_METHOD_RECURSIVE));
    assert(!projected.ring);
    for (slong j = 0; j <= n; j++) unified_mpoly_zero(&modified[n][j]);
    assert(compute_fq_det_mq_projected_rect(&zero, modified, n+1, rows, 5, cols, 7));
    assert(dr_mpoly_length(&zero) == 0); dr_mpoly_clear(&zero);
    for (slong i = 0; i <= n; i++) {
        dr_mpoly_clear(polys+i);
        for (slong j = 0; j <= n; j++) {
            dr_mpoly_clear(&matrix[i][j]); dr_mpoly_clear(&modified[i][j]);
        }
        flint_free(matrix[i]); flint_free(modified[i]);
    }
    flint_free(matrix); flint_free(modified); dr_mpoly_clear(&full); fq_nmod_clear(c, ctx);
}

/* Independent signed permutation expansion, also valid in characteristic two. */
static void expand(fq_nmod_poly_t out, const fq_nmod_poly_mat_t m, slong row,
    ulong used, const fq_nmod_poly_t product, const fq_nmod_ctx_t ctx)
{
    if (row == m->r) { fq_nmod_poly_add(out, out, product, ctx); return; }
    fq_nmod_poly_t tmp; fq_nmod_poly_init(tmp, ctx);
    slong position = 0;
    for (slong c = 0; c < m->c; c++) if (!(used & (UWORD(1)<<c))) {
        fq_nmod_poly_mul(tmp, product, fq_nmod_poly_mat_entry(m,row,c), ctx);
        if (position++ & 1) fq_nmod_poly_neg(tmp,tmp,ctx);
        expand(out, m, row+1, used | (UWORD(1)<<c), tmp, ctx);
    }
    fq_nmod_poly_clear(tmp, ctx);
}

static void copy_matrix(fq_nmod_poly_mat_t a, const fq_nmod_poly_mat_t b, const fq_nmod_ctx_t ctx)
{
    for (slong i = 0; i < b->r; i++) for (slong j = 0; j < b->c; j++)
        fq_nmod_poly_set(fq_nmod_poly_mat_entry(a,i,j), fq_nmod_poly_mat_entry(b,i,j), ctx);
}
static int equal_matrix(const fq_nmod_poly_mat_t a, const fq_nmod_poly_mat_t b, const fq_nmod_ctx_t ctx)
{
    for (slong i = 0; i < b->r; i++) for (slong j = 0; j < b->c; j++)
        if (!fq_nmod_poly_equal(fq_nmod_poly_mat_entry(a,i,j), fq_nmod_poly_mat_entry(b,i,j), ctx)) return 0;
    return 1;
}

static void check_schur(const fq_nmod_ctx_t ctx, flint_rand_t state)
{
    slong rows[] = {1,0,3,4,2}, cols[] = {4,3,2,1,0};
    slong rd[] = {1,0,1,2,2}, cd[] = {0,1,3,2,2};
    dixon_mq_step4_profile p = {5,2,4,rows,cols,rd,cd,1};
    fq_nmod_poly_mat_t mat, copy;
    fq_nmod_poly_mat_init(mat, 5, 5, ctx); fq_nmod_poly_mat_init(copy, 5, 5, ctx);
    fq_nmod_poly_t actual, expected, one;
    fq_nmod_poly_init(actual, ctx); fq_nmod_poly_init(expected, ctx);
    fq_nmod_poly_init(one, ctx); fq_nmod_poly_one(one, ctx);
    fq_nmod_t constant; fq_nmod_init(constant,ctx);
    for (int trial = 0; trial < 6; trial++) {
        /* Vary external permutation parity and the complement sort parity
         * independently of pivot swaps. Nonunit pivots exercise scaling. */
        rows[0] = (trial & 1) ? 0 : 1; rows[1] = (trial & 1) ? 1 : 0;
        rd[2] = trial >= 2 ? 2 : 1; rd[3] = trial >= 2 ? 1 : 2;
        cd[2] = trial >= 2 ? 2 : 3; cd[3] = trial >= 2 ? 3 : 2;
        if (trial == 4) { slong swap = cd[3]; cd[3] = cd[4]; cd[4] = swap; }
        else cd[4] = 2;
        p.odd = dixon_mq_step4_permutation_odd(rows,5) ^ dixon_mq_step4_permutation_odd(cols,5);
        for (slong i = 0; i < 5; i++) for (slong j = 0; j < 5; j++) {
            fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(mat,rows[i],cols[j]);
            fq_nmod_poly_randtest(entry, state, FLINT_MAX(0,5-rd[i]-cd[j]), ctx);
            if (i >= 2 && j >= 2 && rd[i]+cd[j] == 4) {
                fq_nmod_poly_zero(entry, ctx);
                slong ri = 0, ci = 0, count = 0;
                for (slong k = 2; k < 5; k++) {
                    if (rd[k] == rd[i]) { if (k < i) ri++; count++; }
                    if (cd[k] == cd[j] && k < j) ci++;
                }
                if (ri == ((trial & 1) ? count-1-ci : ci)) {
                    fq_nmod_randtest_not_zero(constant,state,ctx);
                    fq_nmod_poly_set_coeff(entry,0,constant,ctx);
                }
            }
        }
        if (trial == 1) for (slong j = 0; j < 5; j++)
            fq_nmod_poly_zero(fq_nmod_poly_mat_entry(mat,rows[0],j), ctx);
        copy_matrix(copy, mat, ctx);
        fq_nmod_poly_zero(expected, ctx);
        expand(expected, mat, 0, 0, one, ctx);
        assert(dixon_mq_step4_try(actual, mat, &p, ctx));
        assert(fq_nmod_poly_equal(expected, actual, ctx));
        assert(equal_matrix(mat, copy, ctx));
        p.sigma = 0;
        assert(!dixon_mq_step4_try(actual, mat, &p, ctx));
        assert(fq_nmod_poly_equal(expected, actual, ctx));
        assert(equal_matrix(mat, copy, ctx)); p.sigma = 4;
        for (slong j = 2; j < 5; j++) fq_nmod_poly_zero(fq_nmod_poly_mat_entry(mat,rows[4],cols[j]), ctx);
        copy_matrix(copy, mat, ctx);
        assert(!dixon_mq_step4_try(actual, mat, &p, ctx));
        assert(fq_nmod_poly_equal(expected, actual, ctx));
        assert(equal_matrix(mat, copy, ctx));
    }
    fq_nmod_clear(constant,ctx);
    fq_nmod_poly_clear(actual,ctx); fq_nmod_poly_clear(expected,ctx); fq_nmod_poly_clear(one,ctx);
    fq_nmod_poly_mat_clear(mat,ctx); fq_nmod_poly_mat_clear(copy,ctx);
}

int main(void)
{
    const slong degrees[] = {4,8,16,32,64,128,8,8};
    const uint64_t moduli[] = {0x13,0x11d,0x1002d,0x100008299,UINT64_C(0x247f43cb7),0x87,0x11b,0x11d};
    g_dixon_verbose_level = 0; g_dixon_det_cache_limit = 1024;
    g_dixon_mq_step1_filter = 1; g_dixon_mq_step1_shared = 1;
    flint_rand_t state; flint_rand_init(state); flint_rand_set_seed(state, 123, 456);
    for (int f = 0; f < 8; f++) {
        nmod_poly_t modulus; nmod_poly_init(modulus, 2);
        nmod_poly_set_coeff_ui(modulus,degrees[f],1);
        for (slong i = 0; i < 64 && i < degrees[f]; i++)
            if ((moduli[f]>>i)&1) nmod_poly_set_coeff_ui(modulus,i,1);
        fq_nmod_ctx_t ctx; fq_nmod_ctx_init_modulus(ctx,modulus,"a"); nmod_poly_clear(modulus);
        check_projection(ctx,state); check_schur(ctx,state);
        printf("MQ GF(2^%ld): exact projection, rank verification, Schur, fallback PASS\n",degrees[f]);
        fq_nmod_ctx_clear(ctx);
    }
    flint_rand_clear(state); flint_cleanup_master(); return 0;
}

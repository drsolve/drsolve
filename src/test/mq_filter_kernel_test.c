/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Exercise private packing/proof helpers without exposing test library APIs. */
#include "../determinant/fq_mpoly_mat_det.c"
#include <assert.h>

static void check_packing(flint_rand_t state, slong n, ordering_t order)
{
    slong rows[10 * 8], cols[10 * 8];
    for (slong i = 0; i < 8 * n; i++) {
        rows[i] = n_randint(state, 2);
        cols[i] = n_randint(state, 2);
    }
    rows[0] = cols[n - 1] = 3;
    mq_det_filter filter;
    assert(mq_filter_init(&filter, n, rows, cols, 8));
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, 2 * n + 1, order, 101);
    mq_filter_prepare_packed(&filter, ctx, n + 2);
    flint_bitcnt_t native_bits = mpoly_fix_bits(1 + FLINT_BIT_COUNT((ulong) n + 2), ctx->minfo);
    if (order == ORD_LEX && n <= (FLINT_BITS - 1) / native_bits)
        assert(filter.packed_bits != 0);
    nmod_mpoly_t input, actual, expected;
    nmod_mpoly_init(input, ctx); nmod_mpoly_init(actual, ctx); nmod_mpoly_init(expected, ctx);
    ulong exp[21];
    for (slong t = 0; t < 256; t++) {
        slong r = n_randint(state, 8), c = n_randint(state, 8);
        for (slong v = 0; v < n; v++) {
            exp[v] = (t % 3 == 0) ? rows[r * n + v] : n_randint(state, rows[r * n + v] + 1);
            exp[n + v] = (t % 3 == 0) ? cols[c * n + v] : n_randint(state, cols[c * n + v] + 1);
        }
        if (t % 4 == 0) exp[t % n] = 4; /* Would alias if high coordinate bits were lost. */
        if (t % 4 == 1) exp[n + t % n] = 4;
        exp[2 * n] = t % 97; /* Parameter exponents must never enter the masks. */
        nmod_mpoly_push_term_ui_ui(input, 1 + t % 100, exp, ctx);
    }
    nmod_mpoly_sort_terms(input, ctx); nmod_mpoly_combine_like_terms(input, ctx);
    flint_bitcnt_t widths[] = {MPOLY_MIN_BITS, native_bits, 16, 32, FLINT_BITS, 2 * FLINT_BITS};
    for (int target = 0; target <= 1; target++) {
        nmod_mpoly_zero(expected, ctx);
        for (slong t = 0; t < input->length; t++) {
            nmod_mpoly_get_term_exp_ui(exp, input, t, ctx);
            if (mq_filter_accepts(&filter, exp, target))
                nmod_mpoly_push_term_ui_ui(expected, nmod_mpoly_get_term_coeff_ui(input, t, ctx), exp, ctx);
        }
        assert(expected->length > 0 && expected->length < input->length);
        for (slong k = 0; k < 6; k++) {
            assert(nmod_mpoly_repack_bits(actual, input, widths[k], ctx));
            mq_filter_poly(actual, ctx, &filter, target);
            assert(nmod_mpoly_is_canonical(actual, ctx));
            assert(nmod_mpoly_equal(actual, expected, ctx));
        }
        /* Also exercise an empty polynomial without reading exponent data. */
        nmod_mpoly_zero(actual, ctx);
        mq_filter_poly(actual, ctx, &filter, target);
        assert(nmod_mpoly_is_zero(actual, ctx));
    }
    nmod_mpoly_clear(input, ctx); nmod_mpoly_clear(actual, ctx); nmod_mpoly_clear(expected, ctx);
    nmod_mpoly_ctx_clear(ctx); mq_filter_clear(&filter);
}

static void check_layer_certificate(void)
{
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx, 101, 1, "a");
    fq_mvpoly_t **m = flint_malloc(4 * sizeof(*m));
    fq_nmod_t one;
    fq_nmod_init(one, ctx); fq_nmod_one(one, ctx);
    slong exp[6] = {0}, par = 0;
    for (slong row = 0; row < 4; row++) {
        m[row] = flint_malloc(4 * sizeof(**m));
        for (slong col = 0; col < 4; col++) {
            fq_mvpoly_init(&m[row][col], 6, 1, ctx);
            fq_mvpoly_add_term_fast(&m[row][col], exp, &par, one);
            if (row == 0) continue;
            slong v = 3 - row;
            exp[v] = 1;
            fq_mvpoly_add_term_fast(&m[row][col], exp, &par, one); exp[v] = 0;
            exp[3 + v] = 1;
            fq_mvpoly_add_term_fast(&m[row][col], exp, &par, one); exp[3 + v] = 0;
        }
    }
    /* Closure includes x0*x1, but not x0*x2. Exactly the last two rows are safe. */
    slong two[] = {2,0,0, 1,1,0, 0,2,0, 0,0,1};
    /* All pure squares are present, but the mixed term is missing. */
    slong one_layer[] = {2,0,0, 0,2,0, 0,0,2};
    slong zero[] = {0,0,0};
    const slong *sets[] = {zero, one_layer, two};
    slong sizes[] = {1,3,4};
    for (slong i = 0; i < 3; i++) {
        mq_det_filter f;
        assert(mq_filter_init(&f, 3, sets[i], sets[i], sizes[i]));
        assert(mq_safe_axis_layers(&f, &f.rows, m, 4, 0) == i);
        assert(mq_safe_axis_layers(&f, &f.cols, m, 4, 1) == i);
        mq_filter_clear(&f);
    }
    /* Full degree-three simplex certifies all LINEAR rows, never row zero. */
    slong all[30], count = 0;
    for (slong a = 0; a <= 3; a++) for (slong b = 0; b <= 3 - a; b++) {
        all[3 * count] = a; all[3 * count + 1] = b;
        all[3 * count++ + 2] = 3 - a - b;
    }
    mq_det_filter f;
    assert(mq_filter_init(&f, 3, all, all, count));
    assert(mq_safe_axis_layers(&f, &f.rows, m, 4, 0) == 3);
    mq_filter_clear(&f);
    for (slong row = 0; row < 4; row++) {
        for (slong col = 0; col < 4; col++) fq_mvpoly_clear(&m[row][col]);
        flint_free(m[row]);
    }
    flint_free(m); fq_nmod_clear(one, ctx); fq_nmod_ctx_clear(ctx);
}

int main(void)
{
    flint_rand_t state;
    flint_rand_init(state); flint_rand_set_seed(state, 995, 192);
    slong axes[] = {1,3,6,7,8,10};
    for (int order = ORD_LEX; order <= ORD_DEGREVLEX; order++)
        for (slong i = 0; i < 6; i++) check_packing(state, axes[i], order);
    check_layer_certificate();
    flint_rand_clear(state); flint_cleanup_master();
    puts("MQ packed filtering and safe-layer certificates passed");
    return 0;
}

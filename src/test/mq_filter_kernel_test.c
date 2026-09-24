/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Exercise private packing/proof helpers without exposing test library APIs. */
#define DRSOLVE_MQ_SUM_TEST 1
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
    assert(mq_filter_init(&filter, n, rows, 8, cols, 8));
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
        assert(mq_filter_init(&f, 3, sets[i], sizes[i], sets[i], sizes[i]));
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
    assert(mq_filter_init(&f, 3, all, count, all, count));
    assert(mq_safe_axis_layers(&f, &f.rows, m, 4, 0) == 3);
    mq_filter_clear(&f);
    for (slong row = 0; row < 4; row++) {
        for (slong col = 0; col < 4; col++) fq_mvpoly_clear(&m[row][col]);
        flint_free(m[row]);
    }
    flint_free(m); fq_nmod_clear(one, ctx); fq_nmod_ctx_clear(ctx);
}

static void check_linear_kernel(flint_rand_t state, slong n, ulong q)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx,2*n+1,ORD_LEX,q);
    slong target[16]; for(slong v=0;v<n;v++)target[v]=1;
    mq_det_filter f; assert(mq_filter_init(&f,n,target,1,target,1));
    mq_filter_prepare_packed(&f,ctx,n+2);
    assert(f.arithmetic_bits);
    nmod_mpoly_t a,b,x,y;
    nmod_mpoly_init(a,ctx);nmod_mpoly_init(b,ctx);nmod_mpoly_init(x,ctx);nmod_mpoly_init(y,ctx);
    ulong exp[33]={0};
    nmod_mpoly_push_term_ui_ui(a,1,exp,ctx);
    for(slong v=0;v<2*n+1;v++) {
        exp[v]=1;nmod_mpoly_push_term_ui_ui(a,1,exp,ctx);exp[v]=0;
    }
    nmod_mpoly_sort_terms(a,ctx);nmod_mpoly_combine_like_terms(a,ctx);
    nmod_mpoly_randtest_bound(b,state,2000,3,ctx);
    assert(nmod_mpoly_repack_bits_inplace(a,f.arithmetic_bits,ctx));
    assert(nmod_mpoly_repack_bits_inplace(b,f.arithmetic_bits,ctx));
    for(int truncate=0;truncate<2;truncate++) {
        nmod_mpoly_mul(y,a,b,ctx);if(truncate)mq_filter_poly(y,ctx,&f,0);
        assert(mq_linear_mul(x,a,b,ctx,&f,truncate));
        assert(nmod_mpoly_equal(x,y,ctx));
        assert(nmod_mpoly_is_canonical(x,ctx));
    }
    nmod_mpoly_t extra, product;
    nmod_mpoly_init(extra,ctx);nmod_mpoly_init(product,ctx);
    nmod_mpoly_randtest_bound(extra,state,1800,3,ctx);
    assert(nmod_mpoly_repack_bits_inplace(extra,f.arithmetic_bits,ctx));
    nmod_mpoly_struct *aa[] = {a,a,a,a};
    nmod_mpoly_struct *bb[] = {b,extra,extra,b};
    for(mq_sum_experiment_mode=1;mq_sum_experiment_mode<=2;mq_sum_experiment_mode++)
    for(slong count=2;count<=4;count++) for(int truncate=0;truncate<2;truncate++) {
        nmod_mpoly_zero(y,ctx);
        for(slong j=0;j<count;j++) {
            nmod_mpoly_mul(product,aa[j],bb[j],ctx);
            if(j&1)nmod_mpoly_sub(y,y,product,ctx);
            else nmod_mpoly_add(y,y,product,ctx);
        }
        if(truncate)mq_filter_poly(y,ctx,&f,0);
        assert(mq_linear_sum(x,aa,bb,count,ctx,&f,truncate));
        assert(nmod_mpoly_equal(x,y,ctx));
        assert(nmod_mpoly_is_canonical(x,ctx));
        if(count==4)assert(nmod_mpoly_is_zero(x,ctx));
    }
    assert(!mq_linear_sum(a,aa,bb,4,ctx,&f,0));
    nmod_mpoly_clear(extra,ctx);nmod_mpoly_clear(product,ctx);
    nmod_mpoly_clear(a,ctx);nmod_mpoly_clear(b,ctx);nmod_mpoly_clear(x,ctx);nmod_mpoly_clear(y,ctx);
    mq_filter_clear(&f);nmod_mpoly_ctx_clear(ctx);
}

static void check_shared_shards(void)
{
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx, 19, ORD_LEX, 65537);
    mq_shared_support prev, shifts;
    mq_shared_support_init(&prev, ctx, 16); /* Native three-word fallback. */
    assert(prev.words == 3);
    mq_shared_empty(&shifts, prev.bits, prev.words);
    ulong exp[64] = {0}, key[3];
    for (unsigned mask = 0; mask < 65536; mask++) {
        for (slong v = 0; v < 16; v++) exp[v] = (mask >> v) & 1;
        mpoly_set_monomial_ui(key, exp, prev.bits, ctx->minfo);
        mq_shared_insert(&prev, key);
    }
    memset(exp, 0, sizeof(exp)); memset(key, 0, sizeof(key)); mq_shared_insert(&shifts, key);
    for (slong v = 0; v < 16; v++) {
        exp[v] = 1; mpoly_set_monomial_ui(key, exp, prev.bits, ctx->minfo);
        mq_shared_insert(&shifts, key); exp[v] = 0;
    }
    slong target[9]; for (slong v = 0; v < 9; v++) target[v] = 1;
    mq_det_filter filter; assert(mq_filter_init(&filter, 9, target, 1, target, 1));
    mq_filter_prepare_packed(&filter, ctx, 16);
    int saved = omp_get_max_threads();
    for (int projected = 0; projected < 2; projected++) {
        const mq_det_filter *f = projected ? &filter : NULL;
        mq_shared_support reference;
        uint32_t *expected = mq_shared_next(&reference, &prev, &shifts, f, ctx, 1, 0);
        int threads[] = {3, 4, 16};
        for (unsigned t = 0; t < 3; t++) {
            omp_set_num_threads(threads[t]);
            mq_shared_support actual;
            uint32_t *map = mq_shared_next(&actual, &prev, &shifts, f, ctx, 1, 1);
            assert(actual.count == reference.count && actual.table == NULL);
            for (size_t i = 0; i < prev.count*shifts.count; i++) {
                assert((map[i] == UINT32_MAX) == (expected[i] == UINT32_MAX));
                if (map[i] != UINT32_MAX)
                    assert(!memcmp(actual.keys+(size_t)map[i]*actual.words,
                                   reference.keys+(size_t)expected[i]*reference.words,
                                   actual.words*sizeof(ulong)));
            }
            flint_free(map); mq_shared_support_clear(&actual);
        }
        flint_free(expected); mq_shared_support_clear(&reference);
    }
    omp_set_num_threads(saved);
    mq_filter_clear(&filter); mq_shared_support_clear(&prev); mq_shared_support_clear(&shifts);
    nmod_mpoly_ctx_clear(ctx);
}

static void check_shared_admission(void)
{
    for (slong n = 8; n <= 9; n++) {
        nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx, 2*n-1, ORD_LEX, 65537);
        nmod_mpoly_t **a = flint_malloc(n*sizeof(*a));
        ulong choose[FLINT_BITS][FLINT_BITS] = {{0}}, e[FLINT_BITS] = {0};
        for (slong i = 0; i <= n; i++) {
            choose[i][0] = 1;
            for (slong j = 1; j <= i; j++) choose[i][j] = choose[i-1][j-1]+choose[i-1][j];
        }
        for (slong r = 0; r < n; r++) {
            a[r] = flint_malloc(n*sizeof(**a));
            for (slong c = 0; c < n; c++) {
                nmod_mpoly_init(a[r][c], ctx);
                for (slong v = 0; v < 2*n-1; v++) {
                    e[v] = r == 0 ? 2 : 1;
                    nmod_mpoly_push_term_ui_ui(a[r][c], 1+c+v, e, ctx); e[v] = 0;
                }
                nmod_mpoly_sort_terms(a[r][c], ctx);
            }
        }
        assert(mq_shared_admit(a, n, ctx, choose)); /* n=9 no longer budget-rejected. */
        if (n == 8) {
            nmod_mpoly_t out; nmod_mpoly_init(out, ctx); nmod_mpoly_one(out, ctx);
            e[0] = 3; nmod_mpoly_push_term_ui_ui(a[0][0], 1, e, ctx); e[0] = 0;
            nmod_mpoly_sort_terms(a[0][0], ctx);
            assert(!mq_shared_det(out, a, n, ctx, 0, NULL, choose));
            assert(nmod_mpoly_is_one(out, ctx)); /* Rejection preserves output. */
            for (slong r = 0; r < n; r++) for (slong c = 0; c < n; c++)
                nmod_mpoly_zero(a[r][c], ctx);
            assert(mq_shared_det(out, a, n, ctx, 0, NULL, choose));
            assert(nmod_mpoly_is_zero(out, ctx));
            nmod_mpoly_clear(out, ctx);
        }
        for (slong r = 0; r < n; r++) {
            for (slong c = 0; c < n; c++) nmod_mpoly_clear(a[r][c], ctx);
            flint_free(a[r]);
        }
        flint_free(a); nmod_mpoly_ctx_clear(ctx);
    }
}

int main(void)
{
    flint_rand_t state;
    flint_rand_init(state); flint_rand_set_seed(state, 995, 192);
    slong axes[] = {1,3,6,7,8,10};
    for (int order = ORD_LEX; order <= ORD_DEGREVLEX; order++)
        for (slong i = 0; i < 6; i++) check_packing(state, axes[i], order);
    check_layer_certificate();
    check_shared_admission();
    check_shared_shards();
    for(slong n=3;n<=9;n+=2) { check_linear_kernel(state,n,2);check_linear_kernel(state,n,65537); }
    flint_rand_clear(state); flint_cleanup_master();
    puts("MQ packed filtering and safe-layer certificates passed");
    return 0;
}

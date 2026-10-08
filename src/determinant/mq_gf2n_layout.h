/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Shared exponent plans, separate native characteristic-two arithmetic. */
#define MQGF_TYPE uint8_t
#define MQGF_NAME mq_shared_gf24
#define MQGF_MEMBER gf24
#define MQGF_ZERO(a) ((a) == 0)
#define MQGF_ADD(a,b) ((a) ^ (b))
#define MQGF_PREPARE(a) const uint8_t *scalar_row = gf24_get_scalar_row(a)
#define MQGF_MUL(a,b) scalar_row[b]
#include "mq_gf2n_template.h"

#define MQGF_TYPE uint8_t
#define MQGF_NAME mq_shared_gf28
#define MQGF_MEMBER gf28
#define MQGF_ZERO(a) ((a) == 0)
#define MQGF_ADD(a,b) ((a) ^ (b))
#define MQGF_PREPARE(a) const uint8_t *scalar_row = gf28_get_scalar_row(a)
#define MQGF_MUL(a,b) scalar_row[b]
#include "mq_gf2n_template.h"

#define MQGF_TYPE uint16_t
#define MQGF_NAME mq_shared_gf216
#define MQGF_MEMBER gf216
#define MQGF_ZERO(a) ((a) == 0)
#define MQGF_ADD(a,b) ((a) ^ (b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf216_mul(a,b)
#include "mq_gf2n_template.h"

#define MQGF_TYPE gf232_t
#define MQGF_NAME mq_shared_gf232
#define MQGF_MEMBER gf232
#define MQGF_ZERO(a) gf232_is_zero(&(a))
#define MQGF_ADD(a,b) gf232_add(&(a), &(b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf232_mul(&(a), &(b))
#include "mq_gf2n_template.h"

#define MQGF_TYPE gf264_t
#define MQGF_NAME mq_shared_gf264
#define MQGF_MEMBER gf264
#define MQGF_ZERO(a) gf264_is_zero(&(a))
#define MQGF_ADD(a,b) gf264_add(&(a), &(b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf264_mul(&(a), &(b))
#include "mq_gf2n_template.h"

#define MQGF_TYPE gf2128_t
#define MQGF_NAME mq_shared_gf2128
#define MQGF_MEMBER gf2128
#define MQGF_ZERO(a) gf2128_is_zero(&(a))
#define MQGF_ADD(a,b) gf2128_add(&(a), &(b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf2128_mul(&(a), &(b))
#include "mq_gf2n_template.h"

#include "mq_fq_layout.h"

static int mq_extension_projected(unified_mpoly_struct *result,
    unified_mpoly_struct **matrix, slong n, mq_det_filter *filter)
{
    if (FLINT_BITS != 64 || n < 4 || 2*n-1 >= FLINT_BITS ||
        !g_dixon_mq_step1_shared || g_field_equation_reduction ||
        g_dixon_det_cache_limit <= 0) return 0;
    field_ctx_t field;
    field_ctx_init(&field, matrix[0][0].ctx);
    int native = field.field_id >= FIELD_ID_GF24 && field.field_id <= FIELD_ID_GF2128;
    slong degree = fq_nmod_ctx_degree(matrix[0][0].ctx);
    if (degree <= 1 || (size_t)degree > SIZE_MAX/sizeof(ulong)) {
        field_ctx_clear(&field); return 0;
    }
    size_t coefficient_bytes = native ? field.elem_size : (size_t)degree*sizeof(ulong);
    ulong choose[FLINT_BITS][FLINT_BITS] = {{0}};
    ulong limit = g_dixon_det_cache_limit;
    int ok = 1;
    for (slong i = 0; i <= n; i++) {
        choose[i][0] = 1;
        for (slong j = 1; j <= i; j++) {
            ulong a = choose[i-1][j-1], b = choose[i-1][j];
            choose[i][j] = a > limit || b > limit-a ? limit+1 : a+b;
        }
    }
    for (slong k = 1; k <= n; k++) {
        ulong count = k == n ? n : choose[n][k], prev = k == 1 ? 0 : choose[n][k-1];
        if (count > limit || prev > limit-count) ok = 0;
    }
    if (!ok) { field_ctx_clear(&field); return 0; }
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, 2*n-1, ORD_LEX, 2);
    nmod_mpoly_t **supports = flint_malloc(n*sizeof(*supports));
    for (slong i = 0; i < n; i++) {
        supports[i] = flint_malloc(n*sizeof(**supports));
        for (slong j = 0; j < n; j++) {
            nmod_mpoly_init(supports[i][j], ctx);
            const unified_mpoly_struct *p = &matrix[i][j];
            for (slong t = 0; t < dr_mpoly_length(p); t++) {
                DR_MPOLY_TERM(term, p, t);
                ulong exp[FLINT_BITS];
                for (slong v = 0; v < 2*n-2; v++) exp[v] = term.var_exp[v];
                exp[2*n-2] = term.par_exp[0];
                nmod_mpoly_set_coeff_ui_ui(supports[i][j], 1, exp, ctx);
            }
        }
    }
    ok = mq_shared_admit_coeffs(supports, n, ctx, choose, coefficient_bytes);
    if (ok) {
        filter->safe_linear_layers = FLINT_MIN(
            mq_safe_axis_layers(filter, &filter->rows, matrix, n, 0),
            mq_safe_axis_layers(filter, &filter->cols, matrix, n, 1));
        mq_filter_prepare_packed(filter, ctx, n+1);
        if (g_dixon_verbose_level >= 2)
            printf("  MQ %s shared-index DP: %zu-byte coefficients\n",
                   native ? "native GF(2^n)" : "fq_nmod", coefficient_bytes);
        int parallel = n >= PARALLEL_THRESHOLD && omp_get_max_threads() > 1;
        switch (field.field_id) {
        case FIELD_ID_GF24: ok = mq_shared_gf24(result, matrix, supports, n, ctx, parallel, filter, choose, &field); break;
        case FIELD_ID_GF28: ok = mq_shared_gf28(result, matrix, supports, n, ctx, parallel, filter, choose, &field); break;
        case FIELD_ID_GF216: ok = mq_shared_gf216(result, matrix, supports, n, ctx, parallel, filter, choose, &field); break;
        case FIELD_ID_GF232: ok = mq_shared_gf232(result, matrix, supports, n, ctx, parallel, filter, choose, &field); break;
        case FIELD_ID_GF264: ok = mq_shared_gf264(result, matrix, supports, n, ctx, parallel, filter, choose, &field); break;
        case FIELD_ID_GF2128: ok = mq_shared_gf2128(result, matrix, supports, n, ctx, parallel, filter, choose, &field); break;
        default: ok = mq_shared_fq(result, matrix, supports, n, ctx, parallel, filter, choose, matrix[0][0].ctx); break;
        }
    }
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) nmod_mpoly_clear(supports[i][j], ctx);
        flint_free(supports[i]);
    }
    flint_free(supports); nmod_mpoly_ctx_clear(ctx); field_ctx_clear(&field);
    return ok;
}

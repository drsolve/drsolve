/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Shared monomial plans with contiguous fq_nmod polynomial-basis coefficients. */
#include "mq_fq_coeffs.h"
static int mq_shared_fq(unified_mpoly_struct *result, unified_mpoly_struct **matrix,
    nmod_mpoly_t **supports, slong n, const nmod_mpoly_ctx_t ctx,
    int parallel, const mq_det_filter *filter, ulong choose[FLINT_BITS][FLINT_BITS],
    const fq_nmod_ctx_t fq)
{
    slong d = fq_nmod_ctx_degree(fq);
    nmod_t mod = fq_nmod_ctx_modulus(fq)->mod;
    mq_shared_support prev;
    mq_shared_support_init(&prev, ctx, n+1);
    ulong *values = flint_calloc(d,sizeof(ulong)); values[0] = 1;
    for (slong k = 1; k <= n; k++) {
        const mq_det_filter *f = k > filter->safe_linear_layers ? filter : NULL;
        slong count = k == n ? n : (slong)choose[n][k];
        mq_shared_support shifts, next;
        mq_shared_row(&shifts, supports[n-k], n, &prev, ctx);
        uint32_t *map = NULL;
        size_t previous_count = k == 1 ? 1 : choose[n][k-1];
        int ranked = mq_rank_enabled() && mq_rank_next(&next, &map, &prev, &shifts,
            filter, ctx, n, k, count, previous_count, parallel);
        if (!ranked) map = mq_shared_next(&next, &prev, &shifts, f, ctx, 1, parallel);
        if (g_dixon_verbose_level >= 2)
            printf("  MQ fq_nmod %s layer %ld: %zu monomials\n",
                   ranked ? "direct-index" : "shared-index", k, next.count);
        if (!next.count) {
            flint_free(map); flint_free(values);
            mq_shared_support_clear(&prev); mq_shared_support_clear(&shifts);
            mq_shared_support_clear(&next);
            dr_mpoly_init(result, 2*n-2, 1, matrix[0][0].ctx);
            return 1;
        }
        ulong *factors = flint_calloc(n*shifts.count*d, sizeof(ulong));
        for (slong c = 0; c < n; c++) {
            const unified_mpoly_struct *p = &matrix[n-k][c];
            for (slong t = 0; t < dr_mpoly_length(p); t++) {
                DR_MPOLY_TERM(term, p, t);
                ulong exp[FLINT_BITS], key[3];
                for (slong v = 0; v < 2*n-2; v++) exp[v] = term.var_exp[v];
                exp[2*n-2] = term.par_exp[0];
                mpoly_set_monomial_ui(key, exp, prev.bits, ctx->minfo);
                uint32_t id = mq_shared_find(&shifts, key);
                FLINT_ASSERT(id != UINT32_MAX);
                mq_fq_store(factors+(c*shifts.count+id)*d, term.coeff, d);
            }
        }
        ulong *output = flint_calloc((size_t)count*next.count*d, sizeof(ulong));
        #pragma omp parallel for if(parallel) schedule(static)
        for (slong index = 0; index < count; index++) {
            slong cols[FLINT_BITS], col = n-1;
            ulong rank = index, prefix[FLINT_BITS], suffix[FLINT_BITS];
            for (slong j = k == n ? 0 : k; j > 0; j--) {
                while (choose[col][j] > rank) col--;
                cols[j-1] = col; rank -= choose[col][j]; col--;
            }
            prefix[0] = 0;
            for (slong j = 0; k != n && j < k; j++) prefix[j+1] = prefix[j]+choose[cols[j]][j+1];
            suffix[k] = 0;
            for (slong j = k-1; k != n && j > 0; j--) suffix[j] = suffix[j+1]+choose[cols[j]][j];
            if (k == n) { cols[0] = index; prefix[0] = n-1-index; suffix[1] = 0; }
            ulong *dest = output+(size_t)index*next.count*d;
            fq_nmod_t scratch; fq_nmod_init(scratch,fq);
            for (slong j = 0; j < (k == n ? 1 : k); j++) {
                const ulong *source = values+(prefix[j]+suffix[j+1])*prev.count*d;
                for (size_t a = 0; a < shifts.count; a++) {
                    const ulong *scalar = factors+(cols[j]*shifts.count+a)*d;
                    if (!mq_fq_length(scalar,d)) continue;
                    const uint32_t *indices = map+a*prev.count;
                    mq_fq_operator op; mq_fq_operator_prepare(&op,scalar,scratch,fq);
                    int subtract = (k == n ? index : j) & 1;
                    for (size_t b = 0; b < prev.count; b++) if (indices[b] != UINT32_MAX)
                        mq_fq_addmul(dest+(size_t)indices[b]*d,&op,source+b*d,subtract,scratch,fq);
                }
            }
            fq_nmod_clear(scratch,fq);
        }
        if (k == n) {
            #pragma omp parallel for if(parallel) schedule(static)
            for (size_t i = 0; i < next.count*d; i++)
                for (slong j = 1; j < n; j++)
                    output[i] = nmod_add(output[i],output[j*next.count*d+i],mod);
        }
        flint_free(map); flint_free(factors); flint_free(values);
        mq_shared_support_clear(&prev); mq_shared_support_clear(&shifts);
        prev = next; values = output;
        if (k == n) {
            dr_mpoly_init(result, 2*n-2, 1, matrix[0][0].ctx);
            ulong exp[FLINT_BITS];
            for (size_t i = 0; i < prev.count; i++) if (mq_fq_length(values+i*d,d)) {
                mpoly_get_monomial_ui(exp, prev.keys+i*prev.words, prev.bits, ctx->minfo);
                if (!mq_filter_accepts(filter, exp, 1)) continue;
                fq_nmod_struct coefficient = mq_fq_view(values+i*d,fq);
                dr_mpoly_add_term_fast(result,(slong *)exp,(slong *)exp+2*n-2,&coefficient);
            }
            flint_free(values); mq_shared_support_clear(&prev);
            return 1;
        }
    }
    return 0;
}
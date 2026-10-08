/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Instantiated for each native binary field. The nmod objects passed here
 * contain ONLY monomial supports with coefficient one; all coefficient
 * arithmetic below is in MQGF_TYPE, never modulo the characteristic. */
static int MQGF_NAME(unified_mpoly_struct *result, unified_mpoly_struct **matrix,
    nmod_mpoly_t **supports, slong n, const nmod_mpoly_ctx_t ctx,
    int parallel, const mq_det_filter *filter, ulong choose[FLINT_BITS][FLINT_BITS],
    const field_ctx_t *field)
{
    mq_shared_support prev;
    mq_shared_support_init(&prev, ctx, n+1);
    MQGF_TYPE *values = flint_malloc(sizeof(MQGF_TYPE));
    field_elem_u one = {0};
    field_set_one(&one, field->field_id, field->ctx.fq_ctx);
    values[0] = one.MQGF_MEMBER;
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
            printf("  MQ GF(2^n) %s layer %ld: %zu monomials\n",
                   ranked ? "direct-index" : "shared-index", k, next.count);
        if (!next.count) {
            flint_free(map); flint_free(values);
            mq_shared_support_clear(&prev); mq_shared_support_clear(&shifts);
            mq_shared_support_clear(&next);
            dr_mpoly_init(result, 2*n-2, 1, matrix[0][0].ctx);
            return 1;
        }
        MQGF_TYPE *factors = flint_calloc(n*shifts.count, sizeof(MQGF_TYPE));
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
                field_elem_u coefficient;
                fq_nmod_to_field_elem(&coefficient, term.coeff, field);
                factors[c*shifts.count+id] = coefficient.MQGF_MEMBER;
            }
        }
        MQGF_TYPE *output = flint_calloc((size_t)count*next.count, sizeof(MQGF_TYPE));
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
            MQGF_TYPE *dest = output+(size_t)index*next.count;
            for (slong j = 0; j < (k == n ? 1 : k); j++) {
                const MQGF_TYPE *source = values+(prefix[j]+suffix[j+1])*prev.count;
                for (size_t a = 0; a < shifts.count; a++) {
                    MQGF_TYPE scalar = factors[cols[j]*shifts.count+a];
                    if (MQGF_ZERO(scalar)) continue;
                    const uint32_t *indices = map+a*prev.count;
                    MQGF_PREPARE(scalar);
                    for (size_t b = 0; b < prev.count; b++) if (indices[b] != UINT32_MAX) {
                        MQGF_TYPE product = MQGF_MUL(scalar, source[b]);
                        dest[indices[b]] = MQGF_ADD(dest[indices[b]], product);
                    }
                }
            }
        }
        if (k == n) {
            #pragma omp parallel for if(parallel) schedule(static)
            for (size_t i = 0; i < next.count; i++)
                for (slong j = 1; j < n; j++)
                    output[i] = MQGF_ADD(output[i], output[j*next.count+i]);
        }
        flint_free(map); flint_free(factors); flint_free(values);
        mq_shared_support_clear(&prev); mq_shared_support_clear(&shifts);
        prev = next; values = output;
        if (k == n) {
            dr_mpoly_init(result, 2*n-2, 1, matrix[0][0].ctx);
            fq_nmod_t coefficient; fq_nmod_init(coefficient, result->ctx);
            ulong exp[FLINT_BITS];
            for (size_t i = 0; i < prev.count; i++) if (!MQGF_ZERO(values[i])) {
                mpoly_get_monomial_ui(exp, prev.keys+i*prev.words, prev.bits, ctx->minfo);
                if (!mq_filter_accepts(filter, exp, 1)) continue;
                field_elem_u value; value.MQGF_MEMBER = values[i];
                field_elem_to_fq_nmod(coefficient, &value, field);
                dr_mpoly_add_term_fast(result, (slong *)exp, (slong *)exp+2*n-2, coefficient);
            }
            fq_nmod_clear(coefficient, result->ctx);
            flint_free(values); mq_shared_support_clear(&prev);
            return 1;
        }
    }
    return 0;
}
#undef MQGF_TYPE
#undef MQGF_NAME
#undef MQGF_MEMBER
#undef MQGF_ZERO
#undef MQGF_PREPARE
#undef MQGF_MUL
#undef MQGF_ADD

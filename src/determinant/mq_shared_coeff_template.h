/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Instantiated for limb and uint32_t coefficient arrays. The caller proves
 * admission and, for the narrow variant, that every unreduced sum fits.
 * No include guard: the two instantiations share the same support logic. */
static void MQ_SHARED_PACK(nmod_mpoly_t out, const MQ_SHARED_COEFF *coeffs,
                           const mq_shared_support *s, const nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(s->bits, ctx->minfo);
    slong words = mpoly_words_per_exp(bits, ctx->minfo);
    nmod_mpoly_fit_length_reset_bits(out, s->count, bits, ctx);
    slong len = 0;
    for (size_t i = 0; i < s->count; i++) if (coeffs[i]) {
        out->coeffs[len] = coeffs[i];
        memcpy(out->exps+len*words, s->keys+i*s->words, s->words*sizeof(ulong));
        len++;
    }
    _nmod_mpoly_set_length(out, len, ctx);
    if (!s->sorted) nmod_mpoly_sort_terms(out, ctx);
}

static int MQ_SHARED_DET(nmod_mpoly_t result, nmod_mpoly_t **matrix, slong n,
                          const nmod_mpoly_ctx_t ctx, int parallel,
                          const mq_det_filter *filter, ulong choose[FLINT_BITS][FLINT_BITS])
{
#ifdef DRSOLVE_MQ_LAYOUT_TEST
    mq_shared_test_used++;
#endif
    if (g_dixon_verbose_level >= 2) {
        if (DRSOLVE_MQ_SHARED_WORKSPACE_BYTES == SIZE_MAX)
            printf("  MQ shared-index DP: size=%ld, no fixed workspace cap\n", n);
        else
            printf("  MQ shared-index DP: size=%ld, workspace cap=%zu MiB\n", n,
                   (size_t)DRSOLVE_MQ_SHARED_WORKSPACE_BYTES/(1024*1024));
    }
    mq_shared_support prev; mq_shared_support_init(&prev, ctx, n+1);
    MQ_SHARED_COEFF *values = flint_malloc(sizeof(MQ_SHARED_COEFF)); values[0] = 1;
    for (slong k = 1; k <= n; k++) {
        const mq_det_filter *f = filter && k > filter->safe_linear_layers ? filter : NULL;
        slong count = k == n ? n : (slong)choose[n][k];
        mq_shared_support shifts, next;
        mq_shared_row(&shifts, matrix[n-k], n, &prev, ctx);
#ifdef DRSOLVE_MQ_LAYOUT_TEST
        double plan_start = omp_get_wtime();
#endif
        uint32_t *map = NULL;
        size_t previous_count = k == 1 ? 1 : choose[n][k-1];
        int ranked = mq_rank_enabled() && mq_rank_next(&next, &map, &prev, &shifts,
                         filter, ctx, n, k, count, previous_count, parallel);
        if (ranked && g_dixon_verbose_level >= 2)
            printf("  MQ direct-index layer %ld: %zu monomials\n", k, next.count);
        if (!ranked) map = mq_shared_next(&next, &prev, &shifts, f, ctx, 1, parallel);
#ifdef DRSOLVE_MQ_LAYOUT_TEST
        mq_shared_test_rank_layers += ranked;
#endif
#ifdef DRSOLVE_MQ_LAYOUT_TEST
        if (mq_shared_variant() & 4)
            mq_shared_order(&next, map, shifts.count*prev.count, parallel);
        mq_shared_test_plan_seconds += omp_get_wtime()-plan_start;
        double arithmetic_start = omp_get_wtime();
#endif
        if (!next.count) {
            flint_free(map); flint_free(values);
            mq_shared_support_clear(&prev); mq_shared_support_clear(&shifts);
            mq_shared_support_clear(&next); nmod_mpoly_zero(result, ctx);
            return 1;
        }
        ulong *factors = flint_calloc(n*shifts.count, sizeof(ulong));
        nmod_mpoly_t tmp; nmod_mpoly_init(tmp, ctx);
        for (slong c = 0; c < n; c++) {
            int packed = nmod_mpoly_repack_bits(tmp, matrix[n-k][c], mpoly_fix_bits(prev.bits, ctx->minfo), ctx);
            FLINT_ASSERT(packed); (void)packed;
            for (slong t = 0; t < tmp->length; t++) {
                const ulong *key = tmp->exps+t*prev.words;
                uint32_t id = mq_shared_find(&shifts, key);
                FLINT_ASSERT(id != UINT32_MAX); factors[c*shifts.count+id] = tmp->coeffs[t];
            }
        }
        nmod_mpoly_clear(tmp, ctx);
#ifdef DRSOLVE_MQ_DAG_TEST
        void *pending=NULL;
        if(k==n-1 && mq_dag_checkpoint(values,MQ_SHARED_NARROW,prev.count,
                factors,shifts.count,map,next.count,n,&pending)) {
            flint_free(map);mq_shared_support_clear(&prev);
            mq_shared_support_clear(&shifts);prev=next;values=pending;continue;
        }
        /* Keep the shared child DAG alive while the experimental consumer
         * solves its complement and adaptively contracts the pending root. */
        if(k==n && mq_dag_root(values,MQ_SHARED_NARROW,prev.count,factors,
                shifts.count,map,next.keys,next.count,next.words,next.bits,n,ctx)) {
            flint_free(map);flint_free(factors);flint_free(values);
            mq_shared_support_clear(&prev);mq_shared_support_clear(&shifts);
            mq_shared_support_clear(&next);nmod_mpoly_zero(result,ctx);return 1;
        }
#endif
        /* The root needs only the sum of its column contributions. Reuse
         * a bounded number of scratch rows instead of retaining all n roots.
         * Keep buffers <= count, preserving the admission memory bound. */
        slong lanes=count;
        if(k==n) {
            slong workers=1;
#ifdef _OPENMP
            if(parallel)workers=omp_get_max_threads();
#endif
            lanes=FLINT_MIN(count-1,FLINT_MIN(4,workers));
        }
        MQ_SHARED_COEFF *output = flint_calloc(
            (size_t)(k==n ? lanes+1 : count)*next.count,sizeof(MQ_SHARED_COEFF));
        /* Each shift is injective: one output receives at most k*#shifts
         * products, not #child-monomials times that number. Negative factors
         * are represented by their canonical residues. Division checks avoid
         * overflow in the proof itself, including full-word primes. */
        ulong contributions = (k == n ? 1 : k)*shifts.count;
        ulong largest = ctx->mod.n-1;
        int lazy = MQ_SHARED_NARROW || (contributions &&
                   largest <= (UWORD_MAX/contributions)/largest);
        for(slong batch=0;batch<count;batch+=lanes) {
            slong end=FLINT_MIN(count,batch+lanes);
            #pragma omp parallel for if(parallel) schedule(static)
            for (slong index = batch; index < end; index++) {
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
                MQ_SHARED_COEFF *dest = output+(size_t)(k==n ? index-batch+1 : index)*next.count;
                if(k==n)memset(dest,0,next.count*sizeof(*dest));
                for (slong j = 0; j < (k == n ? 1 : k); j++) {
                    const MQ_SHARED_COEFF *source = values+(prefix[j]+suffix[j+1])*prev.count;
                    for (size_t a = 0; a < shifts.count; a++) {
                        MQ_SHARED_COEFF scalar = factors[cols[j]*shifts.count+a];
                        if ((k == n ? index : j) & 1) scalar = nmod_neg(scalar, ctx->mod);
                        if (!scalar) continue;
                        const uint32_t *indices = map+a*prev.count;
                        if (lazy) {
                            for (size_t b = 0; b < prev.count; b++) if (indices[b] != UINT32_MAX)
                                dest[indices[b]] += scalar*source[b];
                        } else {
                            for (size_t b = 0; b < prev.count; b++) if (indices[b] != UINT32_MAX) {
                                ulong product = nmod_mul(scalar, source[b], ctx->mod);
                                dest[indices[b]] = nmod_add(dest[indices[b]], product, ctx->mod);
                            }
                        }
                    }
                }
                if (lazy) for (size_t i = 0; i < next.count; i++) {
                    /* FLINT's reduction macros require limb-width operands even
                     * when the accumulator storage is narrower. */
                    ulong z = dest[i];
                    NMOD_RED(z, z, ctx->mod);
                    dest[i] = z;
                }
            }
            if(k==n) {
                #pragma omp parallel for if(parallel) schedule(static)
                for(size_t i=0;i<next.count;i++)for(slong j=1;j<=end-batch;j++)
                    output[i]=nmod_add(output[i],output[j*next.count+i],ctx->mod);
            }
        }
        if (k == n) {
            MQ_SHARED_PACK(result, output, &next, ctx);
        }
        flint_free(map); flint_free(factors); flint_free(values);
#ifdef DRSOLVE_MQ_LAYOUT_TEST
        mq_shared_test_arithmetic_seconds += omp_get_wtime()-arithmetic_start;
#endif
        mq_shared_support_clear(&prev); mq_shared_support_clear(&shifts);
        prev = next; values = output;
        if (k == n) { flint_free(values); mq_shared_support_clear(&prev); return 1; }
    }
    return 0;
}

#undef MQ_SHARED_COEFF
#undef MQ_SHARED_PACK
#undef MQ_SHARED_DET
#undef MQ_SHARED_NARROW

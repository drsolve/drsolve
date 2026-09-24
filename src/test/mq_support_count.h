/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Research-only Boolean support graph. No minor coefficients are computed.
 * Row unions deliberately relax column exclusivity and cancellation. */
int mq_support_count(fq_mvpoly_t **matrix, slong n, const slong *rows,
                     const slong *cols, slong targets)
{
    mq_det_filter filter;
    if (!mq_filter_init(&filter, n-1, rows, targets, cols, targets)) return 0;
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx, 2*n-1, ORD_LEX,
                        fq_nmod_ctx_modulus(matrix[0][0].ctx)->mod.n);
    mq_filter_prepare_packed(&filter, ctx, n+1);
    nmod_mpoly_t **m = flint_malloc(n*sizeof(*m));
    for (slong i=0; i<n; i++) m[i] = flint_malloc(n*sizeof(**m));
    fq_matrix_mvpoly_to_nmod_mpoly(m, matrix, n, ctx);
    mq_shared_support *support = flint_calloc(n+1, sizeof(*support));
    mq_shared_support *shifts = flint_calloc(n+1, sizeof(*shifts));
    unsigned char **live = flint_calloc(n+1, sizeof(*live));
    mq_shared_support_init(support, ctx, n+1);
    for (slong k=1; k<=n; k++) {
        mq_shared_row(shifts+k, m[n-k], n, support+k-1, ctx);
        /* Dense generic columns have the same support. This makes the
         * C(n,k)*k edge weighting exact for the support-level DP model. */
        for(slong col=0;col<n;col++) {
            int packed=nmod_mpoly_repack_bits_inplace(m[n-k][col],support[0].bits,ctx);
            FLINT_ASSERT(packed); (void)packed;
            FLINT_ASSERT((size_t)m[n-k][col]->length==shifts[k].count);
            for(slong t=0;t<m[n-k][col]->length;t++)
                FLINT_ASSERT(mq_shared_find(shifts+k,m[n-k][col]->exps+t*support[0].words)!=UINT32_MAX);
        }
        mq_shared_next(support+k, support+k-1, shifts+k, &filter, ctx, 0, 0);
        fprintf(stderr, "n=%ld forward layer=%ld support=%zu\n", n, k, support[k].count);
    }
    /* Export the exact support-only problem for an independent tuple oracle. */
    printf("{\"kind\":\"model\",\"n\":%ld,\"shifts\":[",n);
    ulong exported[FLINT_BITS];
    for(slong k=1;k<=n;k++) {
        printf("%s[",k==1?"":",");
        for(size_t a=0;a<shifts[k].count;a++) {
            mpoly_get_monomial_ui(exported,shifts[k].keys+a*shifts[k].words,
                                  shifts[k].bits,ctx->minfo);
            printf("%s[",a?",":"");
            for(slong v=0;v<2*n-1;v++) printf("%s%lu",v?",":"",exported[v]);
            printf("]");
        }
        printf("]");
    }
    const slong *axes[2]={rows,cols};
    for(int axis=0;axis<2;axis++) {
        printf(axis?"],\"cols\":[": "],\"rows\":[");
        for(slong i=0;i<targets;i++) {
            printf("%s[",i?",":"");
            for(slong v=0;v<n-1;v++) printf("%s%ld",v?",":"",axes[axis][i*(n-1)+v]);
            printf("]");
        }
    }
    printf("]}\n");
    for (slong k=0; k<=n; k++) live[k] = flint_calloc(support[k].count, 1);
    ulong exp[FLINT_BITS];
    size_t terminal = 0;
    for (size_t b=0; b<support[n].count; b++) {
        mpoly_get_monomial_ui(exp, support[n].keys+b*support[n].words,
                              support[n].bits, ctx->minfo);
        terminal += live[n][b] = mq_filter_accepts(&filter, exp, 1);
    }
    double baseline_work=0, retained_work=0;
    for (slong k=n; k>=1; k--) {
        mq_shared_support *prev=support+k-1, *next=support+k, *shift=shifts+k;
        size_t edges=0, kept_edges=0, kept=0;
        for (size_t b=0; b<next->count; b++) kept += live[k][b];
        for (size_t a=0; a<shift->count; a++) for (size_t b=0; b<prev->count; b++) {
            ulong key[3];
            for (slong w=0; w<prev->words; w++)
                key[w]=shift->keys[a*prev->words+w]+prev->keys[b*prev->words+w];
            uint32_t id=mq_shared_find(next,key);
            if (id==UINT32_MAX) continue;
            edges++;
            if (live[k][id]) { live[k-1][b]=1; kept_edges++; }
        }
        /* C(n,k)*k equally-used columns in the generic row-union model.
         * Root expansion has n cofactors, hence the same weight n. */
        size_t minors=mq_shared_binomial(n,k,SIZE_MAX/2), weight=minors*k;
        baseline_work+=(double)weight*edges;
        retained_work+=(double)weight*kept_edges;
        printf("{\"kind\":\"layer\",\"n\":%ld,\"k\":%ld,\"support\":%zu,\"live\":%zu,\"shifts\":%zu,\"edges\":%zu,\"live_edges\":%zu,\"weight\":%zu}\n",
               n,k,next->count,kept,shift->count,edges,kept_edges,weight);
    }
    printf("{\"kind\":\"summary\",\"n\":%ld,\"axis_targets\":%ld,\"terminal\":%zu,\"baseline_work\":%.0f,\"retained_work\":%.0f,\"work_ratio\":%.9f}\n",
           n,targets,terminal,baseline_work,retained_work,
           baseline_work?retained_work/baseline_work:0);
    for (slong k=0; k<=n; k++) {
        mq_shared_support_clear(support+k);
        if(k) mq_shared_support_clear(shifts+k);
        flint_free(live[k]);
    }
    flint_free(live); flint_free(support); flint_free(shifts);
    for (slong i=0; i<n; i++) {
        for(slong j=0; j<n; j++) nmod_mpoly_clear(m[i][j],ctx);
        flint_free(m[i]);
    }
    flint_free(m); mq_filter_clear(&filter); nmod_mpoly_ctx_clear(ctx);
    return 1;
}

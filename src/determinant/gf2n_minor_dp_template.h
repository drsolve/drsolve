/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Shared compact minor DP and characteristic-two merge. */
static void DP_NAME(swap)(DP_POLY *a, DP_POLY *b)
{
    DP_POLY t = *a; *a = *b; *b = t;
}

/* Nonaliasing merge into reusable output storage. */
static void DP_NAME(add)(DP_POLY *out,
    const DP_POLY *a, const DP_POLY *b,
    const DP_CTX_T ctx)
{
    FLINT_ASSERT(out != a && out != b);
    if (!a->length) { DP_SET(out, b, ctx); return; }
    if (!b->length) { DP_SET(out, a, ctx); return; }
    flint_bitcnt_t bits = FLINT_MAX(a->bits, b->bits);
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    ulong *ae = NULL, *be = NULL;
    const ulong *ax = a->exps, *bx = b->exps;
    if (a->bits != bits) {
        ae = flint_malloc(sizeof(ulong) * N * a->length);
        mpoly_repack_monomials(ae, bits, ax, a->bits, a->length, ctx->minfo);
        ax = ae;
    }
    if (b->bits != bits) {
        be = flint_malloc(sizeof(ulong) * N * b->length);
        mpoly_repack_monomials(be, bits, bx, b->bits, b->length, ctx->minfo);
        bx = be;
    }
    slong capacity = a->length + b->length;
    if (out->coeffs_alloc < capacity) {
        out->coeffs = flint_realloc(out->coeffs, capacity*sizeof(DP_E));
        out->coeffs_alloc = capacity;
    }
    if (out->exps_alloc < N * capacity) {
        out->exps = flint_realloc(out->exps, sizeof(ulong) * N * capacity);
        out->exps_alloc = N * capacity;
    }
    ulong *mask = flint_malloc(sizeof(ulong) * N);
    mpoly_get_cmpmask(mask, N, bits, ctx->minfo);
    slong i = 0, j = 0, len = 0;
    while (i < a->length || j < b->length) {
        int cmp = i == a->length ? -1 : j == b->length ? 1 :
            mpoly_monomial_cmp(ax + N*i, bx + N*j, N, mask);
        const ulong *exp = cmp >= 0 ? ax + N*i : bx + N*j;
        DP_E c;
        if (cmp > 0) c = a->coeffs[i++];
        else if (cmp < 0) c = b->coeffs[j++];
        else { c = DP_ADD(a->coeffs[i],b->coeffs[j]); ++i; ++j; }
        if (!DP_ZERO(c)) {
            out->coeffs[len] = c;
            memcpy(out->exps + N*len++, exp, sizeof(ulong) * N);
        }
    }
    out->length = len; out->bits = bits;
    flint_free(mask); flint_free(ae); flint_free(be);
}

/* Called only after the common DP entry-budget and sparsity preflight.
   Convert the small input matrix once; convert only the final determinant
   back to fq_nmod. All live minors retain packed coefficients. */
static int DP_NAME(det)(unified_mpoly_t result,
    unified_mpoly_t **matrix, slong size, unified_mpoly_ctx_t ctx,
    int use_parallel, const ulong choose[FLINT_BITS][FLINT_BITS])
{
    int failed=0;
    DP_CTX_T nc;
    DP_CTX_INIT(nc, ctx->nvars, ctx->ord);
    DP_POLY *entries = flint_calloc(size * size, sizeof(*entries));
    for (slong i = 0; i < size; ++i)
        for (slong j = 0; j < size; ++j) {
            DP_INIT(entries + i*size+j, nc);
            DP_IMPORT(entries + i*size+j, GET_FQ_POLY(matrix[i][j]), ctx, nc);
        }
    DP_POLY *previous = NULL;
    slong previous_count = 0;
    for (slong k = 1; k <= size; ++k) {
        slong count = k == size ? size : (slong)choose[size][k];
        DP_POLY *current = flint_malloc(count * sizeof(*current));
        for (slong i = 0; i < count; ++i) DP_INIT(current + i, nc);
        if (k == 1) {
            for (slong i = 0; i < count; ++i)
                DP_SET(current + i, entries + (size-1)*size+i, nc);
        } else {
#ifdef _OPENMP
            #pragma omp parallel if(use_parallel && count > 1 && !omp_in_parallel()) num_threads(FLINT_MIN(count, omp_get_max_threads())) reduction(|:failed)
#endif
            {
                DP_POLY_T product, sum;
                DP_INIT(product, nc); DP_INIT(sum, nc);
#ifdef _OPENMP
                #pragma omp for schedule(dynamic, 1)
#endif
                for (slong index = 0; index < count; ++index) {
                    if (k == size) {
                        slong child = size - 1 - index;
                        if (entries[index].length && previous[child].length) {
#ifdef DRSOLVE_DET_TESTING
                            drsolve_det_test_event(1, size);
#endif
                            if (!DP_MUL(current + index, entries + index, previous + child, nc)) failed=1;
                        }
                        /* Each root child has exactly one reader. */
                        DP_CLEAR(previous + child, nc);
                        DP_INIT(previous + child, nc);
                        continue;
                    }
                    slong cols[FLINT_BITS], col = size - 1;
                    ulong prefix[FLINT_BITS], suffix[FLINT_BITS], rank = index;
                    for (slong j = k; j > 0; --j) {
                        while (choose[col][j] > rank) --col;
                        cols[j-1] = col; rank -= choose[col][j]; --col;
                    }
                    prefix[0] = 0; suffix[k] = 0;
                    for (slong j = 0; j < k; ++j)
                        prefix[j+1] = prefix[j] + choose[cols[j]][j+1];
                    for (slong j = k-1; j > 0; --j)
                        suffix[j] = suffix[j+1] + choose[cols[j]][j];
                    for (slong j = 0; j < k; ++j) {
                        ulong child = prefix[j] + suffix[j+1];
                        DP_POLY *entry = entries + (size-k)*size + cols[j];
                        if (!entry->length || !previous[child].length) continue;
                        if (!DP_MUL(product, entry, previous + child, nc)) { failed=1; continue; }
                        if (!current[index].length) DP_NAME(swap)(current + index, product);
                        else {
                            DP_NAME(add)(sum, current + index, product, nc);
                            DP_NAME(swap)(current + index, sum);
                        }
                    }
                }
                if (k == size) {
                    for (slong stride = 1; stride < count; stride *= 2) {
#ifdef _OPENMP
                        #pragma omp for schedule(static)
#endif
                        for (slong left = 0; left < count; left += 2*stride) {
                            if (left + stride >= count) continue;
                            DP_NAME(add)(sum, current + left, current + left+stride, nc);
                            DP_NAME(swap)(current + left, sum);
                            DP_CLEAR(current + left+stride, nc);
                            DP_INIT(current + left+stride, nc);
                        }
                    }
                }
                DP_CLEAR(product, nc); DP_CLEAR(sum, nc);
            }
        }
        for (slong i = 0; i < previous_count; ++i) DP_CLEAR(previous+i, nc);
        flint_free(previous); previous = current; previous_count = count;
        if (failed) break;
    }
    if (!failed) DP_EXPORT(GET_FQ_POLY(result), previous, ctx, nc);
    for (slong i = 0; i < previous_count; ++i) DP_CLEAR(previous+i, nc);
    for (slong i = 0; i < size*size; ++i) DP_CLEAR(entries+i, nc);
    flint_free(previous); flint_free(entries); DP_CTX_CLEAR(nc);
    return !failed;
}

#undef DP_NAME
#undef DP_CTX_T
#undef DP_POLY
#undef DP_POLY_T
#undef DP_E
#undef DP_SET
#undef DP_INIT
#undef DP_CLEAR
#undef DP_CTX_INIT
#undef DP_CTX_CLEAR
#undef DP_IMPORT
#undef DP_EXPORT
#undef DP_ADD
#undef DP_ZERO
#undef DP_MUL

/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_GF28_MINOR_DP_H
#define DRSOLVE_GF28_MINOR_DP_H

/* Byte coefficients for the whole DP, including additions. GF(16) and
   GF(256) share their polynomial layout and characteristic-two addition. */
static void byte_mpoly_swap(gf28_mpoly_struct *a, gf28_mpoly_struct *b)
{
    gf28_mpoly_struct t = *a; *a = *b; *b = t;
}

/* Nonaliasing merge into reusable output storage. */
static void byte_mpoly_add(gf28_mpoly_struct *out,
    const gf28_mpoly_struct *a, const gf28_mpoly_struct *b,
    const gf28_mpoly_ctx_t ctx)
{
    FLINT_ASSERT(out != a && out != b);
    if (!a->length) { gf28_mpoly_set(out, b, ctx); return; }
    if (!b->length) { gf28_mpoly_set(out, a, ctx); return; }
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
        out->coeffs = flint_realloc(out->coeffs, capacity);
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
        uint8_t c;
        if (cmp > 0) c = a->coeffs[i++];
        else if (cmp < 0) c = b->coeffs[j++];
        else { c = a->coeffs[i++] ^ b->coeffs[j++]; }
        if (c) {
            out->coeffs[len] = c;
            memcpy(out->exps + N*len++, exp, sizeof(ulong) * N);
        }
    }
    out->length = len; out->bits = bits;
    flint_free(mask); flint_free(ae); flint_free(be);
}

/* Called only after the common DP entry-budget and sparsity preflight.
   Convert the small input matrix once; convert only the final determinant
   back to fq_nmod. All live minors have one-byte coefficients. */
static void compute_byte_mpoly_det_dp(unified_mpoly_t result,
    unified_mpoly_t **matrix, slong size, unified_mpoly_ctx_t ctx,
    int use_parallel, const ulong choose[FLINT_BITS][FLINT_BITS])
{
    int is_gf24 = result->field_id == FIELD_ID_GF24;
    int (*mul)(gf28_mpoly_struct *, const gf28_mpoly_struct *,
               const gf28_mpoly_struct *, const gf28_mpoly_ctx_struct *) =
        is_gf24 ? gf24_mpoly_mul : gf28_mpoly_mul;
    gf28_mpoly_ctx_t nc;
    gf28_mpoly_ctx_init(nc, ctx->nvars, ctx->ord);
    gf28_mpoly_struct *entries = flint_calloc(size * size, sizeof(*entries));
    for (slong i = 0; i < size; ++i)
        for (slong j = 0; j < size; ++j) {
            gf28_mpoly_init(entries + i*size+j, nc);
            if (is_gf24)
                fq_nmod_mpoly_to_gf24_mpoly(entries + i*size+j,
                    GET_FQ_POLY(matrix[i][j]), ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
            else
                fq_nmod_mpoly_to_gf28_mpoly(entries + i*size+j,
                    GET_FQ_POLY(matrix[i][j]), ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
        }
    gf28_mpoly_struct *previous = NULL;
    slong previous_count = 0;
    for (slong k = 1; k <= size; ++k) {
        slong count = k == size ? size : (slong)choose[size][k];
        gf28_mpoly_struct *current = flint_malloc(count * sizeof(*current));
        for (slong i = 0; i < count; ++i) gf28_mpoly_init(current + i, nc);
        if (k == 1) {
            for (slong i = 0; i < count; ++i)
                gf28_mpoly_set(current + i, entries + (size-1)*size+i, nc);
        } else {
#ifdef _OPENMP
            #pragma omp parallel if(use_parallel && count > 1 && !omp_in_parallel()) num_threads(FLINT_MIN(count, omp_get_max_threads()))
#endif
            {
                gf28_mpoly_t product, sum;
                gf28_mpoly_init(product, nc); gf28_mpoly_init(sum, nc);
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
                            mul(current + index, entries + index, previous + child, nc);
                        }
                        /* Each root child has exactly one reader. */
                        gf28_mpoly_clear(previous + child, nc);
                        gf28_mpoly_init(previous + child, nc);
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
                        gf28_mpoly_struct *entry = entries + (size-k)*size + cols[j];
                        if (!entry->length || !previous[child].length) continue;
                        mul(product, entry, previous + child, nc);
                        if (!current[index].length) byte_mpoly_swap(current + index, product);
                        else {
                            byte_mpoly_add(sum, current + index, product, nc);
                            byte_mpoly_swap(current + index, sum);
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
                            byte_mpoly_add(sum, current + left, current + left+stride, nc);
                            byte_mpoly_swap(current + left, sum);
                            gf28_mpoly_clear(current + left+stride, nc);
                            gf28_mpoly_init(current + left+stride, nc);
                        }
                    }
                }
                gf28_mpoly_clear(product, nc); gf28_mpoly_clear(sum, nc);
            }
        }
        for (slong i = 0; i < previous_count; ++i) gf28_mpoly_clear(previous+i, nc);
        flint_free(previous); previous = current; previous_count = count;
    }
    if (is_gf24)
        gf24_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(result), previous,
                                  ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
    else
        gf28_mpoly_to_fq_nmod_mpoly(GET_FQ_POLY(result), previous,
                                  ctx->field_ctx->ctx.fq_ctx, GET_FQ_CTX(ctx));
    for (slong i = 0; i < previous_count; ++i) gf28_mpoly_clear(previous+i, nc);
    for (slong i = 0; i < size*size; ++i) gf28_mpoly_clear(entries+i, nc);
    flint_free(previous); flint_free(entries); gf28_mpoly_ctx_clear(nc);
}
#endif

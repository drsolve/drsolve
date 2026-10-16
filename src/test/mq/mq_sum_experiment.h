/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Experimental Laplace sum kernels. Included only by explicit test builds;
 * paired measurements found both slower than the production accumulation. */
/* A single heap for the entire alternating Laplace sum. Each stream is a
 * monomial multiple of one child minor; cancellation and projection happen
 * before the sole output write. The previous DP layer is immutable. */
typedef struct
{
    const nmod_mpoly_struct *poly;
    const ulong *shift;
    ulong coeff;
    slong term;
    ulong exp[3];
} mq_sum_stream;

static int mq_direct_sum(nmod_mpoly_t out, nmod_mpoly_struct *const *a,
                         nmod_mpoly_struct *const *b, slong count,
                         const nmod_mpoly_ctx_t ctx, const mq_det_filter *f, int truncate)
{
    if (!f || !f->arithmetic_bits || ctx->minfo->ord != ORD_LEX || count < 2)
        return 0;
    slong words = mpoly_words_per_exp(f->arithmetic_bits, ctx->minfo);
    if (words > 3) return 0;
    slong capacity = 0, reserve = 0;
    for (slong j = 0; j < count; j++) {
        if (!a[j]->length || !b[j]->length) continue;
        if (a[j]->bits != f->arithmetic_bits || b[j]->bits != f->arithmetic_bits ||
            a[j]->length > 2*f->nvars+2 || out == a[j] || out == b[j]) return 0;
        capacity += a[j]->length;
        reserve = FLINT_MAX(reserve, b[j]->length);
    }
    if (reserve < 256 || capacity < 2) return 0;
    mq_sum_stream *streams = flint_malloc((size_t) capacity * sizeof(*streams));
    slong *heap = flint_malloc((size_t) capacity * sizeof(*heap));
    slong heaplen = 0;
    for (slong j = 0; j < count; j++) {
        if (!b[j]->length) continue;
        for (slong i = 0; i < a[j]->length; i++) {
            mq_sum_stream next;
            next.poly = b[j]; next.term = 0;
            next.shift = a[j]->exps + i*words;
            next.coeff = j & 1 ? nmod_neg(a[j]->coeffs[i], ctx->mod) : a[j]->coeffs[i];
            for (slong w = 0; w < words; w++) next.exp[w] = next.shift[w] + b[j]->exps[w];
            slong id = heaplen, pos = heaplen++;
            streams[id] = next;
            while (pos && mq_linear_cmp(next.exp, streams[heap[(pos-1)/2]].exp, words) > 0) {
                heap[pos] = heap[(pos-1)/2]; pos = (pos-1)/2;
            }
            heap[pos] = id;
        }
    }
    slong offset[2], shift[2], low[2], width = f->nvars * f->arithmetic_bits;
    slong vo[2 * FLINT_BITS], vs[2 * FLINT_BITS];
    if (!f->packed_bits)
        for (slong v = 0; v < 2 * f->nvars; v++)
            mpoly_gen_offset_shift_sp(vo + v, vs + v, v, f->arithmetic_bits, ctx->minfo);
    for (slong axis = 0; f->packed_bits && axis < 2; axis++)
    {
        mpoly_gen_offset_shift_sp(offset + axis, shift + axis, (axis + 1) * f->nvars - 1,
                                  f->packed_bits, ctx->minfo);
        low[axis] =
            FLINT_MIN(width, ((FLINT_BITS - shift[axis]) / f->packed_bits) * f->packed_bits);
    }
    nmod_mpoly_fit_length_reset_bits(out, reserve, f->arithmetic_bits, ctx);
    slong len = 0;
    while (heaplen)
    {
        ulong key[3];
        for (slong w = 0; w < words; w++)
            key[w] = streams[heap[0]].exp[w];
        ulong value = 0;
        do
        {
            slong id = heap[0];
            mq_sum_stream *next = streams + id;
            value = nmod_add(value, nmod_mul(next->coeff, next->poly->coeffs[next->term], ctx->mod), ctx->mod);
            if (++next->term == next->poly->length)
                id = heap[--heaplen];
            else
                for (slong w = 0; w < words; w++)
                    next->exp[w] = next->shift[w] + next->poly->exps[next->term * words + w];
            if (heaplen)
            {
                slong pos = 0;
                for (;;)
                {
                    slong child = 2 * pos + 1;
                    if (child >= heaplen)
                        break;
                    if (child + 1 < heaplen &&
                        mq_linear_cmp(streams[heap[child + 1]].exp, streams[heap[child]].exp, words) > 0)
                        child++;
                    if (mq_linear_cmp(streams[id].exp, streams[heap[child]].exp, words) >= 0)
                        break;
                    heap[pos] = heap[child];
                    pos = child;
                }
                heap[pos] = id;
            }
        } while (heaplen && mq_linear_cmp(streams[heap[0]].exp, key, words) == 0);
        if (!value)
            continue;
        if (truncate)
        {
            ulong r = 0, c = 0;
            if (f->packed_bits)
            {
                r = mq_packed_axis(key, offset[0], shift[0], low[0], width);
                c = mq_packed_axis(key, offset[1], shift[1], low[1], width);
                if (!mq_monom_contains(f->packed, r) || !mq_monom_contains(f->packed + 1, c))
                    continue;
            }
            else
            {
                ulong mask = UWORD_MAX >> (FLINT_BITS - f->arithmetic_bits);
                int keep = 1;
                for (slong v = 0; v < f->nvars; v++)
                {
                    ulong a = (key[vo[v]] >> vs[v]) & mask;
                    ulong b = (key[vo[f->nvars + v]] >> vs[f->nvars + v]) & mask;
                    if (a > f->digit_mask || b > f->digit_mask)
                    {
                        keep = 0;
                        break;
                    }
                    r |= a << (v * f->bits);
                    c |= b << (v * f->bits);
                }
                if (!keep || !mq_monom_contains(&f->rows, r) || !mq_monom_contains(&f->cols, c))
                    continue;
            }
        }
        nmod_mpoly_fit_length(out, len + 1, ctx);
        out->coeffs[len] = value;
        for (slong w = 0; w < words; w++)
            out->exps[len * words + w] = key[w];
        len++;
    }
    flint_free(heap);
    flint_free(streams);
    _nmod_mpoly_set_length(out, len, ctx);
    return 1;
}

static int mq_product_sum(nmod_mpoly_t out, nmod_mpoly_struct *const *a,
                         nmod_mpoly_struct *const *b, slong count,
                         const nmod_mpoly_ctx_t ctx, const mq_det_filter *f, int truncate)
{
    if (!f || !f->arithmetic_bits || ctx->minfo->ord != ORD_LEX || count < 2)
        return 0;
    slong words = mpoly_words_per_exp(f->arithmetic_bits, ctx->minfo);
    if (words > 3) return 0;
    slong capacity = 0, reserve = 0;
    for (slong j = 0; j < count; j++) {
        if (!a[j]->length || !b[j]->length) continue;
        if (a[j]->bits != f->arithmetic_bits || b[j]->bits != f->arithmetic_bits ||
            a[j]->length > 2*f->nvars+2 || out == a[j] || out == b[j]) return 0;
        capacity += a[j]->length;
        reserve = FLINT_MAX(reserve, b[j]->length);
    }
    if (reserve < 256 || capacity < 2) return 0;
    /* Materialize each product once, then merge all products together. */
    nmod_mpoly_struct *products = flint_malloc((size_t)count*sizeof(*products));
    mq_sum_stream *streams = flint_malloc((size_t)count*sizeof(*streams));
    slong *heap = flint_malloc((size_t)count*sizeof(*heap));
    slong heaplen = 0;
    ulong zero[3] = {0};
    for (slong j = 0; j < count; j++) {
        nmod_mpoly_init(products+j,ctx);
        if (!mq_linear_mul(products+j,a[j],b[j],ctx,f,truncate)) {
            nmod_mpoly_mul(products+j,a[j],b[j],ctx);
            if (truncate) mq_filter_poly(products+j,ctx,f,0);
        }
        if (!products[j].length) continue;
        /* Zero results may retain their initial packing; nonzero products
         * from generic multiplication must be aligned to the common width. */
        if (!nmod_mpoly_repack_bits_inplace(products+j,f->arithmetic_bits,ctx)) {
            for(slong i=0;i<=j;i++) nmod_mpoly_clear(products+i,ctx);
            flint_free(products);flint_free(streams);flint_free(heap);return 0;
        }
        mq_sum_stream next;
        next.poly=products+j;next.term=0;next.shift=zero;
        next.coeff=j&1 ? ctx->mod.n-1 : 1;
        for(slong w=0;w<words;w++)next.exp[w]=products[j].exps[w];
        slong id=heaplen,pos=heaplen++;streams[id]=next;
        while(pos && mq_linear_cmp(next.exp,streams[heap[(pos-1)/2]].exp,words)>0) {
            heap[pos]=heap[(pos-1)/2];pos=(pos-1)/2;
        }
        heap[pos]=id;
    }
    nmod_mpoly_fit_length_reset_bits(out, reserve, f->arithmetic_bits, ctx);
    slong len = 0;
    while (heaplen)
    {
        ulong key[3];
        for (slong w = 0; w < words; w++)
            key[w] = streams[heap[0]].exp[w];
        ulong value = 0;
        do
        {
            slong id = heap[0];
            mq_sum_stream *next = streams + id;
            value = next->coeff == 1 ? nmod_add(value,next->poly->coeffs[next->term],ctx->mod)
                                    : nmod_sub(value,next->poly->coeffs[next->term],ctx->mod);
            if (++next->term == next->poly->length)
                id = heap[--heaplen];
            else
                for (slong w = 0; w < words; w++)
                    next->exp[w] = next->shift[w] + next->poly->exps[next->term * words + w];
            if (heaplen)
            {
                slong pos = 0;
                for (;;)
                {
                    slong child = 2 * pos + 1;
                    if (child >= heaplen)
                        break;
                    if (child + 1 < heaplen &&
                        mq_linear_cmp(streams[heap[child + 1]].exp, streams[heap[child]].exp, words) > 0)
                        child++;
                    if (mq_linear_cmp(streams[id].exp, streams[heap[child]].exp, words) >= 0)
                        break;
                    heap[pos] = heap[child];
                    pos = child;
                }
                heap[pos] = id;
            }
        } while (heaplen && mq_linear_cmp(streams[heap[0]].exp, key, words) == 0);
        if (!value)
            continue;
        nmod_mpoly_fit_length(out, len + 1, ctx);
        out->coeffs[len] = value;
        for (slong w = 0; w < words; w++)
            out->exps[len * words + w] = key[w];
        len++;
    }
    for(slong j=0;j<count;j++) nmod_mpoly_clear(products+j,ctx);
    flint_free(products);
    flint_free(heap);
    flint_free(streams);
    _nmod_mpoly_set_length(out, len, ctx);
    return 1;
}


static int mq_sum_experiment_mode = DRSOLVE_MQ_SUM_TEST;
static int mq_linear_sum(nmod_mpoly_t out, nmod_mpoly_struct *const *a,
                        nmod_mpoly_struct *const *b, slong count,
                        const nmod_mpoly_ctx_t ctx, const mq_det_filter *f, int truncate)
{
    return mq_sum_experiment_mode == 1 ? mq_direct_sum(out,a,b,count,ctx,f,truncate)
                                       : mq_product_sum(out,a,b,count,ctx,f,truncate);
}

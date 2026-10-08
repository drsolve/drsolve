/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Internal, coefficient-independent planning and heap operations. */
#ifndef GF2N_MPOLY_SPARSE_H
#define GF2N_MPOLY_SPARSE_H

/* FLINT 3.6 has no mpoly_monomial_add_any_bits helper. Product planning
   reserves a guard bit per exponent, so these additions cannot overflow
   into the next packed field. Multi-limb exponents still require carry. */
static inline void gf2n_monomial_add(ulong *out, const ulong *a,
    const ulong *b, slong N, flint_bitcnt_t bits)
{
    if (bits <= FLINT_BITS)
        mpoly_monomial_add(out, a, b, N);
    else
        mpoly_monomial_add_mp(out, a, b, N);
}

static flint_bitcnt_t gf2n_product_bits(const ulong *a, slong alen,
    flint_bitcnt_t abits, const ulong *b, slong blen, flint_bitcnt_t bbits,
    const mpoly_ctx_t ctx, ulong *array_size, ulong *main_size)
{
    slong nf = ctx->nfields;
    fmpz *am = flint_calloc(2 * nf, sizeof(fmpz));
    fmpz *bm = am + nf;
    flint_bitcnt_t bits = MPOLY_MIN_BITS;
    *array_size = 1;
    *main_size = 1;
    mpoly_max_fields_fmpz(am, a, alen, abits, ctx);
    mpoly_max_fields_fmpz(bm, b, blen, bbits, ctx);
    for (slong i = 0; i < nf; ++i) {
        fmpz_add(am + i, am + i, bm + i);
        bits = FLINT_MAX(bits, fmpz_bits(am + i) + 1);
        /* The array path uses signed, machine-sized bounds. */
        if (!fmpz_fits_si(am + i) || fmpz_cmp_ui(am + i, WORD_MAX) >= 0) {
            *array_size = UWORD_MAX;
        } else {
            ulong width = fmpz_get_ui(am + i) + 1;
            if (i == nf - 1) *main_size = width;
            else if (*array_size > UWORD_MAX / width) *array_size = UWORD_MAX;
            else *array_size *= width;
        }
        fmpz_clear(am + i);
        fmpz_clear(bm + i);
    }
    flint_free(am);
    return mpoly_fix_bits(bits, ctx);
}

static int gf2n_prefer_array(const ulong *a, slong alen, flint_bitcnt_t abits,
    const ulong *b, slong blen, flint_bitcnt_t bbits,
    const mpoly_ctx_t ctx, ulong default_limit)
{
    if (!alen || !blen) return 1;
    if (ctx->ord != ORD_LEX || ctx->nvars < 1 || !abits || !bbits ||
        abits > FLINT_BITS || bbits > FLINT_BITS ||
        mpoly_words_per_exp(abits, ctx) > 2 ||
        mpoly_words_per_exp(bbits, ctx) > 2) return 0;
    ulong size, main_size;
    flint_bitcnt_t bits = gf2n_product_bits(a, alen, abits, b, blen, bbits,
                                           ctx, &size, &main_size);
    ulong limit = gf2n_mpoly_array_limit_k >= 0
        ? UWORD(1) << gf2n_mpoly_array_limit_k : default_limit;
    if (size > limit || bits > FLINT_BITS ||
        mpoly_words_per_exp(bits, ctx) > 2) return 0;
    /* Clearing and scanning every main-variable slice should not dwarf
       the number of coefficient products. This also catches sparse gaps. */
    return (double)size * main_size <=
        FLINT_MAX(4096.0, 2.0 * (double)alen * (double)blen);
}

typedef struct {
    slong row, col;
    ulong *exp;
} gf2n_product_heap_entry;

static void gf2n_heap_sift_down(gf2n_product_heap_entry *heap, slong len,
                              slong pos, slong N, const ulong *mask)
{
    gf2n_product_heap_entry item = heap[pos];
    while (pos < len / 2) {
        slong child = 2 * pos + 1;
        if (child + 1 < len && mpoly_monomial_cmp(
                heap[child + 1].exp, heap[child].exp, N, mask) > 0) ++child;
        if (mpoly_monomial_cmp(item.exp, heap[child].exp, N, mask) >= 0) break;
        heap[pos] = heap[child];
        pos = child;
    }
    heap[pos] = item;
}

/* Each row is a sorted stream a_i * b. Keep only its next product in the
   heap, merge equal exponents, and emit canonical terms. Auxiliary storage
   is O(N*(alen+blen)); there is no allocation proportional to the degree box.
   All coefficient operations are specialized at compile time. */
#define DEFINE_GF2N_SPARSE_MUL(FIELD, COEFF_T, IS_ZERO, ADD, MUL, INIT) \
int FIELD##_mpoly_mul_sparse(FIELD##_mpoly_t res, \
    const FIELD##_mpoly_t a, const FIELD##_mpoly_t b, \
    const FIELD##_mpoly_ctx_t ctx) \
{ \
    if (!a->length || !b->length) { FIELD##_mpoly_zero(res, ctx); return 1; } \
    if (!a->bits || !b->bits || ctx->minfo->nvars < 1) return 0; \
    INIT; \
    const FIELD##_mpoly_struct *x = a, *y = b; \
    if (x->length > y->length) { x = b; y = a; } \
    ulong unused_size, unused_main; \
    flint_bitcnt_t bits = gf2n_product_bits(x->exps, x->length, x->bits, \
        y->exps, y->length, y->bits, ctx->minfo, &unused_size, &unused_main); \
    slong N = mpoly_words_per_exp(bits, ctx->minfo); \
    ulong *xcopy = NULL, *ycopy = NULL; \
    const ulong *xe = x->exps, *ye = y->exps; \
    int ok = 1; \
    if (x->bits != bits) { \
        xcopy = flint_malloc(sizeof(ulong) * N * x->length); \
        ok = mpoly_repack_monomials(xcopy, bits, xe, x->bits, x->length, ctx->minfo); \
        xe = xcopy; \
    } \
    if (ok && y->bits != bits) { \
        ycopy = flint_malloc(sizeof(ulong) * N * y->length); \
        ok = mpoly_repack_monomials(ycopy, bits, ye, y->bits, y->length, ctx->minfo); \
        ye = ycopy; \
    } \
    if (!ok) { flint_free(xcopy); flint_free(ycopy); return 0; } \
    ulong *keys = flint_malloc(sizeof(ulong) * N * x->length); \
    ulong *mask = flint_malloc(sizeof(ulong) * N); \
    ulong *exp = flint_malloc(sizeof(ulong) * N); \
    mpoly_get_cmpmask(mask, N, bits, ctx->minfo); \
    gf2n_product_heap_entry *heap = flint_malloc(sizeof(*heap) * x->length); \
    slong len = x->length; \
    for (slong i = 0; i < len; ++i) { \
        heap[i] = (gf2n_product_heap_entry){i, 0, keys + N*i}; \
        gf2n_monomial_add(heap[i].exp, xe + N*i, ye, N, bits); \
    } \
    for (slong i = len / 2; i-- > 0;) gf2n_heap_sift_down(heap, len, i, N, mask); \
    FIELD##_mpoly_t out; \
    FIELD##_mpoly_init(out, ctx); \
    out->bits = bits; \
    while (len) { \
        memcpy(exp, heap[0].exp, sizeof(ulong) * N); \
        COEFF_T coeff = {0}; \
        do { \
            slong i = heap[0].row, j = heap[0].col; \
            COEFF_T product = MUL(x->coeffs[i], y->coeffs[j]); \
            coeff = ADD(coeff, product); \
            if (++j < y->length) { \
                heap[0].col = j; \
                gf2n_monomial_add(heap[0].exp, xe + N*i, ye + N*j, N, bits); \
            } else { heap[0] = heap[--len]; } \
            if (len) gf2n_heap_sift_down(heap, len, 0, N, mask); \
        } while (len && mpoly_monomial_equal(exp, heap[0].exp, N)); \
        if (!IS_ZERO(coeff)) { \
            slong k = out->length; \
            if (k == out->coeffs_alloc) { \
                slong alloc = FLINT_MAX(16, 2 * out->coeffs_alloc); \
                out->coeffs = flint_realloc(out->coeffs, sizeof(COEFF_T) * alloc); \
                out->exps = flint_realloc(out->exps, sizeof(ulong) * N * alloc); \
                out->coeffs_alloc = alloc; out->exps_alloc = N * alloc; \
            } \
            out->coeffs[k] = coeff; \
            memcpy(out->exps + N*k, exp, sizeof(ulong) * N); \
            out->length++; \
        } \
    } \
    FIELD##_mpoly_struct swap = *res; *res = *out; *out = swap; \
    FIELD##_mpoly_clear(out, ctx); \
    flint_free(heap); flint_free(exp); flint_free(mask); flint_free(keys); \
    flint_free(ycopy); flint_free(xcopy); \
    return 1; \
}

#endif

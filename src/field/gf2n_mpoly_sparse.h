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
/* Single-word exponents use FLINT's chained heap: equal monomials share a
 * heap node, so collisions do not each require a full sift-down. */
typedef struct {
    mpoly_heap1_s *heap;
    mpoly_heap_t *rows, *pending;
    const ulong *a, *b;
    slong length, next_loc, columns;
    ulong mask, exponent;
} gf2n_product_stream;

static void gf2n_product_stream_init(gf2n_product_stream *s,
    const ulong *a, slong rows, const ulong *b, slong columns, ulong mask)
{
    s->heap=flint_malloc((rows+1)*sizeof(*s->heap));
    s->rows=flint_malloc(rows*sizeof(*s->rows));
    s->pending=NULL;s->a=a;s->b=b;s->columns=columns;
    s->length=1;s->next_loc=rows+1;s->mask=mask;
    for(slong i=0;i<rows;i++) {
        s->rows[i]=(mpoly_heap_t){i,0,NULL};
        _mpoly_heap_insert1(s->heap,a[i]+b[0],s->rows+i,
                           &s->next_loc,&s->length,mask);
    }
}
static inline int gf2n_product_stream_active(const gf2n_product_stream *s)
{ return s->pending!=NULL || s->length>1; }
static inline void gf2n_product_stream_pop(gf2n_product_stream *s,slong *i,slong *j)
{
    if(!s->pending) {
        s->exponent=s->heap[1].exp;
        s->pending=_mpoly_heap_pop1(s->heap,&s->length,s->mask);
    }
    mpoly_heap_t *row=s->pending;
    s->pending=row->next;*i=row->i;*j=row->j;
    if(++row->j<(ulong)s->columns) {
        row->next=NULL;
        _mpoly_heap_insert1(s->heap,s->a[row->i]+s->b[row->j],row,
                           &s->next_loc,&s->length,s->mask);
    }
}
static void gf2n_product_stream_clear(gf2n_product_stream *s)
{ flint_free(s->heap);flint_free(s->rows); }

/* Common product iterator for fixed and arbitrary moduli. Multi-word
 * exponents retain the general heap and FLINT's comparison mask. */
typedef struct {
    gf2n_product_stream stream;
    gf2n_product_heap_entry *heap;
    ulong *keys;
    const ulong *a, *b, *mask;
    slong length, columns, N;
    flint_bitcnt_t bits;
} gf2n_product_iterator;
static void gf2n_product_iterator_init(gf2n_product_iterator *it,
    const ulong *a,slong rows,const ulong *b,slong columns,
    slong N,flint_bitcnt_t bits,const ulong *mask)
{
    it->N=N;it->bits=bits;it->mask=mask;
    it->a=a;it->b=b;it->columns=columns;it->length=rows;
    if(N==1) {
        gf2n_product_stream_init(&it->stream,a,rows,b,columns,mask[0]);
        return;
    }
    it->keys=flint_malloc(N*rows*sizeof(ulong));
    it->heap=flint_malloc(rows*sizeof(*it->heap));
    for(slong i=0;i<rows;i++) {
        it->heap[i]=(gf2n_product_heap_entry){i,0,it->keys+N*i};
        gf2n_monomial_add(it->heap[i].exp,a+N*i,b,N,bits);
    }
    for(slong i=rows/2;i-- >0;) gf2n_heap_sift_down(it->heap,rows,i,N,mask);
}
static inline int gf2n_product_iterator_active(const gf2n_product_iterator *it)
{ return it->N==1 ? gf2n_product_stream_active(&it->stream) : it->length!=0; }
static inline const ulong *gf2n_product_iterator_key(const gf2n_product_iterator *it)
{
    if(it->N!=1) return it->heap[0].exp;
    return it->stream.pending ? &it->stream.exponent : &it->stream.heap[1].exp;
}
static inline void gf2n_product_iterator_pop(gf2n_product_iterator *it,slong *i,slong *j)
{
    if(it->N==1) { gf2n_product_stream_pop(&it->stream,i,j);return; }
    *i=it->heap[0].row;*j=it->heap[0].col;
    if(*j+1<it->columns) {
        it->heap[0].col=*j+1;
        gf2n_monomial_add(it->heap[0].exp,it->a+it->N*(*i),
                         it->b+it->N*(*j+1),it->N,it->bits);
    } else it->heap[0]=it->heap[--it->length];
    if(it->length) gf2n_heap_sift_down(it->heap,it->length,0,it->N,it->mask);
}
static void gf2n_product_iterator_clear(gf2n_product_iterator *it)
{
    if(it->N==1) gf2n_product_stream_clear(&it->stream);
    else { flint_free(it->keys);flint_free(it->heap); }
}

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
    ulong *mask = flint_malloc(sizeof(ulong) * N); \
    ulong *exp = flint_malloc(sizeof(ulong) * N); \
    mpoly_get_cmpmask(mask, N, bits, ctx->minfo); \
    gf2n_product_iterator iterator; \
    gf2n_product_iterator_init(&iterator,xe,x->length,ye,y->length,N,bits,mask); \
    FIELD##_mpoly_t out; \
    FIELD##_mpoly_init(out, ctx); \
    out->bits = bits; \
    while (gf2n_product_iterator_active(&iterator)) { \
        memcpy(exp,gf2n_product_iterator_key(&iterator),sizeof(ulong)*N); \
        COEFF_T coeff = {0}; \
        do { \
            slong i,j;gf2n_product_iterator_pop(&iterator,&i,&j); \
            COEFF_T product = MUL(x->coeffs[i], y->coeffs[j]); \
            coeff = ADD(coeff, product); \
        } while (gf2n_product_iterator_active(&iterator) && \
                 mpoly_monomial_equal(exp,gf2n_product_iterator_key(&iterator),N)); \
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
    gf2n_product_iterator_clear(&iterator); \
    flint_free(exp); flint_free(mask); \
    flint_free(ycopy); flint_free(xcopy); \
    return 1; \
}

#endif

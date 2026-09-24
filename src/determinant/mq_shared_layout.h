/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Private MQ coefficient-array backend, included after the filter helpers. */
#include <stdint.h>

#ifdef DRSOLVE_MQ_LAYOUT_TEST
static int mq_shared_test_used;
#endif

typedef struct {
    ulong *keys;
    uint32_t *table;
    size_t count, capacity, buckets;
    slong words;
    flint_bitcnt_t bits;
} mq_shared_support;

static size_t mq_shared_hash(const ulong *key, slong words)
{
    ulong h = 0x9e3779b9UL;
    for (slong w = 0; w < words; w++) h = mq_monom_hash(h ^ key[w]);
    return h;
}
static void mq_shared_rehash(mq_shared_support *s, size_t buckets)
{
    flint_free(s->table);
    s->buckets = buckets;
    s->table = flint_calloc(buckets, sizeof(uint32_t));
    for (size_t i = 0; i < s->count; i++) {
        size_t pos = mq_shared_hash(s->keys+i*s->words, s->words) & (buckets-1);
        while (s->table[pos]) pos = (pos+1) & (buckets-1);
        s->table[pos] = (uint32_t)i+1;
    }
}
static uint32_t mq_shared_find(const mq_shared_support *s, const ulong *key)
{
    size_t pos = mq_shared_hash(key, s->words) & (s->buckets-1);
    while (s->table[pos]) {
        uint32_t id = s->table[pos]-1;
        if (!memcmp(key, s->keys+(size_t)id*s->words, s->words*sizeof(ulong))) return id;
        pos = (pos+1) & (s->buckets-1);
    }
    return UINT32_MAX;
}
static uint32_t mq_shared_insert(mq_shared_support *s, const ulong *key)
{
    if (2*(s->count+1) >= s->buckets) mq_shared_rehash(s, 2*s->buckets);
    size_t pos = mq_shared_hash(key, s->words) & (s->buckets-1);
    while (s->table[pos]) {
        uint32_t id = s->table[pos]-1;
        if (!memcmp(key, s->keys+(size_t)id*s->words, s->words*sizeof(ulong))) return id;
        pos = (pos+1) & (s->buckets-1);
    }
    FLINT_ASSERT(s->count < UINT32_MAX);
    if (s->count == s->capacity) {
        s->capacity = s->capacity ? 2*s->capacity : 16;
        s->keys = flint_realloc(s->keys, s->capacity*s->words*sizeof(ulong));
    }
    memcpy(s->keys+s->count*s->words, key, s->words*sizeof(ulong));
    s->table[pos] = (uint32_t)s->count+1;
    return (uint32_t)s->count++;
}
static void mq_shared_empty(mq_shared_support *s, flint_bitcnt_t bits, slong words)
{
    memset(s, 0, sizeof(*s)); s->bits = bits; s->words = words;
    mq_shared_rehash(s, 32);
}
static void mq_shared_support_init(mq_shared_support *s, const nmod_mpoly_ctx_t ctx, slong degree)
{
    flint_bitcnt_t bits = mpoly_fix_bits(1+FLINT_BIT_COUNT((ulong)degree), ctx->minfo);
    slong words = mpoly_words_per_exp(bits, ctx->minfo);
    FLINT_ASSERT(words <= 3);
    mq_shared_empty(s, bits, words);
    ulong zero[3] = {0}; mq_shared_insert(s, zero);
}
static void mq_shared_support_clear(mq_shared_support *s)
{
    flint_free(s->keys); flint_free(s->table); memset(s, 0, sizeof(*s));
}
static size_t mq_shared_support_bytes(const mq_shared_support *s)
{
    return s->capacity*s->words*sizeof(ulong)+s->buckets*sizeof(uint32_t);
}
static void mq_shared_row(mq_shared_support *shifts, nmod_mpoly_t *row,
                          slong n, const mq_shared_support *previous,
                          const nmod_mpoly_ctx_t ctx)
{
    mq_shared_empty(shifts, previous->bits, previous->words);
    nmod_mpoly_t p; nmod_mpoly_init(p, ctx);
    for (slong col = 0; col < n; col++) {
        int packed = nmod_mpoly_repack_bits(p, row[col], previous->bits, ctx);
        FLINT_ASSERT(packed); (void)packed;
        for (slong t = 0; t < p->length; t++)
            mq_shared_insert(shifts, p->exps+t*previous->words);
    }
    nmod_mpoly_clear(p, ctx);
}
static int mq_shared_keep(const ulong *key, const mq_det_filter *f,
                          const nmod_mpoly_ctx_t ctx, flint_bitcnt_t bits)
{
    if (!f) return 1;
    if (f->packed_bits == bits) {
        slong width = f->nvars*bits;
        for (slong axis = 0; axis < 2; axis++) {
            slong word, shift;
            mpoly_gen_offset_shift_sp(&word, &shift, (axis+1)*f->nvars-1, bits, ctx->minfo);
            slong low = FLINT_MIN(width, ((FLINT_BITS-shift)/bits)*bits);
            ulong code = mq_packed_axis(key, word, shift, low, width);
            if (!mq_monom_contains(&f->packed[axis], code)) return 0;
        }
        return 1;
    }
    ulong exp[64];
    mpoly_get_monomial_ui(exp, key, bits, ctx->minfo);
    return mq_filter_accepts(f, exp, 0);
}
static uint32_t *mq_shared_next(mq_shared_support *next, const mq_shared_support *prev,
                                const mq_shared_support *shifts, const mq_det_filter *f,
                                const nmod_mpoly_ctx_t ctx, int maps)
{
    mq_shared_empty(next, prev->bits, prev->words);
    uint32_t *map = maps ? flint_malloc(shifts->count*prev->count*sizeof(uint32_t)) : NULL;
    for (size_t a = 0; a < shifts->count; a++) for (size_t b = 0; b < prev->count; b++) {
        ulong key[3];
        for (slong w = 0; w < prev->words; w++)
            key[w] = shifts->keys[a*prev->words+w]+prev->keys[b*prev->words+w];
        uint32_t id = mq_shared_keep(key, f, ctx, prev->bits) ? mq_shared_insert(next, key) : UINT32_MAX;
        if (maps) map[a*prev->count+b] = id;
    }
    return map;
}

/* Saturating total-degree simplex size. Only small admitted sizes reach the
 * allocator; rejection leaves the caller's output untouched. */
static size_t mq_shared_binomial(slong n, slong k, size_t cap)
{
    size_t value = 1;
    k = FLINT_MIN(k, n-k);
    for (slong j = 1; j <= k; j++) {
        if (value > SIZE_MAX/(n-k+j)) return cap+1;
        value = value*(n-k+j)/j;
        if (value > cap) return cap+1;
    }
    return value;
}

static size_t mq_shared_support_bound(size_t count, slong words)
{
    size_t capacity = 16, buckets = 32;
    while (capacity < count) capacity *= 2;
    while (buckets <= 2*(count+1)) buckets *= 2;
    /* Include simultaneous old/new key storage during realloc growth. */
    return 2*capacity*words*sizeof(ulong)+buckets*sizeof(uint32_t);
}

#ifndef DRSOLVE_MQ_SHARED_WORKSPACE_BYTES
#define DRSOLVE_MQ_SHARED_WORKSPACE_BYTES (256UL * 1024 * 1024)
#endif

static int mq_shared_admit(nmod_mpoly_t **matrix, slong n,
                           const nmod_mpoly_ctx_t ctx,
                           ulong choose[FLINT_BITS][FLINT_BITS])
{
    const size_t budget = DRSOLVE_MQ_SHARED_WORKSPACE_BYTES;
    const size_t cap = FLINT_MIN(budget/sizeof(uint32_t), UINT32_MAX/4);
    slong nv = ctx->minfo->nvars;
    if (FLINT_BITS != 64 || n < 4 || nv != 2*n-1 || nv >= FLINT_BITS || ctx->minfo->ord != ORD_LEX)
        return 0;
    flint_bitcnt_t bits = mpoly_fix_bits(1+FLINT_BIT_COUNT((ulong)n+1), ctx->minfo);
    slong words = mpoly_words_per_exp(bits, ctx->minfo);
    if (words > 3) return 0;
    ulong active = 0, exp[FLINT_BITS];
    slong degree = 0;
    size_t previous = 1, previous_count = 1;
    mq_shared_support packing;
    mq_shared_support_init(&packing, ctx, n+1);
    int ok = 1;
    for (slong k = 1; ok && k <= n; k++) {
        slong row = n-k, row_degree = 0;
        for (slong c = 0; ok && c < n; c++) {
            const nmod_mpoly_struct *p = matrix[row][c];
            if (!nmod_mpoly_degrees_fit_si(p, ctx)) { ok = 0; break; }
            for (slong t = 0; ok && t < p->length; t++) {
                slong d = 0;
                nmod_mpoly_get_term_exp_ui(exp, p, t, ctx);
                for (slong v = 0; v < nv; v++) {
                    if (exp[v] > (ulong)(row == 0 ? 2 : 1)-d) { ok = 0; break; }
                    d += exp[v];
                    if (exp[v]) active |= UWORD(1) << v;
                }
                row_degree = FLINT_MAX(row_degree, d);
            }
        }
        if (!ok) break;
        degree += row_degree;
        slong variables = 0;
        for (ulong mask = active; mask; mask &= mask-1) variables++;
        size_t current = mq_shared_binomial(variables+degree, degree, cap);
        if (current > cap) { ok = 0; break; }
        mq_shared_support shifts;
        mq_shared_row(&shifts, matrix[row], n, &packing, ctx);
        size_t count = k == n ? n : choose[n][k];
        /* All products below fit size_t: current <= cap <= UINT32_MAX/4,
         * count is bounded by the existing DP entry limit, and the n/packing
         * checks restrict n <= 12 on the supported 64-bit layout. Check the
         * coefficient product explicitly before adding the small buffers. */
        if (count > budget/sizeof(ulong)/current ||
            previous_count > budget/sizeof(ulong)/previous) ok = 0;
        if (ok) {
            size_t bytes = (count*current+previous_count*previous+n*shifts.count)*sizeof(ulong)
                +shifts.count*previous*sizeof(uint32_t)
                +mq_shared_support_bound(current, words)
                +mq_shared_support_bound(previous, words)
                +mq_shared_support_bytes(&shifts);
            /* Output packing allowance; FLINT sort scratch, input storage and
             * allocator overhead are not a process-RSS guarantee. */
            if (k == n) bytes += current*(words+1)*sizeof(ulong);
            if (bytes > budget) ok = 0;
        }
        mq_shared_support_clear(&shifts);
        previous = current; previous_count = count;
    }
    mq_shared_support_clear(&packing);
    return ok;
}

static void mq_shared_pack(nmod_mpoly_t out, const ulong *coeffs,
                           const mq_shared_support *s, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_fit_length_reset_bits(out, s->count, s->bits, ctx);
    slong len = 0;
    for (size_t i = 0; i < s->count; i++) if (coeffs[i]) {
        out->coeffs[len] = coeffs[i];
        memcpy(out->exps+len*s->words, s->keys+i*s->words, s->words*sizeof(ulong)); len++;
    }
    _nmod_mpoly_set_length(out, len, ctx);
    nmod_mpoly_sort_terms(out, ctx);
}

/* Shared support and coefficients for all layers, including the root.
 * Preflight proves the degree/packing invariants and bounds owned workspace. */
static int mq_shared_det(nmod_mpoly_t result, nmod_mpoly_t **matrix, slong n,
                          const nmod_mpoly_ctx_t ctx, int parallel,
                          const mq_det_filter *filter, ulong choose[FLINT_BITS][FLINT_BITS])
{
    if (!mq_shared_admit(matrix, n, ctx, choose)) return 0;
#ifdef DRSOLVE_MQ_LAYOUT_TEST
    mq_shared_test_used++;
#endif
    if (g_dixon_verbose_level >= 2)
        printf("  MQ shared-index DP: size=%ld, workspace cap=%zu MiB\n", n,
               (size_t)DRSOLVE_MQ_SHARED_WORKSPACE_BYTES/(1024*1024));
    mq_shared_support prev; mq_shared_support_init(&prev, ctx, n+1);
    ulong *values = flint_malloc(sizeof(ulong)); values[0] = 1;
    for (slong k = 1; k <= n; k++) {
        const mq_det_filter *f = filter && k > filter->safe_linear_layers ? filter : NULL;
        slong count = k == n ? n : (slong)choose[n][k];
        mq_shared_support shifts, next;
        mq_shared_row(&shifts, matrix[n-k], n, &prev, ctx);
        uint32_t *map = mq_shared_next(&next, &prev, &shifts, f, ctx, 1);
        if (!next.count) {
            flint_free(map); flint_free(values);
            mq_shared_support_clear(&prev); mq_shared_support_clear(&shifts);
            mq_shared_support_clear(&next); nmod_mpoly_zero(result, ctx);
            return 1;
        }
        ulong *factors = flint_calloc(n*shifts.count, sizeof(ulong));
        nmod_mpoly_t tmp; nmod_mpoly_init(tmp, ctx);
        for (slong c = 0; c < n; c++) {
            int packed = nmod_mpoly_repack_bits(tmp, matrix[n-k][c], prev.bits, ctx);
            FLINT_ASSERT(packed); (void)packed;
            for (slong t = 0; t < tmp->length; t++) {
                uint32_t id = mq_shared_find(&shifts, tmp->exps+t*prev.words);
                FLINT_ASSERT(id != UINT32_MAX); factors[c*shifts.count+id] = tmp->coeffs[t];
            }
        }
        nmod_mpoly_clear(tmp, ctx);
        ulong *output = flint_calloc((size_t)count*next.count, sizeof(ulong));
        /* Each shift is injective: one output receives at most k*#shifts
         * products, not #child-monomials times that number. Negative factors
         * are represented by their canonical residues. Division checks avoid
         * overflow in the proof itself, including full-word primes. */
        ulong contributions = (k == n ? 1 : k)*shifts.count;
        ulong largest = ctx->mod.n-1;
        int lazy = contributions &&
                   largest <= (UWORD_MAX/contributions)/largest;
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
            ulong *dest = output+(size_t)index*next.count;
            for (slong j = 0; j < (k == n ? 1 : k); j++) {
                const ulong *source = values+(prefix[j]+suffix[j+1])*prev.count;
                for (size_t a = 0; a < shifts.count; a++) {
                    ulong scalar = factors[cols[j]*shifts.count+a];
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
            if (lazy) for (size_t i = 0; i < next.count; i++)
                NMOD_RED(dest[i], dest[i], ctx->mod);
        }
        if (k == n) {
            #pragma omp parallel for if(parallel) schedule(static)
            for (size_t i = 0; i < next.count; i++)
                for (slong j = 1; j < n; j++)
                    output[i] = nmod_add(output[i], output[j*next.count+i], ctx->mod);
            mq_shared_pack(result, output, &next, ctx);
        }
        flint_free(map); flint_free(factors); flint_free(values);
        mq_shared_support_clear(&prev); mq_shared_support_clear(&shifts);
        prev = next; values = output;
        if (k == n) { flint_free(values); mq_shared_support_clear(&prev); return 1; }
    }
    return 0;
}

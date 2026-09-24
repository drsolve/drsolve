/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Private, opt-in research backend. Included only in the layout benchmark.
 * Structural supports are Minkowski sums of row support unions, intersected
 * with the existing downward ideal. They never depend on sampled minor zeros.
 * Coefficients and transition indices are separate from shared packed keys. */
#include <stdint.h>
#include <assert.h>

static int mq_layout_enabled, mq_layout_profile, mq_layout_rotate, mq_layout_rotating;
void mq_layout_configure(int enabled, int profile, int rotate)
{
    mq_shared_test_used = 0;
    mq_shared_test_rank_layers = 0;
    mq_shared_test_variant = enabled >= 4 ? enabled-4 : 0;
    mq_shared_test_plan_seconds = mq_shared_test_arithmetic_seconds = 0;
    g_dixon_mq_step1_shared = enabled >= 4;
    mq_layout_enabled = enabled >= 4 ? 0 : enabled; mq_layout_profile = profile; mq_layout_rotate = rotate;
}
int mq_layout_rank_layers(void) { return mq_shared_test_rank_layers; }
int mq_layout_shared_calls(void) { return mq_shared_test_used; }
double mq_layout_plan_seconds(void) { return mq_shared_test_plan_seconds; }
double mq_layout_arithmetic_seconds(void) { return mq_shared_test_arithmetic_seconds; }

typedef struct {
    ulong *keys;
    uint32_t *table;
    size_t count, capacity, buckets;
    slong words;
    flint_bitcnt_t bits;
} mq_layout_support;

static size_t mq_layout_hash(const ulong *key, slong words)
{
    ulong h = 0x9e3779b9UL;
    for (slong w = 0; w < words; w++) h = mq_monom_hash(h ^ key[w]);
    return h;
}
static void mq_layout_rehash(mq_layout_support *s, size_t buckets)
{
    flint_free(s->table);
    s->buckets = buckets;
    s->table = flint_calloc(buckets, sizeof(uint32_t));
    for (size_t i = 0; i < s->count; i++) {
        size_t pos = mq_layout_hash(s->keys+i*s->words, s->words) & (buckets-1);
        while (s->table[pos]) pos = (pos+1) & (buckets-1);
        s->table[pos] = (uint32_t)i+1;
    }
}
static uint32_t mq_layout_find(const mq_layout_support *s, const ulong *key)
{
    size_t pos = mq_layout_hash(key, s->words) & (s->buckets-1);
    while (s->table[pos]) {
        uint32_t id = s->table[pos]-1;
        if (!memcmp(key, s->keys+(size_t)id*s->words, s->words*sizeof(ulong))) return id;
        pos = (pos+1) & (s->buckets-1);
    }
    return UINT32_MAX;
}
static uint32_t mq_layout_insert(mq_layout_support *s, const ulong *key)
{
    if (2*(s->count+1) >= s->buckets) mq_layout_rehash(s, 2*s->buckets);
    size_t pos = mq_layout_hash(key, s->words) & (s->buckets-1);
    while (s->table[pos]) {
        uint32_t id = s->table[pos]-1;
        if (!memcmp(key, s->keys+(size_t)id*s->words, s->words*sizeof(ulong))) return id;
        pos = (pos+1) & (s->buckets-1);
    }
    assert(s->count < UINT32_MAX);
    if (s->count == s->capacity) {
        s->capacity = s->capacity ? 2*s->capacity : 16;
        s->keys = flint_realloc(s->keys, s->capacity*s->words*sizeof(ulong));
    }
    memcpy(s->keys+s->count*s->words, key, s->words*sizeof(ulong));
    s->table[pos] = (uint32_t)s->count+1;
    return (uint32_t)s->count++;
}
static void mq_layout_empty(mq_layout_support *s, flint_bitcnt_t bits, slong words)
{
    memset(s, 0, sizeof(*s)); s->bits = bits; s->words = words;
    mq_layout_rehash(s, 32);
}
static void mq_layout_support_init(mq_layout_support *s, const nmod_mpoly_ctx_t ctx, slong degree)
{
    flint_bitcnt_t bits = mpoly_fix_bits(1+FLINT_BIT_COUNT((ulong)degree), ctx->minfo);
    slong words = mpoly_words_per_exp(bits, ctx->minfo);
    assert(words <= 3);
    mq_layout_empty(s, bits, words);
    ulong zero[3] = {0}; mq_layout_insert(s, zero);
}
static void mq_layout_support_clear(mq_layout_support *s)
{
    flint_free(s->keys); flint_free(s->table); memset(s, 0, sizeof(*s));
}
static size_t mq_layout_support_bytes(const mq_layout_support *s)
{
    return s->capacity*s->words*sizeof(ulong)+s->buckets*sizeof(uint32_t);
}
static void mq_layout_row(mq_layout_support *shifts, nmod_mpoly_t *row,
                          slong n, const mq_layout_support *previous,
                          const nmod_mpoly_ctx_t ctx)
{
    mq_layout_empty(shifts, previous->bits, previous->words);
    nmod_mpoly_t p; nmod_mpoly_init(p, ctx);
    for (slong col = 0; col < n; col++) {
        assert(nmod_mpoly_repack_bits(p, row[col], previous->bits, ctx));
        for (slong t = 0; t < p->length; t++)
            mq_layout_insert(shifts, p->exps+t*previous->words);
    }
    nmod_mpoly_clear(p, ctx);
}
static int mq_layout_keep(const ulong *key, const mq_det_filter *f,
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
static uint32_t *mq_layout_next(mq_layout_support *next, const mq_layout_support *prev,
                                const mq_layout_support *shifts, const mq_det_filter *f,
                                const nmod_mpoly_ctx_t ctx, int maps)
{
    mq_layout_empty(next, prev->bits, prev->words);
    uint32_t *map = maps ? flint_malloc(shifts->count*prev->count*sizeof(uint32_t)) : NULL;
    for (size_t a = 0; a < shifts->count; a++) for (size_t b = 0; b < prev->count; b++) {
        ulong key[3];
        for (slong w = 0; w < prev->words; w++)
            key[w] = shifts->keys[a*prev->words+w]+prev->keys[b*prev->words+w];
        uint32_t id = mq_layout_keep(key, f, ctx, prev->bits) ? mq_layout_insert(next, key) : UINT32_MAX;
        if (maps) map[a*prev->count+b] = id;
    }
    return map;
}

/* Profiling is deliberately excluded from each reported layer time. It
 * includes both an a-priori structural envelope and actual union occupancy. */
static void mq_layout_observe(mq_layout_support *prev, nmod_mpoly_t *polys, slong live,
                              nmod_mpoly_t *row, slong n, slong k,
                              const nmod_mpoly_ctx_t ctx, const mq_det_filter *filter,
                              double seconds)
{
    mq_layout_support shifts, next;
    mq_layout_row(&shifts, row, n, prev, ctx);
    mq_layout_next(&next, prev, &shifts, filter, ctx, 0);
    unsigned char *seen = flint_calloc(next.count, 1);
    size_t terms = 0, used = 0, allocated = 0;
    nmod_mpoly_t p; nmod_mpoly_init(p, ctx);
    for (slong j = 0; j < live; j++) {
        terms += polys[j]->length;
        allocated += (polys[j]->coeffs_alloc+polys[j]->exps_alloc)*sizeof(ulong);
        assert(nmod_mpoly_repack_bits(p, polys[j], prev->bits, ctx));
        for (slong t = 0; t < p->length; t++) {
            uint32_t id = mq_layout_find(&next, p->exps+t*prev->words);
            assert(id != UINT32_MAX);
            if (!seen[id]) { seen[id] = 1; used++; }
        }
    }
    printf("{\"kind\":\"layer\",\"backend\":\"sparse\",\"n\":%ld,\"quadratic_first\":%d,\"projected\":%d,\"k\":%ld,\"live\":%ld,\"seconds\":%.9f,\"terms\":%zu,\"structural_support\":%zu,\"actual_union\":%zu,\"occupancy\":%.9f,\"sparse_alloc_bytes\":%zu,\"dense_coeff_bytes\":%zu}\n",
           n, mq_layout_rotate, filter != NULL, k, live, seconds, terms, next.count, used,
           next.count ? (double)terms/(live*next.count) : 0, allocated, live*next.count*sizeof(ulong));
    nmod_mpoly_clear(p, ctx); flint_free(seen);
    mq_layout_support_clear(prev); mq_layout_support_clear(&shifts); *prev = next;
}

static void mq_layout_pack(nmod_mpoly_t out, const ulong *coeffs,
                           const mq_layout_support *s, const nmod_mpoly_ctx_t ctx)
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

/* Mode 1 retains the production root multiplication/reduction, mode 2 also
 * shares root indices, mode 3 additionally delays modular reduction when a
 * proven bound fits one limb. No table is cached across determinant calls. */
static int mq_layout_det(nmod_mpoly_t result, nmod_mpoly_t **matrix, slong n,
                          const nmod_mpoly_ctx_t ctx, int parallel,
                          const mq_det_filter *filter, ulong choose[FLINT_BITS][FLINT_BITS])
{
    mq_layout_support prev; mq_layout_support_init(&prev, ctx, n+1);
    ulong *values = flint_malloc(sizeof(ulong)); values[0] = 1;
    slong previous_count = 1;
    for (slong k = 1; k < n+(mq_layout_enabled >= 2); k++) {
        double start = omp_get_wtime();
        const mq_det_filter *f = filter && k > filter->safe_linear_layers ? filter : NULL;
        slong count = k == n ? n : (slong)choose[n][k];
        mq_layout_support shifts, next;
        mq_layout_row(&shifts, matrix[n-k], n, &prev, ctx);
        uint32_t *map = mq_layout_next(&next, &prev, &shifts, f, ctx, 1);
        ulong *factors = flint_calloc(n*shifts.count, sizeof(ulong));
        nmod_mpoly_t tmp; nmod_mpoly_init(tmp, ctx);
        for (slong c = 0; c < n; c++) {
            assert(nmod_mpoly_repack_bits(tmp, matrix[n-k][c], prev.bits, ctx));
            for (slong t = 0; t < tmp->length; t++) {
                uint32_t id = mq_layout_find(&shifts, tmp->exps+t*prev.words);
                assert(id != UINT32_MAX); factors[c*shifts.count+id] = tmp->coeffs[t];
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
        int lazy = mq_layout_enabled >= 3 && contributions &&
                   largest <= (UWORD_MAX/contributions)/largest;
        double setup = omp_get_wtime()-start;
        double arithmetic_start = omp_get_wtime();
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
        double pack = 0;
        if (k == n) {
            #pragma omp parallel for if(parallel) schedule(static)
            for (size_t i = 0; i < next.count; i++)
                for (slong j = 1; j < n; j++)
                    output[i] = nmod_add(output[i], output[j*next.count+i], ctx->mod);
            double phase = omp_get_wtime();
            mq_layout_pack(result, output, &next, ctx);
            pack = omp_get_wtime()-phase;
        }
        double arithmetic = omp_get_wtime()-arithmetic_start;
        if (mq_layout_profile) {
            size_t terms = 0;
            slong live = k == n ? 1 : count;
            for (size_t i = 0; i < (size_t)live*next.count; i++) terms += output[i] != 0;
            size_t bytes = ((size_t)previous_count*prev.count+(size_t)count*next.count+n*shifts.count)*sizeof(ulong)
                + shifts.count*prev.count*sizeof(uint32_t)+mq_layout_support_bytes(&prev)
                +mq_layout_support_bytes(&next)+mq_layout_support_bytes(&shifts);
            printf("{\"kind\":\"layer\",\"backend\":\"shared\",\"n\":%ld,\"quadratic_first\":%d,\"k\":%ld,\"live\":%ld,\"setup\":%.9f,\"arithmetic\":%.9f,\"pack\":%.9f,\"lazy\":%d,\"seconds\":%.9f,\"terms\":%zu,\"structural_support\":%zu,\"occupancy\":%.9f,\"workspace_bytes\":%zu}\n",
                n, mq_layout_rotate, k, live, setup, arithmetic-pack, pack, lazy, setup+arithmetic, terms, next.count,
                next.count ? (double)terms/(live*next.count) : 0, bytes);
        }
        flint_free(map); flint_free(factors); flint_free(values);
        mq_layout_support_clear(&prev); mq_layout_support_clear(&shifts);
        prev = next; values = output; previous_count = count;
        if (k == n) { flint_free(values); mq_layout_support_clear(&prev); return 1; }
    }
    double start = omp_get_wtime();
    nmod_mpoly_t *products = flint_malloc(n*sizeof(*products));
    for (slong j = 0; j < n; j++) nmod_mpoly_init(products[j], ctx);
    double pack_seconds = 0, mul_seconds = 0;
    #pragma omp parallel if(parallel) reduction(+:pack_seconds,mul_seconds)
    {
        nmod_mpoly_t child; nmod_mpoly_init(child, ctx);
        #pragma omp for schedule(static)
        for (slong j = 0; j < n; j++) {
            double phase = omp_get_wtime();
            mq_layout_pack(child, values+(size_t)(n-1-j)*prev.count, &prev, ctx);
            pack_seconds += omp_get_wtime()-phase;
            phase = omp_get_wtime();
            mq_filtered_mul(products[j], matrix[0][j], child, ctx, filter);
            if (j & 1) nmod_mpoly_neg(products[j], products[j], ctx);
            mul_seconds += omp_get_wtime()-phase;
        }
        for (slong stride = 1; stride < n; stride *= 2) {
            #pragma omp for schedule(static)
            for (slong j = 0; j < n; j += 2*stride)
                if (j+stride < n) nmod_mpoly_add(products[j], products[j], products[j+stride], ctx);
        }
        nmod_mpoly_clear(child, ctx);
    }
    nmod_mpoly_swap(result, products[0], ctx);
    if (mq_layout_profile)
        printf("{\"kind\":\"layer\",\"backend\":\"shared\",\"n\":%ld,\"quadratic_first\":%d,\"k\":%ld,\"seconds\":%.9f,\"pack_worker_seconds\":%.9f,\"mul_worker_seconds\":%.9f,\"terms\":%ld}\n",
            n, mq_layout_rotate, n, omp_get_wtime()-start, pack_seconds, mul_seconds, result->length);
    for (slong j = 0; j < n; j++) nmod_mpoly_clear(products[j], ctx);
    flint_free(products); flint_free(values); mq_layout_support_clear(&prev);
    return 1;
}

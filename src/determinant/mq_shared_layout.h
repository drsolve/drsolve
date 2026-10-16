/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Private MQ coefficient-array backend, included after the filter helpers. */
#include <stdint.h>

#ifdef DRSOLVE_MQ_LAYOUT_TEST
static int mq_shared_test_used;
static int mq_shared_test_rank_layers;
static int mq_shared_test_variant;
static double mq_shared_test_plan_seconds, mq_shared_test_arithmetic_seconds;
#endif

/* Test variants use flags 2 (sharding) and 4 (ordering); bit 1 is retired.
 * Variant 8 (benchmark mode 12) follows the retained production policy. */
static int mq_shared_variant(void)
{
#ifdef DRSOLVE_MQ_LAYOUT_TEST
    if (mq_shared_test_variant < 8) return mq_shared_test_variant;
#endif
    /* Sharding works with native multiword keys at every admitted size. */
    return 2;
}
static size_t mq_shared_parallel_threshold(void)
{
#ifdef DRSOLVE_MQ_LAYOUT_TEST
    if (mq_shared_test_variant < 8) return 32768;
#endif
    return 1000000;
}

typedef struct {
    ulong *keys;
    uint32_t *table; /* May be NULL after construction; later layers only need keys. */
    size_t count, capacity, buckets;
    slong words;
    flint_bitcnt_t bits;
    unsigned hash_shift;
    int sorted;
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
        size_t pos = (mq_shared_hash(s->keys+i*s->words, s->words) >> s->hash_shift) & (buckets-1);
        while (s->table[pos]) pos = (pos+1) & (buckets-1);
        s->table[pos] = (uint32_t)i+1;
    }
}
static uint32_t mq_shared_find(const mq_shared_support *s, const ulong *key)
{
    size_t pos = (mq_shared_hash(key, s->words) >> s->hash_shift) & (s->buckets-1);
    while (s->table[pos]) {
        uint32_t id = s->table[pos]-1;
        if (!memcmp(key, s->keys+(size_t)id*s->words, s->words*sizeof(ulong))) return id;
        pos = (pos+1) & (s->buckets-1);
    }
    return UINT32_MAX;
}
static uint32_t mq_shared_insert_hashed(mq_shared_support *s, const ulong *key, size_t hash)
{
    if (2*(s->count+1) >= s->buckets) mq_shared_rehash(s, 2*s->buckets);
    size_t pos = (hash >> s->hash_shift) & (s->buckets-1);
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
static uint32_t mq_shared_insert(mq_shared_support *s, const ulong *key)
{
    return mq_shared_insert_hashed(s, key, mq_shared_hash(key, s->words));
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
        int packed = nmod_mpoly_repack_bits(p, row[col], mpoly_fix_bits(previous->bits, ctx->minfo), ctx);
        FLINT_ASSERT(packed); (void)packed;
        for (slong t = 0; t < p->length; t++) {
            mq_shared_insert(shifts, p->exps+t*previous->words);
        }
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
                                const nmod_mpoly_ctx_t ctx, int maps, int parallel)
{
    mq_shared_empty(next, prev->bits, prev->words);
    uint32_t *map = maps ? flint_malloc(shifts->count*prev->count*sizeof(uint32_t)) : NULL;
    unsigned workers = 1, hash_shift = 0;
#ifdef _OPENMP
    if ((mq_shared_variant() & 2) && parallel && maps &&
        shifts->count*prev->count >= mq_shared_parallel_threshold() && !omp_in_parallel()) {
        unsigned limit = FLINT_MIN(8, omp_get_max_threads());
        while (2*workers <= limit) { workers *= 2; hash_shift++; }
    }
#endif
    if (workers > 1 && shifts->count*prev->count <= UINT32_MAX/workers) {
        mq_shared_support shards[8];
        uint32_t offsets[8];
        for (unsigned owner = 0; owner < workers; owner++) {
            mq_shared_empty(&shards[owner], prev->bits, prev->words);
            shards[owner].hash_shift = hash_shift;
        }
        /* A key belongs to exactly one shard. No locks and no duplicate keys
         * across shards. Each map cell also has exactly one writer. Scanning
         * sums per shard trades cheap packed additions for synchronization. */
        #pragma omp parallel for num_threads(workers) schedule(static)
        for (unsigned owner = 0; owner < workers; owner++) {
            for (size_t a = 0; a < shifts->count; a++) for (size_t b = 0; b < prev->count; b++) {
                ulong key[3];
                for (slong w = 0; w < prev->words; w++)
                    key[w] = shifts->keys[a*prev->words+w]+prev->keys[b*prev->words+w];
                size_t hash = mq_shared_hash(key, prev->words);
                if ((hash & (workers-1)) != owner) continue;
                uint32_t id = UINT32_MAX;
                if (mq_shared_keep(key, f, ctx, prev->bits))
                    id = mq_shared_insert_hashed(&shards[owner], key, hash)*workers+owner;
                map[a*prev->count+b] = id;
            }
        }
        /* Previous layers only require keys, not their construction hash
         * tables. Concatenate disjoint supports and remap encoded local IDs. */
        flint_free(next->table); next->table = NULL; next->buckets = 0;
        for (unsigned owner = 0; owner < workers; owner++) {
            offsets[owner] = next->count; next->count += shards[owner].count;
        }
        next->capacity = next->count;
        next->keys = flint_malloc(next->count*next->words*sizeof(ulong));
        for (unsigned owner = 0; owner < workers; owner++) {
            if (shards[owner].count)
                memcpy(next->keys+(size_t)offsets[owner]*next->words, shards[owner].keys,
                       shards[owner].count*next->words*sizeof(ulong));
            mq_shared_support_clear(&shards[owner]);
        }
        #pragma omp parallel for num_threads(workers) schedule(static)
        for (size_t i = 0; i < shifts->count*prev->count; i++)
            if (map[i] != UINT32_MAX) map[i] = offsets[map[i] & (workers-1)]+(map[i] >> hash_shift);
        return map;
    }
    for (size_t a = 0; a < shifts->count; a++) for (size_t b = 0; b < prev->count; b++) {
        ulong key[3];
        for (slong w = 0; w < prev->words; w++)
            key[w] = shifts->keys[a*prev->words+w]+prev->keys[b*prev->words+w];
        uint32_t id = mq_shared_keep(key, f, ctx, prev->bits) ? mq_shared_insert(next, key) : UINT32_MAX;
        if (maps) map[a*prev->count+b] = id;
    }
    return map;
}

#ifdef DRSOLVE_MQ_LAYOUT_TEST
/* Stable byte radix sort of indices. Sorting the common basis makes every
 * fixed-monomial transition monotone, improving coefficient-array locality.
 * Construction hash tables are no longer needed after transitions exist. */
static void mq_shared_order(mq_shared_support *s, uint32_t *map, size_t entries, int parallel)
{
    if (!s->count) return;
    uint32_t *order = flint_malloc(s->count*sizeof(uint32_t));
    uint32_t *scratch = flint_malloc(s->count*sizeof(uint32_t));
    uint32_t *inverse = flint_malloc(s->count*sizeof(uint32_t));
    ulong used[3] = {0};
    for (size_t i = 0; i < s->count; i++) {
        order[i] = i;
        for (slong w = 0; w < s->words; w++) used[w] |= s->keys[i*s->words+w];
    }
    for (slong w = 0; w < s->words; w++) for (unsigned shift = 0; shift < FLINT_BITS; shift += 8) {
        if (!((used[w] >> shift) & 255)) continue;
        size_t counts[256] = {0}, positions[256];
        for (size_t i = 0; i < s->count; i++)
            counts[(s->keys[(size_t)order[i]*s->words+w] >> shift) & 255]++;
        positions[0] = 0;
        for (unsigned b = 1; b < 256; b++) positions[b] = positions[b-1]+counts[b-1];
        for (size_t i = 0; i < s->count; i++) {
            unsigned b = (s->keys[(size_t)order[i]*s->words+w] >> shift) & 255;
            scratch[positions[b]++] = order[i];
        }
        uint32_t *swap = order; order = scratch; scratch = swap;
    }
    ulong *keys = flint_malloc(s->count*s->words*sizeof(ulong));
    for (size_t i = 0; i < s->count; i++) {
        uint32_t old = order[s->count-1-i];
        inverse[old] = i;
        memcpy(keys+i*s->words, s->keys+(size_t)old*s->words, s->words*sizeof(ulong));
    }
    #pragma omp parallel for if(parallel && entries >= 32768) schedule(static)
    for (size_t i = 0; i < entries; i++)
        if (map[i] != UINT32_MAX) map[i] = inverse[map[i]];
    flint_free(order); flint_free(scratch); flint_free(inverse);
    flint_free(s->keys); flint_free(s->table);
    s->keys = keys; s->capacity = s->count; s->table = NULL; s->buckets = 0; s->sorted = 1;
}
#endif

/* Saturating total-degree simplex size. Index/addressability rejection
 * leaves the caller's output untouched. */
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
/* No fixed memory budget by default. A build may still override this.
 * SIZE_MAX retains addressability checks; uint32 index limits remain below. */
#define DRSOLVE_MQ_SHARED_WORKSPACE_BYTES SIZE_MAX
#endif

static int mq_shared_admit_coeffs(nmod_mpoly_t **matrix, slong n,
                           const nmod_mpoly_ctx_t ctx,
                           ulong choose[FLINT_BITS][FLINT_BITS], size_t coefficient_bytes)
{
    /* At least a limb per slot also covers the parallel plan construction. */
    coefficient_bytes = FLINT_MAX(coefficient_bytes, sizeof(ulong));
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
        if (count > budget/coefficient_bytes/current ||
            previous_count > budget/coefficient_bytes/previous) ok = 0;
        if (ok) {
            /* The plan is built BEFORE current coefficients exist. Disjoint
             * shard capacities/hash tables plus concatenation use at most
             * 3*current*words limbs + 4*current uint32 slots (and metadata).
             * The serial support bound plus current coefficient allowance
             * dominates this: count>=4 and words<=3. Reserve 64 KiB for the
             * <=8 shards' small tables/capacity rounding and stack metadata.
             * Compact keys only decrease this native-packing estimate. */
            long double bytes = 65536.0L+((long double)count*current
                +(long double)previous_count*previous+(long double)n*shifts.count)*coefficient_bytes
                +(long double)shifts.count*previous*sizeof(uint32_t)
                +mq_shared_support_bound(current, words)
                +mq_shared_support_bound(previous, words)
                +mq_shared_support_bytes(&shifts);
            /* Output packing allowance; FLINT sort scratch, input storage and
             * allocator overhead are not a process-RSS guarantee. */
            if (k == n) bytes += (long double)current*(words*sizeof(ulong)+coefficient_bytes);
            if (bytes > budget) ok = 0;
        }
        mq_shared_support_clear(&shifts);
        previous = current; previous_count = count;
    }
    mq_shared_support_clear(&packing);
    return ok;
}

/* Retain the independent limb-sized nmod policy. */
static int mq_shared_admit(nmod_mpoly_t **matrix, slong n,
    const nmod_mpoly_ctx_t ctx, ulong choose[FLINT_BITS][FLINT_BITS])
{
    return mq_shared_admit_coeffs(matrix, n, ctx, choose, sizeof(ulong));
}

#include "mq_rank_layout.h"

/* Admission bounds n and the variable count before this predicate is used.
 * A degree <= 2 row has at most (nv+1)(nv+2)/2 distinct shifts, each injective
 * on the previous support. At most n children contribute to one minor.
 * All source coefficients and signed factors are canonical residues, so this
 * conservative bound certifies the entire sum, not only one multiplication. */
static int mq_shared_narrow_safe(slong n, ulong prime)
{
    if (FLINT_BITS != 64 || n < 4 || n > (FLINT_BITS-1)/2 || prime < 2) return 0;
    ulong nv = 2*n-1, largest = prime-1;
    ulong products = (ulong)n*(nv+1)*(nv+2)/2;
    return largest <= (UINT32_MAX/products)/largest;
}

#define MQ_SHARED_COEFF ulong
#define MQ_SHARED_PACK mq_shared_pack_word
#define MQ_SHARED_DET mq_shared_det_word
#define MQ_SHARED_NARROW 0
#include "mq_shared_coeff_template.h"

#define MQ_SHARED_COEFF uint32_t
#define MQ_SHARED_PACK mq_shared_pack_u32
#define MQ_SHARED_DET mq_shared_det_u32
#define MQ_SHARED_NARROW 1
#include "mq_shared_coeff_template.h"

static int mq_shared_det(nmod_mpoly_t result, nmod_mpoly_t **matrix, slong n,
                          const nmod_mpoly_ctx_t ctx, int parallel,
                          const mq_det_filter *filter, ulong choose[FLINT_BITS][FLINT_BITS])
{
    if (!mq_shared_admit(matrix, n, ctx, choose)) return 0;
    if (mq_shared_narrow_safe(n, ctx->mod.n))
        return mq_shared_det_u32(result, matrix, n, ctx, parallel, filter, choose);
    return mq_shared_det_word(result, matrix, n, ctx, parallel, filter, choose);
}

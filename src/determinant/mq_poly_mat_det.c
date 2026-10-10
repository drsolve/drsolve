/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_poly_mat_det.h"
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Polynomial products use coarse independent row panels, not a parallel
 * region for every scalar pivot. All views have one writer. */
static void mq_block_mul(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B)
{
    slong panels = 1;
#ifdef _OPENMP
    if (!omp_in_parallel())
        panels = FLINT_MIN(omp_get_max_threads(), FLINT_MAX(1, A->r / 32));
#pragma omp parallel for schedule(static) if (panels > 1) num_threads(panels)
#endif
    for (slong p = 0; p < panels; p++)
    {
        slong begin = A->r * p / panels, end = A->r * (p + 1) / panels;
        nmod_poly_mat_t a, c;
        nmod_poly_mat_window_init(a, A, begin, 0, end, A->c);
        nmod_poly_mat_window_init(c, C, begin, 0, end, C->c);
        nmod_poly_mat_mul(c, a, B);
        nmod_poly_mat_window_clear(c);
        nmod_poly_mat_window_clear(a);
    }
}

static int mq_permutation_odd(const slong *p, slong n, slong offset)
{
    unsigned char *seen = flint_calloc((size_t)n, 1);
    int odd = 0;
    for (slong i = 0; i < n; i++)
        if (!seen[i])
        {
            slong j = i, len = 0;
            do
            {
                seen[j] = 1;
                len++;
                j = p[j] - offset;
            } while (!seen[j]);
            odd ^= (len - 1) & 1;
        }
    flint_free(seen);
    return odd;
}

typedef struct
{
    slong begin, size;
    nmod_mat_t lu;
    slong *perm;
} mq_constant_block;

static int mq_schur_valid(const nmod_poly_mat_t core, const ulong *factor,
                          const nmod_poly_mat_t B, const slong *rd,
                          const slong *cd, slong h, slong sigma)
{
    slong n = B->r, budget = 0;
    if (!factor || !rd || !cd || n != B->c || h < 0 || h > n || sigma < 0 || core == B ||
        core->r != h || core->c != h || core->modulus != B->modulus)
        return 0;
    for (slong i = 0; i < n; i++)
        if (rd[i] < 0 || cd[i] < 0 || rd[i] > sigma || cd[i] > sigma)
            return 0;
    for (slong i = 0; i < n; i++) {
        if (i >= h) budget += sigma - rd[i] - cd[i];
        for (slong j = 0; j < n; j++)
            if (!nmod_poly_is_zero(nmod_poly_mat_entry(B, i, j)) &&
                nmod_poly_degree(nmod_poly_mat_entry(B, i, j)) > sigma - rd[i] - cd[j])
                return 0;
    }
    return budget == 0;
}

int nmod_poly_mat_mq_schur(nmod_poly_mat_t core, ulong *factor, const nmod_poly_mat_t B,
                           const slong *rd, const slong *cd, slong h, slong sigma)
{
    slong n = B->r;
    if (!mq_schur_valid(core, factor, B, rd, cd, h, sigma)) return 0;
    slong e = n - h;

    slong *rows = flint_malloc((size_t)e * sizeof(slong));
    slong *cols = flint_malloc((size_t)e * sizeof(slong));
    /* At most e nonempty groups. No full polynomial matrix is copied. */
    mq_constant_block *blocks = flint_malloc((size_t)e * sizeof(*blocks));
    slong nr = 0, nc = 0, groups = 0;
    int ok = 1;
    nmod_t mod;
    nmod_init(&mod, B->modulus);
    ulong scale = 1;
    for (slong d = 0; d <= sigma && nr < e; d++)
    {
        slong begin = nr, cb = nc;
        for (slong i = h; i < n; i++)
            if (rd[i] == d)
                rows[nr++] = i;
        for (slong j = h; j < n; j++)
            if (cd[j] == sigma - d)
                cols[nc++] = j;
        if (nr - begin != nc - cb)
        {
            ok = 0;
            break;
        }
        slong k = nr - begin;
        if (!k)
            continue;
        mq_constant_block *g = blocks + groups++;
        g->begin = begin;
        g->size = k;
        g->perm = flint_malloc((size_t)k * sizeof(slong));
        nmod_mat_init(g->lu, k, k, B->modulus);
        for (slong i = 0; i < k; i++)
            for (slong j = 0; j < k; j++)
                nmod_mat_entry(g->lu, i, j) = nmod_poly_get_coeff_ui(
                    nmod_poly_mat_entry(B, rows[begin + i], cols[begin + j]), 0);
        /* Certify every constant diagonal block before polynomial work. */
        if (nmod_mat_lu(g->perm, g->lu, 1) != k)
        {
            ok = 0;
            break;
        }
        for (slong i = 0; i < k; i++)
            scale = nmod_mul(scale, nmod_mat_entry(g->lu, i, i), mod);
        if (mq_permutation_odd(g->perm, k, 0))
            scale = nmod_neg(scale, mod);
    }
    if (nr != e || nc != e)
        ok = 0;
    if (ok)
    {
        if (mq_permutation_odd(rows, e, h) ^ mq_permutation_odd(cols, e, h))
            scale = nmod_neg(scale, mod);
        /* E X = V, solved from the highest degree block down. X has only
         * e*h entries. E and the input B remain immutable throughout. */
        nmod_poly_mat_t X;
        nmod_poly_mat_init(X, e, h, B->modulus);
        for (slong b = groups; b-- > 0;)
        {
            mq_constant_block *g = blocks + b;
            slong start = g->begin, k = g->size, end = start + k;
            nmod_poly_mat_t rhs;
            nmod_poly_mat_init(rhs, k, h, B->modulus);
            if (end < e)
            {
                nmod_poly_mat_t off, tail;
                nmod_poly_mat_init(off, k, e - end, B->modulus);
                for (slong i = 0; i < k; i++)
                    for (slong j = end; j < e; j++)
                        nmod_poly_set(nmod_poly_mat_entry(off, i, j - end),
                                      nmod_poly_mat_entry(B, rows[start + i], cols[j]));
                nmod_poly_mat_window_init(tail, X, end, 0, e, h);
                mq_block_mul(rhs, off, tail);
                nmod_poly_mat_window_clear(tail);
                nmod_poly_mat_clear(off);
            }
            slong length = 0;
            for (slong i = 0; i < k; i++)
                for (slong j = 0; j < h; j++)
                {
                    nmod_poly_struct *p = nmod_poly_mat_entry(rhs, i, j);
                    nmod_poly_sub(p, nmod_poly_mat_entry(B, rows[start + i], j), p);
                    length = FLINT_MAX(length, p->length);
                }
            if (length)
            {
                /* Reuse the same LU for all parameter coefficients. Pack
                 * them as multiple right-hand sides for blocked field solves. */
                nmod_mat_t packed, lower, solved;
                nmod_mat_init(packed, k, h * length, B->modulus);
                nmod_mat_init(lower, k, h * length, B->modulus);
                nmod_mat_init(solved, k, h * length, B->modulus);
                for (slong i = 0; i < k; i++)
                    for (slong j = 0; j < h; j++)
                    {
                        const nmod_poly_struct *p = nmod_poly_mat_entry(rhs, g->perm[i], j);
                        for (slong t = 0; t < p->length; t++)
                            nmod_mat_entry(packed, i, t * h + j) = p->coeffs[t];
                    }
                nmod_mat_solve_tril(lower, g->lu, packed, 1);
                nmod_mat_solve_triu(solved, g->lu, lower, 0);
                for (slong i = 0; i < k; i++)
                    for (slong j = 0; j < h; j++)
                        for (slong t = 0; t < length; t++)
                            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(X, start + i, j), t,
                                                   nmod_mat_entry(solved, i, t * h + j));
                nmod_mat_clear(solved);
                nmod_mat_clear(lower);
                nmod_mat_clear(packed);
            }
            nmod_poly_mat_clear(rhs);
        }
        nmod_poly_mat_t U, product;
        nmod_poly_mat_init(U, h, e, B->modulus);
        nmod_poly_mat_init(product, h, h, B->modulus);
        for (slong i = 0; i < h; i++)
            for (slong j = 0; j < e; j++)
                nmod_poly_set(nmod_poly_mat_entry(U, i, j), nmod_poly_mat_entry(B, i, cols[j]));
        mq_block_mul(product, U, X);
        for (slong i = 0; i < h; i++)
            for (slong j = 0; j < h; j++)
                nmod_poly_sub(nmod_poly_mat_entry(core, i, j), nmod_poly_mat_entry(B, i, j),
                              nmod_poly_mat_entry(product, i, j));
        nmod_poly_mat_clear(product);
        nmod_poly_mat_clear(U);
        nmod_poly_mat_clear(X);
        *factor = scale;
    }
    for (slong b = 0; b < groups; b++)
    {
        nmod_mat_clear(blocks[b].lu);
        flint_free(blocks[b].perm);
    }
    flint_free(blocks);
    flint_free(cols);
    flint_free(rows);
    return ok;
}

/* Assign pivot indices to complement slots and unused indices to core slots.
 * Only indices of this exact degree move, so all entry degree bounds survive. */
static void mq_place_degree_pivots(slong *map, const slong *weights, slong n,
    slong h, slong degree, const slong *indices, slong count,
    const slong *pivots, slong k)
{
    unsigned char *used = flint_calloc((size_t)count, 1);
    for (slong i = 0; i < k; i++) used[pivots[i]] = 1;
    slong a = 0, b = 0;
    for (slong i = 0; i < n; i++) if (weights[i] == degree) {
        if (i >= h) map[i] = indices[pivots[b++]];
        else {
            while (a < count && used[a]) a++;
            map[i] = indices[a++];
        }
    }
    flint_free(used);
}

int nmod_poly_mat_mq_schur_repair(nmod_poly_mat_t core, ulong *factor,
    const nmod_poly_mat_t B, const slong *rd, const slong *cd, slong h, slong sigma)
{
    if (!mq_schur_valid(core, factor, B, rd, cd, h, sigma)) return 0;
    slong n = B->r;
    slong *rows = flint_malloc((size_t)n * sizeof(slong));
    slong *cols = flint_malloc((size_t)n * sizeof(slong));
    slong *ri = flint_malloc((size_t)n * sizeof(slong));
    slong *ci = flint_malloc((size_t)n * sizeof(slong));
    slong *perm = flint_malloc((size_t)n * sizeof(slong));
    slong *pivots = flint_malloc((size_t)n * sizeof(slong));
    int ok = 1;
    for (slong i = 0; i < n; i++) rows[i] = cols[i] = i;
    for (slong d = 0; d <= sigma && ok; d++) {
        slong nr = 0, nc = 0, kr = 0, kc = 0;
        for (slong i = 0; i < n; i++) {
            if (rd[i] == d) { ri[nr++] = i; kr += i >= h; }
            if (cd[i] == sigma - d) { ci[nc++] = i; kc += i >= h; }
        }
        if (kr != kc) { ok = 0; break; }
        if (!kr) continue;
        nmod_mat_t block;
        nmod_mat_init(block, nr, nc, B->modulus);
        for (slong i = 0; i < nr; i++) for (slong j = 0; j < nc; j++)
            nmod_mat_entry(block, i, j) = nmod_poly_get_coeff_ui(
                nmod_poly_mat_entry(B, ri[i], ci[j]), 0);
        slong rank = nmod_mat_lu(perm, block, 0);
        if (rank < kr) ok = 0;
        else {
            slong next = 0;
            for (slong i = 0; i < kr; i++) {
                while (next < nc && !nmod_mat_entry(block, i, next)) next++;
                pivots[i] = next++;
            }
            mq_place_degree_pivots(rows, rd, n, h, d, ri, nr, perm, kr);
            mq_place_degree_pivots(cols, cd, n, h, sigma - d, ci, nc, pivots, kr);
        }
        nmod_mat_clear(block);
    }
    if (ok) {
        nmod_poly_mat_t reordered;
        nmod_poly_mat_init(reordered, n, n, B->modulus);
        for (slong i = 0; i < n; i++) for (slong j = 0; j < n; j++)
            nmod_poly_set(nmod_poly_mat_entry(reordered, i, j),
                          nmod_poly_mat_entry(B, rows[i], cols[j]));
        ulong scale;
        ok = nmod_poly_mat_mq_schur(core, &scale, reordered, rd, cd, h, sigma);
        if (ok) {
            if (mq_permutation_odd(rows, n, 0) ^ mq_permutation_odd(cols, n, 0))
                scale = B->modulus - scale;
            *factor = scale;
        }
        nmod_poly_mat_clear(reordered);
    }
    flint_free(pivots); flint_free(perm); flint_free(ci); flint_free(ri);
    flint_free(cols); flint_free(rows);
    return ok;
}

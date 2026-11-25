/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "dixon_exterior.h"
#include <flint/nmod_mpoly.h>
#include <flint/ulong_extras.h>
#include <stdlib.h>
#include <string.h>

#define PRIME 65537

typedef struct {
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_struct *diff;
    slong n, r;
    flint_rand_t rng;
} exterior_t;

/* At z=0 a monomial contributes to exactly one divided difference. */
static void top_terms(exterior_t *e, slong col, ulong *exp, slong pos, slong left)
{
    if (pos + 1 < e->r) {
        for (slong a = 0; a <= left; ++a) {
            exp[pos] = a;
            top_terms(e, col, exp, pos + 1, left - a);
        }
        return;
    }
    exp[pos] = left;
    while (!exp[pos])
        --pos;
    --exp[pos];
    nmod_mpoly_set_coeff_ui_ui(e->diff + pos * e->n + col, 1 + n_randint(e->rng, PRIME - 1), exp,
                               e->ctx);
    ++exp[pos];
}

static void minor(nmod_mpoly_t out, exterior_t *e, slong row, slong *cols, slong k)
{
    if (k == 1) {
        nmod_mpoly_set(out, e->diff + row * e->n + cols[0], e->ctx);
        return;
    }
    nmod_mpoly_t tmp;
    nmod_mpoly_init(tmp, e->ctx);
    slong rest[k - 1];
    for (slong j = 0; j < k; ++j) {
        if (nmod_mpoly_is_zero(e->diff + row * e->n + cols[j], e->ctx))
            continue;
        for (slong a = 0, b = 0; a < k; ++a)
            if (a != j)
                rest[b++] = cols[a];
        nmod_mpoly_zero(tmp, e->ctx);
        minor(tmp, e, row + 1, rest, k - 1);
        nmod_mpoly_mul(tmp, tmp, e->diff + row * e->n + cols[j], e->ctx);
        if (j % 2)
            nmod_mpoly_sub(out, out, tmp, e->ctx);
        else
            nmod_mpoly_add(out, out, tmp, e->ctx);
    }
    nmod_mpoly_clear(tmp, e->ctx);
}

/* Compress the monomials actually present into consecutive column indices. */
static slong top_rank(exterior_t *e, const slong *d, slong s)
{
    slong r = e->r, n = e->n, I[r], rows = 0;
    ulong exp[r];
    nmod_mpoly_struct *polys = NULL;
    nmod_mpoly_t support;
    nmod_mpoly_init(support, e->ctx);
    for (slong j = 0; j < r; ++j)
        I[j] = j;
    do {
        slong sum = 0;
        for (slong j = 0; j < r; ++j)
            sum += d[I[j]] - 1;
        if (sum != s)
            continue;
        polys = flint_realloc(polys, (rows + 1) * sizeof(*polys));
        nmod_mpoly_init(polys + rows, e->ctx);
        minor(polys + rows, e, 0, I, r);
        for (slong t = 0; t < nmod_mpoly_length(polys + rows, e->ctx); ++t) {
            nmod_mpoly_get_term_exp_ui(exp, polys + rows, t, e->ctx);
            nmod_mpoly_set_coeff_ui_ui(support, 1, exp, e->ctx);
        }
        ++rows;
    } while (next_subset(I, r, n));
    slong cols = nmod_mpoly_length(support, e->ctx);
    nmod_mat_t mat;
    nmod_mat_init(mat, rows, cols, PRIME);
    for (slong j = 0; j < cols; ++j) {
        nmod_mpoly_get_term_exp_ui(exp, support, j, e->ctx);
        for (slong i = 0; i < rows; ++i)
            nmod_mat_entry(mat, i, j) = nmod_mpoly_get_coeff_ui_ui(polys + i, exp, e->ctx);
    }
    slong rank = nmod_mat_rank(mat);
    nmod_mat_clear(mat);
    for (slong i = 0; i < rows; ++i)
        nmod_mpoly_clear(polys + i, e->ctx);
    flint_free(polys);
    nmod_mpoly_clear(support, e->ctx);
    return rank;
}

/* Column group J is wedged with a two-row minor. Each row group I=J U {i,j}
   is indexed by the exponent of z1, each column by the exponent of w1. */
static slong wedge_rank(exterior_t *e, const slong *d, slong r, slong s)
{
    slong n = e->n, I[r], J[r - 2], rows = 0, cols = 0;
    for (slong j = 0; j < r; ++j)
        I[j] = j;
    do {
        slong rem = -s;
        for (slong j = 0; j < r; ++j)
            rem += d[I[j]] - 1;
        if (rem >= 0)
            rows += rem + 1;
    } while (next_subset(I, r, n));
    for (slong j = 0; j < r - 2; ++j)
        J[j] = j;
    do {
        slong w = s;
        for (slong j = 0; j < r - 2; ++j)
            w -= d[J[j]] - 1;
        if (w >= 0)
            cols += w + 1;
    } while (next_subset(J, r - 2, n));
    nmod_mat_t mat;
    nmod_mat_init(mat, rows, cols, PRIME);
    /* Each pair occurs with many different J. Expand its minor only once. */
    nmod_mpoly_struct *pairs = flint_calloc(n, sizeof(*pairs) * n);
    for (slong i = 0; i < n; ++i) {
        for (slong j = i + 1; j < n; ++j) {
            slong ij[] = {i, j};
            nmod_mpoly_init(pairs + i * n + j, e->ctx);
            minor(pairs + i * n + j, e, 0, ij, 2);
        }
    }
    slong row = 0;
    for (slong j = 0; j < r; ++j)
        I[j] = j;
    do {
        slong rem = -s;
        for (slong j = 0; j < r; ++j)
            rem += d[I[j]] - 1;
        if (rem < 0)
            continue;
        slong col = 0;
        for (slong j = 0; j < r - 2; ++j)
            J[j] = j;
        do {
            slong w = s, pair[2], missing = 0, parity = 0;
            for (slong j = 0; j < r - 2; ++j)
                w -= d[J[j]] - 1;
            if (w < 0)
                continue;
            for (slong a = 0, b = 0; a < r; ++a) {
                if (b < r - 2 && I[a] == J[b]) {
                    ++b;
                    continue;
                }
                if (missing < 2)
                    pair[missing] = I[a];
                ++missing;
                parity ^= b % 2;
            }
            if (missing == 2) {
                nmod_mpoly_struct *poly = pairs + pair[0] * n + pair[1];
                for (slong t = 0; t < nmod_mpoly_length(poly, e->ctx); ++t) {
                    ulong exp[4];
                    nmod_mpoly_get_term_exp_ui(exp, poly, t, e->ctx);
                    if (exp[2] + exp[3] != (ulong)w)
                        continue;
                    ulong c = nmod_mpoly_get_term_coeff_ui(poly, t, e->ctx);
                    nmod_mat_entry(mat, row + exp[0], col + exp[2]) = parity ? PRIME - c : c;
                }
            }
            col += w + 1;
        } while (next_subset(J, r - 2, n));
        row += rem + 1;
    } while (next_subset(I, r, n));
    for (slong i = 0; i < n; ++i)
        for (slong j = i + 1; j < n; ++j)
            nmod_mpoly_clear(pairs + i * n + j, e->ctx);
    flint_free(pairs);
    slong rank = nmod_mat_rank(mat);
    nmod_mat_clear(mat);
    return rank;
}

slong dixon_exterior_rank(const slong *d, slong n, slong r, slong s, int top)
{
    exterior_t e;
    e.r = top ? r : 2;
    e.n = n;
    flint_rand_init(e.rng);
    nmod_mpoly_ctx_init(e.ctx, top ? r : 4, ORD_LEX, PRIME);
    e.diff = flint_malloc(n * e.r * sizeof(*e.diff));
    for (slong j = 0; j < n * e.r; ++j)
        nmod_mpoly_init(e.diff + j, e.ctx);
    for (slong i = 0; i < n; ++i) {
        if (top) {
            ulong exp[r];
            top_terms(&e, i, exp, 0, d[i]);
        } else {
            for (slong k = 0; k <= d[i]; ++k) {
                ulong c = 1 + n_randint(e.rng, PRIME - 1);
                for (slong a = 0; a < k; ++a) {
                    ulong exp[] = {a, d[i] - k, k - 1 - a, 0};
                    nmod_mpoly_set_coeff_ui_ui(e.diff + i, c, exp, e.ctx);
                }
                for (slong a = 0; a < d[i] - k; ++a) {
                    ulong exp[] = {0, a, k, d[i] - k - 1 - a};
                    nmod_mpoly_set_coeff_ui_ui(e.diff + n + i, c, exp, e.ctx);
                }
            }
        }
    }
    slong rank = top ? top_rank(&e, d, s) : wedge_rank(&e, d, r, s);
    for (slong j = 0; j < n * e.r; ++j)
        nmod_mpoly_clear(e.diff + j, e.ctx);
    flint_free(e.diff);
    nmod_mpoly_ctx_clear(e.ctx);
    flint_rand_clear(e.rng);
    return rank;
}

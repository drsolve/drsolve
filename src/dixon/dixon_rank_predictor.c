/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Generic Dixon rank prediction by midpoint reduction and boundary filtration. */
#include "dixon_exterior.h"
#include "dixon_rank_predictor.h"
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct system_data {
    slong n, N, max, *d, *cap, *block, **bounds;
    fmpz *A, *C;
    struct system_data *sub;
} system_t_data;

static size_t checked_mul(size_t a, size_t b)
{
    if (b && a > SIZE_MAX / b) {
        flint_throw(FLINT_OVERFLOW, "Dixon rank predictor allocation size overflow");
    }
    return a * b;
}

static slong min_s(slong a, slong b) { return a < b ? a : b; }
static slong max_s(slong a, slong b) { return a > b ? a : b; }
static slong sat_add(slong a, slong b) { return a > WORD_MAX - b ? WORD_MAX : a + b; }
static int compare_long(const void *a, const void *b)
{
    slong x = *(const slong *)a, y = *(const slong *)b;
    return (x > y) - (x < y);
}

static system_t_data *system_new(slong *d, slong n)
{
    system_t_data *h = flint_calloc(1, sizeof(*h));
    h->d = d;
    h->n = n;
    h->cap = flint_calloc(n, sizeof(slong));
    for (slong j = 1; j < n; ++j)
        h->cap[j] = h->cap[j - 1] + d[n - j] - 1;
    h->max = h->cap[n - 1];
    h->N = h->max + (d[0] - 1);
    h->C = _fmpz_vec_init(h->N + 1);
    fmpz_one(h->C);
    slong deg = 0;
    for (slong i = 0; i < n; ++i) {
        deg += d[i] - 1;
        /* Multiply by (1-x^d)/(1-x), in place. */
        for (slong t = deg; t >= d[i]; --t)
            fmpz_sub(h->C + t, h->C + t, h->C + t - d[i]);
        for (slong t = 1; t <= deg; ++t)
            fmpz_add(h->C + t, h->C + t, h->C + t - 1);
    }
    h->A = _fmpz_vec_init(h->max + 1);
    fmpz_one(h->A);
    /* One prefix sum per suffix length gives all R[q] simultaneously. */
    for (slong j = 1; j < n; ++j)
        for (slong t = 1; t <= h->cap[j]; ++t)
            fmpz_add(h->A + t, h->A + t, h->A + t - 1);
    for (slong t = 1; t <= h->max; ++t)
        fmpz_add(h->A + t, h->A + t, h->A + t - 1);
    h->block = flint_malloc((h->N / 2 + 2) * sizeof(slong));
    for (slong q = 0; q <= ((h->N + 1) / 2); ++q)
        h->block[q] = -1;
    h->bounds = flint_calloc(n, sizeof(*h->bounds));
    return h;
}

static void system_free(system_t_data *h)
{
    if (h->sub)
        system_free(h->sub);
    for (slong r = 1; r < h->n; ++r)
        flint_free(h->bounds[r]);
    flint_free(h->bounds);
    flint_free(h->block);
    flint_free(h->cap);
    _fmpz_vec_clear(h->A, h->max + 1);
    _fmpz_vec_clear(h->C, h->N + 1);
    flint_free(h);
}

static void binom(fmpz_t out, slong n, slong k)
{
    if (k < 0 || k > n)
        fmpz_zero(out);
    else
        fmpz_bin_uiui(out, n, k);
}

static void uniform_drop(fmpz_t drop, slong n, slong r, slong d, slong ell)
{
    if (ell >= d)
        return;
    slong t = r - 1, Np = n - 1, Q = (ell + 1) / t, R = (ell + 1) % t;
    fmpz_t a, b;
    fmpz_init(a);
    fmpz_init(b);
    if (Q >= 1 + (d - 1) / Np) {
        binom(a, Np, t);
        binom(b, n - 2, r - 2);
        fmpz_mul_ui(drop, a, ell + 1);
        fmpz_submul_ui(drop, b, d);
    } else if (d / Np <= Q) {
        slong d0 = d - Q * Np;
        for (slong j = 0; j < R; ++j) {
            binom(a, d0, j);
            binom(b, Np - d0, t - j);
            fmpz_mul(a, a, b);
            fmpz_addmul_ui(drop, a, R - j);
        }
    }
    fmpz_clear(a);
    fmpz_clear(b);
}

/* Convert bounds only after exact subtraction. No layer can have more
   survivors than the number of physically allocated prefix records. */
static slong bound_word(const fmpz_t x) { return fmpz_fits_si(x) ? fmpz_get_si(x) : WORD_MAX; }

static slong full_bound(system_t_data *h, slong r, slong s)
{
    slong n = h->n, *d = h->d;
    if (!h->bounds[r]) {
        h->bounds[r] = flint_malloc((h->cap[r] + 1) * sizeof(slong));
        for (slong j = 0; j <= h->cap[r]; ++j)
            h->bounds[r][j] = -1;
    }
    slong *result = h->bounds[r] + s;
    if (*result >= 0)
        return *result;
    if (r == 1) {
        *result = 0;
        for (slong i = 0; i < n; ++i)
            *result += d[i] - 1 >= s;
        return *result;
    }
    int uniform = d[0] == d[n - 1];
    if (!uniform && r >= 3 && s == h->cap[r])
        return *result = dixon_exterior_rank(d, n, r, s, 1);

    fmpz_t total, term, drop, z2;
    fmpz_init(total);
    fmpz_init(term);
    fmpz_init(drop);
    fmpz_init(z2);
    if (uniform) {
        slong ell = h->cap[r] - s;
        binom(term, ell + r - 1, r - 1);
        binom(total, n - 1, r - 1);
        fmpz_mul(total, total, term);
        if (r >= 3)
            uniform_drop(drop, n, r, d[0], ell);
    } else {
        slong I[r - 1];
        for (slong j = 0; j < r - 1; ++j)
            I[j] = j;
        do {
            slong ell = d[n - 1] - 1 - s;
            for (slong j = 0; j < r - 1; ++j)
                ell += d[I[j]] - 1;
            if (ell < 0)
                continue;
            binom(term, ell + r - 1, r - 1);
            fmpz_add(total, total, term);
            fmpz_add_ui(z2, z2, ell + 1);
        } while (next_subset(I, r - 1, n - 1));
        if (r >= 3 && !fmpz_is_zero(z2)) {
            slong rank = dixon_exterior_rank(d, n, r, s, 0);
            if (rank)
                fmpz_sub_ui(drop, z2, rank);
            if (fmpz_sgn(drop) < 0)
                fmpz_zero(drop);
        }
    }
    fmpz_sub(total, total, drop);
    *result = bound_word(total);
    fmpz_clear(total);
    fmpz_clear(term);
    fmpz_clear(drop);
    fmpz_clear(z2);
    return *result;
}

static slong prefix_bound(system_t_data *h, slong k, const slong *g, slong s)
{
    slong n = h->n, r = n - 1 - k, req[k], acc = 0, cap = 0;
    int unrestricted = 1;
    for (slong j = 0; j < k; ++j) {
        req[j] = (acc += g[k - 1 - j]);
        cap += h->d[k - j] - 1;
        if (acc > cap)
            unrestricted = 0;
    }
    if (unrestricted)
        return full_bound(h, r, s);
    if (h->d[0] == h->d[n - 1])
        return 0;

    /* At each residual degree u, count allowed subsets containing the
       largest active vertex. A subset contributes through ell inclusive. */
    slong levels = h->cap[r] - s + 1, I[r];
    slong *vertex = flint_malloc(checked_mul(levels, sizeof(slong)));
    slong *count = flint_calloc(levels, sizeof(slong));
    for (slong u = 0; u < levels; ++u) {
        vertex[u] = -1;
        count[u] = 0;
    }
    for (slong j = 0; j < r; ++j)
        I[j] = j;
    slong allowed = 0;
    do {
        slong ell = -s;
        for (slong j = 0; j < r; ++j)
            ell += h->d[I[j]] - 1;
        if (ell < 0)
            continue;
        slong a = r - 1, b = 0, sum = 0;
        for (slong i = n - 1; b < k; --i) {
            if (a >= 0 && I[a] == i) {
                --a;
                continue;
            }
            sum += h->d[i] - 1;
            if (sum < req[b])
                break;
            ++b;
        }
        if (b != k)
            continue;
        ++allowed;
        if (r == 1)
            continue;
        for (slong u = 0; u <= ell; ++u) {
            if (vertex[u] < I[r - 1]) {
                vertex[u] = I[r - 1];
                count[u] = 0;
            }
            if (vertex[u] == I[r - 1])
                ++count[u];
        }
    } while (next_subset(I, r, n));
    if (r == 1) {
        flint_free(vertex);
        flint_free(count);
        return allowed;
    }
    fmpz_t total, term;
    fmpz_init(total);
    fmpz_init(term);
    for (slong u = 0; u < levels; ++u) {
        binom(term, u + r - 2, r - 2);
        fmpz_addmul_ui(total, term, count[u]);
    }
    flint_free(vertex);
    flint_free(count);
    slong val = bound_word(total);
    fmpz_clear(total);
    fmpz_clear(term);
    return val;
}

typedef struct {
    slong *g, sum, parent_sum, val, bound;
} prefix_t;
typedef struct {
    prefix_t *items;
    slong len, capacity, k;
    slong *exps;
} layer_t;

/* Generate in lexical order, keeping complete prefix subtrees consecutive. */
static void generate(layer_t *L, system_t_data *h, slong *g, slong j, slong left, slong q)
{
    if (j < L->k) {
        for (slong i = max_s(0, left - h->cap[L->k - j]); i <= left; ++i) {
            g[j] = i;
            generate(L, h, g, j + 1, left - i, q);
        }
        return;
    }
    if (left > h->cap[1])
        return;
    if (L->len == L->capacity) {
        if (L->capacity > WORD_MAX / 2)
            flint_throw(FLINT_OVERFLOW, "Dixon rank predictor layer size overflow");
        L->capacity = L->capacity ? 2 * L->capacity : 32;
        L->items = flint_realloc(L->items, checked_mul(L->capacity, sizeof(prefix_t)));
        L->exps =
            flint_realloc(L->exps, checked_mul(checked_mul(L->capacity, L->k), sizeof(slong)));
    }
    memcpy(L->exps + L->len * L->k, g, L->k * sizeof(slong));
    L->items[L->len++] = (prefix_t){.sum = q - left, .val = 1};
}

static int group_order(const void *a, const void *b)
{
    prefix_t *x = *(prefix_t *const *)a, *y = *(prefix_t *const *)b;
    if (x->sum != y->sum)
        return compare_long(&x->sum, &y->sum);
    if (x->parent_sum != y->parent_sum)
        return compare_long(&x->parent_sum, &y->parent_sum);
    return (x > y) - (x < y);
}

static void window_rank(fmpz_t out, system_t_data *h, slong low, slong high);

static slong block_rank(system_t_data *h, slong q)
{
    if (h->block[q] >= 0)
        return h->block[q];
    slong n = h->n, k0 = n - 2, g[k0], dn = h->d[n - 1], prev = h->d[n - 2];
    layer_t L = {.k = k0};
    generate(&L, h, g, 0, q, q);
    for (slong i = 0; i < L.len; ++i)
        L.items[i].g = L.exps + i * k0;
    prefix_t **order = flint_malloc(L.len * sizeof(*order));
    for (slong k = k0; k > 0; --k) {
        for (slong i = 0; i < L.len; ++i) {
            prefix_t *p = L.items + i;
            p->parent_sum = p->sum - p->g[k - 1];
            p->bound = prefix_bound(h, k, p->g, q - p->sum);
            order[i] = p;
        }
        qsort(order, L.len, sizeof(*order), group_order);
        slong surv = 0, ub = 0;
        for (slong start = 0, end; start < L.len; start = end) {
            slong t = order[start]->sum;
            for (end = start; end < L.len && order[end]->sum == t; ++end) {
                surv += order[end]->val;
                ub = sat_add(ub, order[end]->bound);
            }
            int sublevel = k == k0 && q - t >= prev;
            if (sublevel && n >= 4 && q - prev > (h->N - (dn - 1)) / 2) {
                if (!h->sub)
                    h->sub = system_new(h->d, n - 1);
                fmpz_t v;
                fmpz_init(v);
                window_rank(v, h->sub, max_s(0, q - (dn - 1)), t);
                ub = min_s(ub, bound_word(v));
                fmpz_clear(v);
            }
            slong excess = surv - ub;
            if (excess <= 0)
                continue;
            surv = ub;
            for (slong j = start; j < end && excess; ++j) {
                prefix_t *p = order[j];
                slong cut = min_s(excess, max_s(0, p->val - p->bound));
                p->val -= cut;
                excess -= cut;
            }
            for (slong j = start; j < end && excess; ++j) {
                prefix_t *p = order[sublevel ? end - 1 - (j - start) : j];
                slong cut = min_s(excess, p->val);
                p->val -= cut;
                excess -= cut;
            }
        }
        if (k > 1) {
            slong len = 0;
            for (slong i = 0; i < L.len; ++i) {
                prefix_t p = L.items[i];
                if (len && !memcmp(L.items[len - 1].g, p.g, (k - 1) * sizeof(slong)))
                    L.items[len - 1].val += p.val;
                else {
                    p.sum = p.parent_sum;
                    L.items[len++] = p;
                }
            }
            L.len = len;
        }
    }
    slong rank = 0;
    for (slong i = 0; i < L.len; ++i)
        rank += L.items[i].val;
    flint_free(order);
    flint_free(L.exps);
    flint_free(L.items);
    return h->block[q] = rank;
}

static void deficits(fmpz_t total, system_t_data *h, slong low, slong high, FILE *details)
{
    fmpz_zero(total);
    if (h->n == 2)
        return;
    fmpz_t delta;
    fmpz_init(delta);
    for (slong q = h->N + 1 - h->max; q <= (h->N + 1) / 2; ++q) {
        slong p = h->N + 1 - q;
        int mult = (low <= q && q <= high) + (p != q && low <= p && p <= high);
        if (!mult)
            continue;
        fmpz_sub(delta, h->A + q, h->A + q - 1);
        fmpz_sub_ui(delta, delta, block_rank(h, q));
        if (fmpz_sgn(delta) <= 0)
            continue;
        fmpz_addmul_ui(total, delta, mult);
        if (details) {
            flint_fprintf(details, "  boundary_deficit (%wd,%wd) = %{fmpz}", p, q, delta);
            if (p != q)
                flint_fprintf(details, " (also (%wd,%wd))", q, p);
            fputc('\n', details);
        }
    }
    fmpz_clear(delta);
}

static void midpoint(fmpz_t rho, system_t_data *h)
{
    slong a = h->N / 2, b = (h->N + 1) / 2;
    fmpz_add(rho, h->A + a, h->A + b);
    fmpz_sub(rho, rho, h->C + a);
}

static void window_rank(fmpz_t out, system_t_data *h, slong low, slong high)
{
    high = min_s(high, h->max);
    if (low > high) {
        fmpz_zero(out);
        return;
    }
    if (h->n == 2) {
        fmpz_set_ui(out, high - low + 1);
        return;
    }
    fmpz_t tmp, base, def;
    fmpz_init(tmp);
    fmpz_init(base);
    fmpz_init(def);
    if (low == 0 && high == h->max)
        midpoint(out, h);
    else {
        if (low)
            fmpz_set(base, h->A + low - 1);
        fmpz_sub(out, h->A + high, base);
        for (slong s = low - 1; s < high; ++s) {
            slong j = h->N - s;
            fmpz_set(tmp, h->A + min_s(j, h->max));
            if (j <= h->N)
                fmpz_sub(tmp, tmp, h->C + j);
            if (s >= 0)
                fmpz_add(tmp, tmp, h->A + s);
            fmpz_sub(tmp, tmp, base);
            if (fmpz_cmp(tmp, out) < 0)
                fmpz_set(out, tmp);
        }
    }
    deficits(def, h, low, high, NULL);
    fmpz_sub(out, out, def);
    if (low || high != h->max) {
        deficits(def, h, 0, h->max, NULL);
        if (!fmpz_is_zero(def)) {
            midpoint(tmp, h);
            fmpz_sub(tmp, tmp, def);
            if (fmpz_cmp(tmp, out) < 0)
                fmpz_set(out, tmp);
            if (fmpz_sgn(out) < 0)
                fmpz_zero(out);
        }
    }
    fmpz_clear(tmp);
    fmpz_clear(base);
    fmpz_clear(def);
}

int dixon_generic_rank_from_degrees(fmpz_t rank, const long *degrees, slong n, FILE *details)
{
    if (!degrees || n < 2)
        return 0;
    if (n == 2) {
        if (degrees[0] <= 0 || degrees[1] <= 0)
            return 0;
        fmpz_set_si(rank, max_s(degrees[0], degrees[1]));
        if (details)
            flint_fprintf(details, "n=2 M=%{fmpz}\n  rho = %{fmpz}\n"
                          "  predicted_deficit = 0\n  predicted_rank = %{fmpz}\n",
                          rank, rank, rank);
        return 1;
    }
    slong total = 0;
    for (slong i = 0; i < n; ++i) {
        if (degrees[i] <= 0 || degrees[i] - 1 > WORD_MAX - n - 2 - total)
            return 0;
        total += degrees[i] - 1;
    }
    if ((ulong)total + 1 > SIZE_MAX / sizeof(fmpz))
        return 0;
    slong *sorted = flint_malloc(checked_mul(n, sizeof(*sorted)));
    for (slong i = 0; i < n; ++i)
        sorted[i] = degrees[i];
    qsort(sorted, n, sizeof(*sorted), compare_long);
    system_t_data *system = system_new(sorted, n);
    fmpz_t rho, deficit;
    fmpz_init(rho);
    fmpz_init(deficit);
    midpoint(rho, system);
    if (details)
        flint_fprintf(details, "n=%wd M=%{fmpz}\n", n, system->A + system->max);
    deficits(deficit, system, 0, system->max, details);
    fmpz_sub(rank, rho, deficit);
    if (details)
        flint_fprintf(details, "  rho = %{fmpz}\n  predicted_deficit = %{fmpz}\n"
                      "  predicted_rank = %{fmpz}\n", rho, deficit, rank);
    fmpz_clear(rho);
    fmpz_clear(deficit);
    system_free(system);
    flint_free(sorted);
    return 1;
}

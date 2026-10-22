/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_dixon.h"
#include "unified_mpoly_resultant.h"
#include "unified_mpoly_det.h"
#include "unified_mat_poly.h"
#include <flint/fmpz_mod_mat.h>
#include <flint/nmod_mat.h>

void unified_mpoly_reindex(unified_mpoly_t out, const unified_mpoly_t in, const slong *map)
{
    switch (in->field_id) {
    case FIELD_ID_NMOD:
        nmod_mpoly_compose_nmod_mpoly_gen(GET_NMOD_POLY(out), GET_NMOD_POLY(in), map,
                                          GET_NMOD_CTX(in->ctx_ptr), GET_NMOD_CTX(out->ctx_ptr));
        break;
    case FIELD_ID_FMPZ_MOD:
        fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(GET_BIG_POLY(out), GET_BIG_POLY(in), map,
                                                  GET_BIG_CTX(in->ctx_ptr),
                                                  GET_BIG_CTX(out->ctx_ptr));
        break;
    case FIELD_ID_FQ_ZECH: {
        /* FLINT declares, but does not export, the Zech composition helper. */
        unified_mpoly_t tmp = unified_mpoly_init(out->ctx_ptr);
        ulong *e = flint_calloc(FLINT_MAX(1, in->ctx_ptr->nvars), sizeof(ulong));
        ulong *d = flint_calloc(FLINT_MAX(1, out->ctx_ptr->nvars), sizeof(ulong));
        field_elem_u c, old;
        field_ctx_t *f = in->ctx_ptr->field_ctx;
        const void *fc = field_ctx_data(f);
        field_init_elem(&c, f->field_id, fc);
        field_init_elem(&old, f->field_id, fc);
        for (slong t = 0; t < unified_mpoly_length(in); t++) {
            unified_mpoly_get_term_exp_ui(e, in, t, in->ctx_ptr);
            memset(d, 0, out->ctx_ptr->nvars * sizeof(ulong));
            for (slong j = 0; j < in->ctx_ptr->nvars; j++)
                d[map[j]] += e[j];
            unified_mpoly_get_term_coeff_ui(&c, in, t, in->ctx_ptr);
            unified_mpoly_get_coeff_ui(&old, tmp, d);
            field_add(&c, &c, &old, f->field_id, fc);
            unified_mpoly_set_coeff_ui(tmp, &c, d);
        }
        unified_mpoly_swap(out, tmp);
        unified_mpoly_clear(tmp);
        field_clear_elem(&c, f->field_id, fc);
        field_clear_elem(&old, f->field_id, fc);
        flint_free(e);
        flint_free(d);
        break;
    }
    default:
        fq_nmod_mpoly_compose_fq_nmod_mpoly_gen(GET_FQ_POLY(out), GET_FQ_POLY(in), map,
                                                GET_FQ_CTX(in->ctx_ptr), GET_FQ_CTX(out->ctx_ptr));
        break;
    }
}

static field_elem_u *elements(slong count, field_ctx_t *f)
{
    field_elem_u *a = flint_malloc(FLINT_MAX(1, count) * sizeof(*a));
    for (slong i = 0; i < count; i++)
        field_init_elem(a + i, f->field_id, field_ctx_data(f));
    return a;
}
static void clear_elements(field_elem_u *a, slong count, field_ctx_t *f)
{
    for (slong i = 0; i < count; i++)
        field_clear_elem(a + i, f->field_id, field_ctx_data(f));
    flint_free(a);
}

/* Full pivoting records the actual row/column minor, including rectangular matrices. */
static slong select_minor(field_elem_u *a, slong nr, slong nc, slong *rows, slong *cols,
                          field_ctx_t *f)
{
    const void *fc = field_ctx_data(f);
    field_id_t id = f->field_id;
    field_elem_u *tmp = elements(3, f);
    for (slong i = 0; i < nr; i++)
        rows[i] = i;
    for (slong j = 0; j < nc; j++)
        cols[j] = j;
    slong k;
    for (k = 0; k < FLINT_MIN(nr, nc); k++) {
        slong r = k, c = k;
        for (; r < nr; r++) {
            for (c = k; c < nc; c++)
                if (!field_is_zero(a + r * nc + c, id, fc))
                    break;
            if (c < nc)
                break;
        }
        if (r == nr)
            break;
        slong t = rows[k];
        rows[k] = rows[r];
        rows[r] = t;
        t = cols[k];
        cols[k] = cols[c];
        cols[c] = t;
        /* Ownership swaps, not copies: every element remains owned exactly once. */
        for (slong j = 0; j < nc; j++) {
            field_elem_u v = a[k * nc + j];
            a[k * nc + j] = a[r * nc + j];
            a[r * nc + j] = v;
        }
        for (slong i = 0; i < nr; i++) {
            field_elem_u v = a[i * nc + k];
            a[i * nc + k] = a[i * nc + c];
            a[i * nc + c] = v;
        }
        field_inv(tmp, a + k * nc + k, id, fc);
        for (slong i = k + 1; i < nr; i++) {
            if (field_is_zero(a + i * nc + k, id, fc))
                continue;
            field_mul(tmp + 1, a + i * nc + k, tmp, id, fc);
            for (slong j = k + 1; j < nc; j++) {
                field_mul(tmp + 2, tmp + 1, a + k * nc + j, id, fc);
                field_sub(a + i * nc + j, a + i * nc + j, tmp + 2, id, fc);
            }
            field_set_zero(a + i * nc + k, id, fc);
        }
    }
    clear_elements(tmp, 3, f);
    return k;
}

static slong monomial_index(ulong *basis, slong *count, const ulong *e, slong n)
{
    for (slong i = 0; i < *count; i++)
        if (!memcmp(basis + i * n, e, n * sizeof(ulong)))
            return i;
    memcpy(basis + (*count) * n, e, n * sizeof(ulong));
    return (*count)++;
}

/* Greedy low-degree bases avoid unnecessary parameter factors in the minor. */
static void degree_order(slong *order, slong *degree, slong count)
{
    for (slong i = 0; i < count; i++) {
        order[i] = i;
        slong j = i;
        while (j > 0 && degree[order[j - 1]] > degree[order[j]]) {
            slong t = order[j];
            order[j] = order[j - 1];
            order[j - 1] = t;
            j--;
        }
    }
}

/* The same weighted-degree certificate used by the nmod Schur path:
 * deg M[i,j] <= sigma-rd[i]-cd[j]. Equal-weight pairs are constant;
 * negative bounds force structural zeros. */
static int schur_degrees(unified_mpoly_t **m, slong nr, slong nc, const slong *rd, const slong *cd,
                         slong sigma)
{
    if (!rd || !cd || sigma < 0)
        return 0;
    for (slong i = 0; i < nr; i++)
        for (slong j = 0; j < nc; j++)
            if (!unified_mpoly_is_zero(m[i][j]) &&
                unified_mpoly_total_degree_si(m[i][j]) > sigma - rd[i] - cd[j])
                return 0;
    return 1;
}

/* Collect independent constant pivots in each degree layer, just as the
 * nmod complement-first selector. Their union is block triangular. */
static slong schur_seed(slong *rows, slong *cols, unified_mpoly_t **m, slong nr, slong nc,
                        const slong *rd, const slong *cd, slong sigma)
{
    field_ctx_t *f = m[0][0]->ctx_ptr->field_ctx;
    slong *ri = flint_malloc(nr * sizeof(slong)), *ci = flint_malloc(nc * sizeof(slong));
    slong *rp = flint_malloc(nr * sizeof(slong)), *cp = flint_malloc(nc * sizeof(slong));
    unsigned char *ur = flint_calloc(nr, 1), *uc = flint_calloc(nc, 1);
    field_elem_u *a = elements(nr * nc, f);
    slong count = 0;
    for (slong i = 0; i < nr; i++) {
        slong d = rd[i], r = 0, c = 0;
        int seen = 0;
        for (slong j = 0; j < i; j++)
            if (rd[j] == d) {
                seen = 1;
                break;
            }
        if (seen)
            continue;
        for (slong j = 0; j < nr; j++)
            if (rd[j] == d)
                ri[r++] = j;
        for (slong j = 0; j < nc; j++)
            if (cd[j] == sigma - d)
                ci[c++] = j;
        for (slong u = 0; u < r; u++)
            for (slong v = 0; v < c; v++) {
                field_set_zero(a + u * c + v, f->field_id, field_ctx_data(f));
                if (!unified_mpoly_is_zero(m[ri[u]][ci[v]]))
                    unified_mpoly_get_term_coeff_ui(a + u * c + v, m[ri[u]][ci[v]], 0,
                                                    m[0][0]->ctx_ptr);
            }
        slong rank = select_minor(a, r, c, rp, cp, f);
        for (slong j = 0; j < rank; j++) {
            rows[count] = ri[rp[j]];
            cols[count] = ci[cp[j]];
            ur[rows[count]] = uc[cols[count]] = 1;
            count++;
        }
    }
    slong r = count, c = count;
    for (slong i = 0; i < nr; i++)
        if (!ur[i])
            rows[r++] = i;
    for (slong i = 0; i < nc; i++)
        if (!uc[i])
            cols[c++] = i;
    clear_elements(a, nr * nc, f);
    flint_free(ri);
    flint_free(ci);
    flint_free(rp);
    flint_free(cp);
    flint_free(ur);
    flint_free(uc);
    return count;
}

slong unified_dixon_select_minor(slong *best_r, slong *best_c, unified_mpoly_t **m, slong nr,
                                 slong nc, const slong *rd, const slong *cd, slong sigma)
{
    if (!nr || !nc)
        return 0;
    field_ctx_t *f = m[0][0]->ctx_ptr->field_ctx;
    slong np = m[0][0]->ctx_ptr->nvars, rank = 0, seed = 0;
    slong cap = FLINT_MAX(nr, nc);
    slong *ro = flint_malloc(cap * sizeof(slong)), *co = flint_malloc(cap * sizeof(slong));
    slong *rows = flint_malloc(cap * sizeof(slong)), *cols = flint_malloc(cap * sizeof(slong));
    slong *degree = flint_malloc(cap * sizeof(slong));
    slong *cache = flint_malloc(nr * nc * sizeof(slong));
    ulong *points = flint_malloc(FLINT_MAX(1, np) * sizeof(ulong));
    field_elem_u *values = elements(nr * nc, f), *work = elements(nr * nc, f);
    for (slong i = 0; i < nr; i++)
        for (slong j = 0; j < nc; j++)
            cache[i * nc + j] = unified_mpoly_total_degree_si(m[i][j]);
    if (np == 1 && g_dixon_mq_step4_schur && g_dixon_schur_seed &&
        schur_degrees(m, nr, nc, rd, cd, sigma))
        seed = schur_seed(ro, co, m, nr, nc, rd, cd, sigma);
    int best_verified = 0;
    for (slong trial = 0; trial < 3; trial++) {
        for (slong j = 0; j < np; j++)
            points[j] = 13 + 7 * (trial + 1) * (j + 1);
        for (slong i = 0; i < nr; i++)
            for (slong j = 0; j < nc; j++)
                unified_mpoly_evaluate_all_ui(values + i * nc + j, m[i][j], points);
        if (!seed) {
            for (slong i = 0; i < nr; i++) {
                degree[i] = -1;
                for (slong j = 0; j < nc; j++)
                    degree[i] = FLINT_MAX(degree[i], cache[i * nc + j]);
            }
            degree_order(ro, degree, nr);
            for (slong j = 0; j < nc; j++)
                co[j] = j;
        }
        for (slong i = 0; i < nr; i++)
            for (slong j = 0; j < nc; j++)
                field_set_elem(work + i * nc + j, values + ro[i] * nc + co[j], f->field_id,
                               field_ctx_data(f));
        slong r = select_minor(work, nr, nc, rows, cols, f);
        for (slong i = 0; i < r; i++) {
            rows[i] = ro[rows[i]];
            cols[i] = co[cols[i]];
        }
        /* Without a certified seed, alternate greedy low-degree column and
         * row bases. Degrees are measured only against the opposite basis. */
        for (slong iteration = 0; !seed && r && iteration < 3; iteration++) {
            for (slong side = 0; side < 2; side++) {
                slong count = side ? nr : nc;
                slong *selected = side ? rows : cols, *opposite = side ? cols : rows;
                for (slong i = 0; i < count; i++) {
                    degree[i] = -1;
                    for (slong j = 0; j < r; j++)
                        degree[i] = FLINT_MAX(degree[i], side ? cache[i * nc + opposite[j]]
                                                              : cache[opposite[j] * nc + i]);
                }
                slong *order = flint_malloc(count * sizeof(slong));
                degree_order(order, degree, count);
                for (slong i = 0; i < count; i++)
                    for (slong j = 0; j < r; j++)
                        field_set_elem(work + i * r + j,
                                       values + (side ? order[i] * nc + opposite[j]
                                                      : opposite[j] * nc + order[i]),
                                       f->field_id, field_ctx_data(f));
                select_minor(work, count, r, ro, co, f);
                for (slong i = 0; i < r; i++)
                    selected[i] = order[ro[i]];
                flint_free(order);
            }
        }
        int verified = 1;
        if (g_dixon_step3_second_verification && r) {
            for (slong j = 0; j < np; j++)
                points[j] = 13 + 7 * (trial + 2) * (j + 1);
            for (slong i = 0; i < r; i++)
                for (slong j = 0; j < r; j++)
                    unified_mpoly_evaluate_all_ui(work + i * r + j, m[rows[i]][cols[j]], points);
            /* Keep seed ordering intact across retries. */
            slong *rp = flint_malloc(r * sizeof(slong)), *cp = flint_malloc(r * sizeof(slong));
            verified = select_minor(work, r, r, rp, cp, f) == r;
            flint_free(rp);
            flint_free(cp);
        }
        if (r > rank || (r == rank && verified > best_verified)) {
            rank = r;
            best_verified = verified;
            memcpy(best_r, rows, r * sizeof(slong));
            memcpy(best_c, cols, r * sizeof(slong));
        }
        if (best_verified && rank == FLINT_MIN(nr, nc))
            break;
    }
    if (seed && g_dixon_verbose_level >= 2)
        dixon_info_log("  Schur-first selection: seed=%ld, selected=%ld, retained<=%ld\n", seed,
                       rank, rank - seed);
    if (!seed)
        dixon_info_log("  Degree-aware row/column basis selection\n");
    if (rank && !best_verified)
        dixon_info_log("  Warning: selected minor failed second specialization check\n");
    clear_elements(values, nr * nc, f);
    clear_elements(work, nr * nc, f);
    flint_free(ro);
    flint_free(co);
    flint_free(rows);
    flint_free(cols);
    flint_free(degree);
    flint_free(cache);
    flint_free(points);
    return rank;
}

/* Exact adaptive Schur elimination. Only certified constant degree-layer
 * pivots are removed, so updates stay polynomial and respect the same bounds.
 * The leading h x h block and factor satisfy det(input) = factor * det(core).
 * Swaps include inactive cells so the caller can clear the original allocation. */
slong unified_dixon_schur(unified_mpoly_t **m, slong n, unified_mpoly_t factor, slong *rd,
                          slong *cd, slong sigma)
{
    unified_mpoly_one(factor);
    if (!schur_degrees(m, n, n, rd, cd, sigma))
        return n;
    unified_mpoly_ctx_struct *ctx = factor->ctx_ptr;
    field_ctx_t *f = ctx->field_ctx;
    field_elem_u *v = elements(2, f);
    unified_mpoly_t inv = unified_mpoly_init(ctx), tmp = unified_mpoly_init(ctx);
    ulong *zero = flint_calloc(FLINT_MAX(1, ctx->nvars), sizeof(ulong));
    slong h = n;
    while (h) {
        slong r = -1, c = -1;
        /* Bottom degree group first, matching the existing Schur kernels. */
        for (slong i = 0; i < h; i++)
            for (slong j = 0; j < h; j++)
                if (rd[i] + cd[j] == sigma && !unified_mpoly_is_zero(m[i][j]) &&
                    (r < 0 || rd[i] > rd[r])) {
                    r = i;
                    c = j;
                }
        if (r < 0)
            break;
        --h;
        if (r != h) {
            for (slong j = 0; j < n; j++)
                unified_mpoly_swap(m[r][j], m[h][j]);
            slong d = rd[r];
            rd[r] = rd[h];
            rd[h] = d;
            unified_mpoly_neg(factor, factor);
        }
        if (c != h) {
            for (slong i = 0; i < n; i++)
                unified_mpoly_swap(m[i][c], m[i][h]);
            slong d = cd[c];
            cd[c] = cd[h];
            cd[h] = d;
            unified_mpoly_neg(factor, factor);
        }
        unified_mpoly_mul(factor, factor, m[h][h]);
        unified_mpoly_get_term_coeff_ui(v, m[h][h], 0, ctx);
        field_inv(v + 1, v, f->field_id, field_ctx_data(f));
        unified_mpoly_set_coeff_ui(inv, v + 1, zero);
        for (slong j = 0; j < h; j++)
            unified_mpoly_mul(m[h][j], m[h][j], inv);
        for (slong i = 0; i < h; i++)
            if (!unified_mpoly_is_zero(m[i][h]))
                for (slong j = 0; j < h; j++)
                    if (!unified_mpoly_is_zero(m[h][j])) {
                        unified_mpoly_mul(tmp, m[i][h], m[h][j]);
                        unified_mpoly_sub(m[i][j], m[i][j], tmp);
                    }
    }
    flint_free(zero);
    clear_elements(v, 2, f);
    unified_mpoly_clear(inv);
    unified_mpoly_clear(tmp);
    return h;
}

/* Use the same HNF implementation as the extension-field path. */
static int polynomial_minor_det(unified_mpoly_t out, unified_mpoly_t **m, slong n)
{
    field_ctx_t *f = out->ctx_ptr->field_ctx;
    unified_poly_mat_t polys;
    unified_poly_mat_init(polys, n, n, f);
    field_elem_u c;
    field_init_elem(&c, f->field_id, field_ctx_data(f));
    for (slong i = 0; i < n; i++)
        for (slong j = 0; j < n; j++) {
            unified_poly_struct *p = unified_poly_mat_entry(polys, i, j);
            for (slong t = 0; t < unified_mpoly_length(m[i][j]); t++) {
                ulong e;
                unified_mpoly_get_term_exp_ui(&e, m[i][j], t, out->ctx_ptr);
                if (e > WORD_MAX) {
                    field_clear_elem(&c, f->field_id, field_ctx_data(f));
                    unified_poly_mat_clear(polys);
                    return 0;
                }
                unified_mpoly_get_term_coeff_ui(&c, m[i][j], t, out->ctx_ptr);
                unified_poly_set_coeff(p, (slong)e, &c);
            }
        }
    unified_mat_poly_t matrix;
    unified_mat_poly_init(matrix, n, n, f);
    unified_mat_poly_from_unified(matrix, polys);
    unified_poly_mat_clear(polys);
    unified_poly_t det;
    unified_poly_init(det, f);
    int ok = unified_mat_poly_det_hnf(det, matrix);
    if (ok) {
        unified_mpoly_zero(out);
        for (slong i = 0; i < det->length; i++) {
            ulong e = i;
            unified_poly_get_coeff(&c, det, i);
            unified_mpoly_set_coeff_ui(out, &c, &e);
        }
    }
    field_clear_elem(&c, f->field_id, field_ctx_data(f));
    unified_poly_clear(det);
    unified_mat_poly_clear(matrix);
    return ok;
}

/* For a univariate parameter, interpolate using a rigorous row/column degree bound.
   Singular evaluation points contribute determinant zero and are never skipped. */
static int interpolate_minor(unified_mpoly_t out, unified_mpoly_t **m, slong n)
{
    field_ctx_t *f = out->ctx_ptr->field_ctx;
    const void *fc = field_ctx_data(f);
    field_id_t id = f->field_id;
    if (id != FIELD_ID_FMPZ_MOD && id != FIELD_ID_NMOD)
        return 0;
    slong row_bound = 0, col_bound = 0;
    for (slong i = 0; i < n; i++) {
        slong rd = 0, cd = 0;
        for (slong j = 0; j < n; j++) {
            rd = FLINT_MAX(rd, unified_mpoly_total_degree_si(m[i][j]));
            cd = FLINT_MAX(cd, unified_mpoly_total_degree_si(m[j][i]));
        }
        if (rd > WORD_MAX - row_bound || cd > WORD_MAX - col_bound)
            return 0;
        row_bound += rd;
        col_bound += cd;
    }
    slong bound = FLINT_MIN(row_bound, col_bound);
    if (bound == WORD_MAX || (id == FIELD_ID_NMOD && (ulong)bound >= f->ctx.nmod_ctx.n))
        return 0;
    field_elem_u *v = elements(5, f);
    unified_mpoly_t basis = unified_mpoly_init(out->ctx_ptr);
    unified_mpoly_t term = unified_mpoly_init(out->ctx_ptr);
    unified_mpoly_t factor = unified_mpoly_init(out->ctx_ptr);
    unified_mpoly_t constant = unified_mpoly_init(out->ctx_ptr);
    unified_mpoly_one(basis);
    unified_mpoly_zero(out);
    fmpz_mod_mat_t big;
    nmod_mat_t small;
    if (id == FIELD_ID_FMPZ_MOD)
        fmpz_mod_mat_init(big, n, n, f->ctx.big_ctx);
    else
        nmod_mat_init(small, n, n, f->ctx.nmod_ctx.n);
    fmpz_t x;
    fmpz_init(x);
    for (slong k = 0; k <= bound; k++) {
        ulong point = (ulong)k;
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++) {
                unified_mpoly_evaluate_all_ui(v, m[i][j], &point);
                if (id == FIELD_ID_FMPZ_MOD)
                    fmpz_set(fmpz_mod_mat_entry(big, i, j), &v->big);
                else
                    nmod_mat_entry(small, i, j) = v->nmod;
            }
        if (id == FIELD_ID_FMPZ_MOD)
            fmpz_mod_mat_det(&v->big, big, f->ctx.big_ctx);
        else
            v->nmod = nmod_mat_det(small);
        unified_mpoly_evaluate_all_ui(v + 1, out, &point);
        unified_mpoly_evaluate_all_ui(v + 2, basis, &point);
        field_sub(v, v, v + 1, id, fc);
        field_inv(v + 3, v + 2, id, fc);
        field_mul(v, v, v + 3, id, fc);
        ulong zero = 0;
        unified_mpoly_zero(constant);
        unified_mpoly_set_coeff_ui(constant, v, &zero);
        unified_mpoly_mul(term, basis, constant);
        unified_mpoly_add(out, out, term);
        unified_mpoly_gen(factor, 0, out->ctx_ptr);
        fmpz_set_si(x, k);
        unified_mpoly_set_fmpz(constant, x);
        unified_mpoly_sub(factor, factor, constant);
        unified_mpoly_mul(basis, basis, factor);
    }
    fmpz_clear(x);
    if (id == FIELD_ID_FMPZ_MOD)
        fmpz_mod_mat_clear(big, f->ctx.big_ctx);
    else
        nmod_mat_clear(small);
    unified_mpoly_clear(basis);
    unified_mpoly_clear(term);
    unified_mpoly_clear(factor);
    unified_mpoly_clear(constant);
    clear_elements(v, 5, f);
    return 1;
}

int unified_dixon_resultant(unified_mpoly_t out, unified_mpoly_t *polys, slong n)
{
    field_ctx_t *f = out->ctx_ptr->field_ctx;
    if (n < 0 || (f->field_id != FIELD_ID_NMOD && f->field_id != FIELD_ID_FMPZ_MOD))
        return 0;
    dixon_info_log("\nStep 1: Build Dixon polynomial\n");
    double step_start = get_wall_time();
    slong np = out->ctx_ptr->nvars, dim = 2 * n + np;
    unified_mpoly_ctx_t ring = unified_mpoly_ctx_init(dim, ORD_LEX, f);
    unified_mpoly_t **cancel = unified_mpoly_mat_init(n + 1, n + 1, ring);
    slong *map = flint_malloc(FLINT_MAX(1, n + np) * sizeof(slong));
    for (slong i = 0; i <= n; i++) {
        for (slong j = 0; j < n; j++)
            map[j] = j < i ? n + j : j;
        for (slong j = 0; j < np; j++)
            map[n + j] = 2 * n + j;
        for (slong j = 0; j <= n; j++)
            unified_mpoly_reindex(cancel[i][j], polys[j], map);
    }
    flint_free(map);
    unified_mpoly_t divisor = unified_mpoly_init(ring), tmp = unified_mpoly_init(ring);
    int ok = 1;
    for (slong i = n; i > 0 && ok; i--) {
        unified_mpoly_gen(divisor, i - 1, ring);
        unified_mpoly_gen(tmp, n + i - 1, ring);
        unified_mpoly_sub(divisor, divisor, tmp);
        for (slong j = 0; j <= n; j++) {
            unified_mpoly_sub(cancel[i][j], cancel[i][j], cancel[i - 1][j]);
            if (!unified_mpoly_divides(cancel[i][j], cancel[i][j], divisor)) {
                ok = 0;
                break;
            }
        }
    }
    unified_mpoly_t dixon = unified_mpoly_init(ring);
    if (ok) {
        int method = dixon_global_method_step1;
        if (method != DET_METHOD_KRONECKER_NMOD && method != DET_METHOD_BALANCED_SPLIT)
            method = DET_METHOD_RECURSIVE;
        dixon_info_log("  Determinant method: %s\n", dixon_det_method_name(method));
        compute_unified_mpoly_det_with_method(dixon, cancel, n + 1, ring, 0, method);
    }
    unified_mpoly_mat_clear(cancel, n + 1, n + 1);
    unified_mpoly_clear(divisor);
    unified_mpoly_clear(tmp);
    slong len = unified_mpoly_length(dixon), nr = 0, nc = 0;
    dixon_info_log("  Dixon polynomial: %ld terms\n", len);
    dixon_maybe_print_step_time("Step 1", get_wall_time() - step_start);

    dixon_info_log("\nStep 2: Construct Dixon matrix\n");
    step_start = get_wall_time();
    ulong *xb = flint_calloc(FLINT_MAX(1, len * n), sizeof(ulong));
    ulong *yb = flint_calloc(FLINT_MAX(1, len * n), sizeof(ulong));
    ulong *exp = flint_calloc(FLINT_MAX(1, dim), sizeof(ulong));
    slong *ri = flint_malloc(FLINT_MAX(1, len) * sizeof(slong));
    slong *ci = flint_malloc(FLINT_MAX(1, len) * sizeof(slong));
    for (slong t = 0; t < len; t++) {
        unified_mpoly_get_term_exp_ui(exp, dixon, t, ring);
        ri[t] = monomial_index(xb, &nr, exp, n);
        ci[t] = monomial_index(yb, &nc, exp + n, n);
    }
    unified_mpoly_t **matrix = unified_mpoly_mat_init(nr, nc, out->ctx_ptr);
    field_elem_u *v = elements(2, f);
    for (slong t = 0; t < len; t++) {
        unified_mpoly_get_term_exp_ui(exp, dixon, t, ring);
        unified_mpoly_get_term_coeff_ui(v, dixon, t, ring);
        unified_mpoly_t entry = matrix[ri[t]][ci[t]];
        unified_mpoly_get_coeff_ui(v + 1, entry, exp + 2 * n);
        field_add(v, v, v + 1, f->field_id, field_ctx_data(f));
        unified_mpoly_set_coeff_ui(entry, v, exp + 2 * n);
    }
    clear_elements(v, 2, f);
    unified_mpoly_clear(dixon);
    unified_mpoly_ctx_clear(ring);
    flint_free(exp);
    flint_free(ri);
    flint_free(ci);
    dixon_info_log("  Dixon matrix size: %ld x %ld\n", nr, nc);
    dixon_maybe_print_step_time("Step 2", get_wall_time() - step_start);

    dixon_info_log("\nStep 3: Extract maximal-rank submatrix\n");
    step_start = get_wall_time();
    slong *best_r = flint_malloc(FLINT_MAX(1, nr) * sizeof(slong));
    slong *best_c = flint_malloc(FLINT_MAX(1, nc) * sizeof(slong));
    slong *rd = flint_calloc(FLINT_MAX(1, nr), sizeof(slong));
    slong *cd = flint_calloc(FLINT_MAX(1, nc), sizeof(slong));
    slong sigma = 1;
    for (slong i = 0; i <= n; i++)
        sigma += unified_mpoly_total_degree_si(polys[i]) - 1;
    for (slong i = 0; i < nr; i++)
        for (slong j = 0; j < n; j++)
            rd[i] += xb[i * n + j];
    for (slong i = 0; i < nc; i++)
        for (slong j = 0; j < n; j++)
            cd[i] += yb[i * n + j];
    slong rank = unified_dixon_select_minor(best_r, best_c, matrix, nr, nc, rd, cd, sigma);
    dixon_info_log("  Submatrix size: %ld x %ld\n", rank, rank);
    dixon_maybe_print_step_time("Step 3", get_wall_time() - step_start);

    dixon_info_log("\nStep 4: Compute resultant\n");
    step_start = get_wall_time();
    if (rank && ok) {
        unified_mpoly_t **minor = unified_mpoly_mat_init(rank, rank, out->ctx_ptr);
        for (slong i = 0; i < rank; i++)
            for (slong j = 0; j < rank; j++)
                unified_mpoly_set(minor[i][j], matrix[best_r[i]][best_c[j]]);
        unified_mpoly_t factor = unified_mpoly_init(out->ctx_ptr);
        unified_mpoly_one(factor);
        slong size = rank;
        if (np == 1 && g_dixon_mq_step4_schur) {
            slong *rdeg = flint_malloc(rank * sizeof(slong)),
                  *cdeg = flint_malloc(rank * sizeof(slong));
            for (slong i = 0; i < rank; i++) {
                rdeg[i] = rd[best_r[i]];
                cdeg[i] = cd[best_c[i]];
            }
            size = unified_dixon_schur(minor, rank, factor, rdeg, cdeg, sigma);
            flint_free(rdeg);
            flint_free(cdeg);
            if (size < rank)
                dixon_info_log("  Step 4 adaptive Schur: %ld -> %ld\n", rank, size);
        }
        int method = dixon_global_method_step4;
        if (!size) {
            unified_mpoly_one(out);
        } else if (method == DET_METHOD_RECURSIVE || method == DET_METHOD_BALANCED_SPLIT) {
            dixon_info_log("  Determinant method: %s\n", dixon_det_method_name(method));
            compute_unified_mpoly_det_with_method(out, minor, size, out->ctx_ptr, 0, method);
        } else {
            int done = 0;
            if (np == 1 && (method < 0 || method == DET_METHOD_KRONECKER)) {
                dixon_info_log("  Determinant method: %s\n",
                               dixon_det_method_name(DET_METHOD_KRONECKER));
                done = polynomial_minor_det(out, minor, size);
                if (!done)
                    dixon_info_log("  HNF unavailable for this matrix; trying interpolation.\n");
            }
            if (np == 1 && !done && method != DET_METHOD_KRONECKER_NMOD) {
                dixon_info_log("  Determinant method: %s\n",
                               dixon_det_method_name(DET_METHOD_INTERPOLATION));
                done = interpolate_minor(out, minor, size);
            }
            if (!done) {
                dixon_info_log("  Determinant method: %s\n",
                               dixon_det_method_name(DET_METHOD_KRONECKER_NMOD));
                ok = compute_unified_mpoly_det_bareiss(out, minor, size, out->ctx_ptr);
            }
        }
        if (ok)
            unified_mpoly_mul(out, out, factor);
        unified_mpoly_clear(factor);
        unified_mpoly_mat_clear(minor, rank, rank);
    } else {
        unified_mpoly_zero(out);
        dixon_info_log("  %s\n", ok ? "Empty coefficient matrix, resultant is 0"
                                    : "Dixon polynomial construction failed");
    }
    if (ok)
        dixon_info_log("  Resultant summary: %ld terms, total degree %ld\n",
                       unified_mpoly_length(out), unified_mpoly_total_degree_si(out));
    dixon_maybe_print_step_time("Step 4", get_wall_time() - step_start);
    unified_mpoly_mat_clear(matrix, nr, nc);
    flint_free(best_r);
    flint_free(best_c);
    flint_free(rd);
    flint_free(cd);
    flint_free(xb);
    flint_free(yb);
    return ok;
}

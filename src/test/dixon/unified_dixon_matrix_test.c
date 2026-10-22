/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_dixon.h"
#include "unified_mpoly_det.h"
#include <assert.h>

static void parse(unified_mpoly_t p, const char *s)
{
    const char *names[] = {"t"};
    int status =
        p->field_id == FIELD_ID_NMOD
            ? nmod_mpoly_set_str_pretty(GET_NMOD_POLY(p), s, names, GET_NMOD_CTX(p->ctx_ptr))
            : fmpz_mod_mpoly_set_str_pretty(GET_BIG_POLY(p), s, names, GET_BIG_CTX(p->ctx_ptr));
    assert(!status);
}

static void compression(unified_mpoly_ctx_t ctx, const char **entries, const slong *weights,
                        slong sigma, slong expected)
{
    unified_mpoly_t **a = unified_mpoly_mat_init(3, 3, ctx);
    unified_mpoly_t before = unified_mpoly_init(ctx), after = unified_mpoly_init(ctx);
    unified_mpoly_t factor = unified_mpoly_init(ctx), diff = unified_mpoly_init(ctx);
    slong rd[3], cd[3];
    for (slong i = 0; i < 3; i++) {
        rd[i] = weights[i];
        cd[i] = weights[i + 3];
        for (slong j = 0; j < 3; j++)
            parse(a[i][j], entries[3 * i + j]);
    }
    assert(compute_unified_mpoly_det_bareiss(before, a, 3, ctx));
    slong h = unified_dixon_schur(a, 3, factor, rd, cd, sigma);
    assert(h == expected);
    if (h)
        assert(compute_unified_mpoly_det_bareiss(after, a, h, ctx));
    else
        unified_mpoly_one(after);
    unified_mpoly_mul(after, after, factor);
    unified_mpoly_sub(diff, after, before);
    assert(unified_mpoly_is_zero(diff));
    unified_mpoly_clear(before);
    unified_mpoly_clear(after);
    unified_mpoly_clear(factor);
    unified_mpoly_clear(diff);
    unified_mpoly_mat_clear(a, 3, 3);
}

static void selection(unified_mpoly_ctx_t ctx)
{
    unified_mpoly_t **a = unified_mpoly_mat_init(3, 3, ctx);
    const char *entries[] = {"t^8", "0", "t^20", "0", "1", "0", "1", "0", "t^12"};
    for (slong i = 0; i < 3; i++)
        for (slong j = 0; j < 3; j++)
            parse(a[i][j], entries[3 * i + j]);
    slong rows[3], cols[3];
    slong rank = unified_dixon_select_minor(rows, cols, a, 3, 3, NULL, NULL, -1);
    assert(rank == 2);
    for (slong i = 0; i < rank; i++) {
        assert(rows[i] != 0);
        assert(cols[i] != 2);
    }
    /* Two unlucky specializations must not hide a higher rank. */
    const char *bad[] = {"(t-20)*(t-27)", "0", "0", "0", "1", "0", "0", "0", "1"};
    for (slong i = 0; i < 3; i++)
        for (slong j = 0; j < 3; j++)
            parse(a[i][j], bad[3 * i + j]);
    assert(unified_dixon_select_minor(rows, cols, a, 3, 3, NULL, NULL, -1) == 3);
    /* Rank-deficient rectangular matrix: preserve the constant complement
     * while completing the basis, including on repeated specializations. */
    unified_mpoly_mat_clear(a, 3, 3);
    a = unified_mpoly_mat_init(3, 4, ctx);
    const char *seeded[] = {"t^2", "t", "1", "2", "t", "1", "0", "0", "2*t", "2", "0", "0"};
    for (slong i = 0; i < 3; i++)
        for (slong j = 0; j < 4; j++)
            parse(a[i][j], seeded[4 * i + j]);
    slong sr[3], sc[4], rd[] = {0, 1, 1}, cd[] = {0, 1, 2, 2};
    assert(unified_dixon_select_minor(sr, sc, a, 3, 4, rd, cd, 2) == 2);
    unified_mpoly_t **minor = unified_mpoly_mat_init(2, 2, ctx);
    slong mr[2], mc[2];
    for (slong i = 0; i < 2; i++) {
        mr[i] = rd[sr[i]];
        mc[i] = cd[sc[i]];
        for (slong j = 0; j < 2; j++)
            unified_mpoly_set(minor[i][j], a[sr[i]][sc[j]]);
    }
    unified_mpoly_t factor = unified_mpoly_init(ctx);
    assert(unified_dixon_schur(minor, 2, factor, mr, mc, 2) == 0);
    assert(!unified_mpoly_is_zero(factor));
    unified_mpoly_clear(factor);
    unified_mpoly_mat_clear(minor, 2, 2);
    unified_mpoly_mat_clear(a, 3, 4);
}

int main(void)
{
    const char *primes[] = {"257", "18446744073709551629",
                            "170141183460469231731687303715884105727"};
    for (int i = 0; i < 3; i++) {
        fmpz_t p;
        fmpz_init(p);
        fmpz_set_str(p, primes[i], 10);
        field_ctx_t f;
        field_ctx_init_prime(&f, p);
        unified_mpoly_ctx_t ctx = unified_mpoly_ctx_init(1, ORD_LEX, &f);
        const slong w[] = {0, 1, 1, 0, 1, 1};
        const char *partial[] = {"t^2+1", "t", "2*t+1", "t+2", "0", "3", "2*t", "5", "0"};
        compression(ctx, partial, w, 2, 1); /* off-diagonal pivot, nonunit factors */
        const char *singular[] = {"t^2+1", "t", "2*t", "t", "1", "2", "2*t", "2", "4"};
        compression(ctx, singular, w, 2, 2); /* retain singular complement directions */
        const char *invalid[] = {"t^3", "t", "1", "t", "1", "0", "1", "0", "1"};
        compression(ctx, invalid, w, 2, 3); /* degree certificate fails, unchanged */
        const slong triangular[] = {0, 1, 2, 0, 1, 2};
        const char *all[] = {"t^2", "t", "7", "t+1", "3", "0", "5", "0", "0"};
        compression(ctx, all, triangular, 2, 0); /* nonconstant triangular complement */
        selection(ctx);
        unified_mpoly_ctx_clear(ctx);
        field_ctx_clear(&f);
        fmpz_clear(p);
    }
    puts("Unified Dixon selection and exact Schur determinant checks passed.");
    return 0;
}

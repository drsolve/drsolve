/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Include the implementation to exercise the private repair routine without
   adding a test-only symbol to the library ABI. */
#include "../dixon/dixon_flint.c"
#include <assert.h>

static int lu_backend;
static slong reserve_sigma = -1;

static void check_repair(const nmod_mat_t values, slong target,
                         const fq_nmod_ctx_t ctx, int expected, int degree_case)
{
    slong nr = values->r, nc = values->c;
    fq_mvpoly_t ***matrix = flint_malloc((size_t) nr * sizeof(*matrix));
    fq_nmod_t coeff, params[1];
    fq_nmod_init(coeff, ctx);
    fq_nmod_init(params[0], ctx);
    fq_nmod_one(params[0], ctx);
    slong zero = 0;
    for (slong i = 0; i < nr; i++) {
        matrix[i] = flint_calloc((size_t) nc, sizeof(**matrix));
        for (slong j = 0; j < nc; j++) {
            if (!nmod_mat_entry(values, i, j)) continue;
            matrix[i][j] = flint_malloc(sizeof(fq_mvpoly_t));
            fq_mvpoly_init(matrix[i][j], 1, 1, ctx);
            fq_nmod_set_ui(coeff, nmod_mat_entry(values, i, j), ctx);
            slong degree = ((degree_case == 1 && i == 1) ||
                            (degree_case == 2 && j == 1)) ? 8 : 0;
            fq_mvpoly_add_term(matrix[i][j], &zero, &degree, coeff);
        }
    }
    slong *rows = flint_malloc((size_t) target * sizeof(slong));
    slong *cols = flint_malloc((size_t) target * sizeof(slong));
    slong *perm = flint_malloc((size_t) target * sizeof(slong));
    fq_index_degree_pair *ro = flint_malloc((size_t) nr * sizeof(*ro));
    fq_index_degree_pair *co = flint_malloc((size_t) nc * sizeof(*co));
    for (slong i = 0; i < nr; i++) { ro[i].index = i; ro[i].degree = i; }
    for (slong i = 0; i < nc; i++) { co[i].index = i; co[i].degree = i; }
    nmod_mat_t lu, result;
    nmod_mat_init(lu, target, target, values->mod.n);
    nmod_mat_init(result, target, target, values->mod.n);
    for (slong i = 0; i < target; i++) {
        rows[i] = cols[i] = i;
        for (slong j = 0; j < target; j++)
            nmod_mat_entry(lu, i, j) = nmod_mat_entry(values, i, j);
    }
    slong rank = lu_backend == 1 ? nmod_mat_lu_classical(perm, lu, 0)
               : lu_backend == 2 ? nmod_mat_lu_recursive(perm, lu, 0)
               : nmod_mat_lu(perm, lu, 0);
    assert(rank < target);
    int ok = dixon_repair_predicted_minor(matrix, nr, nc, rows, cols, target,
                                          lu, perm, rank, ro, co, reserve_sigma, params, ctx);
    assert(ok == expected);
    if (ok) {
        for (slong i = 0; i < target; i++) {
            assert(rows[i] >= 0 && rows[i] < nr && cols[i] >= 0 && cols[i] < nc);
            for (slong j = 0; j < i; j++) {
                assert(rows[i] != rows[j]);
                assert(cols[i] != cols[j]);
            }
            for (slong j = 0; j < target; j++)
                nmod_mat_entry(result, i, j) = nmod_mat_entry(values, rows[i], cols[j]);
        }
        assert(nmod_mat_rank(result) == target);
        if (degree_case) {
            assert(degree_case == 1 ? rows[1] == 2 : cols[1] == 2);
            assert(dixon_minor_degree_bound(matrix, rows, cols, target) == 0);
        }
    } else {
        for (slong i = 0; i < target; i++) assert(rows[i] == i && cols[i] == i);
    }
    nmod_mat_clear(result);
    nmod_mat_clear(lu);
    flint_free(co); flint_free(ro); flint_free(perm); flint_free(cols); flint_free(rows);
    for (slong i = 0; i < nr; i++) {
        for (slong j = 0; j < nc; j++) if (matrix[i][j]) {
            fq_mvpoly_clear(matrix[i][j]);
            flint_free(matrix[i][j]);
        }
        flint_free(matrix[i]);
    }
    flint_free(matrix);
    fq_nmod_clear(params[0], ctx);
    fq_nmod_clear(coeff, ctx);
}

int main(void)
{
    fq_nmod_ctx_t ctx;
    fmpz_t prime;
    fmpz_init_set_ui(prime, 65537);
    fq_nmod_ctx_init(ctx, prime, 1, "a");
    flint_rand_t state;
    flint_rand_init(state);
    flint_rand_set_seed(state, 31297, 87123);
    /* Random rank-r matrices whose leading r rows span only s directions.
       A zero leading column and row force skipped U columns and row swaps. */
    for (slong trial = 0; trial < 100; trial++) {
        slong target = 4 + trial % 13, s = target - 1 - trial % 3;
        slong nr = target + 10, nc = target + 12;
        nmod_mat_t left, right, values;
        nmod_mat_init(left, nr, target, 65537);
        nmod_mat_init(right, target, nc, 65537);
        nmod_mat_init(values, nr, nc, 65537);
        lu_backend = trial % 3;
        for (slong i = 0; i < nr; i++)
            for (slong j = 0; j < target; j++)
                nmod_mat_entry(left, i, j) = n_randint(state, 65537);
        for (slong i = 0; i < target; i++)
            for (slong j = 0; j < nc; j++)
                nmod_mat_entry(right, i, j) = n_randint(state, 65537);
        /* Plant independent directions so the expected full rank is certain. */
        for (slong i = 0; i < target; i++) {
            slong row = i < s ? i + 1 : target + i - s;
            slong col = i < s ? i + 1 : target + i - s;
            for (slong j = 0; j < target; j++) {
                nmod_mat_entry(left, row, j) = i == j;
                nmod_mat_entry(right, j, col) = i == j;
            }
        }
        for (slong i = 0; i < target; i++) {
            nmod_mat_entry(right, i, 0) = 0;
            for (slong j = s; j < target; j++) nmod_mat_entry(left, i, j) = 0;
        }
        for (slong j = 0; j < target; j++) nmod_mat_entry(left, 0, j) = 0;
        nmod_mat_mul(values, left, right);
        assert(nmod_mat_rank(values) == target);
        check_repair(values, target, ctx, 1, 0);
        nmod_mat_clear(values); nmod_mat_clear(right); nmod_mat_clear(left);
    }
    /* The useful direction arrives after several expansions; one dimension
       reaches its limit first.  Transposition covers the other orientation. */
    nmod_mat_t a, at;
    nmod_mat_init(a, 24, 7, 65537);
    nmod_mat_entry(a, 0, 0) = nmod_mat_entry(a, 1, 1) = 1;
    nmod_mat_entry(a, 20, 6) = 1;
    check_repair(a, 3, ctx, 1, 0);
    nmod_mat_init(at, 7, 24, 65537);
    nmod_mat_transpose(at, a);
    check_repair(at, 3, ctx, 1, 0);
    nmod_mat_entry(a, 20, 6) = 0;
    check_repair(a, 3, ctx, 0, 0); /* Overprediction: preserve fallback inputs. */
    nmod_mat_zero(a);
    check_repair(a, 3, ctx, 0, 0); /* Rank-zero specialization. */
    nmod_mat_clear(at); nmod_mat_clear(a);
    nmod_mat_init(a, 4, 4, 65537);
    nmod_mat_entry(a, 0, 0) = nmod_mat_entry(a, 1, 2) = nmod_mat_entry(a, 2, 2) = 1;
    check_repair(a, 2, ctx, 1, 1); /* Replace the added degree-8 row with degree 0. */
    nmod_mat_init(at, 4, 4, 65537);
    nmod_mat_transpose(at, a);
    check_repair(at, 2, ctx, 1, 2); /* Column degree improvement. */
    nmod_mat_clear(at);
    nmod_mat_clear(a);
    nmod_mat_init(a, 24, 8, 65537);
    nmod_mat_entry(a, 0, 0) = nmod_mat_entry(a, 1, 1) = 1;
    nmod_mat_entry(a, 4, 4) = nmod_mat_entry(a, 20, 7) = 1;
    check_repair(a, 4, ctx, 1, 0); /* Partial completion before further growth. */
    nmod_mat_clear(a);
    /* A direction outside the bounded pool must trigger fallback. */
    nmod_mat_init(a, 40, 40, 65537);
    nmod_mat_entry(a, 0, 0) = nmod_mat_entry(a, 39, 39) = 1;
    check_repair(a, 2, ctx, 0, 0);
    nmod_mat_zero(a);
    nmod_mat_entry(a, 0, 0) = nmod_mat_entry(a, 1, 39) = 1;
    reserve_sigma = 40;
    check_repair(a, 2, ctx, 1, 0); /* Complementary pair beyond the column prefix. */
    nmod_mat_clear(a);
    flint_rand_clear(state);
    fq_nmod_ctx_clear(ctx);
    fmpz_clear(prime);
    flint_cleanup();
    puts("Schur repair: 109 cases passed");
    return 0;
}

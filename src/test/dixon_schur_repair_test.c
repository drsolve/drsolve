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
                            (degree_case == 2 && j == 1) ||
                            (degree_case == 3 && i == 0) ||
                            (degree_case == 4 && j == 0)) ? 8 : 0;
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
            if (degree_case <= 2) assert(degree_case == 1 ? rows[1] == 2 : cols[1] == 2);
            else assert(degree_case == 3 ? rows[0] == 39 : cols[0] == 39);
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

/* Compare each exchange sweep with the original sorted greedy selection on
   exactly the same fixed opposite set and specialization.  This checks the
   objective, index ties, inverse updates in BOTH orientations and rebasing. */
static void check_exchange_axes(flint_rand_t state, const fq_nmod_ctx_t ctx)
{
    slong total_exchanges = 0, total_rebases = 0;
    for (slong trial = 0; trial < 100; trial++) {
        slong r = 1 + trial % 12, nr = 2 * r + 3, nc = 2 * r + 4;
        nmod_mat_t left, right, values;
        nmod_mat_init(left, nr, r, 65537);
        nmod_mat_init(right, r, nc, 65537);
        nmod_mat_init(values, nr, nc, 65537);
        for (slong i = 0; i < nr; i++) for (slong j = 0; j < r; j++)
            nmod_mat_entry(left, i, j) = i < nr - r ? n_randint(state, 65537) : i - (nr - r) == j;
        for (slong i = 0; i < r; i++) for (slong j = 0; j < nc; j++)
            nmod_mat_entry(right, i, j) = j < nc - r ? n_randint(state, 65537) : j - (nc - r) == i;
        nmod_mat_mul(values, left, right);
        fq_mvpoly_t ***matrix = flint_malloc((size_t) nr * sizeof(*matrix));
        fq_nmod_t coeff, params[1];
        fq_nmod_init(coeff, ctx); fq_nmod_init(params[0], ctx); fq_nmod_one(params[0], ctx);
        slong zero = 0;
        for (slong i = 0; i < nr; i++) {
            matrix[i] = flint_calloc((size_t) nc, sizeof(**matrix));
            for (slong j = 0; j < nc; j++) if (nmod_mat_entry(values, i, j)) {
                matrix[i][j] = flint_malloc(sizeof(fq_mvpoly_t));
                fq_mvpoly_init(matrix[i][j], 1, 1, ctx);
                fq_nmod_set_ui(coeff, nmod_mat_entry(values, i, j), ctx);
                slong degree = (trial % 2) ? n_randint(state, 4) + 4 * (i >= nr-r) + 4 * (j >= nc-r) : 0;
                fq_mvpoly_add_term(matrix[i][j], &zero, &degree, coeff);
            }
        }
        dixon_eval_cache_t cache;
        dixon_eval_cache_init(&cache, matrix, nc, params, ctx);
        dixon_exchange_solver_t solver;
        dixon_exchange_solver_init(&solver, r, 65537);
        nmod_mat_one(solver.lu); nmod_mat_one(solver.transposed);
        slong *rows = flint_malloc((size_t) r * sizeof(slong));
        slong *cols = flint_malloc((size_t) r * sizeof(slong));
        mp_limb_t *vector = flint_malloc((size_t) r * sizeof(mp_limb_t));
        mp_limb_t *rhs = flint_malloc((size_t) r * sizeof(mp_limb_t));
        for (slong i = 0; i < r; i++) { rows[i] = nr-r+i; cols[i] = nc-r+i; solver.perm[i] = i; }
        for (int sweep = 0; sweep < 8; sweep++) {
            int columns = sweep % 2;
            slong count = columns ? nc : nr;
            fq_index_degree_pair *order = flint_malloc((size_t) count * sizeof(*order));
            mp_limb_t *vectors = flint_malloc((size_t) count * r * sizeof(*vectors));
            for (slong i = 0; i < count; i++) {
                order[i].index = i;
                order[i].degree = columns
                    ? compute_fq_selected_rows_col_max_total_degree(matrix, rows, r, i, 1)
                    : compute_fq_selected_cols_row_max_total_degree(matrix, i, cols, r, 1);
                for (slong j = 0; j < r; j++)
                    vectors[i*r+j] = columns ? nmod_mat_entry(values, rows[j], i)
                                             : nmod_mat_entry(values, i, cols[j]);
            }
            qsort(order, (size_t) count, sizeof(*order), compare_fq_degrees);
            nmod_row_basis_tracker_t reference;
            nmod_row_basis_tracker_init(&reference, r, r, &values->mod);
            for (slong i = 0; i < count && reference.current_rank < r; i++)
                nmod_try_add_row_to_basis_raw(&reference, vectors, order[i].index);
            assert(reference.current_rank == r);
            slong *stream_rows = flint_malloc((size_t) r * sizeof(slong));
            slong *stream_cols = flint_malloc((size_t) r * sizeof(slong));
            memcpy(stream_rows, rows, (size_t) r * sizeof(slong));
            memcpy(stream_cols, cols, (size_t) r * sizeof(slong));
            dixon_degree_stream_axis(matrix, nr, nc, r, stream_rows, stream_cols,
                                      columns, params, ctx, 1 + trial % 7);
            slong *stream_selected = columns ? stream_cols : stream_rows;
            /* Pivot columns of the transposed panels must preserve greedy ORDER,
               including dependent panels, skipped pivots and row permutations. */
            for (slong i = 0; i < r; i++)
                assert(stream_selected[i] == reference.selected_indices[i]);
            flint_free(stream_cols); flint_free(stream_rows);
            total_exchanges += dixon_exchange_axis(&solver, &cache, nr, nc, rows, cols, columns);
            slong *selected = columns ? cols : rows;
            for (slong i = 0; i < r; i++) {
                int found = 0;
                for (slong j = 0; j < r; j++) if (selected[i] == reference.selected_indices[j]) found = 1;
                assert(found);
            }
            /* Independent residual check after mixed row and column updates. */
            for (int orientation = 0; orientation < 2; orientation++) {
                for (slong i = 0; i < r; i++) rhs[i] = vector[i] = n_randint(state, 65537);
                dixon_exchange_solve(&solver, vector, orientation);
                for (slong i = 0; i < r; i++) {
                    mp_limb_t sum = 0;
                    for (slong j = 0; j < r; j++)
                        sum = nmod_add(sum, nmod_mul(vector[j], orientation
                            ? nmod_mat_entry(values, rows[i], cols[j])
                            : nmod_mat_entry(values, rows[j], cols[i]), values->mod), values->mod);
                    assert(sum == rhs[i]);
                }
            }
            nmod_row_basis_tracker_clear(&reference);
            flint_free(vectors); flint_free(order);
        }
        total_rebases += solver.rebases;
        /* Force hash growth, zero caching and collision handling independently. */
        for (slong i = 0; i < nr; i++) for (slong j = 0; j < nc; j++)
            assert(dixon_eval_cached(&cache, i, j) == nmod_mat_entry(values, i, j));
        assert(cache.count == (size_t) nr * nc);
        size_t cached = cache.count;
        assert(dixon_eval_cached(&cache, 0, 0) == nmod_mat_entry(values, 0, 0));
        assert(cache.count == cached);
        dixon_exchange_solver_clear(&solver); dixon_eval_cache_clear(&cache);
        flint_free(rhs); flint_free(vector); flint_free(cols); flint_free(rows);
        for (slong i = 0; i < nr; i++) {
            for (slong j = 0; j < nc; j++) if (matrix[i][j]) {
                fq_mvpoly_clear(matrix[i][j]); flint_free(matrix[i][j]);
            }
            flint_free(matrix[i]);
        }
        flint_free(matrix);
        fq_nmod_clear(params[0], ctx); fq_nmod_clear(coeff, ctx);
        nmod_mat_clear(values); nmod_mat_clear(right); nmod_mat_clear(left);
    }
    assert(total_exchanges > 100 && total_rebases > 0);
    printf("Degree-aware exchange and blocked streaming: 800 greedy comparisons each passed (%ld exchanges, %ld rebases)\n",
           total_exchanges, total_rebases);
}

/* Exercise automatic large-minor dispatch without building a Dixon polynomial.
   Shared monomials keep this rank-256 regression below a few MiB of input. */
static void check_large_repair_and_cache_limit(const fq_nmod_ctx_t ctx)
{
    slong r = DIXON_DEGREE_BLOCK_THRESHOLD, n = 2 * r + 4, zero = 0;
    fq_mvpoly_t entries[3];
    fq_nmod_t one, params[1];
    fq_nmod_init(one, ctx); fq_nmod_one(one, ctx);
    fq_nmod_init(params[0], ctx); fq_nmod_one(params[0], ctx);
    for (int i = 0; i < 3; i++) {
        slong degree = 8 * i;
        fq_mvpoly_init(&entries[i], 1, 1, ctx);
        fq_mvpoly_add_term(&entries[i], &zero, &degree, one);
    }
    fq_mvpoly_t ***matrix = flint_malloc((size_t) n * sizeof(*matrix));
    fq_index_degree_pair *order = flint_malloc((size_t) n * sizeof(*order));
    slong *rows = flint_malloc((size_t) r * sizeof(slong));
    slong *cols = flint_malloc((size_t) r * sizeof(slong));
    slong *perm = flint_malloc((size_t) r * sizeof(slong));
    for (slong i = 0; i < n; i++) {
        order[i].index = order[i].degree = i;
        matrix[i] = flint_calloc((size_t) n, sizeof(**matrix));
        for (slong j = 0; j < n; j++)
            if (i != r-1 && (i == r ? r-1 : i % r) == j % r)
                matrix[i][j] = &entries[(i < r) + (j < r)];
    }
    nmod_mat_t lu;
    nmod_mat_init(lu, r, r, fq_nmod_ctx_prime(ctx));
    for (slong i = 0; i < r; i++) {
        rows[i] = cols[i] = i;
        if (i < r-1) nmod_mat_entry(lu, i, i) = 1;
    }
    slong rank = nmod_mat_lu(perm, lu, 0);
    assert(rank == r-1);
    assert(dixon_repair_predicted_minor(matrix, n, n, rows, cols, r, lu, perm,
                                       rank, order, order, -1, params, ctx));
    assert(dixon_minor_degree_bound(matrix, rows, cols, r) == 0);
    for (slong i = 0; i < r; i++) for (slong j = 0; j < r; j++)
        nmod_mat_entry(lu, i, j) = matrix[rows[i]][cols[j]] != NULL;
    assert(nmod_mat_rank(lu) == r);
    nmod_mat_clear(lu);
    flint_free(perm); flint_free(cols); flint_free(rows); flint_free(order);
    for (slong i = 0; i < n; i++) flint_free(matrix[i]);
    flint_free(matrix);

    /* Two entirely dependent panels after a row-swapping first pivot. */
    fq_mvpoly_t *sparse[12][4] = {{NULL}}, *transpose[4][12];
    fq_mvpoly_t **small_rows[12], **small_cols[4];
    for (slong i = 0; i < 12; i++) {
        sparse[i][i < 9 ? 3 : i-9] = &entries[0];
        small_rows[i] = sparse[i];
        for (slong j = 0; j < 4; j++) transpose[j][i] = sparse[i][j];
    }
    for (slong j = 0; j < 4; j++) small_cols[j] = transpose[j];
    slong sr[4] = {0,9,10,11}, sc[4] = {0,1,2,3};
    assert(dixon_degree_stream_axis(small_rows,12,4,4,sr,sc,0,params,ctx,3) == 0);
    assert(sr[0] == 0 && sr[1] == 9 && sr[2] == 10 && sr[3] == 11);
    assert(dixon_degree_stream_axis(small_cols,4,12,4,sc,sr,1,params,ctx,3) == 0);
    assert(sr[0] == 0 && sr[1] == 9 && sr[2] == 10 && sr[3] == 11);

    const size_t limit = ((size_t) 16 << 20) / sizeof(dixon_eval_slot_t);
    slong count = (slong) (limit / 2 + 1024);
    fq_mvpoly_t **line = flint_calloc((size_t) count, sizeof(*line));
    line[count-1] = &entries[0];
    dixon_eval_cache_t cache;
    dixon_eval_cache_init(&cache, &line, count, params, ctx);
    for (slong j = 0; j < count; j++)
        assert(dixon_eval_cached(&cache, 0, j) == (j == count-1));
    assert(cache.capacity == limit && cache.count == limit / 2);
    assert(dixon_eval_cached(&cache, 0, count-1) == 1);
    assert(cache.count == limit / 2);
    dixon_eval_cache_clear(&cache); flint_free(line);
    for (int i = 0; i < 3; i++) fq_mvpoly_clear(&entries[i]);
    fq_nmod_clear(params[0], ctx); fq_nmod_clear(one, ctx);
    puts("Large repair dispatch and bounded cache passed");
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
    reserve_sigma = -1;
    nmod_mat_init(a, 40, 4, 65537);
    nmod_mat_entry(a, 0, 0) = nmod_mat_entry(a, 2, 2) = nmod_mat_entry(a, 39, 0) = 1;
    check_repair(a, 2, ctx, 1, 3); /* Replace a degree-8 CORE row outside reserve pool. */
    nmod_mat_init(at, 4, 40, 65537);
    nmod_mat_transpose(at, a);
    check_repair(at, 2, ctx, 1, 4); /* Same regression for a core column. */
    nmod_mat_clear(at); nmod_mat_clear(a);
    check_exchange_axes(state, ctx);
    check_large_repair_and_cache_limit(ctx);
    flint_rand_clear(state);
    fq_nmod_ctx_clear(ctx);
    fmpz_clear(prime);
    flint_cleanup();
    puts("Schur repair: 111 cases passed");
    return 0;
}

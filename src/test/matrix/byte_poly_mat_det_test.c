/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "fq_poly_mat_det.h"
#include <assert.h>

/* Independent permutation expansion (signs are all + in characteristic 2). */
static void expand(fq_nmod_poly_t result, const fq_nmod_poly_mat_t m,
    slong row, ulong used, const fq_nmod_poly_t product, const fq_nmod_ctx_t ctx)
{
    if (row == m->r) { fq_nmod_poly_add(result, result, product, ctx); return; }
    fq_nmod_poly_t term;
    fq_nmod_poly_init(term, ctx);
    for (slong col = 0; col < m->c; ++col) {
        if ((used >> col) & 1) continue;
        fq_nmod_poly_mul(term, product, fq_nmod_poly_mat_entry(m, row, col), ctx);
        expand(result, m, row+1, used | (UWORD(1) << col), term, ctx);
    }
    fq_nmod_poly_clear(term, ctx);
}

static void check(ulong modulus_bits, slong degree)
{
    nmod_poly_t modulus;
    nmod_poly_init(modulus, 2);
    for (slong i = 0; i <= degree; ++i)
        nmod_poly_set_coeff_ui(modulus, i, (modulus_bits >> i) & 1);
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_modulus(ctx, modulus, "a");
    nmod_poly_clear(modulus);
    flint_rand_t state;
    flint_rand_init(state); flint_rand_set_seed(state, modulus_bits, 73);
    for (slong n = 0; n <= 5; ++n) {
        fq_nmod_poly_mat_t mat, copy;
        fq_nmod_poly_mat_init(mat, n, n, ctx);
        fq_nmod_poly_mat_init(copy, n, n, ctx);
        fq_nmod_poly_t expected, actual, one;
        fq_nmod_poly_init(expected, ctx); fq_nmod_poly_init(actual, ctx);
        fq_nmod_poly_init(one, ctx); fq_nmod_poly_one(one, ctx);
        for (int trial = 0; trial < 12; ++trial) {
            for (slong i = 0; i < n; ++i)
                for (slong j = 0; j < n; ++j) {
                    fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(mat, i, j);
                    fq_nmod_poly_randtest(entry, state, 1 + trial % 5, ctx);
                    if (trial == 0 || (trial == 1 && i > j)) fq_nmod_poly_zero(entry, ctx);
                    if (trial == 9 && ((i+j) & 1)) fq_nmod_poly_shift_left(entry, entry, 17, ctx);
                }
            if (trial == 2 && n > 1)
                for (slong j = 0; j < n; ++j)
                    fq_nmod_poly_set(fq_nmod_poly_mat_entry(mat, 1, j),
                                     fq_nmod_poly_mat_entry(mat, 0, j), ctx);
            for (slong i = 0; i < n; ++i)
                for (slong j = 0; j < n; ++j)
                    fq_nmod_poly_set(fq_nmod_poly_mat_entry(copy, i, j),
                                     fq_nmod_poly_mat_entry(mat, i, j), ctx);
            fq_nmod_poly_zero(expected, ctx);
            expand(expected, mat, 0, 0, one, ctx);
            fq_nmod_poly_mat_det_iter(actual, mat, ctx);
            assert(fq_nmod_poly_equal(actual, expected, ctx));
            for (slong i = 0; i < n; ++i)
                for (slong j = 0; j < n; ++j)
                    assert(fq_nmod_poly_equal(fq_nmod_poly_mat_entry(copy, i, j),
                                              fq_nmod_poly_mat_entry(mat, i, j), ctx));
        }
        fq_nmod_poly_clear(expected, ctx); fq_nmod_poly_clear(actual, ctx);
        fq_nmod_poly_clear(one, ctx);
        fq_nmod_poly_mat_clear(mat, ctx); fq_nmod_poly_mat_clear(copy, ctx);
    }
    flint_rand_clear(state); fq_nmod_ctx_clear(ctx);
}

int main(void)
{
    check(0x13, 4); check(0x11d, 8); check(0x11b, 8); check(0x11d, 8);
    puts("Byte polynomial matrix determinants match permutation expansion; inputs unchanged.");
    flint_cleanup_master();
    return 0;
}

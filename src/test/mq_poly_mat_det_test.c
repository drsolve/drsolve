/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_poly_mat_det.h"
#include <assert.h>
#include <stdio.h>

/* Complement [[a(t),b],[c,0]]: nonconstant E with constant nonzero det.
 * This checks exact scaling/sign, including characteristic two, as well as
 * safe rejection without corrupting the caller's fallback input/output. */
int main(void)
{
    ulong primes[] = {2, 3, 257, 65537};
    for (slong z = 0; z < 4; z++)
    {
        ulong p = primes[z], factor = 0;
        slong d[] = {1, 1, 2};
        nmod_poly_mat_t B, saved, core;
        nmod_poly_mat_init(B, 3, 3, p);
        nmod_poly_mat_init(saved, 3, 3, p);
        nmod_poly_mat_init(core, 1, 1, p);
        for (slong i = 0; i < 3; i++)
            for (slong j = 0; j < 3; j++)
            {
                slong degree = 3 - d[i] - d[j];
                for (slong k = 0; k <= degree; k++)
                    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, i, j), k,
                                           (i + 2 * j + k + 1) % p);
            }
        nmod_poly_one(nmod_poly_mat_entry(B, 1, 2));
        nmod_poly_one(nmod_poly_mat_entry(B, 2, 1));
        nmod_poly_mat_set(saved, B);
        nmod_poly_t expected, got;
        nmod_poly_init(expected, p);
        nmod_poly_init(got, p);
        nmod_poly_mat_det(expected, B);
        assert(nmod_poly_mat_mq_schur(core, &factor, B, d, d, 1, 3));
        nmod_poly_mat_det(got, core);
        nmod_poly_scalar_mul_nmod(got, got, factor);
        assert(nmod_poly_equal(expected, got));
        assert(nmod_poly_mat_equal(B, saved));
        assert(factor == p - 1);
        nmod_poly_zero(nmod_poly_mat_entry(B, 2, 1));
        nmod_poly_mat_set(saved, B);
        nmod_poly_one(nmod_poly_mat_entry(core, 0, 0));
        factor = 42;
        assert(!nmod_poly_mat_mq_schur(core, &factor, B, d, d, 1, 3));
        assert(factor == 42 && nmod_poly_is_one(nmod_poly_mat_entry(core, 0, 0)) &&
               nmod_poly_mat_equal(B, saved));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, 2, 2), 1, 1);
        assert(!nmod_poly_mat_mq_schur(core, &factor, B, d, d, 1, 3));
        assert(factor == 42 && nmod_poly_is_one(nmod_poly_mat_entry(core, 0, 0)));
        assert(!nmod_poly_mat_mq_schur(B, &factor, B, d, d, 3, 3));
        /* B can be nonsingular even when this specified E is singular. */
        nmod_poly_mat_zero(B);
        nmod_poly_one(nmod_poly_mat_entry(B, 0, 2));
        nmod_poly_one(nmod_poly_mat_entry(B, 2, 0));
        nmod_poly_set_coeff_ui(nmod_poly_mat_entry(B, 1, 1), 1, 1);
        nmod_poly_mat_det(expected, B);
        assert(nmod_poly_degree(expected) == 1);
        nmod_poly_mat_set(saved, B);
        assert(!nmod_poly_mat_mq_schur(core, &factor, B, d, d, 1, 3));
        assert(nmod_poly_mat_equal(B, saved) && factor == 42);
        nmod_poly_clear(expected);
        nmod_poly_clear(got);
        nmod_poly_mat_clear(B);
        nmod_poly_mat_clear(saved);
        nmod_poly_mat_clear(core);
    }
    puts("MQ Schur: exact determinant/scaling, characteristic two, rejection and input "
         "preservation PASS");
    flint_cleanup_master();
    return 0;
}

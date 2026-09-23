/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_poly_mat_det.h"
#include <flint/nmod.h>
#include <flint/ulong_extras.h>

int nmod_poly_mat_mq_schur(nmod_poly_mat_t core, ulong *factor, const nmod_poly_mat_t B,
                          const slong *rd, const slong *cd, slong h, slong sigma)
{
    slong n = B->r, e = n - h;
    if (!factor || !rd || !cd || n != B->c || h < 0 || h > n || sigma < 0 || core == B ||
        core->r != h || core->c != h || core->modulus != B->modulus)
        return 0;
    slong budget = 0;
    for (slong i = 0; i < n; i++)
    {
        if (rd[i] < 0 || cd[i] < 0 || rd[i] > sigma || cd[i] > sigma)
            return 0;
        if (i >= h)
            budget += sigma - rd[i] - cd[i];
        for (slong j = 0; j < n; j++)
            if (!nmod_poly_is_zero(nmod_poly_mat_entry(B, i, j)) &&
                nmod_poly_degree(nmod_poly_mat_entry(B, i, j)) > sigma - rd[i] - cd[j])
                return 0;
    }
    if (budget != 0)
        return 0;
    nmod_poly_mat_t work;
    nmod_poly_mat_init(work, n, n, B->modulus);
    slong *r = flint_malloc(n * sizeof(slong)), *c = flint_malloc(n * sizeof(slong));
    /* Move E first. Applying the same permutation on both axes has sign +1. */
    for (slong i = 0; i < n; i++)
    {
        slong ii = i < e ? h + i : i - e;
        r[i] = rd[ii];
        c[i] = cd[ii];
        for (slong j = 0; j < n; j++)
        {
            slong jj = j < e ? h + j : j - e;
            nmod_poly_set(nmod_poly_mat_entry(work, i, j), nmod_poly_mat_entry(B, ii, jj));
        }
    }
    nmod_t mod;
    nmod_init(&mod, B->modulus);
    ulong scale = 1;
    int ok = 1;
    for (slong k = 0; k < e; k++)
    {
        slong pr = -1, pc = -1;
        /* Highest-degree boundary row first: it has the narrowest support. */
        for (slong i = k; i < e; i++)
            for (slong j = k; j < e; j++)
                if (r[i] + c[j] == sigma && !nmod_poly_is_zero(nmod_poly_mat_entry(work, i, j)) &&
                    (pr < 0 || r[i] > r[pr]))
                {
                    pr = i;
                    pc = j;
                }
        if (pr < 0)
        {
            ok = 0;
            break;
        }
        if (pr != k)
        {
            for (slong j = 0; j < n; j++)
                nmod_poly_swap(nmod_poly_mat_entry(work, pr, j), nmod_poly_mat_entry(work, k, j));
            slong t = r[pr];
            r[pr] = r[k];
            r[k] = t;
            scale = nmod_neg(scale, mod);
        }
        if (pc != k)
        {
            for (slong i = 0; i < n; i++)
                nmod_poly_swap(nmod_poly_mat_entry(work, i, pc), nmod_poly_mat_entry(work, i, k));
            slong t = c[pc];
            c[pc] = c[k];
            c[k] = t;
            scale = nmod_neg(scale, mod);
        }
        const nmod_poly_struct *pivot = nmod_poly_mat_entry(work, k, k);
        FLINT_ASSERT(nmod_poly_degree(pivot) == 0);
        ulong p = nmod_poly_get_coeff_ui(pivot, 0), inv = n_invmod(p, B->modulus);
        scale = nmod_mul(scale, p, mod);
        /* No polynomial division and no explicit inverse of E. Each worker
         * owns a row and its scratch polynomials; the pivot row is read-only. */
#pragma omp parallel if (n - k >= 64)
        {
            nmod_poly_t multiplier, product;
            nmod_poly_init(multiplier, B->modulus);
            nmod_poly_init(product, B->modulus);
#pragma omp for schedule(static)
            for (slong i = k + 1; i < n; i++)
            {
                if (nmod_poly_is_zero(nmod_poly_mat_entry(work, i, k)))
                    continue;
                nmod_poly_scalar_mul_nmod(multiplier, nmod_poly_mat_entry(work, i, k), inv);
                for (slong j = k + 1; j < n; j++)
                {
                    if (nmod_poly_is_zero(nmod_poly_mat_entry(work, k, j)))
                        continue;
                    nmod_poly_mul(product, multiplier, nmod_poly_mat_entry(work, k, j));
                    nmod_poly_sub(nmod_poly_mat_entry(work, i, j), nmod_poly_mat_entry(work, i, j),
                                  product);
                }
                nmod_poly_zero(nmod_poly_mat_entry(work, i, k));
            }
            nmod_poly_clear(multiplier);
            nmod_poly_clear(product);
        }
    }
    if (ok)
    {
        for (slong i = 0; i < h; i++)
            for (slong j = 0; j < h; j++)
                nmod_poly_set(nmod_poly_mat_entry(core, i, j),
                              nmod_poly_mat_entry(work, e + i, e + j));
        *factor = scale;
    }
    nmod_poly_mat_clear(work);
    flint_free(r);
    flint_free(c);
    return ok;
}

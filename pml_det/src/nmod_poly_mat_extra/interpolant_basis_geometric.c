/*
    Copyright (C) 2026 Gilles Villard

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
/* Optional, per-call profiling; no global state across concurrent calls. */
static double profile_now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9*t.tv_nsec;
}
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_vec.h>

#include "nmod_extra.h"
/* Bundled from PML 22bc6b8166854052370831f4fe4cd5fbdad69eb3.
 * Local adaptation: use drsolve's existing matrix multiplication dispatcher. */
#include "nmod_poly_mat_multiply.h"
#include "nmod_poly_mat_interpolant.h"
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"

/* Same D&C shape as nmod_poly_mat_pmintbasis (interpolant_basis.c), but
   points are rho_start * rho^k (k=0..d-1), rho = r^2, using one shared,
   already-built G (see nmod_poly_mat_pmintbasis_geometric below).
   rho_start is passed down as a single field element rather than an
   index, so no point array needs to be materialized except at the base
   case. */
static void _pmintbasis_geometric_rec(nmod_poly_mat_t intbas,
                                      slong * shift,
                                      const nmod_mat_struct * E,
                                      slong n,
                                      ulong r,
                                      ulong rho_start,
                                      slong d,
                                      nmod_t mod,
                                      const nmod_geometric_progression_t G, double *prof)
{
    double t = prof ? profile_now() : 0;
    ulong rho = nmod_mul(r, r, mod);
    const slong m = intbas->r;

    if (d <= PMINTBASIS_THRES)
    {
        ulong * pts = (ulong *) flint_malloc(d * sizeof(ulong));
        ulong cur = rho_start;
        for (slong k = 0; k < d; k++)
        {
            pts[k] = cur;
            cur = nmod_mul(cur, rho, mod);
        }
        nmod_poly_mat_mintbasis(intbas, shift, pts, E, d);
        if (prof) prof[0] += profile_now()-t;
        flint_free(pts);
        return;
    }

    const slong d1 = (d + 1) / 2;
    const slong d2 = d - d1;

    /* No truncated copy of E needed for the first half -- E may be longer
       than d1, extra entries simply unused. */
    nmod_poly_mat_t P1;
    nmod_poly_mat_init(P1, m, m, intbas->modulus);
    _pmintbasis_geometric_rec(P1, shift, E, n, r, rho_start, d1, mod, G, prof);
    if (prof) t = profile_now();

     /* R_i = P1(zeta*rho^i) * E_{d1+i}, i=0..d2-1, zeta = rho_start*rho^d1
       (the absolute start of the second half's points). Evaluating P1's
       entries at zeta*rho^i is done by rescaling each entry's coefficient
       k by zeta^k (Q(x) := entry(zeta*x)) and then evaluating the rescaled
       polynomial at the *first* d2 points of the shared G -- no offset
       handling needed inside G itself. */
    ulong zeta = nmod_mul(rho_start, nmod_pow_ui(rho, d1, mod), mod);

    /* a row of an order-d1 interpolant basis can have degree up to d1
       (length up to d1+1): the trivial witness M(x)*e_i, where M is the
       degree-d1 polynomial vanishing at all d1 points, is always a valid
       (if non-minimal) solution, so minimality never needs to exceed this
       -- same role X^d1*e_i plays for approximants, just a different
       degree-d1 witness since X^d1 alone isn't in the interpolant module
       (evaluation doesn't vanish automatically the way truncation does). */
    ulong * zeta_pow = (ulong *) flint_malloc((d1 + 1) * sizeof(ulong));
    zeta_pow[0] = 1;
    for (slong k = 1; k <= d1; k++)
        zeta_pow[k] = nmod_mul(zeta_pow[k - 1], zeta, mod);

    ulong * P1vals = (ulong *) flint_malloc(m * m * d2 * sizeof(ulong));
    /* Entry evaluations are independent. Keep small recursion nodes serial. */
#ifdef _OPENMP
    #pragma omp parallel if(!omp_in_parallel() && d2 >= 256 && m*m >= 8)
#endif
    {
        ulong * scaled = (ulong *) flint_malloc((d1 + 1) * sizeof(ulong));
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (slong ij = 0; ij < m*m; ij++)
        {
            const nmod_poly_struct * entry = nmod_poly_mat_entry(P1, ij/m, ij%m);
            slong len = nmod_poly_length(entry);
            /* Write directly into the final buffer; no per-entry copy. */
            ulong *vs = P1vals + ij*d2;
            if (len <= 1)
            {
                ulong c = len ? entry->coeffs[0] : 0;
                for (slong idx=0; idx<d2; idx++) vs[idx]=c;
            }
            else
            {
                for (slong k = 0; k < len; k++)
                    scaled[k] = nmod_mul(entry->coeffs[k], zeta_pow[k], mod);
                _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(vs, scaled, len, G, d2, mod);
            }
        }
        flint_free(scaled);
    }
    flint_free(zeta_pow);

    if (prof) { prof[1] += profile_now()-t; t=profile_now(); }
    /* R, the "new E" for the second half's recursive call, is a plain
       array too -- nmod_mat_mul writes each R_i directly. */
    nmod_mat_struct * R = (nmod_mat_struct *) flint_malloc(FLINT_MAX(d2, 1) * sizeof(nmod_mat_struct));
#ifdef _OPENMP
    #pragma omp parallel if(!omp_in_parallel() && d2 >= 256)
#endif
    {
        nmod_mat_t evalP1;
        nmod_mat_init(evalP1, m, m, intbas->modulus);
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (slong idx = 0; idx < d2; idx++)
        {
            for (slong i = 0; i < m; i++)
                for (slong j = 0; j < m; j++)
                    nmod_mat_entry(evalP1, i, j) = P1vals[(i * m + j) * d2 + idx];
            nmod_mat_init(R + idx, m, n, intbas->modulus);
            nmod_mat_mul(R + idx, evalP1, E + (d1 + idx));
        }
        nmod_mat_clear(evalP1);
    }
    flint_free(P1vals);

    if (prof) prof[2] += profile_now()-t;
    nmod_poly_mat_t P2;
    nmod_poly_mat_init(P2, m, m, intbas->modulus);
    _pmintbasis_geometric_rec(P2, shift, R, n, r, zeta, d2, mod, G, prof);
    for (slong idx = 0; idx < d2; idx++)
        nmod_mat_clear(R + idx);
    flint_free(R);

    if (prof) t=profile_now();
    drsolve_nmod_poly_mat_mul(intbas, P2, P1);
    if (prof) prof[3] += profile_now()-t;

    nmod_poly_mat_clear(P1);
    nmod_poly_mat_clear(P2);
}

void nmod_poly_mat_pmintbasis_geometric(nmod_poly_mat_t intbas,
                                        slong * shift,
                                        ulong * pts,
                                        const nmod_mat_struct * E,
                                        ulong r,
                                        slong d)
{
    /* d=0: no points, output is the identity, shift unchanged -- matching
       the same convention as every other function in this project TO SEE. */
    if (d == 0)
    {
        nmod_poly_mat_one(intbas);
        return;
    }

    if (r == 0)
    flint_throw(FLINT_ERROR,
                "Exception (nmod_poly_mat_pmintbasis_geometric). "
                "r must be nonzero (r=0 degenerates FLINT's geometric-"
                "progression setup, which computes 1/r unconditionally).\n");

    const slong n = E[0].c;

    nmod_t mod;
    nmod_init(&mod, intbas->modulus);

    nmod_geometric_progression_t G;
    _nmod_geometric_progression_init_function(G, r, d, mod, 1); /* 1 = evaluation only */

    double prof[4]={0};
    int profiling = getenv("DRSOLVE_PROFILE_INTERPOLANT") != NULL;
    _pmintbasis_geometric_rec(intbas, shift, E, n, r, 1, d, mod, G, profiling ? prof : NULL);
    if (profiling)
        fprintf(stderr, "interpolant: base %.6f evaluation %.6f residual %.6f multiply %.6f\n",
                prof[0], prof[1], prof[2], prof[3]);

    nmod_geometric_progression_clear(G);

    if (pts)
    {
        ulong rho = nmod_mul(r, r, mod);
        ulong cur = 1;
        for (slong k = 0; k < d; k++)
        {
            pts[k] = cur;
            cur = nmod_mul(cur, rho, mod);
        }
    }
}

/* Local boundary fix: accept every mathematically available length,
 * including d=(p-1)/2 and the single-point case in characteristic two. */
void nmod_poly_mat_pmintbasis_geometric_auto(nmod_poly_mat_t basis,
    slong *shift, ulong *points, const nmod_mat_struct *values, slong count)
{
    if(count==0){nmod_poly_mat_one(basis);return;}
    if(count==1){nmod_poly_mat_pmintbasis_geometric(basis,shift,points,values,1,1);return;}
    if((ulong)count>(basis->modulus-1)/2)
        flint_throw(FLINT_ERROR,"Not enough distinct square-geometric points.\n");
    ulong root=n_primitive_root_prime(basis->modulus);
    nmod_poly_mat_pmintbasis_geometric(basis,shift,points,values,root,count);
}

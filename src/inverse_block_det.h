/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_INVERSE_BLOCK_DET_H
#define DRSOLVE_INVERSE_BLOCK_DET_H
#include <flint/nmod_poly_mat.h>
#include <flint/nmod_poly.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    slong degree_bound, block, points, attempts, checks;
    double setup, evaluation, inverse, interpolation, determinant, verification;
} dr_inverse_block_stats;

/* A priori heuristic; threads is the number of available sampling workers.
 * Returns a block in [1,n] (0 for invalid dimensions). No benchmark is run. */
slong dr_inverse_block_choose(slong n, slong degree_bound, int threads);

/* Maximum assignment weight of the entry degrees (zeros padded with degree
 * zero). This is a certified upper bound, not an estimate of det's degree.
 * Returns -1 for invalid dimensions or arithmetic overflow. */
slong dr_nmod_poly_mat_det_degree_bound(const nmod_poly_mat_t matrix);

/* Random inverse projection and matrix fraction reconstruction over a word
 * prime field. block=0 selects automatically; otherwise 1<=block<=rows.
 * threads<=0 uses the current OpenMP setting. Inputs remain unchanged and
 * failure leaves out unchanged, allowing the caller to fall back to HNF.
 * Success returns the exact scalar normalization, probabilistically checked
 * at fresh uniform points. At most two candidates, with combined false
 * acceptance probability <= 2^-64 under independent uniform randomness.
 * Needs PML, enough square-geometric sample points, and p>2*degree_bound for
 * base-field verification; unsupported cases return 0. This is not a proof
 * certificate. No genericity assumption is required for the fallback.
 * state supplies all projections, dilations, and verification points.
 * stats may be NULL; all timings are wall time and accumulate over retries. */
int dr_nmod_poly_mat_det_inverse_block_try(nmod_poly_t out,
    const nmod_poly_mat_t matrix, slong block, int threads,
    flint_rand_t state, dr_inverse_block_stats *stats);
#ifdef __cplusplus
}
#endif
#endif

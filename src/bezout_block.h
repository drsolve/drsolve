/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_BEZOUT_BLOCK_H
#define DRSOLVE_BEZOUT_BLOCK_H
#include <flint/nmod_mpoly.h>
#include <flint/nmod_poly.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Process-wide setting: positive = initial point count, 0 = automatic,
 * negative = consult DRSOLVE_BEZOUT_INITIAL_POINTS (the default).
 * Configure before starting computations; does not disable verification. */
void dr_bezout_block_set_initial_points(slong points);

typedef struct {
    slong points, block, attempts, verification_extension_degree, verification_checks;
    double evaluation, inverse, interpolation, determinant, verification;
} dr_bezout_block_stats;
/* Dixon/Bezout inverse-Hankel projection with verified reconstruction.
 * Equal total degrees d: constant, nonzero x^d coefficients are required;
 * block=0 selects max(min(6,d),floor(cuberoot(d))). The contiguous projection
 * starts with O(d) points, doubling up to the total-degree reconstruction
 * budget. DRSOLVE_BEZOUT_INITIAL_POINTS optionally overrides the equal-degree
 * starting count (bounded by that budget); invalid values use the default.
 * DRSOLVE_PROFILE_BEZOUT=1 reports each attempt to stderr.
 * Every candidate is independently verified before it is accepted.
 * Unequal total degrees A>B: the smaller polynomial must have nonzero
 * constant x^B coefficient. Reduce the larger modulo its monic normalization
 * and use a spaced projection of the B-by-B Bezout inverse. block=0 selects
 * floor(cuberoot(B)); otherwise 1 <= block <= min(A,B).
 * var selects x; result is in the other variable, with its exact scalar
 * (including input-order sign and leading-coefficient factors).
 * Returns 1 after probabilistic verification, 0 to request fallback; failure
 * leaves out unchanged. Inputs are never modified. Prime nmod fields only.
 * For p >= 2*(A*B+1), use native prime-field checks with the least k such
 * that ((A*B+1)/p)^k <= 2^-80. Otherwise use eight fresh checks in an extension
 * of size >= 1024*(A*B+1). The scalar-recovery sample is not counted.
 * For uniform independent samples, the false acceptance bound per candidate
 * is <= 2^-80. This is not a certificate. Parallel checks use thread-local
 * scratch and samples drawn before entering the parallel region.
 * state supplies the random samples; stats may be NULL. The evaluation timing
 * includes polynomial degree reduction for unequal total degrees. */
int dr_nmod_mpoly_resultant_bezout_block_try(nmod_poly_t out,
    const nmod_mpoly_t a,const nmod_mpoly_t b,slong var,
    const nmod_mpoly_ctx_t ctx,slong block,flint_rand_t state,
    dr_bezout_block_stats *stats);
#ifdef __cplusplus
}
#endif
#endif

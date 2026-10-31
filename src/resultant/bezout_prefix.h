/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_BEZOUT_PREFIX_H
#define DRSOLVE_BEZOUT_PREFIX_H
#include <flint/nmod_poly.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Word-size prime fields; no NTL/PML dependency. */
typedef struct {
    slong euclidean_steps, hgcd_blocks, full_inverse_fallbacks;
    double setup, hgcd, division, update, fallback, cleanup, series;
} dr_bezout_prefix_stats;
/* Return 1 on success, 0 if noninvertible, -1 for invalid arguments.
 * Output is reverse(g^{-1} mod f, deg(f)-1) mod x^precision.
 * cutoff=0 selects classical Euclid; positive cutoff enables half-GCD.
 * The cutoff selects an algorithm; it does not restrict input sizes.
 * Output is unchanged on failure; aliases with f/g are supported. */
int dr_bezout_inverse_prefix(nmod_poly_t out, const nmod_poly_t f,
    const nmod_poly_t g, slong precision, slong cutoff, dr_bezout_prefix_stats *stats);
/* Coefficients h[i+j] give entries of the monomial Bezout inverse. */
int dr_bezout_left_moments(nmod_poly_t out, const nmod_poly_t f,
    const nmod_poly_t g, slong precision, slong cutoff, dr_bezout_prefix_stats *stats);
/* Internal sampled profiling: timed enables per-operation clocks. */
int dr_bezout_left_moments_profile(nmod_poly_t out,const nmod_poly_t f,
    const nmod_poly_t g,slong precision,slong cutoff,dr_bezout_prefix_stats *stats,int timed);
#ifdef __cplusplus
}
#endif
#endif

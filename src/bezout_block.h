/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_BEZOUT_BLOCK_H
#define DRSOLVE_BEZOUT_BLOCK_H
#include <flint/nmod_mpoly.h>
#include <flint/nmod_poly.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct {
    slong points, block, attempts, verification_extension_degree;
    double evaluation, inverse, interpolation, determinant, verification;
} dr_bezout_block_stats;
/* Experimental Dixon/Bezout inverse-Hankel projection and shifted interpolation.
 * Currently requires two total-degree-d bivariates with constant, nonzero x^d
 * coefficients. var selects x; result is in the other variable, with its exact
 * scalar (not monic). block=0 selects min(6,d); otherwise 1 <= block <= d.
 * Returns 1 after probabilistic verification, 0 to request fallback; failure
 * leaves out unchanged. Inputs are never modified. Prime nmod fields only.
 * Eight fresh checks use a field of size >= 256*(d^2+1), extending the field
 * for verification only if needed. For uniform independent samples, the false
 * acceptance bound per candidate is <= 2^-64. This is not a certificate.
 * state supplies the random samples; stats may be NULL. */
int dr_nmod_mpoly_resultant_bezout_block_try(nmod_poly_t out,
    const nmod_mpoly_t a,const nmod_mpoly_t b,slong var,
    const nmod_mpoly_ctx_t ctx,slong block,flint_rand_t state,
    dr_bezout_block_stats *stats);
#ifdef __cplusplus
}
#endif
#endif

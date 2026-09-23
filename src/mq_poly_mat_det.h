/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_MQ_POLY_MAT_DET_H
#define DRSOLVE_MQ_POLY_MAT_DET_H
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Blocked prime-field MQ compression. B has its h retained indices first
 * on BOTH axes; the remaining complementary block must have degree budget 0.
 * core is an initialized h by h matrix over B's prime field and must not
 * overlap B's storage (including matrix windows). On success,
 * det(B) = *factor * det(core), with factor nonzero. On failure, B, core and
 * *factor are unchanged. Input degrees include all eliminated monomial weights.
 * No generic-rank assumption is used. Caller must supply a prime modulus.
 * Constant degree-diagonal blocks are factored once; block back substitution
 * computes E X = V and core = A - U X without a full copy or inverse of E.
 * The solver attempts this by default when eligible; --no-mq-step4-schur disables it.
 */
int nmod_poly_mat_mq_schur(nmod_poly_mat_t core, ulong *factor, const nmod_poly_mat_t B,
                          const slong *row_degree, const slong *col_degree, slong h, slong sigma);
#ifdef __cplusplus
}
#endif
#endif

/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_QUOTIENT_SOLVER_H
#define DRSOLVE_QUOTIENT_SOLVER_H

#include "polynomial_system_solver.h"

/* Certified finite-field backend; prime fields retain the nmod fast path.
 * On failure returns 0 and sets sols->error_message; it never reports a
 * degree/resource limit as absence of solutions or positive dimension.
 * On success returns every point rational over sols->ctx, without multiplicity.
 * A zero max_degree or memory_mb disables that optional limit.
 */
int solve_by_quotient_closure(char **polys, slong count,
                            variable_info_t *vars, slong n,
                            polynomial_solutions_t *sols,
                            slong max_degree, slong memory_mb);

/* Extension backend, using the existing native binary / Zech / fq dispatch. */
int solve_by_quotient_closure_extension(char **polys, slong count,
                            variable_info_t *vars, slong n,
                            polynomial_solutions_t *sols,
                            slong max_degree, slong memory_mb);

/* Exact, bounded factor cancellation using the quotient relation space.
 * All variables (including parameters) belong to ctx; resultant must involve
 * parameters only, as checked by the caller. Each accepted candidate lies in
 * <polys>, even for positive-dimensional ideals. No point enumeration,
 * field equations, specialization, or radicalization is used.
 * Returns 1 after a bounded search, 0 on invalid input/factorization failure.
 * Unproved cancellations are retained: the output need not be minimal.
 * max_degree must be positive. There is no memory budget; only addressable
 * matrix dimensions are checked. Outputs may alias
 * neither each other nor the input. used_degree is the last constructed degree.
 */
#include <flint/nmod_mpoly.h>
int quotient_sift_resultant(nmod_mpoly_t sifted, nmod_mpoly_t removed,
                           const nmod_mpoly_t resultant,
                           nmod_mpoly_struct *polys, slong count,
                           const nmod_mpoly_ctx_t ctx,
                           slong max_degree,
                           slong *used_degree);

/* Same bounded ideal-membership certificate for extension coefficients.
 * The original fq context, including its defining polynomial, is preserved.
 * Closure arithmetic uses the native binary / general field dispatch.
 */
#include <flint/fq_nmod_mpoly.h>
int quotient_sift_resultant_extension(fq_nmod_mpoly_t sifted,
                           fq_nmod_mpoly_t removed,
                           const fq_nmod_mpoly_t resultant,
                           fq_nmod_mpoly_struct *polys, slong count,
                           const fq_nmod_mpoly_ctx_t ctx,
                           slong max_degree, slong *used_degree);

#endif

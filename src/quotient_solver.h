/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_QUOTIENT_SOLVER_H
#define DRSOLVE_QUOTIENT_SOLVER_H

#include "polynomial_system_solver.h"

/* Prime-field backend. Does not print or change global arithmetic settings.
 * On failure returns 0 and sets sols->error_message; it never reports a
 * degree/resource limit as absence of solutions or positive dimension.
 * On success returns every F_p-rational point, without multiplicity.
 * A zero max_degree or memory_mb disables that optional limit.
 */
int solve_by_quotient_closure(char **polys, slong count,
                            variable_info_t *vars, slong n,
                            polynomial_solutions_t *sols,
                            slong max_degree, slong memory_mb);

#endif

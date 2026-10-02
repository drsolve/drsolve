/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef MACAULAY_FLINT_H
#define MACAULAY_FLINT_H

#include <flint/flint.h>
#include "dr_mpoly.h"
#include "dixon_flint.h"

void fq_macaulay_resultant(unified_mpoly_struct *result, unified_mpoly_struct *polys, slong nvars,
                           slong npars);

void fq_macaulay_resultant_with_names(unified_mpoly_struct *result, unified_mpoly_struct *polys,
                                      slong nvars, slong npars, char **var_names, char **par_names,
                                      const char *gen_name);

#endif

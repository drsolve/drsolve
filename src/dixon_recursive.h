/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DIXON_RECURSIVE_H
#define DIXON_RECURSIVE_H

#include <flint/flint.h>

#include <stdarg.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif
#include "dr_mpoly.h"
#include "dixon_flint.h"

void fq_dixon_fast_resultant(unified_mpoly_struct *result, unified_mpoly_struct *polys, slong nvars,
                             slong npars);

void fq_dixon_fast_resultant_with_names(unified_mpoly_struct *result, unified_mpoly_struct *polys,
                                        slong nvars, slong npars, char **var_names,
                                        char **par_names, const char *gen_name);

#endif

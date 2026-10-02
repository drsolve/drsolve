/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DR_MPOLY_H
#define DR_MPOLY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <unistd.h>
#include <ctype.h>

#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_poly.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_mpoly.h>
#include <flint/nmod_poly.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/profiler.h>

#include "unified_mpoly_interface.h"
#include "fq_unified_interface.h"

/* ============================================================================
 * Data Structures for Multivariate Polynomials over Finite Extension Fields
 * ============================================================================ */

/* Stack-only view of one native term. Never free its coefficient or retain
 * the view across mutation of the polynomial. The exponent buffer is local. */
typedef struct {
    slong *var_exp, *par_exp;
    fq_nmod_t coeff;
} dr_mpoly_term_view;

static inline slong dr_mpoly_length(const unified_mpoly_struct *p)
{
    if (!p->ctx_ptr)
        return 0;
    if (p->field_id == FIELD_ID_NMOD)
        return p->data.nmod_poly.length;
    if (p->field_id == FIELD_ID_FQ_ZECH)
        return p->data.zech_poly.length;
    return p->data.fq_poly.length;
}

void dr_mpoly_read_term(dr_mpoly_term_view *view, slong *exp, ulong *coefficient,
                        const unified_mpoly_struct *p, slong i);
#define DR_MPOLY_TERM(name, poly, index)                                                           \
    slong name##_exponents[FLINT_MAX(1, (poly)->nvars + (poly)->npars)];                           \
    ulong name##_coefficient[FLINT_MAX(1, fq_nmod_ctx_degree((poly)->ctx))];                       \
    dr_mpoly_term_view name;                                                                       \
    dr_mpoly_read_term(&name, name##_exponents, name##_coefficient, (poly), (index))

/* Solver polynomials use unified_mpoly_struct's native FLINT storage. There is
 * no separately allocated term representation. Prime-field kernels may borrow
 * GET_NMOD_POLY directly; extension arithmetic retains unified's optimizations.
 *
 * Objects start zero-initialized. init resets an object, arithmetic replaces its
 * output (aliasing is supported), and clear is safe on a zero object. The caller's
 * fq_nmod context must outlive its polynomials. move transfers ownership.
 *
 * Fast insertion and exponent mutation leave terms unnormalized. Normalize once
 * after construction, before sharing a polynomial among parallel readers or
 * passing its native storage to a FLINT operation requiring canonical input.
 * Generic arithmetic normalizes unfinished inputs lazily. */
void dr_mpoly_normalize(unified_mpoly_struct *p);
void dr_mpoly_fit_length(unified_mpoly_struct *p, slong length);
void dr_mpoly_move(unified_mpoly_struct *out, unified_mpoly_struct *in);
void dr_mpoly_set_term_exp(unified_mpoly_struct *p, slong i, const slong *exp);
void dr_mpoly_neg(unified_mpoly_struct *out, const unified_mpoly_struct *in);

// Comparison structure for sorting degrees
typedef struct {
    slong index;
    slong degree;
} fq_index_degree_pair;

/* ============================================================================
 * Field Equation Reduction
 * ============================================================================ */

/* Global flag: when non-zero, dr_mpoly_mul reduces each variable modulo
 * x^q - x after every multiplication (working in F_q[vars] / <x^q - x, ...>). */
extern int g_field_equation_reduction;

/* Global flag: when non-zero, only final output/resultants are reduced modulo
 * x^q - x, while intermediate arithmetic stays in the ordinary polynomial ring. */
extern int g_field_equation_final_only;

/* Enable or disable field-equation reduction mode. */
void dr_mpoly_set_field_equation_reduction(int enable);

/* Enable or disable final-result-only field-equation reduction mode. */
void dr_mpoly_set_field_equation_final_only(int enable);

/* Reduce all variable exponents of poly in-place using x^q = x. */
void dr_mpoly_reduce_field_equation(unified_mpoly_struct *poly);

/* ============================================================================
 * Basic Operations
 * ============================================================================ */

// Initialize a multivariate polynomial
void dr_mpoly_init(unified_mpoly_struct *p, slong nvars, slong npars, const fq_nmod_ctx_t ctx);

// Clear a multivariate polynomial and free memory
void dr_mpoly_clear(unified_mpoly_struct *p);

// Copy one polynomial to another
void dr_mpoly_copy(unified_mpoly_struct *dest, const unified_mpoly_struct *src);

/* ============================================================================
 * Term Management
 * ============================================================================ */

// Add a term to the polynomial with duplicate checking
void dr_mpoly_add_term(unified_mpoly_struct *p, const slong *var_exp, const slong *par_exp,
                       const fq_nmod_t coeff);

// Add a term to the polynomial without duplicate checking (faster)
void dr_mpoly_add_term_fast(unified_mpoly_struct *p, const slong *var_exp, const slong *par_exp,
                            const fq_nmod_t coeff);

/* ============================================================================
 * Display Functions
 * ============================================================================ */

// Print polynomial in standard format
void dr_mpoly_print(const unified_mpoly_struct *p, const char *name);

// Print polynomial with expanded dual variable notation
void dr_mpoly_print_expanded(const unified_mpoly_struct *p, const char *name, int use_dual);
void dr_mpoly_print_with_names(const unified_mpoly_struct *poly, const char *poly_name,
                               char **var_names, char **par_names, const char *gen_name,
                               int expanded_format);
/* ============================================================================
 * Arithmetic Operations
 * ============================================================================ */

// Multiply two multivariate polynomials
void dr_mpoly_mul(unified_mpoly_struct *result, const unified_mpoly_struct *a,
                  const unified_mpoly_struct *b);

// Compute polynomial to a power
void dr_mpoly_pow(unified_mpoly_struct *result, const unified_mpoly_struct *base, slong power);

// Multiply polynomial by scalar
void dr_mpoly_scalar_mul(unified_mpoly_struct *result, const unified_mpoly_struct *p,
                         const fq_nmod_t scalar);

// Add two polynomials
void dr_mpoly_add(unified_mpoly_struct *result, const unified_mpoly_struct *a,
                  const unified_mpoly_struct *b);

// Subtract two polynomials
void dr_mpoly_sub(unified_mpoly_struct *result, const unified_mpoly_struct *a,
                  const unified_mpoly_struct *b);

// Make polynomial monic (leading coefficient = 1)
void dr_mpoly_make_monic(unified_mpoly_struct *poly);

/* ============================================================================
 * Conversion Functions
 * ============================================================================ */

// Convert dr_mpoly to FLINT's fq_nmod_mpoly format
void dr_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t mpoly, const unified_mpoly_struct *poly,
                               fq_nmod_mpoly_ctx_t mpoly_ctx);

// Convert FLINT's fq_nmod_mpoly to dr_mpoly format
void fq_nmod_mpoly_to_dr_mpoly(unified_mpoly_struct *poly, const fq_nmod_mpoly_t mpoly, slong nvars,
                               slong npars, fq_nmod_mpoly_ctx_t mpoly_ctx, const fq_nmod_ctx_t ctx);

// Convert dr_mpoly to nmod_mpoly for prime fields
void dr_mpoly_to_nmod_mpoly(nmod_mpoly_t mpoly, const unified_mpoly_struct *poly,
                            nmod_mpoly_ctx_t mpoly_ctx);

// Convert nmod_mpoly to dr_mpoly
void nmod_mpoly_to_dr_mpoly(unified_mpoly_struct *result, const nmod_mpoly_t poly, slong nvars,
                            slong npars, const nmod_mpoly_ctx_t mpoly_ctx,
                            const fq_nmod_ctx_t field_ctx);

/* Transfer a native result without copying its coefficient/exponent arrays. */
void dr_mpoly_take_nmod(unified_mpoly_struct *out, nmod_mpoly_t in, slong nvars, slong npars,
                        const nmod_mpoly_ctx_t ctx, const fq_nmod_ctx_t field);

/* ============================================================================
 * Division Operations
 * ============================================================================ */

// Divide by linear factor using FLINT's built-in functions
void dr_mpoly_divide_difference(unified_mpoly_struct *quotient,
                                const unified_mpoly_struct *dividend, slong var_idx, slong nvars,
                                slong npars);

/* ============================================================================
 * Evaluation Functions
 * ============================================================================ */

// Evaluate polynomial at given parameter values
void evaluate_dr_mpoly_at_params(fq_nmod_t result, const unified_mpoly_struct *poly,
                                 const fq_nmod_t *param_vals);

/* ============================================================================
 * Matrix Display and Analysis Functions
 * ============================================================================ */

// Print matrix of polynomials with optional details
void dr_mpoly_matrix_print(unified_mpoly_struct **matrix, slong nrows, slong ncols,
                           const char *matrix_name, int show_details);

/* ============================================================================
 * Comparison Function
 * ============================================================================ */

// Comparison function for sorting degrees
int compare_fq_degrees(const void *a, const void *b);

#endif /* DR_MPOLY_H */

/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_UNIFIED_MAT_POLY_H
#define DRSOLVE_UNIFIED_MAT_POLY_H

#include "fq_unified_interface.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Polynomial with matrix coefficients. Coefficient k is a row-major r*c
 * block. Native binary fields use compact coefficients; other fields use field_elem_u.
 * The context is borrowed and must outlive the matrix. Do not access data
 * through unified_poly_struct::coeffs: these are different storage formats. */
typedef struct {
    void *data;
    slong r, c, length, alloc;
    field_ctx_t *ctx;
} unified_mat_poly_struct;
typedef unified_mat_poly_struct unified_mat_poly_t[1];

void unified_mat_poly_init(unified_mat_poly_t a, slong r, slong c, field_ctx_t *ctx);
void unified_mat_poly_clear(unified_mat_poly_t a);
int unified_mat_poly_is_byte(const unified_mat_poly_t a);
void unified_mat_poly_from_fq(unified_mat_poly_t a, const fq_nmod_poly_mat_t b,
                              const fq_nmod_ctx_t ctx);
void unified_mat_poly_to_fq(fq_nmod_poly_mat_t a, const unified_mat_poly_t b,
                            const fq_nmod_ctx_t ctx);
void unified_mat_poly_from_unified(unified_mat_poly_t a, const unified_poly_mat_t b);
void unified_mat_poly_to_unified(unified_poly_mat_t a, const unified_mat_poly_t b);
/* Outputs must have the indicated shape and the same context. Aliasing is
 * supported. A window returns coefficients [start,start+length), shifted
 * down by start. No field-equation reduction is performed. */
void unified_mat_poly_mul(unified_mat_poly_t out, const unified_mat_poly_t a,
                          const unified_mat_poly_t b);
void unified_mat_poly_mul_window(unified_mat_poly_t out, const unified_mat_poly_t a,
                                 const unified_mat_poly_t b, slong start, slong length);
/* Left approximant basis; shift is updated to shifted row degrees.
 * Output shape is a->r by a->r. order >= 0. */
void unified_mat_poly_pmbasis(unified_mat_poly_t out, slong *shift,
                              const unified_mat_poly_t a, slong order);
/* Exact left kernel basis for a full-column-rank tall matrix. Output shape
 * (a->r-a->c) by a->r. Returns 0 for rank-deficient/unsupported input;
 * callers must not interpret failure as a zero determinant. */
int unified_mat_poly_kernel(unified_mat_poly_t out, const unified_mat_poly_t a);
/* Exact recursive HNF-style determinant, no random evaluation shortcuts.
 * Returns 0 on a structural failure, leaving det unchanged. */
int unified_mat_poly_det_hnf(unified_poly_t det, const unified_mat_poly_t a);
int unified_fq_poly_mat_det_hnf(fq_nmod_poly_t det, const fq_nmod_poly_mat_t a,
                               const fq_nmod_ctx_t ctx);
#ifdef __cplusplus
}
#endif
#endif

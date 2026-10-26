/*
 * Minimal public header for the drsolve-bundled PML determinant subset.
 *
 * This header intentionally exposes only the polynomial-matrix determinant
 * entry points and finite-field interpolation bases used by drsolve.
 */
#ifndef DRSOLVE_PML_DET_H
#define DRSOLVE_PML_DET_H

#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Optional FLINT FFT backend; disable to retain the original dispatcher. */
void drsolve_pml_set_fft_enabled(int enabled);
int drsolve_pml_fft_available(void);


void nmod_poly_mat_det_iter(nmod_poly_t det, nmod_poly_mat_t mat);

int nmod_poly_mat_det_hnf(nmod_poly_t det,
                          const nmod_poly_mat_t mat);

int nmod_poly_mat_det_hnf_knowing_degree(nmod_poly_t det,
                                         const nmod_poly_mat_t mat,
                                         slong degree);

/* Input values[k] is an r-by-c matrix at points[k]; basis is r-by-r.
 * shift has r entries, updated to the output shifted row degrees.
 * For the geometric variant, points[k] = ratio^(2*k); points may be NULL.
 * Inputs are not modified; the caller supplies valid dimensions and a prime
 * modulus. Distinct points give the usual sum-of-ranks determinant degree. */
void nmod_poly_mat_mintbasis(nmod_poly_mat_t basis, slong *shift,
    const ulong *points, const nmod_mat_struct *values, slong count);
void nmod_poly_mat_pmintbasis(nmod_poly_mat_t basis, slong *shift,
    const ulong *points, const nmod_mat_struct *values, slong count);
void nmod_poly_mat_pmintbasis_geometric(nmod_poly_mat_t basis, slong *shift,
    ulong *points, const nmod_mat_struct *values, ulong ratio, slong count);

/* Chooses a primitive-root ratio; supports count <= (p-1)/2, or count <= 1. */
void nmod_poly_mat_pmintbasis_geometric_auto(nmod_poly_mat_t basis, slong *shift,
    ulong *points, const nmod_mat_struct *values, slong count);

#ifdef __cplusplus
}
#endif

#endif /* DRSOLVE_PML_DET_H */

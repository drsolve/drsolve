#ifndef DRSOLVE_SD_FFT_BACKEND_H
#define DRSOLVE_SD_FFT_BACKEND_H
#include <flint/nmod_poly_mat.h>
void drsolve_nmod_poly_mat_mul_sd_fft_direct(nmod_poly_mat_t C,
    const nmod_poly_mat_t A, const nmod_poly_mat_t B);
void drsolve_nmod_poly_mat_mulmid_sd_fft_direct_bounded(nmod_poly_mat_t C,
    const nmod_poly_mat_t A, slong alen, const nmod_poly_mat_t B, slong blen,
    slong lo, slong hi, ulong unused);
void drsolve_nmod_poly_mat_mulmid_sd_fft_direct(nmod_poly_mat_t C,
    const nmod_poly_mat_t A, slong alen, const nmod_poly_mat_t B, slong blen,
    slong lo, slong hi);
#endif

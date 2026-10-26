/* Compile AND link this probe: development snapshots can share a version. */
#include <flint/flint.h>
#include <flint/fft_small.h>
#include <flint/machine_vectors.h>
#if __FLINT_RELEASE < 30700 || !defined(FLINT_HAVE_FFT_SMALL)
#error New FLINT FFT API required
#endif
static void (*volatile required[])(void) = {
    (void (*)(void)) fft_small_plan_init_nmod,
    (void (*)(void)) fft_small_plan_clear,
    (void (*)(void)) _fft_small_plan_set_window_nmod,
    (void (*)(void)) fft_small_op_init_borrowed,
    (void (*)(void)) fft_small_op_clear,
    (void (*)(void)) fft_small_op_sizeof_data,
    (void (*)(void)) fft_small_fft_nmod,
    (void (*)(void)) fft_small_fft_tr,
    (void (*)(void)) fft_small_ifft,
    (void (*)(void)) fft_small_ifft_tr_nmod,
    (void (*)(void)) fft_small_export_nmod_range
};
int main(void) {
    vec4d v = vec4d_set_d(0.0);
    (void) v;
    for (unsigned i=0; i<sizeof(required)/sizeof(required[0]); ++i)
        if (!required[i]) return 1;
    return 0;
}

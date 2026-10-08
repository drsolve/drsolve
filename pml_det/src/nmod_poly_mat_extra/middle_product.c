/*
    Copyright (C) 2025 Vincent Neiger, Éric Schost

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <flint/ulong_extras.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_extra.h"
#include "nmod_poly_mat_multiply.h"

/* A small prime-field NTT backend. Transforms are per entry; point values
 * are packed into FLINT matrices so the expensive inner products share them.
 * The cyclic transform is sufficient for a middle window: if L >= hi and
 * deg(A_trunc B_trunc) < L + start, no wrapped coefficient reaches it. */
static void drsolve_ntt(ulong *a, slong n, const ulong *roots, nmod_t mod)
{
    for (slong i = 1, j = 0; i < n; i++) {
        slong bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) { ulong t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    for (slong len = 2; len <= n; len <<= 1) {
        slong half = len >> 1, step = n / len;
        for (slong i = 0; i < n; i += len)
            for (slong j = 0; j < half; j++) {
                ulong u = a[i+j], v = nmod_mul(a[i+j+half], roots[j*step], mod);
                a[i+j] = nmod_add(u, v, mod);
                a[i+j+half] = nmod_sub(u, v, mod);
            }
    }
}

int drsolve_nmod_poly_mat_mul_window_ntt(nmod_poly_mat_t C,
    const nmod_poly_mat_t A, const nmod_poly_mat_t B, slong start, slong count)
{
    const slong m = A->r, k = A->c, n = B->c;
    if (start < 0 || count < 0 || start > WORD_MAX-count || k != B->r ||
        C->r != m || C->c != n || A->modulus != B->modulus || C->modulus != A->modulus)
        return 0;
    if (!count || !m || !k || !n) { nmod_poly_mat_zero(C); return 1; }
    const slong hi = start + count;
    slong alen = FLINT_MIN(hi, nmod_poly_mat_max_length(A));
    slong blen = FLINT_MIN(hi, nmod_poly_mat_max_length(B));
    if (!alen || !blen) { nmod_poly_mat_zero(C); return 1; }
    if (alen > WORD_MAX-blen) return 0;
    slong plen = alen + blen - 1;
    if (start >= plen) { nmod_poly_mat_zero(C); return 1; }
    slong need = FLINT_MAX(FLINT_MIN(hi, plen), plen-start), len = 1;
    /* Bound both transform length and total workspace before allocating.
     * This backend is optional; unsupported fields use FLINT unchanged. */
    while (len < need && len < 65536) len <<= 1;
    if (len < need || A->modulus < 3 || (A->modulus-1) % len) return 0;
    size_t limit = (size_t)256*1024*1024/sizeof(ulong);
    if ((size_t)m > limit/(size_t)k || (size_t)k > limit/(size_t)n ||
        (size_t)m > limit/(size_t)n) return 0;
    size_t na = (size_t)m*k, nb = (size_t)k*n, nc = (size_t)m*n;
    if (na+nb+nc+3 > limit/(size_t)(len+1)) return 0;
    /* The caller uses prime fields. Checking here also makes this standalone
     * entry point reject composite moduli instead of searching for a root. */
    if (!n_is_prime(A->modulus)) return 0;
    nmod_t mod; nmod_init(&mod, A->modulus);
    ulong root = 1;
    if (len > 1) {
        /* For prime p, a^((p-1)/L) has order L iff its L/2 power is -1. */
        for (ulong a = 2; a < mod.n; a++) {
            root = nmod_pow_ui(a, (mod.n-1)/len, mod);
            if (nmod_pow_ui(root, len/2, mod) == mod.n-1) break;
        }
    }
    ulong *roots = flint_malloc((size_t)len*sizeof(ulong));
    ulong *inverse = flint_malloc((size_t)len*sizeof(ulong));
    ulong *av = flint_malloc(na*len*sizeof(ulong));
    ulong *bv = flint_malloc(nb*len*sizeof(ulong));
    ulong *cv = flint_malloc(nc*len*sizeof(ulong));
    roots[0] = inverse[0] = 1;
    ulong invroot = nmod_inv(root, mod);
    for (slong i = 1; i < len; i++) {
        roots[i] = nmod_mul(roots[i-1],root,mod);
        inverse[i] = nmod_mul(inverse[i-1],invroot,mod);
    }
    /* Bound per-worker scratch as well as the shared transforms. Small
     * products and calls inside an existing team stay serial. */
#ifdef _OPENMP
    int workers = 1;
    if (!omp_in_parallel()) {
        /* Amortize team overhead with at least 32K coefficient slots
         * per worker; pointwise products add more work on wide matrices. */
        size_t jobs = (na+nb+nc)*(size_t)len / 32768;
        size_t scratch = na+nb+nc+(size_t)len;
        size_t available = limit-(na+nb+nc+2)*(size_t)len;
        workers = (int) FLINT_MAX(1, FLINT_MIN((size_t)omp_get_max_threads(),
            FLINT_MIN((size_t)len, FLINT_MIN(jobs, available/scratch))));
    }
#endif
    ulong scale = nmod_inv((ulong)len,mod);
    slong outlen = FLINT_MIN(count,plen-start);
#ifdef _OPENMP
    #pragma omp parallel num_threads(workers) if(workers > 1)
#endif
    {
        ulong *tmp = flint_malloc((size_t)len*sizeof(ulong));
        nmod_mat_t am, bm, cm;
        nmod_mat_init(am,m,k,mod.n);
        nmod_mat_init(bm,k,n,mod.n);
        nmod_mat_init(cm,m,n,mod.n);
        for (int which = 0; which < 2; which++) {
            const nmod_poly_mat_struct *M = which ? B : A;
            ulong *values = which ? bv : av;
            size_t stride = which ? nb : na;
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            for (slong entry = 0; entry < (slong)stride; entry++) {
                const nmod_poly_struct *f = nmod_poly_mat_entry(M,entry/M->c,entry%M->c);
                slong used = FLINT_MIN(f->length,hi);
                if (used) memcpy(tmp,f->coeffs,(size_t)used*sizeof(ulong));
                memset(tmp+used,0,(size_t)(len-used)*sizeof(ulong));
                drsolve_ntt(tmp,len,roots,mod);
                for (slong t = 0; t < len; t++) values[(size_t)t*stride+entry] = tmp[t];
            }
        }
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (slong t = 0; t < len; t++) {
            for (slong i = 0; i < m; i++)
                memcpy(nmod_mat_entry_ptr(am,i,0),av+(size_t)t*na+(size_t)i*k,(size_t)k*sizeof(ulong));
            for (slong i = 0; i < k; i++)
                memcpy(nmod_mat_entry_ptr(bm,i,0),bv+(size_t)t*nb+(size_t)i*n,(size_t)n*sizeof(ulong));
            nmod_mat_mul(cm,am,bm);
            for (slong i = 0; i < m; i++)
                memcpy(cv+(size_t)t*nc+(size_t)i*n,nmod_mat_entry_ptr(cm,i,0),(size_t)n*sizeof(ulong));
        }
        /* The barriers above finish all input reads before any output is
         * changed, including when C aliases A/B or a window of either. */
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (slong entry = 0; entry < (slong)nc; entry++) {
            for (slong t = 0; t < len; t++) tmp[t] = cv[(size_t)t*nc+entry];
            drsolve_ntt(tmp,len,inverse,mod);
            nmod_poly_struct *f = nmod_poly_mat_entry(C,entry/n,entry%n);
            nmod_poly_fit_length(f,outlen);
            for (slong t = 0; t < outlen; t++) f->coeffs[t] = nmod_mul(tmp[start+t],scale,mod);
            f->length = outlen; _nmod_poly_normalise(f);
        }
        nmod_mat_clear(cm); nmod_mat_clear(bm); nmod_mat_clear(am);
        flint_free(tmp);
    }
    flint_free(bv); flint_free(av); flint_free(cv);
    flint_free(inverse); flint_free(roots);
    return 1;
}

/* legacy: original full products; flint: bounded inputs without NTT;
 * ntt: force a transform attempt; auto: conservative shape/field dispatch. */
static int drsolve_mul_mode(void)
{
    const char *mode = getenv("DRSOLVE_PML_MUL");
    if (mode && strcmp(mode,"legacy") == 0) return -1;
    if (mode && strcmp(mode,"flint") == 0) return 0;
    if (mode && strcmp(mode,"ntt") == 0) return 1;
    return 2;
}

static int drsolve_try_ntt(const nmod_poly_mat_t A, const nmod_poly_mat_t B, int mode)
{
    if (mode == 1) return 1;
    /* Dispatch depends on shape; the backend checks root availability and
     * workspace uniformly for every prime. */
    return mode == 2 &&
        FLINT_MIN(A->r,FLINT_MIN(A->c,B->c)) >= 2;
}

void drsolve_nmod_poly_mat_mul(nmod_poly_mat_t C,
    const nmod_poly_mat_t A, const nmod_poly_mat_t B)
{
    if (drsolve_try_ntt(A,B,drsolve_mul_mode())) {
        slong a = nmod_poly_mat_max_length(A), b = nmod_poly_mat_max_length(B);
        if (a && b && a <= WORD_MAX-b &&
            drsolve_nmod_poly_mat_mul_window_ntt(C,A,B,0,a+b-1)) return;
    }
    nmod_poly_mat_mul(C,A,B);
}

void drsolve_nmod_poly_mat_middle_product(nmod_poly_mat_t C,
    const nmod_poly_mat_t A, const nmod_poly_mat_t B, slong start, slong count)
{
    int mode = drsolve_mul_mode();
    if (drsolve_try_ntt(A,B,mode) &&
        drsolve_nmod_poly_mat_mul_window_ntt(C,A,B,start,count)) return;
    /* Coefficients >= start+count cannot contribute. In the first PMBasis
     * branch B can still be the much higher-degree input of an ancestor. */
    nmod_poly_mat_t at,bt;
    const nmod_poly_mat_struct *a=A, *b=B;
    if (mode >= 0 && start >= 0 && count > 0 && start <= WORD_MAX-count) {
        slong hi = start+count;
        if (nmod_poly_mat_max_length(A) > hi) {
            nmod_poly_mat_init(at,A->r,A->c,A->modulus);
            nmod_poly_mat_set_trunc(at,A,hi); a=at;
        }
        if (nmod_poly_mat_max_length(B) > hi) {
            nmod_poly_mat_init(bt,B->r,B->c,B->modulus);
            nmod_poly_mat_set_trunc(bt,B,hi); b=bt;
        }
    }
    nmod_poly_mat_mul(C,a,b);
    nmod_poly_mat_shift_right(C,C,start);
    nmod_poly_mat_truncate(C,count);
    if (b != B) nmod_poly_mat_clear(bt);
    if (a != A) nmod_poly_mat_clear(at);
}

/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1), assuming deg(A) <= dA and deg(B) <= dA + dB
 *  output can alias input
 *  naive implementation (multiply, shift, truncate)
 */
void nmod_poly_mat_middle_product_naive(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                        const ulong dA, const ulong dB)
{
    nmod_poly_mat_mul(C, A, B);
    nmod_poly_mat_shift_right(C, C, dA);
    nmod_poly_mat_truncate(C, dB + 1);
}


/** Middle product for polynomial matrices
 *  sets C = ((A * B) div x^dA) mod x^(dB+1)
 *  output can alias input
 *  ASSUME: deg(A) <= dA and deg(B) <= dA + dB
 *  ASSUME: existence of primitive root
 *  uses geometric evaluation and interpolation
 */
void nmod_poly_mat_middle_product_geometric(nmod_poly_mat_t C, const nmod_poly_mat_t A, const nmod_poly_mat_t B,
                                            const ulong dA, const ulong dB)
{
    if (dA > WORD_MAX || dB >= WORD_MAX || dA > WORD_MAX-dB-1 ||
        nmod_poly_mat_max_length(A) > (slong)dA+1 ||
        nmod_poly_mat_max_length(B) > (slong)(dA+dB+1)) {
        nmod_poly_mat_middle_product_naive(C,A,B,dA,dB);
        return;
    }
    nmod_mat_t *mod_A, *mod_B, *mod_C;
    ulong ellC, order;
    ulong i, j, ell, m, k, n, u;
    long v;
    ulong p, w;
    nn_ptr val, val2;
    nmod_t mod;
    nmod_geometric_progression_t F;
    nmod_poly_t tmp_poly;

    m = A->r;
    k = A->c;
    n = B->c;
    p = A->modulus;

    if (m < 1 || n < 1 || k < 1)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    if (C == A || C == B)
    {
        nmod_poly_mat_t T;
        nmod_poly_mat_init(T, m, n, p);
        nmod_poly_mat_middle_product_geometric(T, A, B, dA, dB);
        nmod_poly_mat_swap_entrywise(C, T);
        nmod_poly_mat_clear(T);
        return;
    }

    // length = 0 iff matrix is zero
    if (nmod_poly_mat_max_length(A) == 0 || nmod_poly_mat_max_length(B) == 0)
    {
        nmod_poly_mat_zero(C);
        return;
    }

    nmod_init(&mod, p);

    ellC = dA + dB + 1;  // length(C) = length(A) + length(B) - 1
    order = ellC;
    nmod_init(&mod, p);
    w = nmod_find_root(2*order, mod);  /* TODO check necessary order */
    nmod_geometric_progression_init(F, w, order, mod);

    mod_A = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_B = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    mod_C = FLINT_ARRAY_ALLOC(ellC, nmod_mat_t);
    val = _nmod_vec_init(ellC);
    val2 = _nmod_vec_init(ellC);
    nmod_poly_init2(tmp_poly, mod.n, ellC);

#ifdef DIRTY_ALLOC_MATRIX
    // we alloc the memory for all matrices at once
    nn_ptr tmp = (nn_ptr) flint_malloc((m*k + k*n + m*n) * ellC * sizeof(ulong));
    nn_ptr bak;

    bak = tmp;
    for (i = 0; i < ellC; i++)
    {
        mod_A[i]->entries = tmp + i*m*k;
        mod_A[i]->stride = k;
        mod_A[i]->r = m;
        mod_A[i]->c = k;
        mod_A[i]->mod.n = mod.n;
        mod_A[i]->mod.norm = mod.norm;
        mod_A[i]->mod.ninv = mod.ninv;
    }
    tmp += ellC*m*k;

    for (i = 0; i < ellC; i++)
    {
        mod_B[i]->entries = tmp + i*k*n;
        mod_B[i]->stride = n;
        mod_B[i]->r = k;
        mod_B[i]->c = n;
        mod_B[i]->mod.n = mod.n;
        mod_B[i]->mod.norm = mod.norm;
        mod_B[i]->mod.ninv = mod.ninv;
    }
    tmp += ellC*k*n;

    for (i = 0; i < ellC; i++)
    {
        mod_C[i]->entries = tmp + i*m*n;
        mod_C[i]->stride = n;
        mod_C[i]->r = m;
        mod_C[i]->c = n;
        mod_C[i]->mod.n = mod.n;
        mod_C[i]->mod.norm = mod.norm;
        mod_C[i]->mod.ninv = mod.ninv;
    }
    tmp = bak;
#else
    for (i = 0; i < ellC; i++)
    {
        nmod_mat_init(mod_A[i], m, k, p);
        nmod_mat_init(mod_B[i], k, n, p);
        nmod_mat_init(mod_C[i], m, n, p);
    }
#endif

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < k; j++)
        {
            ulong len = nmod_poly_mat_entry(A, i, j)->length;
            nn_ptr src = nmod_poly_mat_entry(A, i, j)->coeffs;
            nn_ptr dest = tmp_poly->coeffs;
            v = dA;
            for (u = 0; u < len; u++, v--)
                dest[v] = src[u];
            for (; v >= 0; v--)
                dest[v] = 0;
            tmp_poly->length = dA + 1;
            _nmod_poly_normalise(tmp_poly);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val, tmp_poly->coeffs, tmp_poly->length, F, ellC, F->mod);
            for (ell = 0; ell < ellC; ell++)
                nmod_mat_entry(mod_A[ell], i, j) = val[ell];
        }
    }


    for (i = 0; i < k; i++)
    {
        for (j = 0; j < n; j++)
        {
            ulong len = nmod_poly_mat_entry(B, i, j)->length;
            nn_ptr src = nmod_poly_mat_entry(B, i, j)->coeffs;
            for (u = 0; u < len; u++)
                val2[u] = src[u];
            for (; u < ellC; u++)
                val2[u] = 0;
            nmod_poly_interpolate_geometric_nmod_vec_fast_precomp(tmp_poly, val2, F, ellC);
            len = tmp_poly->length;
            src = tmp_poly->coeffs;
            for (ell = 0; ell < len; ell++)
                nmod_mat_entry(mod_B[ell], i, j) = src[ell];
            for (; ell < ellC; ell++)
                nmod_mat_entry(mod_B[ell], i, j) = 0;
        }
    }

    for (ell = 0; ell < ellC; ell++)
        nmod_mat_mul(mod_C[ell], mod_A[ell], mod_B[ell]);

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            nn_ptr dest = tmp_poly->coeffs;
            for (ell = 0; ell < ellC; ell++)
                dest[ell] = nmod_mat_entry(mod_C[ell], i, j);
            tmp_poly->length = ellC;
            _nmod_poly_normalise(tmp_poly);

            _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(val2, tmp_poly->coeffs, tmp_poly->length, F, ellC, F->mod);

            nmod_poly_realloc(nmod_poly_mat_entry(C, i, j), dB + 1);
            nmod_poly_mat_entry(C, i, j)->length = dB + 1;
            dest = nmod_poly_mat_entry(C, i, j)->coeffs;
            for (u = 0; u <= dB; u++)
                dest[u] = val2[u];
            _nmod_poly_normalise(nmod_poly_mat_entry(C, i, j));
        }
    }

#ifdef DIRTY_ALLOC_MATRIX
    flint_free(tmp);
#else
    for (i = 0; i < ellC; i++)
    {
        nmod_mat_clear(mod_A[i]);
        nmod_mat_clear(mod_B[i]);
        nmod_mat_clear(mod_C[i]);
    }
#endif

    flint_free(mod_A);
    flint_free(mod_B);
    flint_free(mod_C);
    nmod_poly_clear(tmp_poly);
    _nmod_vec_clear(val2);
    _nmod_vec_clear(val);
    nmod_geometric_progression_clear(F);
}

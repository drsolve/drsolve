/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_GF28_POLY_MAT_DET_H
#define DRSOLVE_GF28_POLY_MAT_DET_H

/* Fused dst += scalar*x^shift*src. Distinct rows make this nonaliasing;
   one table lookup and XOR per coefficient, without field_elem_u temporaries. */
static void byte_poly_addmul_shift(gf28_poly_struct *dst,
    const gf28_poly_struct *src, slong shift, const uint8_t *scalar_row)
{
    if (!src->length) return;
    slong len = FLINT_MAX(dst->length, src->length + shift);
    gf28_poly_fit_length(dst, len);
    if (len > dst->length) memset(dst->coeffs + dst->length, 0, len - dst->length);
    for (slong i = 0; i < src->length; ++i)
        dst->coeffs[i + shift] ^= scalar_row[src->coeffs[i]];
    dst->length = len;
    gf28_poly_normalise(dst);
}

/* The same iterative weak-Popov elimination as the unified backend, with
   byte polynomial storage. Row additions and permutations preserve the
   determinant in characteristic two. The caller's matrix is read-only. */
static void byte_poly_mat_det(fq_nmod_poly_t det, const fq_nmod_poly_mat_t mat,
                              const fq_nmod_ctx_t ctx, int is_gf24)
{
    slong n = mat->r;
    gf28_poly_struct *entries = flint_malloc(n*n * sizeof(*entries));
    gf28_poly_struct **rows = flint_malloc(n * sizeof(*rows));
    gf28_poly_struct **ordered = flint_malloc(n * sizeof(*ordered));
    slong *pivots = flint_malloc(n * sizeof(*pivots));
    fq_nmod_t coeff;
    fq_nmod_init(coeff, ctx);
    for (slong i = 0; i < n; ++i) {
        rows[i] = entries + i*n;
        for (slong j = 0; j < n; ++j) {
            gf28_poly_struct *p = rows[i] + j;
            const fq_nmod_poly_struct *src = fq_nmod_poly_mat_entry(mat, i, j);
            gf28_poly_init(p);
            gf28_poly_fit_length(p, src->length);
            p->length = src->length;
            for (slong k = 0; k < src->length; ++k) {
                fq_nmod_poly_get_coeff(coeff, src, k, ctx);
                p->coeffs[k] = is_gf24 ? fq_nmod_to_gf24_elem(coeff, ctx)
                                         : fq_nmod_to_gf28_elem(coeff, ctx);
            }
        }
    }
    for (slong d = n-1; d >= 1; --d) {
        for (slong j = 0; j < d; ++j) pivots[j] = -1;
        slong rank = 0, zeros = 0;
        while (rank + zeros < d+1 && zeros < 2) {
            slong pivot = -1, length = 0;
            for (slong j = 0; j < d; ++j) {
                if (rows[rank][j].length > length) {
                    length = rows[rank][j].length; pivot = j;
                }
            }
            if (pivot < 0) {
                gf28_poly_struct *zero = rows[rank];
                for (slong i = rank; i < d; ++i) rows[i] = rows[i+1];
                rows[d] = zero; ++zeros;
            } else if (pivots[pivot] < 0) {
                pivots[pivot] = rank++;
            } else {
                slong owner = pivots[pivot];
                if (length < rows[owner][pivot].length) {
                    gf28_poly_struct *tmp = rows[owner];
                    rows[owner] = rows[rank]; rows[rank] = tmp;
                }
                const gf28_poly_struct *a = rows[rank] + pivot;
                const gf28_poly_struct *b = rows[owner] + pivot;
                slong shift = a->length - b->length;
                uint8_t lc = a->coeffs[a->length-1], divisor = b->coeffs[b->length-1];
                uint8_t scalar = is_gf24 ? gf24_mul(lc, gf24_inv(divisor))
                                          : gf28_mul(lc, gf28_inv(divisor));
                const uint8_t *table = is_gf24 ? gf24_get_scalar_row(scalar)
                                               : gf28_get_scalar_row(scalar);
                for (slong j = 0; j <= d; ++j)
                    byte_poly_addmul_shift(rows[rank] + j, rows[owner] + j, shift, table);
            }
        }
        if (rank < d || !rows[d][d].length) {
            fq_nmod_poly_zero(det, ctx);
            goto cleanup;
        }
        for (slong j = 0; j < d; ++j) ordered[j] = rows[pivots[j]];
        memcpy(rows, ordered, d * sizeof(*rows));
    }
    {
        fq_nmod_poly_t diagonal;
        fq_nmod_poly_init(diagonal, ctx);
        fq_nmod_poly_one(det, ctx);
        for (slong i = 0; i < n; ++i) {
            const gf28_poly_struct *p = rows[i] + i;
            fq_nmod_poly_zero(diagonal, ctx);
            for (slong k = 0; k < p->length; ++k) {
                if (is_gf24) gf24_elem_to_fq_nmod(coeff, p->coeffs[k], ctx);
                else gf28_elem_to_fq_nmod(coeff, p->coeffs[k], ctx);
                fq_nmod_poly_set_coeff(diagonal, k, coeff, ctx);
            }
            fq_nmod_poly_mul(det, det, diagonal, ctx);
        }
        fq_nmod_poly_clear(diagonal, ctx);
    }
cleanup:
    fq_nmod_clear(coeff, ctx);
    for (slong i = 0; i < n*n; ++i) gf28_poly_clear(entries+i);
    flint_free(pivots); flint_free(ordered); flint_free(rows); flint_free(entries);
}
#endif

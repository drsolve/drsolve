/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Shared exponent plans, separate native characteristic-two arithmetic. */
#define MQGF_TYPE uint8_t
#define MQGF_NAME dixon_schur_gf24
#define MQGF_MEMBER gf24
#define MQGF_ZERO(a) ((a) == 0)
#define MQGF_ADD(a,b) ((a) ^ (b))
#define MQGF_PREPARE(a) const uint8_t *scalar_row = gf24_get_scalar_row(a)
#define MQGF_MUL(a,b) scalar_row[b]
#include "dixon_gf2n_schur_template.h"

#define MQGF_TYPE uint8_t
#define MQGF_NAME dixon_schur_gf28
#define MQGF_MEMBER gf28
#define MQGF_ZERO(a) ((a) == 0)
#define MQGF_ADD(a,b) ((a) ^ (b))
#define MQGF_PREPARE(a) const uint8_t *scalar_row = gf28_get_scalar_row(a)
#define MQGF_MUL(a,b) scalar_row[b]
#include "dixon_gf2n_schur_template.h"

#define MQGF_TYPE uint16_t
#define MQGF_NAME dixon_schur_gf216
#define MQGF_MEMBER gf216
#define MQGF_ZERO(a) ((a) == 0)
#define MQGF_ADD(a,b) ((a) ^ (b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf216_mul(a,b)
#include "dixon_gf2n_schur_template.h"

#define MQGF_TYPE gf232_t
#define MQGF_NAME dixon_schur_gf232
#define MQGF_MEMBER gf232
#define MQGF_ZERO(a) gf232_is_zero(&(a))
#define MQGF_ADD(a,b) gf232_add(&(a), &(b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf232_mul(&(a), &(b))
#include "dixon_gf2n_schur_template.h"

#define MQGF_TYPE gf264_t
#define MQGF_NAME dixon_schur_gf264
#define MQGF_MEMBER gf264
#define MQGF_ZERO(a) gf264_is_zero(&(a))
#define MQGF_ADD(a,b) gf264_add(&(a), &(b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf264_mul(&(a), &(b))
#include "dixon_gf2n_schur_template.h"

#define MQGF_TYPE gf2128_t
#define MQGF_NAME dixon_schur_gf2128
#define MQGF_MEMBER gf2128
#define MQGF_ZERO(a) gf2128_is_zero(&(a))
#define MQGF_ADD(a,b) gf2128_add(&(a), &(b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf2128_mul(&(a), &(b))
#include "dixon_gf2n_schur_template.h"

/* Validate the weighted degree certificate before allocating native storage.
 * Sort only the complement; the core keeps the selected row/column order. */
#include "dixon_fq_schur.h"

static int dixon_extension_schur(fq_nmod_poly_mat_t core, fq_nmod_t factor,
    const fq_nmod_poly_mat_t matrix, const dixon_mq_step4_profile *p,
    const fq_nmod_ctx_t ctx)
{
    slong n = p->size, h = p->h, budget = 0;
    if (!n || matrix->r != n || matrix->c != n || h <= 0 || h >= n || p->sigma < 0) return 0;
    for (slong i = 0; i < n; i++) {
        if (p->rd[i] < 0 || p->cd[i] < 0 || p->rd[i] > p->sigma || p->cd[i] > p->sigma) return 0;
        if (i >= h) budget += p->sigma-p->rd[i]-p->cd[i];
        for (slong j = 0; j < n; j++) {
            const fq_nmod_poly_struct *entry = fq_nmod_poly_mat_entry(matrix, p->rows[i], p->cols[j]);
            if (entry->length && entry->length-1 > p->sigma-p->rd[i]-p->cd[j]) return 0;
        }
    }
    if (budget) return 0;
    slong *rows = flint_malloc(n*sizeof(slong)), *cols = flint_malloc(n*sizeof(slong));
    for (slong i = 0; i < n; i++) rows[i] = cols[i] = i;
    for (slong i = h+1; i < n; i++) {
        slong j = i, r = rows[i], c = cols[i];
        while (j > h && p->rd[rows[j-1]] > p->rd[r]) { rows[j] = rows[j-1]; j--; }
        rows[j] = r; j = i;
        while (j > h && p->cd[cols[j-1]] < p->cd[c]) { cols[j] = cols[j-1]; j--; }
        cols[j] = c;
    }
    int ok = 1;
    for (slong i = h; i < n; i++) if (p->rd[rows[i]]+p->cd[cols[i]] != p->sigma) ok = 0;
    field_ctx_t field; field_ctx_init(&field, ctx);
    if (ok) switch (field.field_id) {
    case FIELD_ID_GF24: ok = dixon_schur_gf24(core, factor, matrix, p, rows, cols, &field); break;
    case FIELD_ID_GF28: ok = dixon_schur_gf28(core, factor, matrix, p, rows, cols, &field); break;
    case FIELD_ID_GF216: ok = dixon_schur_gf216(core, factor, matrix, p, rows, cols, &field); break;
    case FIELD_ID_GF232: ok = dixon_schur_gf232(core, factor, matrix, p, rows, cols, &field); break;
    case FIELD_ID_GF264: ok = dixon_schur_gf264(core, factor, matrix, p, rows, cols, &field); break;
    case FIELD_ID_GF2128: ok = dixon_schur_gf2128(core, factor, matrix, p, rows, cols, &field); break;
    default: ok = dixon_schur_fq(core, factor, matrix, p, rows, cols, ctx); break;
    }
    field_ctx_clear(&field); flint_free(rows); flint_free(cols);
    return ok;
}

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
    if (!n || matrix->r != n || matrix->c != n || h < 0 || h >= n || p->sigma < 0) return 0;
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

/* Exact degree-layer pivot selection also works when the predicted Hilbert
 * quotas fail. Negative degree bounds make the selected complement block
 * triangular, with nonsingular constant diagonal blocks certified by LU.
 * Keep original matrix indices so the existing native kernels need no copy. */
static int dixon_extension_schur_adaptive(fq_nmod_poly_mat_t core,
    fq_nmod_t factor, const fq_nmod_poly_mat_t matrix,
    const dixon_mq_step4_profile *p, const fq_nmod_ctx_t ctx)
{
    slong n = matrix->r, sigma = p->size ? p->sigma : p->adaptive_sigma;
    const slong *rd = p->size ? p->rd : p->adaptive_rd;
    const slong *cd = p->size ? p->cd : p->adaptive_cd;
    if (!n || n != matrix->c || !rd || !cd || sigma < 0) return 0;
    for (slong i = 0; i < n; i++) {
        if (rd[i] < 0 || cd[i] < 0 || rd[i] > sigma || cd[i] > sigma) return 0;
        for (slong j = 0; j < n; j++) {
            const fq_nmod_poly_struct *a = fq_nmod_poly_mat_entry(matrix,
                p->size ? p->rows[i] : i, p->size ? p->cols[j] : j);
            if (a->length && a->length-1 > sigma-rd[i]-cd[j]) return 0;
        }
    }
    slong *ri = flint_malloc(n*sizeof(slong)), *ci = flint_malloc(n*sizeof(slong));
    slong *perm = flint_malloc(n*sizeof(slong));
    slong *rows = flint_malloc(n*sizeof(slong)), *cols = flint_malloc(n*sizeof(slong));
    unsigned char *usedr = flint_calloc(n,1), *usedc = flint_calloc(n,1);
    slong count = 0;
    for (slong degree = 0; degree <= sigma; degree++) {
        slong nr = 0, nc = 0;
        for (slong i = 0; i < n; i++) {
            if (rd[i] == degree) ri[nr++] = i;
            if (cd[i] == sigma-degree) ci[nc++] = i;
        }
        if (!nr || !nc) continue;
        fq_nmod_mat_t block;
        fq_nmod_mat_init(block,nr,nc,ctx);
        for (slong i = 0; i < nr; i++) for (slong j = 0; j < nc; j++)
            fq_nmod_poly_get_coeff(fq_nmod_mat_entry(block,i,j),
                fq_nmod_poly_mat_entry(matrix,p->size ? p->rows[ri[i]] : ri[i],
                    p->size ? p->cols[ci[j]] : ci[j]),0,ctx);
        slong rank = fq_nmod_mat_lu(perm,block,0,ctx), next = 0;
        for (slong i = 0; i < rank; i++) {
            while (next < nc && fq_nmod_is_zero(fq_nmod_mat_entry(block,i,next),ctx)) next++;
            rows[count] = ri[perm[i]]; cols[count] = ci[next++];
            usedr[rows[count]] = usedc[cols[count]] = 1; count++;
        }
        fq_nmod_mat_clear(block,ctx);
    }
    int ok = 0;
    if (count) {
        slong h = n-count;
        memmove(rows+h,rows,count*sizeof(slong));
        memmove(cols+h,cols,count*sizeof(slong));
        slong r = 0, c = 0;
        for (slong i = 0; i < n; i++) {
            if (!usedr[i]) rows[r++] = i;
            if (!usedc[i]) cols[c++] = i;
        }
        dixon_mq_step4_profile selected = {.size=n, .h=h, .sigma=sigma};
        selected.rows = flint_malloc(n*sizeof(slong));
        selected.cols = flint_malloc(n*sizeof(slong));
        selected.rd = flint_malloc(n*sizeof(slong));
        selected.cd = flint_malloc(n*sizeof(slong));
        for (slong i = 0; i < n; i++) {
            selected.rd[i] = rd[rows[i]]; selected.cd[i] = cd[cols[i]];
            selected.rows[i] = p->size ? p->rows[rows[i]] : rows[i];
            selected.cols[i] = p->size ? p->cols[cols[i]] : cols[i];
        }
        selected.odd = dixon_mq_step4_permutation_odd(selected.rows,n) ^
                       dixon_mq_step4_permutation_odd(selected.cols,n);
        fq_nmod_poly_mat_t candidate;
        fq_nmod_poly_mat_init(candidate,h,h,ctx);
        ok = dixon_extension_schur(candidate,factor,matrix,&selected,ctx);
        if (ok) {
            fq_nmod_poly_mat_struct swap = *core;
            *core = *candidate; *candidate = swap;
        }
        fq_nmod_poly_mat_clear(candidate,ctx);
        dixon_mq_step4_profile_clear(&selected);
    }
    flint_free(usedr); flint_free(usedc); flint_free(rows); flint_free(cols);
    flint_free(ri); flint_free(ci); flint_free(perm);
    return ok;
}

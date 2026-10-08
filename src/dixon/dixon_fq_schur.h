/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_fq_coeffs.h"

/* The caller certifies weighted degrees and pairs complement degree groups.
 * Store all coefficient vectors in one buffer, eliminating constant pivots
 * from the bottom degree group upwards. All permutation signs are retained. */
static int dixon_schur_fq(fq_nmod_poly_mat_t core, fq_nmod_t factor,
    const fq_nmod_poly_mat_t matrix, const dixon_mq_step4_profile *p,
    const slong *rows, const slong *cols, const fq_nmod_ctx_t ctx)
{
    slong n = p->size, h = p->h, d = fq_nmod_ctx_degree(ctx);
    if ((size_t)n > SIZE_MAX/(size_t)n/sizeof(ulong *)) return 0;
    size_t coefficients = 0, limit = SIZE_MAX/sizeof(ulong)/(size_t)d;
    for (slong i = 0; i < n; i++) for (slong j = 0; j < n; j++) {
        slong bound = p->sigma-p->rd[rows[i]]-p->cd[cols[j]];
        size_t capacity = bound < 0 ? 0 : (size_t)bound+1;
        if (capacity > limit-coefficients) return 0;
        coefficients += capacity;
    }
    /* Allocate only the certified degree capacity of each entry. Entries
     * with negative degree bounds consume no coefficient storage. */
    ulong *entries = flint_calloc(coefficients*d,sizeof(ulong));
    ulong **cells = flint_malloc((size_t)n*n*sizeof(*cells));
    ulong ***work = flint_malloc(n*sizeof(*work));
    slong *length = flint_calloc((size_t)n*n,sizeof(slong));
    slong **len = flint_malloc(n*sizeof(*len));
    size_t offset = 0;
    for (slong i = 0; i < n; i++) {
        work[i] = cells+(size_t)i*n; len[i] = length+(size_t)i*n;
        for (slong j = 0; j < n; j++) {
            work[i][j] = entries+offset*d;
            slong bound = p->sigma-p->rd[rows[i]]-p->cd[cols[j]];
            if (bound >= 0) offset += (size_t)bound+1;
            const fq_nmod_poly_struct *source = fq_nmod_poly_mat_entry(matrix,p->rows[rows[i]],p->cols[cols[j]]);
            len[i][j] = source->length;
            for (slong t = 0; t < source->length; t++)
                mq_fq_store(work[i][j]+t*d,source->coeffs+t,d);
        }
    }
    int odd = p->odd ^ dixon_mq_step4_permutation_odd(rows,n) ^ dixon_mq_step4_permutation_odd(cols,n);
    fq_nmod_t scale, scratch;
    fq_nmod_init(scale,ctx); fq_nmod_one(scale,ctx); fq_nmod_init(scratch,ctx);
    ulong *inverse = flint_malloc(d*sizeof(ulong)), *temp = flint_malloc(d*sizeof(ulong));
    int ok = 1;
    for (slong k = n; k-- > h;) {
        slong pivot = k;
        while (pivot >= h && p->rd[rows[pivot]] == p->rd[rows[k]] &&
               !mq_fq_length(work[pivot][k],d)) pivot--;
        if (pivot < h || p->rd[rows[pivot]] != p->rd[rows[k]]) { ok = 0; break; }
        if (pivot != k) {
            ulong **swap = work[pivot]; work[pivot] = work[k]; work[k] = swap;
            slong *lswap = len[pivot]; len[pivot] = len[k]; len[k] = lswap;
            odd ^= 1;
        }
        FLINT_ASSERT(len[k][k] == 1);
        fq_nmod_struct value = mq_fq_view(work[k][k],ctx);
        fq_nmod_mul(scale,scale,&value,ctx);
        fq_nmod_inv(scratch,&value,ctx); mq_fq_store(inverse,scratch,d);
        mq_fq_operator inv_op; mq_fq_operator_prepare(&inv_op,inverse,scratch,ctx);
        for (slong j = 0; j < k; j++) for (slong t = 0; t < len[k][j]; t++) {
            ulong *entry = work[k][j]+t*d;
            memset(temp,0,d*sizeof(ulong)); mq_fq_addmul(temp,&inv_op,entry,0,scratch,ctx);
            memcpy(entry,temp,d*sizeof(ulong));
        }
        #pragma omp parallel if(k > 128)
        {
            fq_nmod_t product; fq_nmod_init(product,ctx);
            #pragma omp for schedule(static)
            for (slong i = 0; i < k; i++) {
                slong alen = len[i][k];
                const ulong *a = work[i][k];
                for (slong ad = 0; ad < alen; ad++) {
                    if (!mq_fq_length(a+ad*d,d)) continue;
                    mq_fq_operator op; mq_fq_operator_prepare(&op,a+ad*d,product,ctx);
                    for (slong j = 0; j < k; j++) {
                        slong blen = len[k][j];
                        if (!blen) continue;
                        FLINT_ASSERT(alen+blen-2 <= p->sigma-p->rd[rows[i]]-p->cd[cols[j]]);
                        for (slong bd = 0; bd < blen; bd++)
                            mq_fq_addmul(work[i][j]+(ad+bd)*d,&op,
                                         work[k][j]+bd*d,1,product,ctx);
                        len[i][j] = FLINT_MAX(len[i][j],ad+blen);
                    }
                }
                for (slong j = 0; j < k; j++)
                    while (len[i][j] && !mq_fq_length(work[i][j]+(len[i][j]-1)*d,d)) len[i][j]--;
            }
            fq_nmod_clear(product,ctx);
        }
    }
    if (ok) {
        for (slong i = 0; i < h; i++) for (slong j = 0; j < h; j++) {
            fq_nmod_poly_struct *dest = fq_nmod_poly_mat_entry(core,i,j);
            fq_nmod_poly_zero(dest,ctx);
            for (slong t = 0; t < len[i][j]; t++) {
                fq_nmod_struct view = mq_fq_view(work[i][j]+t*d,ctx);
                fq_nmod_poly_set_coeff(dest,t,&view,ctx);
            }
        }
        if (odd) fq_nmod_neg(scale,scale,ctx);
        fq_nmod_set(factor,scale,ctx);
    }
    flint_free(temp); flint_free(inverse); fq_nmod_clear(scale,ctx); fq_nmod_clear(scratch,ctx);
    flint_free(len); flint_free(length); flint_free(work); flint_free(cells); flint_free(entries);
    return ok;
}

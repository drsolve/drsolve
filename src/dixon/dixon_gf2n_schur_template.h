/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Degree-graded Schur elimination over native GF(2^n). Every eliminated
 * pivot is a nonzero constant. Weighted degree bounds certify that all
 * updates fit the fixed polynomial stride; signs vanish in characteristic 2.
 * The input is immutable, including on singular-complement rejection. */
static int MQGF_NAME(fq_nmod_poly_mat_t core, fq_nmod_t factor,
    const fq_nmod_poly_mat_t matrix, const dixon_mq_step4_profile *p,
    const slong *rows, const slong *cols, const field_ctx_t *field)
{
    const fq_nmod_ctx_struct *ctx = field->ctx.fq_ctx;
    slong n = p->size, h = p->h, stride = p->sigma+1;
    if (stride <= 0 || (size_t)n > SIZE_MAX/(size_t)n/sizeof(MQGF_TYPE)/(size_t)stride) return 0;
    MQGF_TYPE *entries = flint_calloc((size_t)n*n*stride, sizeof(MQGF_TYPE));
    MQGF_TYPE **work = flint_malloc(n*sizeof(*work));
    slong *length = flint_calloc(n*n, sizeof(*length));
    slong **len = flint_malloc(n*sizeof(*len));
    field_elem_u value;
    field_set_one(&value, field->field_id, ctx);
    MQGF_TYPE scale = value.MQGF_MEMBER;
    for (slong i = 0; i < n; i++) {
        work[i] = entries+(size_t)i*n*stride; len[i] = length+i*n;
        for (slong j = 0; j < n; j++) {
            const fq_nmod_poly_struct *source = fq_nmod_poly_mat_entry(matrix, p->rows[rows[i]], p->cols[cols[j]]);
            len[i][j] = source->length;
            for (slong d = 0; d < source->length; d++) {
                fq_nmod_to_field_elem(&value, source->coeffs+d, field);
                work[i][j*stride+d] = value.MQGF_MEMBER;
            }
        }
    }
    int ok = 1;
    for (slong k = n; k-- > h;) {
        slong pivot = k;
        while (pivot >= h && p->rd[rows[pivot]] == p->rd[rows[k]] &&
               MQGF_ZERO(work[pivot][k*stride])) pivot--;
        if (pivot < h || p->rd[rows[pivot]] != p->rd[rows[k]]) { ok = 0; break; }
        if (pivot != k) {
            MQGF_TYPE *swap = work[pivot]; work[pivot] = work[k]; work[k] = swap;
            slong *lswap = len[pivot]; len[pivot] = len[k]; len[k] = lswap;
        }
        FLINT_ASSERT(len[k][k] == 1);
        value.MQGF_MEMBER = work[k][k*stride];
        field_elem_u inverse;
        field_inv(&inverse, &value, field->field_id, ctx);
        { MQGF_TYPE scalar = value.MQGF_MEMBER; MQGF_PREPARE(scalar);
          scale = MQGF_MUL(scalar, scale); }
        /* Normalize the pivot row once. Its complement entries preceding
         * the degree block vanish; its core entries can be polynomials. */
        { MQGF_TYPE scalar = inverse.MQGF_MEMBER; MQGF_PREPARE(scalar);
          for (slong j = 0; j < k; j++) for (slong d = 0; d < len[k][j]; d++)
              work[k][j*stride+d] = MQGF_MUL(scalar, work[k][j*stride+d]); }
        #pragma omp parallel for if(k > 128) schedule(static)
        for (slong i = 0; i < k; i++) {
            const MQGF_TYPE *a = work[i]+k*stride;
            slong alen = len[i][k];
            if (!alen) continue;
            for (slong j = 0; j < k; j++) {
                slong blen = len[k][j];
                if (!blen) continue;
                slong total = alen+blen-1;
                FLINT_ASSERT(total <= stride);
                MQGF_TYPE *dest = work[i]+j*stride;
                for (slong ad = 0; ad < alen; ad++) {
                    MQGF_TYPE scalar = a[ad];
                    if (MQGF_ZERO(scalar)) continue;
                    MQGF_PREPARE(scalar);
                    for (slong bd = 0; bd < blen; bd++) {
                        MQGF_TYPE product = MQGF_MUL(scalar, work[k][j*stride+bd]);
                        dest[ad+bd] = MQGF_ADD(dest[ad+bd], product);
                    }
                }
                total = FLINT_MAX(total, len[i][j]);
                while (total && MQGF_ZERO(dest[total-1])) total--;
                len[i][j] = total;
            }
        }
    }
    if (ok) {
        fq_nmod_t coefficient; fq_nmod_init(coefficient, ctx);
        for (slong i = 0; i < h; i++) for (slong j = 0; j < h; j++) {
            fq_nmod_poly_struct *dest = fq_nmod_poly_mat_entry(core, i, j);
            fq_nmod_poly_zero(dest, ctx);
            for (slong d = 0; d < len[i][j]; d++) {
                value.MQGF_MEMBER = work[i][j*stride+d];
                field_elem_to_fq_nmod(coefficient, &value, field);
                fq_nmod_poly_set_coeff(dest, d, coefficient, ctx);
            }
        }
        value.MQGF_MEMBER = scale; field_elem_to_fq_nmod(factor, &value, field);
        fq_nmod_clear(coefficient, ctx);
    }
    flint_free(len); flint_free(length); flint_free(work); flint_free(entries);
    return ok;
}
#undef MQGF_TYPE
#undef MQGF_NAME
#undef MQGF_MEMBER
#undef MQGF_ZERO
#undef MQGF_PREPARE
#undef MQGF_MUL
#undef MQGF_ADD

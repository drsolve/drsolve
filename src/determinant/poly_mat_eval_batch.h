/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_POLY_MAT_EVAL_BATCH_H
#define DRSOLVE_POLY_MAT_EVAL_BATCH_H
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>

/* Flatten entries into coefficient columns once. Every worker then evaluates a
 * bounded panel of points by [1, x, ...] * C. No roots of unity, prime-specific
 * arithmetic, or determinant-degree assumptions are needed for evaluation. */
typedef struct {
    nmod_mat_t coefficients;
    slong entries, length, batch, step;
    slong *indices;
    slong offsets[9];
} dixon_poly_mat_eval_plan;

static int dixon_poly_mat_eval_prepare(dixon_poly_mat_eval_plan *plan,
    const nmod_poly_mat_t mat, slong count, slong workers, int force)
{
    if (mat->r <= 0 || mat->c <= 0 || count <= 0 || workers <= 0 ||
        mat->r > WORD_MAX/mat->c) return 0;
    slong entries = mat->r*mat->c, length = nmod_poly_mat_max_length(mat);
    if (length <= 0 || (!force && (entries < 16 || count < 16 || length < 2))) return 0;
    const size_t coeff_limit = (size_t)64*1024*1024/sizeof(ulong);
    const size_t work_limit = (size_t)192*1024*1024/sizeof(ulong);
    if ((size_t)entries > coeff_limit/((size_t)length+1)) return 0;
    size_t per_worker = work_limit/(size_t)workers;
    /* Include an evaluated matrix, FLINT's determinant copy, and entry pointers. */
    if ((size_t)entries >= per_worker/3) return 0;
    size_t available = (per_worker-3*(size_t)entries)/((size_t)entries+(size_t)length);
    slong batch = FLINT_MIN(count,128);
    if (available < (size_t)batch) batch = (slong)available;
    if (batch < (force ? 1 : 8)) return 0;
    plan->entries=entries; plan->length=length; plan->batch=batch;
    /* At most eight degree groups. Columns are packed by their actual length,
     * avoiding multiplication by long zero tails in graded/sparse matrices. */
    plan->step=length/8+(length%8!=0);
    slong counts[8]={0}, next[8];
    for (slong i=0;i<mat->r;i++) for (slong j=0;j<mat->c;j++) {
        slong len=nmod_poly_mat_entry(mat,i,j)->length;
        counts[len ? (len-1)/plan->step : 0]++;
    }
    plan->offsets[0]=0;
    for (slong g=0;g<8;g++) {
        next[g]=plan->offsets[g];
        plan->offsets[g+1]=next[g]+counts[g];
    }
    nmod_mat_init(plan->coefficients,length,entries,mat->modulus);
    plan->indices=flint_malloc((size_t)entries*sizeof(slong));
    for (slong i=0;i<mat->r;i++) for (slong j=0;j<mat->c;j++) {
        const nmod_poly_struct *f=nmod_poly_mat_entry(mat,i,j);
        slong g=f->length ? (f->length-1)/plan->step : 0, row=next[g]++;
        plan->indices[row]=i*mat->c+j;
        for (slong d=0;d<f->length;d++)
            nmod_mat_entry(plan->coefficients,d,row)=f->coeffs[d];
    }
    return 1;
}

/* powers is batch x length and values is batch x entries. Each output row
 * is one point; column c represents flattened entry plan->indices[c]. Keeping
 * a whole evaluated matrix contiguous avoids a strided gather per point. */
static void dixon_poly_mat_eval_panel(nmod_mat_t values, nmod_mat_t powers,
    const dixon_poly_mat_eval_plan *plan, const ulong *xs, slong count)
{
    for (slong j=0;j<count;j++) {
        nmod_mat_entry(powers,j,0)=1;
        for (slong d=1;d<plan->length;d++)
            nmod_mat_entry(powers,j,d)=nmod_mul(nmod_mat_entry(powers,j,d-1),xs[j],powers->mod);
    }
    for (slong g=0;g<8;g++) {
        slong begin=plan->offsets[g], end=plan->offsets[g+1];
        if (begin==end) continue;
        slong length=FLINT_MIN(plan->length,(g+1)*plan->step);
        nmod_mat_t cv,pv,vv;
        nmod_mat_window_init(cv,plan->coefficients,0,begin,length,end);
        nmod_mat_window_init(pv,powers,0,0,count,length);
        nmod_mat_window_init(vv,values,0,begin,count,end);
        nmod_mat_mul(vv,pv,cv);
        nmod_mat_window_clear(vv); nmod_mat_window_clear(pv); nmod_mat_window_clear(cv);
    }
}

static void dixon_poly_mat_eval_clear(dixon_poly_mat_eval_plan *plan)
{
    flint_free(plan->indices);
    nmod_mat_clear(plan->coefficients);
}
#endif

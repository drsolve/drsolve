/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_mat_poly.h"

int unified_fq_poly_mat_det_hnf(fq_nmod_poly_t det, const fq_nmod_poly_mat_t a,
                               const fq_nmod_ctx_t ctx)
{
    if(a->r!=a->c) return 0;
    field_ctx_t field; field_ctx_init(&field,ctx);
    unified_mat_poly_t native;
    unified_mat_poly_init(native,a->r,a->c,&field);
    unified_mat_poly_from_fq(native,a,ctx);
    unified_poly_t result; unified_poly_init(result,&field);
    int ok=unified_mat_poly_det_hnf(result,native);
    if(ok) unified_to_fq_nmod_poly(det,result,ctx,&field);
    unified_poly_clear(result);unified_mat_poly_clear(native);field_ctx_clear(&field);
    return ok;
}

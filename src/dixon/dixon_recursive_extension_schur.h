/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Decode the recursive construction's raw rectangular indices before any
 * parameter-content extraction; these weights certify the original matrix. */
static int dixon_recursive_extension_schur(unified_mpoly_struct *result,
    unified_mpoly_struct **matrix, slong size, const unified_mpoly_struct *polys,
    slong nvars, const dixon_recursive_det_options *recursive)
{
    const fq_nmod_ctx_struct *ctx = polys[0].ctx;
    monom_t *rm = flint_malloc(size*sizeof(*rm)), *cm = flint_malloc(size*sizeof(*cm));
    slong *exps = flint_calloc(2*size*nvars,sizeof(slong));
    slong *indices = flint_malloc(size*sizeof(slong));
    for (slong i = 0; i < size; i++) {
        slong r = recursive->rows[i], c = recursive->cols[i];
        rm[i].exp = exps+2*i*nvars; cm[i].exp = rm[i].exp+nvars;
        rm[i].idx = cm[i].idx = indices[i] = i;
        for (slong v = nvars; v-- > 0;) {
            slong nr = (v+1)*recursive->variable_degrees[v];
            slong nc = (nvars-v)*recursive->variable_degrees[v];
            rm[i].exp[v] = nr ? r%nr : 0; if (nr) r /= nr;
            cm[i].exp[v] = nc ? c%nc : 0; if (nc) c /= nc;
        }
    }
    long *degrees = dixon_polynomial_degrees(polys,nvars+1,nvars);
    dixon_mq_step4_profile profile = {0};
    dixon_mq_step4_prepare(&profile,rm,cm,indices,indices,size,nvars,degrees);
    flint_free(degrees); flint_free(indices); flint_free(exps); flint_free(cm); flint_free(rm);
    fq_nmod_poly_mat_t a;
    fq_nmod_poly_mat_init(a,size,size,ctx);
    for (slong i = 0; i < size; i++) for (slong j = 0; j < size; j++) {
        const unified_mpoly_struct *source = &matrix[i][j];
        for (slong t = 0; t < dr_mpoly_length(source); t++) {
            DR_MPOLY_TERM(term,source,t);
            slong degree = term.par_exp ? term.par_exp[0] : 0;
            fq_nmod_poly_set_coeff(fq_nmod_poly_mat_entry(a,i,j),degree,term.coeff,ctx);
        }
    }
    fq_nmod_poly_t det; fq_nmod_poly_init(det,ctx);
    int ok = dixon_mq_step4_try(det,a,&profile,ctx);
    if (ok) {
        dr_mpoly_init(result,0,1,ctx);
        for (slong d = 0; d < det->length; d++)
            if (!fq_nmod_is_zero(det->coeffs+d,ctx))
                dr_mpoly_add_term_fast(result,NULL,&d,det->coeffs+d);
    }
    fq_nmod_poly_clear(det,ctx); fq_nmod_poly_mat_clear(a,ctx);
    dixon_mq_step4_profile_clear(&profile);
    return ok;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Certify the candidate in the actual extension field. In characteristic two,
 * the integer 2 is zero, so the third evaluation uses the field generator. */
static int dixon_gf2n_verify_projection(const unified_mpoly_struct *poly,
    slong nvars, slong rank, hash_entry_t **ri, slong rhs,
    hash_entry_t **ci, slong chs, const slong *rmap, const slong *cmap)
{
    field_ctx_t field; field_ctx_init(&field, poly->ctx);
    field_id_t id = field.field_id;
    const void *ctx = field.ctx.fq_ctx;
    field_elem_u *values = flint_calloc(rank*rank, sizeof(*values));
    field_elem_u powers[FLINT_BITS+2], point, coefficient, product, inverse;
    fq_nmod_t generator; fq_nmod_init(generator, poly->ctx);
    int ok = 0;
    for (int attempt = 0; attempt < 3 && !ok; attempt++) {
        memset(values, 0, rank*rank*sizeof(*values));
        if (attempt == 0) fq_nmod_one(generator, poly->ctx);
        else if (attempt == 1) fq_nmod_zero(generator, poly->ctx);
        else fq_nmod_gen(generator, poly->ctx);
        fq_nmod_to_field_elem(&point, generator, &field);
        field_set_one(powers, id, ctx);
        for (slong d = 1; d <= nvars+2; d++) field_mul(powers+d, powers+d-1, &point, id, ctx);
        for (slong t = 0; t < dr_mpoly_length(poly); t++) {
            DR_MPOLY_TERM(term, poly, t);
            slong r = lookup_monom_index(ri, rhs, term.var_exp, nvars);
            slong c = lookup_monom_index(ci, chs, term.var_exp+nvars, nvars);
            FLINT_ASSERT(r >= 0 && c >= 0 && rmap[r] >= 0 && cmap[c] >= 0);
            slong d = term.par_exp[0];
            fq_nmod_to_field_elem(&coefficient, term.coeff, &field);
            field_mul(&product, &coefficient, powers+d, id, ctx);
            field_elem_u *entry = values+rmap[r]*rank+cmap[c];
            field_add(entry, entry, &product, id, ctx);
        }
        ok = 1;
        for (slong col = 0; col < rank; col++) {
            slong pivot = col;
            while (pivot < rank && field_is_zero(values+pivot*rank+col, id, ctx)) pivot++;
            if (pivot == rank) { ok = 0; break; }
            if (pivot != col) for (slong j = col; j < rank; j++) {
                field_elem_u swap = values[col*rank+j];
                values[col*rank+j] = values[pivot*rank+j]; values[pivot*rank+j] = swap;
            }
            field_inv(&inverse, values+col*rank+col, id, ctx);
            for (slong j = col+1; j < rank; j++)
                field_mul(values+col*rank+j, values+col*rank+j, &inverse, id, ctx);
            #pragma omp parallel for if(rank-col > 128) schedule(static)
            for (slong r = col+1; r < rank; r++) {
                field_elem_u scalar = values[r*rank+col], prod;
                if (field_is_zero(&scalar, id, ctx)) continue;
                for (slong j = col+1; j < rank; j++) {
                    field_mul(&prod, &scalar, values+col*rank+j, id, ctx);
                    field_add(values+r*rank+j, values+r*rank+j, &prod, id, ctx);
                }
            }
        }
    }
    fq_nmod_clear(generator, poly->ctx);
    flint_free(values); field_ctx_clear(&field);
    return ok;
}

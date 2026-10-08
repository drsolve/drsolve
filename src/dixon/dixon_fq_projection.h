/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_fq_coeffs.h"

static int dixon_fq_verify_projection(const unified_mpoly_struct *poly,
    slong nvars, slong rank, hash_entry_t **ri, slong rhs,
    hash_entry_t **ci, slong chs, const slong *rmap, const slong *cmap)
{
    const fq_nmod_ctx_struct *ctx = poly->ctx;
    slong d = fq_nmod_ctx_degree(ctx);
    if (rank <= 0 || (size_t)rank > SIZE_MAX/(size_t)rank/(size_t)d/sizeof(ulong)) return 0;
    ulong *values = flint_calloc((size_t)rank*rank*d,sizeof(ulong));
    ulong *inverse = flint_malloc(d*sizeof(ulong)), *temp = flint_malloc(d*sizeof(ulong));
    fq_nmod_struct powers[FLINT_BITS+2];
    for (slong i = 0; i <= nvars+2; i++) fq_nmod_init(powers+i,ctx);
    fq_nmod_t point, product;
    fq_nmod_init(point,ctx); fq_nmod_init(product,ctx);
    int ok = 0;
    for (int attempt = 0; attempt < 3 && !ok; attempt++) {
        memset(values,0,(size_t)rank*rank*d*sizeof(ulong));
        if (attempt == 0) fq_nmod_one(point,ctx);
        else if (attempt == 1) fq_nmod_zero(point,ctx);
        else fq_nmod_gen(point,ctx);
        fq_nmod_one(powers,ctx);
        for (slong i = 1; i <= nvars+2; i++) fq_nmod_mul(powers+i,powers+i-1,point,ctx);
        for (slong t = 0; t < dr_mpoly_length(poly); t++) {
            DR_MPOLY_TERM(term,poly,t);
            slong r = lookup_monom_index(ri,rhs,term.var_exp,nvars);
            slong c = lookup_monom_index(ci,chs,term.var_exp+nvars,nvars);
            FLINT_ASSERT(r >= 0 && c >= 0 && rmap[r] >= 0 && cmap[c] >= 0);
            fq_nmod_mul(product,term.coeff,powers+term.par_exp[0],ctx);
            ulong *entry = values+((size_t)rmap[r]*rank+cmap[c])*d;
            for (slong j = 0; j < product->length; j++)
                entry[j] = nmod_add(entry[j],product->coeffs[j],fq_nmod_ctx_modulus(ctx)->mod);
        }
        ok = 1;
        for (slong col = 0; col < rank; col++) {
            slong pivot = col;
            while (pivot < rank && !mq_fq_length(values+((size_t)pivot*rank+col)*d,d)) pivot++;
            if (pivot == rank) { ok = 0; break; }
            if (pivot != col) for (size_t j = (size_t)col*d; j < (size_t)rank*d; j++) {
                ulong swap = values[(size_t)col*rank*d+j];
                values[(size_t)col*rank*d+j] = values[(size_t)pivot*rank*d+j];
                values[(size_t)pivot*rank*d+j] = swap;
            }
            fq_nmod_struct view = mq_fq_view(values+((size_t)col*rank+col)*d,ctx);
            fq_nmod_inv(product,&view,ctx); mq_fq_store(inverse,product,d);
            mq_fq_operator op; mq_fq_operator_prepare(&op,inverse,product,ctx);
            for (slong j = col+1; j < rank; j++) {
                ulong *entry = values+((size_t)col*rank+j)*d;
                memset(temp,0,d*sizeof(ulong)); mq_fq_addmul(temp,&op,entry,0,product,ctx);
                memcpy(entry,temp,d*sizeof(ulong));
            }
            #pragma omp parallel if(rank-col > 128)
            {
                fq_nmod_t scratch; fq_nmod_init(scratch,ctx);
                #pragma omp for schedule(static)
                for (slong r = col+1; r < rank; r++) {
                    const ulong *scalar = values+((size_t)r*rank+col)*d;
                    if (!mq_fq_length(scalar,d)) continue;
                    mq_fq_operator row_op; mq_fq_operator_prepare(&row_op,scalar,scratch,ctx);
                    for (slong j = col+1; j < rank; j++)
                        mq_fq_addmul(values+((size_t)r*rank+j)*d,&row_op,
                            values+((size_t)col*rank+j)*d,1,scratch,ctx);
                }
                fq_nmod_clear(scratch,ctx);
            }
        }
    }
    fq_nmod_clear(point,ctx); fq_nmod_clear(product,ctx);
    for (slong i = 0; i <= nvars+2; i++) fq_nmod_clear(powers+i,ctx);
    flint_free(values); flint_free(inverse); flint_free(temp);
    return ok;
}

/* A deficient specialization is only a lower bound for polynomial rank.
 * Before abandoning an extension-field candidate, retry distinct points in
 * the actual field. Preserve the original point if no certificate succeeds. */
static slong dixon_fq_candidate_rank(unified_mpoly_struct ***matrix,
    const slong *rows, const slong *cols, slong size, slong npars,
    fq_nmod_t *params, const fq_nmod_ctx_t ctx)
{
    fq_nmod_mat_t values; fq_nmod_mat_init(values,size,size,ctx);
    fq_nmod_t original; fq_nmod_init(original,ctx);
    if (npars == 1) fq_nmod_set(original,params[0],ctx);
    slong best = 0;
    for (int attempt = 0; attempt < (npars == 1 ? 4 : 1); attempt++) {
        if (attempt == 1) fq_nmod_one(params[0],ctx);
        else if (attempt == 2) fq_nmod_zero(params[0],ctx);
        else if (attempt == 3) fq_nmod_gen(params[0],ctx);
        if (attempt && fq_nmod_equal(params[0],original,ctx)) continue;
        for (slong i = 0; i < size; i++) for (slong j = 0; j < size; j++) {
            unified_mpoly_struct *entry = matrix[rows[i]][cols[j]];
            if (entry && dr_mpoly_length(entry))
                evaluate_dr_mpoly_at_params(fq_nmod_mat_entry(values,i,j),entry,params);
            else fq_nmod_zero(fq_nmod_mat_entry(values,i,j),ctx);
        }
        slong rank = fq_nmod_mat_rank(values,ctx);
        best = FLINT_MAX(best,rank);
        if (rank == size) break;
    }
    if (npars == 1 && best < size) fq_nmod_set(params[0],original,ctx);
    fq_nmod_clear(original,ctx); fq_nmod_mat_clear(values,ctx);
    return best;
}

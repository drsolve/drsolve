/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Product divided by each linear factor; valid in every characteristic. */
static fq_nmod_poly_struct *interpolation_weights(const dq_coeff *nodes,slong count,
                                                 const fq_nmod_ctx_t ctx)
{
    fq_nmod_poly_struct *weights=flint_malloc(count*sizeof(*weights));
    fq_nmod_poly_t product,linear,temp,remainder;
    fq_nmod_poly_init(product,ctx);fq_nmod_poly_init(linear,ctx);
    fq_nmod_poly_init(temp,ctx);fq_nmod_poly_init(remainder,ctx);
    fq_nmod_t node,denominator,one;
    fq_nmod_init(node,ctx);fq_nmod_init(denominator,ctx);fq_nmod_init(one,ctx);
    fq_nmod_one(one,ctx);fq_nmod_poly_one(product,ctx);
    fq_nmod_poly_set_coeff(linear,1,one,ctx);
    for(slong i=0;i<count;i++) {
        field_elem_to_fq_nmod(node,nodes+i,&extension_field);fq_nmod_neg(node,node,ctx);
        fq_nmod_poly_set_coeff(linear,0,node,ctx);
        fq_nmod_poly_mul(temp,product,linear,ctx);fq_nmod_poly_swap(temp,product,ctx);
    }
    for(slong i=0;i<count;i++) {
        fq_nmod_poly_init(weights+i,ctx);
        field_elem_to_fq_nmod(node,nodes+i,&extension_field);fq_nmod_neg(node,node,ctx);
        fq_nmod_poly_set_coeff(linear,0,node,ctx);
        fq_nmod_poly_divrem(weights+i,remainder,product,linear,ctx);
        fq_nmod_neg(node,node,ctx);
        fq_nmod_poly_evaluate_fq_nmod(denominator,weights+i,node,ctx);
        assert(!fq_nmod_is_zero(denominator,ctx));fq_nmod_inv(denominator,denominator,ctx);
        fq_nmod_poly_scalar_mul_fq_nmod(weights+i,weights+i,denominator,ctx);
    }
    fq_nmod_clear(node,ctx);fq_nmod_clear(denominator,ctx);fq_nmod_clear(one,ctx);
    fq_nmod_poly_clear(product,ctx);fq_nmod_poly_clear(linear,ctx);
    fq_nmod_poly_clear(temp,ctx);fq_nmod_poly_clear(remainder,ctx);return weights;
}


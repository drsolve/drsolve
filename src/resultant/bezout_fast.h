/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Exact bivariate Dixon/Bezout resultant by scalar resultants and interpolation.
 * Extension adapters borrow the original field and its defining modulus.
 * No monic normalization here: unified_mpoly_resultant returns the exact scalar.
 * Included only by unified_mpoly_resultant.c. */
#include <flint/gr.h>
#include <flint/gr_poly.h>
#include <flint/gr_vec.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "dixon_flint.h"
#include "zech_mpoly_access.h"

static void bezout_get_coeff(gr_ptr c, const unified_mpoly_t a, slong i)
{
    unified_mpoly_ctx_t ctx=a->ctx_ptr;
    switch(a->field_id) {
    case FIELD_ID_NMOD:
        *(ulong *)c=nmod_mpoly_get_term_coeff_ui(GET_NMOD_POLY(a),i,GET_NMOD_CTX(ctx));break;
    case FIELD_ID_FMPZ_MOD:
        fmpz_mod_mpoly_get_term_coeff_fmpz(c,GET_BIG_POLY(a),i,GET_BIG_CTX(ctx));break;
    case FIELD_ID_FQ_ZECH:
        fq_zech_mpoly_get_term_coeff_fq_zech(c,GET_ZECH_POLY(a),i,GET_ZECH_CTX(ctx));break;
    default:
        fq_nmod_mpoly_get_term_coeff_fq_nmod(c,GET_FQ_POLY(a),i,GET_FQ_CTX(ctx));break;
    }
}

static void bezout_set_coeff(unified_mpoly_t a, gr_srcptr c, const ulong *exps)
{
    unified_mpoly_ctx_t ctx=a->ctx_ptr;
    switch(a->field_id) {
    case FIELD_ID_NMOD:
        nmod_mpoly_set_coeff_ui_ui(GET_NMOD_POLY(a),*(const ulong *)c,exps,GET_NMOD_CTX(ctx));break;
    case FIELD_ID_FMPZ_MOD:
        fmpz_mod_mpoly_set_coeff_fmpz_ui(GET_BIG_POLY(a),c,exps,GET_BIG_CTX(ctx));break;
    case FIELD_ID_FQ_ZECH:
        dr_zech_set_coeff(GET_ZECH_POLY(a),c,exps,GET_ZECH_CTX(ctx));break;
    default:
        fq_nmod_mpoly_set_coeff_fq_nmod_ui(GET_FQ_POLY(a),c,exps,GET_FQ_CTX(ctx));break;
    }
}

/* Enumerate distinct extension-field points using base-p digits, not integer
 * coercion (which would repeat after p points). index is known to be < |K|. */
static int bezout_point(gr_ptr out, ulong index, field_ctx_t *field, gr_ctx_t ctx)
{
    if(field->field_id==FIELD_ID_NMOD || field->field_id==FIELD_ID_FMPZ_MOD)
        return gr_set_ui(out,index,ctx);
    const fq_nmod_ctx_struct *fq=field->field_id==FIELD_ID_FQ_ZECH
        ? field->ctx.zech_ctx->fq_nmod_ctx : field->ctx.fq_ctx;
    ulong p=fq_nmod_ctx_prime(fq);
    fq_nmod_t value;fq_nmod_init(value,fq);
    for(slong j=0;index;j++,index/=p) nmod_poly_set_coeff_ui(value,j,index%p);
    if(field->field_id==FIELD_ID_FQ_ZECH)
        fq_zech_set_fq_nmod(out,value,field->ctx.zech_ctx);
    else fq_nmod_set(out,value,fq);
    fq_nmod_clear(value,fq);
    return GR_SUCCESS;
}

static int bezout_ntt(gr_ptr a, slong n, gr_srcptr roots, gr_ctx_t ctx)
{
    slong sz=ctx->sizeof_elem;
    for(slong i=1,j=0;i<n;i++) {
        slong bit=n>>1;
        for(;j&bit;bit>>=1)j^=bit;
        j^=bit;
        if(i<j)gr_swap(GR_ENTRY(a,i,sz),GR_ENTRY(a,j,sz),ctx);
    }
    int status=GR_SUCCESS;
    gr_ptr u=gr_heap_init(ctx),v=gr_heap_init(ctx);
    for(slong len=2;len<=n;len<<=1)
        for(slong i=0;i<n;i+=len)
            for(slong j=0;j<len/2;j++) {
                gr_ptr x=GR_ENTRY(a,i+j,sz),y=GR_ENTRY(a,i+j+len/2,sz);
                status|=gr_set(u,x,ctx);
                status|=gr_mul(v,y,GR_ENTRY(roots,j*(n/len),sz),ctx);
                status|=gr_add(x,u,v,ctx);
                status|=gr_sub(y,u,v,ctx);
            }
    gr_heap_clear(u,ctx);gr_heap_clear(v,ctx);
    return status;
}

/* At the prime-field nodes 0,...,n-1, the inverse derivative of the
 * node polynomial is (-1)^(n-1-i)/(i! (n-1-i)!). n <= p guarantees
 * every factorial is invertible, including the boundary n == p.
 * One inversion replaces derivative evaluation and n separate inversions. */
static int bezout_grid_weights(gr_ptr weights, slong n, gr_ctx_t ctx)
{
    slong sz=ctx->sizeof_elem;
    gr_ptr inverse_factorials=gr_heap_init_vec(n,ctx);
    gr_ptr last=GR_ENTRY(inverse_factorials,n-1,sz);
    int status=gr_one(last,ctx);
    for(slong i=2;i<n;i++)status|=gr_mul_ui(last,last,i,ctx);
    status|=gr_inv(last,last,ctx);
    for(slong i=n-1;i>0;i--)
        status|=gr_mul_ui(GR_ENTRY(inverse_factorials,i-1,sz),
                          GR_ENTRY(inverse_factorials,i,sz),i,ctx);
    for(slong i=0;i<n;i++) {
        gr_ptr w=GR_ENTRY(weights,i,sz);
        status|=gr_mul(w,GR_ENTRY(inverse_factorials,i,sz),
                        GR_ENTRY(inverse_factorials,n-1-i,sz),ctx);
        if((n-1-i)&1)status|=gr_neg(w,w,ctx);
    }
    gr_heap_clear_vec(inverse_factorials,n,ctx);
    return status;
}

int unified_mpoly_resultant_bezout_try(unified_mpoly_t out,
    const unified_mpoly_t a, const unified_mpoly_t b, slong var, unified_mpoly_ctx_t ring)
{
    if(!g_dixon_bezout_fast || g_field_equation_reduction || ring->nvars!=2 ||
       var<0 || var>1 || a->ctx_ptr!=ring || b->ctx_ptr!=ring || out->ctx_ptr!=ring)
        return 0;
    unified_mpoly_t inputs[2]={a,b};
    slong dx[2]={0},total[2]={0};ulong exps[2];
    for(int k=0;k<2;k++) {
        slong degrees[2];unified_mpoly_degrees_si(degrees,inputs[k]);
        if(degrees[0]<0 || degrees[1]<0) return 0;
        dx[k]=degrees[var];total[k]=unified_mpoly_total_degree_si(inputs[k]);
        if(dx[k]<1 || total[k]<dx[k]) return 0;
        for(slong i=0;i<unified_mpoly_length(inputs[k]);i++) {
            unified_mpoly_get_term_exp_ui(exps,inputs[k],i,ring);
            if(exps[var]==(ulong)dx[k] && exps[1-var]) return 0;
        }
    }
    /* Only representation/overflow checks, no degree or memory budget. */
    if(dx[1]>(WORD_MAX-1)/total[0])return 0;
    slong bound=dx[1]*total[0];
    if(total[1]-dx[1]>(WORD_MAX-1-bound)/dx[0])return 0;
    bound+=dx[0]*(total[1]-dx[1]);
    if(dx[1]>WORD_MAX-2 || dx[0]>WORD_MAX-2-dx[1])return 0;
    slong n=1,rows=dx[0]+dx[1]+2;
    while(n<=bound) {
        if(n>WORD_MAX/2){n=0;break;} /* Multipoint may still fit. */
        n*=2;
    }
    field_ctx_t *field=ring->field_ctx;
    gr_ctx_t ctx;fmpz_t order,qm1,exponent;
    fmpz_init(order);fmpz_init(qm1);fmpz_init(exponent);
    if(field->field_id==FIELD_ID_NMOD) {
        _gr_ctx_init_nmod(ctx,&field->ctx.nmod_ctx);
        fmpz_set_ui(order,field->ctx.nmod_ctx.n);
    } else if(field->field_id==FIELD_ID_FMPZ_MOD) {
        fmpz_set(order,fmpz_mod_ctx_modulus(field->ctx.big_ctx));
        gr_ctx_init_fmpz_mod(ctx,order);
    } else if(field->field_id==FIELD_ID_FQ_ZECH) {
        _gr_ctx_init_fq_zech_from_ref(ctx,field->ctx.zech_ctx);
        fq_zech_ctx_order(order,field->ctx.zech_ctx);
    } else {
        _gr_ctx_init_fq_nmod_from_ref(ctx,field->ctx.fq_ctx);
        fq_nmod_ctx_order(order,field->ctx.fq_ctx);
    }
    int ok=0,status=GR_SUCCESS,ntt=0;
    if(field->field_id==FIELD_ID_NMOD || field->field_id==FIELD_ID_FMPZ_MOD)
        status=gr_ctx_set_is_field(ctx,T_TRUE);
    if(status!=GR_SUCCESS) goto cleanup_ctx;
    if(fmpz_cmp_ui(order,bound+1)<0) goto cleanup_ctx;
    /* Validate both possible transform layouts before allocating scratch. */
    slong count=n ? n : bound+1;
    if(rows>WORD_MAX/count || (size_t)rows>SIZE_MAX/(size_t)ctx->sizeof_elem/(size_t)count ||
       (size_t)rows>SIZE_MAX/sizeof(gr_poly_struct)) goto cleanup_ctx;
    gr_ptr root=gr_heap_init(ctx),temp=gr_heap_init(ctx),candidate=gr_heap_init(ctx);
    fmpz_sub_ui(qm1,order,1);
    if(n && fmpz_fdiv_ui(qm1,n)==0) {
        fmpz_divexact_ui(exponent,qm1,n);
        for(ulong i=2;i<130;i++) {
            if(i<66 || field->field_id==FIELD_ID_NMOD || field->field_id==FIELD_ID_FMPZ_MOD)
                status=gr_set_ui(candidate,i,ctx);
            else {
                status=gr_gen(candidate,ctx);
                status|=gr_add_ui(candidate,candidate,i-66,ctx);
            }
            status|=gr_pow_fmpz(root,candidate,exponent,ctx);
            status|=gr_pow_ui(temp,root,n/2,ctx);
            if(status==GR_SUCCESS && gr_is_neg_one(temp,ctx)==T_TRUE) {ntt=1;break;}
        }
    }
    status=GR_SUCCESS;
    if(!ntt)n=bound+1;
    double start=get_wall_time();
    gr_vec_t points,samples,values;
    gr_vec_init(points,n,ctx);gr_vec_init(samples,n,ctx);gr_vec_init(values,rows*n,ctx);
    gr_poly_struct *coeffs=flint_malloc(rows*sizeof(*coeffs));
    for(slong i=0;i<rows;i++)gr_poly_init(coeffs+i,ctx);
    slong sz=ctx->sizeof_elem;
    for(int k=0;k<2;k++) for(slong i=0;i<unified_mpoly_length(inputs[k]);i++) {
        unified_mpoly_get_term_exp_ui(exps,inputs[k],i,ring);
        bezout_get_coeff(temp,inputs[k],i);
        slong row=exps[var]+(k ? dx[0]+1 : 0);
        status|=gr_poly_set_coeff_scalar(coeffs+row,exps[1-var],temp,ctx);
    }
    gr_ptr *tree=NULL;
    if(ntt) {
        status|=gr_one(points->entries,ctx);
        for(slong i=1;i<n;i++)status|=gr_mul(GR_ENTRY(points->entries,i,sz),
            GR_ENTRY(points->entries,i-1,sz),root,ctx);
    } else {
        for(slong i=0;i<n;i++)status|=bezout_point(GR_ENTRY(points->entries,i,sz),i,field,ctx);
        tree=_gr_poly_tree_alloc(n,ctx);
        status|=_gr_poly_tree_build(tree,points->entries,n,ctx);
    }
    if(status==GR_SUCCESS) {
#ifdef _OPENMP
        #pragma omp parallel for schedule(static) reduction(|:status) if(!omp_in_parallel() && rows>=8)
#endif
        for(slong i=0;i<rows;i++) {
            gr_ptr target=GR_ENTRY(values->entries,i*n,sz);
            if(ntt) {
                for(slong j=0;j<n;j++)status|=gr_zero(GR_ENTRY(target,j,sz),ctx);
                for(slong j=0;j<coeffs[i].length;j++)
                    status|=gr_set(GR_ENTRY(target,j,sz),GR_ENTRY(coeffs[i].coeffs,j,sz),ctx);
                status|=bezout_ntt(target,n,points->entries,ctx);
            } else status|=_gr_poly_evaluate_vec_fast_precomp(target,coeffs[i].coeffs,
                coeffs[i].length,tree,n,ctx);
        }
    }
    double evaluated=get_wall_time();
    if(status==GR_SUCCESS) {
#ifdef _OPENMP
        #pragma omp parallel reduction(|:status) if(!omp_in_parallel() && n>=32)
#endif
        {
            gr_poly_t f,g;gr_poly_init(f,ctx);gr_poly_init(g,ctx);
            gr_poly_fit_length(f,dx[0]+1,ctx);gr_poly_fit_length(g,dx[1]+1,ctx);
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            for(slong i=0;i<n;i++) {
                for(slong j=0;j<=dx[0];j++)status|=gr_set(GR_ENTRY(f->coeffs,j,sz),
                    GR_ENTRY(values->entries,j*n+i,sz),ctx);
                for(slong j=0;j<=dx[1];j++)status|=gr_set(GR_ENTRY(g->coeffs,j,sz),
                    GR_ENTRY(values->entries,(dx[0]+1+j)*n+i,sz),ctx);
                f->length=dx[0]+1;g->length=dx[1]+1;
                status|=gr_poly_resultant(GR_ENTRY(samples->entries,i,sz),f,g,ctx);
            }
            gr_poly_clear(f,ctx);gr_poly_clear(g,ctx);
        }
    }
    double computed=get_wall_time();
    gr_poly_t result;gr_poly_init(result,ctx);
    if(status==GR_SUCCESS && ntt) {
        /* Reverse the nonconstant root powers to obtain the inverse NTT. */
        for(slong i=1;i<n/2;i++)gr_swap(GR_ENTRY(points->entries,i,sz),GR_ENTRY(points->entries,n-i,sz),ctx);
        status|=bezout_ntt(samples->entries,n,points->entries,ctx);
        status|=gr_set_ui(temp,n,ctx);status|=gr_inv(temp,temp,ctx);
        for(slong i=0;i<n;i++) {
            gr_ptr c=GR_ENTRY(samples->entries,i,sz);
            status|=gr_mul(c,c,temp,ctx);
            if(i>bound && gr_is_zero(c,ctx)!=T_TRUE)status|=GR_UNABLE;
            status|=gr_poly_set_coeff_scalar(result,i,c,ctx);
        }
    } else if(status==GR_SUCCESS) {
        gr_ptr weights=gr_heap_init_vec(n,ctx);
        if(field->field_id==FIELD_ID_NMOD || field->field_id==FIELD_ID_FMPZ_MOD)
            status|=bezout_grid_weights(weights,n,ctx);
        else status|=_gr_poly_interpolation_weights(weights,tree,n,ctx);
        gr_poly_fit_length(result,n,ctx);
        if(status==GR_SUCCESS) {
            status|=_gr_poly_interpolate_fast_precomp(result->coeffs,samples->entries,tree,weights,n,ctx);
            result->length=n;
            _gr_poly_normalise(result,ctx);
        }
        gr_heap_clear_vec(weights,n,ctx);
    }
    if(status==GR_SUCCESS) {
        /* Publish only after all computations succeeded; aliases are safe. */
        unified_mpoly_t answer=unified_mpoly_init(ring);
        exps[0]=exps[1]=0;
        /* Descending parameter powers append in canonical monomial order;
         * ascending insertion repeatedly shifts the already built output. */
        for(slong i=result->length;i-- >0;) {
            exps[1-var]=i;
            bezout_set_coeff(answer,GR_ENTRY(result->coeffs,i,sz),exps);
        }
        unified_mpoly_swap(out,answer);unified_mpoly_clear(answer);ok=1;
        dixon_info_log("Dixon/Bezout fast resultant: %s, %ld points, degree bound %ld\n",
                       ntt ? "NTT" : "multipoint",n,bound);
        dixon_info_log("  Evaluation %.6fs | scalar resultants %.6fs | interpolation/output %.6fs\n",
            evaluated-start,computed-evaluated,get_wall_time()-computed);
        dixon_info_log("Dixon/Bezout fast time: Wall time: %.6f seconds\n",get_wall_time()-start);
    }
    gr_poly_clear(result,ctx);
    if(tree)_gr_poly_tree_free(tree,n,ctx);
    for(slong i=0;i<rows;i++)gr_poly_clear(coeffs+i,ctx);
    flint_free(coeffs);gr_vec_clear(values,ctx);gr_vec_clear(points,ctx);gr_vec_clear(samples,ctx);
    gr_heap_clear(root,ctx);gr_heap_clear(temp,ctx);gr_heap_clear(candidate,ctx);
cleanup_ctx:
    /* FLINT's from_ref extension adapters allocate no context storage, but
     * gr_ctx_clear would destroy/free the borrowed field. Only clear owned
     * prime-field contexts here. All temporary elements were cleared above. */
    if(field->field_id==FIELD_ID_NMOD || field->field_id==FIELD_ID_FMPZ_MOD)gr_ctx_clear(ctx);
    fmpz_clear(order);fmpz_clear(qm1);fmpz_clear(exponent);
    return ok;
}

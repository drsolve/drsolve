/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Certified Macaulay closure over extension fields. Relations stay in the
 * original ideal. Commuting, cyclic multiplication operators annihilating
 * the input certify the full quotient, including nonreduced algebras. */
#include "quotient_solver.h"
#include "quotient_basis.h"
#include "fq_mat_det.h"
#include <flint/fq_nmod_mpoly.h>
#include <flint/fq_nmod_poly_factor.h>
#include <flint/fq_nmod_mpoly_factor.h>
#include "quotient_byte.h"

#include "quotient_closure.h"

static void extension_init(quotient_field_t *e, const fq_nmod_ctx_t fq)
{
    memset(e,0,sizeof(*e)); e->fq=fq; field_ctx_init(&e->field,fq);
    e->arithmetic=field_ctx_data(&e->field);
    if (e->field.field_id==FIELD_ID_GF2N8)
        memcpy(e->byte_tables,((const gf2n_native_ctx *)e->arithmetic)->byte_simd,sizeof(e->byte_tables));
    if (e->field.field_id==FIELD_ID_GF24)
        quotient_byte_tables(e->byte_tables,16,gf24_get_scalar_row);
    else if (e->field.field_id==FIELD_ID_GF28)
        quotient_byte_tables(e->byte_tables,256,gf28_get_scalar_row);
}

static slong extension_input_length(const quotient_input_t *input, slong i)
{
    return ((const fq_nmod_mpoly_struct *)input->polys)[i].length;
}

static void extension_input_term(field_elem_u *c, ulong *exp,
                                const quotient_input_t *input, slong i, slong t,
                                const quotient_field_t *e)
{
    const fq_nmod_mpoly_struct *f=(const fq_nmod_mpoly_struct *)input->polys+i;
    const fq_nmod_mpoly_ctx_struct *ctx=input->ctx;
    fq_nmod_t a; fq_nmod_init(a,e->fq);
    fq_nmod_mpoly_get_term_coeff_fq_nmod(a,f,t,ctx);
    fq_nmod_to_field_elem(c,a,&e->field);
    fq_nmod_mpoly_get_term_exp_ui(exp,f,t,ctx);
    fq_nmod_clear(a,e->fq);
}

/* Root extraction uses FLINT over the original fq context. Only the small
 * certified operators cross this boundary; closure stays in native fields.
 * Joint eigenspaces handle repeated roots and nonseparating coordinates. */
static int extension_points(polynomial_solutions_t *sols, fq_nmod_mat_struct *T,
              const fq_nmod_mat_t S, slong variable, fq_nmod_struct *point,
              fq_nmod_mpoly_struct *fs, slong m, const fq_nmod_mpoly_ctx_t ctx)
{
    const fq_nmod_ctx_struct *fq=sols->ctx;
    slong n=sols->num_variables,dim=S->r,k=S->c;
    if (variable==n) {
        if (k!=1 || sols->num_solution_sets>=T[0].r) return 0;
        fq_nmod_struct *values[n]; for (slong j=0;j<n;j++) values[j]=point+j;
        fq_nmod_t ev; fq_nmod_init(ev,fq);
        int ok=1;
        for (slong i=0;i<m && ok;i++) {
            fq_nmod_mpoly_evaluate_all_fq_nmod(ev,fs+i,values,ctx);
            ok=fq_nmod_is_zero(ev,fq);
        }
        fq_nmod_clear(ev,fq); if (!ok) return 0;
        slong row=sols->num_solution_sets++;
        sols->solution_sets[row]=calloc(n,sizeof(fq_nmod_t *));
        for (slong j=0;j<n;j++) {
            sols->solution_sets[row][j]=malloc(sizeof(fq_nmod_t));
            fq_nmod_init(sols->solution_sets[row][j][0],fq);
            fq_nmod_set(sols->solution_sets[row][j][0],point+j,fq);
            sols->solutions_per_var[row*n+j]=1;
        }
        return 1;
    }
    fq_nmod_mat_t ST,SB,inv,TS,rows,A;
    fq_nmod_mat_init(ST,k,dim,fq); fq_nmod_mat_transpose(ST,S,fq);
    slong rank=fq_nmod_mat_rref(ST,ST,fq);
    if (rank!=k) { fq_nmod_mat_clear(ST,fq); return 0; }
    slong *rr=flint_malloc(k*sizeof(slong));
    for (slong i=0;i<k;i++) {
        slong j=0; while (j<dim && fq_nmod_is_zero(fq_nmod_mat_entry(ST,i,j),fq)) j++;
        rr[i]=j;
    }
    fq_nmod_mat_clear(ST,fq);
    fq_nmod_mat_init(SB,k,k,fq); fq_nmod_mat_init(inv,k,k,fq);
    for (slong i=0;i<k;i++) for (slong j=0;j<k;j++)
        fq_nmod_set(fq_nmod_mat_entry(SB,i,j),fq_nmod_mat_entry(S,rr[i],j),fq);
    int ok=fq_nmod_mat_inv(inv,SB,fq); fq_nmod_mat_clear(SB,fq);
    if (!ok) { flint_free(rr); fq_nmod_mat_clear(inv,fq); return 0; }
    fq_nmod_mat_init(TS,dim,k,fq); fq_nmod_mat_mul(TS,T+variable,S,fq);
    fq_nmod_mat_init(rows,k,k,fq);
    for (slong i=0;i<k;i++) for (slong j=0;j<k;j++)
        fq_nmod_set(fq_nmod_mat_entry(rows,i,j),fq_nmod_mat_entry(TS,rr[i],j),fq);
    flint_free(rr); fq_nmod_mat_clear(TS,fq);
    fq_nmod_mat_init(A,k,k,fq); fq_nmod_mat_mul(A,inv,rows,fq);
    fq_nmod_mat_clear(inv,fq); fq_nmod_mat_clear(rows,fq);
    fq_nmod_poly_t mu; fq_nmod_poly_init(mu,fq); fq_nmod_mat_minpoly(mu,A,fq);
    fq_nmod_poly_factor_t roots; fq_nmod_poly_factor_init(roots,fq);
    fq_nmod_poly_roots(roots,mu,0,fq); fq_nmod_poly_clear(mu,fq);
    fq_nmod_t a; fq_nmod_init(a,fq);
    for (slong r=0;r<roots->num && ok;r++) {
        fq_nmod_poly_get_coeff(a,roots->poly+r,0,fq); fq_nmod_neg(a,a,fq);
        fq_nmod_set(point+variable,a,fq);
        fq_nmod_mat_t shifted,K,basis,next;
        fq_nmod_mat_init_set(shifted,A,fq);
        for (slong j=0;j<k;j++) fq_nmod_sub(fq_nmod_mat_entry(shifted,j,j),fq_nmod_mat_entry(shifted,j,j),a,fq);
        fq_nmod_mat_init(K,k,k,fq);
        slong nullity=fq_nmod_mat_nullspace(K,shifted,fq); fq_nmod_mat_clear(shifted,fq);
        if (!nullity) { fq_nmod_mat_clear(K,fq); ok=0; break; }
        fq_nmod_mat_window_init(basis,K,0,0,k,nullity,fq);
        fq_nmod_mat_init(next,dim,nullity,fq); fq_nmod_mat_mul(next,S,basis,fq);
        fq_nmod_mat_window_clear(basis,fq); fq_nmod_mat_clear(K,fq);
        ok=extension_points(sols,T,next,variable+1,point,fs,m,ctx);
        fq_nmod_mat_clear(next,fq);
    }
    fq_nmod_clear(a,fq); fq_nmod_poly_factor_clear(roots,fq); fq_nmod_mat_clear(A,fq);
    return ok;
}

static int extension_degree_fits(slong n, slong d, slong memory_mb,
                                 const quotient_field_t *e, slong *columns)
{
    slong N=monomial_count_bounded(n,d,WORD_MAX/(slong)sizeof(field_elem_u));
    if (!N) return -1;
    /* Echelon relations, normal forms, operators, certification products,
     * and the recursive FLINT eigenspaces. Account for fq limb allocations. */
    long double scalar=sizeof(field_elem_u)+3.0L*fq_nmod_ctx_degree(e->fq)*sizeof(ulong);
    long double bytes=(16.0L+4.0L*n)*N*N*scalar;
    /* Use the exact exclusive limit, including when long double is binary64. */
    if (bytes >= (long double)(UWORD(1) << (FLINT_BITS - 1))) return -1;
    if (memory_mb && bytes>(long double)memory_mb*1024*1024) return 0;
    *columns=N; return 1;
}

int solve_by_quotient_closure_extension(char **polys, slong count,
                            variable_info_t *vars, slong n,
                            polynomial_solutions_t *sols,
                            slong max_degree, slong memory_mb)
{
    if (max_degree<0 || memory_mb<0 || n<=0 || count<n) {
        sols->error_message=strdup("Quotient closure requires non-negative limits and at least as many equations as variables");
        return 0;
    }
    if (fq_nmod_ctx_degree(sols->ctx)==1)
        return solve_by_quotient_closure(polys,count,vars,n,sols,max_degree,memory_mb);
    quotient_field_t e; extension_init(&e,sols->ctx);
    fq_nmod_mpoly_ctx_t ctx; fq_nmod_mpoly_ctx_init(ctx,n,ORD_DEGLEX,e.fq);
    const char **names=flint_malloc(n*sizeof(char *));
    for (slong i=0;i<n;i++) names[i]=vars[i].name;
    fq_nmod_mpoly_struct *fs=flint_malloc(count*sizeof(*fs));
    quotient_input_t input={fs,ctx,extension_input_length,extension_input_term};
    slong *degrees=flint_malloc(count*sizeof(slong)),start=1,initialized=0;
    slong used=0,columns=0,dim=0;
    int success=0,certified=0;
    quotient_closure_t previous;
    int have_previous=0;
    unified_mat_struct *operators=NULL;
    const char *reason="degree limit reached";
    fq_nmod_poly_t mu; fq_nmod_poly_init(mu,e.fq);
    for (slong i=0;i<count;i++) {
        fq_nmod_mpoly_init(fs+i,ctx); initialized++;
        if (fq_nmod_mpoly_set_str_pretty(fs+i,polys[i],names,ctx) ||
                !fq_nmod_mpoly_total_degree_fits_si(fs+i,ctx)) {
            sols->error_message=strdup("Could not parse polynomial for quotient closure (or degree overflow)"); goto done;
        }
        degrees[i]=fs[i].length ? fq_nmod_mpoly_total_degree_si(fs+i,ctx) : -1;
        start=FLINT_MAX(start,degrees[i]);
        if (degrees[i]==0) certified=1;
    }
    for (slong d=start;!certified && (!max_degree || d<=max_degree);d++) {
        int fits=extension_degree_fits(n,d,memory_mb,&e,&columns);
        if (fits!=1) { reason=fits<0?"matrix dimensions exceed addressable range":"memory budget reached"; break; }
        used=d;
        quotient_closure_t cl; quotient_closure_init(&cl,n,d);
        if (have_previous) {
            quotient_closure_seed(&cl,&previous,&e);
            quotient_closure_clear(&previous,&e); have_previous=0;
        } else insert_generators(&cl,&input,count,degrees,&e);
        extend_relations(&cl,&e);
        dim=cl.basis.count-cl.rank;
        if (!dim) certified=1;
        else {
            operators=build_operators(&cl,&input,count,&e);
            certified=operators!=NULL;
        }
        if (certified) quotient_closure_clear(&cl,&e);
        else { previous=cl; have_previous=1; }
        if (d==WORD_MAX) { reason="degree exceeds integer range"; break; }
    }
    if (!certified) {
        char error[384]; snprintf(error,sizeof(error),
            "Quotient closure incomplete: %s (degree limit %ld, memory budget %ld MiB; 0 = unlimited; last attempted degree %ld). No conclusion about existence or dimension of solutions.",
            reason,max_degree,memory_mb,used);
        sols->error_message=strdup(error); goto done;
    }
    if (dim) {
        fq_nmod_mat_struct *T=flint_malloc(n*sizeof(*T));
        for (slong v=0;v<n;v++) {
            fq_nmod_mat_init(T+v,dim,dim,e.fq);
            for (slong i=0;i<dim;i++) for (slong j=0;j<dim;j++)
                field_elem_to_fq_nmod(fq_nmod_mat_entry(T+v,i,j),ENTRY(operators+v,i,j),&e.field);
        }
        fq_nmod_mat_minpoly(mu,T+n-1,e.fq);
        sols->solution_sets=calloc(dim,sizeof(fq_nmod_t **));
        sols->solutions_per_var=calloc(dim*n,sizeof(slong));
        fq_nmod_struct *point=flint_malloc(n*sizeof(*point));
        for (slong i=0;i<n;i++) fq_nmod_init(point+i,e.fq);
        fq_nmod_mat_t S; fq_nmod_mat_init(S,dim,dim,e.fq); fq_nmod_mat_one(S,e.fq);
        success=extension_points(sols,T,S,0,point,fs,count,ctx);
        fq_nmod_mat_clear(S,e.fq);
        for (slong i=0;i<n;i++) { fq_nmod_clear(point+i,e.fq); fq_nmod_mat_clear(T+i,e.fq); }
        flint_free(point); flint_free(T);
        if (!success) { sols->error_message=strdup("Quotient multiplication representation failed during rational-point extraction"); goto done; }
    } else fq_nmod_poly_one(mu,e.fq);
    success=1; sols->has_no_solutions=sols->num_solution_sets==0;
    sols->checked_solution_sets=sols->verified_solution_sets=sols->num_solution_sets;
    char summary[384]; snprintf(summary,sizeof(summary),
        "Quotient algebra over F_%lu^%ld (%s): certified dimension %ld, degree %ld, %ld monomials; %ld rational point(s)",
        fq_nmod_ctx_prime(e.fq),fq_nmod_ctx_degree(e.fq),field_ctx_get_description(&e.field),dim,used,columns,sols->num_solution_sets);
    sols->elimination_summary=strdup(summary);
    char *text=fq_nmod_poly_get_str_pretty(mu,names[n-1],e.fq);
    size_t len=strlen(text)+strlen(names[n-1])+96;
    sols->resultant_steps=calloc(1,sizeof(char *)); sols->resultant_steps[0]=malloc(len);
    snprintf(sols->resultant_steps[0],len,"Certified elimination polynomial in %s (minimal polynomial): %s",names[n-1],text);
    sols->num_resultant_steps=sols->resultant_steps_cap=1; flint_free(text);
 done:
    if (have_previous) quotient_closure_clear(&previous,&e);
    if (operators) { for (slong i=0;i<n;i++) unified_mat_clear(operators+i); flint_free(operators); }
    fq_nmod_poly_clear(mu,e.fq);
    for (slong i=0;i<initialized;i++) fq_nmod_mpoly_clear(fs+i,ctx);
    flint_free(fs); flint_free(degrees); flint_free(names);
    fq_nmod_mpoly_ctx_clear(ctx); field_ctx_clear(&e.field);
    return success;
}

/* Select a small set of actual ideal relations for polynomial division.
 * The set need not be a Groebner basis: only a zero remainder is a proof.
 * This also proves membership for candidates above the closure degree. */
static slong extension_divisors(fq_nmod_mpoly_struct ***out,
                                const quotient_closure_t *cl,
                                const fq_nmod_mpoly_ctx_t ctx,
                                const quotient_field_t *e)
{
    slong N=cl->basis.count,n=cl->basis.n,nr=0;
    slong *selected=flint_malloc(cl->rank*sizeof(slong));
    for (slong i=N;i-- > 0;) if (cl->rows[i]) {
        const ulong *exp=cl->basis.exps+i*n;
        int redundant=0;
        for (slong j=0;j<nr && !redundant;j++) {
            const ulong *q=cl->basis.exps+selected[j]*n;
            int divides=1;
            for (slong v=0;v<n;v++) if (q[v]>exp[v]) { divides=0; break; }
            redundant=divides;
        }
        if (!redundant) selected[nr++]=i;
    }
    fq_nmod_mpoly_struct **basis=flint_malloc(nr*sizeof(*basis));
    fq_nmod_t a; fq_nmod_init(a,e->fq);
    for (slong j=0;j<nr;j++) {
        basis[j]=flint_malloc(sizeof(**basis)); fq_nmod_mpoly_init(basis[j],ctx);
        slong i=selected[j];
        for (slong k=i;k<N;k++) if (!field_is_zero(cl->rows[i]+k,FIELD(e))) {
            field_elem_to_fq_nmod(a,cl->rows[i]+k,&e->field);
            fq_nmod_mpoly_push_term_fq_nmod_ui(basis[j],a,cl->basis.exps+k*n,ctx);
        }
        fq_nmod_mpoly_sort_terms(basis[j],ctx);
    }
    fq_nmod_clear(a,e->fq); flint_free(selected); *out=basis; return nr;
}

/* Bounded-degree echelon reduction uses every relation, including those
 * omitted by the minimal-divisor selection. */
static int extension_relation_contains(const quotient_closure_t *cl,
                         const fq_nmod_mpoly_t f,
                         const fq_nmod_mpoly_ctx_t ctx, const quotient_field_t *e)
{
    if (!fq_nmod_mpoly_total_degree_fits_si(f,ctx) ||
        fq_nmod_mpoly_total_degree_si(f,ctx)>cl->basis.degree) return 0;
    slong N=cl->basis.count,n=cl->basis.n;
    field_elem_u *row=vector_init(N,e),a;
    field_init_elem(&a,FIELD(e));
    fq_nmod_t coeff; fq_nmod_init(coeff,e->fq);
    for (slong i=0;i<f->length;i++) {
        ulong exp[n]; fq_nmod_mpoly_get_term_exp_ui(exp,f,i,ctx);
        slong col=basis_find(&cl->basis,exp); assert(col>=0);
        fq_nmod_mpoly_get_term_coeff_fq_nmod(coeff,f,i,ctx);
        fq_nmod_to_field_elem(row+col,coeff,&e->field);
    }
    int zero=1;
    for (slong j=0;j<N;j++) if (!field_is_zero(row+j,FIELD(e))) {
        if (!cl->rows[j]) { zero=0; break; }
        field_neg(&a,row+j,FIELD(e));
        vector_addmul(row+j,cl->rows[j]+j,N-j,&a,e);
    }
    fq_nmod_clear(coeff,e->fq); field_clear_elem(&a,FIELD(e)); vector_clear(row,N,e);
    return zero;
}

int quotient_sift_resultant_extension(fq_nmod_mpoly_t sifted,
                           fq_nmod_mpoly_t removed,
                           const fq_nmod_mpoly_t resultant,
                           fq_nmod_mpoly_struct *polys, slong count,
                           const fq_nmod_mpoly_ctx_t ctx,
                           slong max_degree, slong *used_degree)
{
    *used_degree=0;
    fq_nmod_mpoly_set(sifted,resultant,ctx); fq_nmod_mpoly_one(removed,ctx);
    slong n=ctx->minfo->nvars,start=1;
    if (max_degree<=0 || count<=0 || n<=0) return 0;
    if (fq_nmod_mpoly_is_fq_nmod(resultant,ctx)) return 1;
    if (!fq_nmod_mpoly_total_degree_fits_si(resultant,ctx)) return 0;
    slong *degrees=flint_malloc(count*sizeof(slong));
    for (slong i=0;i<count;i++) {
        if (!fq_nmod_mpoly_total_degree_fits_si(polys+i,ctx)) { flint_free(degrees); return 0; }
        degrees[i]=fq_nmod_mpoly_total_degree_si(polys+i,ctx);
        start=FLINT_MAX(start,degrees[i]);
    }
    quotient_field_t e; extension_init(&e,ctx->fqctx);
    quotient_input_t input={polys,ctx,extension_input_length,extension_input_term};
    slong columns;
    if (start>max_degree || extension_degree_fits(n,start,0,&e,&columns)!=1) {
        flint_free(degrees); field_ctx_clear(&e.field); return 1;
    }
    fq_nmod_mpoly_factor_t factors; fq_nmod_mpoly_factor_init(factors,ctx);
    if (!fq_nmod_mpoly_factor(factors,resultant,ctx)) {
        fq_nmod_mpoly_factor_clear(factors,ctx); flint_free(degrees);
        field_ctx_clear(&e.field); return 0;
    }
    slong nf=fq_nmod_mpoly_factor_length(factors,ctx);
    slong *left=flint_malloc(nf*sizeof(slong));
    for (slong i=0;i<nf;i++) left[i]=fq_nmod_mpoly_factor_get_exp_si(factors,i,ctx);
    fq_nmod_mpoly_t candidate,factor,remainder;
    fq_nmod_mpoly_init(candidate,ctx); fq_nmod_mpoly_init(factor,ctx);
    fq_nmod_mpoly_init(remainder,ctx);
    for (slong d=start;d<=max_degree;d++) {
        if (extension_degree_fits(n,d,0,&e,&columns)!=1) break;
        quotient_closure_t cl; quotient_closure_init(&cl,n,d);
        insert_generators(&cl,&input,count,degrees,&e); extend_relations(&cl,&e);
        *used_degree=d;
        fq_nmod_mpoly_struct **basis;
        slong nr=extension_divisors(&basis,&cl,ctx,&e);
        fq_nmod_mpoly_struct **quotients=flint_malloc(nr*sizeof(*quotients));
        for (slong j=0;j<nr;j++) {
            quotients[j]=flint_malloc(sizeof(**quotients)); fq_nmod_mpoly_init(quotients[j],ctx);
        }
        for (slong i=0;i<nf;i++) {
            fq_nmod_mpoly_factor_get_base(factor,factors,i,ctx);
            while (left[i]>0) {
                int divides=fq_nmod_mpoly_divides(candidate,sifted,factor,ctx);
                assert(divides);
                if (!divides) break;
                int contains=extension_relation_contains(&cl,candidate,ctx,&e);
                if (!contains && nr) {
                    fq_nmod_mpoly_divrem_ideal(quotients,remainder,candidate,basis,nr,ctx);
                    contains=fq_nmod_mpoly_is_zero(remainder,ctx);
                }
                if (!contains) break;
                fq_nmod_mpoly_swap(sifted,candidate,ctx);
                fq_nmod_mpoly_mul(removed,removed,factor,ctx); left[i]--;
            }
        }
        for (slong j=0;j<nr;j++) {
            fq_nmod_mpoly_clear(basis[j],ctx); flint_free(basis[j]);
            fq_nmod_mpoly_clear(quotients[j],ctx); flint_free(quotients[j]);
        }
        flint_free(basis); flint_free(quotients); quotient_closure_clear(&cl,&e);
        if (fq_nmod_mpoly_is_fq_nmod(sifted,ctx) || d==WORD_MAX) break;
    }
    fq_nmod_mpoly_clear(candidate,ctx); fq_nmod_mpoly_clear(factor,ctx);
    fq_nmod_mpoly_clear(remainder,ctx); fq_nmod_mpoly_factor_clear(factors,ctx);
    flint_free(left); flint_free(degrees); field_ctx_clear(&e.field); return 1;
}

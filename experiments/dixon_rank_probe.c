/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Research diagnostic, NOT a replacement for the production rank model.
 * Usage: build/dixon_rank_probe equations degree [prime [seed]]
 * Random dense affine inputs represent a generic specialization of the hidden
 * parameter. All ranks are computed exactly over the supplied prime field.
 * Sampling gives evidence about generic rank, not a proof of a generic upper
 * bound. The four-equation correction below remains a conjecture.
 */
#include "../src/dixon/dixon_flint.c"
#include <assert.h>

static void dense_input(unified_mpoly_struct *p, slong *exp, slong axis,
                        slong remaining, flint_rand_t state)
{
    if (axis == p->nvars) {
        fq_nmod_t c;
        fq_nmod_init(c, p->ctx);
        fq_nmod_randtest_not_zero(c, state, p->ctx);
        slong par = 0;
        dr_mpoly_add_term_fast(p, exp, &par, c);
        fq_nmod_clear(c, p->ctx);
        return;
    }
    for (slong e = 0; e <= remaining; e++) {
        exp[axis] = e;
        dense_input(p, exp, axis + 1, remaining - e, state);
    }
}

static slong four_equation_defect(slong d, slong p, slong q)
{
    slong defect = 0;
    /* Tail of [(1-z^d)^4/(1-z)^2]_+, starting at max(p,q)-d+1. */
    for (slong a = FLINT_MAX(p,q)-d+1; a < 2*d; a++)
        if (a >= d) defect += FLINT_MAX(0,4*d-3*a-3);
    return defect;
}

int main(int argc, char **argv)
{
    if (argc < 3 || argc > 5) {
        fprintf(stderr,"Usage: %s equations degree [prime [seed]]\n",argv[0]);
        return 2;
    }
    slong n = atol(argv[1]), d = atol(argv[2]), m = n-1;
    ulong prime = argc > 3 ? strtoul(argv[3],NULL,10) : 65537;
    ulong seed = argc > 4 ? strtoul(argv[4],NULL,10) : 12345;
    if (n < 3 || n > 8 || d < 2 || !n_is_prime(prime)) return 2;
    omp_set_num_threads(4);
    g_dixon_verbose_level = 0;
    flint_rand_t state;
    flint_rand_init(state); flint_rand_set_seed(state,seed,808);
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx,prime,1,"a");
    unified_mpoly_struct *polys = flint_calloc(n,sizeof(*polys));
    slong *exp = flint_calloc(m,sizeof(*exp));
    long *degrees = flint_malloc(n*sizeof(*degrees));
    for (slong i = 0; i < n; i++) {
        degrees[i] = d;
        dr_mpoly_init(polys+i,m,1,ctx);
        dense_input(polys+i,exp,0,d,state);
    }
    unified_mpoly_struct **matrix, **a, full = {0};
    build_fq_cancellation_matrix(&matrix,polys,m,1);
    perform_fq_matrix_row_operations(&a,&matrix,m,1);
    compute_fq_cancel_matrix_det(&full,a,m,1,DET_METHOD_RECURSIVE);
    monom_t *rm = NULL, *cm = NULL;
    slong nr=0,nc=0,rc=0,cc=0,rhs=16,chs=16;
    hash_entry_t **ri=flint_calloc(rhs,sizeof(*ri)), **ci=flint_calloc(chs,sizeof(*ci));
    slong nt=dr_mpoly_length(&full);
    slong *rows=flint_malloc(nt*sizeof(*rows)), *cols=flint_malloc(nt*sizeof(*cols));
    for (slong t=0; t<nt; t++) {
        DR_MPOLY_TERM(term,&full,t);
        rows[t]=dixon_intern_monom(&rm,&nr,&rc,&ri,&rhs,term.var_exp,m);
        cols[t]=dixon_intern_monom(&cm,&nc,&cc,&ci,&chs,term.var_exp+m,m);
    }
    slong *R,*H,rl,hl,sigma,rho;
    assert(dixon_rank_profile_from_degrees(&R,&rl,&H,&hl,&sigma,&rho,degrees,n,m));
    slong *rd=flint_calloc(nr,sizeof(*rd)), *cd=flint_calloc(nc,sizeof(*cd));
    slong h=0,boundary_rank=0;
    for(slong i=0;i<nr;i++)for(slong v=0;v<m;v++)rd[i]+=rm[i].exp[v];
    for(slong j=0;j<nc;j++)for(slong v=0;v<m;v++)cd[j]+=cm[j].exp[v];
    for(slong i=0;i<hl;i++)h+=H[i];
    nmod_mat_t B,copy;
    nmod_mat_init(B,nr,nc,prime);
    for(slong t=0;t<nt;t++) {
        DR_MPOLY_TERM(term,&full,t);
        assert(!term.par_exp || !term.par_exp[0]);
        nmod_mat_entry(B,rows[t],cols[t])=nmod_poly_get_coeff_ui(term.coeff,0);
    }
    nmod_mat_init_set(copy,B);
    slong full_rank=nmod_mat_rank(copy);
    nmod_mat_clear(copy);
    printf("{\"n\":%ld,\"d\":%ld,\"prime\":%lu,\"seed\":%lu,"
           "\"rows\":%ld,\"cols\":%ld,\"predicted\":%ld,\"rank\":%ld,"
           "\"h\":%ld,\"sigma\":%ld,\"blocks\":[",n,d,prime,seed,nr,nc,rho,full_rank,h,sigma);
    int first=1;
    for(slong p=0;p<rl;p++) {
        slong q=sigma-p;
        if(q<0 || q>=rl)continue;
        slong rn=0,cn=0;
        for(slong i=0;i<nr;i++)rn+=rd[i]==p;
        for(slong j=0;j<nc;j++)cn+=cd[j]==q;
        nmod_mat_t block,pattern;
        nmod_mat_init(block,rn,cn,prime); nmod_mat_init(pattern,rn,cn,prime);
        slong r=0;
        for(slong i=0;i<nr;i++)if(rd[i]==p) {
            slong c=0;
            for(slong j=0;j<nc;j++)if(cd[j]==q) {
                ulong value=nmod_mat_entry(B,i,j);
                nmod_mat_entry(block,r,c)=value;
                if(value)nmod_mat_entry(pattern,r,c)=1+n_randint(state,prime-1);
                c++;
            }
            r++;
        }
        slong rank=nmod_mat_rank(block), pattern_rank=nmod_mat_rank(pattern);
        slong k=FLINT_MIN(p,q), predicted=FLINT_MAX(0,FLINT_MIN(R[p],R[q])-(k<hl?H[k]:0));
        boundary_rank+=rank;
        printf("%s{\"p\":%ld,\"q\":%ld,\"rows\":%ld,\"cols\":%ld,"
               "\"predicted\":%ld,\"rank\":%ld,\"randomized_pattern_rank\":%ld",
               first?"":",",p,q,rn,cn,predicted,rank,pattern_rank);
        if(n==4)printf(",\"candidate_rank\":%ld",predicted-four_equation_defect(d,p,q));
        printf("}");first=0;
        nmod_mat_clear(block);nmod_mat_clear(pattern);
    }
    printf("],\"boundary_rank\":%ld,\"residual_rank\":%ld",boundary_rank,full_rank-boundary_rank);
    if(n==4) {
        slong k=(d-1)/3,delta=k*(k+1)*(d-2*k-1);
        printf(",\"candidate_rank\":%ld,\"candidate_defect\":%ld",(5*d*d*d-2*d)/3-delta,delta);
    }
    printf("}\n");
    nmod_mat_clear(B);
    flint_free(rd);flint_free(cd);flint_free(R);flint_free(H);
    free_monom_index(ri,rhs);free_monom_index(ci,chs);
    for(slong i=0;i<nr;i++)flint_free(rm[i].exp);
    for(slong j=0;j<nc;j++)flint_free(cm[j].exp);
    flint_free(rm);flint_free(cm);flint_free(rows);flint_free(cols);
    dr_mpoly_clear(&full);
    for(slong i=0;i<n;i++) {
        dr_mpoly_clear(polys+i);
        for(slong j=0;j<n;j++) {
            dr_mpoly_clear(&matrix[i][j]);dr_mpoly_clear(&a[i][j]);
        }
        flint_free(matrix[i]);flint_free(a[i]);
    }
    flint_free(matrix);flint_free(a);flint_free(polys);flint_free(exp);flint_free(degrees);
    fq_nmod_ctx_clear(ctx);flint_rand_clear(state);flint_cleanup_master();
    return 0;
}

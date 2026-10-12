/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Adaptive Macaulay closure in shrinking quotient coordinates over F_p.
 * Derived from the independently audited standalone closure experiments.
 * All inserted relations lie in I. A successful algebra certificate proves
 * that the resulting multiplication representation is exactly F_p[x]/I.
 * No finite-field equations or squarefree/radical assumptions are imposed.
 */
#include "quotient_solver.h"
#include <flint/nmod_mpoly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <flint/nmod_mpoly_factor.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>

#include "quotient_basis.h"

static void matrix_init(nmod_mat_t M, slong r, slong c, ulong p)
{
    assert(r>=0 && c>=0); /* Dimensions are checked before each degree. */
    nmod_mat_init(M,r,c,p);
}

static void select_matrix(nmod_mat_t dst, const nmod_mat_t src,
                          const slong *rows, slong nr, const slong *cols, slong nc)
{
    matrix_init(dst,nr,nc,src->mod.n);
    for (slong i=0;i<nr;i++) for (slong j=0;j<nc;j++)
        nmod_mat_entry(dst,i,j)=nmod_mat_entry(src,rows?rows[i]:i,cols?cols[j]:j);
}

static slong *pivots(const nmod_mat_t E)
{
    slong *p=flint_malloc(E->r*sizeof(slong));
    for (slong i=0;i<E->r;i++) {
        slong j=0; while (j<E->c && !nmod_mat_entry(E,i,j)) j++;
        assert(j<E->c && nmod_mat_entry(E,i,j)==1); p[i]=j;
    }
    return p;
}

/* RREF with implicit identity pivot columns. p and f are the sorted
 * pivot and free columns; only the nontrivial tails are stored. */
typedef struct {
    slong r, c;
    slong *p, *f;
    nmod_mat_t tail;
    nmod_t mod;
} relations_struct;
typedef relations_struct relations_t[1];

static void relations_init(relations_t E, slong nc, ulong prime)
{
    E->r=0; E->c=nc; E->p=NULL;
    E->f=flint_malloc(nc*sizeof(slong));
    for (slong j=0;j<nc;j++) E->f[j]=j;
    matrix_init(E->tail,0,nc,prime); E->mod=E->tail->mod;
}

static void relations_clear(relations_t E)
{
    flint_free(E->p); flint_free(E->f);
    nmod_mat_clear(E->tail);
}

static void relations_projection(nmod_mat_t W, const relations_t E)
{
    slong nf=E->c-E->r;
    matrix_init(W,E->c,nf,E->mod.n);
    for (slong j=0;j<nf;j++) {
        nmod_mat_entry(W,E->f[j],j)=1;
        for (slong i=0;i<E->r;i++)
            nmod_mat_entry(W,E->p[i],j)=nmod_neg(nmod_mat_entry(E->tail,i,j),E->mod);
    }
}

typedef struct {
    basis_t basis;
    relations_t E;
    /* Only pivots are needed to identify a complement of the already
     * processed subspace; its full matrix need not be retained. */
    slong *processed, nprocessed;
    slong *shift, submitted, rref_calls;
    int original_only, processed_ready, inherited;
    uint64_t rref_cells, matmul_work;
} closure_t;

static void compact_rref(nmod_mat_t M, closure_t *cl)
{
    if (cl) { cl->rref_cells+=(uint64_t)M->r*M->c; cl->rref_calls++; }
    slong rank=nmod_mat_rref(M);
    nmod_mat_t small; select_matrix(small,M,NULL,rank,NULL,M->c);
    nmod_mat_swap(M,small); nmod_mat_clear(small);
}

static void closure_init(closure_t *cl, slong n, slong degree, ulong p)
{
    memset(cl,0,sizeof(*cl)); cl->original_only=1; basis_init(&cl->basis,n,degree);
    slong nc=cl->basis.count;
    relations_init(cl->E,nc,p);
    cl->shift=flint_malloc(nc*n*sizeof(slong));
    for (slong j=0;j<nc;j++) for (slong v=0;v<n;v++) {
        ulong e[n]; memcpy(e,cl->basis.exps+j*n,n*sizeof(ulong)); e[v]++;
        cl->shift[j*n+v]=basis_find(&cl->basis,e);
    }
}

static void closure_clear(closure_t *cl)
{
    basis_clear(&cl->basis); relations_clear(cl->E); flint_free(cl->processed); flint_free(cl->shift);
}

/* Retain the previous degree's relations and one variable multiple per
 * distinct leading monomial. Multiplication preserves the monomial order,
 * so these rows are already unit upper triangular in their pivot columns.
 * Only solve for the free columns; a fresh dense LU would discard this
 * structure. Every seed row still belongs to the original ideal. */
static void closure_seed(closure_t *cl, const closure_t *prev)
{
    slong nc=cl->basis.count, pc=prev->basis.count, n=cl->basis.n;
    slong offset=nc-pc, rank=0, nf=0;
    cl->inherited=1;
    cl->original_only=prev->original_only;
    /* If no closure relation has been added, the previous space is exactly
     * the Macaulay space at degree d-1. All its variable multiples are
     * covered by insert_original_rows at degree d, so reuse its RREF. */
    if (prev->original_only) {
        cl->nprocessed=prev->E->r;
        cl->processed=flint_malloc(cl->nprocessed*sizeof(slong));
        for (slong i=0;i<cl->nprocessed;i++) cl->processed[i]=prev->E->p[i]+offset;
        cl->processed_ready=1;
    }
    const slong *p=prev->E->p;
    slong *source=flint_malloc(nc*sizeof(slong));
    slong *variable=flint_malloc(nc*sizeof(slong));
    slong *index=flint_malloc(nc*sizeof(slong));
    slong *cols=flint_malloc(nc*sizeof(slong));
    for (slong j=0;j<nc;j++) source[j]=-1;
    for (slong i=0;i<prev->E->r;i++) {
        slong j=p[i]+offset; source[j]=i; variable[j]=-1;
    }
    for (slong i=0;i<prev->E->r;i++) for (slong v=0;v<n;v++) {
        slong j=cl->shift[(p[i]+offset)*n+v];
        assert(j>=0);
        if (source[j]<0) { source[j]=i; variable[j]=v; }
    }
    for (slong j=0;j<nc;j++) {
        if (source[j]>=0) index[j]=rank++;
        else { cols[nf]=j; index[j]=-(++nf); }
    }
    /* The old RREF has nonzero tails only in its free columns. Back
     * substitution through those sparse tails avoids a rank-by-rank U. */
    slong pn=pc-prev->E->r;
    const slong *pf=prev->E->f;
    nmod_mat_t V;
    matrix_init(V,rank,nf,cl->E->mod.n);
    for (slong j=nc-1;nf && j>=0;j--) if (source[j]>=0) {
        slong row=index[j], src=source[j], v=variable[j];
        for (slong t=0;t<pn;t++) {
            ulong a=nmod_mat_entry(prev->E->tail,src,t);
            if (!a) continue;
            slong dest=v<0?pf[t]+offset:cl->shift[(pf[t]+offset)*n+v];
            assert(dest>j);
            if (index[dest]>=0) {
                _nmod_vec_scalar_addmul_nmod(nmod_mat_entry_ptr(V,row,0),
                    nmod_mat_entry_ptr(V,index[dest],0),nf,nmod_neg(a,V->mod),V->mod);
                cl->matmul_work+=nf;
            } else {
                slong col=-index[dest]-1;
                nmod_mat_entry(V,row,col)=nmod_add(nmod_mat_entry(V,row,col),a,V->mod);
            }
        }
    }
    cl->E->r=rank;
    cl->E->p=flint_realloc(cl->E->p,rank*sizeof(slong));
    for (slong j=0;j<nc;j++) if (source[j]>=0) cl->E->p[index[j]]=j;
    memcpy(cl->E->f,cols,nf*sizeof(slong));
    nmod_mat_swap(cl->E->tail,V); nmod_mat_clear(V);
    flint_free(source); flint_free(variable); flint_free(index); flint_free(cols);
}

/* Merge relations expressed only in the current free columns. */
static void insert_reduced_rows(closure_t *cl, nmod_mat_t S)
{
    relations_struct *E=cl->E;
    slong rank=E->r, nf=S->c;
    assert(nf==E->c-rank);
    compact_rref(S,cl);
    if (!S->r) return;
    slong *np=pivots(S), *rest=flint_malloc(nf*sizeof(slong));
    slong nr=0, pi=0;
    for (slong j=0;j<nf;j++) {
        if (pi<S->r && np[pi]==j) pi++;
        else rest[nr++]=j;
    }
    nmod_mat_t tail,En,product,updated,joined;
    select_matrix(tail,S,NULL,S->r,rest,nr);
    select_matrix(En,E->tail,NULL,rank,np,S->r);
    matrix_init(product,rank,nr,S->mod.n); nmod_mat_mul(product,En,tail);
    select_matrix(updated,E->tail,NULL,rank,rest,nr); nmod_mat_sub(updated,updated,product);
    cl->matmul_work+=(uint64_t)rank*S->r*nr;
    for (slong i=0;i<S->r;i++) np[i]=E->f[np[i]];
    slong *p=flint_malloc((rank+S->r)*sizeof(slong));
    matrix_init(joined,rank+S->r,nr,S->mod.n);
    slong old=0,add=0;
    for (slong i=0;i<joined->r;i++) {
        int use_old=old<rank && (add==S->r || E->p[old]<np[add]);
        p[i]=use_old?E->p[old]:np[add];
        for (slong j=0;j<nr;j++)
            nmod_mat_entry(joined,i,j)=use_old?nmod_mat_entry(updated,old,j):nmod_mat_entry(tail,add,j);
        if (use_old) old++; else add++;
    }
    for (slong i=0;i<nr;i++) E->f[i]=E->f[rest[i]];
    flint_free(E->p); E->p=p; E->r+=S->r;
    nmod_mat_swap(E->tail,joined);
    nmod_mat_clear(joined); nmod_mat_clear(tail); nmod_mat_clear(En);
    nmod_mat_clear(product); nmod_mat_clear(updated); flint_free(np); flint_free(rest);
}

static void insert_rows(closure_t *cl, const nmod_mat_t M)
{
    relations_struct *E=cl->E;
    slong rank=E->r, nf=E->c-rank; cl->submitted+=M->r;
    if (!M->r || !nf) return;
    nmod_mat_t S;
    matrix_init(S,M->r,nf,M->mod.n);
    /* Bound the pivot-column scratch space even for tall batches. */
    for (slong start=0;start<M->r;start+=256) {
        slong rows=FLINT_MIN(256,M->r-start);
        nmod_mat_t Mp,block;
        matrix_init(Mp,rows,rank,M->mod.n);
        for (slong i=0;i<rows;i++) for (slong j=0;j<rank;j++)
            nmod_mat_entry(Mp,i,j)=nmod_mat_entry(M,start+i,E->p[j]);
        nmod_mat_window_init(block,S,start,0,start+rows,nf);
        nmod_mat_mul(block,Mp,E->tail);
        for (slong i=0;i<rows;i++) for (slong j=0;j<nf;j++)
            nmod_mat_entry(block,i,j)=nmod_sub(nmod_mat_entry(M,start+i,E->f[j]),nmod_mat_entry(block,i,j),M->mod);
        nmod_mat_window_clear(block); nmod_mat_clear(Mp);
    }
    cl->matmul_work+=(uint64_t)M->r*rank*nf;
    insert_reduced_rows(cl,S); nmod_mat_clear(S);
}

static void closure_run(closure_t *cl, const nmod_mat_t initial,
                        const nmod_mpoly_ctx_t ctx)
{
    slong n=cl->basis.n,nc=cl->basis.count;
    insert_rows(cl,initial);
    for (;;) {
        if (cl->E->r==nc) break;
        slong *low=flint_malloc(cl->E->r*sizeof(slong));
        slong *delta=flint_malloc(cl->E->r*sizeof(slong));
        slong nl=0,nd=0,np=0;
        /* processed is a subspace of the low-degree space. Its pivots
         * are a subset of the latter's; rows with new pivots form a
         * complement even when old pivot rows have changed tails. */
        for (slong i=0;i<cl->E->r;i++) {
            slong pivot=cl->E->p[i];
            if (cl->basis.degrees[pivot]>=cl->basis.degree) continue;
            low[nl++]=pivot;
            if (np<cl->nprocessed && cl->processed[np]==pivot) np++;
            else delta[nd++]=i;
        }
        assert(np==cl->nprocessed);
        flint_free(cl->processed); cl->processed=low; cl->nprocessed=nl;
        if (!nd) { flint_free(delta); break; }
        nmod_mat_t extra; matrix_init(extra,nd*n,nc,ctx->mod.n);
        for (slong i=0;i<nd;i++) {
            slong row=delta[i];
            for (slong v=0;v<n;v++)
                nmod_mat_entry(extra,i*n+v,cl->shift[cl->E->p[row]*n+v])=1;
            for (slong j=0;j<cl->E->tail->c;j++) {
                ulong a=nmod_mat_entry(cl->E->tail,row,j);
                if (!a) continue;
                for (slong v=0;v<n;v++) {
                    slong dest=cl->shift[cl->E->f[j]*n+v]; assert(dest>=0);
                    nmod_mat_entry(extra,i*n+v,dest)=a;
                }
            }
        }
        slong oldrank=cl->E->r; insert_rows(cl,extra);
        nmod_mat_clear(extra); flint_free(delta);
        if (cl->E->r==oldrank) break;
        cl->original_only=0;
    }
}

static void original_rows(nmod_mat_t M, nmod_mpoly_struct *fs, slong m,
                          const nmod_mpoly_ctx_t ctx, const basis_t *basis, int all, slong bound)
{
    basis_t multipliers; basis_init(&multipliers,basis->n,all?basis->degree:0);
    slong count=0;
    for (slong i=0;i<m;i++) if (fs[i].length) {
        slong d=nmod_mpoly_total_degree_si(fs+i,ctx); assert(d<=basis->degree);
        for (slong j=0;j<multipliers.count;j++) if (!all || multipliers.degrees[j]<=bound-d) count++;
    }
    matrix_init(M,count,basis->count,ctx->mod.n); slong row=0;
    for (slong i=0;i<m;i++) if (fs[i].length) {
        slong d=nmod_mpoly_total_degree_si(fs+i,ctx);
        for (slong j=0;j<multipliers.count;j++) if (!all || multipliers.degrees[j]<=bound-d) {
            for (slong a=0;a<fs[i].length;a++) {
                ulong e[basis->n]; nmod_mpoly_get_term_exp_ui(e,fs+i,a,ctx);
                for (slong v=0;v<basis->n;v++) e[v]+=multipliers.exps[j*basis->n+v];
                slong col=basis_find(basis,e); assert(col>=0);
                nmod_mat_entry(M,row,col)=fs[i].coeffs[a];
            }
            row++;
        }
    }
    assert(row==count); basis_clear(&multipliers);
}

typedef struct {
    int certified;
    slong dimension;
    nmod_mat_struct *operators;
    nmod_poly_t eliminant;
} quotient_result;

static void projection(nmod_mat_t W, const nmod_mat_t E, slong **free_out)
{
    slong *p=pivots(E), *f=flint_malloc(E->c*sizeof(slong)), nf=0, pi=0;
    for (slong j=0;j<E->c;j++) {
        if (pi<E->r && p[pi]==j) pi++; else f[nf++]=j;
    }
    matrix_init(W,E->c,nf,E->mod.n);
    for (slong j=0;j<nf;j++) {
        nmod_mat_entry(W,f[j],j)=1;
        for (slong i=0;i<E->r;i++) nmod_mat_entry(W,p[i],j)=nmod_neg(nmod_mat_entry(E,i,f[j]),E->mod);
    }
    flint_free(p); *free_out=f;
}

/* Construct Macaulay relations directly in quotient coordinates. Each
 * input multiple has only fs[i].length nonzeros, even when the ambient
 * monomial basis has thousands of columns. No dense ambient batch is needed.
 * After seeding, all multiples of lower degree have already been inherited. */
static void insert_original_rows(closure_t *cl, nmod_mpoly_struct *fs, slong m,
                                  const nmod_mpoly_ctx_t ctx)
{
    if (cl->E->r==cl->basis.count) return;
    if (!cl->E->r) {
        nmod_mat_t M; original_rows(M,fs,m,ctx,&cl->basis,1,cl->basis.degree);
        insert_rows(cl,M); nmod_mat_clear(M); return;
    }
    slong n=cl->basis.n, degree=cl->basis.degree, md=0, count=0;
    for (slong i=0;i<m;i++) if (fs[i].length)
        md=FLINT_MAX(md,degree-nmod_mpoly_total_degree_si(fs+i,ctx));
    basis_t multipliers; basis_init(&multipliers,n,md);
    for (slong i=0;i<m;i++) if (fs[i].length) {
        slong d=nmod_mpoly_total_degree_si(fs+i,ctx);
        for (slong j=0;j<multipliers.count;j++) if (multipliers.degrees[j]<=degree-d &&
            (!cl->inherited || multipliers.degrees[j]==degree-d)) count++;
    }
    nmod_mat_t W,S; relations_projection(W,cl->E);
    matrix_init(S,count,W->c,ctx->mod.n);
    slong row=0;
    for (slong i=0;i<m;i++) if (fs[i].length) {
        slong d=nmod_mpoly_total_degree_si(fs+i,ctx);
        ulong *exps=flint_malloc(fs[i].length*n*sizeof(ulong));
        for (slong a=0;a<fs[i].length;a++) nmod_mpoly_get_term_exp_ui(exps+a*n,fs+i,a,ctx);
        for (slong j=0;j<multipliers.count;j++) if (multipliers.degrees[j]<=degree-d &&
            (!cl->inherited || multipliers.degrees[j]==degree-d)) {
            for (slong a=0;a<fs[i].length;a++) {
                ulong e[n];
                for (slong v=0;v<n;v++) e[v]=exps[a*n+v]+multipliers.exps[j*n+v];
                slong col=basis_find(&cl->basis,e); assert(col>=0);
                _nmod_vec_scalar_addmul_nmod(nmod_mat_entry_ptr(S,row,0),
                    nmod_mat_entry_ptr(W,col,0),W->c,fs[i].coeffs[a],ctx->mod);
            }
            row++;
        }
        flint_free(exps);
    }
    assert(row==count); cl->submitted+=count;
    nmod_mat_clear(W); basis_clear(&multipliers);
    insert_reduced_rows(cl,S);
    nmod_mat_clear(S);
}

static void row_times(ulong *out, const ulong *in, const nmod_mat_t T)
{
    memset(out,0,T->c*sizeof(ulong));
    for (slong i=0;i<T->r;i++) if (in[i]) for (slong j=0;j<T->c;j++)
        out[j]=nmod_add(out[j],nmod_mul(in[i],nmod_mat_entry(T,i,j),T->mod),T->mod);
}

static void act_monomial(ulong *out, const ulong *one, const ulong *e,
                         nmod_mat_struct *T, slong n, slong dim)
{
    ulong *tmp=flint_malloc(dim*sizeof(ulong)); memcpy(out,one,dim*sizeof(ulong));
    for (slong i=0;i<n;i++) for (ulong k=0;k<e[i];k++) {
        row_times(tmp,out,T+i); memcpy(out,tmp,dim*sizeof(ulong));
    }
    flint_free(tmp);
}

static int certify_operators(quotient_result *result, nmod_mat_struct *T,
    const ulong *one, const slong *B, const basis_t *basis,
    nmod_mpoly_struct *fs, slong m, const nmod_mpoly_ctx_t ctx)
{
    slong dim=T[0].r,n=basis->n; int ok=1;
    nmod_mat_t AB,BA; matrix_init(AB,dim,dim,ctx->mod.n); matrix_init(BA,dim,dim,ctx->mod.n);
    for (slong i=0;i<n && ok;i++) for (slong j=0;j<i && ok;j++) {
        nmod_mat_mul(AB,T+i,T+j); nmod_mat_mul(BA,T+j,T+i); ok=nmod_mat_equal(AB,BA);
    }
    ulong *v=flint_malloc(dim*sizeof(ulong)),*sum=flint_malloc(dim*sizeof(ulong));
    for (slong i=0;i<dim && ok;i++) {
        act_monomial(v,one,basis->exps+B[i]*n,T,n,dim);
        for (slong j=0;j<dim;j++) if (v[j]!=(ulong)(i==j)) ok=0;
    }
    for (slong i=0;i<m && ok;i++) {
        memset(sum,0,dim*sizeof(ulong));
        for (slong j=0;j<fs[i].length;j++) {
            ulong e[n]; nmod_mpoly_get_term_exp_ui(e,fs+i,j,ctx);
            act_monomial(v,one,e,T,n,dim);
            for (slong k=0;k<dim;k++) sum[k]=nmod_add(sum[k],nmod_mul(fs[i].coeffs[j],v[k],ctx->mod),ctx->mod);
        }
        for (slong j=0;j<dim;j++) if (sum[j]) ok=0;
    }
    if (ok) {
        result->certified=1;result->dimension=dim;
        result->operators=flint_malloc(n*sizeof(*result->operators));
        for(slong i=0;i<n;i++) nmod_mat_init_set(result->operators+i,T+i);
        nmod_mat_minpoly(result->eliminant,T+n-1);
    }
    flint_free(v);flint_free(sum);nmod_mat_clear(AB);nmod_mat_clear(BA);return ok;
}

static void quotient_close(quotient_result *result, closure_t *cl,
                            nmod_mpoly_struct *fs, slong m, const nmod_mpoly_ctx_t ctx)
{
    slong n=cl->basis.n,nc=cl->basis.count;slong *freecols;
    nmod_mat_t W;relations_projection(W,cl->E);
    slong *low=flint_malloc(nc*sizeof(slong)),nl=0;
    for (slong j=nc-1;j>=0;j--) if (cl->basis.degrees[j]<cl->basis.degree) low[nl++]=j;
    slong *shift=flint_malloc(n*nl*sizeof(slong));
    for (slong i=0;i<n;i++) for (slong j=0;j<nl;j++) shift[i*nl+j]=cl->shift[low[j]*n+i];
    for (;;) {
        slong dim=W->c;
        int one_nonzero=0;
        for(slong j=0;j<dim;j++)one_nonzero|=nmod_mat_entry(W,nc-1,j)!=0;
        if (!dim || !one_nonzero) {
            result->certified=1;result->dimension=0;nmod_poly_one(result->eliminant);break;
        }
        nmod_mat_t WL,WT;select_matrix(WL,W,low,nl,NULL,dim);matrix_init(WT,dim,nl,ctx->mod.n);
        nmod_mat_transpose(WT,WL);compact_rref(WT,NULL);
        if (WT->r!=dim) {nmod_mat_clear(WL);nmod_mat_clear(WT);break;}
        slong *local=pivots(WT),*B=flint_malloc(dim*sizeof(slong));
        for(slong j=0;j<dim;j++)B[j]=low[local[j]];
        nmod_mat_clear(WT);
        nmod_mat_t WB,inv,newW;select_matrix(WB,W,B,dim,NULL,dim);matrix_init(inv,dim,dim,ctx->mod.n);
        int invertible=nmod_mat_inv(inv,WB);
        if(!invertible) {
            nmod_mat_clear(inv); nmod_mat_clear(WB); nmod_mat_clear(WL);
            flint_free(local); flint_free(B); break;
        }
        matrix_init(newW,nc,dim,ctx->mod.n);nmod_mat_mul(newW,W,inv);
        nmod_mat_swap(W,newW);nmod_mat_clear(newW);nmod_mat_clear(inv);nmod_mat_clear(WB);
        nmod_mat_clear(WL);select_matrix(WL,W,low,nl,NULL,dim);
        nmod_mat_struct *T=flint_malloc(n*sizeof(*T));
        /* Build all operators, but stream consistency relations. A nonzero
         * block already gives valid ideal relations and strictly lowers
         * the dimension. Restart in the smaller coordinates immediately.
         * Certification still requires a complete pass with no relations. */
        for(slong i=0;i<n;i++) {
            slong *sb=flint_malloc(dim*sizeof(slong));
            for(slong j=0;j<dim;j++)sb[j]=shift[i*nl+local[j]];
            select_matrix(T+i,W,sb,dim,NULL,dim);flint_free(sb);
        }
        nmod_mat_t C; matrix_init(C,0,dim,ctx->mod.n);
        for(slong i=0;i<n && !C->r;i++) {
            for(slong start=0;start<nl;start+=256) {
                slong rows=FLINT_MIN(256,nl-start);
                nmod_mat_t block,product;
                nmod_mat_window_init(block,WL,start,0,start+rows,dim);
                matrix_init(product,rows,dim,ctx->mod.n);
                nmod_mat_mul(product,block,T+i);
                nmod_mat_window_clear(block);
                nmod_mat_clear(C); matrix_init(C,rows,dim,ctx->mod.n);
                for(slong a=0;a<rows;a++)for(slong b=0;b<dim;b++)
                    nmod_mat_entry(C,a,b)=nmod_sub(nmod_mat_entry(W,shift[i*nl+start+a],b),
                                                  nmod_mat_entry(product,a,b),ctx->mod);
                nmod_mat_clear(product); cl->matmul_work+=(uint64_t)rows*dim*dim;
                compact_rref(C,cl);
                if(C->r) break;
            }
        }
        nmod_mat_clear(WL);
        int stop=0;
        if (!C->r) {
            ulong *one=flint_malloc(dim*sizeof(ulong));
            for(slong j=0;j<dim;j++)one[j]=nmod_mat_entry(W,nc-1,j);
            certify_operators(result,T,one,B,&cl->basis,fs,m,ctx);flint_free(one);stop=1;
        } else {
            nmod_mat_t project,next;projection(project,C,&freecols);flint_free(freecols);
            matrix_init(next,nc,project->c,ctx->mod.n);nmod_mat_mul(next,W,project);
            cl->matmul_work+=(uint64_t)nc*dim*project->c;
            nmod_mat_swap(W,next);nmod_mat_clear(next);nmod_mat_clear(project);
        }
        for(slong i=0;i<n;i++)nmod_mat_clear(T+i);
        flint_free(T);flint_free(local);flint_free(B);nmod_mat_clear(C);
        if(stop)break;
    }
    flint_free(low);flint_free(shift);nmod_mat_clear(W);
}


/* Enumerate joint eigencharacters, restricting to an eigenspace at each
 * variable. This handles nonseparating coordinates, nilpotents, and fibers
 * with no F_p-rational points; no random separating form is required. */
static int extract_points(polynomial_solutions_t *sols, nmod_mat_struct *T,
                          const nmod_mat_t S, slong variable, ulong *point,
                          nmod_mpoly_struct *fs, slong m,
                          const nmod_mpoly_ctx_t ctx)
{
    slong n=sols->num_variables, dim=S->r, k=S->c;
    if (variable==n) {
        if(k!=1) return 0; /* Cyclicity implies a one-dimensional character. */
        for(slong i=0;i<m;i++)
            if(nmod_mpoly_evaluate_all_ui(fs+i,point,ctx)) return 0;
        if(sols->num_solution_sets>=T[0].r) return 0;
        slong row=sols->num_solution_sets++;
        sols->solution_sets[row]=calloc(n,sizeof(fq_nmod_t *));
        for(slong j=0;j<n;j++) {
            sols->solution_sets[row][j]=malloc(sizeof(fq_nmod_t));
            fq_nmod_init(sols->solution_sets[row][j][0],sols->ctx);
            fq_nmod_set_ui(sols->solution_sets[row][j][0],point[j],sols->ctx);
            sols->solutions_per_var[row*n+j]=1;
        }
        return 1;
    }
    nmod_mat_t ST,SB,inv,TS,rows,A;
    matrix_init(ST,k,dim,ctx->mod.n); nmod_mat_transpose(ST,S); compact_rref(ST,NULL);
    if(ST->r!=k) { nmod_mat_clear(ST); return 0; }
    slong *rr=pivots(ST); nmod_mat_clear(ST);
    select_matrix(SB,S,rr,k,NULL,k); matrix_init(inv,k,k,ctx->mod.n);
    int ok=nmod_mat_inv(inv,SB); nmod_mat_clear(SB);
    if(!ok) { flint_free(rr); nmod_mat_clear(inv); return 0; }
    matrix_init(TS,dim,k,ctx->mod.n); nmod_mat_mul(TS,T+variable,S);
    select_matrix(rows,TS,rr,k,NULL,k); flint_free(rr); nmod_mat_clear(TS);
    matrix_init(A,k,k,ctx->mod.n); nmod_mat_mul(A,inv,rows);
    nmod_mat_clear(inv); nmod_mat_clear(rows);
    nmod_poly_t mu; nmod_poly_init(mu,ctx->mod.n); nmod_mat_minpoly(mu,A);
    nmod_poly_factor_t roots; nmod_poly_factor_init(roots);
    nmod_poly_roots(roots,mu,0); nmod_poly_clear(mu);
    for(slong r=0;r<roots->num && ok;r++) {
        ulong a=nmod_neg(nmod_poly_get_coeff_ui(roots->p+r,0),ctx->mod);
        /* FLINT returns monic linear factors. */
        point[variable]=a;
        nmod_mat_t shifted,K,basis,next;
        nmod_mat_init_set(shifted,A);
        for(slong j=0;j<k;j++)
            nmod_mat_entry(shifted,j,j)=nmod_sub(nmod_mat_entry(shifted,j,j),a,ctx->mod);
        matrix_init(K,k,k,ctx->mod.n);
        slong nullity=nmod_mat_nullspace(K,shifted); nmod_mat_clear(shifted);
        if(!nullity) { nmod_mat_clear(K); ok=0; break; }
        select_matrix(basis,K,NULL,k,NULL,nullity); nmod_mat_clear(K);
        matrix_init(next,dim,nullity,ctx->mod.n); nmod_mat_mul(next,S,basis);
        nmod_mat_clear(basis);
        ok=extract_points(sols,T,next,variable+1,point,fs,m,ctx);
        nmod_mat_clear(next);
    }
    nmod_poly_factor_clear(roots); nmod_mat_clear(A); return ok;
}

/* Conservative scratch-space allowance: full RREF, projections, closure
 * batches, all operators and recursive eigencharacter restrictions. Input
 * polynomial storage and FLINT's allocator overhead are additional. */
static int degree_fits(slong n, slong degree, const slong *degrees,
                       slong m, slong memory_mb, slong *columns)
{
    slong N=monomial_count_bounded(n,degree,WORD_MAX/(slong)sizeof(ulong));
    if(!N) return -1;
    long double rows=0;
    for(slong i=0;i<m;i++) if(degrees[i]>=0)
        rows+=monomial_count_bounded(n,degree-degrees[i],N);
    long double cells=(30.0L+4.0L*n)*N*N+2.0L*rows*N;
    if(cells>(long double)WORD_MAX/sizeof(ulong)) return -1;
    if(memory_mb>0 && cells*sizeof(ulong)>(long double)memory_mb*1024*1024) return 0;
    *columns=N; return 1;
}

int solve_by_quotient_closure(char **polys, slong count,
                            variable_info_t *vars, slong n,
                            polynomial_solutions_t *sols,
                            slong max_degree, slong memory_mb)
{
    if(max_degree<0 || memory_mb<0) {
        sols->error_message=strdup("Quotient limits must be non-negative (0 = unlimited)");
        return 0;
    }
    if(fq_nmod_ctx_degree(sols->ctx)!=1) {
        return solve_by_quotient_closure_extension(polys, count, vars, n, sols,
                                                    max_degree, memory_mb);
    }
    if(n<=0 || count<n) {
        sols->error_message=strdup("Quotient closure requires at least as many equations as variables, and at least one variable");
        return 0;
    }
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_ctx_init(ctx,n,ORD_DEGLEX,fq_nmod_ctx_prime(sols->ctx));
    const char **names=flint_malloc(n*sizeof(char *));
    for(slong i=0;i<n;i++) names[i]=vars[i].name;
    nmod_mpoly_struct *fs=flint_malloc(count*sizeof(*fs));
    slong *degrees=flint_malloc(count*sizeof(slong)),start_degree=1;
    slong initialized=0,used_degree=0,columns=0;
    closure_t previous;
    int have_previous=0;
    int success=0;
    const char *stop_reason="degree limit reached";
    quotient_result result={0}; nmod_poly_init(result.eliminant,ctx->mod.n);
    for(slong i=0;i<count;i++) {
        nmod_mpoly_init(fs+i,ctx); initialized++;
        if(nmod_mpoly_set_str_pretty(fs+i,polys[i],names,ctx) ||
           !nmod_mpoly_total_degree_fits_si(fs+i,ctx)) {
            sols->error_message=strdup("Could not parse polynomial for quotient closure (or degree overflow)");
            goto done;
        }
        degrees[i]=fs[i].length ? nmod_mpoly_total_degree_si(fs+i,ctx) : -1;
        start_degree=FLINT_MAX(start_degree,degrees[i]);
    }
    for(slong i=0;i<count;i++) if(degrees[i]==0) {
        result.certified=1; result.dimension=0; nmod_poly_one(result.eliminant);
    }
    for(slong d=start_degree;!result.certified && (!max_degree || d<=max_degree);d++) {
        int fits=degree_fits(n,d,degrees,count,memory_mb,&columns);
        if(fits!=1) {
            stop_reason=fits<0 ? "matrix dimensions exceed addressable range" : "memory budget reached";
            break;
        }
        used_degree=d;
        closure_t cl; closure_init(&cl,n,d,ctx->mod.n);
        if (have_previous) {
            closure_seed(&cl,&previous);
            closure_clear(&previous); have_previous=0;
        }
        insert_original_rows(&cl,fs,count,ctx);
        quotient_close(&result,&cl,fs,count,ctx);
        if(!result.certified) {
            if (!cl.processed_ready) {
                nmod_mat_t processed;
                original_rows(processed,fs,count,ctx,&cl.basis,1,d-1);
                compact_rref(processed,NULL);
                cl.processed=pivots(processed); cl.nprocessed=processed->r;
                nmod_mat_clear(processed);
            }
            nmod_mat_t empty; matrix_init(empty,0,cl.basis.count,ctx->mod.n);
            slong oldrank=cl.E->r;
            closure_run(&cl,empty,ctx); nmod_mat_clear(empty);
            if(cl.E->r>oldrank) quotient_close(&result,&cl,fs,count,ctx);
        }
        if(result.certified) { closure_clear(&cl); break; }
        previous=cl; have_previous=1;
        if(d==WORD_MAX) { stop_reason="degree exceeds integer range"; break; }
    }
    if(!result.certified) {
        char error[384];
        snprintf(error,sizeof(error),
            "Quotient closure incomplete: %s (degree limit %ld, memory budget %ld MiB; 0 = unlimited; last attempted degree %ld). No conclusion about existence or dimension of solutions.",
            stop_reason,max_degree,memory_mb,used_degree);
        sols->error_message=strdup(error); goto done;
    }
    if(result.dimension) {
        sols->solution_sets=calloc(result.dimension,sizeof(fq_nmod_t **));
        sols->solutions_per_var=calloc(result.dimension*n,sizeof(slong));
        ulong *point=flint_calloc(n,sizeof(ulong));
        nmod_mat_t S; matrix_init(S,result.dimension,result.dimension,ctx->mod.n); nmod_mat_one(S);
        success=extract_points(sols,result.operators,S,0,point,fs,count,ctx);
        nmod_mat_clear(S); flint_free(point);
        if(!success) {
            sols->error_message=strdup("Quotient multiplication representation failed during rational-point extraction");
            goto done;
        }
    }
    success=1; sols->has_no_solutions=sols->num_solution_sets==0;
    sols->checked_solution_sets=sols->verified_solution_sets=sols->num_solution_sets;
    char summary[256];
    snprintf(summary,sizeof(summary),"Quotient algebra over F_%lu: certified dimension %ld, degree %ld, %ld monomials; %ld rational point(s)",
             ctx->mod.n,result.dimension,used_degree,columns,sols->num_solution_sets);
    sols->elimination_summary=strdup(summary);
    char *mu=nmod_poly_get_str_pretty(result.eliminant,names[n-1]);
    size_t len=strlen(mu)+strlen(names[n-1])+96;
    sols->resultant_steps=calloc(1,sizeof(char *));
    sols->resultant_steps[0]=malloc(len);
    snprintf(sols->resultant_steps[0],len,"Certified elimination polynomial in %s (minimal polynomial): %s",names[n-1],mu);
    sols->num_resultant_steps=sols->resultant_steps_cap=1; flint_free(mu);
 done:
    if(have_previous) closure_clear(&previous);
    if(result.operators) {
        for(slong i=0;i<n;i++) nmod_mat_clear(result.operators+i);
        flint_free(result.operators);
    }
    nmod_poly_clear(result.eliminant);
    for(slong i=0;i<initialized;i++) nmod_mpoly_clear(fs+i,ctx);
    flint_free(fs); flint_free(degrees); flint_free(names); nmod_mpoly_ctx_clear(ctx);
    return success;
}

/* Membership in a subspace of I is a sufficient certificate. Failure to
 * reduce to zero is deliberately NOT interpreted as nonmembership in I. */
static int sift_relation_contains(const closure_t *cl, const nmod_mpoly_t f,
                                  const nmod_mpoly_ctx_t ctx)
{
    if (!nmod_mpoly_total_degree_fits_si(f,ctx) ||
        nmod_mpoly_total_degree_si(f,ctx)>cl->basis.degree) return 0;
    slong nc=cl->basis.count, n=cl->basis.n;
    ulong *row=flint_calloc(nc,sizeof(ulong));
    ulong e[n];
    for (slong j=0;j<f->length;j++) {
        nmod_mpoly_get_term_exp_ui(e,f,j,ctx);
        slong col=basis_find(&cl->basis,e);
        assert(col>=0); row[col]=f->coeffs[j];
    }
    const relations_struct *E=cl->E;
    for (slong i=0;i<E->r;i++) {
        ulong a=row[E->p[i]];
        if (!a) continue;
        for (slong j=0;j<E->tail->c;j++)
            row[E->f[j]]=nmod_sub(row[E->f[j]],
                nmod_mul(a,nmod_mat_entry(E->tail,i,j),E->mod),E->mod);
    }
    int zero=1;
    for (slong j=0;j<E->tail->c;j++) if (row[E->f[j]]) { zero=0; break; }
    flint_free(row); return zero;
}

/* Multiples of certified relations are still in I, even when the tested
 * polynomial exceeds the closure degree. This is ordinary division by a
 * (not necessarily Groebner) set: only a zero remainder proves membership. */
static int sift_reduce_contains(const closure_t *cl, const nmod_mpoly_t f,
                                const nmod_mpoly_ctx_t ctx)
{
    const relations_struct *E=cl->E;
    slong n=cl->basis.n, nr=0;
    slong *selected=flint_malloc(E->r*sizeof(slong));
    /* Keep the minimal leading monomials to avoid thousands of redundant
     * monomial multiples in FLINT's division heap. */
    for (slong i=E->r;i-- > 0;) {
        const ulong *e=cl->basis.exps+E->p[i]*n;
        int redundant=0;
        for (slong j=0;j<nr && !redundant;j++) {
            const ulong *q=cl->basis.exps+E->p[selected[j]]*n;
            int divides=1;
            for (slong v=0;v<n;v++) if (q[v]>e[v]) { divides=0; break; }
            redundant=divides;
        }
        if (!redundant) selected[nr++]=i;
    }
    if (!nr) { flint_free(selected); return nmod_mpoly_is_zero(f,ctx); }
    nmod_mpoly_struct **basis=flint_malloc(nr*sizeof(*basis));
    nmod_mpoly_struct **quotients=flint_malloc(nr*sizeof(*quotients));
    for (slong j=0;j<nr;j++) {
        basis[j]=flint_malloc(sizeof(**basis)); nmod_mpoly_init(basis[j],ctx);
        quotients[j]=flint_malloc(sizeof(**quotients)); nmod_mpoly_init(quotients[j],ctx);
        slong i=selected[j];
        nmod_mpoly_push_term_ui_ui(basis[j],1,cl->basis.exps+E->p[i]*n,ctx);
        for (slong k=0;k<E->tail->c;k++) {
            ulong a=nmod_mat_entry(E->tail,i,k);
            if (a) nmod_mpoly_push_term_ui_ui(basis[j],a,cl->basis.exps+E->f[k]*n,ctx);
        }
        nmod_mpoly_sort_terms(basis[j],ctx);
    }
    nmod_mpoly_t remainder; nmod_mpoly_init(remainder,ctx);
    nmod_mpoly_divrem_ideal(quotients,remainder,f,basis,nr,ctx);
    int zero=nmod_mpoly_is_zero(remainder,ctx);
    nmod_mpoly_clear(remainder,ctx);
    for (slong j=0;j<nr;j++) {
        nmod_mpoly_clear(basis[j],ctx); flint_free(basis[j]);
        nmod_mpoly_clear(quotients[j],ctx); flint_free(quotients[j]);
    }
    flint_free(basis); flint_free(quotients); flint_free(selected);
    return zero;
}

int quotient_sift_resultant(nmod_mpoly_t sifted, nmod_mpoly_t removed,
                           const nmod_mpoly_t resultant,
                           nmod_mpoly_struct *polys, slong count,
                           const nmod_mpoly_ctx_t ctx,
                           slong max_degree,
                           slong *used_degree)
{
    *used_degree=0;
    nmod_mpoly_set(sifted,resultant,ctx);
    nmod_mpoly_one(removed,ctx);
    if (max_degree<=0 || count<=0 ||
        ctx->minfo->nvars<=0) return 0;
    if (nmod_mpoly_is_zero(resultant,ctx) || nmod_mpoly_is_ui(resultant,ctx)) return 1;
    if (!nmod_mpoly_total_degree_fits_si(resultant,ctx)) return 0;
    slong n=ctx->minfo->nvars, start=1;
    slong *degrees=flint_malloc(count*sizeof(slong));
    for (slong i=0;i<count;i++) {
        if (!nmod_mpoly_total_degree_fits_si(polys+i,ctx)) { flint_free(degrees); return 0; }
        degrees[i]=nmod_mpoly_total_degree_si(polys+i,ctx);
        start=FLINT_MAX(start,degrees[i]);
    }
    /* Check the degree bound and addressability before factorization. */
    slong columns;
    if (start>max_degree || degree_fits(n,start,degrees,count,0,&columns)!=1) {
        flint_free(degrees); return 1;
    }
    nmod_mpoly_factor_t factors;
    nmod_mpoly_factor_init(factors,ctx);
    if (!nmod_mpoly_factor(factors,resultant,ctx)) {
        nmod_mpoly_factor_clear(factors,ctx); flint_free(degrees); return 0;
    }
    slong nf=nmod_mpoly_factor_length(factors,ctx);
    slong *left=flint_malloc(nf*sizeof(slong));
    for (slong i=0;i<nf;i++) left[i]=nmod_mpoly_factor_get_exp_si(factors,i,ctx);
    nmod_mpoly_t candidate, factor;
    nmod_mpoly_init(candidate,ctx); nmod_mpoly_init(factor,ctx);
    closure_t previous;
    int have_previous=0;
    for (slong d=start;d<=max_degree;d++) {
        if (degree_fits(n,d,degrees,count,0,&columns)!=1) break;
        closure_t cl; closure_init(&cl,n,d,ctx->mod.n);
        if (have_previous) { closure_seed(&cl,&previous); closure_clear(&previous); }
        insert_original_rows(&cl,polys,count,ctx);
        nmod_mat_t empty; matrix_init(empty,0,cl.basis.count,ctx->mod.n);
        closure_run(&cl,empty,ctx); nmod_mat_clear(empty);
        *used_degree=d;
        for (slong i=0;i<nf;i++) {
            nmod_mpoly_factor_get_base(factor,factors,i,ctx);
            while (left[i]>0) {
                int divides=nmod_mpoly_divides(candidate,sifted,factor,ctx);
                assert(divides);
                if (!divides || (!sift_relation_contains(&cl,candidate,ctx) &&
                                 !sift_reduce_contains(&cl,candidate,ctx))) break;
                nmod_mpoly_swap(sifted,candidate,ctx);
                nmod_mpoly_mul(removed,removed,factor,ctx);
                left[i]--;
            }
        }
        previous=cl; have_previous=1;
        if (nmod_mpoly_is_ui(sifted,ctx) || d==WORD_MAX) break;
    }
    if (have_previous) closure_clear(&previous);
    nmod_mpoly_clear(candidate,ctx); nmod_mpoly_clear(factor,ctx);
    nmod_mpoly_factor_clear(factors,ctx);
    flint_free(left); flint_free(degrees);
    return 1;
}

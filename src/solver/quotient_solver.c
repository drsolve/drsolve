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
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>

typedef struct {
    slong n, degree, count, capacity, slots;
    ulong *exps;
    slong *degrees, *lookup;
} basis_t;

static uint64_t code(const ulong *e, slong n)
{
    uint64_t h = UINT64_C(1469598103934665603);
    for (slong i = 0; i < n; i++) { h ^= e[i]; h *= UINT64_C(1099511628211); }
    return h;
}

static slong hash_slot(uint64_t key, slong slots)
{
    key ^= key>>30; key *= UINT64_C(0xbf58476d1ce4e5b9);
    key ^= key>>27; key *= UINT64_C(0x94d049bb133111eb);
    return (slong)((key^(key>>31)) & (uint64_t)(slots-1));
}

static void basis_enumerate(basis_t *b, ulong *e, slong pos, slong left, slong degree)
{
    if (pos==b->n-1) {
        assert(b->count<b->capacity); e[pos]=left;
        memcpy(b->exps+b->count*b->n,e,b->n*sizeof(ulong));
        b->degrees[b->count++]=degree; return;
    }
    for (slong i=left;i>=0;i--) { e[pos]=i; basis_enumerate(b,e,pos+1,left-i,degree); }
}

static void basis_init(basis_t *b, slong n, slong degree)
{
    memset(b,0,sizeof(*b)); b->n=n; b->degree=degree;
    slong count=1;
    for (slong i=1;i<=degree;i++) count=count*(n+i)/i;
    b->capacity=count;
    b->exps=flint_malloc(count*n*sizeof(ulong)); b->degrees=flint_malloc(count*sizeof(slong));
    ulong e[n]; memset(e,0,sizeof(e));
    for (slong d=degree;d>=0;d--) basis_enumerate(b,e,0,d,d);
    assert(b->count==count); b->slots=1;
    while (b->slots<2*count) b->slots*=2;
    b->lookup=flint_malloc(b->slots*sizeof(slong));
    for (slong j=0;j<b->slots;j++) b->lookup[j]=-1;
    for (slong j=0;j<count;j++) {
        uint64_t k=code(b->exps+j*n,n); slong slot=hash_slot(k,b->slots);
        while (b->lookup[slot]>=0) slot=(slot+1)&(b->slots-1);
        b->lookup[slot]=j;
    }
}

static slong basis_find(const basis_t *b, const ulong *e)
{
    uint64_t k=code(e,b->n); slong slot=hash_slot(k,b->slots);
    while (b->lookup[slot]>=0) {
        if (!memcmp(b->exps+b->lookup[slot]*b->n,e,b->n*sizeof(ulong))) return b->lookup[slot];
        slot=(slot+1)&(b->slots-1);
    }
    return -1;
}

static void basis_clear(basis_t *b)
{
    flint_free(b->exps); flint_free(b->degrees); flint_free(b->lookup);
}

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

typedef struct {
    basis_t basis;
    nmod_mat_t E, processed;
    slong *shift, submitted, rref_calls;
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
    memset(cl,0,sizeof(*cl)); basis_init(&cl->basis,n,degree);
    slong nc=cl->basis.count;
    matrix_init(cl->E,0,nc,p); matrix_init(cl->processed,0,nc,p);
    cl->shift=flint_malloc(nc*n*sizeof(slong));
    for (slong j=0;j<nc;j++) for (slong v=0;v<n;v++) {
        ulong e[n]; memcpy(e,cl->basis.exps+j*n,n*sizeof(ulong)); e[v]++;
        cl->shift[j*n+v]=basis_find(&cl->basis,e);
    }
}

static void closure_clear(closure_t *cl)
{
    basis_clear(&cl->basis); nmod_mat_clear(cl->E); nmod_mat_clear(cl->processed); flint_free(cl->shift);
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
    slong *p=pivots(prev->E);
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
    nmod_mat_t U,V;
    matrix_init(U,rank,rank,cl->E->mod.n);
    matrix_init(V,rank,nf,cl->E->mod.n);
    for (slong j=0;j<nc;j++) if (source[j]>=0) {
        slong row=index[j], src=source[j], v=variable[j];
        for (slong k=p[src];k<pc;k++) {
            ulong a=nmod_mat_entry(prev->E,src,k);
            if (!a) continue;
            slong dest=v<0?k+offset:cl->shift[(k+offset)*n+v];
            assert(dest>=j);
            if (index[dest]>=0) nmod_mat_entry(U,row,index[dest])=a;
            else nmod_mat_entry(V,row,-index[dest]-1)=a;
        }
    }
    nmod_mat_solve_triu(V,U,V,1);
    cl->matmul_work+=(uint64_t)rank*rank*nf;
    nmod_mat_clear(cl->E); matrix_init(cl->E,rank,nc,U->mod.n);
    for (slong j=0;j<nc;j++) if (source[j]>=0) {
        slong row=index[j]; nmod_mat_entry(cl->E,row,j)=1;
        for (slong k=0;k<nf;k++) nmod_mat_entry(cl->E,row,cols[k])=nmod_mat_entry(V,row,k);
    }
    nmod_mat_clear(U); nmod_mat_clear(V);
    flint_free(p); flint_free(source); flint_free(variable); flint_free(index); flint_free(cols);
}

static void insert_rows(closure_t *cl, const nmod_mat_t M)
{
    slong nc=cl->basis.count, rank=cl->E->r; cl->submitted+=M->r;
    if (!M->r || rank==nc) return;
    if (!rank) {
        nmod_mat_clear(cl->E); nmod_mat_init_set(cl->E,M); compact_rref(cl->E,cl); return;
    }
    slong *p=pivots(cl->E), *freecols=flint_malloc(nc*sizeof(slong)), nf=0, pi=0;
    for (slong j=0;j<nc;j++) { if (pi<rank && p[pi]==j) pi++; else freecols[nf++]=j; }
    nmod_mat_t Ef,S;
    select_matrix(Ef,cl->E,NULL,rank,freecols,nf);
    matrix_init(S,M->r,nf,M->mod.n);
    /* Bound the pivot-column scratch space even for tall Macaulay batches. */
    for (slong start=0;start<M->r;start+=256) {
        slong rows=FLINT_MIN(256,M->r-start);
        nmod_mat_t Mp,block;
        matrix_init(Mp,rows,rank,M->mod.n);
        for (slong i=0;i<rows;i++) for (slong j=0;j<rank;j++)
            nmod_mat_entry(Mp,i,j)=nmod_mat_entry(M,start+i,p[j]);
        nmod_mat_window_init(block,S,start,0,start+rows,nf);
        nmod_mat_mul(block,Mp,Ef);
        for (slong i=0;i<rows;i++) for (slong j=0;j<nf;j++)
            nmod_mat_entry(block,i,j)=nmod_sub(nmod_mat_entry(M,start+i,freecols[j]),nmod_mat_entry(block,i,j),M->mod);
        nmod_mat_window_clear(block); nmod_mat_clear(Mp);
    }
    cl->matmul_work+=(uint64_t)M->r*rank*nf;
    nmod_mat_clear(Ef); compact_rref(S,cl);
    if (S->r) {
        slong *np=pivots(S);
        slong *rest=flint_malloc(nf*sizeof(slong)), *global=flint_malloc(nf*sizeof(slong)), nr=0;
        pi=0;
        for (slong j=0;j<nf;j++) {
            if (pi<S->r && np[pi]==j) pi++;
            else { rest[nr]=j; global[nr++]=freecols[j]; }
        }
        nmod_mat_t tail,En,product,updated,joined;
        select_matrix(tail,S,NULL,S->r,rest,nr);
        for (slong i=0;i<S->r;i++) np[i]=freecols[np[i]];
        select_matrix(En,cl->E,NULL,rank,np,S->r);
        matrix_init(product,rank,nr,M->mod.n); nmod_mat_mul(product,En,tail);
        select_matrix(updated,cl->E,NULL,rank,global,nr); nmod_mat_sub(updated,updated,product);
        cl->matmul_work+=(uint64_t)rank*S->r*nr;
        matrix_init(joined,rank+S->r,nc,M->mod.n);
        slong old=0,add=0;
        for (slong i=0;i<joined->r;i++) {
            int use_old=old<rank && (add==S->r || p[old]<np[add]);
            nmod_mat_entry(joined,i,use_old?p[old]:np[add])=1;
            for (slong j=0;j<nr;j++) nmod_mat_entry(joined,i,global[j])=use_old?nmod_mat_entry(updated,old,j):nmod_mat_entry(tail,add,j);
            if (use_old) old++; else add++;
        }
        nmod_mat_swap(cl->E,joined);
        nmod_mat_clear(joined); nmod_mat_clear(tail); nmod_mat_clear(En);
        nmod_mat_clear(product); nmod_mat_clear(updated); flint_free(np);
        flint_free(rest); flint_free(global);
    }
    nmod_mat_clear(S); flint_free(p); flint_free(freecols);
}

static void closure_run(closure_t *cl, const nmod_mat_t initial,
                        const nmod_mpoly_ctx_t ctx)
{
    slong n=cl->basis.n,nc=cl->basis.count;
    insert_rows(cl,initial);
    for (;;) {
        if (cl->E->r==nc) break;
        slong *p=pivots(cl->E);
        slong *lowidx=flint_malloc(cl->E->r*sizeof(slong)),nl=0;
        for (slong i=0;i<cl->E->r;i++) if (cl->basis.degrees[p[i]]<cl->basis.degree) lowidx[nl++]=i;
        nmod_mat_t low,delta;
        select_matrix(low,cl->E,lowidx,nl,NULL,nc); nmod_mat_init_set(delta,low);
        if (cl->processed->r) {
            slong *pp=pivots(cl->processed);
            nmod_mat_t Lp,product; select_matrix(Lp,low,NULL,nl,pp,cl->processed->r);
            matrix_init(product,nl,nc,ctx->mod.n); nmod_mat_mul(product,Lp,cl->processed);
            nmod_mat_sub(delta,delta,product);
            cl->matmul_work+=(uint64_t)nl*cl->processed->r*nc;
            nmod_mat_clear(Lp); nmod_mat_clear(product); flint_free(pp);
        }
        nmod_mat_swap(cl->processed,low); nmod_mat_clear(low); compact_rref(delta,cl);
        flint_free(p); flint_free(lowidx);
        if (!delta->r) { nmod_mat_clear(delta); break; }
        nmod_mat_t extra; matrix_init(extra,delta->r*n,nc,ctx->mod.n);
        for (slong i=0;i<delta->r;i++) for (slong j=0;j<nc;j++) if (nmod_mat_entry(delta,i,j))
            for (slong v=0;v<n;v++) {
                slong dest=cl->shift[j*n+v]; assert(dest>=0);
                nmod_mat_entry(extra,i*n+v,dest)=nmod_mat_entry(delta,i,j);
            }
        slong oldrank=cl->E->r; insert_rows(cl,extra);
        nmod_mat_clear(extra); nmod_mat_clear(delta);
        if (cl->E->r==oldrank) break;
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
    nmod_mat_t W;projection(W,cl->E,&freecols);flint_free(freecols);
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
        nmod_mat_t C,product,shifted;matrix_init(C,n*nl,dim,ctx->mod.n);matrix_init(product,nl,dim,ctx->mod.n);
        for(slong i=0;i<n;i++) {
            slong *sb=flint_malloc(dim*sizeof(slong));
            for(slong j=0;j<dim;j++)sb[j]=shift[i*nl+local[j]];
            select_matrix(T+i,W,sb,dim,NULL,dim);flint_free(sb);
            nmod_mat_mul(product,WL,T+i);select_matrix(shifted,W,shift+i*nl,nl,NULL,dim);
            for(slong a=0;a<nl;a++)for(slong b=0;b<dim;b++)
                nmod_mat_entry(C,i*nl+a,b)=nmod_sub(nmod_mat_entry(shifted,a,b),nmod_mat_entry(product,a,b),ctx->mod);
            nmod_mat_clear(shifted);cl->matmul_work+=(uint64_t)nl*dim*dim;
        }
        nmod_mat_clear(product);nmod_mat_clear(WL);compact_rref(C,cl);
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

static slong monomial_count_bounded(slong n, slong d, slong limit)
{
    slong count=1;
    for(slong i=1;i<=d;i++) {
        if(n>WORD_MAX-i || count>WORD_MAX/(n+i)) return 0;
        count=count*(n+i)/i;
        if(count>limit) return 0;
    }
    return count;
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
        sols->error_message=strdup("Quotient closure currently supports machine-word prime fields only");
        return 0;
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
        nmod_mat_t M; original_rows(M,fs,count,ctx,&cl.basis,1,d);
        insert_rows(&cl,M); nmod_mat_clear(M);
        quotient_close(&result,&cl,fs,count,ctx);
        if(!result.certified) {
            nmod_mat_clear(cl.processed);
            original_rows(cl.processed,fs,count,ctx,&cl.basis,1,d-1);
            compact_rref(cl.processed,NULL);
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

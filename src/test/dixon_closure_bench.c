/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Standalone exact experiment: original relation closure vs partial Dixon.
 * Not a production solver dispatch. Native FLINT arithmetic, one thread.
 * Input: prime n m on line 1, followed by m polynomials in x0,...,x(n-1).
 * Usage: dixon_closure_bench input method degree audit [repeats [matrix-dump]]
 * Methods: seed, hybrid3, hybrid4, macaulay, augment3, augment4.
 * Every method uses the same incremental scalar RREF closure. */
#include <flint/nmod_mpoly.h>
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 8
#define MAX_D 6
#define MAX_CELLS 32000000UL

static double now(void)
{
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9*t.tv_nsec;
}

typedef struct {
    slong n, degree, count, capacity, slots;
    ulong *exps;
    slong *degrees, *lookup;
    uint64_t *keys;
} basis_t;

static uint64_t code(const ulong *e, slong n)
{
    uint64_t x=0;
    for (slong i=0;i<n;i++) { assert(e[i]<16); x |= (uint64_t)e[i]<<(4*i); }
    return x;
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
    assert(count<65536); b->capacity=count;
    b->exps=flint_malloc(count*n*sizeof(ulong)); b->degrees=flint_malloc(count*sizeof(slong));
    ulong e[MAX_N]={0};
    for (slong d=degree;d>=0;d--) basis_enumerate(b,e,0,d,d);
    assert(b->count==count); b->slots=1;
    while (b->slots<2*count) b->slots*=2;
    b->keys=flint_calloc(b->slots,sizeof(uint64_t)); b->lookup=flint_malloc(b->slots*sizeof(slong));
    for (slong j=0;j<b->slots;j++) b->lookup[j]=-1;
    for (slong j=0;j<count;j++) {
        uint64_t k=code(b->exps+j*n,n); slong slot=hash_slot(k,b->slots);
        while (b->lookup[slot]>=0) slot=(slot+1)&(b->slots-1);
        b->lookup[slot]=j; b->keys[slot]=k;
    }
}

static slong basis_find(const basis_t *b, const ulong *e)
{
    uint64_t k=code(e,b->n); slong slot=hash_slot(k,b->slots);
    while (b->lookup[slot]>=0) {
        if (b->keys[slot]==k) return b->lookup[slot];
        slot=(slot+1)&(b->slots-1);
    }
    return -1;
}

static void basis_clear(basis_t *b)
{
    flint_free(b->exps); flint_free(b->degrees); flint_free(b->lookup); flint_free(b->keys);
}

static void matrix_init(nmod_mat_t M, slong r, slong c, ulong p)
{
    assert(r>=0 && c>=0 && (!c || (ulong)r<=MAX_CELLS/(ulong)c));
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
    slong *shift, history[128], nhistory, submitted, rref_calls;
    uint64_t rref_cells, matmul_work;
    int solved;
    ulong point[MAX_N];
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
        ulong e[MAX_N]; memcpy(e,cl->basis.exps+j*n,n*sizeof(ulong)); e[v]++;
        cl->shift[j*n+v]=basis_find(&cl->basis,e);
    }
}

static void closure_clear(closure_t *cl)
{
    basis_clear(&cl->basis); nmod_mat_clear(cl->E); nmod_mat_clear(cl->processed); flint_free(cl->shift);
}

static void insert_rows(closure_t *cl, const nmod_mat_t M)
{
    slong nc=cl->basis.count, rank=cl->E->r; cl->submitted+=M->r;
    if (!rank) {
        nmod_mat_clear(cl->E); nmod_mat_init_set(cl->E,M); compact_rref(cl->E,cl); return;
    }
    slong *p=pivots(cl->E), *freecols=flint_malloc(nc*sizeof(slong)), nf=0, pi=0;
    for (slong j=0;j<nc;j++) { if (pi<rank && p[pi]==j) pi++; else freecols[nf++]=j; }
    nmod_mat_t Mp,Ef,S,tmp;
    select_matrix(Mp,M,NULL,M->r,p,rank);
    select_matrix(Ef,cl->E,NULL,rank,freecols,nf);
    select_matrix(S,M,NULL,M->r,freecols,nf);
    matrix_init(tmp,M->r,nf,M->mod.n); nmod_mat_mul(tmp,Mp,Ef); nmod_mat_sub(S,S,tmp);
    cl->matmul_work+=(uint64_t)M->r*rank*nf;
    nmod_mat_clear(Mp); nmod_mat_clear(Ef); nmod_mat_clear(tmp); compact_rref(S,cl);
    if (S->r) {
        slong *np=pivots(S);
        for (slong i=0;i<S->r;i++) np[i]=freecols[np[i]];
        nmod_mat_t lift,En,product,updated,joined;
        matrix_init(lift,S->r,nc,M->mod.n);
        for (slong i=0;i<S->r;i++) for (slong j=0;j<nf;j++)
            nmod_mat_entry(lift,i,freecols[j])=nmod_mat_entry(S,i,j);
        select_matrix(En,cl->E,NULL,rank,np,S->r);
        matrix_init(product,rank,nc,M->mod.n); nmod_mat_mul(product,En,lift);
        matrix_init(updated,rank,nc,M->mod.n); nmod_mat_sub(updated,cl->E,product);
        cl->matmul_work+=(uint64_t)rank*S->r*nc;
        matrix_init(joined,rank+S->r,nc,M->mod.n);
        slong old=0,add=0;
        for (slong i=0;i<joined->r;i++) {
            int use_old=old<rank && (add==S->r || p[old]<np[add]);
            for (slong j=0;j<nc;j++) nmod_mat_entry(joined,i,j)=use_old?nmod_mat_entry(updated,old,j):nmod_mat_entry(lift,add,j);
            if (use_old) old++; else add++;
        }
        nmod_mat_swap(cl->E,joined);
        nmod_mat_clear(joined); nmod_mat_clear(lift); nmod_mat_clear(En);
        nmod_mat_clear(product); nmod_mat_clear(updated); flint_free(np);
    }
    nmod_mat_clear(S); flint_free(p); flint_free(freecols);
}

static void closure_run(closure_t *cl, const nmod_mat_t initial,
                        nmod_mpoly_struct *fs, slong m, const nmod_mpoly_ctx_t ctx)
{
    slong n=cl->basis.n,nc=cl->basis.count;
    insert_rows(cl,initial);
    for (;;) {
        assert(cl->nhistory<128); cl->history[cl->nhistory++]=cl->E->r;
        if (cl->E->r==nc) { cl->solved=-1; break; }
        slong *p=pivots(cl->E);
        if (cl->E->r==nc-1 && p[nc-2]==nc-2) {
            for (slong v=0;v<n;v++) {
                ulong e[MAX_N]={0}; e[v]=1; slong j=basis_find(&cl->basis,e);
                assert(p[j]==j); cl->point[v]=nmod_neg(nmod_mat_entry(cl->E,j,nc-1),cl->E->mod);
            }
            for (slong i=0;i<m;i++) assert(nmod_mpoly_evaluate_all_ui(fs+i,cl->point,ctx)==0);
            cl->solved=1; flint_free(p); break;
        }
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
                ulong e[MAX_N]; nmod_mpoly_get_term_exp_ui(e,fs+i,a,ctx);
                for (slong v=0;v<basis->n;v++) e[v]+=multipliers.exps[j*basis->n+v];
                slong col=basis_find(basis,e); assert(col>=0);
                nmod_mat_entry(M,row,col)=fs[i].coeffs[a];
            }
            row++;
        }
    }
    assert(row==count); basis_clear(&multipliers);
}

typedef struct { slong terms, rows, peak; uint64_t products; double seconds; } dixon_stats_t;

static void filter_auxiliary(nmod_mpoly_t dest, const nmod_mpoly_t src,
                             slong n, slong lo, slong hi, const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_zero(dest,ctx);
    for (slong j=0;j<src->length;j++) {
        ulong e[2*MAX_N]; nmod_mpoly_get_term_exp_ui(e,src,j,ctx); slong d=0;
        for (slong v=n;v<2*n-1;v++) d+=e[v];
        if (lo<=d && d<=hi) nmod_mpoly_push_term_ui_ui(dest,src->coeffs[j],e,ctx);
    }
    /* Filtering a canonical term sequence preserves its order. */
    assert(nmod_mpoly_is_canonical(dest,ctx));
}

static void construct_dixon(nmod_mpoly_t out, nmod_mpoly_struct *fs, slong n,
                            const nmod_mpoly_ctx_t pctx, const nmod_mpoly_ctx_t qctx,
                            slong rd, int full, dixon_stats_t *stats)
{
    slong size=1L<<n,k=n-1;
    nmod_mpoly_struct *rows=flint_malloc(n*n*sizeof(*rows)), *dp=flint_malloc(size*sizeof(*dp));
    for (slong i=0;i<n*n;i++) nmod_mpoly_init(rows+i,qctx);
    for (slong i=0;i<size;i++) nmod_mpoly_init(dp+i,qctx);
    for (slong col=0;col<n;col++) for (slong term=0;term<fs[col].length;term++) {
        ulong e[MAX_N],v[2*MAX_N]={0}; nmod_mpoly_get_term_exp_ui(e,fs+col,term,pctx);
        memcpy(v,e,n*sizeof(ulong)); nmod_mpoly_push_term_ui_ui(rows+k*n+col,fs[col].coeffs[term],v,qctx);
        for (slong r=0;r<k;r++) for (ulong h=0;h<e[r];h++) {
            memset(v,0,sizeof(v));
            for (slong j=0;j<n;j++) v[j<r?n+j:j]=e[j];
            v[r]=e[r]-1-h; v[n+r]=h;
            nmod_mpoly_push_term_ui_ui(rows+r*n+col,fs[col].coeffs[term],v,qctx);
        }
    }
    for (slong i=0;i<n*n;i++) { nmod_mpoly_sort_terms(rows+i,qctx); nmod_mpoly_combine_like_terms(rows+i,qctx); }
    nmod_mpoly_one(dp,qctx);
    nmod_mpoly_t value,product,filtered;
    nmod_mpoly_init(value,qctx); nmod_mpoly_init(product,qctx); nmod_mpoly_init(filtered,qctx);
    slong lo=FLINT_MAX(0,n+1-rd),hi=k-1;
    for (slong mask=1;mask<size;mask++) {
        slong r=__builtin_popcountl((ulong)mask)-1,j=0;
        nmod_mpoly_zero(value,qctx);
        for (slong col=0;col<n;col++) if (mask&(1L<<col)) {
            nmod_mpoly_struct *a=rows+r*n+col,*b=dp+(mask^(1L<<col));
            stats->products+=(uint64_t)a->length*b->length;
            nmod_mpoly_mul(product,a,b,qctx);
            if ((r+j)&1) nmod_mpoly_sub(value,value,product,qctx);
            else nmod_mpoly_add(value,value,product,qctx);
            j++;
        }
        if (full) nmod_mpoly_swap(dp+mask,value,qctx);
        else {
            filter_auxiliary(filtered,value,n,lo-FLINT_MAX(0,k-r-1),hi,qctx);
            nmod_mpoly_swap(dp+mask,filtered,qctx);
        }
        stats->peak=FLINT_MAX(stats->peak,dp[mask].length);
    }
    nmod_mpoly_set(out,dp+size-1,qctx); stats->terms=out->length;
    for (slong i=0;i<n*n;i++) nmod_mpoly_clear(rows+i,qctx);
    for (slong i=0;i<size;i++) nmod_mpoly_clear(dp+i,qctx);
    flint_free(rows); flint_free(dp);
    nmod_mpoly_clear(value,qctx); nmod_mpoly_clear(product,qctx); nmod_mpoly_clear(filtered,qctx);
}

static void dixon_rows(nmod_mat_t G, const nmod_mpoly_t delta,
                       const nmod_mpoly_ctx_t ctx, const basis_t *basis, dixon_stats_t *stats)
{
    slong n=basis->n; basis_t aux; basis_init(&aux,n-1,n-1);
    slong *map=flint_malloc(aux.count*sizeof(slong));
    for (slong i=0;i<aux.count;i++) map[i]=-1;
    for (slong t=0;t<delta->length;t++) {
        ulong e[2*MAX_N]; nmod_mpoly_get_term_exp_ui(e,delta,t,ctx);
        slong row=basis_find(&aux,e+n); assert(row>=0); map[row]=0;
    }
    slong count=0;
    for (slong i=0;i<aux.count;i++) if (map[i]>=0) map[i]=count++;
    matrix_init(G,count,basis->count,ctx->mod.n);
    for (slong t=0;t<delta->length;t++) {
        ulong e[2*MAX_N]; nmod_mpoly_get_term_exp_ui(e,delta,t,ctx);
        slong row=map[basis_find(&aux,e+n)], col=basis_find(basis,e); assert(col>=0 && row>=0);
        nmod_mat_entry(G,row,col)=delta->coeffs[t];
    }
    stats->rows=count; basis_clear(&aux); flint_free(map);
}

static uint64_t fingerprint(const nmod_mat_t M)
{
    uint64_t h=UINT64_C(14695981039346656037);
    for (slong i=0;i<M->r;i++) for (slong j=0;j<M->c;j++) h=(h^nmod_mat_entry(M,i,j))*UINT64_C(1099511628211);
    return h;
}

static void audit_dixon(const nmod_mat_t G, const nmod_mpoly_t delta, slong rd,
                        nmod_mpoly_struct *fs, slong m, const nmod_mpoly_ctx_t pctx,
                        const nmod_mpoly_ctx_t qctx, const basis_t *basis)
{
    nmod_mat_t M,Gp,product;
    original_rows(M,fs,m,pctx,basis,1,basis->degree); compact_rref(M,NULL); slong *p=pivots(M);
    select_matrix(Gp,G,NULL,G->r,p,M->r); matrix_init(product,G->r,G->c,pctx->mod.n);
    nmod_mat_mul(product,Gp,M); assert(nmod_mat_equal(product,G));
    nmod_mat_clear(M); nmod_mat_clear(Gp); nmod_mat_clear(product); flint_free(p);
    if (basis->n<=5 && rd) {
        nmod_mpoly_t full,filtered; nmod_mpoly_init(full,qctx); nmod_mpoly_init(filtered,qctx);
        dixon_stats_t stats={0}; construct_dixon(full,fs,basis->n,pctx,qctx,rd,1,&stats);
        filter_auxiliary(filtered,full,basis->n,FLINT_MAX(0,basis->n+1-rd),basis->n-2,qctx);
        assert(nmod_mpoly_equal(delta,filtered,qctx));
        nmod_mpoly_clear(full,qctx); nmod_mpoly_clear(filtered,qctx);
    }
}

int main(int argc,char **argv)
{
    if (argc<5 || argc>7) { fprintf(stderr,"usage: %s input method degree audit [repeats [matrix-dump]]\n",argv[0]); return 2; }
    const char *method=argv[2]; slong degree=atol(argv[3]),rd=0;
    int all=!strcmp(method,"macaulay") || !strncmp(method,"augment",7), audit=atoi(argv[4]);
    if (!strncmp(method,"hybrid",6) || !strncmp(method,"augment",7)) rd=atol(method+strlen(method)-1);
    if (strcmp(method,"seed") && strcmp(method,"macaulay") && rd!=3 && rd!=4) return 2;
    int repeats=argc>=6?atoi(argv[5]):1;
    FILE *fp=fopen(argv[1],"r"); if (!fp) return 2;
    ulong prime; slong n,m;
    if (fscanf(fp,"%lu %ld %ld\n",&prime,&n,&m)!=3 || !n_is_prime(prime) || n<2 || n>MAX_N || m<n || m>256 || degree<2 || degree>MAX_D || rd>degree || repeats<1) return 2;
    flint_set_num_threads(1);
    nmod_mpoly_ctx_t ctx,qctx; nmod_mpoly_ctx_init(ctx,n,ORD_DEGLEX,prime); nmod_mpoly_ctx_init(qctx,2*n-1,ORD_DEGLEX,prime);
    char *names[MAX_N];
    for (slong i=0;i<n;i++) { names[i]=flint_malloc(24); snprintf(names[i],24,"x%ld",i); }
    nmod_mpoly_struct *fs=flint_malloc(m*sizeof(*fs)); char *line=NULL; size_t capacity=0;
    for (slong i=0;i<m;i++) {
        nmod_mpoly_init(fs+i,ctx);
        if (getline(&line,&capacity,fp)<0) { fprintf(stderr,"missing polynomial %ld\n",i); return 2; }
        line[strcspn(line,"\r\n")]='\0';
        if (nmod_mpoly_set_str_pretty(fs+i,line,(const char**)names,ctx) || nmod_mpoly_total_degree_si(fs+i,ctx)>2) {
            fprintf(stderr,"invalid quadratic polynomial %ld\n",i); return 2;
        }
    }
    free(line); fclose(fp);
    for (int repeat=0;repeat<repeats;repeat++) {
        double start=now(), cpu_start=(double)clock()/CLOCKS_PER_SEC;
        closure_t cl; closure_init(&cl,n,degree,prime);
        nmod_mpoly_t delta; nmod_mpoly_init(delta,qctx); dixon_stats_t ds={0};
        nmod_mat_t G,M,initial;
        double dt=now(), dcpu=(double)clock()/CLOCKS_PER_SEC;
        if (rd) { construct_dixon(delta,fs,n,ctx,qctx,rd,0,&ds); dixon_rows(G,delta,qctx,&cl.basis,&ds); }
        else matrix_init(G,0,cl.basis.count,prime);
        ds.seconds=now()-dt; double dixon_cpu=(double)clock()/CLOCKS_PER_SEC-dcpu;
        original_rows(M,fs,m,ctx,&cl.basis,all,degree);
        if (all) {
            nmod_mat_clear(cl.processed);
            original_rows(cl.processed,fs,m,ctx,&cl.basis,1,degree-1);
            compact_rref(cl.processed,&cl);
        }
        matrix_init(initial,M->r+G->r,M->c,prime); nmod_mat_concat_vertical(initial,M,G);
        double built=now(), cpu_built=(double)clock()/CLOCKS_PER_SEC;
        closure_run(&cl,initial,fs,m,ctx);
        double done=now(), cpu_done=(double)clock()/CLOCKS_PER_SEC;
        if (audit) audit_dixon(G,delta,rd,fs,m,ctx,qctx,&cl.basis);
        if (argc==7) {
            FILE *dump=fopen(argv[6],"w"); if (!dump) return 2;
            fprintf(dump,"%ld %ld\n",cl.E->r,cl.E->c);
            for (slong i=0;i<cl.E->r;i++) {
                for (slong j=0;j<cl.E->c;j++) fprintf(dump,"%s%lu",j?" ":"",nmod_mat_entry(cl.E,i,j));
                fputc('\n',dump);
            }
            fclose(dump);
        }
        printf("{\"method\":\"%s\",\"n\":%ld,\"m\":%ld,\"prime\":%lu,\"degree\":%ld,\"repeat\":%d,"
               "\"seconds\":%.9f,\"setup_seconds\":%.9f,\"closure_seconds\":%.9f,\"dixon_seconds\":%.9f,"
               "\"dixon_rows\":%ld,\"dixon_terms\":%ld,\"candidate_term_products\":%llu,\"peak_polynomial_terms\":%ld,"
               "\"initial_rows\":%ld,\"columns\":%ld,\"submitted_rows\":%ld,\"rref_cells\":%llu,\"rref_calls\":%ld,"
               "\"dense_matmul_multiply_adds\":%llu,\"nullity\":%ld,\"audit\":%s,\"rowspace_hash\":\"%016llx\",\"rank_history\":[",
               method,n,m,prime,degree,repeat,done-start,built-start,done-built,ds.seconds,
               ds.rows,ds.terms,(unsigned long long)ds.products,ds.peak,initial->r,cl.basis.count,
               cl.submitted,(unsigned long long)cl.rref_cells,cl.rref_calls,(unsigned long long)cl.matmul_work,
               cl.basis.count-cl.E->r,audit?"true":"false",(unsigned long long)fingerprint(cl.E));
        for (slong i=0;i<cl.nhistory;i++) printf("%s%ld",i?",":"",cl.history[i]);
        printf("],\"recovered\":");
        if (cl.solved==1) {
            putchar('['); for (slong v=0;v<n;v++) printf("%s%lu",v?",":"",cl.point[v]); putchar(']');
        } else printf("null");
        printf(",\"inconsistent\":%s,\"cpu_seconds\":%.9f,\"cpu_setup_seconds\":%.9f,"
               "\"cpu_closure_seconds\":%.9f,\"cpu_dixon_seconds\":%.9f}\n",
               cl.solved==-1?"true":"false",cpu_done-cpu_start,cpu_built-cpu_start,cpu_done-cpu_built,dixon_cpu);
        fflush(stdout);
        nmod_mat_clear(G); nmod_mat_clear(M); nmod_mat_clear(initial); nmod_mpoly_clear(delta,qctx); closure_clear(&cl);
    }
    for (slong i=0;i<m;i++) nmod_mpoly_clear(fs+i,ctx);
    flint_free(fs); for (slong i=0;i<n;i++) flint_free(names[i]);
    nmod_mpoly_ctx_clear(ctx); nmod_mpoly_ctx_clear(qctx); flint_cleanup_master(); return 0;
}

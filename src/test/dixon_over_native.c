/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Native experiments; reuse the audited scalar basis/RREF and Dixon kernel.
 * This standalone target does not change production solver dispatch. */
#define main closure_baseline_main
#include "dixon_closure_bench.c"
#undef main
#include <flint/nmod_poly.h>

typedef struct {
    int certified;
    slong dimension, dims[128], ndims;
    ulong point[MAX_N];
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
            ulong e[MAX_N]; nmod_mpoly_get_term_exp_ui(e,fs+i,j,ctx);
            act_monomial(v,one,e,T,n,dim);
            for (slong k=0;k<dim;k++) sum[k]=nmod_add(sum[k],nmod_mul(fs[i].coeffs[j],v[k],ctx->mod),ctx->mod);
        }
        for (slong j=0;j<dim;j++) if (sum[j]) ok=0;
    }
    if (ok) {
        result->certified=1;result->dimension=dim;
        if (dim==1) {
            for (slong i=0;i<n;i++) result->point[i]=nmod_mat_entry(T+i,0,0);
            for (slong i=0;i<m;i++) assert(!nmod_mpoly_evaluate_all_ui(fs+i,result->point,ctx));
        }
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
        slong dim=W->c;assert(result->ndims<128);result->dims[result->ndims++]=dim;
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
        assert(nmod_mat_inv(inv,WB));matrix_init(newW,nc,dim,ctx->mod.n);nmod_mat_mul(newW,W,inv);
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

static void make_column(nmod_mpoly_struct *col,const nmod_mpoly_t f,slong n,
                        const nmod_mpoly_ctx_t ctx,const nmod_mpoly_ctx_t qctx)
{
    for(slong r=0;r<n;r++)nmod_mpoly_zero(col+r,qctx);
    for(slong term=0;term<f->length;term++) {
        ulong e[MAX_N],v[2*MAX_N]={0};nmod_mpoly_get_term_exp_ui(e,f,term,ctx);
        memcpy(v,e,n*sizeof(ulong));nmod_mpoly_push_term_ui_ui(col+n-1,f->coeffs[term],v,qctx);
        for(slong r=0;r<n-1;r++)for(ulong h=0;h<e[r];h++) {
            memset(v,0,sizeof(v));for(slong j=0;j<n;j++)v[j<r?n+j:j]=e[j];
            v[r]=e[r]-1-h;v[n+r]=h;nmod_mpoly_push_term_ui_ui(col+r,f->coeffs[term],v,qctx);
        }
    }
    for(slong r=0;r<n;r++){nmod_mpoly_sort_terms(col+r,qctx);nmod_mpoly_combine_like_terms(col+r,qctx);}
}

static nmod_mpoly_struct *cofactor_cache(nmod_mpoly_struct *fs,slong n,slong degree,
    const nmod_mpoly_ctx_t ctx,const nmod_mpoly_ctx_t qctx,dixon_stats_t *stats)
{
    slong size=1L<<n,lo=FLINT_MAX(0,n+1-degree),hi=n-2;
    nmod_mpoly_struct *columns=flint_malloc(n*(n-1)*sizeof(*columns)),*dp=flint_malloc(size*sizeof(*dp));
    for(slong j=0;j<n*(n-1);j++)nmod_mpoly_init(columns+j,qctx);
    for(slong j=0;j<n-1;j++)make_column(columns+j*n,fs+j,n,ctx,qctx);
    for(slong j=0;j<size;j++)nmod_mpoly_init(dp+j,qctx);
    nmod_mpoly_one(dp,qctx);nmod_mpoly_t value,tmp;
    nmod_mpoly_init(value,qctx);nmod_mpoly_init(tmp,qctx);
    for(slong mask=1;mask<size-1;mask++) {
        slong col=__builtin_popcountl((ulong)mask)-1,j=0;nmod_mpoly_zero(value,qctx);
        for(slong r=0;r<n;r++)if(mask&(1L<<r)) {
            nmod_mpoly_struct *a=columns+col*n+r,*b=dp+(mask^(1L<<r));
            stats->products+=(uint64_t)a->length*b->length;nmod_mpoly_mul(tmp,a,b,qctx);
            if((j+col)&1)nmod_mpoly_sub(value,value,tmp,qctx);else nmod_mpoly_add(value,value,tmp,qctx);
            j++;
        }
        filter_auxiliary(dp+mask,value,n,lo-(n-col-1),hi,qctx);
        stats->peak=FLINT_MAX(stats->peak,dp[mask].length);
    }
    nmod_mpoly_struct *cof=flint_malloc(n*sizeof(*cof));
    for(slong r=0;r<n;r++) {
        nmod_mpoly_init(cof+r,qctx);nmod_mpoly_swap(cof+r,dp+((size-1)^(1L<<r)),qctx);
        if((r+n-1)&1)nmod_mpoly_neg(cof+r,cof+r,qctx);
    }
    for(slong j=0;j<n*(n-1);j++)nmod_mpoly_clear(columns+j,qctx);
    for(slong j=0;j<size;j++)nmod_mpoly_clear(dp+j,qctx);
    flint_free(columns);flint_free(dp);nmod_mpoly_clear(value,qctx);nmod_mpoly_clear(tmp,qctx);return cof;
}

static void cached_dixon(nmod_mpoly_t out,nmod_mpoly_struct *cof,const nmod_mpoly_t g,
    slong n,slong degree,const nmod_mpoly_ctx_t ctx,const nmod_mpoly_ctx_t qctx,dixon_stats_t *stats)
{
    nmod_mpoly_struct *col=flint_malloc(n*sizeof(*col));
    for(slong r=0;r<n;r++)nmod_mpoly_init(col+r,qctx);
    make_column(col,g,n,ctx,qctx);nmod_mpoly_t tmp,sum;
    nmod_mpoly_init(tmp,qctx);nmod_mpoly_init(sum,qctx);
    for(slong r=0;r<n;r++) {
        stats->products+=(uint64_t)col[r].length*cof[r].length;
        nmod_mpoly_mul(tmp,col+r,cof+r,qctx);nmod_mpoly_add(sum,sum,tmp,qctx);
    }
    filter_auxiliary(out,sum,n,FLINT_MAX(0,n+1-degree),n-2,qctx);stats->terms+=out->length;
    for(slong r=0;r<n;r++)nmod_mpoly_clear(col+r,qctx);
    flint_free(col);nmod_mpoly_clear(tmp,qctx);nmod_mpoly_clear(sum,qctx);
}

static nmod_mpoly_struct *group_cofactors(nmod_mpoly_struct *cof,slong n,basis_t *aux,
    const nmod_mpoly_ctx_t ctx,const nmod_mpoly_ctx_t qctx)
{
    basis_init(aux,n-1,n-1);slong count=n*aux->count;
    nmod_mpoly_struct *groups=flint_malloc(count*sizeof(*groups));
    for(slong i=0;i<count;i++)nmod_mpoly_init(groups+i,ctx);
    for(slong r=0;r<n;r++)for(slong j=0;j<cof[r].length;j++) {
        ulong e[2*MAX_N];nmod_mpoly_get_term_exp_ui(e,cof+r,j,qctx);
        slong k=basis_find(aux,e+n);assert(k>=0);
        nmod_mpoly_push_term_ui_ui(groups+r*aux->count+k,cof[r].coeffs[j],e,ctx);
    }
    for(slong i=0;i<count;i++){nmod_mpoly_sort_terms(groups+i,ctx);nmod_mpoly_combine_like_terms(groups+i,ctx);}
    return groups;
}

static void cached_dixon_grouped(nmod_mpoly_t out,nmod_mpoly_struct *cof,const basis_t *aux,
    const nmod_mpoly_t g,slong n,slong degree,const nmod_mpoly_ctx_t ctx,
    const nmod_mpoly_ctx_t qctx,dixon_stats_t *stats)
{
    nmod_mpoly_struct *col=flint_malloc(n*sizeof(*col)),*entries=flint_malloc(n*n*sizeof(*entries));
    for(slong r=0;r<n;r++)nmod_mpoly_init(col+r,qctx);
    for(slong r=0;r<n*n;r++)nmod_mpoly_init(entries+r,ctx);
    make_column(col,g,n,ctx,qctx);
    for(slong r=0;r<n;r++)for(slong j=0;j<col[r].length;j++) {
        ulong e[2*MAX_N];nmod_mpoly_get_term_exp_ui(e,col+r,j,qctx);slong which=0;
        for(slong k=0;k<n-1;k++)if(e[n+k]){assert(!which && e[n+k]==1);which=k+1;}
        nmod_mpoly_push_term_ui_ui(entries+r*n+which,col[r].coeffs[j],e,ctx);
    }
    for(slong r=0;r<n*n;r++){nmod_mpoly_sort_terms(entries+r,ctx);nmod_mpoly_combine_like_terms(entries+r,ctx);}
    nmod_mpoly_t value,product;nmod_mpoly_init(value,ctx);nmod_mpoly_init(product,ctx);nmod_mpoly_zero(out,qctx);
    slong lo=FLINT_MAX(0,n+1-degree),hi=n-2;
    for(slong b=0;b<aux->count;b++)if(lo<=aux->degrees[b] && aux->degrees[b]<=hi) {
        nmod_mpoly_zero(value,ctx);
        for(slong r=0;r<n;r++)for(slong j=0;j<n;j++)if(entries[r*n+j].length) {
            ulong beta[MAX_N];memcpy(beta,aux->exps+b*(n-1),(n-1)*sizeof(ulong));
            if(j && !beta[j-1])continue;
            if(j)beta[j-1]--;
            slong k=basis_find(aux,beta);assert(k>=0);
            nmod_mpoly_struct *a=entries+r*n+j,*h=cof+r*aux->count+k;
            if(!h->length)continue;
            stats->products+=(uint64_t)a->length*h->length;nmod_mpoly_mul(product,a,h,ctx);
            nmod_mpoly_add(value,value,product,ctx);
        }
        for(slong j=0;j<value->length;j++) {
            ulong e[2*MAX_N]={0};nmod_mpoly_get_term_exp_ui(e,value,j,ctx);
            memcpy(e+n,aux->exps+b*(n-1),(n-1)*sizeof(ulong));
            nmod_mpoly_push_term_ui_ui(out,value->coeffs[j],e,qctx);
        }
    }
    nmod_mpoly_sort_terms(out,qctx);nmod_mpoly_combine_like_terms(out,qctx);stats->terms+=out->length;
    for(slong r=0;r<n;r++)nmod_mpoly_clear(col+r,qctx);
    for(slong r=0;r<n*n;r++)nmod_mpoly_clear(entries+r,ctx);
    flint_free(col);flint_free(entries);nmod_mpoly_clear(value,ctx);nmod_mpoly_clear(product,ctx);
}

int main(int argc,char **argv)
{
    if(argc<5 || argc>6){fprintf(stderr,"usage: %s input quotient|adaptive|shared degree audit [repeats]\n",argv[0]);return 2;}
    const char *method=argv[2];int shared=!strcmp(method,"shared") || !strcmp(method,"shared_full"),adaptive=!strcmp(method,"adaptive"),oracle=!strcmp(method,"oracle_seed");
    if(!shared && !adaptive && !oracle && strcmp(method,"quotient"))return closure_baseline_main(argc,argv);
    slong degree=atol(argv[3]);int audit=atoi(argv[4]),repeats=argc==6?atoi(argv[5]):1;
    FILE *fp=fopen(argv[1],"r");if(!fp)return 2;
    ulong prime;slong n,m;
    if(fscanf(fp,"%lu %ld %ld\n",&prime,&n,&m)!=3 || !n_is_prime(prime) || n<2 || n>MAX_N || m<n || m>256 || degree<2 || degree>MAX_D || repeats<1)return 2;
    flint_set_num_threads(1);nmod_mpoly_ctx_t ctx,qctx;
    nmod_mpoly_ctx_init(ctx,n,ORD_DEGLEX,prime);nmod_mpoly_ctx_init(qctx,2*n-1,ORD_DEGLEX,prime);
    char *names[MAX_N];for(slong i=0;i<n;i++){names[i]=flint_malloc(24);snprintf(names[i],24,"x%ld",i);}
    nmod_mpoly_struct *fs=flint_malloc(m*sizeof(*fs));char *line=NULL;size_t cap=0;
    for(slong i=0;i<m;i++) {
        nmod_mpoly_init(fs+i,ctx);if(getline(&line,&cap,fp)<0)return 2;line[strcspn(line,"\r\n")]='\0';
        if(nmod_mpoly_set_str_pretty(fs+i,line,(const char**)names,ctx) || nmod_mpoly_total_degree_si(fs+i,ctx)>(oracle?degree:2))return 2;
    }
    free(line);fclose(fp);
    for(int rep=0;rep<repeats;rep++) {
        double start=now(),cpu=(double)clock()/CLOCKS_PER_SEC;
        if(shared) {
            dixon_stats_t stats={0};nmod_mpoly_struct *cof=cofactor_cache(fs,n,degree,ctx,qctx,&stats);
            basis_t aux;nmod_mpoly_struct *groups=NULL;
            if(!strcmp(method,"shared"))groups=group_cofactors(cof,n,&aux,ctx,qctx);
            double cache_seconds=now()-start,cache_cpu=(double)clock()/CLOCKS_PER_SEC-cpu;
            nmod_mpoly_struct *out=flint_malloc((m-n+1)*sizeof(*out));
            for(slong i=n-1;i<m;i++){
                nmod_mpoly_init(out+i-n+1,qctx);
                if(groups)cached_dixon_grouped(out+i-n+1,groups,&aux,fs+i,n,degree,ctx,qctx,&stats);
                else cached_dixon(out+i-n+1,cof,fs+i,n,degree,ctx,qctx,&stats);
            }
            double elapsed=now()-start,cpus=(double)clock()/CLOCKS_PER_SEC-cpu;
            for(slong r=0;r<n;r++)nmod_mpoly_clear(cof+r,qctx);
            flint_free(cof);
            if(groups){for(slong j=0;j<n*aux.count;j++)nmod_mpoly_clear(groups+j,ctx);flint_free(groups);basis_clear(&aux);}
            dixon_stats_t reference={0};double rs=now();
            nmod_mpoly_struct *reference_out=flint_malloc((m-n+1)*sizeof(*reference_out));
            for(slong i=n-1;i<m;i++) {
                nmod_mpoly_struct sub[MAX_N];for(slong j=0;j<n-1;j++)sub[j]=fs[j];sub[n-1]=fs[i];
                nmod_mpoly_init(reference_out+i-n+1,qctx);
                construct_dixon(reference_out+i-n+1,sub,n,ctx,qctx,degree,0,&reference);
            }
            double repeated=now()-rs;
            printf("{\"method\":\"%s\",\"seconds\":%.9f,\"cpu_seconds\":%.9f,\"cache_seconds\":%.9f,\"cache_cpu_seconds\":%.9f,\"repeated_seconds\":%.9f,\"subsystems\":%ld,\"shared_term_products\":%llu,\"repeated_term_products\":%llu,\"audit\":%s}\n",
                method,elapsed,cpus,cache_seconds,cache_cpu,repeated,m-n+1,(unsigned long long)stats.products,(unsigned long long)reference.products,audit?"true":"false");
            for(slong i=0;i<m-n+1;i++) {
                if(audit)assert(nmod_mpoly_equal(reference_out+i,out+i,qctx));
                nmod_mpoly_clear(reference_out+i,qctx);nmod_mpoly_clear(out+i,qctx);
            }
            flint_free(reference_out);flint_free(out);
        } else {
            quotient_result r={0};nmod_poly_init(r.eliminant,prime);uint64_t cells=0,work=0;slong used_degree=degree;
            for(slong d=adaptive?2:degree;d<=degree;d++) {
                closure_t cl;closure_init(&cl,n,d,prime);nmod_mat_t M;
                original_rows(M,fs,m,ctx,&cl.basis,!oracle,d);
                if(oracle)closure_run(&cl,M,fs,m,ctx);else insert_rows(&cl,M);
                if(oracle && cl.solved==1) {
                    r.certified=1;r.dimension=1;memcpy(r.point,cl.point,n*sizeof(ulong));
                    nmod_poly_set_coeff_ui(r.eliminant,0,nmod_neg(cl.point[n-1],ctx->mod));nmod_poly_set_coeff_ui(r.eliminant,1,1);
                }else quotient_close(&r,&cl,fs,m,ctx);
                if(!r.certified) {
                    nmod_mat_clear(cl.processed);original_rows(cl.processed,fs,m,ctx,&cl.basis,1,d-1);
                    compact_rref(cl.processed,&cl);nmod_mat_t empty;matrix_init(empty,0,cl.basis.count,prime);
                    slong oldrank=cl.E->r;closure_run(&cl,empty,fs,m,ctx);nmod_mat_clear(empty);
                    if(cl.E->r>oldrank)quotient_close(&r,&cl,fs,m,ctx);
                }
                cells+=cl.rref_cells;work+=cl.matmul_work;used_degree=d;nmod_mat_clear(M);closure_clear(&cl);
                if(r.certified)break;
            }
            double elapsed=now()-start,cpus=(double)clock()/CLOCKS_PER_SEC-cpu;
            printf("{\"method\":\"%s\",\"seconds\":%.9f,\"cpu_seconds\":%.9f,\"certified\":%s,\"dimension\":%ld,\"degree\":%ld,\"rref_cells\":%llu,\"dense_matmul_multiply_adds\":%llu,\"dimensions\":[",
                method,elapsed,cpus,r.certified?"true":"false",r.dimension,used_degree,(unsigned long long)cells,(unsigned long long)work);
            for(slong j=0;j<r.ndims;j++)printf("%s%ld",j?",":"",r.dims[j]);
            printf("],\"recovered\":");
            if(r.certified && r.dimension==1){putchar('[');for(slong j=0;j<n;j++)printf("%s%lu",j?",":"",r.point[j]);putchar(']');}else printf("null");
            printf(",\"eliminant\":[");for(slong j=0;j<r.eliminant->length;j++)printf("%s%lu",j?",":"",r.eliminant->coeffs[j]);printf("]}\n");
            nmod_poly_clear(r.eliminant);
        }
        fflush(stdout);
    }
    for(slong i=0;i<m;i++)nmod_mpoly_clear(fs+i,ctx);
    flint_free(fs);
    for(slong i=0;i<n;i++)flint_free(names[i]);
    nmod_mpoly_ctx_clear(ctx);nmod_mpoly_ctx_clear(qctx);flint_cleanup_master();return 0;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Experimental shared construction and coefficient-panel Schur compression.
 * No full nmod_poly_mat candidate is needed by either direct variant. */
#define MQ_DIRECT_CORE_NO_MAIN
#include "mq_direct_core_test.c"

typedef struct {
    slong n, sigma;
    const slong *rd,*cd;
    size_t *offset;
    ulong *values;
    size_t slots;
} coefficient_source;

static ulong coefficient(const coefficient_source *s,slong r,slong c,slong d)
{
    slong bound=s->sigma-s->rd[r]-s->cd[c];
    return d>=0 && d<=bound ? s->values[s->offset[r*s->n+c]+d] : 0;
}
static int coefficient_panel(nmod_poly_mat_t out,panel_source *s,const slong *rows,
    slong nr,const slong *cols,slong nc)
{
    coefficient_source *p=s->data;
    for(slong i=0;i<nr;i++)for(slong j=0;j<nc;j++) {
        slong r=rows[i],c=cols[j],len=FLINT_MAX(0,p->sigma-p->rd[r]-p->cd[c]+1);
        nmod_poly_struct *a=nmod_poly_mat_entry(out,i,j);
        nmod_poly_fit_length(a,len);a->length=len;
        if(len)memcpy(a->coeffs,p->values+p->offset[r*p->n+c],len*sizeof(ulong));
        _nmod_poly_normalise(a);
    }
    return 1;
}
static int coefficient_source_init(coefficient_source *out,const panel_source *s,
                                    const fq_mq_compact *p)
{
    coefficient_source a={.n=s->size,.sigma=s->sigma,.rd=s->rd,.cd=s->cd};
    a.offset=flint_malloc((size_t)a.n*a.n*sizeof(size_t));
    for(slong r=0;r<a.n;r++)for(slong c=0;c<a.n;c++) {
        a.offset[r*a.n+c]=a.slots;
        a.slots+=FLINT_MAX(0,a.sigma-a.rd[r]-a.cd[c]+1);
    }
    a.values=flint_calloc(a.slots,sizeof(ulong));
    slong *rm=flint_malloc(p->nrows*sizeof(slong)),*cm=flint_malloc(p->ncols*sizeof(slong));
    for(int axis=0;axis<2;axis++) {
        slong count=axis?p->ncols:p->nrows,*map=axis?cm:rm;
        const slong *input=axis?p->cols:p->rows,*target=axis?s->cx:s->rx;
        for(slong i=0;i<count;i++) {
            map[i]=-1;
            for(slong j=0;j<s->size;j++)if(!memcmp(input+i*p->nvars,target+j*p->nvars,p->nvars*sizeof(slong))) {map[i]=j;break;}
            assert(map[i]>=0);
        }
    }
    for(slong i=0;i<p->nrows;i++)for(slong j=p->offset[i];j<p->offset[i+1];j++) {
        const fq_mq_compact_term *t=p->terms+j;
        slong r=rm[i],c=cm[t->column];
        assert(t->degree<=a.sigma-a.rd[r]-a.cd[c]);
        a.values[a.offset[r*a.n+c]+t->degree]=t->coefficient;
    }
    flint_free(rm);flint_free(cm);*out=a;return 1;
}

/* E_ij(t) X_j(t) is computed by multiplying each coefficient matrix of E_ij
 * by ALL coefficient blocks of X_j at once. The product is shifted into the
 * packed RHS. Only constant diagonal LU and scalar matrix kernels are used. */
static int coefficient_core(nmod_poly_mat_t core,ulong *factor,const coefficient_source *a,
    slong h,ulong prime,double *multiply,double *solve)
{
    slong n=a->n,e=n-h,groups=0,nr=0,nc=0;
    slong *rows=flint_malloc(e*sizeof(slong)),*cols=flint_malloc(e*sizeof(slong));
    mq_constant_block *g=flint_malloc(e*sizeof(*g));
    nmod_mat_struct *x=flint_calloc(e,sizeof(*x));
    nmod_t mod;nmod_init(&mod,prime);ulong scale=1;int ok=1;
    slong minr=WORD_MAX,minc=WORD_MAX;
    for(slong i=0;i<h;i++){minr=FLINT_MIN(minr,a->rd[i]);minc=FLINT_MIN(minc,a->cd[i]);}
    for(slong d=0;d<=a->sigma;d++) {
        slong begin=nr,cb=nc;
        for(slong i=h;i<n;i++)if(a->rd[i]==d)rows[nr++]=i;
        for(slong j=h;j<n;j++)if(a->cd[j]==a->sigma-d)cols[nc++]=j;
        if(nr-begin!=nc-cb){ok=0;break;}
        slong k=nr-begin;if(!k)continue;
        mq_constant_block *b=g+groups++;b->begin=begin;b->size=k;
        b->perm=flint_malloc(k*sizeof(slong));nmod_mat_init(b->lu,k,k,prime);
        for(slong i=0;i<k;i++)for(slong j=0;j<k;j++)
            nmod_mat_entry(b->lu,i,j)=coefficient(a,rows[begin+i],cols[begin+j],0);
        if(nmod_mat_lu(b->perm,b->lu,1)!=k){ok=0;break;}
        for(slong i=0;i<k;i++)scale=nmod_mul(scale,nmod_mat_entry(b->lu,i,i),mod);
        if(mq_permutation_odd(b->perm,k,0))scale=nmod_neg(scale,mod);
    }
    if(nr!=e || nc!=e)ok=0;
    slong length=a->sigma-minr-minc+1;
    nmod_mat_t reduced;nmod_mat_init(reduced,h,h*length,prime);
    if(ok) {
        if(mq_permutation_odd(rows,e,h)^mq_permutation_odd(cols,e,h))scale=nmod_neg(scale,mod);
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)for(slong d=0;d<length;d++)
            nmod_mat_entry(reduced,i,d*h+j)=coefficient(a,i,j,d);
        for(slong b=groups;b-- >0;) {
            slong start=g[b].begin,k=g[b].size,degree=a->rd[rows[start]];
            slong len=a->sigma-degree-minc+1;
            nmod_mat_t rhs,perm,lower;nmod_mat_init(rhs,k,h*len,prime);
            nmod_mat_init(perm,k,h*len,prime);nmod_mat_init(lower,k,h*len,prime);
            nmod_mat_init(x+b,k,h*len,prime);
            for(slong i=0;i<k;i++)for(slong j=0;j<h;j++)for(slong d=0;d<len;d++)
                nmod_mat_entry(rhs,i,d*h+j)=coefficient(a,rows[start+i],j,d);
            double time=get_wall_time();
            for(slong c=b+1;c<groups;c++) {
                slong cs=g[c].begin,ck=g[c].size;
                slong bound=a->sigma-degree-a->cd[cols[cs]];
                nmod_mat_t entry,view;nmod_mat_init(entry,k,ck,prime);
                for(slong d=0;d<=bound;d++) {
                    for(slong i=0;i<k;i++)for(slong j=0;j<ck;j++)
                        nmod_mat_entry(entry,i,j)=coefficient(a,rows[start+i],cols[cs+j],d);
                    nmod_mat_window_init(view,rhs,0,d*h,k,d*h+x[c].c);
                    nmod_mat_submul(view,view,entry,x+c);
                    nmod_mat_window_clear(view);
                }
                nmod_mat_clear(entry);
            }
            *multiply+=get_wall_time()-time;time=get_wall_time();
            for(slong i=0;i<k;i++)for(slong j=0;j<h*len;j++)
                nmod_mat_entry(perm,i,j)=nmod_mat_entry(rhs,g[b].perm[i],j);
            nmod_mat_solve_tril(lower,g[b].lu,perm,1);
            nmod_mat_solve_triu(x+b,g[b].lu,lower,0);
            *solve+=get_wall_time()-time;
            nmod_mat_clear(lower);nmod_mat_clear(perm);nmod_mat_clear(rhs);
            time=get_wall_time();
            slong ubound=a->sigma-minr-a->cd[cols[start]];
            nmod_mat_t u,view;nmod_mat_init(u,h,k,prime);
            for(slong d=0;d<=ubound;d++) {
                for(slong i=0;i<h;i++)for(slong j=0;j<k;j++)
                    nmod_mat_entry(u,i,j)=coefficient(a,i,cols[start+j],d);
                nmod_mat_window_init(view,reduced,0,d*h,h,d*h+x[b].c);
                nmod_mat_submul(view,view,u,x+b);
                nmod_mat_window_clear(view);
            }
            nmod_mat_clear(u);*multiply+=get_wall_time()-time;
        }
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)for(slong d=0;d<length;d++)
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(core,i,j),d,nmod_mat_entry(reduced,i,d*h+j));
        *factor=scale;
    }
    nmod_mat_clear(reduced);
    for(slong b=0;b<groups;b++){nmod_mat_clear(g[b].lu);flint_free(g[b].perm);if(x[b].r)nmod_mat_clear(x+b);}
    flint_free(x);flint_free(g);flint_free(rows);flint_free(cols);return ok;
}


#ifndef MQ_SHARED_CORE_NO_MAIN
int main(int argc,char **argv)
{
    slong n=argc>1?atol(argv[1]):7;int threads=argc>2?atoi(argv[2]):1;
    ulong prime=argc>3?strtoul(argv[3],NULL,10):536870923,seed=argc>4?strtoul(argv[4],NULL,10):12345;
    assert(n>=3 && n<=9);omp_set_num_threads(threads);flint_set_num_threads(threads);
    g_dixon_verbose_level=0;g_dixon_det_cache_limit=100000;
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,seed,941);
    fq_nmod_ctx_t fq;fq_nmod_ctx_init_ui(fq,prime,1,"a");
    unified_mpoly_struct *polys=random_mq(n-1,fq,rng),**matrix,**a;
    build_fq_cancellation_matrix(&matrix,polys,n-1,1);perform_fq_matrix_row_operations(&a,&matrix,n-1,1);
    panel_source source={.n=n,.matrix=a};slong h;make_profile(&source,&h);
    double t=get_wall_time();fq_mq_compact compact={0};
    assert(compute_fq_det_mq_compact(&compact,a,n,source.rx,source.size,source.cx,source.size));
    double generation=get_wall_time()-t;t=get_wall_time();
    coefficient_source cache;assert(coefficient_source_init(&cache,&source,&compact));
    double packing=get_wall_time()-t;slong terms=compact.nterms;fq_mq_compact_clear(&compact);
    source.data=&cache;source.read=coefficient_panel;
    nmod_poly_mat_t core,coeff,expected,B;
    nmod_poly_mat_init(core,h,h,prime);nmod_poly_mat_init(coeff,h,h,prime);nmod_poly_mat_init(expected,h,h,prime);
    ulong factor=17,cf=17,ef=17;
    t=get_wall_time();int ok=direct_core(core,&factor,source.size,prime,source.rd,source.cd,h,source.sigma,&source);
    double direct=get_wall_time()-t,multiply=0,solve=0;t=get_wall_time();
    int ck=coefficient_core(coeff,&cf,&cache,h,prime,&multiply,&solve);double coefficient_time=get_wall_time()-t;
    slong *all=flint_malloc(source.size*sizeof(slong));for(slong i=0;i<source.size;i++)all[i]=i;
    t=get_wall_time();nmod_poly_mat_init(B,source.size,source.size,prime);
    assert(coefficient_panel(B,&source,all,source.size,all,source.size));double matrix_time=get_wall_time()-t;
    t=get_wall_time();int ref=nmod_poly_mat_mq_schur(expected,&ef,B,source.rd,source.cd,h,source.sigma);
    double schur=get_wall_time()-t;
    assert(ok==ref && ck==ref);
    if(ok){assert(factor==ef && cf==ef);assert(nmod_poly_mat_equal(core,expected));assert(nmod_poly_mat_equal(coeff,expected));}
    printf("{\"n\":%ld,\"p\":%lu,\"seed\":%lu,\"threads\":%d,\"size\":%ld,\"h\":%ld,\"ok\":%d,\"generation\":%.6f,\"packing\":%.6f,\"shared_panels\":%.6f,\"coefficient_core\":%.6f,\"multiply\":%.6f,\"solve\":%.6f,\"matrix\":%.6f,\"schur\":%.6f,\"terms\":%ld,\"slots\":%zu,\"equal\":true}\n",
        n,prime,seed,threads,source.size,h,ok,generation,packing,direct,coefficient_time,multiply,solve,matrix_time,schur,terms,cache.slots);
    flint_free(all);flint_free(cache.values);flint_free(cache.offset);
    nmod_poly_mat_clear(B);nmod_poly_mat_clear(core);nmod_poly_mat_clear(coeff);nmod_poly_mat_clear(expected);
    flint_free(source.rx);flint_free(source.cx);flint_free(source.rd);flint_free(source.cd);
    clear_input(polys,matrix,a,n-1);fq_nmod_ctx_clear(fq);flint_rand_clear(rng);flint_cleanup_master();return 0;
}

#endif

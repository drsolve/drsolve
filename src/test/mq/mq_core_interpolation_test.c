/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Research benchmark: specialize t before minor recursion, compute numeric
 * Schur cores, and interpolate ONLY the core. Never used by the solver.
 * The cached-input mode isolates interpolation/linear-algebra overhead.
 * Usage: n threads prime mode seed; mode 0=cached, 1=specialize before DP,
 * 2=two-pass column generation, 3=cached plus finite inverse expansion. */
#define MQ_SHARED_CORE_NO_MAIN
#include "mq_shared_core_test.c"

typedef struct {
    slong n,h,e,ng;
    slong *rows,*cols;
    mq_constant_block *g;
} point_plan;

static int point_plan_init(point_plan *p,const coefficient_source *a,slong h,ulong prime)
{
    *p=(point_plan){.n=a->n,.h=h,.e=a->n-h};
    p->rows=flint_malloc(p->e*sizeof(slong));p->cols=flint_malloc(p->e*sizeof(slong));
    p->g=flint_malloc(p->e*sizeof(*p->g));slong nr=0,nc=0;
    for(slong d=0;d<=a->sigma;d++) {
        slong begin=nr,cb=nc;
        for(slong i=h;i<a->n;i++)if(a->rd[i]==d)p->rows[nr++]=i;
        for(slong j=h;j<a->n;j++)if(a->cd[j]==a->sigma-d)p->cols[nc++]=j;
        assert(nr-begin==nc-cb);slong k=nr-begin;if(!k)continue;
        mq_constant_block *b=p->g+p->ng++;b->begin=begin;b->size=k;
        b->perm=flint_malloc(k*sizeof(slong));nmod_mat_init(b->lu,k,k,prime);
        for(slong i=0;i<k;i++)for(slong j=0;j<k;j++)
            nmod_mat_entry(b->lu,i,j)=coefficient(a,p->rows[begin+i],p->cols[begin+j],0);
        if(nmod_mat_lu(b->perm,b->lu,1)!=k)return 0;
    }
    return nr==p->e && nc==p->e;
}
static void point_plan_clear(point_plan *p)
{
    for(slong b=0;b<p->ng;b++){nmod_mat_clear(p->g[b].lu);flint_free(p->g[b].perm);}
    flint_free(p->g);flint_free(p->rows);flint_free(p->cols);
}
static void point_core(nmod_mat_t core,const nmod_mat_t a,const point_plan *p)
{
    slong h=p->h;ulong prime=a->mod.n;
    nmod_mat_struct *x=flint_calloc(p->ng,sizeof(*x));
    for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)nmod_mat_entry(core,i,j)=nmod_mat_entry(a,i,j);
    for(slong b=p->ng;b-- >0;) {
        const mq_constant_block *g=p->g+b;slong k=g->size,start=g->begin;
        nmod_mat_t rhs,perm,lower;nmod_mat_init(rhs,k,h,prime);
        nmod_mat_init(perm,k,h,prime);nmod_mat_init(lower,k,h,prime);nmod_mat_init(x+b,k,h,prime);
        for(slong i=0;i<k;i++)for(slong j=0;j<h;j++)nmod_mat_entry(rhs,i,j)=nmod_mat_entry(a,p->rows[start+i],j);
        for(slong c=b+1;c<p->ng;c++) {
            slong ck=p->g[c].size,cs=p->g[c].begin;
            nmod_mat_t entry;nmod_mat_init(entry,k,ck,prime);
            for(slong i=0;i<k;i++)for(slong j=0;j<ck;j++)
                nmod_mat_entry(entry,i,j)=nmod_mat_entry(a,p->rows[start+i],p->cols[cs+j]);
            nmod_mat_submul(rhs,rhs,entry,x+c);nmod_mat_clear(entry);
        }
        for(slong i=0;i<k;i++)for(slong j=0;j<h;j++)nmod_mat_entry(perm,i,j)=nmod_mat_entry(rhs,g->perm[i],j);
        nmod_mat_solve_tril(lower,g->lu,perm,1);nmod_mat_solve_triu(x+b,g->lu,lower,0);
        nmod_mat_t u;nmod_mat_init(u,h,k,prime);
        for(slong i=0;i<h;i++)for(slong j=0;j<k;j++)nmod_mat_entry(u,i,j)=nmod_mat_entry(a,i,p->cols[start+j]);
        nmod_mat_submul(core,core,u,x+b);
        nmod_mat_clear(u);nmod_mat_clear(rhs);nmod_mat_clear(perm);nmod_mat_clear(lower);
    }
    for(slong b=0;b<p->ng;b++)nmod_mat_clear(x+b);flint_free(x);
}

/* Test the inverse-free finite expansion E=D+N. The strictly block upper
 * triangular D^-1 N is nilpotent, so the series terminates exactly. This
 * exposes candidate products for earlier fusion but repeats block products. */
static void neumann_core(nmod_mat_t core,const nmod_mat_t a,const point_plan *p)
{
    slong h=p->h;ulong prime=a->mod.n;
    nmod_mat_struct *z=flint_calloc(p->ng,sizeof(*z)),*next=flint_calloc(p->ng,sizeof(*next));
    for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)nmod_mat_entry(core,i,j)=nmod_mat_entry(a,i,j);
    for(slong b=0;b<p->ng;b++) {
        const mq_constant_block *g=p->g+b;slong k=g->size,start=g->begin;
        nmod_mat_t rhs,lower;nmod_mat_init(rhs,k,h,prime);nmod_mat_init(lower,k,h,prime);
        nmod_mat_init(z+b,k,h,prime);nmod_mat_init(next+b,k,h,prime);
        for(slong i=0;i<k;i++)for(slong j=0;j<h;j++)
            nmod_mat_entry(rhs,i,j)=nmod_mat_entry(a,p->rows[start+g->perm[i]],j);
        nmod_mat_solve_tril(lower,g->lu,rhs,1);nmod_mat_solve_triu(z+b,g->lu,lower,0);
        nmod_mat_clear(rhs);nmod_mat_clear(lower);
    }
    for(slong step=0;step<p->ng;step++) {
        for(slong b=0;b<p->ng;b++) {
            const mq_constant_block *g=p->g+b;slong k=g->size,start=g->begin;
            nmod_mat_t u;nmod_mat_init(u,h,k,prime);
            for(slong i=0;i<h;i++)for(slong j=0;j<k;j++)nmod_mat_entry(u,i,j)=nmod_mat_entry(a,i,p->cols[start+j]);
            nmod_mat_submul(core,core,u,z+b);nmod_mat_clear(u);
            nmod_mat_t rhs,perm,lower;nmod_mat_init(rhs,k,h,prime);nmod_mat_init(perm,k,h,prime);nmod_mat_init(lower,k,h,prime);
            for(slong c=b+1;c<p->ng;c++) {
                slong ck=p->g[c].size,cs=p->g[c].begin;
                nmod_mat_t entry;nmod_mat_init(entry,k,ck,prime);
                for(slong i=0;i<k;i++)for(slong j=0;j<ck;j++)
                    nmod_mat_entry(entry,i,j)=nmod_mat_entry(a,p->rows[start+i],p->cols[cs+j]);
                nmod_mat_submul(rhs,rhs,entry,z+c);nmod_mat_clear(entry);
            }
            for(slong i=0;i<k;i++)for(slong j=0;j<h;j++)nmod_mat_entry(perm,i,j)=nmod_mat_entry(rhs,g->perm[i],j);
            nmod_mat_solve_tril(lower,g->lu,perm,1);nmod_mat_solve_triu(next+b,g->lu,lower,0);
            nmod_mat_clear(rhs);nmod_mat_clear(perm);nmod_mat_clear(lower);
        }
        nmod_mat_struct *tmp=z;z=next;next=tmp;
    }
    for(slong b=0;b<p->ng;b++){assert(nmod_mat_is_zero(z+b));nmod_mat_clear(z+b);nmod_mat_clear(next+b);}
    flint_free(z);flint_free(next);
}
static unified_mpoly_struct **specialize(unified_mpoly_struct **a,slong n,ulong value)
{
    unified_mpoly_struct **b=flint_malloc(n*sizeof(*b));
    fq_nmod_t c;fq_nmod_init(c,a[0][0].ctx);slong zero=0;
    nmod_t mod;nmod_init(&mod,fq_nmod_ctx_prime(a[0][0].ctx));
    for(slong i=0;i<n;i++) {
        b[i]=flint_calloc(n,sizeof(**b));
        for(slong j=0;j<n;j++) {
            dr_mpoly_init(&b[i][j],a[i][j].nvars,1,a[i][j].ctx);
            for(slong k=0;k<dr_mpoly_length(&a[i][j]);k++) {
                DR_MPOLY_TERM(term,&a[i][j],k);
                ulong v=nmod_mul(nmod_poly_get_coeff_ui(term.coeff,0),nmod_pow_ui(value,term.par_exp[0],mod),mod);
                fq_nmod_set_ui(c,v,a[i][j].ctx);
                if(v)dr_mpoly_add_term_fast(&b[i][j],term.var_exp,&zero,c);
            }
        }
    }
    fq_nmod_clear(c,a[0][0].ctx);return b;
}
static void free_specialized(unified_mpoly_struct **b,slong n)
{
    for(slong i=0;i<n;i++){for(slong j=0;j<n;j++)dr_mpoly_clear(&b[i][j]);flint_free(b[i]);}flint_free(b);
}

/* Two-pass prerequisite for a frozen left projection [I,-U E^-1]: first
 * generate complementary columns, then the retained columns. This measures
 * the existing pruned recursion, NOT a fused/transposed recurrence. */
static void split_generation(const panel_source *s,const coefficient_source *cache,
                             slong h,slong full_terms,double full_time)
{
    double elapsed[2];slong terms[2];
    for(int pass=0;pass<2;pass++) {
        slong first=pass?0:h,count=pass?h:s->size-h;
        double start=get_wall_time();fq_mq_compact p={0};
        assert(compute_fq_det_mq_compact(&p,s->matrix,s->n,s->rx,s->size,
                                       s->cx+first*(s->n-1),count));
        elapsed[pass]=get_wall_time()-start;terms[pass]=p.nterms;
        slong *rm=flint_malloc(p.nrows*sizeof(slong)),*cm=flint_malloc(p.ncols*sizeof(slong));
        for(int axis=0;axis<2;axis++) {
            slong total=axis?p.ncols:p.nrows,*map=axis?cm:rm;
            const slong *input=axis?p.cols:p.rows,*target=axis?s->cx:s->rx;
            for(slong i=0;i<total;i++) {
                map[i]=-1;
                for(slong j=0;j<s->size;j++)if(!memcmp(input+i*p.nvars,target+j*p.nvars,p.nvars*sizeof(slong))){map[i]=j;break;}
                assert(map[i]>=0);
                if(axis)assert(map[i]>=first && map[i]<first+count);
            }
        }
        for(slong r=0;r<p.nrows;r++)for(slong k=p.offset[r];k<p.offset[r+1];k++) {
            const fq_mq_compact_term *term=p.terms+k;
            assert(term->coefficient==coefficient(cache,rm[r],cm[term->column],term->degree));
        }
        flint_free(rm);flint_free(cm);fq_mq_compact_clear(&p);
    }
    assert(terms[0]+terms[1]==full_terms);
    printf("{\"n\":%ld,\"mode\":\"split_columns\",\"full_generation_packing\":%.6f,\"complement_generation\":%.6f,\"retained_generation\":%.6f,\"complement_terms\":%ld,\"retained_terms\":%ld,\"equal\":true}\n",
           s->n,full_time,elapsed[0],elapsed[1],terms[0],terms[1]);
}

int main(int argc,char **argv)
{
    slong n=argc>1?atol(argv[1]):7;int threads=argc>2?atoi(argv[2]):1;
    ulong prime=argc>3?strtoul(argv[3],NULL,10):536870923;
    int fresh=argc>4?atoi(argv[4]):1;ulong seed=argc>5?strtoul(argv[5],NULL,10):12345;
    assert(n>=3 && n<=9);omp_set_num_threads(threads);flint_set_num_threads(threads);
    g_dixon_verbose_level=0;g_dixon_det_cache_limit=100000;
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,seed,941);
    fq_nmod_ctx_t fq;fq_nmod_ctx_init_ui(fq,prime,1,"a");
    unified_mpoly_struct *polys=random_mq(n-1,fq,rng),**matrix,**a;
    build_fq_cancellation_matrix(&matrix,polys,n-1,1);perform_fq_matrix_row_operations(&a,&matrix,n-1,1);
    panel_source source={.n=n,.matrix=a};slong h;make_profile(&source,&h);
    double t=get_wall_time();fq_mq_compact compact={0};
    assert(compute_fq_det_mq_compact(&compact,a,n,source.rx,source.size,source.cx,source.size));
    coefficient_source cache;assert(coefficient_source_init(&cache,&source,&compact));
    slong terms=compact.nterms;fq_mq_compact_clear(&compact);double generation=get_wall_time()-t;
    if(fresh==2) {
        split_generation(&source,&cache,h,terms,generation);
        flint_free(cache.values);flint_free(cache.offset);
        flint_free(source.rx);flint_free(source.cx);flint_free(source.rd);flint_free(source.cd);
        clear_input(polys,matrix,a,n-1);fq_nmod_ctx_clear(fq);flint_rand_clear(rng);flint_cleanup_master();return 0;
    }
    nmod_poly_mat_t expected,result;nmod_poly_mat_init(expected,h,h,prime);nmod_poly_mat_init(result,h,h,prime);
    ulong factor;double mul=0,solve=0;t=get_wall_time();
    int ok=coefficient_core(expected,&factor,&cache,h,prime,&mul,&solve);double reference=get_wall_time()-t;
    point_plan plan;t=get_wall_time();int pok=point_plan_init(&plan,&cache,h,prime);assert(ok==pok);
    double planning=get_wall_time()-t;
    if(!ok) {printf("{\"n\":%ld,\"prime\":%lu,\"singular\":true}\n",n,prime);goto cleanup;}
    slong minr=WORD_MAX,minc=WORD_MAX;for(slong i=0;i<h;i++){minr=FLINT_MIN(minr,cache.rd[i]);minc=FLINT_MIN(minc,cache.cd[i]);}
    slong points=cache.sigma-minr-minc+1;
    if((ulong)points>prime){printf("{\"n\":%ld,\"prime\":%lu,\"insufficient_points\":%ld}\n",n,prime,points);goto cleanup;}
    ulong *xs=flint_malloc(points*sizeof(ulong)),*ys=flint_malloc(h*h*points*sizeof(ulong));
    nmod_mat_t numeric,core;nmod_mat_init(numeric,source.size,source.size,prime);nmod_mat_init(core,h,h,prime);
    double eval=0,compress=0,interp=0,neumann=0;
    for(slong q=0;q<points;q++) {
        fprintf(stderr,"point %ld/%ld\n",q+1,points);fflush(stderr);
        xs[q]=q;t=get_wall_time();
        if(fresh==1) {
            unified_mpoly_struct **special=specialize(a,n,q);fq_mq_compact sample={0};
            assert(compute_fq_det_mq_compact(&sample,special,n,source.rx,source.size,source.cx,source.size));
            coefficient_source sc;assert(coefficient_source_init(&sc,&source,&sample));
            if(q==0) {
                /* Fresh evaluation needs no oracle diagonal: recover and
                 * factor it from the first sample, then reuse it at all t. */
                point_plan_clear(&plan);double pt=get_wall_time();
                assert(point_plan_init(&plan,&sc,h,prime));planning=get_wall_time()-pt;
                t+=planning;
            }
            for(slong i=0;i<source.size;i++)for(slong j=0;j<source.size;j++)
                nmod_mat_entry(numeric,i,j)=coefficient(&sc,i,j,0);
            flint_free(sc.values);flint_free(sc.offset);fq_mq_compact_clear(&sample);free_specialized(special,n);
        } else {
            for(slong i=0;i<source.size;i++)for(slong j=0;j<source.size;j++) {
                slong len=FLINT_MAX(0,cache.sigma-cache.rd[i]-cache.cd[j]+1);
                nmod_mat_entry(numeric,i,j)=_nmod_poly_evaluate_nmod(cache.values+cache.offset[i*source.size+j],len,q,numeric->mod);
            }
        }
        eval+=get_wall_time()-t;t=get_wall_time();point_core(core,numeric,&plan);compress+=get_wall_time()-t;
        if(fresh==3) {
            nmod_mat_t expanded;nmod_mat_init(expanded,h,h,prime);t=get_wall_time();
            neumann_core(expanded,numeric,&plan);neumann+=get_wall_time()-t;
            assert(nmod_mat_equal(expanded,core));nmod_mat_clear(expanded);
        }
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)ys[(i*h+j)*points+q]=nmod_mat_entry(core,i,j);
    }
    t=get_wall_time();size_t nonzero=0,slots=0,raw_core_nonzero=0,raw_core_slots=0;
    for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
        nmod_poly_struct *v=nmod_poly_mat_entry(result,i,j);
        nmod_poly_interpolate_nmod_vec(v,xs,ys+(i*h+j)*points,points);
        slong bound=FLINT_MAX(0,cache.sigma-cache.rd[i]-cache.cd[j]+1);slots+=bound;
        for(slong d=0;d<bound;d++){nonzero+=nmod_poly_get_coeff_ui(v,d)!=0;raw_core_nonzero+=coefficient(&cache,i,j,d)!=0;}
        raw_core_slots+=bound;
    }
    interp=get_wall_time()-t;assert(nmod_poly_mat_equal(expected,result));
    if(fresh==3)fprintf(stderr,"Neumann groups=%ld seconds=%.6f blocked_seconds=%.6f\n",plan.ng,neumann,compress);
    size_t left=0,right=0;
    for(slong i=0;i<source.size;i++)for(slong j=0;j<source.size;j++) {
        size_t len=FLINT_MAX(0,cache.sigma-cache.rd[i]-cache.cd[j]+1);
        if(j>=h)left+=len;if(i>=h)right+=len;
    }
    printf("{\"n\":%ld,\"prime\":%lu,\"threads\":%d,\"fresh\":%d,\"size\":%ld,\"h\":%ld,\"points\":%ld,\"generation\":%.6f,\"coefficient_core\":%.6f,\"planning\":%.6f,\"evaluation\":%.6f,\"point_core\":%.6f,\"interpolation\":%.6f,\"terms\":%ld,\"slots\":%zu,\"core_nonzero\":%zu,\"core_slots\":%zu,\"raw_core_nonzero\":%zu,\"raw_core_slots\":%zu,\"left_weight_input_slots\":%zu,\"right_weight_input_slots\":%zu,\"equal\":true}\n",n,prime,threads,fresh,source.size,h,points,generation,reference,planning,eval,compress,interp,terms,cache.slots,nonzero,slots,raw_core_nonzero,raw_core_slots,left,right);
    nmod_mat_clear(numeric);nmod_mat_clear(core);flint_free(xs);flint_free(ys);
cleanup:
    point_plan_clear(&plan);flint_free(cache.values);flint_free(cache.offset);
    nmod_poly_mat_clear(expected);nmod_poly_mat_clear(result);
    flint_free(source.rx);flint_free(source.cx);flint_free(source.rd);flint_free(source.cd);
    clear_input(polys,matrix,a,n-1);fq_nmod_ctx_clear(fq);flint_rand_clear(rng);flint_cleanup_master();return 0;
}

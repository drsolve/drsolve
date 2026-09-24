/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_pencil_det.h"
#include "mq_coefficient_filter.h"
#include "fq_mpoly_mat_det.h"
#include <flint/nmod_mat.h>
#include <math.h>

static double pencil_seconds(void)
{
    struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);
    return t.tv_sec+1e-9*t.tv_nsec;
}

static nmod_mpoly_struct *pencil_polys(slong count,const nmod_mpoly_ctx_t ctx)
{
    nmod_mpoly_struct *a=flint_malloc((size_t)count*sizeof(*a));
    for(slong i=0;i<count;i++)nmod_mpoly_init(a+i,ctx);
    return a;
}
static void pencil_clear(nmod_mpoly_struct *a,slong count,const nmod_mpoly_ctx_t ctx)
{
    for(slong i=0;i<count;i++)nmod_mpoly_clear(a+i,ctx);
    flint_free(a);
}
/* Split the quadratic border row once. Recurrence and coefficient buckets
 * remain parameter-free; append parameter exponents only at the end. */
static void pencil_split(nmod_mpoly_struct *out,const nmod_mpoly_t in,
                         const nmod_mpoly_ctx_t ctx)
{
    ulong exp[65];slong t=ctx->minfo->nvars-1;
    for(slong i=0;i<in->length;i++) {
        nmod_mpoly_get_term_exp_ui(exp,in,i,ctx);
        ulong degree=exp[t];exp[t]=0;
        nmod_mpoly_push_term_ui_ui(out+degree,in->coeffs[i],exp,ctx);
    }
    for(slong d=0;d<3;d++)nmod_mpoly_sort_terms(out+d,ctx);
}

static int pencil_compute(fq_mvpoly_t *result,fq_mvpoly_t **matrix,
                           slong n,mq_pencil_stats *stats,
                           const slong *rows,slong nr,const slong *cols,slong nc)
{
    mq_pencil_stats local={0};if(!stats)stats=&local;
    memset(stats,0,sizeof(*stats));stats->size=n;
    double start=pencil_seconds();
    stats->reason="requires prime-field, single-parameter MQ divided differences";
    if(!matrix || n<2 || n>32)return 0;
    const fq_nmod_ctx_struct *fq=matrix[0][0].ctx;
    if(fq_nmod_ctx_degree(fq)!=1)return 0;
    slong m=n-1,nv=2*m;
    ulong prime=fq_nmod_ctx_modulus(fq)->mod.n;
    if(prime<=(ulong)m) {
        stats->reason="characteristic must exceed n-1 for the degree recurrence";return 0;
    }
    for(slong i=0;i<n;i++)for(slong j=0;j<n;j++) {
        const fq_mvpoly_t *a=&matrix[i][j];
        if(a->nvars!=nv || a->npars!=1)return 0;
        for(slong k=0;k<a->nterms;k++) {
            slong degree=a->terms[k].par_exp[0];
            if(degree<0 || degree>(i?1:2))return 0;
            for(slong v=0;v<nv;v++) {
                slong e=a->terms[k].var_exp[v];
                if(e<0 || e>(i?1:2))return 0;
                degree+=e;
            }
            if(degree>(i?1:2))return 0;
        }
    }
    /* This first implementation stores two polynomial matrices. Bound dense
     * coefficient slots before allocation, rather than risking enormous runs.
     * n=8 and n=9 fit; the sparse minor backend handles larger workspaces. */
    double slots=2.0*m*m;
    for(slong k=1;k<=m-1;k++)slots*=((double)nv+k)/k;
    if(slots>33554432.0) {
        stats->reason="degree recurrence workspace exceeds 33554432 coefficient slots";return 0;
    }
    nmod_mat_t C,P,Pi,T;
    nmod_mat_init(C,m,n,prime);nmod_mat_init(P,m,m,prime);
    nmod_mat_init(Pi,m,m,prime);nmod_mat_init(T,n,n,prime);
    for(slong i=0;i<m;i++)for(slong j=0;j<n;j++) {
        const fq_mvpoly_t *a=&matrix[i+1][j];
        for(slong k=0;k<a->nterms;k++)if(a->terms[k].par_exp[0]==1)
            nmod_mat_entry(C,i,j)=nmod_add(nmod_mat_entry(C,i,j),
                nmod_poly_get_coeff_ui(a->terms[k].coeff,0),C->mod);
    }
    slong piv[32];int normalized=0;
    for(slong free_col=0;free_col<n && !normalized;free_col++) {
        slong pos=0;for(slong j=0;j<n;j++)if(j!=free_col)piv[pos++]=j;
        for(slong i=0;i<m;i++)for(slong j=0;j<m;j++)nmod_mat_entry(P,i,j)=nmod_mat_entry(C,i,piv[j]);
        if(!nmod_mat_inv(Pi,P))continue;
        nmod_mat_zero(T);nmod_mat_entry(T,free_col,0)=1;
        for(slong i=0;i<m;i++) {
            ulong c=0;
            for(slong j=0;j<m;j++) {
                ulong z=nmod_mat_entry(Pi,i,j);
                nmod_mat_entry(T,piv[i],j+1)=z;
                c=nmod_add(c,nmod_mul(z,nmod_mat_entry(C,j,free_col),C->mod),C->mod);
            }
            nmod_mat_entry(T,piv[i],0)=nmod_neg(c,C->mod);
        }
        normalized=1;
    }
    if(!normalized) {
        nmod_mat_clear(C);nmod_mat_clear(P);nmod_mat_clear(Pi);nmod_mat_clear(T);
        stats->reason="parameter coefficient matrix is rank deficient";return 0;
    }
    ulong scale=nmod_inv(nmod_mat_det(T),C->mod);
    nmod_mpoly_ctx_t ctx;nmod_mpoly_ctx_init(ctx,nv+1,ORD_LEX,prime);
    void *filter=NULL;
    if(rows) {
        filter=mq_coefficient_filter_create(m,rows,nr,cols,nc,ctx,n+1);
        if(!filter) {
            nmod_mat_clear(C);nmod_mat_clear(P);nmod_mat_clear(Pi);nmod_mat_clear(T);
            nmod_mpoly_ctx_clear(ctx);stats->reason="unsupported coefficient closure";return 0;
        }
    }
    nmod_mpoly_struct *input=pencil_polys(n*n,ctx),*a=pencil_polys(n*n,ctx);
    nmod_mpoly_t product,power,trace,coefficient,answer;
    nmod_mpoly_init(product,ctx);nmod_mpoly_init(power,ctx);
    nmod_mpoly_init(trace,ctx);
    nmod_mpoly_init(coefficient,ctx);nmod_mpoly_init(answer,ctx);
    for(slong i=0;i<n;i++)for(slong j=0;j<n;j++)fq_mvpoly_to_nmod_mpoly(input+i*n+j,&matrix[i][j],ctx);
    for(slong i=0;i<n;i++)for(slong j=0;j<n;j++)for(slong k=0;k<n;k++) {
        ulong c=nmod_mat_entry(T,k,j);if(!c)continue;
        nmod_mpoly_scalar_mul_ui(product,input+i*n+k,c,ctx);
        nmod_mpoly_add(a+i*n+j,a+i*n+j,product,ctx);
    }
    pencil_clear(input,n*n,ctx);
    nmod_mat_clear(C);nmod_mat_clear(P);nmod_mat_clear(Pi);nmod_mat_clear(T);
    /* Lower-right block tI+L -> K=-L. All lower entries are then t-free. */
    nmod_mpoly_gen(power,nv,ctx);
    for(slong i=0;i<m;i++)for(slong j=0;j<m;j++) {
        nmod_mpoly_struct *entry=a+(i+1)*n+j+1;
        if(i==j)nmod_mpoly_sub(entry,entry,power,ctx);
        nmod_mpoly_neg(entry,entry,ctx);
    }
    nmod_mpoly_struct *B=pencil_polys(m*m,ctx),*next=pencil_polys(m*m,ctx);
    nmod_mpoly_struct *parts=pencil_polys(3*m,ctx);
    nmod_mpoly_struct *q=pencil_polys(3*n,ctx),*buckets=pencil_polys(n+2,ctx);
    for(slong i=0;i<n;i++)pencil_split(q+3*i,a+i,ctx);
    for(slong i=0;i<m;i++)nmod_mpoly_one(B+i*m+i,ctx);
    for(slong i=n;i<n*n;i++)mq_coefficient_filter_repack(a+i,ctx,filter);
    for(slong i=0;i<m*m;i++)mq_coefficient_filter_repack(B+i,ctx,filter);
    nmod_mpoly_one(coefficient,ctx);
    for(slong d=0;d<3;d++)nmod_mpoly_set(buckets+m+d,q+d,ctx);
    stats->normalization=pencil_seconds()-start;
    slong threads=1;
#ifdef _OPENMP
    if(!omp_in_parallel())threads=FLINT_MIN(omp_get_max_threads(),m*m);
#endif
    stats->threads=threads;
    /* c_0=1, B_0=I. Before each update, append q B_k u to the border
     * polynomial. Only two matrix layers and the current c_k are retained. */
    for(slong k=0;k<m;k++) {
        double phase=pencil_seconds();
#ifdef _OPENMP
#pragma omp parallel num_threads(threads) if(threads>1)
#endif
        {
            nmod_mpoly_t p,v;nmod_mpoly_init(p,ctx);nmod_mpoly_init(v,ctx);
#ifdef _OPENMP
#pragma omp for schedule(dynamic,1)
#endif
            for(slong i=0;i<m;i++) {
                nmod_mpoly_zero(v,ctx);
                for(slong j=0;j<m;j++) {
                    nmod_mpoly_mul(p,a+(j+1)*n,B+i*m+j,ctx);
                    nmod_mpoly_add(v,v,p,ctx);
                }
                mq_coefficient_filter_apply(v,ctx,filter,0);
                for(slong d=0;d<3;d++) {
                    nmod_mpoly_mul(parts+3*i+d,q+3*(i+1)+d,v,ctx);
                    mq_coefficient_filter_apply(parts+3*i+d,ctx,filter,1);
                }
            }
            nmod_mpoly_clear(p,ctx);nmod_mpoly_clear(v,ctx);
        }
        for(slong d=0;d<3;d++)for(slong i=0;i<m;i++)
            nmod_mpoly_sub(buckets+m-1-k+d,buckets+m-1-k+d,parts+3*i+d,ctx);
        stats->assembly+=pencil_seconds()-phase;phase=pencil_seconds();
        /* The last update only needs the trace: B_m=0 is never materialized. */
        slong count=k==m-1?m:m*m;
#ifdef _OPENMP
#pragma omp parallel num_threads(threads) if(threads>1)
#endif
        {
            nmod_mpoly_t p;nmod_mpoly_init(p,ctx);
#ifdef _OPENMP
#pragma omp for schedule(dynamic,1)
#endif
            for(slong idx=0;idx<count;idx++) {
                slong i=k==m-1?idx:idx/m,j=k==m-1?idx:idx%m;
                nmod_mpoly_struct *entry=next+i*m+j;
                nmod_mpoly_zero(entry,ctx);
                for(slong l=0;l<m;l++) {
                    nmod_mpoly_mul(p,a+(i+1)*n+l+1,B+l*m+j,ctx);
                    nmod_mpoly_add(entry,entry,p,ctx);
                }
                mq_coefficient_filter_apply(entry,ctx,filter,0);
            }
            nmod_mpoly_clear(p,ctx);
        }
        nmod_mpoly_zero(trace,ctx);
        for(slong i=0;i<m;i++)nmod_mpoly_add(trace,trace,next+i*m+i,ctx);
        nmod_mpoly_scalar_mul_ui(coefficient,trace,nmod_neg(nmod_inv(k+1,ctx->mod),ctx->mod),ctx);
        stats->recurrence+=pencil_seconds()-phase;phase=pencil_seconds();
        for(slong d=0;d<3;d++) {
            nmod_mpoly_mul(product,q+d,coefficient,ctx);
            mq_coefficient_filter_apply(product,ctx,filter,1);
            nmod_mpoly_add(buckets+m-1-k+d,buckets+m-1-k+d,product,ctx);
        }
        stats->assembly+=pencil_seconds()-phase;phase=pencil_seconds();
        if(k<m-1)for(slong i=0;i<m;i++)nmod_mpoly_add(next+i*m+i,next+i*m+i,coefficient,ctx);
        for(slong i=0;i<m*m;i++)mq_coefficient_filter_repack(next+i,ctx,filter);
        slong terms=0;
        for(slong i=0;i<m*m;i++)terms+=B[i].length+next[i].length;
        stats->peak_terms=FLINT_MAX(stats->peak_terms,terms);
        nmod_mpoly_struct *swap=B;B=next;next=swap;
        stats->recurrence+=pencil_seconds()-phase;
    }
    double phase=pencil_seconds();
    for(slong d=0;d<n+2;d++)mq_coefficient_filter_apply(buckets+d,ctx,filter,1);
    ulong exp[65];
    for(slong d=0;d<n+2;d++)for(slong i=0;i<buckets[d].length;i++) {
        nmod_mpoly_get_term_exp_ui(exp,buckets+d,i,ctx);exp[nv]=d;
        nmod_mpoly_push_term_ui_ui(answer,nmod_mul(buckets[d].coeffs[i],scale,ctx->mod),exp,ctx);
    }
    nmod_mpoly_sort_terms(answer,ctx);
    nmod_mpoly_to_fq_mvpoly(result,answer,nv,1,ctx,fq);
    stats->assembly+=pencil_seconds()-phase;
    pencil_clear(B,m*m,ctx);pencil_clear(next,m*m,ctx);pencil_clear(parts,3*m,ctx);pencil_clear(q,3*n,ctx);pencil_clear(buckets,n+2,ctx);pencil_clear(a,n*n,ctx);
    nmod_mpoly_clear(product,ctx);nmod_mpoly_clear(power,ctx);
    nmod_mpoly_clear(trace,ctx);
    nmod_mpoly_clear(coefficient,ctx);nmod_mpoly_clear(answer,ctx);nmod_mpoly_ctx_clear(ctx);
    mq_coefficient_filter_destroy(filter);
    stats->total=pencil_seconds()-start;stats->reason=NULL;return 1;
}

int compute_fq_det_mq_pencil(fq_mvpoly_t *result,fq_mvpoly_t **matrix,
                           slong n,mq_pencil_stats *stats)
{
    return pencil_compute(result,matrix,n,stats,NULL,0,NULL,0);
}
int compute_fq_det_mq_pencil_projected(fq_mvpoly_t *result,fq_mvpoly_t **matrix,
    slong n,const slong *rows,slong nr,const slong *cols,slong nc,mq_pencil_stats *stats)
{
    if(!rows || !cols || nr<=0 || nc<=0) {
        if(stats){memset(stats,0,sizeof(*stats));stats->reason="empty projection targets";}
        return 0;
    }
    return pencil_compute(result,matrix,n,stats,rows,nr,cols,nc);
}

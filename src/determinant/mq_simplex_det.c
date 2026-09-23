/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "mq_simplex_det.h"
#include <flint/nmod_mat.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Bounded dense workspace; n=8 has 1,307,504 points, n=9 has 8,436,285.
 * Larger jobs can still use the existing sparse minor-DP path. */
#define MQ_SIMPLEX_MAX_N 16
#define MQ_SIMPLEX_MAX_POINTS (WORD(1) << 24)
#define MQ_SIMPLEX_BINOM (3*MQ_SIMPLEX_MAX_N+1)
#define MQ_SIMPLEX_VARS (2*MQ_SIMPLEX_MAX_N-1)
#define MQ_SIMPLEX_DEG (MQ_SIMPLEX_MAX_N+1)
#define MQ_SIMPLEX_BATCH 128

typedef struct { ulong coeff; int x, y; } mq_eval_term;

static double simplex_seconds(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9*t.tv_nsec;
}

static slong simplex_rank(const uint8_t *e, slong v, slong d,
                          const slong choose[MQ_SIMPLEX_BINOM][MQ_SIMPLEX_BINOM])
{
    slong rank=0;
    for (slong j=0;j<v;j++) {
        slong left=v-j;
        rank+=choose[d+left][left]-choose[d-e[j]+left][left];
        d-=e[j];
    }
    return rank;
}

static void simplex_enumerate(uint8_t *exps, uint8_t *e, slong *count,
                              slong v, slong j, slong d)
{
    if (j==v) { memcpy(exps+(*count)++*v,e,(size_t)v); return; }
    for (slong k=0;k<=d;k++) {
        e[j]=(uint8_t)k;
        simplex_enumerate(exps,e,count,v,j+1,d-k);
    }
}

int compute_fq_det_mq_simplex(fq_mvpoly_t *result, fq_mvpoly_t **matrix,
                            slong n, mq_simplex_stats *stats)
{
    mq_simplex_stats local={0};
    if (!stats) stats=&local;
    memset(stats,0,sizeof(*stats));
    double start=simplex_seconds();
    stats->reason="requires prime-field, single-parameter MQ divided differences";
    if (!matrix || n<2 || n>MQ_SIMPLEX_MAX_N) return 0;
    const fq_nmod_ctx_struct *fq=matrix[0][0].ctx;
    if (fq_nmod_ctx_degree(fq)!=1) return 0;
    ulong prime=fq_nmod_ctx_modulus(fq)->mod.n;
    slong v=2*n-1,d=n+1;
    if (prime <= (ulong)d) {
        stats->reason="base field needs at least n+2 distinct nodes";
        return 0;
    }
    slong total_terms=0;
    for (slong i=0;i<n;i++) for (slong j=0;j<n;j++) {
        const fq_mvpoly_t *p=&matrix[i][j];
        if (p->nvars!=v-1 || p->npars!=1) return 0;
        for (slong k=0;k<p->nterms;k++) {
            slong degree=p->terms[k].par_exp[0];
            if (degree<0 || degree>(i?1:2)) return 0;
            for (slong l=0;l<v-1;l++) {
                slong a=p->terms[k].var_exp[l];
                if (a<0 || a>(i?1:2)) return 0;
                degree+=a;
            }
            if (degree>(i?1:2)) return 0;
        }
        if (p->nterms>WORD_MAX-total_terms) return 0;
        total_terms+=p->nterms;
    }
    slong choose[MQ_SIMPLEX_BINOM][MQ_SIMPLEX_BINOM]={{0}};
    for (slong i=0;i<=v+d;i++) {
        choose[i][0]=1;
        for (slong j=1;j<=i;j++)
            choose[i][j]=FLINT_MIN(MQ_SIMPLEX_MAX_POINTS+1,choose[i-1][j-1]+choose[i-1][j]);
    }
    slong N=choose[v+d][v];
    stats->points=N;
    if (N>MQ_SIMPLEX_MAX_POINTS) {
        stats->reason="simplex dense workspace exceeds 16777216 points";
        return 0;
    }
    slong threads=1;
#ifdef _OPENMP
    if (!omp_in_parallel()) threads=FLINT_MIN(omp_get_max_threads(),FLINT_MAX(1,N/MQ_SIMPLEX_BATCH));
#endif
    stats->threads=threads;
    slong batch=MQ_SIMPLEX_BATCH*threads;
    uint8_t *exps=malloc((size_t)N*v);
    ulong *values=malloc((size_t)N*sizeof(ulong));
    mq_eval_term *terms=malloc((size_t)FLINT_MAX(1,total_terms)*sizeof(*terms));
    nmod_mat_struct *mats=malloc((size_t)batch*sizeof(*mats));
    if (!exps || !values || !terms || !mats) {
        free(exps);free(values);free(terms);free(mats);
        stats->reason="simplex workspace allocation failed";return 0;
    }
    uint8_t e[MQ_SIMPLEX_VARS]; slong count=0;
    simplex_enumerate(exps,e,&count,v,0,d);
    slong offsets[MQ_SIMPLEX_MAX_N*MQ_SIMPLEX_MAX_N+1];
    slong pos=0;
    for (slong i=0;i<n;i++) for (slong j=0;j<n;j++) {
        offsets[i*n+j]=pos;
        const fq_mvpoly_t *p=&matrix[i][j];
        for (slong k=0;k<p->nterms;k++) {
            mq_eval_term *t=terms+pos++;
            t->coeff=nmod_poly_get_coeff_ui(p->terms[k].coeff,0);
            t->x=t->y=-1;
            for (slong l=0;l<v;l++) {
                slong power=l==v-1?p->terms[k].par_exp[0]:p->terms[k].var_exp[l];
                for (slong z=0;z<power;z++) {
                    if (t->x<0)t->x=(int)l;else t->y=(int)l;
                }
            }
        }
    }
    offsets[n*n]=pos;
    for (slong b=0;b<batch;b++)nmod_mat_init(mats+b,n,n,prime);
    nmod_t mod; nmod_init(&mod,prime);
    stats->setup=simplex_seconds()-start;
    double phase=0;
    /* One team, bounded batches of numerical matrices. Separate timed stages
     * include their worksharing barrier; no sum of thread times is reported. */
#ifdef _OPENMP
#pragma omp parallel num_threads(threads) if(threads>1)
#endif
    {
        for (slong base=0;base<N;base+=batch) {
            slong len=FLINT_MIN(batch,N-base);
#ifdef _OPENMP
#pragma omp single
#endif
            phase=simplex_seconds();
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (slong b=0;b<len;b++) {
                const uint8_t *point=exps+(base+b)*v;
                for (slong ij=0;ij<n*n;ij++) {
                    ulong sum=0;
                    for (slong k=offsets[ij];k<offsets[ij+1];k++) {
                        mq_eval_term term=terms[k]; ulong z=term.coeff;
                        if (term.x>=0)z=nmod_mul(z,point[term.x],mod);
                        if (term.y>=0)z=nmod_mul(z,point[term.y],mod);
                        sum=nmod_add(sum,z,mod);
                    }
                    nmod_mat_entry(mats+b,ij/n,ij%n)=sum;
                }
            }
#ifdef _OPENMP
#pragma omp single
#endif
            {stats->entry_eval+=simplex_seconds()-phase;phase=simplex_seconds();}
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (slong b=0;b<len;b++)values[base+b]=nmod_mat_det(mats+b);
#ifdef _OPENMP
#pragma omp single
#endif
            stats->determinants+=simplex_seconds()-phase;
        }
    }
    for (slong b=0;b<batch;b++)nmod_mat_clear(mats+b);
    free(mats);free(terms);
    phase=simplex_seconds();
    ulong inv[MQ_SIMPLEX_DEG+1],basis[MQ_SIMPLEX_DEG+1][MQ_SIMPLEX_DEG+1]={{0}};
    for (slong k=1;k<=d;k++)inv[k]=nmod_inv(k,mod);
    basis[0][0]=1;
    for (slong k=1;k<=d;k++)for(slong j=0;j<=k;j++)
        basis[k][j]=nmod_sub(j?basis[k-1][j-1]:0,nmod_mul(k-1,basis[k-1][j],mod),mod);
    /* Fibers on each axis are disjoint. Barrier before changing axes or
     * converting the Newton basis prevents cross-axis read/write races. */
#ifdef _OPENMP
#pragma omp parallel num_threads(threads) if(threads>1)
#endif
    {
        uint8_t point[MQ_SIMPLEX_VARS];
        slong index[MQ_SIMPLEX_DEG+1];
        ulong line[MQ_SIMPLEX_DEG+1],converted[MQ_SIMPLEX_DEG+1];
        for (int pass=0;pass<2;pass++)for(slong axis=0;axis<v;axis++) {
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
            for (slong i=0;i<N;i++) {
                if (exps[i*v+axis])continue;
                memcpy(point,exps+i*v,(size_t)v);slong degree=0;
                for(slong j=0;j<v;j++)degree+=point[j];
                slong len=d-degree;
                for(slong k=0;k<=len;k++) {
                    point[axis]=(uint8_t)k;index[k]=simplex_rank(point,v,d,choose);
                    line[k]=values[index[k]];converted[k]=0;
                }
                if (!pass) {
                    for(slong k=1;k<=len;k++)for(slong j=len;j>=k;j--)
                        line[j]=nmod_mul(nmod_sub(line[j],line[j-1],mod),inv[k],mod);
                } else {
                    for(slong k=0;k<=len;k++)for(slong j=0;j<=k;j++)
                        converted[j]=nmod_add(converted[j],nmod_mul(line[k],basis[k][j],mod),mod);
                    memcpy(line,converted,(size_t)(len+1)*sizeof(ulong));
                }
                for(slong k=0;k<=len;k++)values[index[k]]=line[k];
            }
        }
    }
    stats->interpolation=simplex_seconds()-phase;
    phase=simplex_seconds();
    fq_mvpoly_init(result,v-1,1,fq);
    fq_nmod_t coefficient;fq_nmod_init(coefficient,fq);
    slong powers[MQ_SIMPLEX_VARS];
    for (slong i=N;i-- >0;)if(values[i]) {
        for(slong j=0;j<v;j++)powers[j]=exps[i*v+j];
        fq_nmod_set_ui(coefficient,values[i],fq);
        fq_mvpoly_add_term_fast(result,powers,powers+v-1,coefficient);
    }
    fq_nmod_clear(coefficient,fq);
    free(values);free(exps);
    stats->packing=simplex_seconds()-phase;
    stats->total=simplex_seconds()-start;stats->reason=NULL;
    return 1;
}

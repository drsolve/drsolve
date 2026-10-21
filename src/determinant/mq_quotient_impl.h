/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Research prototype: graded Macaulay normal forms, followed by Bezoutian
 * minor recursion in (S/I)_a tensor (S/I)_b. Shared by the opt-in solver
 * backend and its independent validation executable.
 * Input: n prime; n triangular quadratic coefficient arrays in variables
 * x0,...,x(n-1),1; optional reference degree and ascending coefficients.
 * Usage: input threads mode [kernel] [shared_threads]; mode 0=one point, 1=point-internal parallelism, 2=parallel samples, 3=shared leading frames plus parallel samples.
 * Kernels: -1=automatic (default), 0=FLINT updates, 1=fused accumulation,
 * 2=fused accumulation with zero-prefix/suffix skipping,
 * 3=packed 32-bit matrices, 4=packed matrices with layer-local maps,
 * 5=register-blocked packed updates (automatic),
 * 6=experimental batches of four parent minors.
 * Modes 1, 2 and 3 interpolate and verify; validation time is separate.
 * Mode 3 shares pure-x leading Macaulay pivot inverses across parameters.
 * Reduction modulo I in each tensor factor commutes with multiplication,
 * so the subset recurrence computes the reduced Bezoutian directly.
 * This is an experiment for dense generic quadrics: rank/basis/degree
 * assertions check internal invariants. The public adapter checks admission
 * and returns failure for the caller to use its existing resultant path.
 * All setup is included in compute_seconds. Extra evaluation points and
 * comparison against full Macaulay RREF are in validation_seconds.
 * The hidden variable is x(n-1), the last homogeneous variable is z. */
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly_mat.h>
#include "pml_det.h"

#ifndef MQ_QUOTIENT_EXTENSION
#include "degree_quotient_prime.h"
#endif

typedef struct {
    slong s,h,*free;
    uint64_t *keys;
    slong *hash,buckets;
    dq_mat_t reduction;
    dq_mat_t variables;
} grade;
typedef struct {slong r,c;uint32_t *data;} packed_matrix;
typedef struct {dq_mat_t x,y;slong *ystart,*yend;packed_matrix px,py;} linear_map;
static slong N,M,K,L,threads,key_bits;
/* Nested sample teams must not spawn another team. OpenMP is optional. */
static slong quotient_recurrence_workers(void)
{
#ifdef _OPENMP
    return omp_in_parallel() ? 1 : threads;
#else
    return 1;
#endif
}

/* Packed exponent keys and subset masks must fit their machine words.
 * There is no table or per-equation-count implementation. */
static int quotient_dimensions(slong n)
{
    if(n<2 || n>=FLINT_BITS-1)return 0;
    N=n;M=N/2;K=N-M;L=K;key_bits=1;
    while((UINT64_C(1)<<key_bits)<=(uint64_t)L)key_bits++;
    return key_bits*N<=64;
}
static int recurrence_kernel,verify_recurrence;
static ulong prime;
static nmod_t mod;
static double normal_time,map_time,dp_time;
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static slong choose(slong n,slong k){if(k<0||k>n)return 0;slong z=1;for(slong i=1;i<=k;i++)z=z*(n-i+1)/i;return z;}
static size_t hash_key(uint64_t x){x^=x>>30;x*=UINT64_C(0xbf58476d1ce4e5b9);x^=x>>27;x*=UINT64_C(0x94d049bb133111eb);return x^(x>>31);}
static void enum_x(grade *g,slong v,slong left,uint64_t key,slong *count)
{
    if(v==N-2){g->keys[(*count)++]=key|((uint64_t)left<<(key_bits*v));return;}
    for(slong a=left;a>=0;a--)enum_x(g,v+1,left-a,key|((uint64_t)a<<(key_bits*v)),count);
}
static slong index_key(const grade *g,uint64_t key)
{
    size_t pos=hash_key(key)&(g->buckets-1);
    while(g->hash[pos]){slong j=g->hash[pos]-1;if(g->keys[j]==key)return j;pos=(pos+1)&(g->buckets-1);}
    assert(0);return -1;
}
static void grade_support(grade *g,slong d)
{
    memset(g,0,sizeof(*g));g->s=choose(N+d-1,d);g->h=choose(N,d);
    g->keys=flint_malloc(g->s*sizeof(uint64_t));slong count=0;
    for(slong z=0;z<=d;z++)enum_x(g,0,d-z,(uint64_t)z<<(key_bits*(N-1)),&count);
    assert(count==g->s);g->buckets=16;while(g->buckets<2*g->s)g->buckets*=2;
    g->hash=flint_calloc(g->buckets,sizeof(slong));
    for(slong i=0;i<g->s;i++){size_t p=hash_key(g->keys[i])&(g->buckets-1);while(g->hash[p])p=(p+1)&(g->buckets-1);g->hash[p]=i+1;}
    g->free=flint_malloc(g->h*sizeof(slong));dq_mat_init(g->reduction,g->s,g->h,prime);
}
#include "mq_quotient_normal_forms.h"
static int normal_forms(grade *g,const dq_coeff *q)
{
    double t=now();
    if(frames_enabled)cached_normal_forms(g,q);
    else for(slong d=0;d<=L;d++) {
        grade_support(g+d,d);
        if(d<2) {
            for(slong i=0;i<g[d].s;i++){g[d].free[i]=i;dq_one(dq_mat_entry(g[d].reduction,i,i));}
            continue;
        }
        dq_mat_t relations;dq_mat_init(relations,N*g[d-2].s,g[d].s,prime);
        for(slong f=0;f<N;f++)for(slong k=0;k<g[d-2].s;k++)
            for(slong i=0;i<N;i++)for(slong j=i;j<N;j++) {
                uint64_t key=g[d-2].keys[k]+(UINT64_C(1)<<(key_bits*i))+(UINT64_C(1)<<(key_bits*j));
                dq_set(dq_mat_entry(relations,f*g[d-2].s+k,index_key(g+d,key)),q[(f*N+i)*N+j]);
            }
        slong rank=dq_mat_rref(relations);
        if(rank!=g[d].s-g[d].h){fprintf(stderr,"non-generic rank at degree %ld: %ld\n",d,rank);dq_mat_clear(relations);return 0;}
        slong *pivot=flint_malloc(rank*sizeof(slong)),r=0,nfree=0;
        for(slong c=0;c<g[d].s;c++) {
            if(r<rank && dq_nonzero(dq_mat_entry(relations,r,c))){pivot[r++]=c;}
            else {assert(nfree<g[d].h);g[d].free[nfree++]=c;}
        }
        assert(r==rank && nfree==g[d].h);
        for(slong j=0;j<g[d].h;j++) {
            dq_one(dq_mat_entry(g[d].reduction,g[d].free[j],j));
            for(slong i=0;i<rank;i++)dq_neg(dq_mat_entry(g[d].reduction,pivot[i],j),dq_mat_entry(relations,i,g[d].free[j]));
        }
        flint_free(pivot);dq_mat_clear(relations);
    }
    double elapsed=now()-t;
    #pragma omp atomic update
    normal_time+=elapsed;
    t=now();
    for(slong d=0;d<L;d++) {
        dq_mat_init(g[d].variables,N,g[d].h*g[d+1].h,prime);
        for(slong v=0;v<N;v++)for(slong i=0;i<g[d].h;i++) {
            slong k=index_key(g+d+1,g[d].keys[g[d].free[i]]+(UINT64_C(1)<<(key_bits*v)));
            dq_vec_set(dq_mat_entry_ptr(g[d].variables,v,i*g[d+1].h),dq_mat_entry_ptr(g[d+1].reduction,k,0),g[d+1].h);
        }
    }
    elapsed=now()-t;
    #pragma omp atomic update
    map_time+=elapsed;
    return 1;
}
static void clear_grades(grade *g)
{
    for(slong d=0;d<=L;d++) {
        if(!g[d].keys)continue;
        flint_free(g[d].keys);flint_free(g[d].hash);flint_free(g[d].free);dq_mat_clear(g[d].reduction);
        if(g[d].variables->r)dq_mat_clear(g[d].variables);
    }
}
#ifndef MQ_QUOTIENT_EXTENSION
static void pack_matrix(packed_matrix *out,dq_mat_t in)
{
    out->r=in->r;out->c=in->c;out->data=flint_malloc(out->r*out->c*sizeof(uint32_t));
    for(slong r=0;r<out->r;r++)for(slong c=0;c<out->c;c++)out->data[r*out->c+c]=(uint32_t)dq_mat_entry(in,r,c);
    dq_mat_clear(in);memset(in,0,sizeof(*in));
}
#endif
static void clear_linear_maps(linear_map *op)
{
    if(!op)return;
    for(slong i=0;i<N*N*L;i++) {
        if(op[i].x->r)dq_mat_clear(op[i].x);if(op[i].y->r)dq_mat_clear(op[i].y);
        flint_free(op[i].ystart);flint_free(op[i].yend);flint_free(op[i].px.data);flint_free(op[i].py.data);
    }
    flint_free(op);
}
static linear_map *make_linear_maps_for_row(const grade *g,const dq_coeff *q,slong only_row)
{
    double t=now();linear_map *op=flint_calloc(N*N*L,sizeof(*op));
    for(slong d=0;d<L;d++) {
        slong *selected=flint_malloc(2*N*N*sizeof(slong)),count=0;
        slong *ix=flint_malloc(N*N*sizeof(slong)),*iy=flint_malloc(N*N*sizeof(slong));
        for(slong f=0;f<N;f++)for(slong j=0;j<N;j++) {
            slong id=f*N+j,other=f-d;ix[id]=iy[id]=-1;
            if(only_row>=0 && f!=only_row)continue;
            if(d<M && other>=0 && other<=K){ix[id]=count;selected[count++]=2*id;}
            if(d<K && other>=0 && other<=M){iy[id]=count;selected[count++]=2*id+1;}
        }
        if(!count){flint_free(selected);flint_free(ix);flint_free(iy);continue;}
        slong width=g[d].h*g[d+1].h;
        dq_mat_t coefficients,combined;
        dq_mat_init(coefficients,count,N,prime);
        for(slong r=0;r<count;r++) {
            slong id=selected[r]/2,f=id/N,j=id%N;int side=selected[r]&1;
            for(slong v=0;v<N;v++)if(side?v<=j:v>=j)
                dq_set(dq_mat_entry(coefficients,r,v),q[(f*N+FLINT_MIN(j,v))*N+FLINT_MAX(j,v)]);
        }
        dq_mat_init(combined,count,width,prime);
        dq_mat_mul(combined,coefficients,g[d].variables);
        for(slong id=0;id<N*N;id++) {
            linear_map *a=op+id*L+d;
            if(ix[id]>=0) {
                dq_mat_init(a->x,g[d+1].h,g[d].h,prime);
                const dq_coeff *x=dq_mat_entry_ptr(combined,ix[id],0);
                for(slong r=0;r<g[d].h;r++)for(slong c=0;c<g[d+1].h;c++)
                    dq_set(dq_mat_entry(a->x,c,r),x[r*g[d+1].h+c]);
            }
            if(iy[id]>=0) {
                dq_mat_init(a->y,g[d].h,g[d+1].h,prime);
                const dq_coeff *y=dq_mat_entry_ptr(combined,iy[id],0);
                for(slong r=0;r<g[d].h;r++)
                    dq_vec_set(dq_mat_entry_ptr(a->y,r,0),y+r*g[d+1].h,g[d+1].h);
                if(recurrence_kernel==2) {
                    a->ystart=flint_malloc(a->y->r*sizeof(slong));a->yend=flint_malloc(a->y->r*sizeof(slong));
                    for(slong r=0;r<a->y->r;r++) {
                        slong lo=0,hi=a->y->c;
                        while(lo<hi && !dq_nonzero(dq_mat_entry(a->y,r,lo)))lo++;
                        while(hi>lo && !dq_nonzero(dq_mat_entry(a->y,r,hi-1)))hi--;
                        a->ystart[r]=lo;a->yend[r]=hi;
                    }
                }
            }
        }
#ifndef MQ_QUOTIENT_EXTENSION
        if(recurrence_kernel>=3)for(slong id=0;id<N*N;id++) {
            linear_map *a=op+id*L+d;
            if(a->x->r)pack_matrix(&a->px,a->x);
            if(a->y->r)pack_matrix(&a->py,a->y);
        }
#endif
        dq_mat_clear(coefficients);dq_mat_clear(combined);
        flint_free(selected);flint_free(ix);flint_free(iy);
    }
    double elapsed=now()-t;
    #pragma omp atomic update
    map_time+=elapsed;
    return op;
}
static linear_map *make_linear_maps(const grade *g,const dq_coeff *q)
{return make_linear_maps_for_row(g,q,-1);}
#ifndef MQ_QUOTIENT_EXTENSION
#include "mq_quotient_accumulate.h"
#include "mq_quotient_packed.h"
#else
#include "mq_quotient_byte.h"
#endif
static void quotient_det(dq_mat_t out,const grade *g,const linear_map *op,const dq_coeff *q)
{
    double t=now();
#ifndef MQ_QUOTIENT_EXTENSION
    if(recurrence_kernel>=3){quotient_det_packed(out,g,op,q);return;}
#endif
#ifdef MQ_QUOTIENT_EXTENSION
    if(mq_byte_enabled()){mq_byte_det(out,g,q);return;}
#endif
    slong masks=1L<<N,workers=quotient_recurrence_workers();
    dq_mat_struct *state=flint_calloc(masks*(M+1),sizeof(*state));
    dq_mat_init(state,1,1,prime);dq_one(dq_mat_entry(state,0,0));
    for(slong k=1;k<=N;k++) {
        linear_map *layer_maps=op?NULL:make_linear_maps_for_row(g,q,k-1);
        const linear_map *active_op=op?op:layer_maps;
        #pragma omp parallel for schedule(dynamic) num_threads(workers) if(workers>1)
        for(slong mask=1;mask<masks;mask++)if(__builtin_popcountl(mask)==k) {
            for(slong a=FLINT_MAX(0,k-K);a<=FLINT_MIN(M,k);a++) {
                slong b=k-a;dq_mat_struct *dest=state+mask*(M+1)+a;
                dq_mat_init(dest,g[a].h,g[b].h,prime);
#ifndef MQ_QUOTIENT_EXTENSION
                if(recurrence_kernel && fused_minor(dest,state,op,mask,k,a,b))continue;
#endif
                slong pos=0;
                for(slong j=0;j<N;j++)if(mask&(1L<<j)) {
                    slong child=mask^(1L<<j);int neg=(k-1+pos++)&1;
                    if(a) {
                        const dq_mat_struct *src=state+child*(M+1)+a-1;
                        const linear_map *map=active_op+((k-1)*N+j)*L+a-1;
                        if(neg)dq_mat_submul(dest,dest,map->x,src);else dq_mat_addmul(dest,dest,map->x,src);
                    }
                    if(b) {
                        const dq_mat_struct *src=state+child*(M+1)+a;
                        const linear_map *map=active_op+((k-1)*N+j)*L+b-1;
                        if(neg)dq_mat_submul(dest,dest,src,map->y);else dq_mat_addmul(dest,dest,src,map->y);
                    }
                }
            }
        }
        clear_linear_maps(layer_maps);
        for(slong mask=0;mask<masks;mask++)if(__builtin_popcountl(mask)==k-1)
            for(slong a=0;a<=M;a++)if(state[mask*(M+1)+a].r)dq_mat_clear(state+mask*(M+1)+a);
    }
    dq_mat_set(out,state+(masks-1)*(M+1)+M);
    dq_mat_clear(state+(masks-1)*(M+1)+M);flint_free(state);
    double elapsed=now()-t;
    #pragma omp atomic update
    dp_time+=elapsed;
}
static void specialize(dq_coeff *q,const dq_coeff *input,dq_coeff value)
{
    slong w=N+1;dq_coeff square;dq_init(square);dq_pow(square,value,2);
    for(slong f=0;f<N;f++) {
        for(slong i=0;i<N-1;i++)for(slong j=i;j<N-1;j++)dq_set(q[(f*N+i)*N+j],input[(f*w+i)*w+j]);
        for(slong i=0;i<N-1;i++) {
            dq_set(q[(f*N+i)*N+N-1],input[(f*w+i)*w+N]);
            dq_addmul(q[(f*N+i)*N+N-1],value,input[(f*w+i)*w+N-1]);
        }
        dq_coeff *c=q+(f*N+N-1)*N+N-1;
        dq_set(*c,input[(f*w+N)*w+N]);
        dq_addmul(*c,value,input[(f*w+N-1)*w+N]);
        dq_addmul(*c,square,input[(f*w+N-1)*w+N-1]);
    }
    dq_clear(square);
}
static void evaluate(dq_mat_t core,const dq_coeff *input,dq_coeff value,uint64_t *basis,uint64_t *basis_y,int first)
{
    dq_coeff *q=dq_vec_new(N*N*N);specialize(q,input,value);
    grade *g=flint_calloc(L+1,sizeof(grade));assert(normal_forms(g,q));
    for(slong j=0;j<g[M].h;j++) {
        uint64_t key=g[M].keys[g[M].free[j]];
        if(first)basis[j]=key;else assert(basis[j]==key);
        key=g[K].keys[g[K].free[j]];
        if(first)basis_y[j]=key;else assert(basis_y[j]==key);
    }
    linear_map *op=recurrence_kernel>=4?NULL:make_linear_maps(g,q);
    if(op && !dq_nonzero(value) && getenv("DRSOLVE_QUOTIENT_PROFILE"))for(slong d=0;d<L;d++) {
        slong nonzero=0,entries=0,single=0,total_rows=0;
        for(slong v=0;v<N;v++)for(slong i=0;i<g[d].h;i++) {
            slong count=0;for(slong j=0;j<g[d+1].h;j++)count+=dq_nonzero(dq_mat_entry(g[d].variables,v,i*g[d+1].h+j));
            single+=count==1;total_rows++;
        }
        for(slong f=0;f<N;f++)for(slong j=0;j<N;j++) {
            const linear_map *m=op+(f*N+j)*L+d;
            for(slong r=0;r<m->y->r;r++)for(slong c=0;c<m->y->c;c++) {
                nonzero+=dq_nonzero(dq_mat_entry(m->y,r,c));entries++;
            }
        }
        fprintf(stderr,"grade %ld: variable one-entry rows %ld/%ld; right map density %.4f\n",d,single,total_rows,(double)nonzero/entries);
    }
    quotient_det(core,g,op,q);
    if(verify_recurrence && recurrence_kernel) {
        dq_mat_t original;dq_mat_init(original,core->r,core->c,prime);
        int saved=recurrence_kernel;recurrence_kernel=0;
        linear_map *reference_op=saved>=3?make_linear_maps(g,q):op;
        quotient_det(original,g,reference_op,q);
        if(saved>=3)clear_linear_maps(reference_op);
        recurrence_kernel=saved;
        assert(dq_mat_equal(original,core));dq_mat_clear(original);
    }
    clear_linear_maps(op);clear_grades(g);flint_free(g);dq_vec_free(q,N*N*N);
}

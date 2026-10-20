/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Shared prime-field general-degree quotient construction.
 * Homogeneous normal forms use full Macaulay RREF at each parameter point.
 * The library adapter uses the grouped Bezoutian recurrence; the standalone
 * experiment also retains term-by-term and full-block comparison paths.
 * No shared leading frames or packed kernels. The library adapter can run
 * independent samples in parallel, with a private context for each worker.
 * Library admission checks representation and allocation arithmetic. Rank/basis rejection
 * returns to its caller; the experiment reports failures as process errors.
 */
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <limits.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly_mat.h>
#include "pml_det.h"

static slong n,max_degree,limit,left_degree,right_degree,socle;
static slong *degrees,*prefix_degree,*hilbert_dims;
static slong expected_degree;
static ulong prime;
static int grouped_kernel=1;
static nmod_t mod;
static double seconds(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
#ifndef DEGREE_QUOTIENT_LIBRARY
static void fail(const char *message) {fprintf(stderr,"General-degree quotient experiment: %s\n",message);exit(2);}
#else
static void fail(const char *message);
#endif
static slong positive_integer(const char *s) {
    char *end;errno=0;long value=strtol(s,&end,10);
    if(errno || !*s || *end || value<=0)fail("expected a positive integer");
    return value;
}
static slong binomial(slong a,slong b) {
    if(b>a-b)b=a-b;
    slong r=1;
    for(slong j=1;j<=b;j++) {
        __uint128_t next=(__uint128_t)r*(a-j+1)/j;
        if(next>LONG_MAX)fail("dimension exceeds machine-word capacity");
        r=(slong)next;
    }
    return r;
}
static slong exponent(uint64_t key,slong v) {return (key>>(8*v))&255;}
static slong key_degree(uint64_t key) {slong d=0;for(slong v=0;v<n;v++)d+=exponent(key,v);return d;}
typedef struct {slong size,h;uint64_t *keys;slong *free; nmod_mat_t reduce;} grade;
/* Dimensions, monomial keys and input coefficients remain read-only during
 * sampling. All mutable normal forms and instrumentation belong to a worker.
 * Recurrence maps/workspaces already have sample-local lifetimes. */
typedef struct {
    grade *grades;
    double setup_seconds,update_seconds;
    slong products;
} degree_worker;
static _Thread_local degree_worker active_worker;
#define spaces (active_worker.grades)
#define operator_setup_seconds (active_worker.setup_seconds)
#define operator_update_seconds (active_worker.update_seconds)
#define matrix_products (active_worker.products)
static void init_dimensions(const char *spec,slong memory_mib) {
    degrees=flint_malloc(n*sizeof(slong));prefix_degree=flint_calloc(n+1,sizeof(slong));
    char *copy=flint_malloc(strlen(spec)+1);strcpy(copy,spec);
    char *cursor=copy;slong count=0;
    do {
        char *comma=strchr(cursor,',');if(comma)*comma='\0';
        if(count==n)fail("too many equation degrees");
        degrees[count++]=positive_integer(cursor);
        if(!comma)break;
        cursor=comma+1;
    } while(1);
    flint_free(copy);
    if(count==1)for(slong i=1;i<n;i++)degrees[i]=degrees[0];
    else if(count!=n)fail("degree list must contain exactly N entries");
    expected_degree=1;
    for(slong i=0;i<n;i++) {
        if(degrees[i]>255 || socle>255-(degrees[i]-1))fail("total Bezoutian degree exceeds packed exponent capacity (255)");
        socle+=degrees[i]-1;prefix_degree[i+1]=socle;
        max_degree=FLINT_MAX(max_degree,degrees[i]);
        if(expected_degree>LONG_MAX/degrees[i])fail("resultant degree exceeds machine-word capacity");
        expected_degree*=degrees[i];
    }
    left_degree=socle/2;right_degree=socle-left_degree;limit=FLINT_MAX(right_degree,max_degree);
    hilbert_dims=flint_calloc(limit+1,sizeof(slong));hilbert_dims[0]=1;
    for(slong i=0;i<n;i++) {
        slong *next=flint_calloc(limit+1,sizeof(slong));
        for(slong j=0;j<=limit;j++)for(slong k=0;k<degrees[i] && j+k<=limit;k++) {
            if(next[j+k]>LONG_MAX-hilbert_dims[j])fail("Hilbert dimension overflow");
            next[j+k]+=hilbert_dims[j];
        }
        memcpy(hilbert_dims,next,(limit+1)*sizeof(slong));flint_free(next);
    }
    /* Conservative storage bound for both comparison paths, their cached
     * multiplication maps, normal forms, interpolation, and Macaulay RREF. */
    long double bytes=0,peak_rel=0,previous=0,peak_states=0;
    slong max_h=0,max_s=0;
    for(slong d=0;d<=limit;d++) {
        slong size=binomial(n+d-1,d),h=hilbert_dims[d];
        bytes+=(long double)size*(sizeof(uint64_t)+(long double)h*sizeof(ulong))+h*sizeof(slong);
        if(d<=right_degree) {
            long double rows=0;
            for(slong f=0;f<n;f++)if(d>=degrees[f])rows+=binomial(n+d-degrees[f]-1,d-degrees[f]);
            if(rows*size>peak_rel)peak_rel=rows*size;
            max_h=FLINT_MAX(max_h,h);max_s=FLINT_MAX(max_s,size);
        }
    }
    for(slong k=0;k<=n;k++) {
        long double layer=0;slong total=prefix_degree[k];
        for(slong a=FLINT_MAX(0,total-right_degree);a<=FLINT_MIN(left_degree,total);a++)
            layer+=(long double)binomial(n+a-1,a)*binomial(n+total-a-1,total-a);
        layer*=binomial(n,k);
        if(previous+layer>peak_states)peak_states=previous+layer;
        previous=layer;
    }
    bytes+=(peak_rel+peak_states+(long double)(socle+3)*max_h*max_h*4)*sizeof(ulong);
    bytes+=(long double)binomial(n+max_degree-1,max_degree-1)*(right_degree+1)*
        ((long double)max_h*max_h*sizeof(ulong)+max_h*sizeof(slong)+128);
    bytes+=(long double)binomial(n+max_degree,max_degree)*(n*sizeof(ulong)+32);
    bytes+=(long double)(1L<<n)*(left_degree+1)*sizeof(nmod_mat_struct);
    if(grouped_kernel) {
        long double peak_operators=0;
        for(slong f=0;f<n;f++) {
            long double layer=0;slong total=prefix_degree[f+1];
            for(slong j=0;j<n;j++)
                for(slong a=FLINT_MAX(0,total-right_degree);a<=FLINT_MIN(left_degree,total);a++) {
                    slong b=total-a;
                    /* Term arrays are allocated at the ungrouped entry size. */
                    layer+=(long double)degrees[f]*binomial(n+degrees[f]-1,degrees[f])*256;
                    for(slong ax=0;ax<=FLINT_MIN(a,degrees[f]-1);ax++) {
                        slong by=degrees[f]-1-ax;if(by>b)continue;
                        slong groups=FLINT_MIN(binomial(j+ax,ax),binomial(n-j+by-1,by));
                        layer+=(long double)groups*((long double)hilbert_dims[a]*hilbert_dims[a-ax]+
                                (long double)hilbert_dims[b-by]*hilbert_dims[b])*sizeof(ulong);
                    }
                }
            if(layer>peak_operators)peak_operators=layer;
        }
        bytes+=2*(long double)max_h*max_s*sizeof(ulong);
        bytes+=peak_operators+2*(long double)(right_degree+1)*(right_degree+1)*max_h*max_h*sizeof(ulong);
    }
    if(bytes>(long double)SIZE_MAX)fail("workspace dimensions exceed machine address capacity");
    if(memory_mib>0 && bytes>(long double)memory_mib*1024*1024) {
        fprintf(stderr,"Estimated storage bound %.1Lf MiB exceeds budget %ld MiB.\n",bytes/1048576,memory_mib);
        fail("increase the optional memory-MiB argument or reduce the problem size");
    }
}
static void enumerate(grade *g,slong v,slong remaining,uint64_t key,slong *count) {
    if(v==n-2){g->keys[(*count)++]=key|((uint64_t)remaining<<(8*v));return;}
    for(slong a=remaining;a>=0;a--)enumerate(g,v+1,remaining-a,key|((uint64_t)a<<(8*v)),count);
}
static slong index_of(const grade *g,uint64_t key) {
    for(slong i=0;i<g->size;i++)if(g->keys[i]==key)return i;
    fprintf(stderr,"missing monomial %llu\n",(unsigned long long)key);abort();
}
static void init_spaces(void) {
    spaces=flint_calloc(limit+1,sizeof(*spaces));
    for(slong d=0;d<=limit;d++) {
        grade *g=spaces+d;g->size=binomial(n+d-1,d);g->h=hilbert_dims[d];
        g->keys=flint_malloc(g->size*sizeof(uint64_t));g->free=flint_malloc(g->h*sizeof(slong));
        slong count=0;
        for(slong z=0;z<=d;z++)enumerate(g,0,d-z,(uint64_t)z<<(8*(n-1)),&count);
        assert(count==g->size);nmod_mat_init(g->reduce,g->size,g->h,prime);
    }
}
static int normal_forms(ulong **q,int first) {
    for(slong d=0;d<=right_degree;d++) {
        grade *g=spaces+d;nmod_mat_zero(g->reduce);slong rows=0;
        for(slong f=0;f<n;f++)if(d>=degrees[f])rows+=spaces[d-degrees[f]].size;
        if(!rows) {
            for(slong i=0;i<g->size;i++){g->free[i]=i;nmod_mat_entry(g->reduce,i,i)=1;}
            continue;
        }
        nmod_mat_t rel;nmod_mat_init(rel,rows,g->size,prime);slong offset=0;
        for(slong f=0;f<n;f++)if(d>=degrees[f]) {
            grade *mult=spaces+d-degrees[f],*terms=spaces+degrees[f];
            for(slong a=0;a<mult->size;a++)for(slong t=0;t<terms->size;t++)
                nmod_mat_entry(rel,offset+a,index_of(g,mult->keys[a]+terms->keys[t]))=q[f][t];
            offset+=mult->size;
        }
        slong rank=nmod_mat_rref(rel);
        if(rank!=g->size-g->h){nmod_mat_clear(rel);return 0;}
        slong *pivot=flint_malloc(rank*sizeof(slong));slong r=0,free_count=0;
        for(slong c=0;c<g->size;c++) {
            if(r<rank && nmod_mat_entry(rel,r,c))pivot[r++]=c;
            else {
                if(!first && g->free[free_count]!=c){flint_free(pivot);nmod_mat_clear(rel);return 0;}
                g->free[free_count++]=c;
            }
        }
        assert(free_count==g->h);
        for(slong j=0;j<g->h;j++) {
            nmod_mat_entry(g->reduce,g->free[j],j)=1;
            for(slong i=0;i<rank;i++)nmod_mat_entry(g->reduce,pivot[i],j)=nmod_neg(nmod_mat_entry(rel,i,g->free[j]),mod);
        }
        flint_free(pivot);nmod_mat_clear(rel);
    }
    return 1;
}
/* Dense affine input; each equation retains its own total-degree bound. */
typedef struct {uint64_t key;ulong *coeff;} input_term;
static input_term *input;static slong input_count;
static void input_enum(slong v,slong remaining,uint64_t key,flint_rand_t rng) {
    if(v==n) {
        input[input_count].key=key;input[input_count].coeff=flint_malloc(n*sizeof(ulong));
        for(slong f=0;f<n;f++)input[input_count].coeff[f]=key_degree(key)<=degrees[f]?n_randint(rng,prime):0;
        input_count++;return;
    }
    for(slong a=0;a<=remaining;a++)input_enum(v+1,remaining-a,key|((uint64_t)a<<(8*v)),rng);
}
static void specialize(ulong **q,ulong value) {
    for(slong f=0;f<n;f++)memset(q[f],0,spaces[degrees[f]].size*sizeof(ulong));
    for(slong t=0;t<input_count;t++) {
        uint64_t key=input[t].key;slong pexp=exponent(key,n-1);
        uint64_t pure=key&((UINT64_C(1)<<(8*(n-1)))-1);
        ulong factor=nmod_pow_ui(value,pexp,mod);
        for(slong f=0;f<n;f++)if(input[t].coeff[f]) {
            slong zexp=degrees[f]-key_degree(pure);
            uint64_t homogeneous=pure|((uint64_t)zexp<<(8*(n-1)));
            slong col=index_of(spaces+degrees[f],homogeneous);
            q[f][col]=nmod_add(q[f][col],nmod_mul(input[t].coeff[f],factor,mod),mod);
        }
    }
}
typedef struct {uint64_t x,y;ulong c;} pair_term;
typedef struct {slong count;pair_term *terms;} entry;
static entry *difference_matrix(ulong **q) {
    entry *out=flint_calloc(n*n,sizeof(*out));
    for(slong f=0;f<n;f++)for(slong j=0;j<n;j++) {
        grade *g=spaces+degrees[f];
        entry *e=out+f*n+j;e->terms=flint_malloc(g->size*degrees[f]*sizeof(pair_term));
        for(slong t=0;t<g->size;t++)if(q[f][t]) {
            uint64_t key=g->keys[t];slong power=exponent(key,j);
            for(slong a=0;a<power;a++) {
                pair_term p={0,0,q[f][t]};
                for(slong v=0;v<n;v++) {
                    if(v<j)p.x|=(uint64_t)exponent(key,v)<<(8*v);
                    else if(v>j)p.y|=(uint64_t)exponent(key,v)<<(8*v);
                    else {p.x|=(uint64_t)a<<(8*v);p.y|=(uint64_t)(power-1-a)<<(8*v);}
                }
                e->terms[e->count++]=p;
            }
        }
    }
    return out;
}
/* Maps from degree a to a+deg(monomial); all cached per sample. */
typedef struct map_node {slong a;uint64_t key;slong *rows;nmod_mat_t matrix;struct map_node *next;} map_node;
static map_node *monomial_rows(map_node **head,slong a,uint64_t key) {
    for(map_node *p=*head;p;p=p->next)if(p->a==a && p->key==key)return p;
    grade *src=spaces+a,*dst=spaces+a+key_degree(key);
    map_node *p=flint_calloc(1,sizeof(*p));p->a=a;p->key=key;p->next=*head;*head=p;
    p->rows=flint_malloc(src->h*sizeof(slong));
    for(slong i=0;i<src->h;i++)p->rows[i]=index_of(dst,src->keys[src->free[i]]+key);
    return p;
}
static const nmod_mat_struct *multiplication(map_node **head,slong a,uint64_t key) {
    map_node *p=monomial_rows(head,a,key);
    grade *src=spaces+a,*dst=spaces+a+key_degree(key);
    if(!p->matrix->r) {
        nmod_mat_init(p->matrix,dst->h,src->h,prime);
        for(slong i=0;i<src->h;i++)for(slong j=0;j<dst->h;j++)
            nmod_mat_entry(p->matrix,j,i)=nmod_mat_entry(dst->reduce,p->rows[i],j);
    }
    return p->matrix;
}
#include "degree_quotient_grouped.h"

static void bezout_core(nmod_mat_t out,const entry *entries,int reduced) {
    if(reduced){operator_setup_seconds=operator_update_seconds=0;matrix_products=0;}
    slong masks=1L<<n,stride=left_degree+1;
    nmod_mat_struct *state=flint_calloc(masks*stride,sizeof(*state));map_node *maps=NULL;
    nmod_mat_init(state,1,1,prime);nmod_mat_entry(state,0,0)=1;
    for(slong k=1;k<=n;k++) {
        grouped_operator *operators=NULL;
        grouped_workspace *temps=NULL,*products=NULL;
        for(slong mask=1;mask<masks;mask++)if(__builtin_popcountl(mask)==k) {
            slong total=prefix_degree[k];
            for(slong a=FLINT_MAX(0,total-right_degree);a<=FLINT_MIN(left_degree,total);a++) {
                slong b=total-a;grade *ga=spaces+a,*gb=spaces+b;
                nmod_mat_struct *dst=state+mask*stride+a;
                nmod_mat_init(dst,reduced?ga->h:ga->size,reduced?gb->h:gb->size,prime);
                slong pos=0;
                for(slong j=0;j<n;j++)if(mask&(1L<<j)) {
                    slong child=mask^(1L<<j);int neg=(k-1+pos++)&1;const entry *e=entries+(k-1)*n+j;
                    if(reduced && grouped_kernel) {
                        double stage=seconds();
                        grouped_operator *op=grouped_build(&operators,&maps,e,a,b);
                        operator_setup_seconds+=seconds()-stage;stage=seconds();
                        grouped_apply(dst,op,state+child*stride,neg,&temps,&products);
                        operator_update_seconds+=seconds()-stage;
                        continue;
                    }
                    for(slong t=0;t<e->count;t++) {
                        pair_term p=e->terms[t];slong ax=key_degree(p.x),by=key_degree(p.y);
                        if(a<ax || b<by)continue;
                        slong ca=a-ax,cb=b-by;nmod_mat_struct *src=state+child*stride+ca;
                        if(!src->r)continue;
                        ulong scalar=neg?nmod_neg(p.c,mod):p.c;
                        if(reduced) {
                            const nmod_mat_struct *mx=multiplication(&maps,ca,p.x),*my=multiplication(&maps,cb,p.y);
                            nmod_mat_t temp,right,product;
                            nmod_mat_init(temp,mx->r,src->c,prime);nmod_mat_init(right,my->c,my->r,prime);
                            nmod_mat_init(product,dst->r,dst->c,prime);nmod_mat_transpose(right,my);
                            nmod_mat_mul(temp,mx,src);nmod_mat_mul(product,temp,right);matrix_products+=2;
                            for(slong r=0;r<dst->r;r++)for(slong c=0;c<dst->c;c++)
                                nmod_mat_entry(dst,r,c)=nmod_add(nmod_mat_entry(dst,r,c),nmod_mul(scalar,nmod_mat_entry(product,r,c),mod),mod);
                            nmod_mat_clear(temp);nmod_mat_clear(right);nmod_mat_clear(product);
                        } else {
                            for(slong r=0;r<src->r;r++) {
                                slong rr=index_of(ga,spaces[ca].keys[r]+p.x);
                                for(slong c=0;c<src->c;c++) {
                                    slong cc=index_of(gb,spaces[cb].keys[c]+p.y);
                                    nmod_mat_entry(dst,rr,cc)=nmod_add(nmod_mat_entry(dst,rr,cc),nmod_mul(scalar,nmod_mat_entry(src,r,c),mod),mod);
                                }
                            }
                        }
                    }
                }
            }
        }
        grouped_clear(operators);grouped_scratch_clear(temps);grouped_scratch_clear(products);
        for(slong mask=0;mask<masks;mask++)if(__builtin_popcountl(mask)==k-1)
            for(slong a=0;a<stride;a++)if(state[mask*stride+a].r){nmod_mat_clear(state+mask*stride+a);memset(state+mask*stride+a,0,sizeof(*state));}
    }
    nmod_mat_struct *last=state+(masks-1)*stride+left_degree;
    if(reduced)nmod_mat_set(out,last);
    else {
        nmod_mat_t transpose,temp;nmod_mat_init(transpose,spaces[left_degree].h,spaces[left_degree].size,prime);
        nmod_mat_init(temp,transpose->r,last->c,prime);nmod_mat_transpose(transpose,spaces[left_degree].reduce);
        nmod_mat_mul(temp,transpose,last);nmod_mat_mul(out,temp,spaces[right_degree].reduce);
        nmod_mat_clear(transpose);nmod_mat_clear(temp);
    }
    nmod_mat_clear(last);flint_free(state);
    while(maps){map_node *next=maps->next;if(maps->matrix->r)nmod_mat_clear(maps->matrix);flint_free(maps->rows);flint_free(maps);maps=next;}
}

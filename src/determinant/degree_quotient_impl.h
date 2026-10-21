/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Shared finite-field general-degree quotient construction.
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
#ifndef DEGREE_QUOTIENT_EXTENSION
#include "degree_quotient_prime.h"
#endif
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
#include "degree_quotient_keys.h"
typedef struct {slong size,h;dq_key *keys;slong *free; dq_mat_t reduce;} grade;
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
    if(n<2 || n>=FLINT_BITS-1)fail("subset workspace exceeds machine address capacity");
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
        if(socle>LONG_MAX-n-degrees[i])fail("total Bezoutian degree exceeds machine address capacity");
        socle+=degrees[i]-1;prefix_degree[i+1]=socle;
        max_degree=FLINT_MAX(max_degree,degrees[i]);
        if(expected_degree>LONG_MAX/degrees[i])fail("resultant degree exceeds machine-word capacity");
        expected_degree*=degrees[i];
    }
    left_degree=socle/2;right_degree=socle-left_degree;limit=FLINT_MAX(right_degree,max_degree);
    slong key_count=binomial(n+limit,limit);
    if((ulong)key_count>SIZE_MAX/((n+2)*sizeof(slong)))
        fail("monomial tables exceed machine address capacity");
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
        bytes+=(long double)size*((n+2)*sizeof(slong)+(long double)h*DQ_COEFF_BYTES)+h*sizeof(slong);
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
    bytes+=(peak_rel+peak_states+(long double)(socle+3)*max_h*max_h*4)*DQ_COEFF_BYTES;
    bytes+=(long double)binomial(n+max_degree-1,max_degree-1)*(right_degree+1)*
        ((long double)max_h*max_h*DQ_COEFF_BYTES+max_h*sizeof(slong)+128);
    bytes+=(long double)binomial(n+max_degree,max_degree)*(n*DQ_COEFF_BYTES+32);
    bytes+=(long double)(1L<<n)*(left_degree+1)*sizeof(dq_mat_struct);
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
                                (long double)hilbert_dims[b-by]*hilbert_dims[b])*DQ_COEFF_BYTES;
                    }
                }
            if(layer>peak_operators)peak_operators=layer;
        }
        bytes+=2*(long double)max_h*max_s*DQ_COEFF_BYTES;
        bytes+=peak_operators+2*(long double)(right_degree+1)*(right_degree+1)*max_h*max_h*DQ_COEFF_BYTES;
    }
    if(bytes>(long double)SIZE_MAX)fail("workspace dimensions exceed machine address capacity");
    if(memory_mib>0 && bytes>(long double)memory_mib*1024*1024) {
        fprintf(stderr,"Estimated storage bound %.1Lf MiB exceeds budget %ld MiB.\n",bytes/1048576,memory_mib);
        fail("increase the optional memory-MiB argument or reduce the problem size");
    }
}
static void enumerate(grade *g,slong d,slong v,slong remaining,slong *e,slong *count) {
    if(v==n-2) {
        e[v]=remaining;dq_key key=key_offsets[d]+*count;
        g->keys[(*count)++]=key;key_degrees[key]=d;
        memcpy(key_exponents+key*n,e,n*sizeof(slong));return;
    }
    for(slong a=remaining;a>=0;a--){e[v]=a;enumerate(g,d,v+1,remaining-a,e,count);}
}
static slong index_of(const grade *g,dq_key key) {
    slong index=key-g->keys[0];assert(index>=0 && index<g->size);return index;
}
static void init_spaces(void) {
    key_tables_init();spaces=flint_calloc(limit+1,sizeof(*spaces));
    for(slong d=0;d<=limit;d++) {
        grade *g=spaces+d;g->size=binomial(n+d-1,d);g->h=hilbert_dims[d];
        g->keys=flint_malloc(g->size*sizeof(dq_key));g->free=flint_malloc(g->h*sizeof(slong));
        slong count=0,e[n];
        for(slong z=0;z<=d;z++){e[n-1]=z;enumerate(g,d,0,d-z,e,&count);}
        assert(count==g->size);dq_mat_init(g->reduce,g->size,g->h,prime);
    }
}
static int normal_forms(dq_coeff **q,int first) {
    for(slong d=0;d<=right_degree;d++) {
        grade *g=spaces+d;dq_mat_zero(g->reduce);slong rows=0;
        for(slong f=0;f<n;f++)if(d>=degrees[f])rows+=spaces[d-degrees[f]].size;
        if(!rows) {
            for(slong i=0;i<g->size;i++){g->free[i]=i;dq_one(dq_mat_entry(g->reduce,i,i));}
            continue;
        }
        dq_mat_t rel;dq_mat_init(rel,rows,g->size,prime);slong offset=0;
        for(slong f=0;f<n;f++)if(d>=degrees[f]) {
            grade *mult=spaces+d-degrees[f],*terms=spaces+degrees[f];
            for(slong a=0;a<mult->size;a++)for(slong t=0;t<terms->size;t++)
                dq_set(dq_mat_entry(rel,offset+a,index_of(g,key_add(mult->keys[a],terms->keys[t]))),q[f][t]);
            offset+=mult->size;
        }
        slong rank=dq_mat_rref(rel);
        if(rank!=g->size-g->h){dq_mat_clear(rel);return 0;}
        slong *pivot=flint_malloc(rank*sizeof(slong));slong r=0,free_count=0;
        for(slong c=0;c<g->size;c++) {
            if(r<rank && dq_nonzero(dq_mat_entry(rel,r,c)))pivot[r++]=c;
            else {
                if(!first && g->free[free_count]!=c){flint_free(pivot);dq_mat_clear(rel);return 0;}
                g->free[free_count++]=c;
            }
        }
        assert(free_count==g->h);
        for(slong j=0;j<g->h;j++) {
            dq_one(dq_mat_entry(g->reduce,g->free[j],j));
            for(slong i=0;i<rank;i++)dq_neg(dq_mat_entry(g->reduce,pivot[i],j),dq_mat_entry(rel,i,g->free[j]));
        }
        flint_free(pivot);dq_mat_clear(rel);
    }
    return 1;
}
/* Dense affine input; each equation retains its own total-degree bound. */
typedef struct {dq_key key;dq_coeff *coeff;} input_term;
static input_term *input;static slong input_count;
#ifndef DEGREE_QUOTIENT_EXTENSION
static void input_enum(slong v,slong remaining,dq_key key,flint_rand_t rng) {
    if(v==n) {
        input[input_count].key=key;input[input_count].coeff=flint_malloc(n*DQ_COEFF_BYTES);
        for(slong f=0;f<n;f++)input[input_count].coeff[f]=key_degree(key)<=degrees[f]?n_randint(rng,prime):0;
        input_count++;return;
    }
    for(slong a=0;a<=remaining;a++)input_enum(v+1,remaining-a,key_set(key,v,a),rng);
}
#endif
static void specialize(dq_coeff **q,dq_coeff value) {
    for(slong f=0;f<n;f++)dq_vec_zero(q[f],spaces[degrees[f]].size);
    for(slong t=0;t<input_count;t++) {
        dq_key key=input[t].key;slong pexp=exponent(key,n-1);
        dq_key pure=key_set(key,n-1,0);
        dq_coeff factor;dq_init(factor);dq_pow(factor,value,pexp);
        for(slong f=0;f<n;f++)if(dq_nonzero(input[t].coeff[f])) {
            slong zexp=degrees[f]-key_degree(pure);
            dq_key homogeneous=key_set(pure,n-1,zexp);
            slong col=index_of(spaces+degrees[f],homogeneous);
            dq_addmul(q[f][col],input[t].coeff[f],factor);
        }
        dq_clear(factor);
    }
}
/* Pair coefficients are borrowed from q for this sample. */
typedef struct {dq_key x,y;dq_coeff c;} pair_term;
typedef struct {slong count;pair_term *terms;} entry;
static entry *difference_matrix(dq_coeff **q) {
    entry *out=flint_calloc(n*n,sizeof(*out));
    for(slong f=0;f<n;f++)for(slong j=0;j<n;j++) {
        grade *g=spaces+degrees[f];
        entry *e=out+f*n+j;e->terms=flint_malloc(g->size*degrees[f]*sizeof(pair_term));
        for(slong t=0;t<g->size;t++)if(dq_nonzero(q[f][t])) {
            dq_key key=g->keys[t];slong power=exponent(key,j);
            for(slong a=0;a<power;a++) {
                pair_term p={0,0,q[f][t]};
                for(slong v=0;v<n;v++) {
                    if(v<j)p.x=key_set(p.x,v,exponent(key,v));
                    else if(v>j)p.y=key_set(p.y,v,exponent(key,v));
                    else {p.x=key_set(p.x,v,a);p.y=key_set(p.y,v,power-1-a);}
                }
                e->terms[e->count++]=p;
            }
        }
    }
    return out;
}
/* Maps from degree a to a+deg(monomial); all cached per sample. */
typedef struct map_node {slong a;dq_key key;slong *rows;dq_mat_t matrix;struct map_node *next;} map_node;
static map_node *monomial_rows(map_node **head,slong a,dq_key key) {
    for(map_node *p=*head;p;p=p->next)if(p->a==a && p->key==key)return p;
    grade *src=spaces+a,*dst=spaces+a+key_degree(key);
    map_node *p=flint_calloc(1,sizeof(*p));p->a=a;p->key=key;p->next=*head;*head=p;
    p->rows=flint_malloc(src->h*sizeof(slong));
    for(slong i=0;i<src->h;i++)p->rows[i]=index_of(dst,key_add(src->keys[src->free[i]],key));
    return p;
}
static const dq_mat_struct *multiplication(map_node **head,slong a,dq_key key) {
    map_node *p=monomial_rows(head,a,key);
    grade *src=spaces+a,*dst=spaces+a+key_degree(key);
    if(!p->matrix->r) {
        dq_mat_init(p->matrix,dst->h,src->h,prime);
        for(slong i=0;i<src->h;i++)for(slong j=0;j<dst->h;j++)
            dq_set(dq_mat_entry(p->matrix,j,i),dq_mat_entry(dst->reduce,p->rows[i],j));
    }
    return p->matrix;
}
#include "degree_quotient_grouped.h"

static void bezout_core(dq_mat_t out,const entry *entries,int reduced) {
    if(reduced){operator_setup_seconds=operator_update_seconds=0;matrix_products=0;}
    slong masks=1L<<n,stride=left_degree+1;
    dq_mat_struct *state=flint_calloc(masks*stride,sizeof(*state));map_node *maps=NULL;
    dq_mat_init(state,1,1,prime);dq_one(dq_mat_entry(state,0,0));
    for(slong k=1;k<=n;k++) {
        grouped_operator *operators=NULL;
        grouped_workspace *temps=NULL,*products=NULL;
        for(slong mask=1;mask<masks;mask++)if(__builtin_popcountl(mask)==k) {
            slong total=prefix_degree[k];
            for(slong a=FLINT_MAX(0,total-right_degree);a<=FLINT_MIN(left_degree,total);a++) {
                slong b=total-a;grade *ga=spaces+a,*gb=spaces+b;
                dq_mat_struct *dst=state+mask*stride+a;
                dq_mat_init(dst,reduced?ga->h:ga->size,reduced?gb->h:gb->size,prime);
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
                        slong ca=a-ax,cb=b-by;dq_mat_struct *src=state+child*stride+ca;
                        if(!src->r)continue;
                        dq_coeff scalar;dq_init(scalar);
                        if(neg)dq_neg(scalar,p.c);else dq_set(scalar,p.c);
                        if(reduced) {
                            const dq_mat_struct *mx=multiplication(&maps,ca,p.x),*my=multiplication(&maps,cb,p.y);
                            dq_mat_t temp,right,product;
                            dq_mat_init(temp,mx->r,src->c,prime);dq_mat_init(right,my->c,my->r,prime);
                            dq_mat_init(product,dst->r,dst->c,prime);dq_mat_transpose(right,my);
                            dq_mat_mul(temp,mx,src);dq_mat_mul(product,temp,right);matrix_products+=2;
                            for(slong r=0;r<dst->r;r++)for(slong c=0;c<dst->c;c++)
                                dq_addmul(dq_mat_entry(dst,r,c),scalar,dq_mat_entry(product,r,c));
                            dq_mat_clear(temp);dq_mat_clear(right);dq_mat_clear(product);
                        } else {
                            for(slong r=0;r<src->r;r++) {
                                slong rr=index_of(ga,key_add(spaces[ca].keys[r],p.x));
                                for(slong c=0;c<src->c;c++) {
                                    slong cc=index_of(gb,key_add(spaces[cb].keys[c],p.y));
                                    dq_addmul(dq_mat_entry(dst,rr,cc),scalar,dq_mat_entry(src,r,c));
                                }
                            }
                        }
                        dq_clear(scalar);
                    }
                }
            }
        }
        grouped_clear(operators);grouped_scratch_clear(temps);grouped_scratch_clear(products);
        for(slong mask=0;mask<masks;mask++)if(__builtin_popcountl(mask)==k-1)
            for(slong a=0;a<stride;a++)if(state[mask*stride+a].r){dq_mat_clear(state+mask*stride+a);memset(state+mask*stride+a,0,sizeof(*state));}
    }
    dq_mat_struct *last=state+(masks-1)*stride+left_degree;
    if(reduced)dq_mat_set(out,last);
    else {
        dq_mat_t transpose,temp;dq_mat_init(transpose,spaces[left_degree].h,spaces[left_degree].size,prime);
        dq_mat_init(temp,transpose->r,last->c,prime);dq_mat_transpose(transpose,spaces[left_degree].reduce);
        dq_mat_mul(temp,transpose,last);dq_mat_mul(out,temp,spaces[right_degree].reduce);
        dq_mat_clear(transpose);dq_mat_clear(temp);
    }
    dq_mat_clear(last);flint_free(state);
    while(maps){map_node *next=maps->next;if(maps->matrix->r)dq_mat_clear(maps->matrix);flint_free(maps->rows);flint_free(maps);maps=next;}
}

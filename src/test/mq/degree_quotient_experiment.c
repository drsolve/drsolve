/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Standalone general-degree experiment, not a production backend.
 * Build homogeneous normal forms by full Macaulay RREF at each parameter
 * point. Compare reduction during the Bezoutian subset recurrence against
 * constructing its full middle coefficient block and reducing afterwards.
 * No shared leading frames, packed kernels, or sample parallelism.
 * Usage: executable N degrees seed output-prefix [prime] [memory-MiB].
 * degrees is one positive integer or an N-entry comma-separated list.
 * N=2..8; monomial exponents must fit in 8-bit packed slots.
 * Outputs input .dr, monic determinant coefficients .coeff, and JSON timings.
 * Generic ranks, fixed bases, extra interpolation points and full-vs-reduced
 * core equality are checked. Rejected ranks/bases and resource limits report an error and exit.
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
static nmod_t mod;
static double seconds(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void fail(const char *message) {fprintf(stderr,"General-degree quotient experiment: %s\n",message);exit(2);}
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
static grade *spaces;
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
    slong max_h=0;
    for(slong d=0;d<=limit;d++) {
        slong size=binomial(n+d-1,d),h=hilbert_dims[d];
        bytes+=(long double)size*(sizeof(uint64_t)+(long double)h*sizeof(ulong))+h*sizeof(slong);
        if(d<=right_degree) {
            long double rows=0;
            for(slong f=0;f<n;f++)if(d>=degrees[f])rows+=binomial(n+d-degrees[f]-1,d-degrees[f]);
            if(rows*size>peak_rel)peak_rel=rows*size;
            max_h=FLINT_MAX(max_h,h);
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
    bytes+=(long double)binomial(n+max_degree-1,max_degree-1)*(right_degree+1)*max_h*max_h*sizeof(ulong);
    bytes+=(long double)binomial(n+max_degree,max_degree)*(n*sizeof(ulong)+32);
    bytes+=(long double)(1L<<n)*(left_degree+1)*sizeof(nmod_mat_struct);
    if(bytes>(long double)memory_mib*1024*1024) {
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
static void normal_forms(ulong **q,int first) {
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
        if(rank!=g->size-g->h)fail("non-generic Macaulay rank; this sample is unsupported");
        slong *pivot=flint_malloc(rank*sizeof(slong));slong r=0,free_count=0;
        for(slong c=0;c<g->size;c++) {
            if(r<rank && nmod_mat_entry(rel,r,c))pivot[r++]=c;
            else {
                if(!first && g->free[free_count]!=c)fail("quotient basis changed between samples");
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
typedef struct map_node {slong a;uint64_t key;nmod_mat_t matrix;struct map_node *next;} map_node;
static const nmod_mat_struct *multiplication(map_node **head,slong a,uint64_t key) {
    for(map_node *p=*head;p;p=p->next)if(p->a==a && p->key==key)return p->matrix;
    grade *src=spaces+a,*dst=spaces+a+key_degree(key);
    map_node *p=flint_malloc(sizeof(*p));p->a=a;p->key=key;p->next=*head;*head=p;
    nmod_mat_init(p->matrix,dst->h,src->h,prime);
    for(slong i=0;i<src->h;i++) {
        slong row=index_of(dst,src->keys[src->free[i]]+key);
        for(slong j=0;j<dst->h;j++)nmod_mat_entry(p->matrix,j,i)=nmod_mat_entry(dst->reduce,row,j);
    }
    return p->matrix;
}
static void bezout_core(nmod_mat_t out,const entry *entries,int reduced) {
    slong masks=1L<<n,stride=left_degree+1;
    nmod_mat_struct *state=flint_calloc(masks*stride,sizeof(*state));map_node *maps=NULL;
    nmod_mat_init(state,1,1,prime);nmod_mat_entry(state,0,0)=1;
    for(slong k=1;k<=n;k++) {
        for(slong mask=1;mask<masks;mask++)if(__builtin_popcountl(mask)==k) {
            slong total=prefix_degree[k];
            for(slong a=FLINT_MAX(0,total-right_degree);a<=FLINT_MIN(left_degree,total);a++) {
                slong b=total-a;grade *ga=spaces+a,*gb=spaces+b;
                nmod_mat_struct *dst=state+mask*stride+a;
                nmod_mat_init(dst,reduced?ga->h:ga->size,reduced?gb->h:gb->size,prime);
                slong pos=0;
                for(slong j=0;j<n;j++)if(mask&(1L<<j)) {
                    slong child=mask^(1L<<j);int neg=(k-1+pos++)&1;const entry *e=entries+(k-1)*n+j;
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
                            nmod_mat_mul(temp,mx,src);nmod_mat_mul(product,temp,right);
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
    while(maps){map_node *next=maps->next;nmod_mat_clear(maps->matrix);flint_free(maps);maps=next;}
}
static void write_input(const char *prefix) {
    char path[1024];snprintf(path,sizeof(path),"%s.dr",prefix);FILE *f=fopen(path,"w");assert(f);
    for(slong v=0;v<n-1;v++)fprintf(f,"%sx%ld",v?",":"",v);fprintf(f,"\n%lu\n",prime);
    for(slong eq=0;eq<n;eq++) {
        int first=1;
        for(slong t=0;t<input_count;t++)if(input[t].coeff[eq]) {
            fprintf(f,"%s%lu",first?"":"+",input[t].coeff[eq]);first=0;
            for(slong v=0;v<n;v++)if(exponent(input[t].key,v))fprintf(f,"*x%ld^%ld",v,exponent(input[t].key,v));
        }
        fprintf(f,"\n");
    }
    fclose(f);
}
int main(int argc,char **argv) {
    if(argc<5 || argc>7)fail("usage: N degrees seed output-prefix [prime] [memory-MiB]");
    n=positive_integer(argv[1]);
    if(n<2 || n>8)fail("N must be between 2 and 8 for packed monomial keys");
    char *end;errno=0;ulong seed=strtoul(argv[3],&end,10);
    if(errno || !*argv[3] || *end || argv[3][0]=='-')fail("invalid seed");
    prime=argc>=6?(ulong)positive_integer(argv[5]):65537;
    if(!n_is_prime(prime))fail("coefficient modulus must be prime");
    init_dimensions(argv[2],argc>=7?positive_integer(argv[6]):512);
    slong samples=socle+3;
    if(prime<=(ulong)(samples+2))fail("prime too small for distinct interpolation and validation points");
    nmod_init(&mod,prime);flint_set_num_threads(1);
    init_spaces();input=flint_calloc(binomial(n+max_degree,max_degree),sizeof(*input));
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,seed,seed+1);input_enum(0,max_degree,0,rng);write_input(argv[4]);
    slong h=spaces[left_degree].h;assert(h==spaces[right_degree].h);
    ulong **q=flint_malloc(n*sizeof(*q));
    for(slong f=0;f<n;f++)q[f]=flint_malloc(spaces[degrees[f]].size*sizeof(ulong));
    ulong *xs=flint_malloc(samples*sizeof(ulong)),*ys=flint_malloc(h*h*samples*sizeof(ulong));
    double normal_seconds=0,reduced_seconds=0,direct_seconds=0;
    nmod_mat_t core,direct;nmod_mat_init(core,h,h,prime);nmod_mat_init(direct,h,h,prime);
    for(slong s=0;s<samples+2;s++) {
        double t=seconds();specialize(q,s);normal_forms(q,s==0);double normal=seconds()-t;
        entry *entries=difference_matrix(q);
        t=seconds();bezout_core(core,entries,1);double reduced_time=seconds()-t;
        t=seconds();bezout_core(direct,entries,0);double direct_time=seconds()-t;
        assert(nmod_mat_equal(core,direct));
        if(s<samples) {
            normal_seconds+=normal;reduced_seconds+=reduced_time;direct_seconds+=direct_time;xs[s]=s;
            for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)ys[(i*h+j)*samples+s]=nmod_mat_entry(core,i,j);
        }
        for(slong i=0;i<n*n;i++)flint_free(entries[i].terms);flint_free(entries);
        if(s>=samples) {
            nmod_poly_t p;nmod_poly_init(p,prime);
            for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
                nmod_poly_interpolate_nmod_vec(p,xs,ys+(i*h+j)*samples,samples);
                if(nmod_poly_evaluate_nmod(p,s)!=nmod_mat_entry(core,i,j))fail("interpolation failed at an extra validation point");
            }
            nmod_poly_clear(p);
        }
    }
    double t=seconds();nmod_poly_mat_t matrix;nmod_poly_mat_init(matrix,h,h,prime);
    for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
        nmod_poly_struct *entry=nmod_poly_mat_entry(matrix,i,j);
        nmod_poly_interpolate_nmod_vec(entry,xs,ys+(i*h+j)*samples,samples);
        slong wi=left_degree-exponent(spaces[left_degree].keys[spaces[left_degree].free[i]],n-1);
        slong wj=right_degree-exponent(spaces[right_degree].keys[spaces[right_degree].free[j]],n-1);
        if(nmod_poly_degree(entry)>socle+1-wi-wj)fail("interpolated core violates the weighted degree bound");
    }
    double interpolation=seconds()-t;t=seconds();nmod_poly_t det;nmod_poly_init(det,prime);
    if(!nmod_poly_mat_det_hnf(det,matrix))fail("HNF determinant failed");
    if(nmod_poly_degree(det)!=expected_degree)fail("resultant degree differs from the product of equation degrees");nmod_poly_make_monic(det,det);double determinant=seconds()-t;
    char path[1024];snprintf(path,sizeof(path),"%s.coeff",argv[4]);FILE *f=fopen(path,"w");assert(f);
    for(slong i=0;i<det->length;i++)fprintf(f,"%lu%c",det->coeffs[i],i+1==det->length?'\n':' ');fclose(f);
    printf("{\"degrees\":[");
    for(slong f=0;f<n;f++)printf("%s%ld",f?",":"",degrees[f]);
    printf("],\"n\":%ld,\"degree\":%ld,\"seed\":%lu,\"core\":%ld,\"full_rows\":%ld,\"full_cols\":%ld,\"samples\":%ld,\"normal\":%.9f,\"reduced\":%.9f,\"direct\":%.9f,\"interpolation\":%.9f,\"determinant\":%.9f,\"resultant_degree\":%ld}\n",n,max_degree,seed,h,spaces[left_degree].size,spaces[right_degree].size,samples,normal_seconds,reduced_seconds,direct_seconds,interpolation,determinant,nmod_poly_degree(det));
    nmod_poly_clear(det);nmod_poly_mat_clear(matrix);nmod_mat_clear(core);nmod_mat_clear(direct);
    for(slong f=0;f<n;f++)flint_free(q[f]);
    flint_free(q);flint_free(xs);flint_free(ys);
    flint_free(degrees);flint_free(prefix_degree);flint_free(hilbert_dims);
    for(slong d=0;d<=limit;d++){flint_free(spaces[d].keys);flint_free(spaces[d].free);nmod_mat_clear(spaces[d].reduce);}flint_free(spaces);
    for(slong i=0;i<input_count;i++)flint_free(input[i].coeff);flint_free(input);flint_rand_clear(rng);flint_cleanup_master();return 0;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Standalone research backend. mode 1 memoizes requested root coefficients;
 * mode 2 pushes U's contraction one level into the root sum; mode 3 performs
 * the same row-wise contraction after aggregating each root, for comparison.
 * No complete Dixon polynomial or matrix is formed by the experimental run.
 * Usage: n threads prime mode seed depth (depth 1 or 2;
 * n<=7 for unaggregated contraction).
 * The full reference run is separate and is not included in dag_total. */
#define MQ_SHARED_CORE_NO_MAIN
#include "mq_shared_core_test.c"
#include "mq_dag_hook.h"
#include <flint/nmod_vec.h>

typedef struct {uint32_t shift,child;} dag_edge;
typedef struct {
    panel_source *profile;
    coefficient_source coeff;
    const void *children;
    size_t child_count,shifts,count;
    int narrow,mode,active,called,ok;
    slong n,h;
    ulong *factors;
    size_t *offset;
    dag_edge *edges;
    uint32_t *ids;
    unsigned char *ready;
    nmod_poly_mat_struct *core;
    ulong factor;
    nmod_t mod;
    size_t requested,root_products,fused_edges,contract_scalars,edges_count;
    double indexing,solve,multiply;
} dag_state;
static dag_state dag;
static struct {
    int enabled,used,narrow;
    slong n;
    void *leaf,*cache;
    ulong *factors;
    unsigned char *ready;
    size_t previous,count,shifts,*offset,requested,products;
    dag_edge *edges;
    slong parent[FLINT_BITS][FLINT_BITS];
    double indexing;
} checkpoint;

static dag_edge *dag_reverse_edges(size_t **offset,size_t count,const uint32_t *map,
                                   size_t shifts,size_t child_count)
{
    size_t *o=flint_calloc(count+1,sizeof(size_t));
    for(size_t a=0;a<shifts;a++)for(size_t b=0;b<child_count;b++) {
        uint32_t id=map[a*child_count+b];if(id!=UINT32_MAX)o[id+1]++;
    }
    for(size_t i=0;i<count;i++)o[i+1]+=o[i];
    dag_edge *edges=flint_malloc(o[count]*sizeof(dag_edge));
    size_t *cursor=flint_malloc(count*sizeof(size_t));memcpy(cursor,o,count*sizeof(size_t));
    for(size_t a=0;a<shifts;a++)for(size_t b=0;b<child_count;b++) {
        uint32_t id=map[a*child_count+b];if(id!=UINT32_MAX)edges[cursor[id]++]=(dag_edge){a,b};
    }
    flint_free(cursor);*offset=o;return edges;
}

int mq_dag_checkpoint(void *children,int narrow,size_t child_count,
    ulong *factors,size_t shifts,const uint32_t *map,size_t count,slong n,void **output)
{
    if(!dag.active || !checkpoint.enabled)return 0;
    double start=get_wall_time();assert(!checkpoint.used++);
    checkpoint.leaf=children;checkpoint.narrow=narrow;checkpoint.previous=child_count;
    checkpoint.n=n;checkpoint.factors=factors;checkpoint.shifts=shifts;checkpoint.count=count;
    checkpoint.edges=dag_reverse_edges(&checkpoint.offset,count,map,shifts,child_count);
    checkpoint.cache=flint_calloc(n*count,narrow?sizeof(uint32_t):sizeof(ulong));
    checkpoint.ready=flint_calloc(n*count,1);*output=checkpoint.cache;
    ulong choose[FLINT_BITS][FLINT_BITS]={{0}};
    for(slong i=0;i<=n;i++){choose[i][0]=1;for(slong j=1;j<=i;j++)choose[i][j]=choose[i-1][j-1]+choose[i-1][j];}
    for(slong omitted=0;omitted<n;omitted++)for(slong j=0;j<n;j++)if(j!=omitted) {
        slong rank=0,k=0;
        for(slong c=0;c<n;c++)if(c!=omitted && c!=j)rank+=choose[c][++k];
        checkpoint.parent[omitted][j]=rank;
    }
    checkpoint.indexing=get_wall_time()-start;return 1;
}
static ulong dag_child(slong omitted,uint32_t child)
{
    size_t pos=(size_t)(dag.n-1-omitted)*dag.child_count+child;
    if(checkpoint.used && !checkpoint.ready[pos]) {
        __uint128_t total=0;
        for(size_t e=checkpoint.offset[child];e<checkpoint.offset[child+1];e++) {
            dag_edge edge=checkpoint.edges[e];
            for(slong c=0;c<dag.n;c++)if(c!=omitted) {
                ulong f=checkpoint.factors[c*checkpoint.shifts+edge.shift];
                if(!f)continue;
                if((c-(c>omitted))&1)f=nmod_neg(f,dag.mod);
                size_t source=(size_t)checkpoint.parent[omitted][c]*checkpoint.previous+edge.child;
                ulong v=checkpoint.narrow?((uint32_t*)checkpoint.leaf)[source]:((ulong*)checkpoint.leaf)[source];
                total+=(__uint128_t)f*v;checkpoint.products++;
            }
        }
        ulong lo=(ulong)total,hi=(ulong)(total>>64),r;NMOD_RED2(r,hi,lo,dag.mod);
        if(checkpoint.narrow)((uint32_t*)checkpoint.cache)[pos]=r;else ((ulong*)checkpoint.cache)[pos]=r;
        checkpoint.ready[pos]=1;checkpoint.requested++;
    }
    return dag.narrow?((const uint32_t*)dag.children)[pos]:((const ulong*)dag.children)[pos];
}

static __uint128_t dag_edge_sum(dag_edge edge)
{
    __uint128_t total=0;
    const ulong *f=dag.factors+(size_t)edge.shift*dag.n;
    for(slong c=0;c<dag.n;c++)if(f[c]) {
        ulong value=dag_child(c,edge.child);
        total+=(__uint128_t)f[c]*value;
        dag.root_products++;
    }
    return total;
}
static ulong dag_reduce(__uint128_t total)
{
    ulong lo=(ulong)total,hi=(ulong)(total>>64),r;NMOD_RED2(r,hi,lo,dag.mod);return r;
}
static ulong dag_coefficient(const coefficient_source *a,slong r,slong c,slong d)
{
    if(d<0 || d>a->sigma-a->rd[r]-a->cd[c])return 0;
    size_t pos=a->offset[(size_t)r*a->n+c]+d;
    if(!dag.ready[pos]) {
        uint32_t id=dag.ids[pos];__uint128_t sum=0;
        if(id!=UINT32_MAX)for(size_t k=dag.offset[id];k<dag.offset[id+1];k++)
            sum+=dag_edge_sum(dag.edges[k]);
        a->values[pos]=dag_reduce(sum);dag.ready[pos]=1;dag.requested++;
    }
    return a->values[pos];
}
static void dag_contract(nmod_mat_t reduced,const coefficient_source *a,
    const nmod_mat_t x,const slong *cols,slong start,slong k,slong h,slong ubound)
{
    for(slong i=0;i<h;i++)for(slong j=0;j<k;j++) {
        slong col=cols[start+j],bound=FLINT_MIN(ubound,a->sigma-a->rd[i]-a->cd[col]);
        for(slong d=0;d<=bound;d++) {
            ulong *dest=nmod_mat_entry_ptr(reduced,i,d*h);
            const ulong *weights=nmod_mat_entry_ptr(x,j,0);
            if(dag.mode==3) {
                ulong value=dag_coefficient(a,i,col,d);
                if(value){_nmod_vec_scalar_addmul_nmod(dest,weights,x->c,dag.mod.n-value,dag.mod);dag.contract_scalars+=x->c;}
            } else {
                uint32_t id=dag.ids[a->offset[(size_t)i*a->n+col]+d];
                if(id==UINT32_MAX)continue;
                /* Distribute the projection over the last multiplication:
                 * no U coefficient is assembled or cached in this path. */
                for(size_t e=dag.offset[id];e<dag.offset[id+1];e++) {
                    ulong value=dag_reduce(dag_edge_sum(dag.edges[e]));dag.fused_edges++;
                    if(value){_nmod_vec_scalar_addmul_nmod(dest,weights,x->c,dag.mod.n-value,dag.mod);dag.contract_scalars+=x->c;}
                }
            }
        }
    }
}

/* The blocked solver below is shared with the coefficient-panel experiment
 * through a generated specialization kept in source for review. */
#include "mq_dag_solve.h"

int mq_dag_root(const void *children,int narrow,size_t child_count,
    const ulong *factors,size_t shifts,const uint32_t *map,
    const ulong *keys,size_t count,slong words,flint_bitcnt_t bits,
    slong n,const nmod_mpoly_ctx_t ctx)
{
    if(!dag.active)return 0;
    assert(!dag.called++ && FLINT_BITS==64 && ctx->mod.n<((ulong)1<<31));
    double t=get_wall_time();dag.children=children;dag.narrow=narrow;dag.child_count=child_count;
    dag.shifts=shifts;dag.count=count;dag.n=n;dag.mod=ctx->mod;
    dag.factors=flint_malloc(n*shifts*sizeof(ulong));
    for(size_t a=0;a<shifts;a++)for(slong c=0;c<n;c++) {
        ulong value=factors[c*shifts+a];dag.factors[a*n+c]=(c&1)?nmod_neg(value,ctx->mod):value;
    }
    dag.edges=dag_reverse_edges(&dag.offset,count,map,shifts,child_count);dag.edges_count=dag.offset[count];
    panel_source *s=dag.profile;
    coefficient_source *cache=&dag.coeff;
    *cache=(coefficient_source){.n=s->size,.sigma=s->sigma,.rd=s->rd,.cd=s->cd};
    cache->offset=flint_malloc((size_t)s->size*s->size*sizeof(size_t));
    for(slong i=0;i<s->size;i++)for(slong j=0;j<s->size;j++) {
        cache->offset[(size_t)i*s->size+j]=cache->slots;
        cache->slots+=FLINT_MAX(0,s->sigma-s->rd[i]-s->cd[j]+1);
    }
    cache->values=flint_calloc(cache->slots,sizeof(ulong));
    dag.ids=flint_malloc(cache->slots*sizeof(uint32_t));memset(dag.ids,255,cache->slots*sizeof(uint32_t));
    dag.ready=flint_calloc(cache->slots,1);
    monom_t *rm=NULL,*cm=NULL;slong rn=0,cn=0,rc=0,cc=0,rh=16,ch=16;
    hash_entry_t **ri=flint_calloc(rh,sizeof(*ri)),**ci=flint_calloc(ch,sizeof(*ci));
    for(slong i=0;i<s->size;i++){dixon_intern_monom(&rm,&rn,&rc,&ri,&rh,s->rx+i*(n-1),n-1);dixon_intern_monom(&cm,&cn,&cc,&ci,&ch,s->cx+i*(n-1),n-1);}
    ulong exp[2*FLINT_BITS+1];slong ex[2*FLINT_BITS+1];
    for(size_t id=0;id<count;id++) {
        mpoly_get_monomial_ui(exp,keys+id*words,bits,ctx->minfo);
        for(slong v=0;v<2*n-1;v++)ex[v]=exp[v];
        slong r=lookup_monom_index(ri,rh,ex,n-1),c=lookup_monom_index(ci,ch,ex+n-1,n-1),d=ex[2*n-2];
        if(r<0 || c<0)continue;
        assert(d<=s->sigma-s->rd[r]-s->cd[c]);
        size_t pos=cache->offset[(size_t)r*s->size+c]+d;assert(dag.ids[pos]==UINT32_MAX);dag.ids[pos]=id;
    }
    free_monom_index(ri,rh);free_monom_index(ci,ch);flint_free(rm);flint_free(cm);
    dag.indexing=get_wall_time()-t;
    dag.ok=dag_core(dag.core,&dag.factor,cache,dag.h,ctx->mod.n,&dag.multiply,&dag.solve);
    flint_free(dag.factors);flint_free(dag.offset);flint_free(dag.edges);flint_free(dag.ids);flint_free(dag.ready);
    flint_free(cache->values);flint_free(cache->offset);
    if(checkpoint.used) {
        flint_free(checkpoint.leaf);flint_free(checkpoint.factors);flint_free(checkpoint.ready);
        flint_free(checkpoint.offset);flint_free(checkpoint.edges);
    }
    return 1;
}

int main(int argc,char **argv)
{
    slong n=argc>1?atol(argv[1]):7;int threads=argc>2?atoi(argv[2]):1;
    ulong prime=argc>3?strtoul(argv[3],NULL,10):536870923;int mode=argc>4?atoi(argv[4]):1;
    ulong seed=argc>5?strtoul(argv[5],NULL,10):12345;
    int depth=argc>6?atoi(argv[6]):1;assert(depth==1 || depth==2);
    assert(n>=5 && n<=9 && prime<((ulong)1<<31) && mode>=1 && mode<=3);
    assert(mode!=2 || n<=7);
    omp_set_num_threads(threads);flint_set_num_threads(threads);g_dixon_verbose_level=0;g_dixon_det_cache_limit=100000;
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,seed,941);
    fq_nmod_ctx_t fq;fq_nmod_ctx_init_ui(fq,prime,1,"a");
    unified_mpoly_struct *polys=random_mq(n-1,fq,rng),**matrix,**a;
    build_fq_cancellation_matrix(&matrix,polys,n-1,1);perform_fq_matrix_row_operations(&a,&matrix,n-1,1);
    panel_source source={.n=n,.matrix=a};slong h;make_profile(&source,&h);
    nmod_poly_mat_t expected,actual;nmod_poly_mat_init(expected,h,h,prime);nmod_poly_mat_init(actual,h,h,prime);
    double t=get_wall_time();fq_mq_compact compact={0};
    assert(compute_fq_det_mq_compact(&compact,a,n,source.rx,source.size,source.cx,source.size));
    coefficient_source cache;assert(coefficient_source_init(&cache,&source,&compact));
    fq_mq_compact_clear(&compact);double generation=get_wall_time()-t,mul=0,solve=0;
    ulong factor=17;t=get_wall_time();int ok=coefficient_core(expected,&factor,&cache,h,prime,&mul,&solve);
    double reference=get_wall_time()-t;flint_free(cache.values);flint_free(cache.offset);
    dag=(dag_state){.profile=&source,.h=h,.mode=mode,.active=1,.core=actual,.factor=17};
    checkpoint.enabled=depth==2;
    t=get_wall_time();assert(compute_fq_det_mq_compact(&compact,a,n,source.rx,source.size,source.cx,source.size));
    double total=get_wall_time()-t;fq_mq_compact_clear(&compact);
    assert(dag.called==1 && dag.ok==ok);if(ok){assert(dag.factor==factor);assert(nmod_poly_mat_equal(actual,expected));}
    if(checkpoint.enabled) {
        assert(checkpoint.used==1);
        fprintf(stderr,"checkpoint requested=%zu total=%zu products=%zu indexing=%.6f\n",
            checkpoint.requested,n*checkpoint.count,checkpoint.products,checkpoint.indexing);
    }
    fprintf(stderr,"depth=%d root_reverse_index_bytes=%zu checkpoint_fraction=%.9f\n",depth,
        (dag.count+1)*sizeof(size_t)+dag.edges_count*sizeof(dag_edge),
        checkpoint.used?(double)checkpoint.requested/(n*checkpoint.count):0.0);
    printf("{\"n\":%ld,\"prime\":%lu,\"seed\":%lu,\"threads\":%d,\"mode\":%d,\"ok\":%d,\"generation\":%.6f,\"reference_core\":%.6f,\"dag_total\":%.6f,\"indexing\":%.6f,\"solve\":%.6f,\"multiply\":%.6f,\"root_support\":%zu,\"root_edges\":%zu,\"requested_coefficients\":%zu,\"root_products\":%zu,\"fused_edges\":%zu,\"contract_scalars\":%zu,\"equal\":true}\n",
        n,prime,seed,threads,mode,ok,generation,reference,total,dag.indexing,dag.solve,dag.multiply,dag.count,dag.edges_count,dag.requested,dag.root_products,dag.fused_edges,dag.contract_scalars);
    nmod_poly_mat_clear(expected);nmod_poly_mat_clear(actual);
    flint_free(source.rx);flint_free(source.cx);flint_free(source.rd);flint_free(source.cd);
    clear_input(polys,matrix,a,n-1);fq_nmod_ctx_clear(fq);flint_rand_clear(rng);flint_cleanup_master();return 0;
}

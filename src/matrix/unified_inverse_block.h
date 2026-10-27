/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Included after the compact matrix dispatch and conversion helpers. */
#include "inverse_block_det.h"
#include <math.h>
#include <flint/ulong_extras.h>

static void ib_cardinality(fmpz_t q, const field_ctx_t *f)
{
    if(f->field_id==FIELD_ID_FMPZ_MOD) fmpz_set(q,fmpz_mod_ctx_modulus(f->ctx.big_ctx));
    else if(f->field_id==FIELD_ID_NMOD) fmpz_set_ui(q,f->ctx.nmod_ctx.n);
    else if(f->field_id==FIELD_ID_FQ_ZECH) fmpz_set_ui(q,f->ctx.zech_ctx->qm1+1);
    else {
        fmpz_set_ui(q,fq_nmod_ctx_prime(f->ctx.fq_ctx));
        fmpz_pow_ui(q,q,fq_nmod_ctx_degree(f->ctx.fq_ctx));
    }
}

/* Uniform over the entire field, including zero (randtest is not uniform). */
static void ib_random(field_elem_u *x, field_ctx_t *f, flint_rand_t rng)
{
    if(f->field_id==FIELD_ID_FMPZ_MOD) fmpz_randm(&x->big,rng,fmpz_mod_ctx_modulus(f->ctx.big_ctx));
    else if(f->field_id==FIELD_ID_NMOD) x->nmod=n_randint(rng,f->ctx.nmod_ctx.n);
    else if(f->field_id==FIELD_ID_FQ_ZECH) fq_zech_rand(&x->fq_zech,rng,f->ctx.zech_ctx);
    else {
        fq_nmod_t v;fq_nmod_init(v,f->ctx.fq_ctx);fq_nmod_rand(v,rng,f->ctx.fq_ctx);
        fq_nmod_to_field_elem(x,v,f);fq_nmod_clear(v,f->ctx.fq_ctx);
    }
}

static void ib_random_matrix(unified_mat_poly_t a, flint_rand_t rng)
{
    field_ctx_t *f=a->ctx;field_elem_u x;
    field_init_elem(&x,f->field_id,field_ctx_data(f));unified_mat_poly_resize(a,1);
    for(slong i=0;i<a->r*a->c;i++) {ib_random(&x,f,rng);unified_mat_poly_put(a,i,&x);}
    field_clear_elem(&x,f->field_id,field_ctx_data(f));
}

static int ib_points(unified_mat_poly_t points, slong count, flint_rand_t rng, int distinct)
{
    field_ctx_t *f=points->ctx;const void *fc=field_ctx_data(f);
    field_elem_u x,y;field_init_elem(&x,f->field_id,fc);field_init_elem(&y,f->field_id,fc);
    unified_mat_poly_resize(points,count);int ok=1;
    for(slong i=0;i<count;i++) {
        int found=0;
        for(int attempt=0;attempt<64 && !found;attempt++) {
            ib_random(&x,f,rng);found=1;
            if(distinct) for(slong j=0;j<i;j++) {
                unified_mat_poly_get(&y,points,j);
                if(field_equal(&x,&y,f->field_id,fc)) {found=0;break;}
            }
        }
        if(!found) {ok=0;break;}
        unified_mat_poly_put(points,i,&x);
    }
    field_clear_elem(&x,f->field_id,fc);field_clear_elem(&y,f->field_id,fc);return ok;
}

static int ib_reconstruct(unified_poly_t candidate, const unified_mat_poly_t a,
    slong m, slong count, int workers, flint_rand_t rng, dr_inverse_block_stats *stats)
{
    field_ctx_t *f=a->ctx;slong n=a->r;
    unified_mat_poly_t points,u,v,basis,den;
    unified_mat_poly_init(points,1,1,f);unified_mat_poly_init(u,m,n,f);
    unified_mat_poly_init(v,n,m,f);unified_mat_poly_init(basis,2*m,2*m,f);
    unified_mat_poly_init(den,m,m,f);
    int ok=ib_points(points,count,rng,1);
    if(!ok) goto done;
    ib_random_matrix(u,rng);ib_random_matrix(v,rng);
    unified_mat_poly_struct *samples=flint_malloc(count*sizeof(*samples));
    int failed=0;double start=get_wall_time();
#ifdef _OPENMP
#pragma omp parallel for num_threads(workers) schedule(dynamic) reduction(|:failed) if(workers>1)
#endif
    for(slong k=0;k<count;k++) {
        int previous=flint_set_num_workers(0),success;
        unified_mat_poly_t det;unified_mat_poly_init(det,1,1,f);
        unified_mat_poly_init(samples+k,2*m,m,f);
        UM_DISPATCH(a,ib_sample,samples+k,det,a,u,v,points,k,&success);
        failed|=!success;unified_mat_poly_clear(det);flint_reset_num_workers(previous);
    }
    stats->inverse+=get_wall_time()-start;ok=0;
    if(!failed) {
        slong *shift=flint_calloc(2*m,sizeof(slong)),*rows=flint_malloc(2*m*sizeof(slong));
        start=get_wall_time();
        UM_DISPATCH(a,ib_pmintbasis,basis,shift,points,0,samples,count);
        for(slong i=0;i<2*m;i++) rows[i]=i;
        for(slong i=1;i<2*m;i++) {
            slong row=rows[i],j=i;
            while(j && shift[rows[j-1]]>shift[row]) {rows[j]=rows[j-1];j--;}
            rows[j]=row;
        }
        UM_DISPATCH(a,select,den,basis,rows,NULL,0);
        stats->interpolation+=get_wall_time()-start;start=get_wall_time();
        ok=unified_mat_poly_det_hnf(candidate,den);
        stats->determinant+=get_wall_time()-start;flint_free(shift);flint_free(rows);
    }
    for(slong k=0;k<count;k++) unified_mat_poly_clear(samples+k);
    flint_free(samples);
done:
    unified_mat_poly_clear(points);unified_mat_poly_clear(u);unified_mat_poly_clear(v);
    unified_mat_poly_clear(basis);unified_mat_poly_clear(den);return ok;
}

static int ib_verify(unified_poly_t candidate, const unified_mat_poly_t a,
    slong bound, slong checks, int workers, flint_rand_t rng, dr_inverse_block_stats *stats)
{
    if(!candidate->length || candidate->length-1>bound) return 0;
    field_ctx_t *f=a->ctx;const void *fc=field_ctx_data(f);
    unified_mat_poly_t points,u,v,scalar;
    unified_mat_poly_init(points,1,1,f);unified_mat_poly_init(u,0,a->r,f);
    unified_mat_poly_init(v,a->r,0,f);unified_mat_poly_init(scalar,1,1,f);
    unified_mat_poly_resize(scalar,candidate->length);
    for(slong i=0;i<candidate->length;i++) unified_mat_poly_put(scalar,i,candidate->coeffs+i);
    field_elem_u expected,actual,scale;
    field_init_elem(&expected,f->field_id,fc);field_init_elem(&actual,f->field_id,fc);
    field_init_elem(&scale,f->field_id,fc);
    int scaled=0;
    for(int attempt=0;attempt<32 && !scaled;attempt++) {
        ib_points(points,1,rng,0);
        unified_mat_poly_t val,det,block;
        unified_mat_poly_init(val,1,1,f);unified_mat_poly_init(det,1,1,f);unified_mat_poly_init(block,0,0,f);
        UM_DISPATCH(a,ib_evaluate,val,scalar,points,0);unified_mat_poly_get(&actual,val,0);
        if(!field_is_zero(&actual,f->field_id,fc)) {
            int success;UM_DISPATCH(a,ib_sample,block,det,a,u,v,points,0,&success);
            if(success) {
                unified_mat_poly_get(&expected,det,0);field_inv(&scale,&actual,f->field_id,fc);
                field_mul(&scale,&scale,&expected,f->field_id,fc);
                unified_poly_scalar_mul(candidate,candidate,&scale);scaled=1;
            }
        }
        unified_mat_poly_clear(val);unified_mat_poly_clear(det);unified_mat_poly_clear(block);
    }
    int failed=!scaled;
    if(scaled) {
        for(slong i=0;i<candidate->length;i++) unified_mat_poly_put(scalar,i,candidate->coeffs+i);
        ib_points(points,checks,rng,0);
#ifdef _OPENMP
#pragma omp parallel for num_threads(workers) schedule(dynamic) reduction(|:failed) if(workers>1)
#endif
        for(slong k=0;k<checks;k++) {
            int previous=flint_set_num_workers(0),success;
            unified_mat_poly_t val,det,block;
            unified_mat_poly_init(val,1,1,f);unified_mat_poly_init(det,1,1,f);unified_mat_poly_init(block,0,0,f);
            UM_DISPATCH(a,ib_evaluate,val,scalar,points,k);
            UM_DISPATCH(a,ib_sample,block,det,a,u,v,points,k,&success);
            field_elem_u x,y;field_init_elem(&x,f->field_id,fc);field_init_elem(&y,f->field_id,fc);
            unified_mat_poly_get(&x,val,0);unified_mat_poly_get(&y,det,0);
            failed|=!field_equal(&x,&y,f->field_id,fc);
            field_clear_elem(&x,f->field_id,fc);field_clear_elem(&y,f->field_id,fc);
            unified_mat_poly_clear(val);unified_mat_poly_clear(det);unified_mat_poly_clear(block);
            flint_reset_num_workers(previous);
        }
        stats->checks+=checks;
    }
    field_clear_elem(&expected,f->field_id,fc);field_clear_elem(&actual,f->field_id,fc);field_clear_elem(&scale,f->field_id,fc);
    unified_mat_poly_clear(points);unified_mat_poly_clear(u);unified_mat_poly_clear(v);unified_mat_poly_clear(scalar);
    return !failed;
}

int dr_unified_mat_poly_det_inverse_block_try(unified_poly_t out,
    const unified_mat_poly_t a, slong block, int threads,
    flint_rand_t rng, dr_inverse_block_stats *stats)
{
    dr_inverse_block_stats empty={0};if(!stats) stats=&empty;*stats=empty;
    slong n=a->r;
    if(n<1 || n!=a->c || out->ctx!=a->ctx || block<0 || block>n ||
        n>WORD_MAX/n || (size_t)n>SIZE_MAX/(size_t)n/sizeof(slong)) return 0;
    int workers=1;
#ifdef _OPENMP
    if(!omp_in_parallel()) workers=threads>0 ? threads : omp_get_max_threads();
#else
    (void)threads;
#endif
    double start=get_wall_time();
    slong *degrees=flint_malloc((size_t)n*n*sizeof(slong));
    UM_DISPATCH(a,ib_degrees,degrees,a);
    slong bound=dr_inverse_block_degree_bound(degrees,n);flint_free(degrees);
    stats->degree_bound=bound;
    if(bound<1 || bound>WORD_MAX/2) return 0;
    fmpz_t q;fmpz_init(q);ib_cardinality(q,a->ctx);
    int enough=fmpz_cmp_ui(q,2*(ulong)bound)>0;
    /* A lower bound for log2(q) makes the verification count conservative. */
    double logq=(double)(fmpz_bits(q)-1);
    if(fmpz_bits(q)<=FLINT_BITS) logq=log2((double)fmpz_get_ui(q));
    if(!enough) {fmpz_clear(q);return 0;}
    slong checks=FLINT_MAX(8,(slong)ceil(66.0/(logq-log2((double)bound))));
    slong m=block ? block : dr_inverse_block_choose(n,bound,workers);
    if(m<1 || m>WORD_MAX/2 || (size_t)m>SIZE_MAX/(size_t)m/4/sizeof(field_elem_u)) {
        fmpz_clear(q);return 0;
    }
    slong half=bound/m+(bound%m!=0),count=2*half+1;
    stats->block=m;stats->setup=get_wall_time()-start;
    unified_poly_t candidate;unified_poly_init(candidate,a->ctx);int ok=0;
    for(int attempt=0;attempt<2;attempt++) {
        if(fmpz_cmp_ui(q,(ulong)count)<0 || (size_t)count>SIZE_MAX/sizeof(unified_mat_poly_struct)) break;
        stats->attempts++;stats->points=count;unified_poly_zero(candidate);
        if(ib_reconstruct(candidate,a,m,count,workers,rng,stats)) {
            start=get_wall_time();ok=ib_verify(candidate,a,bound,checks,workers,rng,stats);
            stats->verification+=get_wall_time()-start;
            if(ok) {unified_poly_set(out,candidate);break;}
        }
        if(count>(WORD_MAX-1)/2) break;
        count=2*count+1;
    }
    unified_poly_clear(candidate);fmpz_clear(q);return ok;
}

int dr_unified_mat_poly_det_inverse_block_logged(unified_poly_t out,
    const unified_mat_poly_t matrix)
{
    extern slong g_dixon_inverse_block_size;
    extern void dixon_info_log(const char *format, ...);
    flint_rand_t rng;flint_rand_init(rng);dr_inverse_block_stats stats;double start=get_wall_time();
    int ok=dr_unified_mat_poly_det_inverse_block_try(out,matrix,g_dixon_inverse_block_size,0,rng,&stats);
    flint_rand_clear(rng);
    if(ok && g_dixon_verbose_level>=2) {
        dixon_info_log("  Dixon inverse-block: %ld x %ld, degree bound %ld, block %ld, points %ld, attempts %ld\n",
            matrix->r,matrix->c,stats.degree_bound,stats.block,stats.points,stats.attempts);
        dixon_info_log("  Dixon inverse-block time: %.6fs (setup %.6f, sampling %.6f, reconstruction %.6f, determinant %.6f, verification %.6f; %ld checks)\n",
            get_wall_time()-start,stats.setup,stats.inverse,stats.interpolation,stats.determinant,stats.verification,stats.checks);
    } else if(!ok) dixon_info_log("  Dixon inverse-block: insufficient field points, unavailable or verification failed; fallback to existing determinant backend\n");
    return ok;
}

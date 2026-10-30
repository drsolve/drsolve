/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "bezout_block.h"

static slong initial_points_setting = -1;
void dr_bezout_block_set_initial_points(slong points)
{
    initial_points_setting = points < 0 ? -1 : points;
}
#ifdef HAVE_PML
#include "bezout_prefix.h"
#include "pml_det.h"
#include <stdint.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include <flint/fmpz.h>
#include <flint/gr.h>
#include <flint/gr_poly.h>
#include <flint/fq_nmod.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static double seconds(void)
{
    struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);
    return t.tv_sec+1e-9*t.tv_nsec;
}
static int eval(gr_ptr out,const nmod_poly_struct *p,gr_srcptr point,gr_ctx_t ctx)
{
    int status=gr_zero(out,ctx);
    for(slong i=p->length;i-->0;) {
        status|=gr_mul(out,out,point,ctx);
        status|=gr_add_ui(out,out,p->coeffs[i],ctx);
    }
    return status;
}
/* Recover the base-field scalar at one point; then check eight independent
 * points. Use an extension when base-field checks have a weak degree bound.
 * 80 bits per candidate leave room for the adaptive-attempt union bound. */
static int verify(nmod_poly_t candidate,const nmod_poly_struct *coeffs,
    slong d,slong bound,flint_rand_t state,slong *extension_degree)
{
    if(!candidate->length || nmod_poly_degree(candidate)>bound)return 0;
    ulong p=candidate->mod.n;slong ext=1;
    fmpz_t q,target;fmpz_init_set_ui(q,p);fmpz_init_set_ui(target,(ulong)bound+1);
    fmpz_mul_2exp(target,target,10);
    while(fmpz_cmp(q,target)<0){fmpz_mul_ui(q,q,p);ext++;}
    fmpz_clear(q);fmpz_clear(target);*extension_degree=ext;
    fq_nmod_ctx_t fq;gr_ctx_t ctx;
    int status=GR_SUCCESS;
    if(ext==1){status=gr_ctx_init_nmod(ctx,p);status|=gr_ctx_set_is_field(ctx,T_TRUE);}
    else {fq_nmod_ctx_init_ui(fq,p,ext,"v");_gr_ctx_init_fq_nmod_from_ref(ctx,fq);}
    gr_ptr point=gr_heap_init(ctx),value=gr_heap_init(ctx),expected=gr_heap_init(ctx);
    gr_ptr actual=gr_heap_init(ctx),scale=gr_heap_init(ctx),product=gr_heap_init(ctx);
    gr_poly_t f,g;gr_poly_init(f,ctx);gr_poly_init(g,ctx);
    int have_scale=0,checks=0,ok=0;ulong scalar=0;
    for(int tries=0;tries<32 && checks<8 && status==GR_SUCCESS;tries++) {
        if(ext==1)*(ulong *)point=n_randint(state,p);
        else fq_nmod_rand(point,state,fq);
        for(int k=0;k<2;k++)for(slong i=0;i<=d;i++) {
            status|=eval(value,coeffs+k*(d+1)+i,point,ctx);
            status|=gr_poly_set_coeff_scalar(k?g:f,i,value,ctx);
        }
        status|=gr_poly_resultant(expected,f,g,ctx);
        status|=eval(actual,candidate,point,ctx);
        if(status!=GR_SUCCESS)break;
        if(!have_scale) {
            if(gr_is_zero(expected,ctx)==T_TRUE)continue;
            if(gr_is_zero(actual,ctx)==T_TRUE)break;
            status|=gr_div(scale,expected,actual,ctx);
            if(ext==1)scalar=*(ulong *)scale;
            else if(((nmod_poly_struct *)scale)->length<=1)
                scalar=nmod_poly_get_coeff_ui((nmod_poly_struct *)scale,0);
            else status|=GR_DOMAIN;
            if(status!=GR_SUCCESS || !scalar)break;
            have_scale=1;continue;
        }
        status|=gr_mul(product,actual,scale,ctx);
        if(status!=GR_SUCCESS || gr_equal(product,expected,ctx)!=T_TRUE)break;
        checks++;
    }
    if(status==GR_SUCCESS && checks==8){nmod_poly_scalar_mul_nmod(candidate,candidate,scalar);ok=1;}
    gr_poly_clear(f,ctx);gr_poly_clear(g,ctx);
    gr_heap_clear(point,ctx);gr_heap_clear(value,ctx);gr_heap_clear(expected,ctx);
    gr_heap_clear(actual,ctx);gr_heap_clear(scale,ctx);gr_heap_clear(product,ctx);
    if(ext==1)gr_ctx_clear(ctx);else fq_nmod_ctx_clear(fq);
    return ok;
}
static int reconstruct(nmod_poly_t result,const nmod_poly_struct *coeffs,
    slong d,slong m,slong count,int shifted,ulong root,ulong dilation,
    const slong *indices,dr_bezout_block_stats *stats)
{
    const slong rows=2*(d+1);const ulong p=result->mod.n;
    if((size_t)count>SIZE_MAX/sizeof(ulong)/(size_t)rows ||
       (size_t)count>SIZE_MAX/sizeof(nmod_mat_struct))return 0;
    double start=seconds();
    ulong *values=flint_malloc((size_t)rows*count*sizeof(ulong));
    nmod_geometric_progression_t geo;
    _nmod_geometric_progression_init_function(geo,root,count,result->mod,1);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(!omp_in_parallel() && rows>=8)
#endif
    for(slong i=0;i<rows;i++) {
        nmod_poly_t scaled;nmod_poly_init(scaled,p);
        nmod_poly_fit_length(scaled,coeffs[i].length);
        ulong power=1;
        for(slong j=0;j<coeffs[i].length;j++) {
            scaled->coeffs[j]=nmod_mul(coeffs[i].coeffs[j],power,result->mod);
            power=nmod_mul(power,dilation,result->mod);
        }
        scaled->length=coeffs[i].length;
        _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(values+(size_t)i*count,
            scaled->coeffs,scaled->length,geo,count,result->mod);
        nmod_poly_clear(scaled);
    }
    nmod_geometric_progression_clear(geo);
    stats->evaluation+=seconds()-start;start=seconds();
    nmod_mat_struct *blocks=flint_malloc(count*sizeof(*blocks));
    int failed=0;
#ifdef _OPENMP
    #pragma omp parallel reduction(|:failed) if(!omp_in_parallel() && count>=32)
#endif
    {
        nmod_poly_t f,g,h;nmod_poly_init(f,p);nmod_poly_init(g,p);nmod_poly_init(h,p);
        nmod_poly_fit_length(f,d+1);nmod_poly_fit_length(g,d+1);
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for(slong k=0;k<count;k++) {
            nmod_mat_init(blocks+k,2*m,m,p);
            for(slong i=0;i<=d;i++) {
                f->coeffs[i]=values[(size_t)i*count+k];
                g->coeffs[i]=values[(size_t)(d+1+i)*count+k];
            }
            f->length=g->length=d+1;
            /* The x^d coefficient of f is a nonzero constant. Cancel the
             * constant x^d coefficient of g (already zero after reduction).
             * Validation uses the pair before this cancellation. */
            ulong scale=nmod_div(g->coeffs[d],f->coeffs[d],result->mod);
            for(slong i=0;i<=d;i++)g->coeffs[i]=nmod_sub(g->coeffs[i],nmod_mul(scale,f->coeffs[i],result->mod),result->mod);
            _nmod_poly_normalise(g);
            slong precision=2*(indices?indices[m-1]:m-1)+1;
            if(dr_bezout_left_moments(h,f,g,precision,1024,NULL)!=1){failed=1;continue;}
            for(slong i=0;i<m;i++) {
                for(slong j=0;j<m;j++)nmod_mat_entry(blocks+k,i,j)=
                    nmod_poly_get_coeff_ui(h,indices?indices[i]+indices[j]:i+j);
                nmod_mat_entry(blocks+k,m+i,i)=1;
            }
        }
        nmod_poly_clear(f);nmod_poly_clear(g);nmod_poly_clear(h);
    }
    flint_free(values);stats->inverse+=seconds()-start;
    int ok=0;
    if(!failed) {
        slong *shift=flint_calloc(2*m,sizeof(slong)),*weights=flint_calloc(2*m,sizeof(slong));
        slong *order=flint_malloc(2*m*sizeof(slong));
        if(shifted)for(slong j=0;j<m;j++){weights[j]=j;weights[m+j]=2*d-1-j;}
        for(slong i=0;i<2*m;i++){shift[i]=weights[i];order[i]=i;}
        nmod_poly_mat_t basis,den;nmod_poly_mat_init(basis,2*m,2*m,p);nmod_poly_mat_init(den,m,m,p);
        start=seconds();
        nmod_poly_mat_pmintbasis_geometric(basis,shift,NULL,blocks,root,count);
        /* Ordered weak Popov is pivot ordered. Select the m lowest shifted
         * row degrees, retaining row order on ties as in the research probe. */
        for(slong i=1;i<2*m;i++) {
            slong row=order[i],j=i;
            while(j>0 && shift[order[j-1]]>shift[row]){order[j]=order[j-1];j--;}
            order[j]=row;
        }
        for(slong i=0;i<m;i++)for(slong j=0;j<m;j++)
            nmod_poly_set(nmod_poly_mat_entry(den,i,j),nmod_poly_mat_entry(basis,order[i],j));
        stats->interpolation+=seconds()-start;start=seconds();
        ok=nmod_poly_mat_det_hnf(result,den);
        stats->determinant+=seconds()-start;
        nmod_poly_mat_clear(basis);nmod_poly_mat_clear(den);
        flint_free(shift);flint_free(weights);flint_free(order);
        if(ok) {
            ulong inv=nmod_inv(dilation,result->mod),power=1;
            for(slong i=0;i<result->length;i++) {
                result->coeffs[i]=nmod_mul(result->coeffs[i],power,result->mod);
                power=nmod_mul(power,inv,result->mod);
            }
        }
    }
    for(slong k=0;k<count;k++)nmod_mat_clear(blocks+k);
    flint_free(blocks);return ok;
}

/* Integer cube root: the direct matrix FFT uses classical matrix products.
 * This is an asymptotic block-size choice, not a cap on the input size. */
static slong default_block(slong d)
{
    slong lo=1,hi=d;
    while(lo<hi) {
        slong mid=lo+(hi-lo+1)/2;
        if(mid<=d/mid/mid)lo=mid;else hi=mid-1;
    }
    return lo;
}

/* Encode x=z^(A+1), y=z. Since Q has constant nonzero x^B coefficient,
 * division in K[y][x] preserves total degree <= A. Its remainder has x-degree
 * < B, so its encoding has degree < B*(A+1). Thus ordinary univariate
 * remainder is exactly the encoding of P mod Q (no coefficient carry).
 * Normalize Q to monic so specialization of Res(Q, remainder) remains valid
 * even when the degree of the remainder drops at a verification point. */
static int prepare_unbalanced(nmod_poly_struct *coeffs,ulong *correction,
    const nmod_mpoly_t a,const nmod_mpoly_t b,slong var,
    const nmod_mpoly_ctx_t ctx,slong A,slong B,int a_is_large)
{
    const nmod_mpoly_struct *P=a_is_large?a:b,*Q=a_is_large?b:a;
    slong base=A+1;
    if(A>(WORD_MAX-1)/base)return 0;
    ulong p=nmod_mpoly_ctx_modulus(ctx),exps[2];
    nmod_poly_t packed_p,packed_q,rem;
    nmod_poly_init(packed_p,p);nmod_poly_init(packed_q,p);nmod_poly_init(rem,p);
    const nmod_mpoly_struct *inputs[2]={P,Q};
    for(int k=0;k<2;k++)for(slong i=0;i<inputs[k]->length;i++) {
        nmod_mpoly_get_term_exp_ui(exps,inputs[k],i,ctx);
        slong index=(slong)exps[var]*base+(slong)exps[1-var];
        nmod_poly_set_coeff_ui(k?packed_q:packed_p,index,
            nmod_mpoly_get_term_coeff_ui(inputs[k],i,ctx));
    }
    ulong leading=nmod_poly_get_coeff_ui(packed_q,B*base);
    nmod_poly_scalar_mul_nmod(packed_q,packed_q,nmod_inv(leading,packed_q->mod));
    nmod_poly_rem(rem,packed_p,packed_q);
    int ok=!nmod_poly_is_zero(rem);
    if(ok) {
        const nmod_poly_struct *reduced[2]={packed_q,rem};
        for(int k=0;k<2;k++)for(slong i=0;i<reduced[k]->length;i++)
            if(reduced[k]->coeffs[i])
                nmod_poly_set_coeff_ui(coeffs+k*(B+1)+i/base,i%base,reduced[k]->coeffs[i]);
        slong degree_p=nmod_mpoly_degree_si(P,var,ctx);
        *correction=nmod_pow_ui(leading,(ulong)degree_p,packed_q->mod);
        /* Res(P,Q)=(-1)^(deg_x(P)*B) lc_x(Q)^deg_x(P) Res(Q_monic,P mod Q). */
        if(a_is_large && (degree_p&1) && (B&1))
            *correction=nmod_neg(*correction,packed_q->mod);
    }
    nmod_poly_clear(packed_p);nmod_poly_clear(packed_q);nmod_poly_clear(rem);
    return ok;
}
int dr_nmod_mpoly_resultant_bezout_block_try(nmod_poly_t out,
    const nmod_mpoly_t a,const nmod_mpoly_t b,slong var,
    const nmod_mpoly_ctx_t ctx,slong block,flint_rand_t state,dr_bezout_block_stats *stats)
{
    dr_bezout_block_stats local={0};if(!stats)stats=&local;*stats=local;
    if(nmod_mpoly_ctx_nvars(ctx)!=2 || var<0 || var>1 || block<0 ||
       out->mod.n!=nmod_mpoly_ctx_modulus(ctx) || !n_is_prime(out->mod.n) ||
       !nmod_mpoly_degrees_fit_si(a,ctx) || !nmod_mpoly_degrees_fit_si(b,ctx) ||
       !nmod_mpoly_total_degree_fits_si(a,ctx) || !nmod_mpoly_total_degree_fits_si(b,ctx))return 0;
    slong da=nmod_mpoly_total_degree_si(a,ctx),db=nmod_mpoly_total_degree_si(b,ctx);
    slong A=FLINT_MAX(da,db),d=FLINT_MIN(da,db);
    int unbalanced=da!=db;
    if(d<1 || A>(WORD_MAX-2)/2 || A>(WORD_MAX-1)/d)return 0;
    if(unbalanced) {
        const nmod_mpoly_struct *small=da<db?a:b;
        if(nmod_mpoly_degree_si(small,var,ctx)!=d)return 0;
    } else if(nmod_mpoly_degree_si(a,var,ctx)!=d ||
              nmod_mpoly_degree_si(b,var,ctx)!=d)return 0;
    slong m=block?block:default_block(d);
    if(!block && !unbalanced)m=FLINT_MAX(m,FLINT_MIN(6,d));
    slong bound=A*d;
    if(m>d || m>WORD_MAX/2 ||
       (size_t)(2*m)>SIZE_MAX/sizeof(nmod_poly_struct)/(size_t)(2*m))return 0;
    slong short_count,full_count;
    if(unbalanced) {
        slong stride=d/m+(d%m!=0);
        if(A>((WORD_MAX-1)/2)/stride)return 0;
        short_count=full_count=2*A*stride+1;
    } else {
        slong half=bound/m+(bound%m!=0);
        if(half>(WORD_MAX-1-2*d)/2)return 0;
        short_count=1+2*half-2*d+2*m;full_count=1+2*half+2*d;
    }
    /* The total-degree bound is often very loose for composed systems.
     * Start with O(d) points, then double until independent resultant checks
     * succeed. Only the full bound is trusted by verification; an undersampled
     * reconstruction is never accepted just because its degree looks small. */
    slong count=unbalanced?short_count:FLINT_MIN(short_count,2*d+1);
    /* Optional starting budget; correctness still requires the same checks.
     * Full_count is the existing theoretical reconstruction budget. */
    const char *initial=getenv("DRSOLVE_BEZOUT_INITIAL_POINTS");
    if(!unbalanced && initial_points_setting < 0 && initial && *initial) {
        char *end;errno=0;
        unsigned long value=strtoul(initial,&end,10);
        if(errno || initial[0]=='-' || end==initial || *end || !value || value>(ulong)WORD_MAX) {
            fprintf(stderr,"Invalid DRSOLVE_BEZOUT_INITIAL_POINTS: expected a positive slong; using automatic budget.\n");
        } else count=FLINT_MIN((slong)value,full_count);
    }
    if(!unbalanced && initial_points_setting > 0)
        count=FLINT_MIN(initial_points_setting,full_count);
    const char *profile_env=getenv("DRSOLVE_PROFILE_BEZOUT");
    int profile=profile_env && *profile_env && strcmp(profile_env,"0");
    if((ulong)count>(out->mod.n-1)/2)return 0;
    slong rows=2*(d+1);
    if((size_t)rows>SIZE_MAX/sizeof(nmod_poly_struct))return 0;
    nmod_poly_struct *coeffs=flint_malloc(rows*sizeof(*coeffs));
    for(slong i=0;i<rows;i++)nmod_poly_init(coeffs+i,out->mod.n);
    ulong correction=1;int prepared=1;
    double prep_start=seconds();
    if(unbalanced)prepared=prepare_unbalanced(coeffs,&correction,a,b,var,ctx,A,d,da>db);
    else {
        const nmod_mpoly_struct *inputs[2]={a,b};ulong exps[2];
        for(int k=0;k<2;k++)for(slong i=0;i<inputs[k]->length;i++) {
            nmod_mpoly_get_term_exp_ui(exps,inputs[k],i,ctx);
            nmod_poly_set_coeff_ui(coeffs+k*(d+1)+exps[var],exps[1-var],
                nmod_mpoly_get_term_coeff_ui(inputs[k],i,ctx));
        }
    }
    if(unbalanced)stats->evaluation+=seconds()-prep_start;
    slong *indices=NULL;
    if(unbalanced) {
        indices=flint_malloc(m*sizeof(*indices));
        slong stride=d/m+(d%m!=0);
        /* The default m=floor(cuberoot(d)) permits a uniform stride. For
         * larger caller-selected blocks, rounded spacing still covers all
         * d basis positions with gaps <= ceil(d/m), with the same degree
         * bound. Neither branch requires a cap on valid block sizes. */
        for(slong i=0;i<m;i++)indices[i]=(m-1)<=(d-1)/stride?i*stride:i*d/m;
    }
    nmod_poly_t candidate;nmod_poly_init(candidate,out->mod.n);
    ulong root=prepared?n_primitive_root_prime(out->mod.n):1;int ok=0;
    stats->block=m;
    for(int attempt=0;prepared;attempt++) {
        if((ulong)count>(out->mod.n-1)/2)break;
        stats->points=count;stats->attempts++;
        ulong dilation=1+n_randint(state,out->mod.n-1);
        nmod_poly_zero(candidate);
        double attempt_start=seconds();
        dr_bezout_block_stats before=*stats;
        int reconstructed=reconstruct(candidate,coeffs,d,m,count,
                       !unbalanced && !attempt && count==short_count,
                       root,dilation,indices,stats);
        if(reconstructed) {
            double start=seconds();
            ok=verify(candidate,coeffs,d,bound,state,&stats->verification_extension_degree);
            stats->verification+=seconds()-start;
        }
        if(profile)fprintf(stderr,
            "Bezout attempt %d: points=%ld status=%s degree=%ld wall=%.6f evaluation=%.6f inverse=%.6f interpolation=%.6f determinant=%.6f verification=%.6f\n",
            attempt+1,count,ok?"accepted":reconstructed?"verification_failed":"reconstruction_failed",
            nmod_poly_degree(candidate),seconds()-attempt_start,
            stats->evaluation-before.evaluation,stats->inverse-before.inverse,
            stats->interpolation-before.interpolation,stats->determinant-before.determinant,
            stats->verification-before.verification);
        if(ok){
            nmod_poly_scalar_mul_nmod(candidate,candidate,correction);
            nmod_poly_swap(out,candidate);break;
        }
        if(unbalanced) {
            if(attempt)break;
            count=full_count;
        } else {
            if(count==full_count)break;
            count=count>(full_count-1)/2?full_count:2*count+1;
        }
    }
    nmod_poly_clear(candidate);flint_free(indices);
    for(slong i=0;i<rows;i++)nmod_poly_clear(coeffs+i);
    flint_free(coeffs);return ok;
}

#else
int dr_nmod_mpoly_resultant_bezout_block_try(nmod_poly_t out,
    const nmod_mpoly_t a,const nmod_mpoly_t b,slong var,
    const nmod_mpoly_ctx_t ctx,slong block,flint_rand_t state,dr_bezout_block_stats *stats)
{
    (void)out;(void)a;(void)b;(void)var;(void)ctx;(void)block;(void)state;
    if(stats){dr_bezout_block_stats empty={0};*stats=empty;}
    return 0;
}
#endif

#include "unified_mpoly_resultant.h"
#include "dixon_flint.h"

int unified_mpoly_resultant_bezout_block_try(unified_mpoly_t out,
    const unified_mpoly_t a,const unified_mpoly_t b,slong var,unified_mpoly_ctx_t ctx)
{
    if(ctx->nvars!=2 || ctx->field_ctx->field_id!=FIELD_ID_NMOD ||
       a->ctx_ptr!=ctx || b->ctx_ptr!=ctx || out->ctx_ptr!=ctx ||
       var<0 || var>1 || g_field_equation_reduction)return 0;
    slong da=0,db=0;
    nmod_poly_t result;nmod_poly_init(result,ctx->field_ctx->ctx.nmod_ctx.n);
    flint_rand_t state;flint_rand_init(state);
    dr_bezout_block_stats stats;
    double start=get_wall_time();
    int ok=dr_nmod_mpoly_resultant_bezout_block_try(result,GET_NMOD_POLY(a),GET_NMOD_POLY(b),
        var,GET_NMOD_CTX(ctx),0,state,&stats);
    if(ok) {
        da=nmod_mpoly_total_degree_si(GET_NMOD_POLY(a),GET_NMOD_CTX(ctx));
        db=nmod_mpoly_total_degree_si(GET_NMOD_POLY(b),GET_NMOD_CTX(ctx));
        unified_mpoly_t answer=unified_mpoly_init(ctx);ulong exp[2]={0,0};
        for(slong i=result->length;i-->0;) {
            exp[1-var]=i;
            nmod_mpoly_set_coeff_ui_ui(GET_NMOD_POLY(answer),result->coeffs[i],exp,GET_NMOD_CTX(ctx));
        }
        unified_mpoly_swap(out,answer);unified_mpoly_clear(answer);
        dixon_info_log("Dixon/Bezout block-Hankel resultant: %s, block %ld, %ld points, %ld attempt(s)\n",
            da==db?"C adaptive interpolation":"C reduced spaced projection",
            stats.block,stats.points,stats.attempts);
        if(da!=db)dixon_info_log("  Total degrees %ld/%ld; reduced Bezout order %ld (preparation included in evaluation)\n",
            da,db,FLINT_MIN(da,db));
        dixon_info_log("  Evaluation %.6fs | inverse %.6fs | interpolation %.6fs | determinant %.6fs | verification %.6fs (field degree %ld)\n",
            stats.evaluation,stats.inverse,stats.interpolation,stats.determinant,stats.verification,stats.verification_extension_degree);
        dixon_info_log("Dixon/Bezout block time: Wall time: %.6f seconds\n",
            get_wall_time()-start);
    }
    flint_rand_clear(state);nmod_poly_clear(result);return ok;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "inverse_block_det.h"
#include "dixon_flint.h"
#include <assert.h>

static void check(field_ctx_t *f, const fq_nmod_ctx_struct *fq, flint_rand_t rng)
{
    const void *fc=field_ctx_data(f);
    unified_poly_mat_t a;unified_poly_mat_init(a,4,4,f);
    unified_mat_poly_t mat;unified_mat_poly_init(mat,4,4,f);
    unified_poly_t ref,got,sentinel;
    unified_poly_init(ref,f);unified_poly_init(got,f);unified_poly_init(sentinel,f);
    field_elem_u c;field_init_elem(&c,f->field_id,fc);
    field_set_one(&c,f->field_id,fc);unified_poly_set_coeff(sentinel,0,&c);
    fq_nmod_t v;if(fq) fq_nmod_init(v,fq);
    int successes=0;
    for(int trial=0;trial<4;trial++) {
        for(slong i=0;i<4;i++) for(slong j=0;j<4;j++) {
            unified_poly_zero(a->rows[i]+j);
            for(slong d=0;d<4;d++) {
                if(fq) {fq_nmod_rand(v,rng,fq);fq_nmod_to_field_elem(&c,v,f);}
                else fmpz_randm(&c.big,rng,fmpz_mod_ctx_modulus(f->ctx.big_ctx));
                unified_poly_set_coeff(a->rows[i]+j,d,&c);
            }
        }
        unified_mat_poly_from_unified(mat,a);assert(unified_mat_poly_det_hnf(ref,mat));
        for(slong block=1;block<=4;block++) {
            unified_poly_set(got,sentinel);dr_inverse_block_stats stats;
            int ok=dr_unified_mat_poly_det_inverse_block_try(got,mat,block,trial&1 ? 2 : 1,rng,&stats);
            if(ok) {
                successes++;assert(got->length==ref->length);
                for(slong k=0;k<got->length;k++) assert(field_equal(got->coeffs+k,ref->coeffs+k,f->field_id,fc));
            } else {
                assert(got->length==1 && field_equal(got->coeffs,sentinel->coeffs,f->field_id,fc));
            }
        }
    }
    assert(successes>8);
    /* Repeated invariant factors cannot be recovered by a scalar projection.
     * A full block does recover multiplicity; a failed small block is inert. */
    field_set_one(&c,f->field_id,fc);
    for(slong i=0;i<4;i++) for(slong j=0;j<4;j++) {
        unified_poly_zero(a->rows[i]+j);
        if(i==j) {
            unified_poly_set_coeff(a->rows[i]+j,0,&c);
            unified_poly_set_coeff(a->rows[i]+j,1,&c);
        }
    }
    unified_mat_poly_from_unified(mat,a);unified_poly_set(got,sentinel);
    assert(!dr_unified_mat_poly_det_inverse_block_try(got,mat,1,1,rng,NULL));
    assert(got->length==1 && field_equal(got->coeffs,sentinel->coeffs,f->field_id,fc));
    assert(unified_mat_poly_det_hnf(ref,mat));
    assert(dr_unified_mat_poly_det_inverse_block_try(got,mat,4,2,rng,NULL));
    assert(got->length==ref->length);
    for(slong k=0;k<got->length;k++) assert(field_equal(got->coeffs+k,ref->coeffs+k,f->field_id,fc));
    /* Identically singular input must never publish a reconstructed answer. */
    for(slong j=0;j<4;j++) unified_poly_set(a->rows[1]+j,a->rows[0]+j);
    unified_mat_poly_from_unified(mat,a);unified_poly_set(got,sentinel);
    assert(!dr_unified_mat_poly_det_inverse_block_try(got,mat,2,1,rng,NULL));
    assert(got->length==1 && field_equal(got->coeffs,sentinel->coeffs,f->field_id,fc));
    assert(!dr_unified_mat_poly_det_inverse_block_try(got,mat,5,1,rng,NULL));
    if(fq) fq_nmod_clear(v,fq);
    field_clear_elem(&c,f->field_id,fc);unified_poly_clear(ref);unified_poly_clear(got);
    unified_poly_clear(sentinel);unified_mat_poly_clear(mat);unified_poly_mat_clear(a);
}
int main(void)
{
    g_dixon_verbose_level=0;g_dixon_debug_mode=0;
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,4321,9123);
    const ulong primes[]={3,101,2,2,2,2,2,2,2};
    const slong degrees[]={6,3,8,16,32,64,128,9,7};
    fmpz_t p;fmpz_init(p);
    for(size_t i=0;i<sizeof(primes)/sizeof(*primes);i++) {
        fmpz_set_ui(p,primes[i]);fq_nmod_ctx_t fq;fq_nmod_ctx_init(fq,p,degrees[i],"t");
        field_ctx_t field;field_ctx_init(&field,fq);check(&field,fq,rng);
        printf("PASS unified inverse-block %lu^%ld (%s)\n",primes[i],degrees[i],field.description);
        field_ctx_clear(&field);fq_nmod_ctx_clear(fq);
    }
    fmpz_set_str(p,"18446744073709551629",10);field_ctx_t big;field_ctx_init_prime(&big,p);
    check(&big,NULL,rng);field_ctx_clear(&big);
    /* Insufficient cardinality: deterministic early fallback, no attempts. */
    fmpz_set_ui(p,2);fq_nmod_ctx_t fq;fq_nmod_ctx_init(fq,p,4,"t");
    field_ctx_t field;field_ctx_init(&field,fq);
    unified_poly_mat_t a;unified_poly_mat_init(a,1,1,&field);
    field_elem_u one;field_init_elem(&one,field.field_id,field_ctx_data(&field));
    field_set_one(&one,field.field_id,field_ctx_data(&field));unified_poly_set_coeff(a->rows[0],12,&one);
    unified_mat_poly_t mat;unified_mat_poly_init(mat,1,1,&field);unified_mat_poly_from_unified(mat,a);
    unified_poly_t got;unified_poly_init(got,&field);unified_poly_set_coeff(got,0,&one);
    dr_inverse_block_stats stats;
    assert(!dr_unified_mat_poly_det_inverse_block_try(got,mat,1,1,rng,&stats));
    assert(!stats.attempts && got->length==1);
    unified_poly_clear(got);unified_mat_poly_clear(mat);unified_poly_mat_clear(a);
    field_clear_elem(&one,field.field_id,field_ctx_data(&field));field_ctx_clear(&field);fq_nmod_ctx_clear(fq);
    fmpz_clear(p);flint_rand_clear(rng);puts("Unified inverse-block exact comparisons and fallback passed");
}

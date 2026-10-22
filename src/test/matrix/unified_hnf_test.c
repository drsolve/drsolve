/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_mat_poly.h"
#include "unified_byte.h"
#include "fq_poly_mat_det.h"
#include <assert.h>
#include <limits.h>
#include <flint/fmpz.h>

static void expansion(fq_nmod_poly_t out,const fq_nmod_poly_mat_t a,slong row,
                       ulong used,const fq_nmod_poly_t term,const fq_nmod_ctx_t ctx)
{
    if(row==a->r) { fq_nmod_poly_add(out,out,term,ctx);return; }
    fq_nmod_poly_t next;fq_nmod_poly_init(next,ctx);
    for(slong j=0;j<a->c;j++) if(!(used&(UWORD(1)<<j))) {
        fq_nmod_poly_mul(next,term,fq_nmod_poly_mat_entry(a,row,j),ctx);
        unsigned parity=0;
        for(slong k=j+1;k<a->c;k++) parity^=(used>>k)&1;
        if(parity) fq_nmod_poly_neg(next,next,ctx);
        expansion(out,a,row+1,used|(UWORD(1)<<j),next,ctx);
    }
    fq_nmod_poly_clear(next,ctx);
}

static void random_matrix(fq_nmod_poly_mat_t a,flint_rand_t rng,slong len,const fq_nmod_ctx_t ctx)
{
    for(slong i=0;i<a->r;i++) for(slong j=0;j<a->c;j++)
        fq_nmod_poly_randtest(fq_nmod_poly_mat_entry(a,i,j),rng,len,ctx);
}

static void check_equal(const fq_nmod_poly_mat_t a,const fq_nmod_poly_mat_t b,const fq_nmod_ctx_t ctx)
{
    for(slong i=0;i<a->r;i++) for(slong j=0;j<a->c;j++)
        assert(fq_nmod_poly_equal(fq_nmod_poly_mat_entry(a,i,j),fq_nmod_poly_mat_entry(b,i,j),ctx));
}

static void product(fq_nmod_poly_mat_t out,const fq_nmod_poly_mat_t a,
                     const fq_nmod_poly_mat_t b,const fq_nmod_ctx_t ctx)
{
    fq_nmod_poly_t t;fq_nmod_poly_init(t,ctx);
    for(slong i=0;i<a->r;i++) for(slong j=0;j<b->c;j++) {
        fq_nmod_poly_struct *p=fq_nmod_poly_mat_entry(out,i,j);fq_nmod_poly_zero(p,ctx);
        for(slong k=0;k<a->c;k++) {
            fq_nmod_poly_mul(t,fq_nmod_poly_mat_entry(a,i,k),fq_nmod_poly_mat_entry(b,k,j),ctx);
            fq_nmod_poly_add(p,p,t,ctx);
        }
    }
    fq_nmod_poly_clear(t,ctx);
}

static void check_operations(field_ctx_t *field,const fq_nmod_ctx_t ctx,flint_rand_t rng)
{
    fq_nmod_poly_mat_t a,b,ref,got,basis,ann;
    fq_nmod_poly_mat_init(a,4,3,ctx);fq_nmod_poly_mat_init(b,3,4,ctx);
    fq_nmod_poly_mat_init(ref,4,4,ctx);fq_nmod_poly_mat_init(got,4,4,ctx);
    fq_nmod_poly_mat_init(basis,4,4,ctx);fq_nmod_poly_mat_init(ann,4,3,ctx);
    random_matrix(a,rng,39,ctx);random_matrix(b,rng,35,ctx);
    unified_mat_poly_t ua,ub,uc,p;
    unified_mat_poly_init(ua,4,3,field);unified_mat_poly_init(ub,3,4,field);
    unified_mat_poly_init(uc,4,4,field);unified_mat_poly_init(p,4,4,field);
    unified_mat_poly_from_fq(ua,a,ctx);unified_mat_poly_from_fq(ub,b,ctx);
    unified_poly_mat_t old;unified_poly_mat_init(old,4,3,field);
    unified_mat_poly_to_unified(old,ua);unified_mat_poly_from_unified(ua,old);
    fq_nmod_poly_mat_t roundtrip;fq_nmod_poly_mat_init(roundtrip,4,3,ctx);
    unified_mat_poly_to_fq(roundtrip,ua,ctx);check_equal(roundtrip,a,ctx);
    product(ref,a,b,ctx);unified_mat_poly_mul(uc,ua,ub);
    unified_mat_poly_to_fq(got,uc,ctx);check_equal(got,ref,ctx);
    for(slong start=0;start<=80;start+=7) {
        unified_mat_poly_mul_window(uc,ua,ub,start,21);
        unified_mat_poly_to_fq(got,uc,ctx);
        fq_nmod_poly_t t;fq_nmod_poly_init(t,ctx);
        for(slong i=0;i<4;i++) for(slong j=0;j<4;j++) {
            fq_nmod_poly_shift_right(t,fq_nmod_poly_mat_entry(ref,i,j),start,ctx);
            fq_nmod_poly_truncate(t,21,ctx);
            assert(fq_nmod_poly_equal(t,fq_nmod_poly_mat_entry(got,i,j),ctx));
        }
        fq_nmod_poly_clear(t,ctx);
    }
    for(slong order=0;order<=23;order+=order<3 ? 1 : 5) {
        slong shift[4]={-2,3,0,1},original[4]={-2,3,0,1};
        unified_mat_poly_pmbasis(p,shift,ua,order);
        unified_mat_poly_to_fq(basis,p,ctx);product(ann,basis,a,ctx);
        for(slong i=0;i<4;i++) {
            slong degree=LONG_MIN;
            for(slong j=0;j<4;j++) {
                slong d=fq_nmod_poly_degree(fq_nmod_poly_mat_entry(basis,i,j),ctx);
                if(d>=0) degree=FLINT_MAX(degree,d+original[j]);
            }
            assert(degree==shift[i]);
            for(slong j=0;j<3;j++) {
                fq_nmod_poly_truncate(fq_nmod_poly_mat_entry(ann,i,j),order,ctx);
                assert(fq_nmod_poly_is_zero(fq_nmod_poly_mat_entry(ann,i,j),ctx));
            }
        }
        fq_nmod_poly_t determinant,one;
        fq_nmod_poly_init(determinant,ctx);fq_nmod_poly_init(one,ctx);fq_nmod_poly_one(one,ctx);
        expansion(determinant,basis,0,0,one,ctx);
        slong degree=0;for(slong i=0;i<4;i++) degree+=shift[i]-original[i];
        fq_nmod_poly_shift_left(one,one,degree,ctx);
        assert(fq_nmod_poly_equal(determinant,one,ctx));
        fq_nmod_poly_clear(one,ctx);fq_nmod_poly_clear(determinant,ctx);
    }
    /* In-place PMBasis and multiplication exercise ownership of FQ objects. */
    unified_mat_poly_from_fq(p,basis,ctx);
    product(ref,basis,basis,ctx);unified_mat_poly_mul(p,p,p);
    unified_mat_poly_to_fq(got,p,ctx);check_equal(got,ref,ctx);
    unified_mat_poly_mul_window(p,p,p,10000,5);assert(p->length==0);
    unified_mat_poly_from_fq(p,basis,ctx);
    slong alias_shift[4]={0,0,0,0};
    unified_mat_poly_pmbasis(p,alias_shift,p,11);
    unified_mat_poly_to_fq(got,p,ctx);product(ref,got,basis,ctx);
    for(slong i=0;i<4;i++) for(slong j=0;j<4;j++) {
        fq_nmod_poly_truncate(fq_nmod_poly_mat_entry(ref,i,j),11,ctx);
        assert(fq_nmod_poly_is_zero(fq_nmod_poly_mat_entry(ref,i,j),ctx));
    }
    unified_poly_mat_clear(old);fq_nmod_poly_mat_clear(roundtrip,ctx);
    unified_mat_poly_clear(p);unified_mat_poly_clear(uc);unified_mat_poly_clear(ub);unified_mat_poly_clear(ua);
    fq_nmod_poly_mat_clear(ann,ctx);fq_nmod_poly_mat_clear(basis,ctx);fq_nmod_poly_mat_clear(got,ctx);
    fq_nmod_poly_mat_clear(ref,ctx);fq_nmod_poly_mat_clear(b,ctx);fq_nmod_poly_mat_clear(a,ctx);
}

static void check_determinants(field_ctx_t *field,const fq_nmod_ctx_t ctx,flint_rand_t rng)
{
    for(slong n=0;n<=7;n++) {
        fq_nmod_poly_mat_t a,copy;fq_nmod_poly_mat_init(a,n,n,ctx);fq_nmod_poly_mat_init(copy,n,n,ctx);
        unified_mat_poly_t native;unified_mat_poly_init(native,n,n,field);
        unified_poly_t result;unified_poly_init(result,field);
        fq_nmod_poly_t expected,actual,one,content;
        fq_nmod_poly_init(expected,ctx);fq_nmod_poly_init(actual,ctx);
        fq_nmod_poly_init(one,ctx);fq_nmod_poly_one(one,ctx);fq_nmod_poly_init(content,ctx);
        fq_nmod_poly_set_coeff(content,0,one->coeffs,ctx);fq_nmod_poly_set_coeff(content,2,one->coeffs,ctx);
        for(int trial=0;trial<14;trial++) {
            random_matrix(a,rng,trial==0 ? 1 : 2+trial%5,ctx);
            if(trial==1) for(slong i=0;i<n;i++) for(slong j=0;j<i;j++) fq_nmod_poly_zero(fq_nmod_poly_mat_entry(a,i,j),ctx);
            if(trial==2 && n>1) for(slong j=0;j<n;j++) fq_nmod_poly_set(fq_nmod_poly_mat_entry(a,1,j),fq_nmod_poly_mat_entry(a,0,j),ctx);
            if(trial==3 && n) for(slong j=0;j<n;j++) fq_nmod_poly_zero(fq_nmod_poly_mat_entry(a,0,j),ctx);
            /* Nonconstant row content defeats a probabilistic saturation shortcut. */
            if(trial==4) for(slong i=0;i<(n+1)/2;i++) for(slong j=0;j<n;j++)
                fq_nmod_poly_mul(fq_nmod_poly_mat_entry(a,i,j),fq_nmod_poly_mat_entry(a,i,j),content,ctx);
            if(trial==5) for(slong i=0;i<n;i++) for(slong j=0;j<n;j++)
                fq_nmod_poly_shift_left(fq_nmod_poly_mat_entry(a,i,j),fq_nmod_poly_mat_entry(a,i,j),i+2*j,ctx);
            unified_mat_poly_from_fq(native,a,ctx);unified_mat_poly_to_fq(copy,native,ctx);
            fq_nmod_poly_zero(expected,ctx);
            if(n<=5) expansion(expected,a,0,0,one,ctx);
            else {
                fq_nmod_poly_mat_det_set_method(FQ_NMOD_POLY_DET_METHOD_ITER);
                fq_nmod_poly_mat_det_iter(expected,a,ctx);
                fq_nmod_poly_mat_det_set_method(FQ_NMOD_POLY_DET_METHOD_AUTO);
            }
            unified_poly_zero(result);
            int ok=unified_mat_poly_det_hnf(result,native);
            if(ok) {
                unified_to_fq_nmod_poly(actual,result,ctx,field);
                if(!fq_nmod_poly_equal(expected,actual,ctx)) {
                    fprintf(stderr,"det mismatch field=%d n=%ld trial=%d\n",field->field_id,n,trial);abort();
                }
            } else assert(fq_nmod_poly_is_zero(expected,ctx));
            /* Automatic dispatcher must also handle singular fallback. */
            fq_nmod_poly_mat_det_iter(actual,a,ctx);
            assert(fq_nmod_poly_equal(expected,actual,ctx));check_equal(a,copy,ctx);
        }
        fq_nmod_poly_clear(content,ctx);fq_nmod_poly_clear(one,ctx);fq_nmod_poly_clear(actual,ctx);fq_nmod_poly_clear(expected,ctx);
        unified_poly_clear(result);unified_mat_poly_clear(native);fq_nmod_poly_mat_clear(copy,ctx);fq_nmod_poly_mat_clear(a,ctx);
    }
}

static void check_kernel(field_ctx_t *field,const fq_nmod_ctx_t ctx,flint_rand_t rng)
{
    fq_nmod_poly_mat_t a,k,zero;
    fq_nmod_poly_mat_init(a,5,2,ctx);fq_nmod_poly_mat_init(k,3,5,ctx);fq_nmod_poly_mat_init(zero,3,2,ctx);
    random_matrix(a,rng,8,ctx);
    /* Identity block guarantees full column rank and a saturated kernel. */
    for(slong i=0;i<2;i++) for(slong j=0;j<2;j++) {
        fq_nmod_poly_zero(fq_nmod_poly_mat_entry(a,i,j),ctx);
        if(i==j) fq_nmod_poly_one(fq_nmod_poly_mat_entry(a,i,j),ctx);
    }
    unified_mat_poly_t ua,uk;unified_mat_poly_init(ua,5,2,field);unified_mat_poly_init(uk,3,5,field);
    unified_mat_poly_from_fq(ua,a,ctx);assert(unified_mat_poly_kernel(uk,ua));
    unified_mat_poly_to_fq(k,uk,ctx);product(zero,k,a,ctx);
    for(slong i=0;i<3;i++) for(slong j=0;j<2;j++) assert(fq_nmod_poly_is_zero(fq_nmod_poly_mat_entry(zero,i,j),ctx));
    for(slong i=0;i<5;i++) fq_nmod_poly_set(fq_nmod_poly_mat_entry(a,i,1),fq_nmod_poly_mat_entry(a,i,0),ctx);
    unified_mat_poly_from_fq(ua,a,ctx);assert(!unified_mat_poly_kernel(uk,ua));
    unified_mat_poly_clear(uk);unified_mat_poly_clear(ua);
    fq_nmod_poly_mat_clear(zero,ctx);fq_nmod_poly_mat_clear(k,ctx);fq_nmod_poly_mat_clear(a,ctx);
}

static void check_byte_tails(field_ctx_t *field)
{
    if(field->field_id!=FIELD_ID_GF24 && field->field_id!=FIELD_ID_GF28) return;
    uint8_t tables[256][32],src[80],got[80],expected[80];
    int order=field->field_id==FIELD_ID_GF24 ? 16 : 256;
    const uint8_t *(*row)(uint8_t)=order==16 ? gf24_get_scalar_row : gf28_get_scalar_row;
    unified_byte_tables(tables,order,row);
    for(int scalar=0;scalar<order;scalar++) for(int n=0;n<=67;n++) for(int offset=0;offset<4;offset++) {
        for(int i=0;i<80;i++) { src[i]=(17*i+scalar)&(order-1);got[i]=expected[i]=i&(order-1); }
        for(int i=offset;i<offset+n;i++) expected[i]^=row(scalar)[src[i]];
        unified_byte_addmul(got+offset,src+offset,n,tables[scalar]);
        assert(!memcmp(got,expected,sizeof(got)));
    }
}

static void check(ulong prime,slong degree,ulong modulus_bits,int force_generic)
{
    fq_nmod_ctx_t ctx;
    if(modulus_bits) {
        nmod_poly_t modulus;nmod_poly_init(modulus,prime);
        for(slong i=0;i<FLINT_BITS && i<=degree;i++) nmod_poly_set_coeff_ui(modulus,i,(modulus_bits>>i)&1);
        nmod_poly_set_coeff_ui(modulus,degree,1);
        fq_nmod_ctx_init_modulus(ctx,modulus,"a");nmod_poly_clear(modulus);
    } else { fmpz_t p;fmpz_init_set_ui(p,prime);fq_nmod_ctx_init(ctx,p,degree,"a");fmpz_clear(p); }
    field_ctx_t field;field_ctx_init_enhanced(&field,ctx,force_generic ? 0 : 65536);
    if(force_generic) assert(field.field_id==FIELD_ID_FQ);
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,12,257);
    check_byte_tails(&field);check_operations(&field,ctx,rng);
    check_kernel(&field,ctx,rng);check_determinants(&field,ctx,rng);
    printf("Unified HNF checked: p=%lu degree=%ld backend=%d\n",prime,degree,field.field_id);
    flint_rand_clear(rng);field_ctx_clear(&field);fq_nmod_ctx_clear(ctx);
}

int main(void)
{
    check(2,4,0x13,0);check(2,8,0x11d,0);check(2,8,0x11b,0);
    check(3,2,0,0);check(3,2,0,1);check(257,1,0,0);
    check(2,16,0x1002d,0);check(2,32,UWORD(0x100008299),0);
    check(2,64,UWORD(0x247f43cb7),0);check(2,128,0x87,0);
    const slong arbitrary[]={2,7,9,15,17,31,33,63,65,95,127};
    for(unsigned i=0;i<sizeof(arbitrary)/sizeof(*arbitrary);i++) check(2,arbitrary[i],0,0);
    check(2,4,0x19,0);check(2,16,0x1002b,0);check(2,64,0x1b,0);
    flint_cleanup_master();return 0;
}

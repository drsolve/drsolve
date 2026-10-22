/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "fq_unified_interface.h"
#include <assert.h>
#include <flint/nmod_poly_factor.h>
#include <omp.h>

static void check(const fq_nmod_ctx_t fq, flint_rand_t rng)
{
    gf2n_native_ctx *c=gf2n_native_create(fq);
    field_ctx_t f; field_ctx_init(&f,fq);
    assert(f.field_id==FIELD_ID_NMOD || field_id_is_gf2n(f.field_id));
    fq_nmod_t a,b,ref,out; fq_nmod_init(a,fq); fq_nmod_init(b,fq);
    fq_nmod_init(ref,fq); fq_nmod_init(out,fq);
    const void *ar=f.field_id==FIELD_ID_NMOD ? (const void *)&f.ctx.nmod_ctx : (const void *)f.ctx.fq_ctx;
    for(int i=0;i<96;i++) {
        fq_nmod_randtest(a,rng,fq); fq_nmod_randtest(b,rng,fq);
        if(i==0) fq_nmod_zero(a,fq);
        if(i==1) fq_nmod_one(b,fq);
        gf2128_t x=gf2n_native_from_fq(a),y=gf2n_native_from_fq(b);
        field_elem_u u,v,w; fq_nmod_to_field_elem(&u,a,&f); fq_nmod_to_field_elem(&v,b,&f);
        fq_nmod_mul(ref,a,b,fq);
        for(int software=0;software<2;software++) {
            c->clmul=software ? 0 : has_pclmulqdq();
            gf2n_native_to_fq(out,gf2n_native_mul(x,y,c),fq);
            assert(fq_nmod_equal(out,ref,fq));
            gf2n_native_to_fq(out,gf2n_native_mul_fast(x,y,c),fq);
            assert(fq_nmod_equal(out,ref,fq));
        }
        field_mul(&w,&u,&v,f.field_id,ar); field_elem_to_fq_nmod(out,&w,&f);
        assert(fq_nmod_equal(out,ref,fq));
        field_mul(&u,&u,&v,f.field_id,ar); field_elem_to_fq_nmod(out,&u,&f);
        assert(fq_nmod_equal(out,ref,fq));
        fq_nmod_to_field_elem(&u,a,&f);
        fq_nmod_add(ref,a,b,fq); field_add(&w,&u,&v,f.field_id,ar);
        field_elem_to_fq_nmod(out,&w,&f); assert(fq_nmod_equal(out,ref,fq));
        c->clmul=(i&1) ? has_pclmulqdq() : 0;
        fq_nmod_sqr(ref,a,fq); gf2n_native_to_fq(out,gf2n_native_sqr(x,c),fq);
        assert(fq_nmod_equal(out,ref,fq));
        if(!fq_nmod_is_zero(a,fq)) {
            fq_nmod_inv(ref,a,fq); gf2n_native_to_fq(out,gf2n_native_inv(x,c),fq);
            assert(fq_nmod_equal(out,ref,fq));
            field_inv(&w,&u,f.field_id,ar); field_elem_to_fq_nmod(out,&w,&f);
            assert(fq_nmod_equal(out,ref,fq));
        }
    }
    fq_nmod_clear(a,fq); fq_nmod_clear(b,fq); fq_nmod_clear(ref,fq); fq_nmod_clear(out,fq);
    field_ctx_clear(&f); gf2n_native_destroy(c);
}
int main(void)
{
    /* Independent random irreducible moduli, including every awkward width.
     * OpenMP also checks that different contexts do not share mutable state. */
    #pragma omp parallel for schedule(dynamic) num_threads(4)
    for(int n=1;n<=128;n++) {
        flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng,311+n,991+n);
        nmod_poly_t p; nmod_poly_init(p,2);
        for(int trial=0;trial<2;trial++) {
            nmod_poly_randtest_monic_irreducible(p,rng,n+1);
            fq_nmod_ctx_t fq; fq_nmod_ctx_init_modulus(fq,p,"t"); check(fq,rng); fq_nmod_ctx_clear(fq);
        }
        nmod_poly_clear(p); flint_rand_clear(rng);
    }
    flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng,1,2);
    const unsigned ns[]={4,8,8,16,32,64,64,128};
    const uint64_t ps[]={0x13,0x11d,0x11b,0x1002d,0x100008299,UINT64_C(0x247f43cb7),0x1b,0x87};
    for(unsigned k=0;k<8;k++) {
        nmod_poly_t p; nmod_poly_init(p,2);
        for(unsigned j=0;j<64 && j<ns[k];j++) if((ps[k]>>j)&1) nmod_poly_set_coeff_ui(p,j,1);
        nmod_poly_set_coeff_ui(p,ns[k],1);
        fq_nmod_ctx_t fq; fq_nmod_ctx_init_modulus(fq,p,"t"); check(fq,rng);
        if(ns[k]>=32) {
            gf232_mul_func m32=gf232_mul; gf264_mul_func m64=gf264_mul; gf2128_mul_func m128=gf2128_mul;
            gf232_mul=gf232_mul_software; gf264_mul=gf264_mul_software; gf2128_mul=gf2128_mul_software;
            check(fq,rng);
            gf232_mul=m32; gf264_mul=m64; gf2128_mul=m128;
        }
        /* Concurrent initialization and use of the same modulus must also
         * work, in addition to the different-modulus workers above. */
        if(ns[k]==8 || ns[k]==16) {
            #pragma omp parallel for num_threads(4)
            for(int j=0;j<4;j++) {
                flint_rand_t local; flint_rand_init(local); flint_rand_set_seed(local,j+55,k+77);
                check(fq,local); flint_rand_clear(local);
            }
        }
        fq_nmod_ctx_clear(fq); nmod_poly_clear(p);
    }
    /* Exercise the two-fold 128-bit reduction with dense lower halves,
     * including bit 63, rather than only the usual sparse 0x87 modulus. */
    for (int trial=0; trial<3; ++trial) {
        nmod_poly_t modulus; nmod_poly_init(modulus,2);
        do {
            nmod_poly_randtest(modulus,rng,64);
            nmod_poly_set_coeff_ui(modulus,0,1);
            nmod_poly_set_coeff_ui(modulus,63,1);
            nmod_poly_set_coeff_ui(modulus,128,1);
        } while (!nmod_poly_is_irreducible(modulus));
        fq_nmod_ctx_t fq; fq_nmod_ctx_init_modulus(fq,modulus,"t");
        check(fq,rng);
        fq_nmod_ctx_clear(fq); nmod_poly_clear(modulus);
    }
    /* Keep two distinct representations of the same degree live together.
     * Matching moduli share immutable tables; distinct moduli never do. */
    fq_nmod_ctx_t pair[2];
    field_ctx_t live[2];
    nmod_poly_t p; nmod_poly_init(p,2);
    for(int k=0;k<2;k++) {
        unsigned bits=k ? 0x11b : 0x11d;
        nmod_poly_zero(p);
        for(unsigned j=0;j<=8;j++) if((bits>>j)&1) nmod_poly_set_coeff_ui(p,j,1);
        fq_nmod_ctx_init_modulus(pair[k],p,"t"); field_ctx_init(live+k,pair[k]);
    }
    #pragma omp parallel for num_threads(4)
    for(int j=0;j<8;j++) {
        flint_rand_t local; flint_rand_init(local); flint_rand_set_seed(local,j+71,991);
        check(pair[j%2],local); flint_rand_clear(local);
    }
    for(int k=0;k<2;k++) { field_ctx_clear(live+k); fq_nmod_ctx_clear(pair[k]); }
    nmod_poly_clear(p);
    fq_nmod_ctx_t large; fq_nmod_ctx_init_ui(large,2,129,"t");
    field_ctx_t fallback; field_ctx_init(&fallback,large); assert(fallback.field_id==FIELD_ID_FQ);
    field_ctx_clear(&fallback); fq_nmod_ctx_clear(large);
    flint_rand_clear(rng); flint_cleanup_master();
    puts("GF(2^n): all 1<=n<=128, two random moduli each, hardware/software and unified arithmetic PASS");
    return 0;
}

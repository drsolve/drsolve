/* Exact unified API resultants, including non-monic scalars and aliases. */
#include <assert.h>
#include "unified_mpoly_resultant.h"
#include "dixon_flint.h"
#include "unified_dixon.h"

static void check(field_ctx_t *field)
{
    const void *fc=field_ctx_data(field);
    unified_mpoly_ctx_t ctx=unified_mpoly_ctx_init(2,ORD_LEX,field);
    unified_mpoly_t a=unified_mpoly_init(ctx),b=unified_mpoly_init(ctx);
    unified_mpoly_t ref=unified_mpoly_init(ctx),got=unified_mpoly_init(ctx);
    unified_mpoly_t tmp=unified_mpoly_init(ctx);
    field_elem_u c;field_init_elem(&c,field->field_id,fc);
    for(slong var=0;var<2;var++) {
        unified_mpoly_zero(a);unified_mpoly_zero(b);
        for(int k=0;k<2;k++) for(ulong x=0;x<=(ulong)(3-k);x++)
            for(ulong y=0;y<=(ulong)(3-k)-x;y++) {
                ulong e[2];e[var]=x;e[1-var]=y;
                ulong value=2+x+3*y+k;
                if(field->field_id==FIELD_ID_NMOD)c.nmod=value%field->ctx.nmod_ctx.n;
                else if(field->field_id==FIELD_ID_FMPZ_MOD)fmpz_set_ui(&c.big,value);
                else {
                    const fq_nmod_ctx_struct *fq=field->field_id==FIELD_ID_FQ_ZECH
                        ? field->ctx.zech_ctx->fq_nmod_ctx : field->ctx.fq_ctx;
                    fq_nmod_t z;fq_nmod_init(z,fq);
                    fq_nmod_set_ui(z,value,fq);
                    nmod_poly_set_coeff_ui(z,1,1);
                    fq_nmod_to_field_elem(&c,z,field);fq_nmod_clear(z,fq);
                }
                unified_mpoly_set_coeff_ui(k ? b : a,&c,e);
            }
        g_dixon_bezout_fast=0;
        assert(unified_mpoly_resultant(ref,a,b,var,ctx));
        g_dixon_bezout_fast=1;
        assert(unified_mpoly_resultant_bezout_try(got,a,b,var,ctx));
        assert(unified_mpoly_equal(got,ref));
        assert(unified_mpoly_resultant(got,a,b,var,ctx));
        assert(unified_mpoly_equal(got,ref));
        if(var==0) {
            unified_mpoly_ctx_t output=unified_mpoly_ctx_init(1,ORD_LEX,field);
            unified_mpoly_t reduced=unified_mpoly_init(output),expected=unified_mpoly_init(output);
            unified_mpoly_t pair[2]={a,b};const slong map[2]={0,0};
            unified_mpoly_reindex(expected,ref,map);
            assert(unified_dixon_resultant(reduced,pair,1));
            assert(unified_mpoly_equal(reduced,expected));
            unified_mpoly_clear(reduced);unified_mpoly_clear(expected);unified_mpoly_ctx_clear(output);
        }
        unified_mpoly_set(tmp,a);
        assert(unified_mpoly_resultant(tmp,tmp,b,var,ctx));
        assert(unified_mpoly_equal(tmp,ref));
        unified_mpoly_set(tmp,b);
        assert(unified_mpoly_resultant(tmp,a,tmp,var,ctx));
        assert(unified_mpoly_equal(tmp,ref));
        /* Shared factor: exact zero, not a nonzero maximal-rank minor. */
        assert(unified_mpoly_resultant(got,a,a,var,ctx));
        assert(unified_mpoly_is_zero(got));
        /* Rejection must leave an existing output untouched. */
        unified_mpoly_set(got,ref);unified_mpoly_zero(tmp);
        assert(!unified_mpoly_resultant_bezout_try(got,tmp,b,var,ctx));
        assert(unified_mpoly_equal(got,ref));
    }
    field_clear_elem(&c,field->field_id,fc);
    unified_mpoly_clear(a);unified_mpoly_clear(b);unified_mpoly_clear(tmp);
    unified_mpoly_clear(ref);unified_mpoly_clear(got);unified_mpoly_ctx_clear(ctx);
}

int main(void)
{
    g_dixon_verbose_level=0;
    fmpz_t p;fmpz_init(p);
    const char *primes[]={"7","65537","101","18446744073709551557",
                         "18446744073709551629","170141183460469231731687303715884105727"};
    for(size_t i=0;i<sizeof(primes)/sizeof(*primes);i++) {
        fmpz_set_str(p,primes[i],10);
        field_ctx_t f;field_ctx_init_prime(&f,p);check(&f);field_ctx_clear(&f);
    }
    /* An arbitrary-precision NTT-friendly prime, p == 1 mod 64. */
    fmpz_one(p);fmpz_mul_2exp(p,p,80);fmpz_add_ui(p,p,1);
    while(!fmpz_is_prime(p))fmpz_add_ui(p,p,64);
    field_ctx_t big;field_ctx_init_prime(&big,p);check(&big);field_ctx_clear(&big);
    const ulong chars[]={2,2,3,5,65537};
    const slong degrees[]={4,8,2,2,2};
    for(size_t i=0;i<sizeof(chars)/sizeof(*chars);i++) {
        fmpz_set_ui(p,chars[i]);fq_nmod_ctx_t fq;fq_nmod_ctx_init(fq,p,degrees[i],"t");
        field_ctx_t f;field_ctx_init(&f,fq);check(&f);field_ctx_clear(&f);
        fq_nmod_ctx_clear(fq);
    }
    /* Exercise the Zech adapter even if automatic Zech selection is off. */
    field_ctx_t zech={0};zech.field_id=FIELD_ID_FQ_ZECH;
    zech.ctx.zech_ctx=flint_malloc(sizeof(fq_zech_ctx_struct));
    fq_zech_ctx_init_ui(zech.ctx.zech_ctx,3,2,"t");
    check(&zech);field_ctx_clear(&zech);
    fmpz_clear(p);
    puts("Unified Dixon/Bezout fast exact coefficients, fields and aliases PASS");
    return 0;
}

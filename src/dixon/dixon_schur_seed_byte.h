/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Compact GF(16)/GF(256) matrices. Windows borrow storage; LU keeps FLINT's
 * packed L convention so seed and completion pivot indices remain identical.
 * Arithmetic tables belong to the call, not a process-global active field. */
#include "../solver/quotient_byte.h"
typedef struct {
    field_ctx_t field;
    uint8_t tables[256][32];
    const uint8_t *(*row)(uint8_t);
    uint8_t (*inverse)(uint8_t);
} seed_byte_field;
typedef struct { slong r,c,stride; uint8_t *data; } seed_byte_mat_struct;
typedef seed_byte_mat_struct seed_byte_mat_t[1];
#define SEED_BYTE_ENTRY(a,i,j) ((a)->data[(i)*(a)->stride+(j)])
static uint8_t seed_byte_coeff(const fq_nmod_t value,const seed_byte_field *byte)
{
    field_elem_u converted;
    fq_nmod_to_field_elem(&converted,value,&byte->field);
    return byte->field.field_id == FIELD_ID_GF24 ? converted.gf24 : converted.gf28;
}
static void seed_byte_init(seed_byte_mat_t a,slong r,slong c)
{
    a->r=r; a->c=a->stride=c;
    a->data=flint_calloc(FLINT_MAX((size_t)r*c,1),1);
}
static void seed_byte_clear(seed_byte_mat_t a,const seed_byte_field *byte)
{ (void)byte; flint_free(a->data); }
static void seed_byte_window_init(seed_byte_mat_t a,const seed_byte_mat_t b,
    slong r0,slong c0,slong r1,slong c1,const seed_byte_field *byte)
{
    (void)byte; a->r=r1-r0; a->c=c1-c0; a->stride=b->stride;
    a->data=a->r && a->c ? b->data+r0*b->stride+c0 : b->data;
}
static void seed_byte_window_clear(seed_byte_mat_t a,const seed_byte_field *byte)
{ (void)a; (void)byte; }
static slong seed_byte_lu(slong *perm,seed_byte_mat_t a,int rank_check,
    const seed_byte_field *byte)
{
    (void)rank_check;
    slong rank=0;
    for(slong i=0;i<a->r;i++)perm[i]=i;
    for(slong col=0;col<a->c && rank<a->r;col++) {
        slong pivot=rank;
        while(pivot<a->r && !SEED_BYTE_ENTRY(a,pivot,col))pivot++;
        if(pivot==a->r)continue;
        if(pivot!=rank) {
            slong swap=perm[pivot];perm[pivot]=perm[rank];perm[rank]=swap;
            for(slong j=0;j<a->c;j++) {
                uint8_t v=SEED_BYTE_ENTRY(a,pivot,j);
                SEED_BYTE_ENTRY(a,pivot,j)=SEED_BYTE_ENTRY(a,rank,j);
                SEED_BYTE_ENTRY(a,rank,j)=v;
            }
        }
        const uint8_t *inv=byte->row(byte->inverse(SEED_BYTE_ENTRY(a,rank,col)));
        for(slong i=rank+1;i<a->r;i++) {
            uint8_t factor=inv[SEED_BYTE_ENTRY(a,i,col)];
            if(factor)quotient_byte_addmul(a->data+i*a->stride+col+1,
                a->data+rank*a->stride+col+1,a->c-col-1,byte->tables[factor]);
            SEED_BYTE_ENTRY(a,i,col)=0;
            SEED_BYTE_ENTRY(a,i,rank)=factor;
        }
        rank++;
    }
    return rank;
}
static int seed_byte_solve(seed_byte_mat_t x,const seed_byte_mat_t a,
    const seed_byte_mat_t b,const seed_byte_field *byte)
{
    seed_byte_mat_t lu;seed_byte_init(lu,a->r,a->c);
    for(slong i=0;i<a->r;i++)memcpy(lu->data+i*lu->stride,a->data+i*a->stride,a->c);
    slong *perm=flint_malloc(FLINT_MAX(a->r,1)*sizeof(slong));
    int ok=seed_byte_lu(perm,lu,0,byte)==a->r;
    if(ok) {
        for(slong i=0;i<x->r;i++)memcpy(x->data+i*x->stride,b->data+perm[i]*b->stride,x->c);
        for(slong i=0;i<x->r;i++)for(slong j=0;j<i;j++) {
            uint8_t scalar=SEED_BYTE_ENTRY(lu,i,j);
            if(scalar)quotient_byte_addmul(x->data+i*x->stride,x->data+j*x->stride,x->c,byte->tables[scalar]);
        }
        for(slong i=x->r;i-- >0;) {
            for(slong j=i+1;j<x->r;j++) {
                uint8_t scalar=SEED_BYTE_ENTRY(lu,i,j);
                if(scalar)quotient_byte_addmul(x->data+i*x->stride,x->data+j*x->stride,x->c,byte->tables[scalar]);
            }
            const uint8_t *inv=byte->row(byte->inverse(SEED_BYTE_ENTRY(lu,i,i)));
            for(slong j=0;j<x->c;j++)SEED_BYTE_ENTRY(x,i,j)=inv[SEED_BYTE_ENTRY(x,i,j)];
        }
    }
    flint_free(perm);seed_byte_clear(lu,byte);return ok;
}
static void seed_byte_mul(seed_byte_mat_t out,const seed_byte_mat_t a,
    const seed_byte_mat_t b,const seed_byte_field *byte)
{
    /* Independent output rows; SIMD handles the contiguous column direction. */
    #pragma omp parallel for schedule(static) if(a->r>128 && b->c>64)
    for(slong i=0;i<out->r;i++) {
        uint8_t *row=out->data+i*out->stride;
        memset(row,0,out->c);
        for(slong k=0;k<a->c;k++) {
            uint8_t scalar=SEED_BYTE_ENTRY(a,i,k);
            if(scalar)quotient_byte_addmul(row,b->data+k*b->stride,b->c,byte->tables[scalar]);
        }
    }
}
static void seed_byte_sub(seed_byte_mat_t out,const seed_byte_mat_t a,
    const seed_byte_mat_t b,const seed_byte_field *byte)
{
    (void)byte;
    for(slong i=0;i<out->r;i++)for(slong j=0;j<out->c;j++)
        SEED_BYTE_ENTRY(out,i,j)=SEED_BYTE_ENTRY(a,i,j)^SEED_BYTE_ENTRY(b,i,j);
}

#define SCHUR_SEED_EXTENSION
#define SCHUR_SEED_EXTRA , const seed_byte_field *byte
#define SCHUR_SEED_NAME dixon_schur_seed_byte
#define SEED_MAT seed_byte_mat_t
#define SEED_INIT(a,r,c) seed_byte_init(a,r,c)
#define SEED_CALL(fn,...) seed_byte_##fn(__VA_ARGS__,byte)
#define SEED_ENTRY(a,i,j) SEED_BYTE_ENTRY(a,i,j)
#define SEED_ZERO(a,i,j) (!SEED_ENTRY(a,i,j))
#define SEED_STORE(a,i,j,v) (SEED_ENTRY(a,i,j)=seed_byte_coeff(v,byte))
#define SEED_ACCUMULATE(a,i,j,v) (SEED_ENTRY(a,i,j)^=seed_byte_coeff(v,byte))
#include "dixon_schur_seed_template.h"
#undef SCHUR_SEED_EXTENSION
#undef SCHUR_SEED_EXTRA
#undef SCHUR_SEED_NAME
#undef SEED_MAT
#undef SEED_INIT
#undef SEED_CALL
#undef SEED_ENTRY
#undef SEED_ZERO
#undef SEED_STORE
#undef SEED_ACCUMULATE
#undef SEED_BYTE_ENTRY

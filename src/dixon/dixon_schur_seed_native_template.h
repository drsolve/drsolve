/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Packed characteristic-two scalar matrices; shared seed algorithm. */
typedef struct { slong r,c,stride; SEED_E *data; } SEED_NAME(mat_struct);
typedef SEED_NAME(mat_struct) SEED_NAME(mat_t)[1];
#define SEED_ENTRY(a,i,j) ((a)->data[(i)*(a)->stride+(j)])
static void SEED_NAME(init)(SEED_NAME(mat_t) a,slong r,slong c)
{
    a->r=r; a->c=a->stride=c;
    a->data=flint_calloc(FLINT_MAX((size_t)r*c,1),sizeof(SEED_E));
}
static void SEED_NAME(clear)(SEED_NAME(mat_t) a,const gf2n_native_ctx *field)
{ (void)field; flint_free(a->data); }
static void SEED_NAME(window_init)(SEED_NAME(mat_t) a,const SEED_NAME(mat_t) b,
    slong r0,slong c0,slong r1,slong c1,const gf2n_native_ctx *field)
{
    (void)field; a->r=r1-r0; a->c=c1-c0; a->stride=b->stride;
    a->data=a->r && a->c ? b->data+r0*b->stride+c0 : b->data;
}
static void SEED_NAME(window_clear)(SEED_NAME(mat_t) a,const gf2n_native_ctx *field)
{ (void)a; (void)field; }
static slong SEED_NAME(lu)(slong *perm,SEED_NAME(mat_t) a,int rank_check,
    const gf2n_native_ctx *field)
{
    (void)rank_check;
    slong rank=0;
    for(slong i=0;i<a->r;i++)perm[i]=i;
    for(slong col=0;col<a->c && rank<a->r;col++) {
        slong pivot=rank;
        while(pivot<a->r && SEED_ISZERO(SEED_ENTRY(a,pivot,col)))pivot++;
        if(pivot==a->r)continue;
        if(pivot!=rank) {
            slong swap=perm[pivot];perm[pivot]=perm[rank];perm[rank]=swap;
            for(slong j=0;j<a->c;j++) {
                SEED_E v=SEED_ENTRY(a,pivot,j);
                SEED_ENTRY(a,pivot,j)=SEED_ENTRY(a,rank,j);
                SEED_ENTRY(a,rank,j)=v;
            }
        }
        SEED_E inv=SEED_INV(SEED_ENTRY(a,rank,col));
        for(slong i=rank+1;i<a->r;i++) {
            SEED_E factor=SEED_MUL(inv,SEED_ENTRY(a,i,col));
            if(!SEED_ISZERO(factor)) for(slong j=col+1;j<a->c;j++)
                SEED_ENTRY(a,i,j)=SEED_ADD(SEED_ENTRY(a,i,j),SEED_MUL(factor,SEED_ENTRY(a,rank,j)));
            SEED_ENTRY(a,i,col)=SEED_FROM_WORD(0);
            SEED_ENTRY(a,i,rank)=factor;
        }
        rank++;
    }
    return rank;
}
static int SEED_NAME(solve)(SEED_NAME(mat_t) x,const SEED_NAME(mat_t) a,
    const SEED_NAME(mat_t) b,const gf2n_native_ctx *field)
{
    SEED_NAME(mat_t) lu;SEED_NAME(init)(lu,a->r,a->c);
    for(slong i=0;i<a->r;i++)memcpy(lu->data+i*lu->stride,a->data+i*a->stride,a->c*sizeof(SEED_E));
    slong *perm=flint_malloc(FLINT_MAX(a->r,1)*sizeof(slong));
    int ok=SEED_NAME(lu)(perm,lu,0,field)==a->r;
    if(ok) {
        for(slong i=0;i<x->r;i++)memcpy(x->data+i*x->stride,b->data+perm[i]*b->stride,x->c*sizeof(SEED_E));
        for(slong i=0;i<x->r;i++)for(slong j=0;j<i;j++) {
            SEED_E scalar=SEED_ENTRY(lu,i,j);
            if(!SEED_ISZERO(scalar)) for(slong k=0;k<x->c;k++)
                SEED_ENTRY(x,i,k)=SEED_ADD(SEED_ENTRY(x,i,k),SEED_MUL(scalar,SEED_ENTRY(x,j,k)));
        }
        for(slong i=x->r;i-- >0;) {
            for(slong j=i+1;j<x->r;j++) {
                SEED_E scalar=SEED_ENTRY(lu,i,j);
                if(!SEED_ISZERO(scalar)) for(slong k=0;k<x->c;k++)
                SEED_ENTRY(x,i,k)=SEED_ADD(SEED_ENTRY(x,i,k),SEED_MUL(scalar,SEED_ENTRY(x,j,k)));
            }
            SEED_E inv=SEED_INV(SEED_ENTRY(lu,i,i));
            for(slong j=0;j<x->c;j++)SEED_ENTRY(x,i,j)=SEED_MUL(inv,SEED_ENTRY(x,i,j));
        }
    }
    flint_free(perm);SEED_NAME(clear)(lu,field);return ok;
}
static void SEED_NAME(mul)(SEED_NAME(mat_t) out,const SEED_NAME(mat_t) a,
    const SEED_NAME(mat_t) b,const gf2n_native_ctx *field)
{
    /* Independent rows, contiguous updates, immutable field context. */
    #pragma omp parallel for schedule(static) if(a->r>128 && b->c>64)
    for(slong i=0;i<out->r;i++) {
        SEED_E *row=out->data+i*out->stride;
        memset(row,0,out->c*sizeof(SEED_E));
        for(slong k=0;k<a->c;k++) {
            SEED_E scalar=SEED_ENTRY(a,i,k);
            if(!SEED_ISZERO(scalar)) for(slong j=0;j<b->c;j++)
                row[j]=SEED_ADD(row[j],SEED_MUL(scalar,SEED_ENTRY(b,k,j)));
        }
    }
}
static void SEED_NAME(sub)(SEED_NAME(mat_t) out,const SEED_NAME(mat_t) a,
    const SEED_NAME(mat_t) b,const gf2n_native_ctx *field)
{
    (void)field;
    for(slong i=0;i<out->r;i++)for(slong j=0;j<out->c;j++)
        SEED_ENTRY(out,i,j)=SEED_ADD(SEED_ENTRY(a,i,j),SEED_ENTRY(b,i,j));
}


#define SEED_NATIVE_EVAL
#define SCHUR_SEED_EXTENSION
#define SCHUR_SEED_EXTRA , const gf2n_native_ctx *field
#define SCHUR_SEED_NAME SEED_NAME(select)
#define SEED_MAT SEED_NAME(mat_t)
#define SEED_INIT(a,r,c) SEED_NAME(init)(a,r,c)
#define SEED_CALL(fn,...) SEED_NAME(fn)(__VA_ARGS__,field)
#define SEED_ZERO(a,i,j) SEED_ISZERO(SEED_ENTRY(a,i,j))
#define SEED_STORE(a,i,j,v) (SEED_ENTRY(a,i,j)=SEED_FROM_FQ(v))
#define SEED_ACCUMULATE(a,i,j,v) (SEED_ENTRY(a,i,j)=SEED_ADD(SEED_ENTRY(a,i,j),SEED_FROM_FQ(v)))
#include "dixon_schur_seed_template.h"
#undef SEED_NATIVE_EVAL
#undef SCHUR_SEED_EXTENSION
#undef SCHUR_SEED_EXTRA
#undef SCHUR_SEED_NAME
#undef SEED_MAT
#undef SEED_INIT
#undef SEED_CALL
#undef SEED_ZERO
#undef SEED_STORE
#undef SEED_ACCUMULATE
#undef SEED_ENTRY
#undef SEED_NAME
#undef SEED_E
#undef SEED_MUL
#undef SEED_INV
#undef SEED_ADD
#undef SEED_ISZERO
#undef SEED_FROM_WORD
#undef SEED_FROM_FQ

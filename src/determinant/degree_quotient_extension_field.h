/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Unified-field instantiation. The context is immutable during sampling;
 * matrices, scalar temporaries and normal forms are owned by each worker. */
#include "fq_unified_interface.h"
#include "solver/quotient_byte.h"
#include <assert.h>

typedef field_elem_u dq_coeff;
static field_ctx_t extension_field;
static const void *extension_arithmetic;
static uint8_t extension_byte_tables[256][32];
static void dq_field_setup(const fq_nmod_ctx_t ctx)
{
    field_ctx_init(&extension_field,ctx);
    extension_arithmetic=extension_field.field_id==FIELD_ID_FQ_ZECH ?
        (const void *)extension_field.ctx.zech_ctx : (const void *)extension_field.ctx.fq_ctx;
    if(extension_field.field_id==FIELD_ID_GF24)
        quotient_byte_tables(extension_byte_tables,16,gf24_get_scalar_row);
    if(extension_field.field_id==FIELD_ID_GF28)
        quotient_byte_tables(extension_byte_tables,256,gf28_get_scalar_row);
}
#define DQ_FIELD extension_field.field_id, extension_arithmetic
/* Account for the separately allocated polynomial representation as well. */
#define DQ_COEFF_BYTES (sizeof(dq_coeff)+(extension_field.field_id==FIELD_ID_FQ ? \
    fq_nmod_ctx_degree(extension_field.ctx.fq_ctx)*sizeof(ulong) : 0))
#define dq_init(a) field_init_elem(&(a),DQ_FIELD)
#define dq_clear(a) field_clear_elem(&(a),DQ_FIELD)
#define dq_set(a,b) field_set_elem(&(a),&(b),DQ_FIELD)
#define dq_one(a) field_set_one(&(a),DQ_FIELD)
#define dq_nonzero(a) (!field_is_zero(&(a),DQ_FIELD))
#define dq_neg(a,b) field_neg(&(a),&(b),DQ_FIELD)
#define dq_add(a,b,c) field_add(&(a),&(b),&(c),DQ_FIELD)

static dq_coeff *dq_vec_new(slong count)
{
    dq_coeff *v=flint_calloc(count,sizeof(*v));
    for(slong i=0;i<count;i++) {
        dq_init(v[i]);field_set_zero(v+i,DQ_FIELD);
    }
    return v;
}
static void dq_vec_free(dq_coeff *v,slong count)
{
    for(slong i=0;i<count;i++)dq_clear(v[i]);
    flint_free(v);
}
static void dq_vec_zero(dq_coeff *v,slong count)
{
    for(slong i=0;i<count;i++)field_set_zero(v+i,DQ_FIELD);
}
static void dq_addmul_impl(dq_coeff *out,const dq_coeff *a,const dq_coeff *b)
{
    dq_coeff t;dq_init(t);
    field_mul(&t,a,b,DQ_FIELD);field_add(out,out,&t,DQ_FIELD);dq_clear(t);
}
#define dq_addmul(out,a,b) dq_addmul_impl(&(out),&(a),&(b))
static void dq_pow_impl(dq_coeff *out,const dq_coeff *base,ulong power)
{
    dq_coeff b;dq_init(b);dq_set(b,*base);dq_one(*out);
    while(power) {
        if(power&1)field_mul(out,out,&b,DQ_FIELD);
        power>>=1;if(power)field_mul(&b,&b,&b,DQ_FIELD);
    }
    dq_clear(b);
}
#define dq_pow(out,a,e) dq_pow_impl(&(out),&(a),(e))

typedef struct {slong r,c;dq_coeff *entries;} dq_mat_struct;
typedef dq_mat_struct dq_mat_t[1];
#define dq_mat_entry(m,i,j) ((m)->entries[(i)*(m)->c+(j)])
static void dq_mat_init(dq_mat_struct *m,slong r,slong c,ulong unused)
{
    (void)unused;m->r=r;m->c=c;m->entries=dq_vec_new(r*c);
}
static void dq_mat_clear(dq_mat_struct *m)
{
    dq_vec_free(m->entries,m->r*m->c);memset(m,0,sizeof(*m));
}
static void dq_mat_zero(dq_mat_struct *m) {dq_vec_zero(m->entries,m->r*m->c);}
static void dq_mat_set(dq_mat_struct *out,const dq_mat_struct *in)
{
    assert(out->r==in->r && out->c==in->c);
    for(slong i=0;i<out->r*out->c;i++)dq_set(out->entries[i],in->entries[i]);
}
static void dq_mat_transpose(dq_mat_struct *out,const dq_mat_struct *in)
{
    assert(out!=in && out->r==in->c && out->c==in->r);
    for(slong i=0;i<out->r;i++)for(slong j=0;j<out->c;j++)
        dq_set(dq_mat_entry(out,i,j),dq_mat_entry(in,j,i));
}
static void dq_mat_add(dq_mat_struct *out,const dq_mat_struct *a,const dq_mat_struct *b)
{
    for(slong i=0;i<out->r*out->c;i++)dq_add(out->entries[i],a->entries[i],b->entries[i]);
}
static void dq_mat_sub(dq_mat_struct *out,const dq_mat_struct *a,const dq_mat_struct *b)
{
    for(slong i=0;i<out->r*out->c;i++)field_sub(out->entries+i,a->entries+i,b->entries+i,DQ_FIELD);
}
/* Use native scalar multiplication tables for byte fields; all larger binary
 * fields retain their existing native arithmetic via field_mul. */
static void dq_row_addmul(dq_coeff *out,const dq_coeff *in,slong count,
                          const dq_coeff *scalar,dq_coeff *temp)
{
    if(!dq_nonzero(*scalar))return;
    if(extension_field.field_id==FIELD_ID_GF28) {
        const uint8_t *row=gf28_get_scalar_row(scalar->gf28);
        for(slong j=0;j<count;j++)out[j].gf28^=row[in[j].gf28];
    } else if(extension_field.field_id==FIELD_ID_GF24) {
        const uint8_t *row=gf24_get_scalar_row(scalar->gf24);
        for(slong j=0;j<count;j++)out[j].gf24^=row[in[j].gf24];
    } else {
        for(slong j=0;j<count;j++) {
            field_mul(temp,scalar,in+j,DQ_FIELD);
            field_add(out+j,out+j,temp,DQ_FIELD);
        }
    }
}
static void dq_mat_mul(dq_mat_struct *out,const dq_mat_struct *a,const dq_mat_struct *b)
{
    assert(out!=a && out!=b && a->c==b->r && out->r==a->r && out->c==b->c);
    if(extension_field.field_id==FIELD_ID_GF24 || extension_field.field_id==FIELD_ID_GF28) {
        int gf16=extension_field.field_id==FIELD_ID_GF24;
        uint8_t *left=flint_malloc(a->r*a->c),*right=flint_malloc(b->r*b->c);
        uint8_t *product=flint_calloc(out->r*out->c,1);
        for(slong i=0;i<a->r*a->c;i++)left[i]=gf16?a->entries[i].gf24:a->entries[i].gf28;
        for(slong i=0;i<b->r*b->c;i++)right[i]=gf16?b->entries[i].gf24:b->entries[i].gf28;
        for(slong i=0;i<a->r;i++)for(slong k=0;k<a->c;k++) {
            uint8_t c=left[i*a->c+k];
            if(c)quotient_byte_addmul(product+i*out->c,right+k*b->c,out->c,extension_byte_tables[c]);
        }
        for(slong i=0;i<out->r*out->c;i++) {
            if(gf16)out->entries[i].gf24=product[i];else out->entries[i].gf28=product[i];
        }
        flint_free(left);flint_free(right);flint_free(product);return;
    }
    dq_mat_zero(out);dq_coeff temp;dq_init(temp);
    for(slong i=0;i<a->r;i++)for(slong k=0;k<a->c;k++)
        dq_row_addmul(out->entries+i*out->c,b->entries+k*b->c,b->c,
                      &dq_mat_entry(a,i,k),&temp);
    dq_clear(temp);
}
static slong dq_mat_rref(dq_mat_struct *m)
{
    slong rank=0;dq_coeff inv,scalar,temp;
    dq_init(inv);dq_init(scalar);dq_init(temp);
    for(slong col=0;col<m->c && rank<m->r;col++) {
        slong pivot=rank;
        while(pivot<m->r && !dq_nonzero(dq_mat_entry(m,pivot,col)))pivot++;
        if(pivot==m->r)continue;
        if(pivot!=rank)for(slong j=0;j<m->c;j++) {
            dq_coeff swap=dq_mat_entry(m,pivot,j);
            dq_mat_entry(m,pivot,j)=dq_mat_entry(m,rank,j);
            dq_mat_entry(m,rank,j)=swap;
        }
        field_inv(&inv,&dq_mat_entry(m,rank,col),DQ_FIELD);
        for(slong j=col;j<m->c;j++)
            field_mul(&dq_mat_entry(m,rank,j),&dq_mat_entry(m,rank,j),&inv,DQ_FIELD);
        for(slong i=0;i<m->r;i++)if(i!=rank) {
            dq_neg(scalar,dq_mat_entry(m,i,col));
            dq_row_addmul(&dq_mat_entry(m,i,col),&dq_mat_entry(m,rank,col),m->c-col,&scalar,&temp);
        }
        rank++;
    }
    dq_clear(inv);dq_clear(scalar);dq_clear(temp);return rank;
}

/* Shared leading-frame MQ operations. LU uses the same packed-L layout as
 * FLINT: multipliers occupy the first rank columns, not the pivot columns. */
#define dq_mat_entry_ptr(m,i,j) (&dq_mat_entry(m,i,j))
static void dq_vec_set(dq_coeff *out,const dq_coeff *in,slong count)
{for(slong i=0;i<count;i++)dq_set(out[i],in[i]);}
static void dq_vec_addmul_impl(dq_coeff *out,const dq_coeff *in,slong count,const dq_coeff *scalar)
{dq_coeff temp;dq_init(temp);dq_row_addmul(out,in,count,scalar,&temp);dq_clear(temp);}
#define dq_vec_addmul(a,b,n,c) dq_vec_addmul_impl((a),(b),(n),&(c))
static void dq_mat_one(dq_mat_struct *m)
{dq_mat_zero(m);for(slong i=0;i<FLINT_MIN(m->r,m->c);i++)dq_one(dq_mat_entry(m,i,i));}
static void dq_mat_neg(dq_mat_struct *out,const dq_mat_struct *in)
{for(slong i=0;i<out->r*out->c;i++)dq_neg(out->entries[i],in->entries[i]);}
static int dq_mat_equal(const dq_mat_struct *a,const dq_mat_struct *b)
{
    if(a->r!=b->r || a->c!=b->c)return 0;
    for(slong i=0;i<a->r*a->c;i++)if(!field_equal(a->entries+i,b->entries+i,DQ_FIELD))return 0;
    return 1;
}
static uint8_t *dq_byte_copy(const dq_mat_struct *m)
{
    uint8_t *out=flint_malloc(m->r*m->c);
    for(slong i=0;i<m->r*m->c;i++)out[i]=extension_field.field_id==FIELD_ID_GF24?m->entries[i].gf24:m->entries[i].gf28;
    return out;
}
static void dq_byte_store(dq_mat_struct *m,const uint8_t *data)
{
    for(slong i=0;i<m->r*m->c;i++) {
        if(extension_field.field_id==FIELD_ID_GF24)m->entries[i].gf24=data[i];else m->entries[i].gf28=data[i];
    }
}
static const uint8_t *dq_byte_inverse_row(uint8_t pivot)
{
    return extension_field.field_id==FIELD_ID_GF24 ?
        gf24_get_scalar_row(gf24_inv(pivot)):gf28_get_scalar_row(gf28_inv(pivot));
}
static slong dq_byte_lu(slong *perm,dq_mat_struct *m)
{
    uint8_t *data=dq_byte_copy(m);slong rank=0;
    for(slong i=0;i<m->r;i++)perm[i]=i;
    for(slong col=0;col<m->c && rank<m->r;col++) {
        slong pivot=rank;while(pivot<m->r && !data[pivot*m->c+col])pivot++;
        if(pivot==m->r)continue;
        if(pivot!=rank) {
            slong p=perm[pivot];perm[pivot]=perm[rank];perm[rank]=p;
            for(slong j=0;j<m->c;j++) {
                uint8_t t=data[pivot*m->c+j];data[pivot*m->c+j]=data[rank*m->c+j];data[rank*m->c+j]=t;
            }
        }
        const uint8_t *inverse=dq_byte_inverse_row(data[rank*m->c+col]);
        for(slong i=rank+1;i<m->r;i++) {
            uint8_t factor=inverse[data[i*m->c+col]];
            if(factor)quotient_byte_addmul(data+i*m->c+col+1,data+rank*m->c+col+1,m->c-col-1,extension_byte_tables[factor]);
            data[i*m->c+col]=0;data[i*m->c+rank]=factor;
        }
        rank++;
    }
    dq_byte_store(m,data);flint_free(data);return rank;
}
static void dq_byte_solve_tri(dq_mat_struct *out,const dq_mat_struct *a,
                               const dq_mat_struct *b,int unit,int upper)
{
    if(!out->c)return;
    uint8_t *matrix=dq_byte_copy(a),*rhs=dq_byte_copy(b);
    for(slong step=0;step<a->r;step++) {
        slong i=upper?a->r-1-step:step;
        for(slong jj=0;jj<step;jj++) {
            slong j=upper?a->r-1-jj:jj;uint8_t scalar=matrix[i*a->c+j];
            if(scalar)quotient_byte_addmul(rhs+i*b->c,rhs+j*b->c,b->c,extension_byte_tables[scalar]);
        }
        if(!unit) {
            const uint8_t *inverse=dq_byte_inverse_row(matrix[i*a->c+i]);
            for(slong j=0;j<b->c;j++)rhs[i*b->c+j]=inverse[rhs[i*b->c+j]];
        }
    }
    dq_byte_store(out,rhs);flint_free(matrix);flint_free(rhs);
}
static slong dq_mat_lu(slong *perm,dq_mat_struct *m,int unused)
{
    if(extension_field.field_id==FIELD_ID_GF24 || extension_field.field_id==FIELD_ID_GF28)
        return dq_byte_lu(perm,m);
    (void)unused;slong rank=0;dq_coeff inv,factor,negative,temp;
    dq_init(inv);dq_init(factor);dq_init(negative);dq_init(temp);
    for(slong i=0;i<m->r;i++)perm[i]=i;
    for(slong col=0;col<m->c && rank<m->r;col++) {
        slong pivot=rank;while(pivot<m->r && !dq_nonzero(dq_mat_entry(m,pivot,col)))pivot++;
        if(pivot==m->r)continue;
        if(pivot!=rank) {
            slong p=perm[pivot];perm[pivot]=perm[rank];perm[rank]=p;
            for(slong j=0;j<m->c;j++) {
                dq_coeff swap=dq_mat_entry(m,pivot,j);
                dq_mat_entry(m,pivot,j)=dq_mat_entry(m,rank,j);dq_mat_entry(m,rank,j)=swap;
            }
        }
        field_inv(&inv,&dq_mat_entry(m,rank,col),DQ_FIELD);
        for(slong i=rank+1;i<m->r;i++) {
            field_mul(&factor,&dq_mat_entry(m,i,col),&inv,DQ_FIELD);dq_neg(negative,factor);
            if(col+1<m->c)dq_row_addmul(&dq_mat_entry(m,i,col+1),&dq_mat_entry(m,rank,col+1),m->c-col-1,&negative,&temp);
            field_set_zero(&dq_mat_entry(m,i,col),DQ_FIELD);
            dq_set(dq_mat_entry(m,i,rank),factor);
        }
        rank++;
    }
    dq_clear(inv);dq_clear(factor);dq_clear(negative);dq_clear(temp);return rank;
}
static void dq_mat_solve_tri(dq_mat_struct *out,const dq_mat_struct *a,
                             const dq_mat_struct *b,int unit,int upper)
{
    if(extension_field.field_id==FIELD_ID_GF24 || extension_field.field_id==FIELD_ID_GF28) {
        dq_byte_solve_tri(out,a,b,unit,upper);return;
    }
    dq_coeff scalar,inv,temp;dq_init(scalar);dq_init(inv);dq_init(temp);
    for(slong step=0;step<a->r;step++) {
        slong i=upper?a->r-1-step:step;
        if(out!=b)dq_vec_set(&dq_mat_entry(out,i,0),&dq_mat_entry(b,i,0),b->c);
        for(slong jj=0;jj<step;jj++) {
            slong j=upper?a->r-1-jj:jj;dq_neg(scalar,dq_mat_entry(a,i,j));
            dq_row_addmul(&dq_mat_entry(out,i,0),&dq_mat_entry(out,j,0),out->c,&scalar,&temp);
        }
        if(!unit) {
            field_inv(&inv,&dq_mat_entry(a,i,i),DQ_FIELD);
            for(slong j=0;j<out->c;j++)field_mul(&dq_mat_entry(out,i,j),&dq_mat_entry(out,i,j),&inv,DQ_FIELD);
        }
    }
    dq_clear(scalar);dq_clear(inv);dq_clear(temp);
}
#define dq_mat_solve_tril(o,a,b,u) dq_mat_solve_tri(o,a,b,u,0)
#define dq_mat_solve_triu(o,a,b,u) dq_mat_solve_tri(o,a,b,u,1)
static void dq_mat_accumulate(dq_mat_struct *out,const dq_mat_struct *base,
                              const dq_mat_struct *a,const dq_mat_struct *b,int negative)
{
    if(out!=base)dq_mat_set(out,base);
    dq_coeff scalar,temp;dq_init(scalar);dq_init(temp);
    for(slong i=0;i<a->r;i++)for(slong k=0;k<a->c;k++) {
        if(negative)dq_neg(scalar,dq_mat_entry(a,i,k));else dq_set(scalar,dq_mat_entry(a,i,k));
        dq_row_addmul(&dq_mat_entry(out,i,0),&dq_mat_entry(b,k,0),out->c,&scalar,&temp);
    }
    dq_clear(scalar);dq_clear(temp);
}
#define dq_mat_addmul(o,c,a,b) dq_mat_accumulate(o,c,a,b,0)
#define dq_mat_submul(o,c,a,b) dq_mat_accumulate(o,c,a,b,1)

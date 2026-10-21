/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Native byte analogue of the prime-field packed MQ kernel. Each subset
 * minor stays in compact native coordinates. Characteristic two removes
 * Laplace signs; nibble shuffles implement exact GF(16)/GF(256) products. */
#include "solver/quotient_byte.h"
typedef struct {slong r,c;uint8_t *data;} mq_byte_matrix;
#define mq_byte_tables extension_byte_tables
static int mq_byte_enabled(void)
{return extension_field.field_id==FIELD_ID_GF24 || extension_field.field_id==FIELD_ID_GF28;}
static mq_byte_matrix mq_byte_pack(const dq_mat_struct *m)
{
    mq_byte_matrix out={m->r,m->c,NULL};
    if(m->r) {
        out.data=flint_malloc(m->r*m->c);
        for(slong i=0;i<m->r*m->c;i++)out.data[i]=extension_field.field_id==FIELD_ID_GF24 ? m->entries[i].gf24:m->entries[i].gf28;
    }
    return out;
}
static void mq_byte_product(mq_byte_matrix *out,const mq_byte_matrix *a,const mq_byte_matrix *b)
{
    for(slong i=0;i<a->r;i++)for(slong k=0;k<a->c;k++) {
        uint8_t c=a->data[i*a->c+k];
        if(c)quotient_byte_addmul(out->data+i*out->c,b->data+k*b->c,out->c,mq_byte_tables[c]);
    }
}
static void mq_byte_det(dq_mat_struct *out,const grade *g,const dq_coeff *q)
{
    slong masks=1L<<N,stride=M+1;
    mq_byte_matrix *state=flint_calloc(masks*stride,sizeof(*state));
    state[0]=(mq_byte_matrix){1,1,flint_malloc(1)};state[0].data[0]=1;
    for(slong k=1;k<=N;k++) {
        /* Only this row's maps are needed by the current subset layer. */
        linear_map *maps=make_linear_maps_for_row(g,q,k-1);
        mq_byte_matrix *left=flint_calloc(N*L,sizeof(*left)),*right=flint_calloc(N*L,sizeof(*right));
        for(slong j=0;j<N;j++)for(slong d=0;d<L;d++) {
            const linear_map *m=maps+((k-1)*N+j)*L+d;
            left[j*L+d]=mq_byte_pack(m->x);right[j*L+d]=mq_byte_pack(m->y);
        }
        clear_linear_maps(maps);
        for(slong mask=1;mask<masks;mask++)if(__builtin_popcountl(mask)==k) {
            for(slong a=FLINT_MAX(0,k-K);a<=FLINT_MIN(M,k);a++) {
                slong b=k-a;mq_byte_matrix *dest=state+mask*stride+a;
                *dest=(mq_byte_matrix){g[a].h,g[b].h,flint_calloc(g[a].h*g[b].h,1)};
                for(slong j=0;j<N;j++)if(mask&(1L<<j)) {
                    slong child=mask^(1L<<j);
                    if(a)mq_byte_product(dest,left+j*L+a-1,state+child*stride+a-1);
                    if(b)mq_byte_product(dest,state+child*stride+a,right+j*L+b-1);
                }
            }
        }
        for(slong i=0;i<N*L;i++){flint_free(left[i].data);flint_free(right[i].data);}
        flint_free(left);flint_free(right);
        for(slong mask=0;mask<masks;mask++)if(__builtin_popcountl(mask)==k-1)
            for(slong a=0;a<stride;a++){flint_free(state[mask*stride+a].data);state[mask*stride+a].data=NULL;}
    }
    const mq_byte_matrix *last=state+(masks-1)*stride+M;
    for(slong i=0;i<out->r*out->c;i++) {
        if(extension_field.field_id==FIELD_ID_GF24)out->entries[i].gf24=last->data[i];
        else out->entries[i].gf28=last->data[i];
    }
    flint_free(last->data);flint_free(state);
}

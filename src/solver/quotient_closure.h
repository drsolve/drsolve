/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Shared certified Macaulay closure for unified field backends. Native input
 * access and rational-point extraction are supplied by the field adapters. */
#ifndef DRSOLVE_QUOTIENT_CLOSURE_H
#define DRSOLVE_QUOTIENT_CLOSURE_H
#include "quotient_basis.h"
#include "fq_mat_det.h"
#include "quotient_byte.h"

typedef struct {
    field_ctx_t field;
    const void *arithmetic;
    const fq_nmod_ctx_struct *fq;
    uint8_t byte_tables[256][32];
} quotient_field_t;

#define FIELD(e) (e)->field.field_id, (e)->arithmetic
#define ENTRY(M,i,j) ((M)->entries+(i)*(M)->c+(j))

typedef struct quotient_input quotient_input_t;
struct quotient_input {
    const void *polys;
    const void *ctx;
    slong (*length)(const quotient_input_t *, slong);
    void (*term)(field_elem_u *, ulong *, const quotient_input_t *, slong, slong,
                 const quotient_field_t *);
};

static field_elem_u *vector_init(slong n, const quotient_field_t *e)
{
    field_elem_u *v=flint_calloc(n,sizeof(*v));
    for (slong i=0;i<n;i++) {
        field_init_elem(v+i,FIELD(e));
        field_set_zero(v+i,FIELD(e)); /* Zech zero is not integer zero. */
    }
    return v;
}

static void vector_clear(field_elem_u *v, slong n, const quotient_field_t *e)
{
    for (slong i=0;i<n;i++) field_clear_elem(v+i,FIELD(e));
    flint_free(v);
}

static void vector_zero(field_elem_u *v, slong n, const quotient_field_t *e)
{
    for (slong i=0;i<n;i++) field_set_zero(v+i,FIELD(e));
}

static void vector_set(field_elem_u *v, const field_elem_u *w, slong n,
                       const quotient_field_t *e)
{
    for (slong i=0;i<n;i++) field_set_elem(v+i,w+i,FIELD(e));
}

/* Hot row kernel: retain native binary arithmetic, including the existing
 * scalar multiplication tables for GF(16) and GF(256). */
static void vector_addmul(field_elem_u *v, const field_elem_u *w, slong n,
                         const field_elem_u *a, const quotient_field_t *e)
{
    if (field_is_zero(a,FIELD(e))) return;
    if (e->field.field_id==FIELD_ID_GF28) {
        const uint8_t *row=gf28_get_scalar_row(a->gf28);
        for (slong j=0;j<n;j++) v[j].gf28^=row[w[j].gf28];
    } else if (e->field.field_id==FIELD_ID_GF24) {
        const uint8_t *row=gf24_get_scalar_row(a->gf24);
        for (slong j=0;j<n;j++) v[j].gf24^=row[w[j].gf24];
    } else {
        field_elem_u t; field_init_elem(&t,FIELD(e));
        for (slong j=0;j<n;j++) {
            field_mul(&t,a,w+j,FIELD(e));
            field_add(v+j,v+j,&t,FIELD(e));
        }
        field_clear_elem(&t,FIELD(e));
    }
}

static void matrix_new(unified_mat_t M, slong r, slong c, quotient_field_t *e)
{
    M->r=r; M->c=c; M->ctx=&e->field;
    M->entries=vector_init(r*c,e);
}

static void matrix_product(unified_mat_t C, const unified_mat_t A,
                           const unified_mat_t B, const quotient_field_t *e)
{
    assert(A->c==B->r && C->r==A->r && C->c==B->c);
    vector_zero(C->entries,C->r*C->c,e);
    for (slong i=0;i<A->r;i++) for (slong k=0;k<A->c;k++)
        vector_addmul(ENTRY(C,i,0),ENTRY(B,k,0),C->c,ENTRY(A,i,k),e);
}

/* Sparse-in-rows echelon storage. A pivot row never changes after insertion,
 * so each low-degree relation needs to be multiplied by each variable once. */
typedef struct {
    basis_t basis;
    field_elem_u **rows;
    /* Packed native tails for the elimination kernel. The unified rows are
     * materialized only when a new pivot is stored, for the shared closure
     * and algebra certificate. */
    void **native_rows;
    void *native_scratch;
    unsigned char *processed;
    slong rank;
} quotient_closure_t;

static void quotient_closure_init(quotient_closure_t *cl, slong n, slong d)
{
    basis_init(&cl->basis,n,d);
    cl->rows=flint_calloc(cl->basis.count,sizeof(*cl->rows)); cl->rank=0;
    cl->native_rows=flint_calloc(cl->basis.count,sizeof(*cl->native_rows));
    cl->native_scratch=NULL;
    cl->processed=flint_calloc(cl->basis.count,1);
}

static void quotient_closure_clear(quotient_closure_t *cl, const quotient_field_t *e)
{
    for (slong j=0;j<cl->basis.count;j++) {
        if (cl->rows[j]) vector_clear(cl->rows[j],cl->basis.count,e);
        flint_free(cl->native_rows[j]);
    }
    flint_free(cl->native_rows); flint_free(cl->native_scratch);
    flint_free(cl->processed);
    flint_free(cl->rows); basis_clear(&cl->basis);
}

#include "quotient_gf2n.h"

static void insert_relation(quotient_closure_t *cl, field_elem_u *row,
                            const quotient_field_t *e)
{
    switch (e->field.field_id) {
    case FIELD_ID_GF24: quotient_insert_gf24(cl,row,e); return;
    case FIELD_ID_GF28: quotient_insert_gf28(cl,row,e); return;
    case FIELD_ID_GF216: quotient_insert_gf216(cl,row,e); return;
    case FIELD_ID_GF232: quotient_insert_gf232(cl,row,e); return;
    case FIELD_ID_GF264: quotient_insert_gf264(cl,row,e); return;
    case FIELD_ID_GF2128: quotient_insert_gf2128(cl,row,e); return;
    case FIELD_ID_GF2N8: quotient_insert_gf2n8(cl,row,e); return;
    case FIELD_ID_GF2N16: quotient_insert_gf2n16(cl,row,e); return;
    case FIELD_ID_GF2N64: quotient_insert_gf2n64(cl,row,e); return;
    case FIELD_ID_GF2N128: quotient_insert_gf2n128(cl,row,e); return;
    default: break;
    }
    slong N=cl->basis.count;
    field_elem_u a; field_init_elem(&a,FIELD(e));
    for (slong j=0;j<N;j++) if (!field_is_zero(row+j,FIELD(e))) {
        if (cl->rows[j]) {
            field_neg(&a,row+j,FIELD(e));
            vector_addmul(row+j,cl->rows[j]+j,N-j,&a,e);
        } else {
            field_inv(&a,row+j,FIELD(e));
            cl->rows[j]=vector_init(N,e);
            for (slong k=j;k<N;k++) field_mul(cl->rows[j]+k,row+k,&a,FIELD(e));
            cl->rank++; break;
        }
    }
    field_clear_elem(&a,FIELD(e));
}

static void extend_relations(quotient_closure_t *cl, const quotient_field_t *e)
{
    slong N=cl->basis.count,n=cl->basis.n;
    unsigned char *processed=cl->processed;
    field_elem_u *row=vector_init(N,e);
    int more=1;
    while (more && cl->rank<N) {
        more=0;
        for (slong j=0;j<N;j++) if (cl->rows[j] && !processed[j] &&
                cl->basis.degrees[j]<cl->basis.degree) {
            processed[j]=1; more=1;
            for (slong v=0;v<n && cl->rank<N;v++) {
                vector_zero(row,N,e);
                for (slong k=j;k<N;k++) if (!field_is_zero(cl->rows[j]+k,FIELD(e))) {
                    ulong exp[n]; memcpy(exp,cl->basis.exps+k*n,sizeof(exp)); exp[v]++;
                    slong dest=basis_find(&cl->basis,exp); assert(dest>=0);
                    field_set_elem(row+dest,cl->rows[j]+k,FIELD(e));
                }
                insert_relation(cl,row,e);
            }
        }
    }
    vector_clear(row,N,e);
}

static void insert_generators(quotient_closure_t *cl, const quotient_input_t *input,
                              slong m, const slong *degrees,
                              const quotient_field_t *e)
{
    slong N=cl->basis.count,n=cl->basis.n;
    field_elem_u *row=vector_init(N,e);
    for (slong i=0;i<m && cl->rank<N;i++) if (degrees[i]>=0) {
        basis_t mult; basis_init(&mult,n,cl->basis.degree-degrees[i]);
        field_elem_u *cs=vector_init(input->length(input,i),e);
        ulong *exps=flint_malloc(input->length(input,i)*n*sizeof(ulong));
        for (slong t=0;t<input->length(input,i);t++) {
            input->term(cs+t,exps+t*n,input,i,t,e);
        }
        for (slong j=0;j<mult.count && cl->rank<N;j++) {
            vector_zero(row,N,e);
            for (slong t=0;t<input->length(input,i);t++) {
                ulong exp[n];
                for (slong v=0;v<n;v++) exp[v]=exps[t*n+v]+mult.exps[j*n+v];
                slong col=basis_find(&cl->basis,exp); assert(col>=0);
                field_set_elem(row+col,cs+t,FIELD(e));
            }
            insert_relation(cl,row,e);
        }
        flint_free(exps); vector_clear(cs,input->length(input,i),e); basis_clear(&mult);
    }
    vector_clear(row,N,e);
}

/* Reuse the preceding degree instead of rebuilding all Macaulay rows.
 * Its relations and their variable multiples span every original multiple
 * needed at the new degree. Distinct leading monomials can be installed
 * directly: multiplication preserves the degree-lex order. Collisions must
 * still be reduced, since their differences can yield new relations. */
static void quotient_closure_seed(quotient_closure_t *cl,
                         const quotient_closure_t *prev, const quotient_field_t *e)
{
    slong N=cl->basis.count,P=prev->basis.count,n=cl->basis.n,offset=N-P;
    slong *shift=flint_malloc(P*n*sizeof(*shift));
    for (slong j=0;j<P;j++) for (slong v=0;v<n;v++) {
        ulong exp[n]; memcpy(exp,prev->basis.exps+j*n,sizeof(exp)); exp[v]++;
        shift[j*n+v]=basis_find(&cl->basis,exp);
        assert(shift[j*n+v]>=0);
    }
    slong *owner=flint_malloc(N*sizeof(*owner));
    for (slong j=0;j<N;j++) owner[j]=-1;
    for (slong j=0;j<P;j++) if (prev->rows[j]) {
        cl->rows[j+offset]=vector_init(N,e);
        vector_set(cl->rows[j+offset]+offset,prev->rows[j],P,e);
        /* All variable multiples of this unchanged row are inserted below. */
        cl->processed[j+offset]=1;
        cl->rank++;
    }
    /* First install one row per new pivot, so duplicate multiples reduce
     * against the entire known triangular block. */
    for (slong j=0;j<P;j++) if (prev->rows[j]) for (slong v=0;v<n;v++) {
        slong pivot=shift[j*n+v];
        if (cl->rows[pivot]) continue;
        field_elem_u *row=cl->rows[pivot]=vector_init(N,e);
        for (slong k=j;k<P;k++)
            field_set_elem(row+shift[k*n+v],prev->rows[j]+k,FIELD(e));
        owner[pivot]=j*n+v; cl->rank++;
    }
    field_elem_u *row=vector_init(N,e);
    for (slong j=0;j<P && cl->rank<N;j++) if (prev->rows[j])
        for (slong v=0;v<n && cl->rank<N;v++) {
            if (owner[shift[j*n+v]]==j*n+v) continue;
            vector_zero(row,N,e);
            for (slong k=j;k<P;k++)
                field_set_elem(row+shift[k*n+v],prev->rows[j]+k,FIELD(e));
            insert_relation(cl,row,e);
        }
    vector_clear(row,N,e); flint_free(owner); flint_free(shift);
}

static void quotient_act_monomial(field_elem_u *out, const field_elem_u *one,
                                   const ulong *exp, unified_mat_struct *T,
                                   slong n, slong dim, const quotient_field_t *e)
{
    field_elem_u *tmp=vector_init(dim,e); vector_set(out,one,dim,e);
    for (slong v=0;v<n;v++) for (ulong k=0;k<exp[v];k++) {
        vector_zero(tmp,dim,e);
        for (slong i=0;i<dim;i++) vector_addmul(tmp,ENTRY(T+v,i,0),dim,out+i,e);
        vector_set(out,tmp,dim,e);
    }
    vector_clear(tmp,dim,e);
}

/* Both directions of the certificate matter: input vanishing and commuting
 * alone only prove a quotient of the desired algebra. The free monomials
 * must also act on 1 as the full coordinate basis. */
static int quotient_certify(unified_mat_struct *T, const field_elem_u *one,
                             const basis_t *basis, const slong *freecols,
                             const quotient_input_t *input, slong m,
                             quotient_field_t *e)
{
    slong n=basis->n,dim=T[0].r; int ok=1;
    unified_mat_t AB,BA; matrix_new(AB,dim,dim,e); matrix_new(BA,dim,dim,e);
    for (slong i=0;i<n && ok;i++) for (slong j=0;j<i && ok;j++) {
        matrix_product(AB,T+i,T+j,e); matrix_product(BA,T+j,T+i,e);
        for (slong k=0;k<dim*dim;k++) if (!field_equal(AB->entries+k,BA->entries+k,FIELD(e))) { ok=0; break; }
    }
    field_elem_u *v=vector_init(dim,e),*sum=vector_init(dim,e),a;
    field_init_elem(&a,FIELD(e));
    for (slong i=0;i<dim && ok;i++) {
        quotient_act_monomial(v,one,basis->exps+freecols[i]*n,T,n,dim,e);
        for (slong j=0;j<dim;j++)
            if (i==j ? !field_is_one(v+j,FIELD(e)) : !field_is_zero(v+j,FIELD(e))) ok=0;
    }
    for (slong i=0;i<m && ok;i++) {
        vector_zero(sum,dim,e);
        for (slong j=0;j<input->length(input,i);j++) {
            ulong exp[n]; input->term(&a,exp,input,i,j,e);
            quotient_act_monomial(v,one,exp,T,n,dim,e);
            vector_addmul(sum,v,dim,&a,e);
        }
        for (slong j=0;j<dim;j++) if (!field_is_zero(sum+j,FIELD(e))) ok=0;
    }
    field_clear_elem(&a,FIELD(e));
    vector_clear(v,dim,e); vector_clear(sum,dim,e);
    unified_mat_clear(AB); unified_mat_clear(BA); return ok;
}

static unified_mat_struct *build_operators(const quotient_closure_t *cl,
                  const quotient_input_t *input, slong m,
                  quotient_field_t *e)
{
    slong N=cl->basis.count,n=cl->basis.n,dim=N-cl->rank,used=0;
    slong *freecols=flint_malloc(dim*sizeof(slong));
    for (slong j=0;j<N;j++) if (!cl->rows[j]) {
        if (cl->basis.degrees[j]>=cl->basis.degree) { flint_free(freecols); return NULL; }
        freecols[used++]=j;
    }
    unified_mat_t W; matrix_new(W,N,dim,e);
    for (slong j=0;j<dim;j++) field_set_one(ENTRY(W,freecols[j],j),FIELD(e));
    field_elem_u a; field_init_elem(&a,FIELD(e));
    for (slong j=N;j-- > 0;) if (cl->rows[j]) {
        for (slong k=j+1;k<N;k++) if (!field_is_zero(cl->rows[j]+k,FIELD(e))) {
            field_neg(&a,cl->rows[j]+k,FIELD(e));
            vector_addmul(ENTRY(W,j,0),ENTRY(W,k,0),dim,&a,e);
        }
    }
    field_clear_elem(&a,FIELD(e));
    unified_mat_struct *T=flint_malloc(n*sizeof(*T));
    for (slong v=0;v<n;v++) {
        matrix_new(T+v,dim,dim,e);
        for (slong j=0;j<dim;j++) {
            ulong exp[n]; memcpy(exp,cl->basis.exps+freecols[j]*n,sizeof(exp)); exp[v]++;
            slong dest=basis_find(&cl->basis,exp); assert(dest>=0);
            vector_set(ENTRY(T+v,j,0),ENTRY(W,dest,0),dim,e);
        }
    }
    int ok=quotient_certify(T,ENTRY(W,N-1,0),&cl->basis,freecols,input,m,e);
    unified_mat_clear(W); flint_free(freecols);
    if (!ok) {
        for (slong v=0;v<n;v++) unified_mat_clear(T+v);
        flint_free(T); return NULL;
    }
    return T;
}


#endif

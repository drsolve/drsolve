/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Instantiated for each native binary field. Only pivot tails are stored;
 * scalar dispatch and union-sized strides stay outside the hot row loop. */
static void MQGF_NAME(quotient_closure_t *cl, const field_elem_u *input,
                       const quotient_field_t *e)
{
    slong N=cl->basis.count;
    if (!cl->native_scratch)
        cl->native_scratch=flint_malloc(N*sizeof(MQGF_TYPE));
    MQGF_TYPE *row=cl->native_scratch;
    for (slong j=0;j<N;j++) row[j]=input[j].MQGF_MEMBER;
    for (slong j=0;j<N;j++) if (!MQGF_ZERO(row[j])) {
        if (cl->rows[j]) {
            MQGF_TYPE *tail=cl->native_rows[j];
            if (!tail) {
                /* Inherited rows are packed on first use. */
                tail=flint_malloc((N-j)*sizeof(*tail));
                for (slong k=j;k<N;k++) tail[k-j]=cl->rows[j][k].MQGF_MEMBER;
                cl->native_rows[j]=tail;
            }
            MQGF_TYPE a=row[j];
#ifdef MQGF_BYTE
            quotient_byte_addmul(row+j+1,tail+1,N-j-1,e->byte_tables[a]);
#else
            MQGF_PREPARE(a);
            for (slong k=j+1;k<N;k++) {
                MQGF_TYPE product=MQGF_MUL(a,tail[k-j]);
                row[k]=MQGF_ADD(row[k],product);
            }
#endif
        } else {
            field_elem_u pivot,inv;
            pivot.MQGF_MEMBER=row[j];
            field_inv(&inv,&pivot,FIELD(e));
            MQGF_TYPE a=inv.MQGF_MEMBER;
            MQGF_PREPARE(a);
            MQGF_TYPE *tail=flint_malloc((N-j)*sizeof(*tail));
            cl->native_rows[j]=tail;
            cl->rows[j]=vector_init(N,e);
            for (slong k=j;k<N;k++) {
                tail[k-j]=MQGF_MUL(a,row[k]);
                cl->rows[j][k].MQGF_MEMBER=tail[k-j];
            }
            cl->rank++; break;
        }
    }
}

#undef MQGF_TYPE
#undef MQGF_NAME
#undef MQGF_MEMBER
#undef MQGF_ZERO
#undef MQGF_ADD
#undef MQGF_PREPARE
#undef MQGF_MUL
#undef MQGF_BYTE

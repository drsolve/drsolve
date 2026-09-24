/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Included by fq_mpoly_mat_det.c. The canonical ORD_LEX output already groups
 * equal original-variable supports, so no term-sized permutation is needed. */
void fq_mq_compact_clear(fq_mq_compact *p)
{
    flint_free(p->rows); flint_free(p->cols); flint_free(p->offset); flint_free(p->terms);
    memset(p,0,sizeof(*p));
}

void fq_mq_compact_materialize(fq_mvpoly_t *out, const fq_mq_compact *p,
                              const fq_nmod_ctx_t ctx)
{
    fq_mvpoly_init(out,2*p->nvars,1,ctx);
    slong *exp=flint_malloc(2*p->nvars*sizeof(slong));
    fq_nmod_t coefficient; fq_nmod_init(coefficient,ctx);
    for(slong r=0;r<p->nrows;r++) {
        memcpy(exp,p->rows+r*p->nvars,p->nvars*sizeof(slong));
        for(slong t=p->offset[r];t<p->offset[r+1];t++) {
            const fq_mq_compact_term *a=p->terms+t;
            memcpy(exp+p->nvars,p->cols+a->column*p->nvars,p->nvars*sizeof(slong));
            fq_nmod_set_ui(coefficient,a->coefficient,ctx);
            fq_mvpoly_add_term_fast(out,exp,&a->degree,coefficient);
        }
    }
    fq_nmod_clear(coefficient,ctx); flint_free(exp);
}

static void mq_compact_from_nmod(fq_mq_compact *out, const nmod_mpoly_t poly,
    const nmod_mpoly_ctx_t ctx, slong n, slong row_capacity, slong col_capacity)
{
    fq_mq_compact p={0};
    p.nvars=n; p.prime=ctx->mod.n; p.nterms=poly->length;
    p.rows=flint_malloc((size_t)row_capacity*n*sizeof(slong));
    p.cols=flint_malloc((size_t)col_capacity*n*sizeof(slong));
    p.offset=flint_calloc(row_capacity+1,sizeof(slong));
    p.terms=flint_malloc((size_t)p.nterms*sizeof(*p.terms));
    slong hash_size=16;
    while(hash_size<2*col_capacity) hash_size*=2;
    slong *table=flint_calloc(hash_size,sizeof(slong));
    ulong *exp=flint_malloc((size_t)(2*n+1)*sizeof(ulong));
    for(slong t=0;t<p.nterms;t++) {
        nmod_mpoly_get_term_exp_ui(exp,poly,t,ctx);
        int same=p.nrows>0;
        for(slong v=0;v<n && same;v++) same=p.rows[(p.nrows-1)*n+v]==(slong)exp[v];
        if(!same) {
            FLINT_ASSERT(p.nrows<row_capacity);
            p.offset[p.nrows]=t;
            for(slong v=0;v<n;v++) p.rows[p.nrows*n+v]=(slong)exp[v];
            p.nrows++;
        }
        /* Hash full exponent vectors; no single-word axis packing limit. */
        ulong hash=0;
        for(slong v=0;v<n;v++) hash=hash*UWORD(1000003)+exp[n+v]+1;
        slong bucket=mq_monom_hash(hash)&(hash_size-1), c;
        for(;;) {
            c=table[bucket]-1;
            if(c<0) {
                FLINT_ASSERT(p.ncols<col_capacity);
                c=p.ncols++; table[bucket]=c+1;
                for(slong v=0;v<n;v++) p.cols[c*n+v]=(slong)exp[n+v];
                break;
            }
            same=1;
            for(slong v=0;v<n && same;v++) same=p.cols[c*n+v]==(slong)exp[n+v];
            if(same) break;
            bucket=(bucket+1)&(hash_size-1);
        }
        p.terms[t]=(fq_mq_compact_term){c,(slong)exp[2*n],poly->coeffs[t]};
        p.degree=FLINT_MAX(p.degree,(slong)exp[2*n]);
    }
    p.offset[p.nrows]=p.nterms;
    flint_free(exp); flint_free(table);
    *out=p;
}

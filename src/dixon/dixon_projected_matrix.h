/* SPDX-License-Identifier: GPL-2.0-or-later */
/* A verified Step 1 projection already is the selected square block. Extract
 * content and degree permutations from terms, then fill the final univariate
 * matrix directly. The input polynomial remains owned by the caller. */
#ifdef DRSOLVE_MQ_STEP2_TEST
static int dixon_step2_test_force_generic;
static int dixon_step2_test_direct_calls;
#endif

static slong dixon_projected_poly_matrix(fq_nmod_poly_mat_t out,
    slong *rows, slong *cols, slong size, const fq_mvpoly_t *poly,
    const slong *term_rows, const slong *term_cols)
{
    slong *rp=flint_malloc(size*sizeof(slong)), *cp=flint_malloc(size*sizeof(slong));
    fq_index_degree_pair *rd=flint_malloc(size*sizeof(*rd)), *cd=flint_malloc(size*sizeof(*cd));
    for(slong i=0;i<size;i++) {
        rp[i]=cp[i]=WORD_MAX;
        rd[i]=(fq_index_degree_pair){i,-1}; cd[i]=(fq_index_degree_pair){i,-1};
    }
    /* Preserve the legacy row-first, then column valuation extraction. */
    for(slong t=0;t<poly->nterms;t++) if(term_rows[t]>=0) {
        slong d=poly->terms[t].par_exp ? poly->terms[t].par_exp[0] : 0;
        rp[term_rows[t]]=FLINT_MIN(rp[term_rows[t]],d);
    }
    for(slong i=0;i<size;i++) if(rp[i]==WORD_MAX) rp[i]=0;
    for(slong t=0;t<poly->nterms;t++) if(term_rows[t]>=0) {
        slong d=(poly->terms[t].par_exp ? poly->terms[t].par_exp[0] : 0)-rp[term_rows[t]];
        cp[term_cols[t]]=FLINT_MIN(cp[term_cols[t]],d);
    }
    slong content=0;
    for(slong i=0;i<size;i++) { if(cp[i]==WORD_MAX) cp[i]=0; content+=rp[i]+cp[i]; }
    for(slong t=0;t<poly->nterms;t++) if(term_rows[t]>=0) {
        slong r=term_rows[t],c=term_cols[t];
        slong d=(poly->terms[t].par_exp ? poly->terms[t].par_exp[0] : 0)-rp[r]-cp[c];
        FLINT_ASSERT(d>=0);
        rd[r].degree=FLINT_MAX(rd[r].degree,d); cd[c].degree=FLINT_MAX(cd[c].degree,d);
    }
    const char *reorder=getenv("DRSOLVE_PREDICT_REORDER");
    if(!reorder || strcmp(reorder,"0")!=0) {
        qsort(rd,size,sizeof(*rd),compare_fq_degrees);
        qsort(cd,size,sizeof(*cd),compare_fq_degrees);
    }
    slong *rmap=flint_malloc(size*sizeof(slong)), *cmap=flint_malloc(size*sizeof(slong));
    for(slong i=0;i<size;i++) {
        rows[i]=rd[i].index; cols[i]=cd[i].index;
        rmap[rows[i]]=i; cmap[cols[i]]=i;
    }
    flint_free(rd); flint_free(cd);
    slong *offset=flint_calloc(size+1,sizeof(slong)), *write=flint_malloc(size*sizeof(slong));
    for(slong t=0;t<poly->nterms;t++) if(term_rows[t]>=0) offset[term_rows[t]+1]++;
    for(slong r=0;r<size;r++) offset[r+1]+=offset[r];
    memcpy(write,offset,size*sizeof(slong));
    slong *terms=flint_malloc(offset[size]*sizeof(slong));
    for(slong t=0;t<poly->nterms;t++) if(term_rows[t]>=0) terms[write[term_rows[t]]++]=t;
    flint_free(write);
    fq_nmod_poly_mat_init(out,size,size,poly->ctx);
#ifdef _OPENMP
    #pragma omp parallel if(size>1)
#endif
    {
        slong *maximum=flint_malloc(size*sizeof(slong));
#ifdef _OPENMP
        #pragma omp for schedule(dynamic,1)
#endif
        for(slong r=0;r<size;r++) {
            for(slong c=0;c<size;c++) maximum[c]=-1;
            for(slong i=offset[r];i<offset[r+1];i++) {
                slong t=terms[i],c=term_cols[t];
                slong d=(poly->terms[t].par_exp ? poly->terms[t].par_exp[0] : 0)-rp[r]-cp[c];
                maximum[c]=FLINT_MAX(maximum[c],d);
            }
            for(slong c=0;c<size;c++) if(maximum[c]>=0)
                fq_nmod_poly_fit_length(fq_nmod_poly_mat_entry(out,rmap[r],cmap[c]),maximum[c]+1,poly->ctx);
            for(slong i=offset[r];i<offset[r+1];i++) {
                slong t=terms[i],c=term_cols[t];
                slong d=(poly->terms[t].par_exp ? poly->terms[t].par_exp[0] : 0)-rp[r]-cp[c];
                fq_nmod_poly_set_coeff(fq_nmod_poly_mat_entry(out,rmap[r],cmap[c]),d,poly->terms[t].coeff,poly->ctx);
            }
        }
        flint_free(maximum);
    }
    flint_free(terms); flint_free(offset);
    flint_free(rmap); flint_free(cmap); flint_free(rp); flint_free(cp);
#ifdef DRSOLVE_MQ_STEP2_TEST
    dixon_step2_test_direct_calls++;
#endif
    return content;
}

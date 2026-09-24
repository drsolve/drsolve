/* SPDX-License-Identifier: GPL-2.0-or-later */
/* A verified Step 1 projection already is the selected square block. Extract
 * content and degree permutations from terms, then fill the final univariate
 * matrix directly. Public/non-prime callers retain their input; the private
 * native path explicitly consumes source terms and term-index arrays. */
#ifdef DRSOLVE_MQ_STEP2_TEST
static int dixon_step2_test_force_generic;
static int dixon_step2_test_direct_calls;
#endif

/* Consume/reorder the native matrix in place. If Schur succeeds, release
 * the large matrix before taking the determinant of the smaller core. */
static void dixon_mq_native_det(fq_nmod_poly_t out, nmod_poly_mat_t matrix,
                              const dixon_mq_step4_profile *p, const fq_nmod_ctx_t ctx)
{
    ulong prime=fq_nmod_ctx_prime(ctx), factor=1;
    if(p->size && p->size==matrix->r && matrix->r==matrix->c) {
        double start=get_wall_time();
        slong n=matrix->r;
        slong *at=flint_malloc(n*sizeof(slong)), *position=flint_malloc(n*sizeof(slong));
        for(int axis=0;axis<2;axis++) {
            for(slong i=0;i<n;i++) at[i]=position[i]=i;
            const slong *order=axis?p->cols:p->rows;
            for(slong i=0;i<n;i++) {
                slong j=position[order[i]];
                if(i==j) continue;
                for(slong v=0;v<n;v++)
                    nmod_poly_swap(nmod_poly_mat_entry(matrix,axis?v:i,axis?i:v),
                                   nmod_poly_mat_entry(matrix,axis?v:j,axis?j:v));
                slong tmp=at[i]; at[i]=at[j]; at[j]=tmp;
                position[at[i]]=i; position[at[j]]=j;
            }
        }
        flint_free(at); flint_free(position);
        nmod_poly_mat_t core; nmod_poly_mat_init(core,p->h,p->h,prime);
        ulong schur_factor=0;
        if(nmod_poly_mat_mq_schur(core,&schur_factor,matrix,p->rd,p->cd,p->h,p->sigma)) {
            factor=schur_factor;
            nmod_poly_mat_swap(matrix,core);
            dixon_info_log("  MQ Step 4 Schur: %ld -> %ld, compression %.3fs\n",
                           n,p->h,get_wall_time()-start);
        } else {
            dixon_info_log("  MQ Step 4 Schur: complement singular or degree check failed; using native determinant backend\n");
        }
        nmod_poly_mat_clear(core);
        if(p->odd) factor=prime-factor;
    }
    nmod_poly_t det; nmod_poly_init(det,prime);
    dixon_nmod_poly_mat_det(det,matrix);
    nmod_poly_scalar_mul_nmod(det,det,factor);
    fq_nmod_t coefficient; fq_nmod_init(coefficient,ctx); fq_nmod_poly_zero(out,ctx);
    for(slong k=0;k<det->length;k++) {
        fq_nmod_set_ui(coefficient,nmod_poly_get_coeff_ui(det,k),ctx);
        fq_nmod_poly_set_coeff(out,k,coefficient,ctx);
    }
    fq_nmod_clear(coefficient,ctx); nmod_poly_clear(det);
}

typedef struct { slong column, degree; ulong coefficient; } dixon_prime_term;

static slong dixon_projected_poly_matrix(fq_nmod_poly_mat_t out,
    nmod_poly_mat_t *prime_out, fq_mvpoly_t *consume,
    slong *rows, slong *cols, slong size, const fq_mvpoly_t *poly,
    slong *term_rows, slong *term_cols)
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
    if(prime_out) {
        FLINT_ASSERT(consume == poly && fq_nmod_ctx_degree(poly->ctx)==1);
        ulong prime=fq_nmod_ctx_prime(poly->ctx);
        dixon_prime_term *packed=flint_malloc(offset[size]*sizeof(*packed));
        /* Monomial labels were packed independently during support collection.
         * No source term is needed after this consuming conversion. */
        for(slong t=0;t<poly->nterms;t++) {
            fq_monomial_t *term=consume->terms+t;
            if(term_rows[t]>=0) {
                slong r=term_rows[t],c=term_cols[t];
                slong d=(term->par_exp?term->par_exp[0]:0)-rp[r]-cp[c];
                packed[write[r]++]=(dixon_prime_term){cmap[c],d,nmod_poly_get_coeff_ui(term->coeff,0)};
            }
            fq_nmod_clear(term->coeff,poly->ctx);
            flint_free(term->var_exp); flint_free(term->par_exp);
        }
        flint_free(consume->terms); consume->terms=NULL; consume->nterms=consume->alloc=0;
        flint_free(write); flint_free(term_rows); flint_free(term_cols);
        dixon_debug_log("  Released source Dixon terms and term maps; allocating native prime-field matrix...\n");
        nmod_poly_mat_init(*prime_out,size,size,prime);
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
                for(slong i=offset[r];i<offset[r+1];i++)
                    maximum[packed[i].column]=FLINT_MAX(maximum[packed[i].column],packed[i].degree);
                for(slong c=0;c<size;c++) if(maximum[c]>=0)
                    nmod_poly_fit_length(nmod_poly_mat_entry(*prime_out,rmap[r],c),maximum[c]+1);
                for(slong i=offset[r];i<offset[r+1];i++)
                    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(*prime_out,rmap[r],packed[i].column),
                                          packed[i].degree,packed[i].coefficient);
            }
            flint_free(maximum);
        }
        flint_free(packed); flint_free(offset);
    } else {
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
    }
    flint_free(rmap); flint_free(cmap); flint_free(rp); flint_free(cp);
#ifdef DRSOLVE_MQ_STEP2_TEST
    dixon_step2_test_direct_calls++;
#endif
    return content;
}

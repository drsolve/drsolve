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
    double phase=get_wall_time();
    slong workers=1, nt=poly->nterms;
#ifdef _OPENMP
    if(prime_out && nt>=65536 && !omp_in_parallel()) workers=omp_get_max_threads();
#endif
    slong *rp=flint_malloc(size*sizeof(slong)), *cp=flint_malloc(size*sizeof(slong));
    slong *offset=flint_calloc(size+1,sizeof(slong));
    slong *counts=flint_calloc((size_t)workers*size,sizeof(slong));
    slong *scratch=flint_malloc((size_t)2*workers*size*sizeof(slong));
    fq_index_degree_pair *rd=flint_malloc(size*sizeof(*rd)), *cd=flint_malloc(size*sizeof(*cd));
    for(slong i=0;i<size;i++) {
        rp[i]=cp[i]=WORD_MAX;
        rd[i]=(fq_index_degree_pair){i,-1}; cd[i]=(fq_index_degree_pair){i,-1};
    }
    /* Fixed contiguous chunks need no atomics and preserve per-row term order.
     * Merge row counting with valuation scanning; reuse counts as pack offsets. */
    #pragma omp parallel for if(workers>1) num_threads(workers) schedule(static)
    for(slong w=0;w<workers;w++) {
        slong *minimum=scratch+(size_t)w*size, *count=counts+(size_t)w*size;
        for(slong r=0;r<size;r++) minimum[r]=WORD_MAX;
        slong begin=(nt/workers)*w+FLINT_MIN(w,nt%workers);
        slong end=begin+nt/workers+(w<nt%workers);
        for(slong t=begin;t<end;t++) if(term_rows[t]>=0) {
            slong r=term_rows[t], d=poly->terms[t].par_exp?poly->terms[t].par_exp[0]:0;
            minimum[r]=FLINT_MIN(minimum[r],d); count[r]++;
        }
    }
    for(slong r=0;r<size;r++) {
        for(slong w=0;w<workers;w++) {
            rp[r]=FLINT_MIN(rp[r],scratch[(size_t)w*size+r]);
            offset[r+1]+=counts[(size_t)w*size+r];
        }
        if(rp[r]==WORD_MAX) rp[r]=0;
        offset[r+1]+=offset[r];
        slong position=offset[r];
        for(slong w=0;w<workers;w++) {
            slong count=counts[(size_t)w*size+r];
            counts[(size_t)w*size+r]=position; position+=count;
        }
    }
    /* Column valuations depend on completed row valuations. */
    #pragma omp parallel for if(workers>1) num_threads(workers) schedule(static)
    for(slong w=0;w<workers;w++) {
        slong *minimum=scratch+(size_t)w*size;
        for(slong c=0;c<size;c++) minimum[c]=WORD_MAX;
        slong begin=(nt/workers)*w+FLINT_MIN(w,nt%workers);
        slong end=begin+nt/workers+(w<nt%workers);
        for(slong t=begin;t<end;t++) if(term_rows[t]>=0) {
            slong c=term_cols[t];
            slong d=(poly->terms[t].par_exp?poly->terms[t].par_exp[0]:0)-rp[term_rows[t]];
            minimum[c]=FLINT_MIN(minimum[c],d);
        }
    }
    slong content=0;
    for(slong c=0;c<size;c++) {
        for(slong w=0;w<workers;w++) cp[c]=FLINT_MIN(cp[c],scratch[(size_t)w*size+c]);
        if(cp[c]==WORD_MAX) cp[c]=0;
        content+=rp[c]+cp[c];
    }
    #pragma omp parallel for if(workers>1) num_threads(workers) schedule(static)
    for(slong w=0;w<workers;w++) {
        slong *rmax=scratch+(size_t)2*w*size,*cmax=rmax+size;
        for(slong i=0;i<size;i++) rmax[i]=cmax[i]=-1;
        slong begin=(nt/workers)*w+FLINT_MIN(w,nt%workers);
        slong end=begin+nt/workers+(w<nt%workers);
        for(slong t=begin;t<end;t++) if(term_rows[t]>=0) {
            slong r=term_rows[t],c=term_cols[t];
            slong d=(poly->terms[t].par_exp?poly->terms[t].par_exp[0]:0)-rp[r]-cp[c];
            FLINT_ASSERT(d>=0);
            rmax[r]=FLINT_MAX(rmax[r],d); cmax[c]=FLINT_MAX(cmax[c],d);
        }
    }
    for(slong i=0;i<size;i++) for(slong w=0;w<workers;w++) {
        rd[i].degree=FLINT_MAX(rd[i].degree,scratch[(size_t)2*w*size+i]);
        cd[i].degree=FLINT_MAX(cd[i].degree,scratch[((size_t)2*w+1)*size+i]);
    }
    flint_free(scratch);
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
    slong *write=flint_malloc(size*sizeof(slong));
    memcpy(write,offset,size*sizeof(slong));
    dixon_debug_log("  Step 2 direct metadata: %.3fs (%ld workers)\n",get_wall_time()-phase,workers);
    phase=get_wall_time();
    if(prime_out) {
        FLINT_ASSERT(consume == poly && fq_nmod_ctx_degree(poly->ctx)==1);
        ulong prime=fq_nmod_ctx_prime(poly->ctx);
        dixon_prime_term *packed=flint_malloc(offset[size]*sizeof(*packed));
        /* Monomial labels were packed independently during support collection.
         * No source term is needed after this consuming conversion. */
        #pragma omp parallel for if(workers>1) num_threads(workers) schedule(static)
        for(slong w=0;w<workers;w++) {
            slong *cursor=counts+(size_t)w*size;
            slong begin=(nt/workers)*w+FLINT_MIN(w,nt%workers);
            slong end=begin+nt/workers+(w<nt%workers);
            for(slong t=begin;t<end;t++) {
                fq_monomial_t *term=consume->terms+t;
                if(term_rows[t]>=0) {
                    slong r=term_rows[t],c=term_cols[t];
                    slong d=(term->par_exp?term->par_exp[0]:0)-rp[r]-cp[c];
                    packed[cursor[r]++]=(dixon_prime_term){cmap[c],d,nmod_poly_get_coeff_ui(term->coeff,0)};
                }
            }
        }
        dixon_debug_log("  Step 2 direct pack: %.3fs\n",get_wall_time()-phase);
        phase=get_wall_time();
        /* Source objects often share an allocator arena. Concurrent frees
         * contend on its locks; keep destruction serial and packing parallel. */
        for(slong t=0;t<nt;t++) {
            fq_monomial_t *term=consume->terms+t;
            fq_nmod_clear(term->coeff,poly->ctx);
            flint_free(term->var_exp); flint_free(term->par_exp);
        }
        flint_free(consume->terms); consume->terms=NULL; consume->nterms=consume->alloc=0;
        flint_free(write); flint_free(term_rows); flint_free(term_cols);
        dixon_debug_log("  Released source Dixon terms and term maps; allocating native prime-field matrix...\n");
        dixon_debug_log("  Step 2 direct source release: %.3fs\n",get_wall_time()-phase);
        phase=get_wall_time();
        nmod_poly_mat_init(*prime_out,size,size,prime);
        dixon_debug_log("  Step 2 direct matrix init: %.3fs\n",get_wall_time()-phase);
        phase=get_wall_time();
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
                for(slong c=0;c<size;c++) if(maximum[c]>=0) {
                    nmod_poly_struct *entry=nmod_poly_mat_entry(*prime_out,rmap[r],c);
                    slong length=maximum[c]+1;
                    nmod_poly_fit_length(entry,length);
                    memset(entry->coeffs,0,length*sizeof(ulong)); entry->length=length;
                }
                /* Fresh, pre-sized arrays: no repeated length/gap handling
                 * for each scalar. Preserve last-write order for duplicates. */
                for(slong i=offset[r];i<offset[r+1];i++)
                    nmod_poly_mat_entry(*prime_out,rmap[r],packed[i].column)->coeffs[packed[i].degree]=packed[i].coefficient;
                for(slong c=0;c<size;c++) if(maximum[c]>=0)
                    _nmod_poly_normalise(nmod_poly_mat_entry(*prime_out,rmap[r],c));
            }
            flint_free(maximum);
        }
        dixon_debug_log("  Step 2 direct matrix fill: %.3fs\n",get_wall_time()-phase);
        phase=get_wall_time();
        flint_free(packed); flint_free(offset);
        dixon_debug_log("  Step 2 direct buffer cleanup: %.3fs\n",get_wall_time()-phase);
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
    flint_free(counts); flint_free(rmap); flint_free(cmap); flint_free(rp); flint_free(cp);
#ifdef DRSOLVE_MQ_STEP2_TEST
    dixon_step2_test_direct_calls++;
#endif
    return content;
}

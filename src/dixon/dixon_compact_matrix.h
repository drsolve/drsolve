/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Private, included after the projected matrix and Schur helpers. */
static int dixon_compact_enabled(slong npars)
{
    const char *setting=getenv("DRSOLVE_MQ_COMPACT");
    return (!setting || strcmp(setting,"0")!=0) && npars==1 &&
        (dixon_global_method_step4==-1 ||
         (g_dixon_mq_step4_schur && dixon_global_method_step4==DET_METHOD_KRONECKER));
}

static void dixon_print_compact(const fq_mq_compact *p,const fq_nmod_ctx_t ctx,
                               char **vars,char **pars,const char *gen)
{
    dixon_info_log("  Dixon polynomial: %ld terms (compact row storage)\n",p->nterms);
    if(g_dixon_verbose_level>=1 && p->nterms<=100) {
        fq_mvpoly_t small; fq_mq_compact_materialize(&small,p,ctx);
        fq_mvpoly_print_with_names(&small,"  DixonPoly",vars,pars,gen,1);
        fq_mvpoly_clear(&small);
    }
    if(g_dixon_debug_mode) {
        slong *degree=flint_calloc(2*p->nvars,sizeof(slong));
        for(slong r=0;r<p->nrows;r++) for(slong v=0;v<p->nvars;v++)
            degree[v]=FLINT_MAX(degree[v],p->rows[r*p->nvars+v]);
        for(slong c=0;c<p->ncols;c++) for(slong v=0;v<p->nvars;v++)
            degree[p->nvars+v]=FLINT_MAX(degree[p->nvars+v],p->cols[c*p->nvars+v]);
        printf("  Dixon polynomial actual degrees:\n");
        print_fq_degree_vector_with_names("original vars",degree,p->nvars,vars,0,0);
        print_fq_degree_vector_with_names("dual vars",degree+p->nvars,p->nvars,vars,1,0);
        print_fq_degree_vector_with_names("parameter vars",&p->degree,1,pars,0,1);
        flint_free(degree);
    }
}

/* Same content removal and degree tie-breaking as the conventional extractor.
 * Input records are already grouped by row; no per-term index or pack buffer. */
static slong dixon_compact_poly_matrix(nmod_poly_mat_t out,
    slong *rows,slong *cols,const fq_mq_compact *p)
{
    slong size=p->nrows;
    FLINT_ASSERT(size==p->ncols && size>0);
    double phase=get_wall_time();
    slong *rp=flint_malloc(size*sizeof(slong)),*cp=flint_malloc(size*sizeof(slong));
    fq_index_degree_pair *rd=flint_malloc(size*sizeof(*rd)),*cd=flint_malloc(size*sizeof(*cd));
    for(slong i=0;i<size;i++) {
        rp[i]=cp[i]=WORD_MAX;
        rd[i]=(fq_index_degree_pair){i,-1}; cd[i]=(fq_index_degree_pair){i,-1};
    }
    #pragma omp parallel for if(size>1) schedule(static)
    for(slong r=0;r<size;r++) {
        for(slong t=p->offset[r];t<p->offset[r+1];t++)
            rp[r]=FLINT_MIN(rp[r],p->terms[t].degree);
        if(rp[r]==WORD_MAX) rp[r]=0;
    }
    for(slong r=0;r<size;r++) for(slong t=p->offset[r];t<p->offset[r+1];t++) {
        const fq_mq_compact_term *a=p->terms+t;
        cp[a->column]=FLINT_MIN(cp[a->column],a->degree-rp[r]);
    }
    slong content=0;
    for(slong c=0;c<size;c++) {
        if(cp[c]==WORD_MAX) cp[c]=0;
        content+=rp[c]+cp[c];
    }
    for(slong r=0;r<size;r++) for(slong t=p->offset[r];t<p->offset[r+1];t++) {
        const fq_mq_compact_term *a=p->terms+t;
        slong d=a->degree-rp[r]-cp[a->column];
        FLINT_ASSERT(d>=0);
        rd[r].degree=FLINT_MAX(rd[r].degree,d);
        cd[a->column].degree=FLINT_MAX(cd[a->column].degree,d);
    }
    const char *reorder=getenv("DRSOLVE_PREDICT_REORDER");
    if(!reorder || strcmp(reorder,"0")!=0) {
        qsort(rd,size,sizeof(*rd),compare_fq_degrees);
        qsort(cd,size,sizeof(*cd),compare_fq_degrees);
    }
    slong *rmap=flint_malloc(size*sizeof(slong)),*cmap=flint_malloc(size*sizeof(slong));
    for(slong i=0;i<size;i++) {
        rows[i]=rd[i].index; cols[i]=cd[i].index;
        rmap[rows[i]]=i; cmap[cols[i]]=i;
    }
    flint_free(rd); flint_free(cd);
    dixon_debug_log("  Step 2 compact metadata: %.3fs\n",get_wall_time()-phase);
    phase=get_wall_time(); nmod_poly_mat_init(out,size,size,p->prime);
    dixon_debug_log("  Step 2 compact matrix init: %.3fs\n",get_wall_time()-phase);
    phase=get_wall_time();
    #pragma omp parallel if(size>1)
    {
        slong *maximum=flint_malloc(size*sizeof(slong));
        #pragma omp for schedule(dynamic,1)
        for(slong r=0;r<size;r++) {
            for(slong c=0;c<size;c++) maximum[c]=-1;
            for(slong t=p->offset[r];t<p->offset[r+1];t++) {
                const fq_mq_compact_term *a=p->terms+t;
                slong d=a->degree-rp[r]-cp[a->column];
                maximum[a->column]=FLINT_MAX(maximum[a->column],d);
            }
            for(slong c=0;c<size;c++) if(maximum[c]>=0) {
                nmod_poly_struct *entry=nmod_poly_mat_entry(out,rmap[r],cmap[c]);
                nmod_poly_fit_length(entry,maximum[c]+1);
                entry->length=maximum[c]+1;
                memset(entry->coeffs,0,entry->length*sizeof(ulong));
            }
            for(slong t=p->offset[r];t<p->offset[r+1];t++) {
                const fq_mq_compact_term *a=p->terms+t;
                nmod_poly_mat_entry(out,rmap[r],cmap[a->column])->coeffs[a->degree-rp[r]-cp[a->column]]=a->coefficient;
            }
            for(slong c=0;c<size;c++) if(maximum[c]>=0)
                _nmod_poly_normalise(nmod_poly_mat_entry(out,rmap[r],cmap[c]));
        }
        flint_free(maximum);
    }
    dixon_debug_log("  Step 2 compact matrix fill: %.3fs\n",get_wall_time()-phase);
    flint_free(rmap); flint_free(cmap); flint_free(rp); flint_free(cp);
    return content;
}

static void dixon_extract_compact_matrix(nmod_poly_mat_t out,slong *size,
    slong *content,fq_mq_compact *p,const long *degrees,
    dixon_mq_step4_profile *profile,char **par_names)
{
    clock_t cpu=clock(); double wall=get_wall_time();
    dixon_info_log("\nStep 2: Construct Dixon matrix\n");
    *size=p->nrows;
    dixon_info_log("  Dixon matrix size: %ld x %ld\n",p->nrows,p->ncols);
    dixon_debug_log("  Using compact Step 1 rows directly; no support scan or term repacking\n");
    slong *rows=flint_malloc(*size*sizeof(slong)),*cols=flint_malloc(*size*sizeof(slong));
    *content=dixon_compact_poly_matrix(out,rows,cols,p);
    if(*content>0) {
        const char *v=par_names && par_names[0] ? par_names[0] : "x";
        dixon_info_log("  Pre-selection full-matrix %s-content: %s^%ld\n",v,v,*content);
    }
    /* Release the only term-sized buffer before Step 3. Support labels live
     * until the profile has captured the degree ordering and sign. */
    double phase=get_wall_time(); flint_free(p->terms); p->terms=NULL;
    dixon_debug_log("  Step 2 compact buffer release: %.3fs\n",get_wall_time()-phase);
    dixon_maybe_print_parallel_step_time("Step 2",(double)(clock()-cpu)/CLOCKS_PER_SEC,get_wall_time()-wall);
    cpu=clock(); wall=get_wall_time();
    dixon_info_log("\nStep 3: Extract maximal-rank submatrix\n");
    dixon_debug_log("  Using MQ candidate verified in Step 1\n");
    dixon_info_log("  Submatrix size: %ld x %ld\n",*size,*size);
    if(profile) {
        monom_t *rm=flint_malloc(*size*sizeof(*rm)),*cm=flint_malloc(*size*sizeof(*cm));
        for(slong i=0;i<*size;i++) {
            rm[i].exp=p->rows+i*p->nvars; rm[i].idx=i;
            cm[i].exp=p->cols+i*p->nvars; cm[i].idx=i;
        }
        dixon_mq_step4_prepare(profile,rm,cm,rows,cols,*size,p->nvars,degrees);
        flint_free(rm); flint_free(cm);
    }
    flint_free(rows); flint_free(cols); fq_mq_compact_clear(p);
    dixon_maybe_print_parallel_step_time("Step 3",(double)(clock()-cpu)/CLOCKS_PER_SEC,get_wall_time()-wall);
}

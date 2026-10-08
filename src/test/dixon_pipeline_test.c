/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Independent cancellation polynomial checks the recursive labels; determinant
 * comparisons cover shared selection/dispatch and construction-specific paths. */
#define main recursive_arithmetic_fixture_main
#include "dixon_recursive_native_test.c"
#undef main

static int pipeline_equal(const unified_mpoly_struct *a, const unified_mpoly_struct *b)
{
    fq_nmod_mpoly_ctx_t ctx;
    fq_nmod_mpoly_ctx_init(ctx,a->nvars+a->npars,ORD_LEX,a->ctx);
    fq_nmod_mpoly_t x,y;
    fq_nmod_mpoly_init(x,ctx);fq_nmod_mpoly_init(y,ctx);
    dr_mpoly_to_fq_nmod_mpoly(x,a,ctx);dr_mpoly_to_fq_nmod_mpoly(y,b,ctx);
    int same=fq_nmod_mpoly_equal(x,y,ctx);
    fq_nmod_mpoly_clear(x,ctx);fq_nmod_mpoly_clear(y,ctx);fq_nmod_mpoly_ctx_clear(ctx);
    return same;
}

static void check_dense_content(void)
{
    fq_nmod_ctx_t ctx;fq_nmod_ctx_init_ui(ctx,101,1,"a");
    unified_mpoly_struct entries[2][2]={0}, *matrix[2]={entries[0],entries[1]};
    unified_mpoly_struct polys[3]={0},result={0},expected={0};
    fq_nmod_t one;fq_nmod_init(one,ctx);fq_nmod_one(one,ctx);
    for(slong i=0;i<2;i++) for(slong j=0;j<2;j++) dr_mpoly_init(&entries[i][j],0,1,ctx);
    for(slong i=0;i<2;i++) {
        slong power=i+1;
        dr_mpoly_add_term_fast(&entries[i][i],NULL,&power,one);
    }
    for(slong i=0;i<3;i++) {
        slong x[2]={1,1},power=1;
        dr_mpoly_init(polys+i,2,1,ctx);dr_mpoly_add_term_fast(polys+i,x,&power,one);
    }
    dixon_global_method_step4=DET_METHOD_RECURSIVE;
    dixon_compute_dense_resultant(&result,matrix,2,polys,2,1,2,NULL,NULL);
    dr_mpoly_init(&expected,0,1,ctx);
    slong power=5;dr_mpoly_add_term_fast(&expected,NULL,&power,one);
    assert(pipeline_equal(&result,&expected));
    dr_mpoly_clear(&result);dr_mpoly_clear(&expected);
    for(slong i=0;i<3;i++) dr_mpoly_clear(polys+i);
    for(slong i=0;i<2;i++) for(slong j=0;j<2;j++) dr_mpoly_clear(&entries[i][j]);
    fq_nmod_clear(one,ctx);fq_nmod_ctx_clear(ctx);
}

static void check_scalar(slong extension)
{
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init_ui(ctx, 101, extension, "a");
    fast_dixon_matrix_t full;
    fast_dixon_matrix_init(&full, 3, 3, 0, ctx);
    const int values[3][3]={{1,1,0},{1,1,0},{0,0,1}};
    fq_nmod_t one; fq_nmod_init(one,ctx); fq_nmod_gen(one,ctx);
    for (slong i=0;i<3;i++) for (slong j=0;j<3;j++)
        if (values[i][j]) dr_mpoly_add_term_fast(&FAST_DIXON_ENTRY(&full,i,j),NULL,NULL,one);
    unified_mpoly_struct ***grid=fast_dixon_build_pointer_grid(&full);
    slong *rows=NULL,*cols=NULL,size=0;
    dixon_select_submatrix(grid,3,3,NULL,NULL,0,0,NULL,0,-1,ctx,&rows,&cols,&size);
    assert(size==2);
    fq_nmod_mat_t selected;
    fq_nmod_mat_init(selected,size,size,ctx);
    for (slong i=0;i<size;i++) for (slong j=0;j<size;j++)
        if (values[rows[i]][cols[j]]) fq_nmod_set(fq_nmod_mat_entry(selected,i,j),one,ctx);
    assert(fq_nmod_mat_rank(selected,ctx)==size);
    fq_nmod_mat_clear(selected,ctx);
    flint_free(rows); flint_free(cols); fast_dixon_free_pointer_grid(grid,3);
    fast_dixon_matrix_clear(&full); fq_nmod_clear(one,ctx); fq_nmod_ctx_clear(ctx);
}

static void check_pipeline(ulong prime, slong extension, slong nvars, slong npars)
{
    fq_nmod_ctx_t ctx; fq_nmod_ctx_init_ui(ctx,prime,extension,"a");
    flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng,1729+nvars,43);
    unified_mpoly_struct p[4]={0}, result[4]={0};
    const unified_mpoly_struct *ptrs[4];
    slong exp[5]={0}, degrees[3];
    fq_nmod_t one; fq_nmod_init(one,ctx); fq_nmod_gen(one,ctx);
    for (slong i=0;i<=nvars;i++) {
        dr_mpoly_init(p+i,nvars,npars,ctx); ptrs[i]=p+i;
        random_terms(p+i,exp,0,2,rng,0);
        for (slong v=0;v<nvars;v++) {
            slong x[3]={0}, par[2]={0}; x[v]=2;
            dr_mpoly_add_term_fast(p+i,x,par,one);
        }
    }
    fast_dixon_compute_degree_bounds(degrees,ptrs,nvars+1,nvars);
    /* Verify all raw labels against an independently constructed Dixon poly. */
    fast_dixon_matrix_t full;
    fast_dixon_profile_reset(nvars); fast_dixon_subproblem_cache_reset_all();
    fast_dixon_build_matrix(&full,ptrs,degrees,nvars,0,npars,ctx);
    unified_mpoly_struct actual={0},expected={0},**cancel,**modified;
    dr_mpoly_init(&actual,2*nvars,npars,ctx);
    slong rc[3],cc[3],labels[6];
    for(slong v=0;v<nvars;v++) {rc[v]=(v+1)*degrees[v];cc[v]=(nvars-v)*degrees[v];}
    for(slong r=0;r<full.rows;r++) for(slong c=0;c<full.cols;c++) {
        fast_dixon_decode_rectangular_index(r,rc,nvars,labels);
        fast_dixon_decode_rectangular_index(c,cc,nvars,labels+nvars);
        const unified_mpoly_struct *entry=&FAST_DIXON_ENTRY(&full,r,c);
        for(slong t=0;t<dr_mpoly_length(entry);t++) {
            DR_MPOLY_TERM(term,entry,t);
            dr_mpoly_add_term_fast(&actual,labels,term.par_exp,term.coeff);
        }
    }
    build_fq_cancellation_matrix(&cancel,p,nvars,npars);
    perform_fq_matrix_row_operations(&modified,&cancel,nvars,npars);
    compute_fq_cancel_matrix_det(&expected,modified,nvars,npars,DET_METHOD_RECURSIVE);
    if(!pipeline_equal(&actual,&expected)) {dr_mpoly_neg(&actual,&actual);assert(pipeline_equal(&actual,&expected));}
    for(slong i=0;i<=nvars;i++) {
        for(slong j=0;j<=nvars;j++) {dr_mpoly_clear(&cancel[i][j]);dr_mpoly_clear(&modified[i][j]);}
        flint_free(cancel[i]);flint_free(modified[i]);
    }
    flint_free(cancel);flint_free(modified);dr_mpoly_clear(&actual);dr_mpoly_clear(&expected);
    fast_dixon_matrix_clear(&full);

    /* Standard MQ/projection, recursive prediction, prediction-disabled fallback,
     * and recursive legacy arithmetic must agree after monic normalization. */
    dixon_global_method_step4=DET_METHOD_KRONECKER;
    setenv("DRSOLVE_PREDICT_MAXRANK","1",1);
    fq_dixon_resultant(result,p,nvars,npars);
    setenv("DRSOLVE_FAST_NATIVE","1",1);
    fq_dixon_fast_resultant(result+1,p,nvars,npars);
    setenv("DRSOLVE_PREDICT_MAXRANK","0",1);
    fq_dixon_fast_resultant(result+2,p,nvars,npars);
    setenv("DRSOLVE_FAST_NATIVE","0",1);
    setenv("DRSOLVE_PREDICT_MAXRANK","1",1);
    fq_dixon_fast_resultant_with_names(result+3,p,nvars,npars,NULL,NULL,NULL);
    for(slong i=1;i<4;i++) assert(pipeline_equal(result,result+i));
    for(slong i=0;i<4;i++) dr_mpoly_clear(result+i);
    for(slong i=0;i<=nvars;i++) dr_mpoly_clear(p+i);
    fq_nmod_clear(one,ctx); flint_rand_clear(rng); fq_nmod_ctx_clear(ctx);
}

int main(void)
{
    g_dixon_verbose_level=getenv("DRSOLVE_PIPELINE_TEST_VERBOSE") ? 3 : 0;
    omp_set_num_threads(2);
    check_scalar(1);check_scalar(2);check_dense_content();
    check_pipeline(65537,1,2,1);
    check_pipeline(101,1,3,1);
    check_pipeline(101,1,2,2);
    check_pipeline(101,2,2,1);
    check_pipeline(101,1,2,0);
    unsetenv("DRSOLVE_FAST_NATIVE");unsetenv("DRSOLVE_PREDICT_MAXRANK");
    dixon_global_method_step4=-1;
    puts("Dixon shared pipeline tests passed");
    return 0;
}

/* SPDX-License-Identifier: GPL-2.0-or-later */
#define DRSOLVE_MQ_STEP2_TEST 1
#define main projected_matrix_fixture_main
#include "dixon_mq_filter_test.c"
#undef main

static void compare_paths(const fq_mvpoly_t *poly, slong m, const long *degrees,
                          int profile, slong expected_content)
{
    fq_mvpoly_t **unused=NULL;
    fq_nmod_poly_mat_t expected, actual;
    slong *r[2],*c[2],size[2],content[2];
    dixon_mq_step4_profile p[2]={{0},{0}};
    for(int method=0;method<2;method++) {
        r[method]=flint_malloc(poly->nterms*sizeof(slong));
        c[method]=flint_malloc(poly->nterms*sizeof(slong));
        dixon_step2_test_force_generic=!method;
        int before=dixon_step2_test_direct_calls;
        extract_fq_coefficient_matrix_from_dixon_impl(&unused,method?&actual:&expected,
            r[method],c[method],size+method,content+method,poly,m,1,NULL,NULL,NULL,
            degrees,m+1,1,profile?p+method:NULL);
        assert(dixon_step2_test_direct_calls-before==method);
        assert(!unused);
    }
    assert(size[0]==size[1] && content[0]==content[1]);
    if(expected_content>=0) assert(content[0]==expected_content);
    assert(!memcmp(r[0],r[1],size[0]*sizeof(slong)));
    assert(!memcmp(c[0],c[1],size[0]*sizeof(slong)));
    for(slong i=0;i<size[0];i++) for(slong j=0;j<size[0];j++)
        assert(fq_nmod_poly_equal(fq_nmod_poly_mat_entry(expected,i,j),
                                  fq_nmod_poly_mat_entry(actual,i,j),poly->ctx));
    assert(p[0].size==p[1].size && p[0].h==p[1].h && p[0].sigma==p[1].sigma && p[0].odd==p[1].odd);
    if(p[0].size) {
        assert(!memcmp(p[0].rows,p[1].rows,p[0].size*sizeof(slong)));
        assert(!memcmp(p[0].cols,p[1].cols,p[0].size*sizeof(slong)));
        assert(!memcmp(p[0].rd,p[1].rd,p[0].size*sizeof(slong)));
        assert(!memcmp(p[0].cd,p[1].cd,p[0].size*sizeof(slong)));
    }
    fq_nmod_poly_mat_clear(expected,poly->ctx); fq_nmod_poly_mat_clear(actual,poly->ctx);
    for(int i=0;i<2;i++) { flint_free(r[i]); flint_free(c[i]); dixon_mq_step4_profile_clear(p+i); }
    dixon_step2_test_force_generic=0;
}
static void synthetic(ulong prime, slong extension)
{
    fq_nmod_ctx_t ctx; fq_nmod_ctx_init_ui(ctx,prime,extension,"a");
    fq_mvpoly_t poly; fq_mvpoly_init(&poly,2,1,ctx);
    fq_nmod_t coeff; fq_nmod_init(coeff,ctx);
    fq_nmod_gen(coeff,ctx); if(fq_nmod_is_zero(coeff,ctx)) fq_nmod_one(coeff,ctx);
    slong powers[2][2]={{5,7},{4,5}};
    for(slong r=0;r<2;r++) for(slong c=0;c<2;c++) {
        slong e[2]={r,c}; fq_mvpoly_add_term_fast(&poly,e,&powers[r][c],coeff);
    }
    /* Original parameter exponents must remain caller-owned and unchanged. */
    compare_paths(&poly,1,NULL,0,10);
    for(slong r=0;r<2;r++) for(slong c=0;c<2;c++) assert(poly.terms[2*r+c].par_exp[0]==powers[r][c]);
    fq_mvpoly_clear(&poly); fq_nmod_clear(coeff,ctx); fq_nmod_ctx_clear(ctx);
}
int main(void)
{
    g_dixon_verbose_level=0; g_dixon_det_cache_limit=100000;
    flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng,132,941);
    for(int reorder=0;reorder<2;reorder++) {
        setenv("DRSOLVE_PREDICT_REORDER",reorder?"1":"0",1);
        for(int threads=1;threads<=4;threads+=3) {
            omp_set_num_threads(threads);
            synthetic(257,1); synthetic(7,2);
            fq_nmod_ctx_t ctx; fq_nmod_ctx_init_ui(ctx,257,1,"a");
            for(slong m=3;m<=5;m++) {
                fq_mvpoly_t *polys=random_mq(m,ctx,rng), **matrix, **a, projected;
                build_fq_cancellation_matrix_mvpoly(&matrix,polys,m,1);
                perform_fq_matrix_row_operations_mvpoly(&a,&matrix,m,1);
                assert(dixon_try_mq_projection(&projected,a,polys,m,1,DET_METHOD_RECURSIVE));
                long degrees[6]; for(slong i=0;i<=m;i++) degrees[i]=2;
                compare_paths(&projected,m,degrees,1,-1);
                fq_mvpoly_clear(&projected); clear_input(polys,matrix,a,m);
            }
            fq_nmod_ctx_clear(ctx);
        }
    }
    unsetenv("DRSOLVE_PREDICT_REORDER");
    flint_rand_clear(rng); flint_cleanup_master();
    puts("Direct projected matrix: 20 exact matrix/content/order/profile comparisons PASS");
    return 0;
}

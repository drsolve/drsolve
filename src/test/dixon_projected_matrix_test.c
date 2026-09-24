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
            degrees,m+1,1,profile?p+method:NULL,NULL,NULL);
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
    if(fq_nmod_ctx_degree(poly->ctx)==1) {
        fq_mvpoly_t owned; fq_mvpoly_init(&owned,poly->nvars,poly->npars,poly->ctx);
        fq_mvpoly_copy(&owned,poly);
        nmod_poly_mat_t native;
        slong *rn=flint_malloc(poly->nterms*sizeof(slong)),*cn=flint_malloc(poly->nterms*sizeof(slong));
        slong ns,ncontent; dixon_mq_step4_profile np={0};
        dixon_step2_test_force_generic=0;
        extract_fq_coefficient_matrix_from_dixon_impl(&unused,NULL,rn,cn,&ns,&ncontent,
            &owned,m,1,NULL,NULL,NULL,degrees,m+1,1,profile?&np:NULL,&native,&owned);
        assert(!owned.terms && !owned.nterms && !owned.alloc);
        fq_mvpoly_clear(&owned); /* Consumed inputs remain safely clearable. */
        assert(ns==size[0] && ncontent==content[0]);
        assert(!memcmp(rn,r[0],ns*sizeof(slong)) && !memcmp(cn,c[0],ns*sizeof(slong)));
        for(slong i=0;i<ns;i++) for(slong j=0;j<ns;j++) {
            nmod_poly_struct *a=nmod_poly_mat_entry(native,i,j);
            fq_nmod_poly_struct *b=fq_nmod_poly_mat_entry(expected,i,j);
            assert(a->length==b->length);
            for(slong d=0;d<a->length;d++) assert(a->coeffs[d]==nmod_poly_get_coeff_ui(b->coeffs+d,0));
        }
        assert(np.size==p[0].size && np.h==p[0].h && np.sigma==p[0].sigma && np.odd==p[0].odd);
        if(np.size) {
            assert(!memcmp(np.rows,p[0].rows,ns*sizeof(slong)));
            assert(!memcmp(np.cols,p[0].cols,ns*sizeof(slong)));
            assert(!memcmp(np.rd,p[0].rd,ns*sizeof(slong)));
            assert(!memcmp(np.cd,p[0].cd,ns*sizeof(slong)));
        }
        for(int method=0;method<3;method++) {
            fq_nmod_poly_mat_det_set_method(method);
            fq_nmod_poly_t want,got; fq_nmod_poly_init(want,poly->ctx); fq_nmod_poly_init(got,poly->ctx);
            fq_nmod_poly_mat_det_iter(want,expected,poly->ctx);
            for(int fallback=0;fallback<2;fallback++) {
                nmod_poly_mat_t work; nmod_poly_mat_init(work,ns,ns,fq_nmod_ctx_prime(poly->ctx));
                nmod_poly_mat_set(work,native);
                slong sigma=np.sigma; if(fallback) np.sigma=0;
                dixon_mq_native_det(got,work,&np,poly->ctx); np.sigma=sigma;
                assert(fq_nmod_poly_equal(want,got,poly->ctx));
                nmod_poly_mat_clear(work);
            }
            fq_nmod_poly_clear(want,poly->ctx); fq_nmod_poly_clear(got,poly->ctx);
        }
        fq_nmod_poly_mat_det_set_method(FQ_NMOD_POLY_DET_METHOD_AUTO);
        nmod_poly_mat_clear(native); dixon_mq_step4_profile_clear(&np); flint_free(rn); flint_free(cn);
    }
    fq_nmod_poly_mat_clear(expected,poly->ctx); fq_nmod_poly_mat_clear(actual,poly->ctx);
    for(int i=0;i<2;i++) { flint_free(r[i]); flint_free(c[i]); dixon_mq_step4_profile_clear(p+i); }
    dixon_step2_test_force_generic=0;
}
static void synthetic(ulong prime, slong extension, slong repeats)
{
    fq_nmod_ctx_t ctx; fq_nmod_ctx_init_ui(ctx,prime,extension,"a");
    fq_mvpoly_t poly; fq_mvpoly_init(&poly,2,1,ctx);
    fq_nmod_t coeff; fq_nmod_init(coeff,ctx);
    fq_nmod_gen(coeff,ctx); if(fq_nmod_is_zero(coeff,ctx)) fq_nmod_one(coeff,ctx);
    slong powers[2][2]={{5,7},{4,5}};
    for(slong repeat=0;repeat<repeats;repeat++) {
        if(repeats>1) fq_nmod_set_ui(coeff,1+(repeat%(prime-1)),ctx);
        for(slong r=0;r<2;r++) for(slong c=0;c<2;c++) {
            slong e[2]={r,c}; fq_mvpoly_add_term_fast(&poly,e,&powers[r][c],coeff);
        }
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
            synthetic(257,1,1); synthetic(7,2,1);
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
    /* Cross the parallel metadata threshold with duplicate terms, checking
     * deterministic last-write ordering and consuming cleanup. The row length
     * is not a multiple of 32, exercising partial staging flushes as well. */
    omp_set_num_threads(4); synthetic(257,1,20001);
    unsetenv("DRSOLVE_PREDICT_REORDER");
    flint_rand_clear(rng); flint_cleanup_master();
    puts("Direct projected matrix: 21 generic/direct comparisons, 17 consuming native comparisons and 102 native determinant checks PASS");
    return 0;
}

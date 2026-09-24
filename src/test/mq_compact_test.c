/* SPDX-License-Identifier: GPL-2.0-or-later */
#define main compact_fixture_main
#include "dixon_mq_filter_test.c"
#undef main

static void same_polynomial(const fq_mvpoly_t *a,const fq_mvpoly_t *b)
{
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx,a->nvars+a->npars,ORD_LEX,fq_nmod_ctx_prime(a->ctx));
    nmod_mpoly_t x,y; nmod_mpoly_init(x,ctx); nmod_mpoly_init(y,ctx);
    fq_mvpoly_to_nmod_mpoly(x,a,ctx); fq_mvpoly_to_nmod_mpoly(y,b,ctx);
    assert(nmod_mpoly_equal(x,y,ctx));
    nmod_mpoly_clear(x,ctx); nmod_mpoly_clear(y,ctx); nmod_mpoly_ctx_clear(ctx);
}

static void same_matrix(fq_mq_compact *compact,fq_mvpoly_t *legacy,slong n,int reorder)
{
    setenv("DRSOLVE_PREDICT_REORDER",reorder ? "1" : "0",1);
    fq_mvpoly_t **unused=NULL, owned;
    fq_mvpoly_init(&owned,legacy->nvars,1,legacy->ctx); fq_mvpoly_copy(&owned,legacy);
    nmod_poly_mat_t expected,actual;
    slong size,content,cs,cc;
    long degrees[8]; for(slong i=0;i<=n;i++) degrees[i]=2;
    dixon_mq_step4_profile a={0},b={0};
    extract_fq_coefficient_matrix_from_dixon_impl(&unused,NULL,NULL,NULL,&size,&content,
        &owned,n,1,NULL,NULL,NULL,degrees,n+1,1,&a,&expected,&owned);
    dixon_extract_compact_matrix(actual,&cs,&cc,compact,degrees,&b,NULL);
    assert(!compact->nvars && !compact->terms && !compact->rows && !compact->offset);
    assert(size==cs && content==cc);
    assert(nmod_poly_mat_equal(expected,actual));
    assert(a.size==b.size && a.h==b.h && a.sigma==b.sigma && a.odd==b.odd);
    if(a.size) {
        assert(!memcmp(a.rows,b.rows,size*sizeof(slong)));
        assert(!memcmp(a.cols,b.cols,size*sizeof(slong)));
        assert(!memcmp(a.rd,b.rd,size*sizeof(slong)));
        assert(!memcmp(a.cd,b.cd,size*sizeof(slong)));
    }
    dixon_mq_step4_profile_clear(&a); dixon_mq_step4_profile_clear(&b);
    nmod_poly_mat_clear(expected); nmod_poly_mat_clear(actual); fq_mvpoly_clear(&owned);
    unsetenv("DRSOLVE_PREDICT_REORDER");
}

int main(void)
{
    g_dixon_verbose_level=0; g_dixon_det_cache_limit=100000;
    flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng,132,941);
    ulong primes[]={2,3,257};
    int native=0,fallback=0,repaired=0;
    for(int threads=1;threads<=4;threads+=3) for(int pi=0;pi<3;pi++) for(slong n=2;n<=5;n++) {
        omp_set_num_threads(threads);
        fq_nmod_ctx_t ctx; fq_nmod_ctx_init_ui(ctx,primes[pi],1,"a");
        fq_mvpoly_t *polys=random_mq(n,ctx,rng),**matrix,**m,old,placeholder;
        build_fq_cancellation_matrix_mvpoly(&matrix,polys,n,1);
        perform_fq_matrix_row_operations_mvpoly(&m,&matrix,n,1);
        /* Rectangular output and its fallback materialization are exact. */
        slong rows[15]={0},cols[15]={0}; rows[n]=2; rows[2*n+1]=1;
        fq_mq_compact compact={0};
        assert(compute_fq_det_mq_projected_rect(&old,m,n+1,rows,3,cols,1));
        assert(compute_fq_det_mq_compact(&compact,m,n+1,rows,3,cols,1));
        fq_mq_compact_materialize(&placeholder,&compact,ctx);
        same_polynomial(&old,&placeholder);
        fq_mvpoly_clear(&old); fq_mvpoly_clear(&placeholder); fq_mq_compact_clear(&compact);
        assert(!compute_fq_det_mq_compact(&compact,m,n+1,rows,0,cols,1));
        assert(!compact.nvars && !compact.terms);
        int want=dixon_try_mq_projection(&old,m,polys,n,1,DET_METHOD_RECURSIVE);
        int got=dixon_try_mq_projection_from_full(&placeholder,m,polys,n,1,DET_METHOD_RECURSIVE,NULL,&compact);
        assert(want==got);
        if(got) {
            if(compact.nvars) {
                assert(placeholder.nterms==0); /* No generic per-term objects. */
                fq_mvpoly_t expanded; fq_mq_compact_materialize(&expanded,&compact,ctx);
                same_polynomial(&old,&expanded); fq_mvpoly_clear(&expanded);
                same_matrix(&compact,&old,n,threads==4); native++;
            } else { same_polynomial(&old,&placeholder); repaired++; }
            fq_mvpoly_clear(&old); fq_mvpoly_clear(&placeholder);
        } else { assert(!compact.nvars && !compact.terms); fallback++; }
        fq_mq_compact_clear(&compact);
        clear_input(polys,matrix,m,n); fq_nmod_ctx_clear(ctx);
    }
    /* A structurally eligible, identically singular system must discard the
     * compact candidate and return a conventional full polynomial. */
    fq_nmod_ctx_t ctx; fq_nmod_ctx_init_ui(ctx,257,1,"a");
    fq_mvpoly_t *polys=random_mq(3,ctx,rng),**matrix,**m,result;
    fq_mvpoly_copy(polys+1,polys);
    build_fq_cancellation_matrix_mvpoly(&matrix,polys,3,1);
    perform_fq_matrix_row_operations_mvpoly(&m,&matrix,3,1);
    fq_mq_compact compact={0};
    assert(!dixon_compute_step1(&result,m,polys,3,1,DET_METHOD_RECURSIVE,&compact));
    assert(!compact.nvars && !compact.terms && !result.nterms);
    fq_mvpoly_clear(&result); fq_mq_compact_clear(&compact);
    clear_input(polys,matrix,m,3); fq_nmod_ctx_clear(ctx);
    /* Positive row/column valuations, coefficient gaps and both reorder modes. */
    fq_nmod_ctx_init_ui(ctx,257,1,"a");
    for(int reorder=0;reorder<2;reorder++) {
        fq_mq_compact p={.nvars=1,.nrows=2,.ncols=2,.nterms=4,.degree=7,.prime=257};
        p.rows=flint_malloc(2*sizeof(slong)); p.cols=flint_malloc(2*sizeof(slong));
        p.offset=flint_malloc(3*sizeof(slong)); p.terms=flint_malloc(4*sizeof(*p.terms));
        p.rows[0]=p.cols[0]=0; p.rows[1]=p.cols[1]=1;
        p.offset[0]=0; p.offset[1]=2; p.offset[2]=4;
        p.terms[0]=(fq_mq_compact_term){0,5,1}; p.terms[1]=(fq_mq_compact_term){1,7,2};
        p.terms[2]=(fq_mq_compact_term){0,4,3}; p.terms[3]=(fq_mq_compact_term){1,5,4};
        fq_mvpoly_t expanded; fq_mq_compact_materialize(&expanded,&p,ctx);
        same_matrix(&p,&expanded,1,reorder); fq_mvpoly_clear(&expanded);
    }
    fq_nmod_ctx_clear(ctx);
    assert(native>0 && fallback>0);
    printf("Compact MQ: %d native matrices/profiles, %d repaired, %d rejected; rectangular and singular fallback PASS\n",native,repaired,fallback);
    flint_rand_clear(rng); flint_cleanup_master();
    return 0;
}

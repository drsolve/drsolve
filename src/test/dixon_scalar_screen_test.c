/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include "dixon_interface_flint.h"
#include "../dixon/dixon_pipeline.h"

static void check_system(const char *polys, const fq_nmod_ctx_t ctx,
                         dixon_scalar_screen_status_t expected)
{
    dixon_scalar_screen_report_t report;
    char *result = dixon_str_with_scalar_screen(polys, "x,y", ctx, &report);
    assert(result);
    assert(strcmp(result, expected == DIXON_SCREEN_NO_COMMON_ZERO ? "1" : "0")==0);
    free(result);
    assert(report.status == expected);
    assert(!dixon_scalar_screen_active());
    if (expected == DIXON_SCREEN_NO_COMMON_ZERO) {
        assert(report.constant_col >= 0);
        assert(report.rank == report.nonconstant_rank + 1);
    } else {
        assert(report.rank == report.nonconstant_rank);
    }
}

/* Column labels are deliberately permuted: the constant is in the middle.
 * Rows annihilate the monomial vector [a,1,a^2], with a a field generator.
 * This also catches accidental projection of extension coefficients to F_p. */
static void check_matrix(const fq_nmod_ctx_t ctx)
{
    unified_mpoly_struct entries[2][3] = {0};
    unified_mpoly_struct *row0[3], *row1[3], **matrix[2] = {row0, row1};
    slong exponents[3] = {1,0,2};
    monom_t labels[3];
    fq_nmod_t a, one, neg;
    fq_nmod_init(a,ctx); fq_nmod_init(one,ctx); fq_nmod_init(neg,ctx);
    fq_nmod_gen(a,ctx); fq_nmod_one(one,ctx);
    for (slong i=0; i<2; i++) for (slong j=0; j<3; j++) {
        dr_mpoly_init(&entries[i][j],0,0,ctx);
        matrix[i][j]=&entries[i][j];
    }
    for (slong j=0; j<3; j++) { labels[j].exp=&exponents[j]; labels[j].idx=j; }
    dr_mpoly_add_term_fast(&entries[0][0],NULL,NULL,one);
    fq_nmod_neg(neg,a,ctx);
    dr_mpoly_add_term_fast(&entries[0][1],NULL,NULL,neg);
    dr_mpoly_add_term_fast(&entries[1][2],NULL,NULL,one);
    fq_nmod_mul(neg,a,a,ctx); fq_nmod_neg(neg,neg,ctx);
    dr_mpoly_add_term_fast(&entries[1][1],NULL,NULL,neg);
    slong *rows=NULL,*cols=NULL,size=-1;
    dixon_scalar_screen_report_t report;
    assert(dixon_scalar_screen_set_report(&report)==NULL);
    dixon_select_submatrix(matrix,2,3,NULL,labels,1,0,NULL,0,0,ctx,&rows,&cols,&size);
    assert(report.status==DIXON_SCREEN_INCONCLUSIVE);
    assert(report.constant_col==1 && report.rank==2 && report.nonconstant_rank==2);
    assert(!rows && !cols && size==0);

    /* Treat the same matrix's columns as nonconstant: no affine certificate. */
    exponents[1]=3;
    dixon_select_submatrix(matrix,2,3,NULL,labels,1,0,NULL,0,0,ctx,&rows,&cols,&size);
    assert(report.status==DIXON_SCREEN_INCONCLUSIVE && report.constant_col==-1);
    assert(!rows && !cols && size==0);
    dixon_select_submatrix(NULL,0,0,NULL,NULL,1,0,NULL,0,-1,ctx,&rows,&cols,&size);
    assert(report.status==DIXON_SCREEN_INCONCLUSIVE && report.rank==0);
    assert(dixon_scalar_screen_set_report(NULL)==&report);
    for (slong i=0; i<2; i++) for (slong j=0; j<3; j++) dr_mpoly_clear(&entries[i][j]);
    fq_nmod_clear(neg,ctx); fq_nmod_clear(one,ctx); fq_nmod_clear(a,ctx);
}

/* Compare the single-LU inference with two independent rank computations.
 * Delayed pivots and dependent rows exercise the packed L/U storage boundary. */
static void check_rank_pair(slong nr, slong nc, slong constant, int mode,
                            const fq_nmod_ctx_t ctx, flint_rand_t rng)
{
    fq_nmod_mat_t full, reduced;
    fq_nmod_mat_init(full,nr,nc,ctx);
    fq_nmod_mat_init(reduced,nr,nc-(constant>=0),ctx);
    for (slong i=0; i<nr; i++) for (slong j=0; j<nc; j++) {
        if (mode==3 || (mode==1 && j%3!=1) || (mode==4 && j!=constant)) continue;
        fq_nmod_rand(fq_nmod_mat_entry(full,i,j),rng,ctx);
    }
    if (mode==2 && constant>=0) {
        for (slong i=0; i<nr; i++) {
            fq_nmod_zero(fq_nmod_mat_entry(full,i,constant),ctx);
            for (slong j=0; j<nc; j++) if (j!=constant)
                fq_nmod_add(fq_nmod_mat_entry(full,i,constant),
                    fq_nmod_mat_entry(full,i,constant),fq_nmod_mat_entry(full,i,j),ctx);
        }
    }
    if (mode==5) for (slong i=1; i<nr; i+=2) for (slong j=0; j<nc; j++)
        fq_nmod_set(fq_nmod_mat_entry(full,i,j),fq_nmod_mat_entry(full,i-1,j),ctx);
    for (slong i=0; i<nr; i++) for (slong j=0, k=0; j<nc; j++) if (j!=constant)
        fq_nmod_set(fq_nmod_mat_entry(reduced,i,k++),fq_nmod_mat_entry(full,i,j),ctx);
    slong rank=fq_nmod_mat_rank(full,ctx), without=fq_nmod_mat_rank(reduced,ctx);

    unified_mpoly_struct ***grid=flint_calloc(nr,sizeof(*grid));
    monom_t *labels=flint_malloc(nc*sizeof(*labels));
    slong *exps=flint_malloc(nc*sizeof(*exps));
    for (slong j=0; j<nc; j++) {
        exps[j]=j==constant ? 0 : j+1; labels[j].exp=exps+j; labels[j].idx=j;
    }
    for (slong i=0; i<nr; i++) {
        grid[i]=flint_calloc(nc,sizeof(*grid[i]));
        for (slong j=0; j<nc; j++) if (!fq_nmod_is_zero(fq_nmod_mat_entry(full,i,j),ctx)) {
            grid[i][j]=flint_calloc(1,sizeof(*grid[i][j]));
            dr_mpoly_init(grid[i][j],0,0,ctx);
            dr_mpoly_add_term_fast(grid[i][j],NULL,NULL,fq_nmod_mat_entry(full,i,j));
        }
    }
    dixon_scalar_screen_report_t report;
    slong *rows=NULL,*cols=NULL,size=-1;
    assert(dixon_scalar_screen_set_report(&report)==NULL);
    dixon_select_submatrix(grid,nr,nc,NULL,labels,1,0,NULL,0,-1,ctx,&rows,&cols,&size);
    assert(report.rank==rank && report.nonconstant_rank==without);
    assert((report.status==DIXON_SCREEN_NO_COMMON_ZERO)==(rank>without));
    assert(!rows && !cols && size==0);
    assert(dixon_scalar_screen_set_report(NULL)==&report);
    for (slong i=0; i<nr; i++) {
        for (slong j=0; j<nc; j++) if (grid[i][j]) {
            dr_mpoly_clear(grid[i][j]); flint_free(grid[i][j]);
        }
        flint_free(grid[i]);
    }
    flint_free(grid); flint_free(labels); flint_free(exps);
    fq_nmod_mat_clear(reduced,ctx); fq_nmod_mat_clear(full,ctx);
}

int main(void)
{
    g_dixon_verbose_level=0;
    flint_rand_t rng;
    flint_rand_init(rng); flint_rand_set_seed(rng,1729,2718);
    for (slong extension=1; extension<=2; extension++) {
        fq_nmod_ctx_t ctx;
        fq_nmod_ctx_init_ui(ctx,7,extension,"t");
        check_matrix(ctx);
        for (slong trial=0; trial<180; trial++) {
            slong nr=n_randint(rng,13), nc=n_randint(rng,13);
            slong constant=nc && trial%4 ? n_randint(rng,nc) : -1;
            check_rank_pair(nr,nc,constant,trial%6,ctx,rng);
        }
        /* Exercise block LU as well as tiny classical elimination. */
        check_rank_pair(160,192,37,1,ctx,rng);
        check_rank_pair(192,160,159,2,ctx,rng);
        for (int recursive=0; recursive<=1; recursive++) {
            g_resultant_method=recursive ? RESULTANT_METHOD_DIXON_RECURSIVE : RESULTANT_METHOD_DIXON;
            g_dixon_fast_use_ksy_precondition=recursive;
            g_dixon_fast_ksy_constant_col=999; /* Must use labels, not this index. */
            check_system("x,y,1",ctx,DIXON_SCREEN_NO_COMMON_ZERO);
            check_system("x-1,y-2,x+y-3",ctx,DIXON_SCREEN_INCONCLUSIVE);
            check_system("x,y,x+y",ctx,DIXON_SCREEN_INCONCLUSIVE);
            check_system("x^2+1,y,x^2+1",ctx,DIXON_SCREEN_INCONCLUSIVE);
            if (extension > 1) {
                check_system("x-t,y-t^2,x*y-t^3",ctx,DIXON_SCREEN_INCONCLUSIVE);
                check_system("x-t,y,t",ctx,DIXON_SCREEN_NO_COMMON_ZERO);
            }
            /* Inconsistent but a zero Dixon polynomial: not a complete solver. */
            check_system("x,x+1,x+2",ctx,DIXON_SCREEN_INCONCLUSIVE);
            dixon_scalar_screen_report_t report;
            char *result=dixon_str_with_scalar_screen("x,y,3","x,y",ctx,&report);
            assert(result && strcmp(result,"1")==0);
            assert(report.status==DIXON_SCREEN_NO_COMMON_ZERO);
            free(result);
            result=dixon_str_with_scalar_screen("x+y+z,x-y,x+1","x,y",ctx,&report);
            assert(result && report.status==DIXON_SCREEN_NOT_RUN);
            assert(!dixon_scalar_screen_active());
            free(result);
            result=dixon_str_with_scalar_screen("x,y","x,y",ctx,&report);
            assert(!result && report.status==DIXON_SCREEN_NOT_RUN);
        }
        fq_nmod_ctx_clear(ctx);
    }
    puts("Dixon scalar screening tests passed");
    flint_rand_clear(rng);
    flint_cleanup_master();
    return 0;
}

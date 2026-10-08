/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Compare the optimized relation spaces with independent full RREFs, and
 * the certified solver's complete point sets with exhaustive enumeration. */
#define solve_by_quotient_closure test_solve_by_quotient_closure
#include "../solver/quotient_solver.c"
#undef solve_by_quotient_closure

static void check_seed(flint_rand_t state, ulong prime)
{
    for (slong n=1;n<=4;n++) for (slong degree=1;degree<=3;degree++) {
        closure_t prev,seeded;
        closure_init(&prev,n,degree,prime);
        closure_init(&seeded,n,degree+1,prime);
        slong pc=prev.basis.count, nc=seeded.basis.count, offset=nc-pc;
        nmod_mat_t initial,all;
        matrix_init(initial,pc/2,pc,prime); nmod_mat_randtest(initial,state);
        insert_rows(&prev,initial); nmod_mat_clear(initial);
        closure_seed(&seeded,&prev);
        matrix_init(all,(n+1)*prev.E->r,nc,prime);
        for (slong i=0;i<prev.E->r;i++) for (slong j=0;j<pc;j++) {
            ulong a=nmod_mat_entry(prev.E,i,j);
            nmod_mat_entry(all,i*(n+1),j+offset)=a;
            for (slong v=0;v<n;v++)
                nmod_mat_entry(all,i*(n+1)+v+1,seeded.shift[(j+offset)*n+v])=a;
        }
        insert_rows(&seeded,all);
        compact_rref(all,NULL);
        assert(nmod_mat_equal(seeded.E,all));
        nmod_mat_clear(all); closure_clear(&prev); closure_clear(&seeded);
    }
}

static void check_batches(flint_rand_t state, ulong prime)
{
    closure_t cl; closure_init(&cl,4,4,prime);
    slong nc=cl.basis.count;
    nmod_mat_t first,batch,reference;
    matrix_init(first,17,nc,prime); nmod_mat_randtest(first,state);
    matrix_init(batch,600,nc,prime); nmod_mat_randtest(batch,state);
    /* Leave free columns as well as nonconsecutive pivot columns. */
    for (slong i=0;i<batch->r;i++) for (slong j=0;j<nc;j++)
        if (j%3==0 || i>=31) nmod_mat_entry(batch,i,j)=0;
    for (slong i=0;i<first->r;i++) for (slong j=0;j<nc;j++)
        if (j%3==0) nmod_mat_entry(first,i,j)=0;
    matrix_init(reference,first->r+batch->r,nc,prime);
    nmod_mat_concat_vertical(reference,first,batch); compact_rref(reference,NULL);
    insert_rows(&cl,first); insert_rows(&cl,batch);
    assert(nmod_mat_equal(cl.E,reference));
    /* A dependent batch must leave the space unchanged. */
    insert_rows(&cl,batch); assert(nmod_mat_equal(cl.E,reference));
    nmod_mat_clear(first); nmod_mat_clear(batch); nmod_mat_clear(reference); closure_clear(&cl);
}

static void check_projected_input(flint_rand_t state, ulong prime)
{
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx,3,ORD_DEGLEX,prime);
    nmod_mpoly_struct fs[5];
    for (slong i=0;i<5;i++) nmod_mpoly_init(fs+i,ctx);
    for (slong trial=0;trial<8;trial++) {
        /* Include mixed degrees, a zero generator, and a constant generator. */
        for (slong i=0;i<3;i++) nmod_mpoly_randtest_bound(fs+i,state,12,2,ctx);
        nmod_mpoly_set_ui(fs+3,trial%2,ctx);
        closure_t sparse,dense;
        closure_init(&sparse,3,4,prime); closure_init(&dense,3,4,prime);
        slong nc=sparse.basis.count, nr=trial==0?0:(trial==1?nc:nc/3);
        nmod_mat_t initial,M;
        matrix_init(initial,nr,nc,prime); nmod_mat_randtest(initial,state);
        insert_rows(&sparse,initial); insert_rows(&dense,initial);
        insert_original_rows(&sparse,fs,5,ctx);
        original_rows(M,fs,5,ctx,&dense.basis,1,4); insert_rows(&dense,M);
        assert(nmod_mat_equal(sparse.E,dense.E));
        nmod_mat_clear(initial); nmod_mat_clear(M);
        closure_clear(&sparse); closure_clear(&dense);
    }
    for (slong i=0;i<5;i++) nmod_mpoly_clear(fs+i,ctx);
    nmod_mpoly_ctx_clear(ctx);
}

static void check_points(char **polys, slong m, slong n, ulong prime,
                         const char *dimension)
{
    const char *names[]={"x","y","z"};
    variable_info_t vars[]={{"x",0,0},{"y",0,1},{"z",0,2}};
    fq_nmod_ctx_t field; fq_nmod_ctx_init_ui(field,prime,1,"t");
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx,n,ORD_DEGLEX,prime);
    nmod_mpoly_struct fs[m];
    for (slong i=0;i<m;i++) {
        nmod_mpoly_init(fs+i,ctx);
        assert(!nmod_mpoly_set_str_pretty(fs+i,polys[i],names,ctx));
    }
    polynomial_solutions_t sols; polynomial_solutions_init(&sols,n,field);
    assert(test_solve_by_quotient_closure(polys,m,vars,n,&sols,8,0));
    assert(!sols.error_message);
    if (dimension) assert(strstr(sols.elimination_summary,dimension));
    ulong total=1; for (slong i=0;i<n;i++) total*=prime;
    unsigned char *seen=calloc(total,1);
    for (slong i=0;i<sols.num_solution_sets;i++) {
        ulong key=0,place=1;
        for (slong j=0;j<n;j++) {
            ulong a=nmod_poly_get_coeff_ui(sols.solution_sets[i][j][0],0);
            assert(a<prime); key+=place*a; place*=prime;
        }
        assert(!seen[key]); seen[key]=1;
    }
    slong count=0;
    for (ulong key=0;key<total;key++) {
        ulong point[n],value=key;
        for (slong j=0;j<n;j++) { point[j]=value%prime; value/=prime; }
        int root=1;
        for (slong i=0;i<m;i++) if (nmod_mpoly_evaluate_all_ui(fs+i,point,ctx)) root=0;
        assert(seen[key]==root); count+=root;
    }
    assert(count==sols.num_solution_sets);
    assert(sols.has_no_solutions==(count==0));
    free(seen); polynomial_solutions_clear(&sols);
    for (slong i=0;i<m;i++) nmod_mpoly_clear(fs+i,ctx);
    nmod_mpoly_ctx_clear(ctx); fq_nmod_ctx_clear(field);
}

static void check_limits(void)
{
    char *polys[]={"x*y","x*y"};
    variable_info_t vars[]={{"x",0,0},{"y",0,1}};
    fq_nmod_ctx_t field; fq_nmod_ctx_init_ui(field,7,1,"t");
    for (int memory=0;memory<2;memory++) {
        polynomial_solutions_t sols; polynomial_solutions_init(&sols,2,field);
        assert(!test_solve_by_quotient_closure(polys,2,vars,2,&sols,memory?0:4,memory?1:0));
        assert(sols.error_message && strstr(sols.error_message,memory?"memory budget reached":"degree limit reached"));
        assert(!sols.has_no_solutions && !sols.num_solution_sets);
        polynomial_solutions_clear(&sols);
    }
    fq_nmod_ctx_clear(field);
}

int main(void)
{
    flint_rand_t state; flint_rand_init(state); flint_rand_set_seed(state,123,456);
    const ulong primes[]={2,3,7,257,65537,
#if FLINT_BITS == 64
        UWORD(18446744073709551557)
#else
        UWORD(4294967291)
#endif
    };
    for (slong p=0;p<6;p++) {
        check_seed(state,primes[p]); check_batches(state,primes[p]);
        check_projected_input(state,primes[p]);
    }
    char *nilpotent[]={"x^2","y^2","x*y"};
    char *fibers[]={"x^2-x","y^2-y","x*y*(x-1)"};
    char *extension[]={"x^2+1","y-x","y^2+1"};
    char *empty[]={"x^2","x^2-1","y"};
    char *mixed[]={"x^3-x","y^2-y","x*y","0"};
    char *univariate[]={"x^3-x","0"};
    check_points(nilpotent,3,2,3,"dimension 3,");
    check_points(fibers,3,2,3,"dimension 4,");
    check_points(extension,3,2,3,"dimension 2,");
    check_points(empty,3,2,3,"dimension 0,");
    check_points(mixed,4,2,3,NULL);
    check_points(univariate,2,1,3,"dimension 3,");
    for (slong p=0;p<3;p++) for (int trial=0;trial<20;trial++) {
        char storage[4][256],*polys[4];
        const char *vars[]={"x","y","z","x*y"};
        for (int i=0;i<4;i++) {
            snprintf(storage[i],sizeof(storage[i]),"(%s)^2+%lu*x+%lu*y+%lu*z+%lu",vars[i],
                     n_randint(state,primes[p]),n_randint(state,primes[p]),
                     n_randint(state,primes[p]),trial%2?n_randint(state,primes[p]):0);
            polys[i]=storage[i];
        }
        check_points(polys,4,3,primes[p],NULL);
    }
    check_limits(); flint_rand_clear(state); flint_cleanup();
    puts("Quotient solver: relation-space and exhaustive point tests passed");
    return 0;
}

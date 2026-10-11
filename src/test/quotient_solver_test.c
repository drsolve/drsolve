/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Compare the optimized relation spaces with independent full RREFs, and
 * the certified solver's complete point sets with exhaustive enumeration. */
#define solve_by_quotient_closure test_solve_by_quotient_closure
#define quotient_sift_resultant test_quotient_sift_resultant
#include "../solver/quotient_solver.c"
#undef solve_by_quotient_closure
#undef quotient_sift_resultant

static ulong relation_entry(const relations_t E, slong row, slong col)
{
    if (E->p[row]==col) return 1;
    for (slong j=0;j<E->tail->c;j++)
        if (E->f[j]==col) return nmod_mat_entry(E->tail,row,j);
    return 0;
}

static void expand_relations(nmod_mat_t M, const relations_t E)
{
    matrix_init(M,E->r,E->c,E->mod.n);
    for (slong i=0;i<E->r;i++) for (slong j=0;j<E->c;j++)
        nmod_mat_entry(M,i,j)=relation_entry(E,i,j);
}

static int relations_equal_dense(const relations_t E, const nmod_mat_t M)
{
    nmod_mat_t dense; expand_relations(dense,E);
    int equal=nmod_mat_equal(dense,M); nmod_mat_clear(dense); return equal;
}

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
            ulong a=relation_entry(prev.E,i,j);
            nmod_mat_entry(all,i*(n+1),j+offset)=a;
            for (slong v=0;v<n;v++)
                nmod_mat_entry(all,i*(n+1)+v+1,seeded.shift[(j+offset)*n+v])=a;
        }
        insert_rows(&seeded,all);
        compact_rref(all,NULL);
        assert(relations_equal_dense(seeded.E,all));
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
    assert(relations_equal_dense(cl.E,reference));
    /* A dependent batch must leave the space unchanged. */
    insert_rows(&cl,batch); assert(relations_equal_dense(cl.E,reference));
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
        nmod_mat_t ref; expand_relations(ref,dense.E);
        assert(relations_equal_dense(sparse.E,ref)); nmod_mat_clear(ref);
        nmod_mat_clear(initial); nmod_mat_clear(M);
        closure_clear(&sparse); closure_clear(&dense);
    }
    for (slong i=0;i<5;i++) nmod_mpoly_clear(fs+i,ctx);
    nmod_mpoly_ctx_clear(ctx);
}

static void full_closure(nmod_mat_t reference, const closure_t *cl)
{
    slong nc=cl->basis.count;
    ulong prime=reference->mod.n;
    for (;;) {
        slong *p=pivots(reference), nl=0;
        for (slong i=0;i<reference->r;i++) if (cl->basis.degrees[p[i]]<cl->basis.degree) nl++;
        nmod_mat_t all; matrix_init(all,reference->r+cl->basis.n*nl,nc,prime);
        slong row=reference->r;
        for (slong i=0;i<reference->r;i++) {
            for (slong j=0;j<nc;j++) nmod_mat_entry(all,i,j)=nmod_mat_entry(reference,i,j);
            if (cl->basis.degrees[p[i]]>=cl->basis.degree) continue;
            for (slong v=0;v<cl->basis.n;v++,row++) for (slong j=0;j<nc;j++) {
                ulong a=nmod_mat_entry(reference,i,j);
                if (a) nmod_mat_entry(all,row,cl->shift[j*cl->basis.n+v])=a;
            }
        }
        flint_free(p); compact_rref(all,NULL);
        int done=all->r==reference->r;
        nmod_mat_swap(reference,all); nmod_mat_clear(all);
        if (done) break;
    }
}

static void check_seeded_input(flint_rand_t state, ulong prime)
{
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx,2,ORD_DEGLEX,prime);
    const char *names[]={"x","y"};
    nmod_mpoly_struct fs[3];
    for (slong i=0;i<3;i++) nmod_mpoly_init(fs+i,ctx);
    for (slong trial=0;trial<10;trial++) {
        if (!trial) {
            assert(!nmod_mpoly_set_str_pretty(fs,"x^2+y",names,ctx));
            assert(!nmod_mpoly_set_str_pretty(fs+1,"x^2+x",names,ctx));
            nmod_mpoly_zero(fs+2,ctx);
        } else {
            for (slong i=0;i<3;i++) nmod_mpoly_randtest_bound(fs+i,state,5,2,ctx);
        }
        closure_t prev,cl;
        closure_init(&prev,2,2,prime); closure_init(&cl,2,3,prime);
        insert_original_rows(&prev,fs,3,ctx);
        if (trial%2==0) {
            nmod_mat_t empty; matrix_init(empty,0,prev.basis.count,prime);
            closure_run(&prev,empty,ctx); nmod_mat_clear(empty);
            if (!trial) assert(!prev.original_only);
        }
        closure_seed(&cl,&prev);
        assert(cl.processed_ready==prev.original_only);
        nmod_mat_t seed,M,all;
        expand_relations(seed,cl.E);
        original_rows(M,fs,3,ctx,&cl.basis,1,3);
        matrix_init(all,seed->r+M->r,cl.basis.count,prime);
        nmod_mat_concat_vertical(all,seed,M); compact_rref(all,NULL);
        insert_original_rows(&cl,fs,3,ctx);
        assert(relations_equal_dense(cl.E,all));
        full_closure(all,&cl);
        nmod_mat_t empty; matrix_init(empty,0,cl.basis.count,prime);
        closure_run(&cl,empty,ctx); nmod_mat_clear(empty);
        assert(relations_equal_dense(cl.E,all));
        nmod_mat_clear(seed); nmod_mat_clear(M); nmod_mat_clear(all);
        closure_clear(&prev); closure_clear(&cl);
    }
    for (slong i=0;i<3;i++) nmod_mpoly_clear(fs+i,ctx);
    nmod_mpoly_ctx_clear(ctx);
}

/* The only nonzero consistency relation occurs in the second variable,
 * beyond the first 256 low monomials. Missing it would leave dimension 2. */
static void check_late_consistency(void)
{
    const ulong prime=7;
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx,2,ORD_DEGLEX,prime);
    const char *names[]={"x","y"};
    nmod_mpoly_struct fs[2];
    for (slong i=0;i<2;i++) {
        nmod_mpoly_init(fs+i,ctx);
        assert(!nmod_mpoly_set_str_pretty(fs+i,names[i],names,ctx));
    }
    closure_t cl; closure_init(&cl,2,24,prime);
    ulong x[]={1,0}, top[]={0,24};
    slong xc=basis_find(&cl.basis,x), tc=basis_find(&cl.basis,top), nc=cl.basis.count;
    nmod_mat_t M; matrix_init(M,nc-2,nc,prime);
    slong row=0;
    for (slong j=0;j<nc-1;j++) if (j!=xc) {
        nmod_mat_entry(M,row,j)=1;
        if (j==tc) nmod_mat_entry(M,row,xc)=prime-1;
        row++;
    }
    insert_rows(&cl,M); nmod_mat_clear(M);
    quotient_result result={0}; nmod_poly_init(result.eliminant,prime);
    quotient_close(&result,&cl,fs,2,ctx);
    assert(result.certified && result.dimension==1);
    assert(cl.rref_calls>=4);
    for (slong i=0;i<2;i++) {
        assert(nmod_mat_is_zero(result.operators+i));
        nmod_mat_clear(result.operators+i); nmod_mpoly_clear(fs+i,ctx);
    }
    flint_free(result.operators); nmod_poly_clear(result.eliminant);
    closure_clear(&cl); nmod_mpoly_ctx_clear(ctx);
}

/* Independent ambient closure: multiply every low-degree relation each
 * round and recompute a full RREF, without processed-space shortcuts. */
static void check_closure(flint_rand_t state, ulong prime)
{
    nmod_mpoly_ctx_t ctx; nmod_mpoly_ctx_init(ctx,3,ORD_DEGLEX,prime);
    for (slong trial=0;trial<16;trial++) {
        closure_t cl; closure_init(&cl,3,3,prime);
        slong nc=cl.basis.count;
        nmod_mat_t initial,reference;
        matrix_init(initial,4+trial,nc,prime);
        nmod_mat_randtest(initial,state);
        /* Force low pivots and changing old tails as closure grows. */
        for (slong i=0;i<initial->r;i++) for (slong j=0;j<nc;j++)
            if (cl.basis.degrees[j]>1+trial%2) nmod_mat_entry(initial,i,j)=0;
        nmod_mat_init_set(reference,initial); compact_rref(reference,NULL);
        full_closure(reference,&cl);
        closure_run(&cl,initial,ctx);
        assert(relations_equal_dense(cl.E,reference));
        nmod_mat_clear(initial); nmod_mat_clear(reference); closure_clear(&cl);
    }
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
        check_projected_input(state,primes[p]); check_closure(state,primes[p]);
        check_seeded_input(state,primes[p]);
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
    check_late_consistency(); check_limits(); flint_rand_clear(state); flint_cleanup();
    puts("Quotient solver: relation-space and exhaustive point tests passed");
    return 0;
}

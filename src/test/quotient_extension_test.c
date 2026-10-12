/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "quotient_solver.h"
#include "fq_unified_interface.h"
#include <flint/fq_nmod_mpoly.h>
#include <assert.h>

static variable_info_t vars[]={{"x",0,0},{"y",0,1}};
static const char *names[]={"x","y"};

static void check(char **polys, slong m, const fq_nmod_ctx_t field,
                  slong expected, const char *dimension, int exhaustive)
{
    polynomial_solutions_t sols; polynomial_solutions_init(&sols,2,field);
    int ok=solve_by_quotient_closure(polys,m,vars,2,&sols,8,0);
    if (!ok) fprintf(stderr,"%s\n",sols.error_message);
    assert(ok && !sols.error_message);
    if (expected>=0) assert(sols.num_solution_sets==expected);
    if (dimension) assert(strstr(sols.elimination_summary,dimension));
    fq_nmod_mpoly_ctx_t ctx; fq_nmod_mpoly_ctx_init(ctx,2,ORD_DEGLEX,field);
    fq_nmod_mpoly_struct fs[m];
    for (slong i=0;i<m;i++) {
        fq_nmod_mpoly_init(fs+i,ctx);
        assert(!fq_nmod_mpoly_set_str_pretty(fs+i,polys[i],names,ctx));
    }
    fq_nmod_t value; fq_nmod_init(value,field);
    for (slong i=0;i<sols.num_solution_sets;i++) {
        fq_nmod_struct *point[]={sols.solution_sets[i][0][0],sols.solution_sets[i][1][0]};
        for (slong j=0;j<m;j++) {
            fq_nmod_mpoly_evaluate_all_fq_nmod(value,fs+j,point,ctx);
            assert(fq_nmod_is_zero(value,field));
        }
        for (slong j=0;j<i;j++)
            assert(!fq_nmod_equal(point[0],sols.solution_sets[j][0][0],field) ||
                   !fq_nmod_equal(point[1],sols.solution_sets[j][1][0],field));
    }
    if (exhaustive) {
        ulong p=fq_nmod_ctx_prime(field),q=1;
        slong d=fq_nmod_ctx_degree(field),found=0;
        for (slong i=0;i<d;i++) q*=p;
        fq_nmod_t x,y; fq_nmod_init(x,field); fq_nmod_init(y,field);
        fq_nmod_struct *point[]={x,y};
        for (ulong a=0;a<q;a++) for (ulong b=0;b<q;b++) {
            ulong aa=a,bb=b;
            for (slong j=0;j<d;j++) {
                nmod_poly_set_coeff_ui(x,j,aa%p); aa/=p;
                nmod_poly_set_coeff_ui(y,j,bb%p); bb/=p;
            }
            int root=1;
            for (slong i=0;i<m && root;i++) {
                fq_nmod_mpoly_evaluate_all_fq_nmod(value,fs+i,point,ctx);
                root=fq_nmod_is_zero(value,field);
            }
            found+=root;
        }
        assert(found==sols.num_solution_sets);
        fq_nmod_clear(x,field); fq_nmod_clear(y,field);
    }
    fq_nmod_clear(value,field);
    for (slong i=0;i<m;i++) fq_nmod_mpoly_clear(fs+i,ctx);
    fq_nmod_mpoly_ctx_clear(ctx); polynomial_solutions_clear(&sols);
}

static void cases(const fq_nmod_ctx_t field, int exhaustive)
{
    char *linear[]={"(t+1)*x+t*y+t", "x+(t+1)*y+1", "x-t*y"};
    char *nilpotent[]={"(x-t)^2", "(y-t-1)^2", "(x-t)*(y-t-1)"};
    char *fibers[]={"(x-t)*(x-t-1)", "y*(y-1)", "0"};
    char *empty[]={"x^2", "x^2-1", "y"};
    char *constant[]={"0", "t"};
    char *closure[]={"x^2+y", "x^2+x", "0"};
    check(linear,3,field,-1,NULL,exhaustive);
    check(nilpotent,3,field,1,"dimension 3,",exhaustive);
    check(fibers,3,field,4,"dimension 4,",exhaustive);
    check(empty,3,field,0,"dimension 0,",exhaustive);
    check(constant,2,field,0,"dimension 0,",exhaustive);
    check(closure,3,field,-1,NULL,exhaustive);
    char *positive[]={"x*y", "0"};
    for (int memory=0;memory<2;memory++) {
        polynomial_solutions_t sols; polynomial_solutions_init(&sols,2,field);
        assert(!solve_by_quotient_closure(positive,2,vars,2,&sols,memory?0:4,memory?1:0));
        assert(strstr(sols.error_message,memory?"memory budget reached":"degree limit reached"));
        assert(!sols.has_no_solutions && !sols.num_solution_sets);
        polynomial_solutions_clear(&sols);
    }
}

int main(void)
{
    const slong ds[]={4,8,16,32,64,128};
    const uint64_t mods[]={0x13,0x11d,0x1002d,0x100008299,UINT64_C(0x247f43cb7),0x87};
    for (int i=0;i<6;i++) {
        nmod_poly_t mod; nmod_poly_init(mod,2); nmod_poly_set_coeff_ui(mod,ds[i],1);
        for (int j=0;j<64;j++) if ((mods[i]>>j)&1) nmod_poly_set_coeff_ui(mod,j,1);
        fq_nmod_ctx_t field; fq_nmod_ctx_init_modulus(field,mod,"t");
        field_ctx_t native; field_ctx_init(&native,field);
        assert(native.field_id==FIELD_ID_GF24+i); field_ctx_clear(&native);
        cases(field,i==0); fq_nmod_ctx_clear(field); nmod_poly_clear(mod);
    }
    /* Change GF(256)'s defining polynomial in the same process. */
    nmod_poly_t mod; nmod_poly_init(mod,2);
    nmod_poly_set_coeff_ui(mod,8,1); nmod_poly_set_coeff_ui(mod,4,1);
    nmod_poly_set_coeff_ui(mod,3,1); nmod_poly_set_coeff_ui(mod,1,1); nmod_poly_set_coeff_ui(mod,0,1);
    fq_nmod_ctx_t field; fq_nmod_ctx_init_modulus(field,mod,"t");
    cases(field,0); fq_nmod_ctx_clear(field); nmod_poly_clear(mod);
    const ulong ps[]={2,3,5,1031};
    const slong degrees[]={3,2,2,2};
    for (int i=0;i<4;i++) {
        fq_nmod_ctx_init_ui(field,ps[i],degrees[i],"t");
        field_ctx_t backend; field_ctx_init(&backend,field);
        assert(backend.field_id==FIELD_ID_FQ); /* Zech is disabled by the project policy. */
        field_ctx_clear(&backend); cases(field,i<3);
        if (i==1) {
            char *no_points[]={"x^2-t", "y"};
            check(no_points,2,field,0,"dimension 2,",1);
            flint_rand_t random; flint_rand_init(random); flint_rand_set_seed(random,17,23);
            for (int trial=0;trial<20;trial++) {
                char strings[3][128],*polys[3];
                for (int j=0;j<3;j++) {
                    snprintf(strings[j],sizeof(strings[j]),"%s+(%lu*t+%lu)*x+(%lu*t+%lu)*y+%lu*t+%lu",
                        j==0?"x^2":j==1?"y^2":"x*y",
                        n_randint(random,3),n_randint(random,3),n_randint(random,3),
                        n_randint(random,3),n_randint(random,3),n_randint(random,3));
                    polys[j]=strings[j];
                }
                check(polys,3,field,-1,NULL,1);
            }
            flint_rand_clear(random);
        }
        fq_nmod_ctx_clear(field);
    }
    flint_cleanup(); puts("Quotient extension fields: native, generic and exhaustive point tests passed");
    return 0;
}

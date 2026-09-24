/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Experimental construction from projected coefficient panels. This does not
 * build the full candidate, but it still computes complementary coefficients;
 * it is NOT an O*(4^n) quotient-algebra construction or a solver backend. */
#define main mq_filter_existing_main
#include "dixon_mq_filter_test.c"
#undef main
/* Reuse the exact blocked solve and permutation primitives in this standalone
 * experiment, without changing the production Schur backend. */
#include "../determinant/mq_poly_mat_det.c"
#include "mq_pencil_det.h"

typedef struct {
    fq_mvpoly_t **matrix;
    slong n,size,sigma;
    slong *rx,*cx,*rd,*cd;
    int pencil,fail_at;
    slong calls,terms,slots;
    double generation,packing;
} panel_source;

static int fetch_panel(nmod_poly_mat_t out,panel_source *s,const slong *rows,
                       slong nr,const slong *cols,slong nc)
{
    if(++s->calls==s->fail_at)return 0;
    double start=get_wall_time();
    slong m=s->n-1;
    slong *rx=flint_malloc(nr*m*sizeof(slong)),*cx=flint_malloc(nc*m*sizeof(slong));
    for(slong i=0;i<nr;i++)memcpy(rx+i*m,s->rx+rows[i]*m,m*sizeof(slong));
    for(slong i=0;i<nc;i++)memcpy(cx+i*m,s->cx+cols[i]*m,m*sizeof(slong));
    fq_mvpoly_t poly;mq_pencil_stats stats;
    int ok=s->pencil
        ? compute_fq_det_mq_pencil_projected(&poly,s->matrix,s->n,rx,nr,cx,nc,&stats)
        : compute_fq_det_mq_projected_rect(&poly,s->matrix,s->n,rx,nr,cx,nc);
    s->generation+=get_wall_time()-start;start=get_wall_time();
    if(ok) {
        monom_t *rm=NULL,*cm=NULL;slong rn=0,cn=0,rc=0,cc=0,rh=16,ch=16;
        hash_entry_t **ri=flint_calloc(16,sizeof(*ri)),**ci=flint_calloc(16,sizeof(*ci));
        for(slong i=0;i<nr;i++)dixon_intern_monom(&rm,&rn,&rc,&ri,&rh,rx+i*m,m);
        for(slong i=0;i<nc;i++)dixon_intern_monom(&cm,&cn,&cc,&ci,&ch,cx+i*m,m);
        nmod_poly_mat_zero(out);
        for(slong t=0;t<poly.nterms;t++) {
            fq_monomial_t *term=poly.terms+t;
            slong i=lookup_monom_index(ri,rh,term->var_exp,m);
            slong j=lookup_monom_index(ci,ch,term->var_exp+m,m);
            assert(i>=0 && j>=0);
            slong d=term->par_exp[0];
            if(d>s->sigma-s->rd[rows[i]]-s->cd[cols[j]]){ok=0;break;}
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(out,i,j),d,nmod_poly_get_coeff_ui(term->coeff,0));
        }
        s->terms+=poly.nterms;
        for(slong i=0;i<nr;i++)for(slong j=0;j<nc;j++)
            s->slots+=FLINT_MAX(0,s->sigma-s->rd[rows[i]]-s->cd[cols[j]]+1);
        free_monom_index(ri,rh);free_monom_index(ci,ch);flint_free(rm);flint_free(cm);
        fq_mvpoly_clear(&poly);
    }
    flint_free(rx);flint_free(cx);s->packing+=get_wall_time()-start;return ok;
}
static int direct_core(nmod_poly_mat_t core, ulong *factor, slong n, ulong prime,
                       const slong *rd,const slong *cd,slong h,slong sigma,
                       panel_source *source)
{
    if(h<=0 || h>=n || core->r!=h || core->c!=h)return 0;
    slong e = n - h, budget = 0;
    for (slong i = 0; i < n; i++)
        if (rd[i] < 0 || cd[i] < 0 || rd[i] > sigma || cd[i] > sigma)
            return 0;
    for (slong i = h; i < n; i++)budget += sigma-rd[i]-cd[i];
    if (budget)
        return 0;

    slong *rows = flint_malloc((size_t)e * sizeof(slong));
    slong *cols = flint_malloc((size_t)e * sizeof(slong));
    /* At most e nonempty groups. No full polynomial matrix is copied. */
    mq_constant_block *blocks = flint_malloc((size_t)e * sizeof(*blocks));
    slong nr = 0, nc = 0, groups = 0;
    int ok = 1;
    nmod_t mod;
    nmod_init(&mod, prime);
    ulong scale = 1;
    for (slong d = 0; d <= sigma && nr < e; d++)
    {
        slong begin = nr, cb = nc;
        for (slong i = h; i < n; i++)
            if (rd[i] == d)
                rows[nr++] = i;
        for (slong j = h; j < n; j++)
            if (cd[j] == sigma - d)
                cols[nc++] = j;
        if (nr - begin != nc - cb)
        {
            ok = 0;
            break;
        }
        slong k = nr - begin;
        if (!k)
            continue;
        mq_constant_block *g = blocks + groups++;
        g->begin = begin;
        g->size = k;
        g->perm = flint_malloc((size_t)k * sizeof(slong));
        nmod_mat_init(g->lu, k, k, prime);
        nmod_poly_mat_t diagonal;
        nmod_poly_mat_init(diagonal,k,k,prime);
        ok=fetch_panel(diagonal,source,rows+begin,k,cols+begin,k);
        if(ok)for(slong i=0;i<k;i++)for(slong j=0;j<k;j++)
            nmod_mat_entry(g->lu,i,j)=nmod_poly_get_coeff_ui(nmod_poly_mat_entry(diagonal,i,j),0);
        nmod_poly_mat_clear(diagonal);
        if(!ok)break;
        /* Certify every constant diagonal block before polynomial work. */
        if (nmod_mat_lu(g->perm, g->lu, 1) != k)
        {
            ok = 0;
            break;
        }
        for (slong i = 0; i < k; i++)
            scale = nmod_mul(scale, nmod_mat_entry(g->lu, i, i), mod);
        if (mq_permutation_odd(g->perm, k, 0))
            scale = nmod_neg(scale, mod);
    }
    if (nr != e || nc != e)
        ok = 0;
    if (ok)
    {
        if (mq_permutation_odd(rows, e, h) ^ mq_permutation_odd(cols, e, h))
            scale = nmod_neg(scale, mod);
        /* E X = V, solved from the highest degree block down. X has only
         * e*h entries. Coefficient panels are requested and discarded on demand. */
        nmod_poly_mat_t X,reduced;
        nmod_poly_mat_init(X,e,h,prime);nmod_poly_mat_init(reduced,h,h,prime);
        slong *retained=flint_malloc((size_t)h*sizeof(slong));
        for(slong i=0;i<h;i++)retained[i]=i;
        ok=fetch_panel(reduced,source,retained,h,retained,h);
        for (slong b = groups; ok && b-- > 0;)
        {
            mq_constant_block *g = blocks + b;
            slong start = g->begin, k = g->size, end = start + k;
            nmod_poly_mat_t rhs,panel;
            nmod_poly_mat_init(rhs,k,h,prime);
            slong width=h+e-end;
            slong *wanted=flint_malloc((size_t)width*sizeof(slong));
            for(slong j=0;j<h;j++)wanted[j]=j;
            for(slong j=end;j<e;j++)wanted[h+j-end]=cols[j];
            nmod_poly_mat_init(panel,k,width,prime);
            ok=fetch_panel(panel,source,rows+start,k,wanted,width);
            flint_free(wanted);
            if(!ok){nmod_poly_mat_clear(panel);nmod_poly_mat_clear(rhs);break;}
            if(end<e) {
                nmod_poly_mat_t off,tail;
                nmod_poly_mat_window_init(off,panel,0,h,k,width);
                nmod_poly_mat_window_init(tail,X,end,0,e,h);
                mq_block_mul(rhs,off,tail);
                nmod_poly_mat_window_clear(tail);nmod_poly_mat_window_clear(off);
            }
            slong length = 0;
            for (slong i = 0; i < k; i++)
                for (slong j = 0; j < h; j++)
                {
                    nmod_poly_struct *p = nmod_poly_mat_entry(rhs, i, j);
                    nmod_poly_sub(p, nmod_poly_mat_entry(panel,i,j), p);
                    length = FLINT_MAX(length, p->length);
                }
            if (length)
            {
                /* Reuse the same LU for all parameter coefficients. Pack
                 * them as multiple right-hand sides for blocked field solves. */
                nmod_mat_t packed, lower, solved;
                nmod_mat_init(packed, k, h * length, prime);
                nmod_mat_init(lower, k, h * length, prime);
                nmod_mat_init(solved, k, h * length, prime);
                for (slong i = 0; i < k; i++)
                    for (slong j = 0; j < h; j++)
                    {
                        const nmod_poly_struct *p = nmod_poly_mat_entry(rhs, g->perm[i], j);
                        for (slong t = 0; t < p->length; t++)
                            nmod_mat_entry(packed, i, t * h + j) = p->coeffs[t];
                    }
                nmod_mat_solve_tril(lower, g->lu, packed, 1);
                nmod_mat_solve_triu(solved, g->lu, lower, 0);
                for (slong i = 0; i < k; i++)
                    for (slong j = 0; j < h; j++)
                        for (slong t = 0; t < length; t++)
                            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(X, start + i, j), t,
                                                   nmod_mat_entry(solved, i, t * h + j));
                nmod_mat_clear(solved);
                nmod_mat_clear(lower);
                nmod_mat_clear(packed);
            }
            nmod_poly_mat_clear(rhs);nmod_poly_mat_clear(panel);
            /* Accumulate this group's U_i X_i immediately. */
            nmod_poly_mat_t u,xi,product;
            nmod_poly_mat_init(u,h,k,prime);nmod_poly_mat_init(product,h,h,prime);
            ok=fetch_panel(u,source,retained,h,cols+start,k);
            if(ok) {
                nmod_poly_mat_window_init(xi,X,start,0,end,h);
                mq_block_mul(product,u,xi);
                nmod_poly_mat_sub(reduced,reduced,product);
                nmod_poly_mat_window_clear(xi);
            }
            nmod_poly_mat_clear(product);nmod_poly_mat_clear(u);
        }
        if(ok){nmod_poly_mat_set(core,reduced);*factor=scale;}
        flint_free(retained);nmod_poly_mat_clear(reduced);nmod_poly_mat_clear(X);
    }
    for (slong b = 0; b < groups; b++)
    {
        nmod_mat_clear(blocks[b].lu);
        flint_free(blocks[b].perm);
    }
    flint_free(blocks);
    flint_free(cols);
    flint_free(rows);
    return ok;
}

/* Exactly the canonical pre-Step-1 ordering used by the production predictor. */
static void make_profile(panel_source *s,slong *h)
{
    slong m=s->n-1,*R,*H,rl,hl,sigma,rho;
    long degrees[16];for(slong i=0;i<s->n;i++)degrees[i]=2;
    assert(dixon_rank_profile_from_degrees(&R,&rl,&H,&hl,&sigma,&rho,degrees,s->n,m));
    slong count=0;for(slong i=0;i<rl;i++)count+=R[i];flint_free(R);
    slong *exps=flint_malloc(count*m*sizeof(slong)),*rev=flint_malloc(count*m*sizeof(slong));
    slong tmp[16]={0},actual=0;
    assert(dixon_mq_support(exps,count,&actual,tmp,m,0,m) && actual==count);
    monom_t *rm=NULL,*cm=NULL;slong rn=0,cn=0,rc=0,cc=0,rh=16,ch=16;
    hash_entry_t **ri=flint_calloc(16,sizeof(*ri)),**ci=flint_calloc(16,sizeof(*ci));
    for(slong excess=0;excess<=m;excess++)for(slong j=0;j<count;j++) {
        slong e=0;for(slong v=0;v<m;v++)e+=FLINT_MAX(0,exps[j*m+v]-1);
        if(e!=excess)continue;
        for(slong v=0;v<m;v++)rev[j*m+v]=exps[j*m+m-1-v];
        dixon_intern_monom(&rm,&rn,&rc,&ri,&rh,rev+j*m,m);
        dixon_intern_monom(&cm,&cn,&cc,&ci,&ch,exps+j*m,m);
    }
    slong *rows,*cols,size;
    assert(dixon_build_predicted_mirror_indices(&rows,&cols,&size,rm,rn,cm,cn,
        ri,rh,ci,ch,m,H,hl,sigma,2,rho));
    dixon_mq_step4_profile p={0};
    dixon_mq_step4_prepare(&p,rm,cm,rows,cols,size,m,degrees);
    assert(p.size==size && p.h>0);
    s->size=size;s->sigma=sigma;*h=p.h;
    s->rx=flint_malloc(size*m*sizeof(slong));s->cx=flint_malloc(size*m*sizeof(slong));
    s->rd=flint_malloc(size*sizeof(slong));s->cd=flint_malloc(size*sizeof(slong));
    for(slong i=0;i<size;i++) {
        memcpy(s->rx+i*m,rm[rows[p.rows[i]]].exp,m*sizeof(slong));
        memcpy(s->cx+i*m,cm[cols[p.cols[i]]].exp,m*sizeof(slong));
        s->rd[i]=p.rd[i];s->cd[i]=p.cd[i];
    }
    dixon_mq_step4_profile_clear(&p);flint_free(rows);flint_free(cols);flint_free(H);
    free_monom_index(ri,rh);free_monom_index(ci,ch);flint_free(rm);flint_free(cm);
    flint_free(exps);flint_free(rev);
}

#ifndef MQ_DIRECT_CORE_NO_MAIN
int main(int argc,char **argv)
{
    slong n=argc>1?atol(argv[1]):5;
    int threads=argc>2?atoi(argv[2]):1,pencil=argc>3?atoi(argv[3]):0;
    ulong seed=argc>4?strtoul(argv[4],NULL,10):12345;
    ulong prime=argc>5?strtoul(argv[5],NULL,10):65537;
    if(n<3 || n>8 || threads<1 || threads>32 || !n_is_prime(prime))return 2;
    omp_set_num_threads(threads);flint_set_num_threads(threads);
    g_dixon_verbose_level=0;g_dixon_det_cache_limit=100000;
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,seed,941);
    fq_nmod_ctx_t fq;fq_nmod_ctx_init_ui(fq,prime,1,"a");
    fq_mvpoly_t *polys=random_mq(n-1,fq,rng),**matrix,**a;
    build_fq_cancellation_matrix_mvpoly(&matrix,polys,n-1,1);
    perform_fq_matrix_row_operations_mvpoly(&a,&matrix,n-1,1);
    panel_source source={0};source.n=n;source.matrix=a;source.pencil=pencil;
    slong h;make_profile(&source,&h);slong size=source.size;
    nmod_poly_mat_t core,expected,B;
    nmod_poly_mat_init(core,h,h,prime);nmod_poly_mat_init(expected,h,h,prime);
    /* Direct path runs BEFORE the full candidate is allocated. */
    nmod_poly_one(nmod_poly_mat_entry(core,0,0));
    nmod_poly_one(nmod_poly_mat_entry(expected,0,0));ulong factor=17,ef=0;
    double start=get_wall_time();
    int ok=direct_core(core,&factor,size,prime,source.rd,source.cd,h,source.sigma,&source);
    double direct=get_wall_time()-start;
    panel_source base=source;base.calls=base.terms=base.slots=0;base.generation=base.packing=0;
    slong *all=flint_malloc(size*sizeof(slong));for(slong i=0;i<size;i++)all[i]=i;
    start=get_wall_time();nmod_poly_mat_init(B,size,size,prime);
    assert(fetch_panel(B,&base,all,size,all,size));
    double baseline_generation=get_wall_time()-start;start=get_wall_time();
    int reference=nmod_poly_mat_mq_schur(expected,&ef,B,source.rd,source.cd,h,source.sigma);
    double schur=get_wall_time()-start;
    assert(ok==reference);
    if(ok) {
        assert(factor==ef && nmod_poly_mat_equal(core,expected));
        assert(source.terms==base.terms && source.slots==base.slots);
        if(n<=6) {
            nmod_poly_t x,y;nmod_poly_init(x,prime);nmod_poly_init(y,prime);
            nmod_poly_mat_det(x,B);nmod_poly_mat_det(y,core);nmod_poly_scalar_mul_nmod(y,y,factor);
            assert(nmod_poly_equal(x,y));nmod_poly_clear(x);nmod_poly_clear(y);
        }
    } else {
        assert(factor==17 && nmod_poly_mat_equal(core,expected));
    }
    printf("{\"n\":%ld,\"q\":%lu,\"threads\":%d,\"pencil\":%d,\"seed\":%lu,\"size\":%ld,\"h\":%ld,\"ok\":%d,\"equal\":true,\"direct\":%.6f,\"generation\":%.6f,\"packing\":%.6f,\"baseline_generation\":%.6f,\"baseline_schur\":%.6f,\"panels\":%ld,\"terms\":%ld,\"baseline_terms\":%ld,\"slots\":%ld,\"baseline_slots\":%ld}\n",
      n,prime,threads,pencil,seed,size,h,ok,direct,source.generation,source.packing,baseline_generation,schur,source.calls,source.terms,base.terms,source.slots,base.slots);
    if(n<=5 && ok) {
        /* Failure after work has begun must not publish a partial core. */
        for(int fail=1;fail<=source.calls;fail++) {
            panel_source broken=source;broken.calls=0;broken.fail_at=fail;
            ulong keep=factor;
            assert(!direct_core(core,&factor,size,prime,source.rd,source.cd,h,source.sigma,&broken));
            assert(factor==keep && nmod_poly_mat_equal(core,expected));
        }
    }
    nmod_poly_mat_clear(core);nmod_poly_mat_clear(expected);nmod_poly_mat_clear(B);
    flint_free(all);flint_free(source.rx);flint_free(source.cx);flint_free(source.rd);flint_free(source.cd);
    clear_input(polys,matrix,a,n-1);fq_nmod_ctx_clear(fq);flint_rand_clear(rng);
    flint_cleanup_master();return 0;
}

#endif

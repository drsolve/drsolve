/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Default complement-first selection; --no-schur-seed disables it. The seed is nonsingular over the
 * polynomial ring: degree bounds make it block triangular with invertible
 * constant diagonal blocks. Complete it over two scalar specializations;
 * as in the existing selector, maximum rank remains a specialization heuristic.
 * All caller outputs are unchanged on rejection. */
static int dixon_schur_seed_select(unified_mpoly_struct ***matrix,
    slong nr, slong nc, const monom_t *rm, const monom_t *cm, slong m,
    const long *degrees, const fq_nmod_ctx_t ctx,
    slong **out_rows, slong **out_cols, slong *out_size)
{
    if (fq_nmod_ctx_degree(ctx) != 1 || nr<=0 || nc<=0) return 0;
    double start=get_wall_time();
    slong sigma=1, cap=FLINT_MIN(nr,nc), count=0;
    for(slong i=0;i<=m;i++)sigma+=degrees[i]-1;
    slong *rd=flint_calloc(nr,sizeof(slong)),*cd=flint_calloc(nc,sizeof(slong));
    for(slong i=0;i<nr;i++)for(slong v=0;v<m;v++)rd[i]+=rm[i].exp[v];
    for(slong j=0;j<nc;j++)for(slong v=0;v<m;v++)cd[j]+=cm[j].exp[v];
    int ok=1;
    /* Check the complete coefficient matrix, including structural zeros. */
    for(slong i=0;i<nr && ok;i++)for(slong j=0;j<nc && ok;j++) {
        unified_mpoly_struct *f=matrix[i][j];
        if(!f)continue;
        for(slong t=0;t<dr_mpoly_length(f);t++) {
            DR_MPOLY_TERM(term,f,t);
            slong d=term.par_exp ? term.par_exp[0] : 0;
            if(d>sigma-rd[i]-cd[j]){ok=0;break;}
        }
    }
    slong *ri=flint_malloc(nr*sizeof(slong)),*ci=flint_malloc(nc*sizeof(slong));
    slong *rows=flint_malloc(nr*sizeof(slong)),*cols=flint_malloc(nc*sizeof(slong));
    slong *perm=flint_malloc(nr*sizeof(slong));
    unsigned char *ur=flint_calloc(nr,1),*uc=flint_calloc(nc,1);
    ulong prime=fq_nmod_ctx_prime(ctx);
    for(slong d=0;ok && d<=sigma;d++) {
        slong a=0,b=0;
        for(slong i=0;i<nr;i++)if(rd[i]==d)ri[a++]=i;
        for(slong j=0;j<nc;j++)if(cd[j]==sigma-d)ci[b++]=j;
        if(!a || !b)continue;
        nmod_mat_t block;nmod_mat_init(block,a,b,prime);
        for(slong i=0;i<a;i++)for(slong j=0;j<b;j++) {
            unified_mpoly_struct *f=matrix[ri[i]][ci[j]];
            if(!f)continue;
            for(slong t=0;t<dr_mpoly_length(f);t++) {
                DR_MPOLY_TERM(term,f,t);
                if(!term.par_exp || term.par_exp[0]==0)
                    nmod_mat_entry(block,i,j)=nmod_add(nmod_mat_entry(block,i,j),
                        nmod_poly_get_coeff_ui(term.coeff,0),block->mod);
            }
        }
        slong rank=nmod_mat_lu(perm,block,0),next=0;
        for(slong i=0;i<rank;i++) {
            while(next<b && !nmod_mat_entry(block,i,next))next++;
            rows[count]=ri[perm[i]];cols[count]=ci[next++];
            ur[rows[count]]=uc[cols[count]]=1;count++;
        }
        nmod_mat_clear(block);
    }
    slong *best_rows=NULL,*best_cols=NULL,best=0;
    if(ok && count) {
        slong a=count,b=count;
        for(slong i=0;i<nr;i++)if(!ur[i])rows[a++]=i;
        for(slong j=0;j<nc;j++)if(!uc[j])cols[b++]=j;
        double seed_time=get_wall_time()-start;
        fq_nmod_t params[1],value;fq_nmod_init(params[0],ctx);fq_nmod_init(value,ctx);
        for(slong trial=0;trial<2;trial++) {
            fq_nmod_set_ui(params[0],(trial+1)%prime,ctx);
            nmod_mat_t full,A,B,C,D,X,S,product;
            nmod_mat_init(full,nr,nc,prime);
            for(slong i=0;i<nr;i++)for(slong j=0;j<nc;j++) {
                unified_mpoly_struct *f=matrix[rows[i]][cols[j]];
                if(f) {
                    evaluate_dr_mpoly_at_params(value,f,params);
                    nmod_mat_entry(full,i,j)=nmod_poly_get_coeff_ui(value,0);
                }
            }
            nmod_mat_window_init(A,full,0,0,count,count);
            nmod_mat_window_init(B,full,0,count,count,nc);
            nmod_mat_window_init(C,full,count,0,nr,count);
            nmod_mat_window_init(D,full,count,count,nr,nc);
            nmod_mat_init(X,count,nc-count,prime);
            nmod_mat_init(S,nr-count,nc-count,prime);
            nmod_mat_init(product,nr-count,nc-count,prime);
            int solved=nmod_mat_solve(X,A,B);
            if(solved) {
                nmod_mat_mul(product,C,X);nmod_mat_sub(S,D,product);
                slong rank=nmod_mat_lu(perm,S,0),next=0;
                if(count+rank>best) {
                    best=count+rank;
                    if(!best_rows) {best_rows=flint_malloc(cap*sizeof(slong));best_cols=flint_malloc(cap*sizeof(slong));}
                    for(slong i=0;i<count;i++){best_rows[i]=rows[i];best_cols[i]=cols[i];}
                    for(slong i=0;i<rank;i++) {
                        while(next<S->c && !nmod_mat_entry(S,i,next))next++;
                        best_rows[count+i]=rows[count+perm[i]];best_cols[count+i]=cols[count+next++];
                    }
                }
            } else ok=0;
            nmod_mat_clear(product);nmod_mat_clear(S);nmod_mat_clear(X);
            nmod_mat_window_clear(D);nmod_mat_window_clear(C);nmod_mat_window_clear(B);nmod_mat_window_clear(A);
            nmod_mat_clear(full);
            if(!ok)break;
        }
        fq_nmod_clear(value,ctx);fq_nmod_clear(params[0],ctx);
        if(ok && best) {
            *out_rows=best_rows;*out_cols=best_cols;*out_size=best;
            best_rows=best_cols=NULL;
            dixon_info_log("  Schur-first selection: seed=%ld, selected=%ld, retained<=%ld; seed %.3fs, completion %.3fs\n",
                count,best,best-count,seed_time,get_wall_time()-start-seed_time);
        }
    }
    flint_free(best_rows);flint_free(best_cols);
    flint_free(ur);flint_free(uc);flint_free(perm);flint_free(ri);flint_free(ci);
    flint_free(rows);flint_free(cols);flint_free(rd);flint_free(cd);
    return ok && best>0;
}

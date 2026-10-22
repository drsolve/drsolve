/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Default complement-first selection; --no-schur-seed disables it. The seed is nonsingular over the
 * polynomial ring: degree bounds make it block triangular with invertible
 * constant diagonal blocks. Complete it over two scalar specializations;
 * as in the existing selector, maximum rank remains a specialization heuristic.
 * All caller outputs are unchanged on rejection. */
static int SCHUR_SEED_NAME(unified_mpoly_struct ***matrix,
    slong nr, slong nc, const monom_t *rm, const monom_t *cm, slong m,
    const long *degrees, const fq_nmod_ctx_t ctx,
    slong **out_rows, slong **out_cols, slong *out_size
#ifdef SCHUR_SEED_EXTRA
    SCHUR_SEED_EXTRA
#endif
    )
{
    if (nr<=0 || nc<=0) return 0;
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
#ifndef SCHUR_SEED_EXTENSION
    ulong prime=fq_nmod_ctx_prime(ctx);
#endif
    for(slong d=0;ok && d<=sigma;d++) {
        slong a=0,b=0;
        for(slong i=0;i<nr;i++)if(rd[i]==d)ri[a++]=i;
        for(slong j=0;j<nc;j++)if(cd[j]==sigma-d)ci[b++]=j;
        if(!a || !b)continue;
        SEED_MAT block;SEED_INIT(block,a,b);
        for(slong i=0;i<a;i++)for(slong j=0;j<b;j++) {
            unified_mpoly_struct *f=matrix[ri[i]][ci[j]];
            if(!f)continue;
            for(slong t=0;t<dr_mpoly_length(f);t++) {
                DR_MPOLY_TERM(term,f,t);
                if(!term.par_exp || term.par_exp[0]==0)
                    SEED_ACCUMULATE(block,i,j,term.coeff);
            }
        }
        slong rank=SEED_CALL(lu,perm,block,0),next=0;
        for(slong i=0;i<rank;i++) {
            while(next<b && SEED_ZERO(block,i,next))next++;
            rows[count]=ri[perm[i]];cols[count]=ci[next++];
            ur[rows[count]]=uc[cols[count]]=1;count++;
        }
        SEED_CALL(clear,block);
    }
    slong *best_rows=NULL,*best_cols=NULL,best=0;
    if(ok && count) {
        slong a=count,b=count;
        for(slong i=0;i<nr;i++)if(!ur[i])rows[a++]=i;
        for(slong j=0;j<nc;j++)if(!uc[j])cols[b++]=j;
        double seed_time=get_wall_time()-start;
        fq_nmod_t params[1],value;fq_nmod_init(params[0],ctx);fq_nmod_init(value,ctx);
        for(slong trial=0;trial<2;trial++) {
#ifdef SCHUR_SEED_EXTENSION
            /* Sample in the extension itself, not just the prime subfield. */
            fq_nmod_gen(params[0],ctx);
            if (trial) { fq_nmod_one(value,ctx); fq_nmod_add(params[0],params[0],value,ctx); }
#else
            fq_nmod_set_ui(params[0],(trial+1)%prime,ctx);
#endif
#ifdef SEED_NATIVE_EVAL
            SEED_E *powers=flint_malloc((sigma+1)*sizeof(*powers));
            powers[0]=SEED_FROM_WORD(1);
            SEED_E point=SEED_FROM_FQ(params[0]);
            for(slong d=1;d<=sigma;d++) powers[d]=SEED_MUL(powers[d-1],point);
#endif
            SEED_MAT full,A,B,C,D,X,S,product;
            SEED_INIT(full,nr,nc);
            for(slong i=0;i<nr;i++)for(slong j=0;j<nc;j++) {
                unified_mpoly_struct *f=matrix[rows[i]][cols[j]];
                if(f) {
#ifdef SEED_NATIVE_EVAL
                    SEED_E v=SEED_FROM_WORD(0);
                    for(slong t=0;t<dr_mpoly_length(f);t++) {
                        DR_MPOLY_TERM(term,f,t);
                        slong d=term.par_exp ? term.par_exp[0] : 0;
                        v=SEED_ADD(v,SEED_MUL(SEED_FROM_FQ(term.coeff),powers[d]));
                    }
                    SEED_ENTRY(full,i,j)=v;
#else
                    evaluate_dr_mpoly_at_params(value,f,params);
                    SEED_STORE(full,i,j,value);
#endif
                }
            }
            SEED_CALL(window_init,A,full,0,0,count,count);
            SEED_CALL(window_init,B,full,0,count,count,nc);
            SEED_CALL(window_init,C,full,count,0,nr,count);
            SEED_CALL(window_init,D,full,count,count,nr,nc);
            SEED_INIT(X,count,nc-count);
            SEED_INIT(S,nr-count,nc-count);
            SEED_INIT(product,nr-count,nc-count);
            int solved=SEED_CALL(solve,X,A,B);
            if(solved) {
                SEED_CALL(mul,product,C,X);SEED_CALL(sub,S,D,product);
                slong rank=SEED_CALL(lu,perm,S,0),next=0;
                if(count+rank>best) {
                    best=count+rank;
                    if(!best_rows) {best_rows=flint_malloc(cap*sizeof(slong));best_cols=flint_malloc(cap*sizeof(slong));}
                    for(slong i=0;i<count;i++){best_rows[i]=rows[i];best_cols[i]=cols[i];}
                    for(slong i=0;i<rank;i++) {
                        while(next<S->c && SEED_ZERO(S,i,next))next++;
                        best_rows[count+i]=rows[count+perm[i]];best_cols[count+i]=cols[count+next++];
                    }
                }
            } else ok=0;
            SEED_CALL(clear,product);SEED_CALL(clear,S);SEED_CALL(clear,X);
            SEED_CALL(window_clear,D);SEED_CALL(window_clear,C);SEED_CALL(window_clear,B);SEED_CALL(window_clear,A);
            SEED_CALL(clear,full);
#ifdef SEED_NATIVE_EVAL
            flint_free(powers);
#endif
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

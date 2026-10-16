/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Experimental specialization of coefficient_core: lazy root lookup and
 * optional distribution of U contraction into root edges. */
static int dag_core(nmod_poly_mat_t core,ulong *factor,const coefficient_source *a,
    slong h,ulong prime,double *multiply,double *solve)
{
    slong n=a->n,e=n-h,groups=0,nr=0,nc=0;
    slong *rows=flint_malloc(e*sizeof(slong)),*cols=flint_malloc(e*sizeof(slong));
    mq_constant_block *g=flint_malloc(e*sizeof(*g));
    nmod_mat_struct *x=flint_calloc(e,sizeof(*x));
    nmod_t mod;nmod_init(&mod,prime);ulong scale=1;int ok=1;
    slong minr=WORD_MAX,minc=WORD_MAX;
    for(slong i=0;i<h;i++){minr=FLINT_MIN(minr,a->rd[i]);minc=FLINT_MIN(minc,a->cd[i]);}
    for(slong d=0;d<=a->sigma;d++) {
        slong begin=nr,cb=nc;
        for(slong i=h;i<n;i++)if(a->rd[i]==d)rows[nr++]=i;
        for(slong j=h;j<n;j++)if(a->cd[j]==a->sigma-d)cols[nc++]=j;
        if(nr-begin!=nc-cb){ok=0;break;}
        slong k=nr-begin;if(!k)continue;
        mq_constant_block *b=g+groups++;b->begin=begin;b->size=k;
        b->perm=flint_malloc(k*sizeof(slong));nmod_mat_init(b->lu,k,k,prime);
        for(slong i=0;i<k;i++)for(slong j=0;j<k;j++)
            nmod_mat_entry(b->lu,i,j)=dag_coefficient(a,rows[begin+i],cols[begin+j],0);
        if(nmod_mat_lu(b->perm,b->lu,1)!=k){ok=0;break;}
        for(slong i=0;i<k;i++)scale=nmod_mul(scale,nmod_mat_entry(b->lu,i,i),mod);
        if(mq_permutation_odd(b->perm,k,0))scale=nmod_neg(scale,mod);
    }
    if(nr!=e || nc!=e)ok=0;
    slong length=a->sigma-minr-minc+1;
    nmod_mat_t reduced;nmod_mat_init(reduced,h,h*length,prime);
    if(ok) {
        if(mq_permutation_odd(rows,e,h)^mq_permutation_odd(cols,e,h))scale=nmod_neg(scale,mod);
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)for(slong d=0;d<length;d++)
            nmod_mat_entry(reduced,i,d*h+j)=dag_coefficient(a,i,j,d);
        for(slong b=groups;b-- >0;) {
            slong start=g[b].begin,k=g[b].size,degree=a->rd[rows[start]];
            slong len=a->sigma-degree-minc+1;
            nmod_mat_t rhs,perm,lower;nmod_mat_init(rhs,k,h*len,prime);
            nmod_mat_init(perm,k,h*len,prime);nmod_mat_init(lower,k,h*len,prime);
            nmod_mat_init(x+b,k,h*len,prime);
            for(slong i=0;i<k;i++)for(slong j=0;j<h;j++)for(slong d=0;d<len;d++)
                nmod_mat_entry(rhs,i,d*h+j)=dag_coefficient(a,rows[start+i],j,d);
            double time=get_wall_time();
            for(slong c=b+1;c<groups;c++) {
                slong cs=g[c].begin,ck=g[c].size;
                slong bound=a->sigma-degree-a->cd[cols[cs]];
                nmod_mat_t entry,view;nmod_mat_init(entry,k,ck,prime);
                for(slong d=0;d<=bound;d++) {
                    for(slong i=0;i<k;i++)for(slong j=0;j<ck;j++)
                        nmod_mat_entry(entry,i,j)=dag_coefficient(a,rows[start+i],cols[cs+j],d);
                    nmod_mat_window_init(view,rhs,0,d*h,k,d*h+x[c].c);
                    nmod_mat_submul(view,view,entry,x+c);
                    nmod_mat_window_clear(view);
                }
                nmod_mat_clear(entry);
            }
            *multiply+=get_wall_time()-time;time=get_wall_time();
            for(slong i=0;i<k;i++)for(slong j=0;j<h*len;j++)
                nmod_mat_entry(perm,i,j)=nmod_mat_entry(rhs,g[b].perm[i],j);
            nmod_mat_solve_tril(lower,g[b].lu,perm,1);
            nmod_mat_solve_triu(x+b,g[b].lu,lower,0);
            *solve+=get_wall_time()-time;
            nmod_mat_clear(lower);nmod_mat_clear(perm);nmod_mat_clear(rhs);
            time=get_wall_time();
            slong ubound=a->sigma-minr-a->cd[cols[start]];
            if(dag.mode>=2)dag_contract(reduced,a,x+b,cols,start,k,h,ubound);
            else {
            nmod_mat_t u,view;nmod_mat_init(u,h,k,prime);
            for(slong d=0;d<=ubound;d++) {
                for(slong i=0;i<h;i++)for(slong j=0;j<k;j++)
                    nmod_mat_entry(u,i,j)=dag_coefficient(a,i,cols[start+j],d);
                nmod_mat_window_init(view,reduced,0,d*h,h,d*h+x[b].c);
                nmod_mat_submul(view,view,u,x+b);
                nmod_mat_window_clear(view);
            }
            nmod_mat_clear(u);
            }
            *multiply+=get_wall_time()-time;
        }
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)for(slong d=0;d<length;d++)
            nmod_poly_set_coeff_ui(nmod_poly_mat_entry(core,i,j),d,nmod_mat_entry(reduced,i,d*h+j));
        *factor=scale;
    }
    nmod_mat_clear(reduced);
    for(slong b=0;b<groups;b++){nmod_mat_clear(g[b].lu);flint_free(g[b].perm);if(x[b].r)nmod_mat_clear(x+b);}
    flint_free(x);flint_free(g);flint_free(rows);flint_free(cols);return ok;
}


/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Included by mq_poly_mat_det.c after the shared constant-block helpers. */
#include <stdint.h>
#include <string.h>
int mq_coeff_mat_init(mq_coeff_mat *out,slong n,const slong *rd,const slong *cd,
                      slong sigma,ulong prime)
{
    if(!out || !rd || !cd || n<=0 || sigma<0 || sigma>=WORD_MAX || prime<2 ||
       (size_t)n>SIZE_MAX/sizeof(size_t)/(size_t)n || n>WORD_MAX/(sigma+1))return 0;
    for(slong i=0;i<n;i++)if(rd[i]<0 || cd[i]<0 || rd[i]>sigma || cd[i]>sigma)return 0;
    mq_coeff_mat a={.n=n,.sigma=sigma,.prime=prime,.rd=rd,.cd=cd};
    a.offset=flint_malloc((size_t)n*n*sizeof(size_t));
    for(slong i=0;i<n;i++)for(slong j=0;j<n;j++) {
        size_t length=FLINT_MAX(0,sigma-rd[i]-cd[j]+1);
        if(length>SIZE_MAX/sizeof(ulong)-a.slots){flint_free(a.offset);return 0;}
        a.offset[(size_t)i*n+j]=a.slots;a.slots+=length;
    }
    a.values=flint_calloc(FLINT_MAX(a.slots,1),sizeof(ulong));*out=a;return 1;
}
void mq_coeff_mat_clear(mq_coeff_mat *a)
{
    flint_free(a->offset);flint_free(a->values);memset(a,0,sizeof(*a));
}
static ulong mq_coeff_get(const mq_coeff_mat *a,slong r,slong c,slong d)
{
    return d>=0 && d<=a->sigma-a->rd[r]-a->cd[c]
        ? a->values[a->offset[(size_t)r*a->n+c]+d] : 0;
}
/* Each panel has one writer; avoid nested parallelism if the caller already
 * configured FLINT's own thread pool or entered an OpenMP region. */
static void mq_coeff_submul(nmod_mat_t dest,const nmod_mat_t left,const nmod_mat_t right)
{
    slong panels=1;
#ifdef _OPENMP
    if(!omp_in_parallel() && flint_get_num_threads()==1)
        panels=FLINT_MIN(omp_get_max_threads(),FLINT_MAX(1,left->r/32));
    #pragma omp parallel for schedule(static) if(panels>1) num_threads(panels)
#endif
    for(slong p=0;p<panels;p++) {
        slong first=left->r*p/panels,last=left->r*(p+1)/panels;
        nmod_mat_t d,l;
        nmod_mat_window_init(d,dest,first,0,last,dest->c);
        nmod_mat_window_init(l,left,first,0,last,left->c);
        nmod_mat_submul(d,d,l,right);
        nmod_mat_window_clear(l);nmod_mat_window_clear(d);
    }
}
/* E_ij(t) X_j(t): multiply each coefficient of E_ij by all coefficients
 * of X_j in one packed matrix, then shift/subtract into the right-hand side.
 * Constant diagonal LU is reused across all right-hand-side coefficients.
 * U_i X_i is accumulated immediately into the core. */
int nmod_poly_mat_mq_schur_coeffs(nmod_poly_mat_t core,ulong *factor,const mq_coeff_mat *a,
    slong h)
{
    if (!a || !factor || !a->rd || !a->cd || !a->values || !a->offset ||
        h<=0 || h>=a->n || core->r!=h || core->c!=h || core->modulus!=a->prime ||
        a->sigma>=WORD_MAX || h>WORD_MAX/(a->sigma+1)) return 0;
    slong n=a->n,e=n-h,groups=0,nr=0,nc=0,budget=0;
    ulong prime=a->prime;
    slong minr=WORD_MAX,minc=WORD_MAX;
    for(slong i=0;i<h;i++){minr=FLINT_MIN(minr,a->rd[i]);minc=FLINT_MIN(minc,a->cd[i]);}
    if(minr>a->sigma-minc || (size_t)n>SIZE_MAX/sizeof(ulong)/(h*(a->sigma+1)))return 0;
    for (slong i=h;i<n;i++) budget+=a->sigma-a->rd[i]-a->cd[i];
    if (budget) return 0;
    slong *rows=flint_malloc(e*sizeof(slong)),*cols=flint_malloc(e*sizeof(slong));
    mq_constant_block *g=flint_malloc(e*sizeof(*g));
    nmod_mat_struct *x=flint_calloc(e,sizeof(*x));
    nmod_t mod;nmod_init(&mod,prime);ulong scale=1;int ok=1;
    for(slong d=0;d<=a->sigma;d++) {
        slong begin=nr,cb=nc;
        for(slong i=h;i<n;i++)if(a->rd[i]==d)rows[nr++]=i;
        for(slong j=h;j<n;j++)if(a->cd[j]==a->sigma-d)cols[nc++]=j;
        if(nr-begin!=nc-cb){ok=0;break;}
        slong k=nr-begin;if(!k)continue;
        mq_constant_block *b=g+groups++;b->begin=begin;b->size=k;
        b->perm=flint_malloc(k*sizeof(slong));nmod_mat_init(b->lu,k,k,prime);
        for(slong i=0;i<k;i++)for(slong j=0;j<k;j++)
            nmod_mat_entry(b->lu,i,j)=mq_coeff_get(a,rows[begin+i],cols[begin+j],0);
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
            nmod_mat_entry(reduced,i,d*h+j)=mq_coeff_get(a,i,j,d);
        for(slong b=groups;b-- >0;) {
            slong start=g[b].begin,k=g[b].size,degree=a->rd[rows[start]];
            slong len=a->sigma-degree-minc+1;
            nmod_mat_t rhs,perm,lower;nmod_mat_init(rhs,k,h*len,prime);
            nmod_mat_init(perm,k,h*len,prime);nmod_mat_init(lower,k,h*len,prime);
            nmod_mat_init(x+b,k,h*len,prime);
            for(slong i=0;i<k;i++)for(slong j=0;j<h;j++)for(slong d=0;d<len;d++)
                nmod_mat_entry(rhs,i,d*h+j)=mq_coeff_get(a,rows[start+i],j,d);
            for(slong c=b+1;c<groups;c++) {
                slong cs=g[c].begin,ck=g[c].size;
                slong bound=a->sigma-degree-a->cd[cols[cs]];
                nmod_mat_t entry,view;nmod_mat_init(entry,k,ck,prime);
                for(slong d=0;d<=bound;d++) {
                    for(slong i=0;i<k;i++)for(slong j=0;j<ck;j++)
                        nmod_mat_entry(entry,i,j)=mq_coeff_get(a,rows[start+i],cols[cs+j],d);
                    nmod_mat_window_init(view,rhs,0,d*h,k,d*h+x[c].c);
                    mq_coeff_submul(view,entry,x+c);
                    nmod_mat_window_clear(view);
                }
                nmod_mat_clear(entry);
            }
            for(slong i=0;i<k;i++)for(slong j=0;j<h*len;j++)
                nmod_mat_entry(perm,i,j)=nmod_mat_entry(rhs,g[b].perm[i],j);
            nmod_mat_solve_tril(lower,g[b].lu,perm,1);
            nmod_mat_solve_triu(x+b,g[b].lu,lower,0);
            nmod_mat_clear(lower);nmod_mat_clear(perm);nmod_mat_clear(rhs);
            slong ubound=a->sigma-minr-a->cd[cols[start]];
            nmod_mat_t u,view;nmod_mat_init(u,h,k,prime);
            for(slong d=0;d<=ubound;d++) {
                for(slong i=0;i<h;i++)for(slong j=0;j<k;j++)
                    nmod_mat_entry(u,i,j)=mq_coeff_get(a,i,cols[start+j],d);
                nmod_mat_window_init(view,reduced,0,d*h,h,d*h+x[b].c);
                mq_coeff_submul(view,u,x+b);
                nmod_mat_window_clear(view);
            }
            nmod_mat_clear(u);
        }
        /* Check the structural degree certificate before publishing anything. */
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)
            for(slong d=FLINT_MAX(0,a->sigma-a->rd[i]-a->cd[j]+1);d<length;d++)
                if(nmod_mat_entry(reduced,i,d*h+j))ok=0;
        if(ok) {
            nmod_poly_mat_zero(core);
            for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)for(slong d=0;d<length;d++)
                nmod_poly_set_coeff_ui(nmod_poly_mat_entry(core,i,j),d,nmod_mat_entry(reduced,i,d*h+j));
            *factor=scale;
        }
    }
    nmod_mat_clear(reduced);
    for(slong b=0;b<groups;b++){nmod_mat_clear(g[b].lu);flint_free(g[b].perm);if(x[b].r)nmod_mat_clear(x+b);}
    flint_free(x);flint_free(g);flint_free(rows);flint_free(cols);return ok;
}


/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Constant leading-x Macaulay frames. A relation of x-degree w has only
 * w-1 and w-2 lower parts. Thus inverse pivot blocks are parameter-free;
 * successive normal forms have polynomial, degree-bounded coefficients.
 * This experiment admits only the expected truncated Hilbert ranks. */
typedef struct {
    grade support;
    slong sx,rank,nfree,lower_s,*row,*pivot,*free,*free_index;
    dq_mat_t inverse,tail;
    slong *lower_index;
} leading_frame;
static leading_frame *frames;
static int frames_enabled;

static int leading_init(const dq_coeff *q)
{
    frames=flint_calloc(L+1,sizeof(leading_frame));
    for(slong w=0;w<=L;w++) {
        leading_frame *f=frames+w;grade_support(&f->support,w);
        f->sx=choose(N+w-2,w);f->nfree=choose(N,w)-choose(N,w-1);
        f->rank=f->sx-f->nfree;
        f->free_index=flint_malloc(f->sx*sizeof(slong));
        for(slong i=0;i<f->sx;i++)f->free_index[i]=-1;
        f->free=flint_malloc(FLINT_MAX(1,f->nfree)*sizeof(slong));
        if(w<2){for(slong i=0;i<f->sx;i++){f->free[i]=i;f->free_index[i]=i;}continue;}
        f->lower_s=frames[w-2].sx;
        dq_mat_t relations;
        dq_mat_init(relations,N*f->lower_s,f->sx,prime);
        for(slong eq=0;eq<N;eq++)for(slong a=0;a<f->lower_s;a++)
            for(slong i=0;i<N-1;i++)for(slong j=i;j<N-1;j++) {
                uint64_t key=frames[w-2].support.keys[a]+(UINT64_C(1)<<(key_bits*i))+(UINT64_C(1)<<(key_bits*j));
                dq_set(dq_mat_entry(relations,eq*f->lower_s+a,index_key(&f->support,key)),q[(eq*N+i)*N+j]);
            }
        /* One rank-revealing LU supplies the selected original rows,
         * pivot columns and both factors of their pivot block B = LU.
         * FLINT packs L below the diagonal in the first rank columns;
         * U keeps its original (possibly nonconsecutive) pivot columns. */
        slong *perm=flint_malloc(relations->r*sizeof(slong));
        slong actual=dq_mat_lu(perm,relations,0);
        if(actual!=f->rank) {
            flint_free(perm);dq_mat_clear(relations);return 0;
        }
        f->row=flint_malloc(f->rank*sizeof(slong));
        f->pivot=flint_malloc(f->rank*sizeof(slong));
        for(slong i=0;i<f->rank;i++)f->row[i]=perm[i];
        flint_free(perm);
        slong col=0,nf=0;
        for(slong r=0;r<f->rank;r++) {
            /* col >= r, so this scan never sees packed L entries. */
            while(col<f->sx && !dq_nonzero(dq_mat_entry(relations,r,col))) {
                assert(nf<f->nfree);
                f->free_index[col]=nf;f->free[nf++]=col++;
            }
            assert(col<f->sx);f->pivot[r]=col++;
        }
        while(col<f->sx) {
            assert(nf<f->nfree);
            f->free_index[col]=nf;f->free[nf++]=col++;
        }
        assert(nf==f->nfree);
        dq_mat_t lower,upper,identity;
        dq_mat_init(lower,f->rank,f->rank,prime);
        dq_mat_init(upper,f->rank,f->rank,prime);
        dq_mat_init(identity,f->rank,f->rank,prime);dq_mat_one(identity);
        dq_mat_init(f->inverse,f->rank,f->rank,prime);
        dq_mat_init(f->tail,f->rank,f->nfree,prime);
        for(slong i=0;i<f->rank;i++) {
            dq_one(dq_mat_entry(lower,i,i));
            for(slong j=0;j<i;j++)dq_set(dq_mat_entry(lower,i,j),dq_mat_entry(relations,i,j));
            for(slong j=i;j<f->rank;j++)dq_set(dq_mat_entry(upper,i,j),dq_mat_entry(relations,i,f->pivot[j]));
            /* U has zero entries before its row's pivot. */
            for(slong j=0;j<f->nfree;j++)if(f->free[j]>f->pivot[i])
                dq_set(dq_mat_entry(f->tail,i,j),dq_mat_entry(relations,i,f->free[j]));
        }
        dq_mat_clear(relations);
        dq_mat_solve_tril(f->inverse,lower,identity,1);
        dq_mat_solve_triu(f->inverse,upper,f->inverse,0);
        /* T = L U_free, hence -B^{-1}T = -U^{-1}U_free. */
        dq_mat_solve_triu(f->tail,upper,f->tail,0);
        dq_mat_neg(f->tail,f->tail);
        dq_mat_clear(lower);dq_mat_clear(upper);dq_mat_clear(identity);
        f->lower_index=flint_malloc(f->rank*N*sizeof(slong));
        for(slong r=0;r<f->rank;r++) {
            uint64_t alpha=frames[w-2].support.keys[f->row[r]%f->lower_s];
            for(slong i=0;i<N;i++)
                f->lower_index[r*N+i]=index_key(&frames[w-1].support,
                    alpha+(UINT64_C(1)<<(key_bits*i)));
        }
    }
    return 1;
}
static void leading_clear(void)
{
    if(!frames)return;
    for(slong w=0;w<=L;w++) {
        leading_frame *f=frames+w;if(!f->support.keys)continue;
        flint_free(f->support.keys);flint_free(f->support.hash);flint_free(f->support.free);dq_mat_clear(f->support.reduction);
        flint_free(f->lower_index);flint_free(f->row);flint_free(f->pivot);flint_free(f->free);flint_free(f->free_index);
        if(f->inverse->r)dq_mat_clear(f->inverse);
        if(f->tail->r)dq_mat_clear(f->tail);
        memset(f,0,sizeof(*f));
    }
    flint_free(frames);frames=NULL;
}
/* Monomials are ordered by increasing z exponent. Consequently the rows
 * after the pure-x block, and the basis after its pure-x free block, are
 * exactly z times the previous grade in the same order. */
static void cached_normal_forms(grade *g,const dq_coeff *q)
{
    for(slong d=0;d<=L;d++) {
        const leading_frame *f=frames+d;
        grade_support(g+d,d);
        for(slong j=0;j<f->nfree;j++) {
            g[d].free[j]=f->free[j];
            dq_one(dq_mat_entry(g[d].reduction,f->free[j],j));
        }
        if(!d)continue;
        assert(g[d].s==f->sx+g[d-1].s);
        assert(g[d].h==f->nfree+g[d-1].h);
        for(slong j=0;j<g[d-1].h;j++)
            g[d].free[f->nfree+j]=f->sx+g[d-1].free[j];
        for(slong i=0;i<g[d-1].s;i++)
            dq_vec_set(dq_mat_entry_ptr(g[d].reduction,f->sx+i,f->nfree),
                   dq_mat_entry_ptr(g[d-1].reduction,i,0),g[d-1].h);
        if(!f->rank)continue;
        dq_mat_t rhs,sol;
        dq_mat_init(rhs,f->rank,g[d-1].h,prime);
        dq_mat_init(sol,f->rank,g[d-1].h,prime);
        for(slong r=0;r<f->rank;r++) {
            slong eq=f->row[r]/f->lower_s;
            for(slong i=0;i<N;i++) {
                const dq_coeff *scalar=q+(eq*N+i)*N+N-1;if(!dq_nonzero(*scalar))continue;
                slong c=f->lower_index[r*N+i];
                dq_vec_addmul(dq_mat_entry_ptr(rhs,r,0),
                    dq_mat_entry_ptr(g[d-1].reduction,c,0),g[d-1].h,*scalar);
            }
        }
        dq_mat_mul(sol,f->inverse,rhs);
        for(slong r=0;r<f->rank;r++) {
            slong c=f->pivot[r];
            for(slong j=0;j<f->nfree;j++)
                dq_set(dq_mat_entry(g[d].reduction,c,j),dq_mat_entry(f->tail,r,j));
            for(slong j=0;j<g[d-1].h;j++)
                dq_neg(dq_mat_entry(g[d].reduction,c,f->nfree+j),dq_mat_entry(sol,r,j));
        }
        dq_mat_clear(rhs);dq_mat_clear(sol);
    }
}
/* Basis extraction needs only the shared pivot/free-column choices. */
static void leading_basis(uint64_t *basis,slong degree)
{
    slong count=0;
    for(slong z=0;z<=degree;z++) {
        const leading_frame *f=frames+degree-z;
        for(slong j=0;j<f->nfree;j++)
            basis[count++]=f->support.keys[f->free[j]]+
                ((uint64_t)z<<(key_bits*(N-1)));
    }
    assert(count==choose(N,degree));
}

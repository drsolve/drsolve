/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Accumulate the complete Laplace update before modular reduction.
 * Canonical residues and an explicit bound prevent word overflow. Large
 * characteristics use the original FLINT recurrence instead. */
static int fused_minor(nmod_mat_t dest,const nmod_mat_struct *state,
                       const linear_map *op,slong mask,slong k,slong a,slong b)
{
    slong inner=0;
    if(a)inner+=choose(N,a-1);
    if(b)inner+=choose(N,b-1);
    /* Divide the bound rather than multiply it: this is safe even for
     * characteristics close to the machine-word limit. */
    if((ulong)(k*inner)>UWORD_MAX/(prime-1)/(prime-1))return 0;
    slong rows=dest->r,cols=dest->c;
    ulong *acc=flint_calloc(rows*cols,sizeof(ulong));
    slong pos=0;
    for(slong j=0;j<N;j++)if(mask&(1L<<j)) {
        slong child=mask^(1L<<j);int neg=(k-1+pos++)&1;
        if(a) {
            const nmod_mat_struct *src=state+child*(M+1)+a-1;
            const nmod_mat_struct *map=op[((k-1)*N+j)*L+a-1].x;
            for(slong i=0;i<rows;i++)for(slong r=0;r<src->r;r++) {
                ulong c=nmod_mat_entry(map,i,r);if(!c)continue;
                if(neg)c=prime-c;
                ulong *restrict out=acc+i*cols;
                const ulong *restrict in=nmod_mat_entry_ptr(src,r,0);
                for(slong v=0;v<cols;v++)out[v]+=c*in[v];
            }
        }
        if(b) {
            const nmod_mat_struct *src=state+child*(M+1)+a;
            const linear_map *operation=op+((k-1)*N+j)*L+b-1;
            const nmod_mat_struct *map=operation->y;
            for(slong i=0;i<rows;i++)for(slong r=0;r<src->c;r++) {
                ulong c=nmod_mat_entry(src,i,r);if(!c)continue;
                if(neg)c=prime-c;
                ulong *restrict out=acc+i*cols;
                const ulong *restrict in=nmod_mat_entry_ptr(map,r,0);
                slong lo=operation->ystart?operation->ystart[r]:0;
                slong hi=operation->yend?operation->yend[r]:cols;
                for(slong v=lo;v<hi;v++)out[v]+=c*in[v];
            }
        }
    }
    for(slong i=0;i<rows;i++)for(slong j=0;j<cols;j++)
        nmod_mat_entry(dest,i,j)=nmod_set_ui(acc[i*cols+j],mod);
    flint_free(acc);return 1;
}

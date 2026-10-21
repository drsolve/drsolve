/* SPDX-License-Identifier: GPL-2.0-or-later */
/* PML-derived construction: Copyright (C) 2025 Vincent Neiger.
 * Redistributed under GPL version 2 or any later version. */
/* Exact two-branch LNZ determinant recursion, adapted from the structure
 * in pml_det/src/nmod_poly_mat_extra/nmod_poly_mat_det.c. Kernel bases are
 * computed through unified PMBasis. No probabilistic saturation shortcut.
 * A = [Au;Ad], N = ker(Au), Vu = leftker(N), Au=B1*Vu.
 * det(A) = det(Vu(0)[:,comp])/det([I[:,comp],N(0)]) * det(B1)*det(Ad*N).
 * A kernel [K1,K2] of [Vu;Au] has K1=-K2*B1 and unimodular K2. */

static int UM_NAME(det)(unified_mat_poly_t out, const unified_mat_poly_t a)
{
    field_ctx_t *ctx=a->ctx;
    slong n=a->r;
    if(!n) { UM_NAME(one)(out); return 1; }
    if(!a->length) { UM_NAME(resize)(out,0); return 1; }
    if(n==1) { UM_NAME(select)(out,a,NULL,NULL,0); return 1; }
    if(a->length==1) {
        UM_NAME(resize)(out,1);
        UM_NAME(const_eliminate)((UM_E *)out->data,NULL,a);
        UM_NAME(normalise)(out); return 1;
    }
    if(n==2) {
        unified_mat_poly_t x,y,z,w,p,q;
        unified_mat_poly_init(x,1,1,ctx);unified_mat_poly_init(y,1,1,ctx);
        unified_mat_poly_init(z,1,1,ctx);unified_mat_poly_init(w,1,1,ctx);
        unified_mat_poly_init(p,1,1,ctx);unified_mat_poly_init(q,1,1,ctx);
        slong zero=0,one=1;
        UM_NAME(select)(x,a,&zero,&zero,0);UM_NAME(select)(y,a,&one,&one,0);
        UM_NAME(select)(z,a,&zero,&one,0);UM_NAME(select)(w,a,&one,&zero,0);
        UM_NAME(mul)(p,x,y); UM_NAME(mul)(q,z,w);
        slong len=FLINT_MAX(p->length,q->length);
        UM_NAME(resize)(p,len); UM_NAME(resize)(q,len);
        for(slong i=0;i<len;i++) UM_SUB((UM_E *)p->data+i,(UM_E *)p->data+i,(UM_E *)q->data+i);
        UM_NAME(normalise)(p); UM_NAME(swap)(out,p);
        UM_NAME(clear)(q);UM_NAME(clear)(p);UM_NAME(clear)(w);
        UM_NAME(clear)(z);UM_NAME(clear)(y);UM_NAME(clear)(x);
        return 1;
    }
    /* Triangular matrices (including zero rows) need no kernel recursion. */
    int upper=1,lower=1;
    for(slong i=0;i<n;i++) for(slong j=0;j<n;j++) if(UM_NAME(entry_degree)(a,i,j)>=0) {
        if(i>j) upper=0;
        if(i<j) lower=0;
    }
    if(upper || lower) {
        unified_mat_poly_t p,x;
        unified_mat_poly_init(p,1,1,ctx);unified_mat_poly_init(x,1,1,ctx);UM_NAME(one)(p);
        for(slong i=0;i<n;i++) { UM_NAME(select)(x,a,&i,&i,0);UM_NAME(mul)(p,p,x); }
        UM_NAME(swap)(out,p);UM_NAME(clear)(p);UM_NAME(clear)(x);return 1;
    }

    slong m=(n+1)/2,k=n-m;
    double started=get_wall_time(),kernel_time=0,popov_time=0,mark;
    int certified=0;
    unified_mat_poly_t au,ad,at,nt,nmat,vu,stack,kmat,k1,k20,b2,d1,d2,ustar,vuc;
    unified_mat_poly_init(au,m,n,ctx); unified_mat_poly_init(ad,k,n,ctx);
    unified_mat_poly_init(at,n,m,ctx); unified_mat_poly_init(nt,k,n,ctx);
    unified_mat_poly_init(nmat,n,k,ctx); unified_mat_poly_init(vu,m,n,ctx);
    unified_mat_poly_init(stack,2*m,m,ctx); unified_mat_poly_init(kmat,m,2*m,ctx);
    unified_mat_poly_init(k1,m,m,ctx); unified_mat_poly_init(k20,m,m,ctx);
    unified_mat_poly_init(b2,k,k,ctx); unified_mat_poly_init(d1,1,1,ctx);
    unified_mat_poly_init(d2,1,1,ctx); unified_mat_poly_init(ustar,n,n,ctx);
    unified_mat_poly_init(vuc,m,m,ctx);
    slong *bottom=flint_malloc(k*sizeof(slong)),*profile=flint_malloc(n*sizeof(slong));
    slong *comp=flint_malloc(m*sizeof(slong)),*right=flint_malloc(m*sizeof(slong));
    unsigned char *used=flint_calloc(n,1);
    UM_E du,dv,c2,inv,scalar;UM_INIT(&du);UM_INIT(&dv);UM_INIT(&c2);UM_INIT(&inv);UM_INIT(&scalar);
    int ok=0;
    for(slong i=0;i<k;i++) bottom[i]=m+i;
    for(slong i=0;i<m;i++) right[i]=m+i;
    UM_NAME(select)(au,a,NULL,NULL,0);UM_NAME(select)(ad,a,bottom,NULL,0);
    mark=get_wall_time();
    slong au_degree=UM_NAME(weak_popov)(au);
    popov_time+=get_wall_time()-mark;
    if(au_degree<0) { UM_NAME(resize)(out,0);ok=1;goto done; }
    UM_NAME(select)(at,au,NULL,NULL,1);
    mark=get_wall_time();
    if(!UM_NAME(kernel)(nt,at)) goto done;
    kernel_time+=get_wall_time()-mark;mark=get_wall_time();
    slong kernel_degree=UM_NAME(weak_popov)(nt);
    popov_time+=get_wall_time()-mark;
    if(kernel_degree<0) goto done;
    UM_NAME(select)(nmat,nt,NULL,NULL,1);
    if(UM_NAME(const_eliminate)(&du,profile,nmat)!=k) goto done;
    for(slong i=0;i<k;i++) used[profile[i]]=1;
    for(slong i=0,j=0;i<n;i++) if(!used[i]) comp[j++]=i;
    UM_NAME(resize)(ustar,1);
    UM_E *u=ustar->data;const UM_E *nd=nmat->data;
    for(slong j=0;j<m;j++) UM_ONE(u+comp[j]*n+j);
    for(slong i=0;i<n;i++) for(slong j=0;j<k;j++) UM_SET(u+i*n+m+j,nd+i*k+j);
    UM_NAME(const_eliminate)(&du,NULL,ustar);
    if(UM_ISZERO(&du)) goto done;

    UM_NAME(mul)(b2,ad,nmat);
    /* Deterministic saturation certificate. For dual minimal bases,
     * sum(rowdeg(Vu)) = sum(coldeg(N)). Since Au is row reduced and
     * Au=B1*Vu, sum(rowdeg(Au)) differs by deg(det(B1)). Equality makes
     * B1 unimodular. Thus Au itself may serve as Vu, without a random
     * evaluation (which is particularly unsafe over small fields). */
    if(au_degree==kernel_degree) {
        certified=1;
        UM_NAME(select)(vuc,au,NULL,comp,0);
        UM_NAME(const_eliminate)(&dv,NULL,vuc);
        if(UM_ISZERO(&dv)) goto done;
        if(!UM_NAME(det)(d2,b2)) goto done;
        UM_INV(&inv,&du);UM_MUL(&scalar,&dv,&inv);
        for(slong i=0;i<d2->length;i++) UM_MUL((UM_E *)d2->data+i,(UM_E *)d2->data+i,&scalar);
        UM_NAME(normalise)(d2);UM_NAME(swap)(out,d2);ok=1;goto done;
    }

    if(!UM_NAME(kernel)(vu,nmat)) goto done;
    UM_NAME(select)(vuc,vu,NULL,comp,0);
    UM_NAME(const_eliminate)(&dv,NULL,vuc);
    if(UM_ISZERO(&dv)) goto done;
    /* Selecting independent columns preserves the left kernel of [Vu;Au].
     * This makes the kernel input full-column-rank, including odd n. */
    UM_NAME(resize)(stack,FLINT_MAX(vu->length,au->length));
    UM_E *sd=stack->data;const UM_E *vd=vu->data,*aud=au->data;
    for(slong d=0;d<stack->length;d++) for(slong i=0;i<m;i++) for(slong j=0;j<m;j++) {
        if(d<vu->length) UM_SET(sd+(d*2*m+i)*m+j,vd+(d*m+i)*n+comp[j]);
        if(d<au->length) UM_SET(sd+(d*2*m+m+i)*m+j,aud+(d*m+i)*n+comp[j]);
    }
    if(!UM_NAME(kernel)(kmat,stack)) goto done;
    UM_NAME(select)(k1,kmat,NULL,NULL,0);UM_NAME(select)(k20,kmat,NULL,right,0);
    UM_NAME(const_eliminate)(&c2,NULL,k20);
    if(UM_ISZERO(&c2)) goto done;
    if(!UM_NAME(det)(d1,k1) || !UM_NAME(det)(d2,b2)) goto done;
    UM_MUL(&scalar,&du,&c2);UM_INV(&inv,&scalar);UM_MUL(&scalar,&dv,&inv);
    if(m&1) UM_NEG(&scalar,&scalar);
    UM_NAME(mul)(d1,d1,d2);
    for(slong i=0;i<d1->length;i++) UM_MUL((UM_E *)d1->data+i,(UM_E *)d1->data+i,&scalar);
    UM_NAME(normalise)(d1);UM_NAME(swap)(out,d1);ok=1;
done:
    if(g_dixon_verbose_level>=3 && n>=16)
        printf("  Unified HNF node: size=%ld degree=%ld kernel=%.3fs popov=%.3fs total=%.3fs certified=%d success=%d\n",
               n,a->length-1,kernel_time,popov_time,get_wall_time()-started,certified,ok);
    UM_CLEAR(&scalar);UM_CLEAR(&inv);UM_CLEAR(&c2);UM_CLEAR(&dv);UM_CLEAR(&du);
    flint_free(used);flint_free(right);flint_free(comp);flint_free(profile);flint_free(bottom);
    UM_NAME(clear)(vuc);UM_NAME(clear)(ustar);UM_NAME(clear)(d2);UM_NAME(clear)(d1);
    UM_NAME(clear)(b2);UM_NAME(clear)(k20);UM_NAME(clear)(k1);UM_NAME(clear)(kmat);
    UM_NAME(clear)(stack);UM_NAME(clear)(vu);UM_NAME(clear)(nmat);UM_NAME(clear)(nt);
    UM_NAME(clear)(at);UM_NAME(clear)(ad);UM_NAME(clear)(au);
    return ok;
}

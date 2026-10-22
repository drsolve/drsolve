/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Instantiated for compact bytes and generic unified field elements.
 * PMBasis follows the divide-and-conquer construction also used in PML
 * (Copyright (C) 2025 Vincent Neiger, GPL-2.0-or-later). */

static void UM_NAME(clear)(unified_mat_poly_t a)
{
    field_ctx_t *ctx=a->ctx;
    UM_E *p=a->data;
    for(slong i=0;i<a->alloc*a->r*a->c;i++) UM_CLEAR(p+i);
    flint_free(p); a->data=NULL; a->alloc=a->length=0;
    (void)ctx;
}

static void UM_NAME(resize)(unified_mat_poly_t a, slong len)
{
    field_ctx_t *ctx=a->ctx;
    const slong size=a->r*a->c;
    if(len>a->alloc) {
        if(size && (size_t)len>SIZE_MAX/(size_t)size/sizeof(UM_E)) flint_abort();
        a->data=flint_realloc(a->data,(size_t)len*size*sizeof(UM_E));
        UM_E *p=a->data;
        for(slong i=a->alloc*size;i<len*size;i++) { UM_INIT(p+i); UM_ZERO(p+i); }
        a->alloc=len;
    }
    UM_E *p=a->data;
    /* Clear both newly exposed and discarded coefficients, including
     * allocated coefficients beyond length which a later resize may reuse. */
    for(slong i=FLINT_MIN(len,a->length)*size;i<FLINT_MAX(len,a->length)*size;i++) UM_ZERO(p+i);
    a->length=len;
    (void)ctx;
}

static void UM_NAME(normalise)(unified_mat_poly_t a)
{
    field_ctx_t *ctx=a->ctx;
    UM_E *p=a->data;
    slong size=a->r*a->c;
    while(a->length) {
        slong i=0;
        for(;i<size;i++) if(!UM_ISZERO(p+(a->length-1)*size+i)) break;
        if(i<size) break;
        --a->length;
    }
    (void)ctx;
}

static void UM_NAME(swap)(unified_mat_poly_t a, unified_mat_poly_t b)
{
    unified_mat_poly_struct t=*a; *a=*b; *b=t;
}

static void UM_NAME(one)(unified_mat_poly_t a)
{
    field_ctx_t *ctx=a->ctx;
    UM_NAME(resize)(a,0); UM_NAME(resize)(a,1);
    UM_E *p=a->data;
    for(slong i=0;i<FLINT_MIN(a->r,a->c);i++) UM_ONE(p+i*a->c+i);
    (void)ctx;
}

/* Copy an indexed submatrix, optionally transposed. NULL indices mean
 * consecutive indices starting at zero. Always safe with aliased input. */
static void UM_NAME(select)(unified_mat_poly_t out, const unified_mat_poly_t a,
                            const slong *rows, const slong *cols, int transpose)
{
    field_ctx_t *ctx=a->ctx;
    unified_mat_poly_t t; unified_mat_poly_init(t,out->r,out->c,ctx);
    UM_NAME(resize)(t,a->length);
    UM_E *d=t->data; const UM_E *s=a->data;
    for(slong k=0;k<a->length;k++) for(slong i=0;i<out->r;i++) for(slong j=0;j<out->c;j++) {
        slong ri=transpose ? j : i, cj=transpose ? i : j;
        ri=rows ? rows[ri] : ri; cj=cols ? cols[cj] : cj;
        UM_SET(d+(k*out->r+i)*out->c+j,s+(k*a->r+ri)*a->c+cj);
    }
    UM_NAME(normalise)(t); UM_NAME(swap)(out,t); UM_NAME(clear)(t);
}

/* Coefficient-window convolution. Byte innermost loops are contiguous
 * row axpys with a single scalar multiplication table. No union accesses. */
static void UM_NAME(mul_window_classical)(unified_mat_poly_t out, const unified_mat_poly_t a,
                                const unified_mat_poly_t b, slong start, slong len)
{
    field_ctx_t *ctx=a->ctx;
    unified_mat_poly_t t; unified_mat_poly_init(t,a->r,b->c,ctx);
    slong total=a->length && b->length ? a->length+b->length-1 : 0;
    len=FLINT_MIN(len,FLINT_MAX(0,total-start));
    UM_NAME(resize)(t,len);
    const UM_E *ad=a->data,*bd=b->data; UM_E *td=t->data;
#if UM_BYTE
    const uint8_t (*tables)[32]=unified_byte_multipliers(ctx);
#endif
#ifdef _OPENMP
    double work=(double)a->r*a->c*b->c*a->length*FLINT_MIN(b->length,len);
    int workers=FLINT_MIN(4,omp_get_max_threads());
#pragma omp parallel for schedule(static) num_threads(workers) if(work>2000000 && a->r>=4 && !omp_in_parallel())
#endif
    for(slong i=0;i<a->r;i++) {
#ifdef UM_RAW_ENABLE
      /* Reduction is linear over GF(2). Accumulate unreduced products in a
       * private output row, then reduce once per coefficient, not per term. */
      if(UM_RAW_ENABLE && a->c>1 && len>0) {
        __m128i *acc=flint_calloc(len*b->c,sizeof(*acc));
        for(slong u=0;u<a->length;u++) {
          slong lo=FLINT_MAX(0,start-u),hi=FLINT_MIN(b->length,start+len-u);
          for(slong v=lo;v<hi;v++) for(slong k=0;k<a->c;k++) {
            const UM_E *scalar=ad+(u*a->r+i)*a->c+k;
            if(UM_ISZERO(scalar)) continue;
            __m128i x=_mm_cvtsi64_si128(UM_RAW_WORD(scalar));
            const UM_E *row=bd+(v*b->r+k)*b->c;
            __m128i *dst=acc+(u+v-start)*b->c;
            for(slong j=0;j<b->c;j++)
              dst[j]=_mm_xor_si128(dst[j],_mm_clmulepi64_si128(x,
                  _mm_cvtsi64_si128(UM_RAW_WORD(row+j)),0));
          }
        }
        for(slong d=0;d<len;d++) for(slong j=0;j<b->c;j++)
          td[(d*a->r+i)*b->c+j]=UM_RAW_REDUCE(acc[d*b->c+j]);
        flint_free(acc);
        continue;
      }
#endif
      UM_E prod; UM_INIT(&prod);
      for(slong u=0;u<a->length;u++) {
        slong lo=FLINT_MAX(0,start-u), hi=FLINT_MIN(b->length,start+len-u);
        for(slong v=lo;v<hi;v++) for(slong k=0;k<a->c;k++) {
            const UM_E *scalar=ad+(u*a->r+i)*a->c+k;
            if(UM_ISZERO(scalar)) continue;
            const UM_E *row=bd+(v*b->r+k)*b->c;
            UM_E *dst=td+((u+v-start)*a->r+i)*b->c;
#if UM_BYTE
            unified_byte_addmul(dst,row,b->c,tables[*scalar]);
#else
            for(slong j=0;j<b->c;j++) {
                UM_MUL(&prod,scalar,row+j); UM_ADD(dst+j,dst+j,&prod);
            }
#endif
        }
      }
      UM_CLEAR(&prod);
    }
    UM_NAME(normalise)(t); UM_NAME(swap)(out,t); UM_NAME(clear)(t);
}

static void UM_NAME(slice)(unified_mat_poly_t out,const unified_mat_poly_t a,slong start,slong len)
{
    field_ctx_t *ctx=a->ctx;
    len=FLINT_MIN(len,FLINT_MAX(0,a->length-start));
    UM_NAME(resize)(out,0);UM_NAME(resize)(out,len);
    UM_E *d=out->data;const UM_E *s=a->data;
    for(slong i=0;i<len*a->r*a->c;i++) UM_SET(d+i,s+start*a->r*a->c+i);
    UM_NAME(normalise)(out);
    (void)ctx;
}

static void UM_NAME(add)(unified_mat_poly_t out,const unified_mat_poly_t a,const unified_mat_poly_t b)
{
    field_ctx_t *ctx=a->ctx;
    UM_NAME(resize)(out,0);UM_NAME(resize)(out,FLINT_MAX(a->length,b->length));
    UM_E *d=out->data;const UM_E *ad=a->data,*bd=b->data;
    slong size=a->r*a->c;
    for(slong i=0;i<out->length*size;i++) {
        if(i<a->length*size && i<b->length*size) UM_ADD(d+i,ad+i,bd+i);
        else if(i<a->length*size) UM_SET(d+i,ad+i);
        else UM_SET(d+i,bd+i);
    }
    UM_NAME(normalise)(out);
    (void)ctx;
}

static void UM_NAME(accumulate_window)(unified_mat_poly_t out,const unified_mat_poly_t a,
                                        slong offset,slong start,int subtract)
{
    field_ctx_t *ctx=a->ctx;
    slong lo=FLINT_MAX(0,start-offset),hi=FLINT_MIN(a->length,start+out->length-offset);
    slong size=a->r*a->c;UM_E *d=out->data;const UM_E *s=a->data;
    for(slong k=lo;k<hi;k++) for(slong i=0;i<size;i++) {
        UM_E *v=d+(k+offset-start)*size+i;
        if(subtract) UM_SUB(v,v,s+k*size+i);else UM_ADD(v,v,s+k*size+i);
    }
    (void)ctx;
}

#ifndef UM_CUSTOM_MUL_WINDOW
/* Karatsuba works over the noncommutative ring of matrix coefficients:
 * preserve operand order in all three products. It avoids quadratic
 * degree growth costs at the smaller, higher-degree HNF recursion nodes.
 * Short windows retain the bounded classical convolution. */
static void UM_NAME(mul_window)(unified_mat_poly_t out,const unified_mat_poly_t a,
                                const unified_mat_poly_t b,slong start,slong len)
{
    slong low=FLINT_MIN(a->length,b->length),high=FLINT_MAX(a->length,b->length);
    if(low<16 || high>2*low || len<low/2) {
        UM_NAME(mul_window_classical)(out,a,b,start,len);return;
    }
    field_ctx_t *ctx=a->ctx;
    slong h=(high+1)/2;
    unified_mat_poly_t a0,a1,b0,b1,as,bs,z0,z1,z2,t;
    unified_mat_poly_init(a0,a->r,a->c,ctx);unified_mat_poly_init(a1,a->r,a->c,ctx);
    unified_mat_poly_init(b0,b->r,b->c,ctx);unified_mat_poly_init(b1,b->r,b->c,ctx);
    unified_mat_poly_init(as,a->r,a->c,ctx);unified_mat_poly_init(bs,b->r,b->c,ctx);
    unified_mat_poly_init(z0,a->r,b->c,ctx);unified_mat_poly_init(z1,a->r,b->c,ctx);
    unified_mat_poly_init(z2,a->r,b->c,ctx);unified_mat_poly_init(t,a->r,b->c,ctx);
    UM_NAME(slice)(a0,a,0,h);UM_NAME(slice)(a1,a,h,high-h);
    UM_NAME(slice)(b0,b,0,h);UM_NAME(slice)(b1,b,h,high-h);
    UM_NAME(add)(as,a0,a1);UM_NAME(add)(bs,b0,b1);
    UM_NAME(mul_window)(z0,a0,b0,0,2*h);
    UM_NAME(mul_window)(z1,as,bs,0,2*h);
    UM_NAME(mul_window)(z2,a1,b1,0,2*h);
    len=FLINT_MIN(len,FLINT_MAX(0,a->length+b->length-1-start));
    UM_NAME(resize)(t,len);
    UM_NAME(accumulate_window)(t,z0,0,start,0);
    UM_NAME(accumulate_window)(t,z0,h,start,1);
    UM_NAME(accumulate_window)(t,z1,h,start,0);
    UM_NAME(accumulate_window)(t,z2,h,start,1);
    UM_NAME(accumulate_window)(t,z2,2*h,start,0);
    UM_NAME(normalise)(t);UM_NAME(swap)(out,t);
    UM_NAME(clear)(t);UM_NAME(clear)(z2);UM_NAME(clear)(z1);UM_NAME(clear)(z0);
    UM_NAME(clear)(bs);UM_NAME(clear)(as);UM_NAME(clear)(b1);UM_NAME(clear)(b0);
    UM_NAME(clear)(a1);UM_NAME(clear)(a0);
}

#endif

static void UM_NAME(mul)(unified_mat_poly_t out, const unified_mat_poly_t a,
                         const unified_mat_poly_t b)
{
    UM_NAME(mul_window)(out,a,b,0,a->length && b->length ? a->length+b->length-1 : 0);
}

/* dst row -= scalar * src row, optionally shifted. */
static void UM_NAME(row_submul)(unified_mat_poly_t a, slong dst, slong src,
                                const UM_E *scalar, slong shift, slong first, slong last)
{
    field_ctx_t *ctx=a->ctx;
    UM_E prod; UM_INIT(&prod);
    UM_E *p=a->data;
#if UM_BYTE
    const uint8_t *table=unified_byte_multipliers(ctx)[*scalar];
#endif
    for(slong d=FLINT_MIN(last,a->length-1);d>=FLINT_MAX(first,shift);d--) {
        UM_E *to=p+(d*a->r+dst)*a->c;
        const UM_E *from=p+((d-shift)*a->r+src)*a->c;
#if UM_BYTE
        unified_byte_addmul(to,from,a->c,table);
#else
        for(slong j=0;j<a->c;j++) {
#ifdef UM_SUBMUL
            UM_SUBMUL(to+j,scalar,from+j);
#else
            UM_MUL(&prod,scalar,from+j); UM_SUB(to+j,to+j,&prod);
#endif
        }
#endif
    }
    UM_CLEAR(&prod);
}

static void UM_NAME(row_shift)(unified_mat_poly_t a, slong row, slong first, slong last)
{
    field_ctx_t *ctx=a->ctx;
    UM_E *p=a->data;
    for(slong d=FLINT_MIN(last,a->length-1);d>first;d--) for(slong j=0;j<a->c;j++)
        UM_SET(p+(d*a->r+row)*a->c+j,p+((d-1)*a->r+row)*a->c+j);
    for(slong j=0;j<a->c;j++) UM_ZERO(p+(first*a->r+row)*a->c+j);
    (void)ctx;
}

static void UM_NAME(mbasis)(unified_mat_poly_t out, slong *shift,
                            const unified_mat_poly_t a, slong order)
{
    field_ctx_t *ctx=a->ctx;
    slong n=a->r;
    unified_mat_poly_t basis,residual;
    unified_mat_poly_init(basis,n,n,ctx); unified_mat_poly_init(residual,n,a->c,ctx);
    UM_NAME(one)(basis);
    UM_NAME(slice)(residual,a,0,order);UM_NAME(resize)(residual,order);
    unsigned char *used=flint_calloc(n,1);
    /* Conservative row-degree bounds avoid updating zero padding. */
    slong *basis_degree=flint_calloc(n,sizeof(slong));
    UM_E inv,scalar; UM_INIT(&inv); UM_INIT(&scalar);
    for(slong d=0;d<order;d++) {
        memset(used,0,n);
        for(slong col=0;col<a->c;col++) {
            UM_E *r=(UM_E *)residual->data+d*n*a->c;
            slong pivot=-1;
            for(slong i=0;i<n;i++) if(!used[i] && !UM_ISZERO(r+i*a->c+col) &&
                (pivot<0 || shift[i]<shift[pivot] || (shift[i]==shift[pivot] && i<pivot))) pivot=i;
            if(pivot<0) continue;
            UM_INV(&inv,r+pivot*a->c+col);
            for(slong i=0;i<n;i++) if(!used[i] && i!=pivot && !UM_ISZERO(r+i*a->c+col)) {
                UM_MUL(&scalar,r+i*a->c+col,&inv);
                basis_degree[i]=FLINT_MAX(basis_degree[i],basis_degree[pivot]);
                UM_NAME(row_submul)(basis,i,pivot,&scalar,0,0,basis_degree[i]);
                /* Earlier residual coefficients are already zero for every row. */
                UM_NAME(row_submul)(residual,i,pivot,&scalar,0,d,residual->length-1);
            }
            used[pivot]=1;
        }
        UM_NAME(resize)(basis,basis->length+1);
        for(slong i=0;i<n;i++) if(used[i]) {
            UM_NAME(row_shift)(basis,i,0,++basis_degree[i]);
            UM_NAME(row_shift)(residual,i,d,residual->length-1);shift[i]++;
        }
        UM_NAME(normalise)(basis);
    }
    UM_CLEAR(&inv); UM_CLEAR(&scalar); flint_free(used); flint_free(basis_degree);
    UM_NAME(swap)(out,basis); UM_NAME(clear)(basis); UM_NAME(clear)(residual);
}

static void UM_NAME(pmbasis)(unified_mat_poly_t out, slong *shift,
                             const unified_mat_poly_t a, slong order)
{
    if(order<=8) { UM_NAME(mbasis)(out,shift,a,order); return; }
    slong half=order/2;
    unified_mat_poly_t first,second,residual;
    unified_mat_poly_init(first,a->r,a->r,a->ctx);
    unified_mat_poly_init(second,a->r,a->r,a->ctx);
    unified_mat_poly_init(residual,a->r,a->c,a->ctx);
    UM_NAME(pmbasis)(first,shift,a,half);
    UM_NAME(mul_window)(residual,first,a,half,order-half);
    UM_NAME(pmbasis)(second,shift,residual,order-half);
    UM_NAME(mul)(out,second,first);
    UM_NAME(clear)(residual); UM_NAME(clear)(second); UM_NAME(clear)(first);
}

static slong UM_NAME(entry_degree)(const unified_mat_poly_t a, slong row, slong col)
{
    field_ctx_t *ctx=a->ctx;
    const UM_E *p=a->data;
    slong d=a->length-1;
    while(d>=0 && UM_ISZERO(p+(d*a->r+row)*a->c+col)) --d;
    (void)ctx;
    return d;
}

/* Unshifted weak Popov reduction, using determinant-preserving row
 * additions only. Distinct rightmost leading positions certify row
 * reducedness. Return the sum of row degrees, or -1 for dependent rows. */
static slong UM_NAME(weak_popov)(unified_mat_poly_t a)
{
    field_ctx_t *ctx=a->ctx;
    slong *degree=flint_malloc(a->r*sizeof(slong));
    slong *pivot=flint_malloc(a->r*sizeof(slong));
    slong *owner=flint_malloc(a->c*sizeof(slong));
    for(slong i=0;i<a->r;i++) {
        degree[i]=-1;pivot[i]=-1;
        for(slong j=0;j<a->c;j++) {
            slong d=UM_NAME(entry_degree)(a,i,j);
            if(d>=0 && d>=degree[i]) { degree[i]=d;pivot[i]=j; }
        }
    }
    UM_E inv,scalar;UM_INIT(&inv);UM_INIT(&scalar);
    slong sum=0;
    for(;;) {
        for(slong j=0;j<a->c;j++) owner[j]=-1;
        slong dst=-1,src=-1;
        for(slong i=0;i<a->r;i++) {
            if(pivot[i]<0) { sum=-1;goto done_popov; }
            slong p=pivot[i];
            if(owner[p]>=0) {
                src=owner[p];dst=i;
                if(degree[dst]<degree[src]) { slong t=dst;dst=src;src=t; }
                break;
            }
            owner[p]=i;
        }
        if(dst<0) break;
        UM_E *data=a->data;
        UM_INV(&inv,data+(degree[src]*a->r+src)*a->c+pivot[src]);
        UM_MUL(&scalar,data+(degree[dst]*a->r+dst)*a->c+pivot[dst],&inv);
        UM_NAME(row_submul)(a,dst,src,&scalar,degree[dst]-degree[src],0,degree[dst]);
        degree[dst]=-1;pivot[dst]=-1;
        for(slong j=0;j<a->c;j++) {
            slong d=UM_NAME(entry_degree)(a,dst,j);
            if(d>=0 && d>=degree[dst]) { degree[dst]=d;pivot[dst]=j; }
        }
    }
    for(slong i=0;i<a->r;i++) sum+=degree[i];
done_popov:
    UM_CLEAR(&scalar);UM_CLEAR(&inv);flint_free(owner);flint_free(pivot);flint_free(degree);
    UM_NAME(normalise)(a);return sum;
}

/* Extract exact zero rows of P*A, not merely approximants. Once there are
 * n-m such rows and A has rank m, they are a saturated kernel basis:
 * P is a polynomial basis of the approximant module, and its remaining
 * m image rows are rationally independent. */
static int UM_NAME(kernel)(unified_mat_poly_t out, const unified_mat_poly_t a)
{
    const slong n=a->r,m=a->c,k=n-m;
    if(!k) { UM_NAME(resize)(out,0); return 1; }
    slong *shift=flint_calloc(n,sizeof(slong)), *original=flint_calloc(n,sizeof(slong));
    slong *rows=flint_malloc(n*sizeof(slong));
    slong maxshift=0,sumshift=0;
    for(slong i=0;i<n;i++) for(slong j=0;j<m;j++)
        original[i]=FLINT_MAX(original[i],UM_NAME(entry_degree)(a,i,j));
    for(slong i=0;i<n;i++) { maxshift=FLINT_MAX(maxshift,original[i]);sumshift+=original[i]; }
    /* Conservative complete order bound; adaptive attempts usually finish
     * near twice the entry degree. No field-size-dependent sampling. */
    if(maxshift && m+2>LONG_MAX/maxshift) { flint_free(rows);flint_free(original);flint_free(shift);return 0; }
    slong bound=(m+2)*maxshift+1;
    slong order=FLINT_MIN(bound,FLINT_MAX(maxshift+1,sumshift/k+1));
    unified_mat_poly_t p,product;
    unified_mat_poly_init(p,n,n,a->ctx); unified_mat_poly_init(product,n,m,a->ctx);
    int ok=0;
    for(;;) {
        memcpy(shift,original,n*sizeof(slong));
        UM_NAME(pmbasis)(p,shift,a,order);
        UM_NAME(mul)(product,p,a);
        slong found=0;
        for(slong i=0;i<n;i++) {
            int zero=1;
            for(slong j=0;j<m && zero;j++) zero=UM_NAME(entry_degree)(product,i,j)<0;
            if(zero) rows[found++]=i;
        }
        if(found==k) { UM_NAME(select)(out,p,rows,NULL,0);ok=1;break; }
        if(found>k || order==bound) break;
        order=order>bound/2 ? bound : order*2;
    }
    UM_NAME(clear)(product); UM_NAME(clear)(p);
    flint_free(rows); flint_free(original); flint_free(shift);
    return ok;
}

/* Constant Gaussian elimination, also returns a row rank profile when
 * requested. The input is copied; determinant is zero if nonsquare/rank
 * deficient. The permutation sign matters for odd characteristic. */
static slong UM_NAME(const_eliminate)(UM_E *det, slong *profile, const unified_mat_poly_t a)
{
    field_ctx_t *ctx=a->ctx;
    unified_mat_poly_t t; unified_mat_poly_init(t,a->r,a->c,ctx);
    UM_NAME(resize)(t,1);
    UM_E *p=t->data; const UM_E *s=a->data;
    if(a->length) for(slong i=0;i<a->r*a->c;i++) UM_SET(p+i,s+i);
    slong *perm=flint_malloc(a->r*sizeof(slong));
    for(slong i=0;i<a->r;i++) perm[i]=i;
    UM_E inv,scalar,tmp; UM_INIT(&inv);UM_INIT(&scalar);UM_INIT(&tmp);UM_ONE(det);
    slong rank=0;
    for(slong j=0;j<a->c && rank<a->r;j++) {
        slong pivot=rank;
        while(pivot<a->r && UM_ISZERO(p+pivot*a->c+j)) pivot++;
        if(pivot==a->r) continue;
        if(pivot!=rank) {
            slong q=perm[pivot];perm[pivot]=perm[rank];perm[rank]=q;
            for(slong c=0;c<a->c;c++) {
                UM_SET(&tmp,p+pivot*a->c+c); UM_SET(p+pivot*a->c+c,p+rank*a->c+c); UM_SET(p+rank*a->c+c,&tmp);
            }
            UM_NEG(det,det);
        }
        UM_MUL(det,det,p+rank*a->c+j); UM_INV(&inv,p+rank*a->c+j);
        for(slong i=rank+1;i<a->r;i++) if(!UM_ISZERO(p+i*a->c+j)) {
            UM_MUL(&scalar,p+i*a->c+j,&inv);
            UM_NAME(row_submul)(t,i,rank,&scalar,0,0,0);
        }
        rank++;
    }
    if(profile) memcpy(profile,perm,a->r*sizeof(slong));
    if(a->r!=a->c || rank!=a->r) UM_ZERO(det);
    UM_CLEAR(&tmp);UM_CLEAR(&scalar);UM_CLEAR(&inv);
    flint_free(perm); UM_NAME(clear)(t); return rank;
}

#include "../determinant/unified_hnf_template.h"

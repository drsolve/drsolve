/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Whole-operation instantiations retain compact binary/fmpz arithmetic. */
static void UM_NAME(ib_evaluate)(unified_mat_poly_t out,
    const unified_mat_poly_t a, const unified_mat_poly_t points, slong index)
{
    field_ctx_t *ctx=a->ctx;
    const UM_E *x=(const UM_E *)points->data+index, *s=a->data;
    UM_NAME(resize)(out,0);UM_NAME(resize)(out,1);
    UM_E *d=out->data;
    for(slong k=a->length;k-- >0;) for(slong i=0;i<a->r*a->c;i++) {
        UM_MUL(d+i,d+i,x);UM_ADD(d+i,d+i,s+k*a->r*a->c+i);
    }
    (void)ctx;
}

/* Solve A X=V by elimination; det is returned even when V has no columns.
 * All temporaries are private to a sampling worker. */
static void UM_NAME(ib_sample)(unified_mat_poly_t block, unified_mat_poly_t determinant,
    const unified_mat_poly_t a, const unified_mat_poly_t u, const unified_mat_poly_t v,
    const unified_mat_poly_t points, slong index, int *ok)
{
    field_ctx_t *ctx=a->ctx;slong n=a->r,m=v->c;
    unified_mat_poly_t val,aug,x,h;
    unified_mat_poly_init(val,n,n,ctx);unified_mat_poly_init(aug,n,n+m,ctx);
    unified_mat_poly_init(x,n,m,ctx);unified_mat_poly_init(h,m,m,ctx);
    UM_NAME(ib_evaluate)(val,a,points,index);UM_NAME(resize)(aug,1);
    UM_NAME(resize)(x,1);UM_NAME(one)(determinant);
    UM_E *p=aug->data,*det=determinant->data,*xd=x->data;
    const UM_E *vd=v->data,*ad=val->data;
    for(slong i=0;i<n;i++) {
        for(slong j=0;j<n;j++) UM_SET(p+i*(n+m)+j,ad+i*n+j);
        for(slong j=0;j<m;j++) UM_SET(p+i*(n+m)+n+j,vd+i*m+j);
    }
    UM_E inv,scalar,tmp;UM_INIT(&inv);UM_INIT(&scalar);UM_INIT(&tmp);
    *ok=1;
    for(slong k=0;k<n;k++) {
        slong pivot=k;
        while(pivot<n && UM_ISZERO(p+pivot*(n+m)+k)) pivot++;
        if(pivot==n) { UM_ZERO(det);*ok=0;break; }
        if(pivot!=k) {
            for(slong j=k;j<n+m;j++) {
                UM_SET(&tmp,p+k*(n+m)+j);UM_SET(p+k*(n+m)+j,p+pivot*(n+m)+j);
                UM_SET(p+pivot*(n+m)+j,&tmp);
            }
            UM_NEG(det,det);
        }
        UM_MUL(det,det,p+k*(n+m)+k);UM_INV(&inv,p+k*(n+m)+k);
        for(slong i=k+1;i<n;i++) if(!UM_ISZERO(p+i*(n+m)+k)) {
            UM_MUL(&scalar,p+i*(n+m)+k,&inv);
            for(slong j=k+1;j<n+m;j++) {
#ifdef UM_SUBMUL
                UM_SUBMUL(p+i*(n+m)+j,&scalar,p+k*(n+m)+j);
#else
                UM_MUL(&tmp,&scalar,p+k*(n+m)+j);UM_SUB(p+i*(n+m)+j,p+i*(n+m)+j,&tmp);
#endif
            }
            UM_ZERO(p+i*(n+m)+k);
        }
    }
    if(*ok && m) {
        for(slong i=n;i-- >0;) {
            UM_INV(&inv,p+i*(n+m)+i);
            for(slong j=0;j<m;j++) {
                UM_SET(&tmp,p+i*(n+m)+n+j);
                for(slong k=i+1;k<n;k++) {
                    UM_MUL(&scalar,p+i*(n+m)+k,xd+k*m+j);UM_SUB(&tmp,&tmp,&scalar);
                }
                UM_MUL(xd+i*m+j,&tmp,&inv);
            }
        }
        UM_NAME(mul)(h,u,x);
        UM_NAME(resize)(block,0);UM_NAME(resize)(block,1);
        UM_E *bd=block->data;
        for(slong i=0;i<m*m;i++) UM_SET(bd+i,(UM_E *)h->data+i);
        for(slong i=0;i<m;i++) UM_ONE(bd+(m+i)*m+i);
    }
    UM_CLEAR(&inv);UM_CLEAR(&scalar);UM_CLEAR(&tmp);
    UM_NAME(clear)(val);UM_NAME(clear)(aug);UM_NAME(clear)(x);UM_NAME(clear)(h);
}

/* Minimal interpolant basis for B(t_i) E_i=0. The pivot order and shifted
 * degrees follow MBasis; a pivot row is multiplied by (x-t_i), not x. */
static void UM_NAME(ib_mintbasis)(unified_mat_poly_t out, slong *shift,
    const unified_mat_poly_t points, slong start,
    const unified_mat_poly_struct *samples, slong count)
{
    field_ctx_t *ctx=out->ctx;slong r=out->r,c=samples->c;
    unified_mat_poly_t value,residual;
    unified_mat_poly_init(value,r,r,ctx);unified_mat_poly_init(residual,r,c,ctx);
    UM_NAME(one)(out);
    unsigned char *used=flint_calloc(r,1);
    UM_E inv,scalar,tmp;UM_INIT(&inv);UM_INIT(&scalar);UM_INIT(&tmp);
    for(slong t=0;t<count;t++) {
        const UM_E *point=(UM_E *)points->data+start+t;
        UM_NAME(ib_evaluate)(value,out,points,start+t);
        UM_NAME(mul)(residual,value,samples+t);
        UM_NAME(resize)(residual,1);
        memset(used,0,r);
        for(slong col=0;col<c;col++) {
            UM_E *e=residual->data;slong pivot=-1;
            for(slong i=0;i<r;i++) if(!used[i] && !UM_ISZERO(e+i*c+col) &&
                (pivot<0 || shift[i]<shift[pivot] || (shift[i]==shift[pivot] && i<pivot))) pivot=i;
            if(pivot<0) continue;
            UM_INV(&inv,e+pivot*c+col);
            for(slong i=0;i<r;i++) if(!used[i] && i!=pivot && !UM_ISZERO(e+i*c+col)) {
                UM_MUL(&scalar,e+i*c+col,&inv);
                UM_NAME(row_submul)(out,i,pivot,&scalar,0,0,out->length-1);
                UM_NAME(row_submul)(residual,i,pivot,&scalar,0,0,0);
            }
            used[pivot]=1;
        }
        slong length=out->length;UM_NAME(resize)(out,length+1);UM_E *b=out->data;
        for(slong i=0;i<r;i++) if(used[i]) {
            for(slong d=length;d>=0;d--) for(slong j=0;j<r;j++) {
                UM_E *dest=b+(d*r+i)*r+j;
                UM_MUL(&tmp,dest,point);
                if(d) UM_SUB(dest,b+((d-1)*r+i)*r+j,&tmp);
                else UM_NEG(dest,&tmp);
            }
            shift[i]++;
        }
        UM_NAME(normalise)(out);
    }
    UM_CLEAR(&inv);UM_CLEAR(&scalar);UM_CLEAR(&tmp);flint_free(used);
    UM_NAME(clear)(value);UM_NAME(clear)(residual);
}

static void UM_NAME(ib_pmintbasis)(unified_mat_poly_t out, slong *shift,
    const unified_mat_poly_t points, slong start,
    const unified_mat_poly_struct *samples, slong count)
{
    if(count<=8) {UM_NAME(ib_mintbasis)(out,shift,points,start,samples,count);return;}
    slong half=count/2,r=out->r,c=samples->c;
    unified_mat_poly_t first,second,value;
    unified_mat_poly_init(first,r,r,out->ctx);unified_mat_poly_init(second,r,r,out->ctx);
    unified_mat_poly_init(value,r,r,out->ctx);
    UM_NAME(ib_pmintbasis)(first,shift,points,start,samples,half);
    unified_mat_poly_struct *residual=flint_malloc((count-half)*sizeof(*residual));
    for(slong i=0;i<count-half;i++) {
        unified_mat_poly_init(residual+i,r,c,out->ctx);
        UM_NAME(ib_evaluate)(value,first,points,start+half+i);
        UM_NAME(mul)(residual+i,value,samples+half+i);
    }
    UM_NAME(ib_pmintbasis)(second,shift,points,start+half,residual,count-half);
    UM_NAME(mul)(out,second,first);
    for(slong i=0;i<count-half;i++) UM_NAME(clear)(residual+i);
    flint_free(residual);UM_NAME(clear)(first);UM_NAME(clear)(second);UM_NAME(clear)(value);
}

static void UM_NAME(ib_degrees)(slong *degrees, const unified_mat_poly_t a)
{
    for(slong i=0;i<a->r;i++) for(slong j=0;j<a->c;j++)
        degrees[i*a->c+j]=UM_NAME(entry_degree)(a,i,j);
}

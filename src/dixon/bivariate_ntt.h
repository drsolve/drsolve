/* Optional bivariate resultant over a word-size prime field.
 * Included by dixon_flint.c. The leading x coefficients must be constant,
 * so every specialization preserves degrees and the resultant commutes
 * with evaluation. Unsupported inputs leave result untouched. */
#include <stdint.h>
static void bivariate_ntt(ulong *a, slong n, const ulong *roots, nmod_t mod)
{
    for (slong i=1,j=0;i<n;i++) {
        slong bit=n>>1;
        for (;j&bit;bit>>=1) j^=bit;
        j^=bit;
        if(i<j){ulong t=a[i];a[i]=a[j];a[j]=t;}
    }
    for(slong len=2;len<=n;len<<=1)
        for(slong i=0;i<n;i+=len)
            for(slong j=0;j<len/2;j++) {
                ulong u=a[i+j],v=nmod_mul(a[i+j+len/2],roots[j*(n/len)],mod);
                a[i+j]=nmod_add(u,v,mod);
                a[i+j+len/2]=nmod_sub(u,v,mod);
            }
}

static int dixon_try_bivariate_ntt(unified_mpoly_struct *result,
    unified_mpoly_struct *polys, slong nvars, slong npars)
{
    if(!g_dixon_bezout_fast) return 0;
    const char *reason="requires two polynomials, one parameter and a prime field";
    if(nvars!=1 || npars!=1 || fq_nmod_ctx_degree(polys[0].ctx)!=1 ||
       g_field_equation_reduction) goto fallback;
    slong dx[2]={0},dt[2]={0};
    reason="degree arithmetic overflow";
    for(int k=0;k<2;k++) for(slong j=0;j<dr_mpoly_length(polys+k);j++) {
        DR_MPOLY_TERM(term,polys+k,j);
        slong x=term.var_exp ? term.var_exp[0] : 0;
        slong t=term.par_exp ? term.par_exp[0] : 0;
        if(x<0 || t<0 || x>WORD_MAX-t) goto fallback;
        dx[k]=FLINT_MAX(dx[k],x);dt[k]=FLINT_MAX(dt[k],x+t);
    }
    reason="requires positive x degrees and constant leading x coefficients";
    if(!dx[0] || !dx[1]) goto fallback;
    for(int k=0;k<2;k++) for(slong j=0;j<dr_mpoly_length(polys+k);j++) {
        DR_MPOLY_TERM(term,polys+k,j);
        if(term.var_exp && term.var_exp[0]==dx[k] &&
           term.par_exp && term.par_exp[0]) goto fallback;
    }
    /* Sylvester determinant weight: n*total(f)+m*total(g)-m*n.
     * Only integer/array representation checks; no degree or memory budget. */
    reason="degree or array length overflow";
    if(dx[1]>(WORD_MAX-1)/dt[0])goto fallback;
    ulong bound=dx[1]*dt[0];
    if((ulong)(dt[1]-dx[1])>((ulong)WORD_MAX-1-bound)/(ulong)dx[0])goto fallback;
    bound+=(ulong)dx[0]*(dt[1]-dx[1]);
    slong length=1;
    while((ulong)length<=bound) {
        if(length>WORD_MAX/2)goto fallback;
        length*=2;
    }
    if(dx[1]>WORD_MAX-2 || dx[0]>WORD_MAX-2-dx[1])goto fallback;
    slong rows=dx[0]+dx[1]+2;
    if((size_t)rows>SIZE_MAX/sizeof(ulong)/(size_t)length)goto fallback;
    ulong prime=fq_nmod_ctx_prime(polys[0].ctx);
    reason="prime field lacks roots of unity of the required order";
    if(prime<3 || (prime-1)%length) goto fallback;
    nmod_t mod;nmod_init(&mod,prime);
    ulong root=1;
    for(ulong a=2;a<prime && a<130;a++) {
        root=nmod_pow_ui(a,(prime-1)/length,mod);
        if(nmod_pow_ui(root,length/2,mod)==prime-1) break;
    }
    if(nmod_pow_ui(root,length/2,mod)!=prime-1) goto fallback;
    double start=get_wall_time();
    ulong *values=flint_calloc((size_t)rows*length,sizeof(ulong));
    ulong *roots=flint_malloc(length*sizeof(ulong));
    ulong *inverse=flint_malloc(length*sizeof(ulong));
    ulong *samples=flint_malloc(length*sizeof(ulong));
    roots[0]=inverse[0]=1;
    ulong invroot=nmod_inv(root,mod);
    for(slong i=1;i<length;i++) {
        roots[i]=nmod_mul(roots[i-1],root,mod);
        inverse[i]=nmod_mul(inverse[i-1],invroot,mod);
    }
    for(int k=0;k<2;k++) for(slong j=0;j<dr_mpoly_length(polys+k);j++) {
        DR_MPOLY_TERM(term,polys+k,j);
        slong x=term.var_exp ? term.var_exp[0] : 0;
        slong t=term.par_exp ? term.par_exp[0] : 0;
        slong row=x+(k ? dx[0]+1 : 0);
        values[(size_t)row*length+t]=nmod_add(values[(size_t)row*length+t],
            nmod_poly_get_coeff_ui(term.coeff,0),mod);
    }
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) if(!omp_in_parallel())
#endif
    for(slong i=0;i<rows;i++) bivariate_ntt(values+(size_t)i*length,length,roots,mod);
    double evaluated=get_wall_time();
#ifdef _OPENMP
    #pragma omp parallel if(!omp_in_parallel())
#endif
    {
        nmod_poly_t f,g;nmod_poly_init(f,prime);nmod_poly_init(g,prime);
        nmod_poly_fit_length(f,dx[0]+1);nmod_poly_fit_length(g,dx[1]+1);
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for(slong i=0;i<length;i++) {
            for(slong j=0;j<=dx[0];j++) f->coeffs[j]=values[(size_t)j*length+i];
            for(slong j=0;j<=dx[1];j++) g->coeffs[j]=values[(size_t)(dx[0]+1+j)*length+i];
            f->length=dx[0]+1;g->length=dx[1]+1;
            samples[i]=nmod_poly_resultant(f,g);
        }
        nmod_poly_clear(f);nmod_poly_clear(g);
    }
    double computed=get_wall_time();
    bivariate_ntt(samples,length,inverse,mod);
    ulong scale=nmod_inv(length,mod);
    /* A deterministic consistency check before publishing the result. */
    int ok=1;
    for(slong i=bound+1;i<length;i++) if(samples[i]) ok=0;
    if(ok) {
        dr_mpoly_init(result,0,1,polys[0].ctx);
        fq_nmod_t coeff;fq_nmod_init(coeff,polys[0].ctx);
        for(slong i=0;i<=(slong)bound;i++) if(samples[i]) {
            fq_nmod_set_ui(coeff,nmod_mul(samples[i],scale,mod),polys[0].ctx);
            dr_mpoly_add_term_fast(result,NULL,&i,coeff);
        }
        fq_nmod_clear(coeff,polys[0].ctx);
        dr_mpoly_make_monic(result);
        dixon_info_log("Dixon/Bezout fast resultant: NTT, %ld points, resultant degree bound %lu\n",length,bound);
        dixon_info_log("  Evaluation %.6fs | scalar resultants %.6fs | interpolation/output %.6fs\n",
            evaluated-start,computed-evaluated,get_wall_time()-computed);
        dixon_info_log("Dixon/Bezout fast time: Wall time: %.6f seconds\n",get_wall_time()-start);
    }
    flint_free(values);flint_free(roots);flint_free(inverse);flint_free(samples);
    if(ok) return 1;
    reason="interpolated coefficients violate the resultant degree bound";
fallback:
    (void)reason; /* The unified path may still support multipoint interpolation. */
    return 0;
}

static int dixon_try_bezout_fast(unified_mpoly_struct *result,
    unified_mpoly_struct *polys, slong nvars, slong npars)
{
    if((!g_dixon_bezout_fast && !g_dixon_bezout_block) || nvars!=1 || npars!=1 || g_field_equation_reduction) return 0;
    /* Parsing and variable permutation can leave native terms unsorted.
     * FLINT degree queries (and the Bezout kernels) require canonical input. */
    dr_mpoly_normalize(polys);
    dr_mpoly_normalize(polys+1);
    unified_mpoly_t temp=unified_mpoly_init(polys[0].ctx_ptr);
    int ok=g_dixon_bezout_block && unified_mpoly_resultant_bezout_block_try(temp,polys,polys+1,0,polys[0].ctx_ptr);
    if(!ok) {
        if(g_dixon_bezout_block)dixon_info_log("Dixon/Bezout block-Hankel: fallback to existing resultant backend\n");
        if(dixon_try_bivariate_ntt(result,polys,nvars,npars)) {unified_mpoly_clear(temp);return 1;}
        ok=unified_mpoly_resultant_bezout_try(temp,polys,polys+1,0,polys[0].ctx_ptr);
    }
    if(ok) {
        dr_mpoly_init(result,0,1,polys[0].ctx);
        field_ctx_t *field=polys[0].ctx_ptr->field_ctx;
        field_elem_u c;field_init_elem(&c,field->field_id,field_ctx_data(field));
        fq_nmod_t coeff;fq_nmod_init(coeff,polys[0].ctx);
        for(slong i=0;i<unified_mpoly_length(temp);i++) {
            ulong exp[2];unified_mpoly_get_term_exp_ui(exp,temp,i,temp->ctx_ptr);
            unified_mpoly_get_term_coeff_ui(&c,temp,i,temp->ctx_ptr);
            if(field->field_id==FIELD_ID_NMOD) fq_nmod_set_ui(coeff,c.nmod,polys[0].ctx);
            else field_elem_to_fq_nmod(coeff,&c,field);
            slong power=exp[1];dr_mpoly_add_term_fast(result,NULL,&power,coeff);
        }
        fq_nmod_clear(coeff,polys[0].ctx);
        field_clear_elem(&c,field->field_id,field_ctx_data(field));
        dr_mpoly_make_monic(result);
    } else dixon_debug_log("Dixon/Bezout fast resultant: unsupported input, insufficient field points or size overflow; using Dixon matrix\n");
    unified_mpoly_clear(temp);
    return ok;
}

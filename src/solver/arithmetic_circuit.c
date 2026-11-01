/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "arithmetic_circuit.h"
#include "fq_nmod_roots.h"
#include "dr_profile.h"
#include <errno.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#include <limits.h>
#include <string.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef enum {C_X,C_Y,C_CONST,C_ADD,C_SUB,C_MUL,C_POW,C_LIN} circuit_op;
typedef struct {
    circuit_op op;
    slong a,b,nargs,*refs,degree,slot;
    ulong constant,exponent,*coeffs,leading;
    int affine;
    ulong low[16];
    int low_exact;
    ulong form[3]; /* exact constant, x, y for affine subcircuits */
} circuit_node;
struct dr_circuit {
    circuit_node *nodes;
    slong count,capacity,outputs[2],slots,degree;
    int eliminate;
    nmod_t mod;
};
static int natural(const char *s,ulong *value)
{
    if(!s || !*s)return 0;
    ulong v=0;
    for(;*s;s++) {
        if(*s<'0' || *s>'9')return 0;
        ulong d=*s-'0';if(v>(UWORD_MAX-d)/10)return 0;v=10*v+d;
    }
    *value=v;return 1;
}
static int coefficient(const char *s,ulong *v,ulong prime)
{
    if(!s || !*s)return 0;
    const char *p=s;if(*p=='-' || *p=='+')p++;
    if(!*p)return 0;
    for(;*p;p++)if(!isdigit((unsigned char)*p))return 0;
    fmpz_t z;fmpz_init(z);
    int ok=!fmpz_set_str(z,*s=='+'?s+1:s,10);
    if(ok)*v=fmpz_fdiv_ui(z,prime);
    fmpz_clear(z);return ok;
}
static int reference(const char *s,slong count,slong *out)
{
    ulong v;if(!natural(s,&v) || v>=(ulong)count)return 0;*out=(slong)v;return 1;
}
static void add_node(dr_circuit *c,circuit_node node)
{
    if(c->count==c->capacity) {
        slong capacity=c->capacity?2*c->capacity:16;
        if(capacity<c->capacity || (ulong)capacity>SIZE_MAX/sizeof(node))flint_abort();
        c->nodes=flint_realloc(c->nodes,capacity*sizeof(node));c->capacity=capacity;
    }
    c->nodes[c->count++]=node;
}
void dr_circuit_clear(dr_circuit *c)
{
    if(!c)return;
    for(slong i=0;i<c->count;i++){flint_free(c->nodes[i].refs);flint_free(c->nodes[i].coeffs);}
    flint_free(c->nodes);flint_free(c);
}
/* Conservative total-degree bound and coefficient of the pure top power.
 * A nonzero top coefficient certifies the exact degree and constant leading
 * coefficient without expanding any bivariate polynomial. */
static int analyze_node(dr_circuit *c,circuit_node *n)
{
    circuit_node *a=c->nodes+n->a,*b=c->nodes+n->b;
    /* Exact affine folding catches cancellation in input parameterizations.
     * It never infers a degree bound from numerical samples. */
    n->affine=1;
    switch(n->op) {
    case C_X:n->form[1]=1;break;
    case C_Y:n->form[2]=1;break;
    case C_CONST:n->form[0]=n->constant;break;
    case C_ADD:case C_SUB:
        n->affine=a->affine && b->affine;
        if(n->affine)for(int j=0;j<3;j++)n->form[j]=n->op==C_ADD?nmod_add(a->form[j],b->form[j],c->mod):nmod_sub(a->form[j],b->form[j],c->mod);
        break;
    case C_MUL:
        n->affine=a->affine && b->affine && (a->degree<=0 || b->degree<=0);
        if(n->affine)for(int j=0;j<3;j++)n->form[j]=a->degree<=0?nmod_mul(a->form[0],b->form[j],c->mod):nmod_mul(b->form[0],a->form[j],c->mod);
        break;
    case C_POW:
        n->affine=!n->exponent || (a->affine && (a->degree<=0 || n->exponent==1));
        if(!n->exponent)n->form[0]=1;
        else if(n->affine) {
            if(n->exponent==1)memcpy(n->form,a->form,sizeof(n->form));
            else n->form[0]=nmod_pow_ui(a->form[0],n->exponent,c->mod);
        }
        break;
    case C_LIN:
        n->form[0]=n->constant;
        for(slong k=0;k<n->nargs;k++) {
            circuit_node *v=c->nodes+n->refs[k];if(!n->coeffs[k])continue;
            if(!v->affine){n->affine=0;break;}
            for(int j=0;j<3;j++)n->form[j]=nmod_add(n->form[j],nmod_mul(n->coeffs[k],v->form[j],c->mod),c->mod);
        }
        break;
    }
    if(n->affine) {
        n->degree=(n->form[1] || n->form[2])?1:n->form[0]?0:-1;
        n->leading=n->degree==1?n->form[c->eliminate?2:1]:n->form[0];
        return 1;
    }
    n->degree=-1;n->leading=0;
    switch(n->op) {
    case C_X:case C_Y:n->degree=1;n->leading=(n->op==C_Y)==c->eliminate;break;
    case C_CONST:if(n->constant){n->degree=0;n->leading=n->constant;}break;
    case C_ADD:case C_SUB:
        n->degree=FLINT_MAX(a->degree,b->degree);
        if(a->degree==n->degree)n->leading=a->leading;
        if(b->degree==n->degree)n->leading=n->op==C_ADD?nmod_add(n->leading,b->leading,c->mod):nmod_sub(n->leading,b->leading,c->mod);
        break;
    case C_MUL:
        if(a->degree<0 || b->degree<0)break;
        if(a->degree>WORD_MAX-b->degree)return 0;
        n->degree=a->degree+b->degree;n->leading=nmod_mul(a->leading,b->leading,c->mod);break;
    case C_POW:
        if(!n->exponent){n->degree=0;n->leading=1;break;}
        if(a->degree<0)break;
        if(a->degree && n->exponent>(ulong)WORD_MAX/(ulong)a->degree)return 0;
        n->degree=(slong)((ulong)a->degree*n->exponent);
        n->leading=nmod_pow_ui(a->leading,n->exponent,c->mod);break;
    case C_LIN:
        if(n->constant){n->degree=0;n->leading=n->constant;}
        for(slong j=0;j<n->nargs;j++) {
            circuit_node *v=c->nodes+n->refs[j];if(!n->coeffs[j] || v->degree<0)continue;
            ulong lead=nmod_mul(v->leading,n->coeffs[j],c->mod);
            if(v->degree>n->degree){n->degree=v->degree;n->leading=lead;}
            else if(v->degree==n->degree)n->leading=nmod_add(n->leading,lead,c->mod);
        }
        break;
    }
    return 1;
}
/* Small-degree constant folding, an optimization rather than an input limit.
 * Above degree three the certified structural bound remains in use. This
 * resolves cancellations in cubic parameterizations without expanding the
 * high-degree system. Index 4*x_degree + y_degree. */
static void low_mul(ulong *out,const ulong *a,const ulong *b,nmod_t mod)
{
    ulong tmp[16]={0};
    for(int i=0;i<=3;i++)for(int j=0;j<=3-i;j++)
        for(int k=0;k<=3-i-j;k++)for(int l=0;l<=3-i-j-k;l++)
            tmp[4*(i+k)+j+l]=nmod_add(tmp[4*(i+k)+j+l],nmod_mul(a[4*i+j],b[4*k+l],mod),mod);
    memcpy(out,tmp,sizeof(tmp));
}
static void fold_low_degree(dr_circuit *c,circuit_node *n)
{
    if(n->degree>3)return;
    circuit_node *a=c->nodes+n->a,*b=c->nodes+n->b;
    n->low_exact=1;
    switch(n->op) {
    case C_X:n->low[4]=1;break;case C_Y:n->low[1]=1;break;
    case C_CONST:n->low[0]=n->constant;break;
    case C_ADD:case C_SUB:
        if(!a->low_exact || !b->low_exact){n->low_exact=0;return;}
        for(int j=0;j<16;j++)n->low[j]=n->op==C_ADD?nmod_add(a->low[j],b->low[j],c->mod):nmod_sub(a->low[j],b->low[j],c->mod);
        break;
    case C_MUL:
        if(n->degree<0)break;
        if(!a->low_exact || !b->low_exact){n->low_exact=0;return;}
        low_mul(n->low,a->low,b->low,c->mod);break;
    case C_POW:
        if(!n->exponent){n->low[0]=1;break;}
        if(a->degree<0)break;
        if(!a->low_exact){n->low_exact=0;return;}
        if(a->degree==0){n->low[0]=nmod_pow_ui(a->low[0],n->exponent,c->mod);break;}
        n->low[0]=1;
        for(ulong j=0;j<n->exponent;j++)low_mul(n->low,n->low,a->low,c->mod);
        break;
    case C_LIN:
        n->low[0]=n->constant;
        for(slong k=0;k<n->nargs;k++) {
            circuit_node *v=c->nodes+n->refs[k];if(!n->coeffs[k])continue;
            if(!v->low_exact){n->low_exact=0;return;}
            for(int j=0;j<16;j++)n->low[j]=nmod_add(n->low[j],nmod_mul(n->coeffs[k],v->low[j],c->mod),c->mod);
        }
        break;
    }
    n->degree=-1;
    for(int i=0;i<=3;i++)for(int j=0;j<=3-i;j++)if(n->low[4*i+j])n->degree=FLINT_MAX(n->degree,i+j);
    n->leading=n->degree<0?0:n->low[(c->eliminate?1:4)*n->degree];
    n->affine=n->degree<=1;
    if(n->affine){n->form[0]=n->low[0];n->form[1]=n->low[4];n->form[2]=n->low[1];}
}
static slong arity(const circuit_node *n) {
    return n->op==C_LIN?n->nargs:n->op==C_POW?1:
        n->op==C_ADD || n->op==C_SUB || n->op==C_MUL?2:0;
}
static slong operand(const circuit_node *n,slong i) {return n->op==C_LIN?n->refs[i]:i?n->b:n->a;}
static void assign_slots(dr_circuit *c)
{
    slong *last=flint_malloc(c->count*sizeof(*last)),*free_slots=flint_malloc(c->count*sizeof(*free_slots)),nf=0;
    for(slong i=0;i<c->count;i++)last[i]=i;
    for(slong i=0;i<c->count;i++)for(slong j=0;j<arity(c->nodes+i);j++)last[operand(c->nodes+i,j)]=i;
    last[c->outputs[0]]=last[c->outputs[1]]=c->count;
    for(slong i=0;i<c->count;i++) {
        c->nodes[i].slot=nf?free_slots[--nf]:c->slots++;
        for(slong j=0;j<arity(c->nodes+i);j++) {
            slong ref=operand(c->nodes+i,j);
            if(last[ref]==i){free_slots[nf++]=c->nodes[ref].slot;last[ref]=-1;}
        }
        if(last[i]==i)free_slots[nf++]=c->nodes[i].slot;
    }
    flint_free(last);flint_free(free_slots);
}
dr_circuit *dr_circuit_read(FILE *input,char *error,size_t error_size)
{
    dr_circuit *c=flint_calloc(1,sizeof(*c));c->eliminate=1;
    char *line=NULL;size_t cap=0;slong lineno=0;int header=0,prime=0,outputs=0;
    const char *why="invalid circuit syntax";
    while(getline(&line,&cap,input)>=0) {
        lineno++;char *comment=strchr(line,'#');if(comment)*comment=0;
        char *save=NULL,*op=strtok_r(line," \t\r\n",&save);if(!op)continue;
        if(outputs){why="unexpected data after outputs";goto fail;}
        char **tokens=NULL;slong nt=0,alloc=0;
        for(char *t=strtok_r(NULL," \t\r\n",&save);t;t=strtok_r(NULL," \t\r\n",&save)) {
            if(nt==alloc){alloc=alloc?2*alloc:8;tokens=flint_realloc(tokens,alloc*sizeof(*tokens));}
            tokens[nt++]=t;
        }
        int ok=1; circuit_node node={0};node.a=node.b=0;
        if(!header) {ok=!strcmp(op,"DRCIRCUIT") && nt==1 && !strcmp(tokens[0],"1");header=ok;}
        else if(!prime) {
            ulong p;ok=!strcmp(op,"prime") && nt==1 && natural(tokens[0],&p) && p>=2 && n_is_prime(p);
            if(ok){nmod_init(&c->mod,p);prime=1;node.op=C_X;add_node(c,node);node.op=C_Y;add_node(c,node);}
        } else if(!strcmp(op,"eliminate")) {
            ok=c->count==2 && nt==1 && (!strcmp(tokens[0],"x") || !strcmp(tokens[0],"y"));
            if(ok)c->eliminate=!strcmp(tokens[0],"y");
        } else if(!strcmp(op,"outputs")) {
            ok=nt==2 && reference(tokens[0],c->count,c->outputs) && reference(tokens[1],c->count,c->outputs+1);outputs=ok;
        } else {
            if(!strcmp(op,"const")){node.op=C_CONST;ok=nt==1 && coefficient(tokens[0],&node.constant,c->mod.n);}
            else if(!strcmp(op,"add") || !strcmp(op,"sub") || !strcmp(op,"mul")) {
                node.op=!strcmp(op,"add")?C_ADD:!strcmp(op,"sub")?C_SUB:C_MUL;
                ok=nt==2 && reference(tokens[0],c->count,&node.a) && reference(tokens[1],c->count,&node.b);
            } else if(!strcmp(op,"pow")) {
                node.op=C_POW;ok=nt==2 && reference(tokens[0],c->count,&node.a) && natural(tokens[1],&node.exponent);
            } else if(!strcmp(op,"lin")) {
                node.op=C_LIN;ok=nt>=1 && (nt&1) && coefficient(tokens[0],&node.constant,c->mod.n);
                if(ok) {
                    node.nargs=(nt-1)/2;node.refs=flint_malloc(node.nargs*sizeof(*node.refs));node.coeffs=flint_malloc(node.nargs*sizeof(*node.coeffs));
                    for(slong j=0;ok && j<node.nargs;j++)ok=coefficient(tokens[1+2*j],node.coeffs+j,c->mod.n) && reference(tokens[2+2*j],c->count,node.refs+j);
                }
            } else ok=0;
            if(ok)add_node(c,node);else{flint_free(node.refs);flint_free(node.coeffs);}
        }
        flint_free(tokens);if(!ok)goto fail;
    }
    if(ferror(input)){why="input read failed";goto fail;}
    if(!outputs){why="missing outputs";goto fail;}
    for(slong i=0;i<c->count;i++) {
        if(!analyze_node(c,c->nodes+i)){why="circuit degree exceeds slong representation";goto fail;}
        fold_low_degree(c,c->nodes+i);
    }
    circuit_node *a=c->nodes+c->outputs[0],*b=c->nodes+c->outputs[1];
    if(a->degree<1 || a->degree!=b->degree || !a->leading || !b->leading) {
        why="unsupported circuit: need certified equal positive total degrees and constant nonzero leading coefficients; try the expanded .dr interface";goto fail;
    }
    c->degree=a->degree;assign_slots(c);free(line);return c;
fail:
    snprintf(error,error_size,"line %ld: %s",lineno,why);free(line);dr_circuit_clear(c);return NULL;
}
ulong dr_circuit_prime(const dr_circuit *c){return c->mod.n;}
slong dr_circuit_degree(const dr_circuit *c){return c->degree;}
static void *workspace_init(const void *data)
{
    const dr_circuit *c=data;nmod_poly_struct *w=flint_malloc(c->slots*sizeof(*w));
    for(slong i=0;i<c->slots;i++)nmod_poly_init(w+i,c->mod.n);return w;
}
static void workspace_clear(void *workspace,const void *data)
{
    const dr_circuit *c=data;nmod_poly_struct *w=workspace;
    for(slong i=0;i<c->slots;i++)nmod_poly_clear(w+i);flint_free(w);
}
static void poly_constant(nmod_poly_t p,ulong value)
{
    nmod_poly_zero(p);nmod_poly_set_coeff_ui(p,0,value);
}
static void specialize(nmod_poly_t f,nmod_poly_t g,ulong point,void *workspace,const void *data)
{
    const dr_circuit *c=data;nmod_poly_struct *w=workspace;
    for(slong i=0;i<c->count;i++) {
        const circuit_node *n=c->nodes+i;nmod_poly_struct *r=w+n->slot,*a=w+c->nodes[n->a].slot,*b=w+c->nodes[n->b].slot;
        switch(n->op) {
        case C_X:case C_Y:
            if((n->op==C_Y)==c->eliminate){nmod_poly_zero(r);nmod_poly_set_coeff_ui(r,1,1);}
            else poly_constant(r,point);break;
        case C_CONST:poly_constant(r,n->constant);break;
        case C_ADD:nmod_poly_add(r,a,b);break;
        case C_SUB:nmod_poly_sub(r,a,b);break;
        case C_MUL:nmod_poly_mul(r,a,b);break;
        case C_POW:nmod_poly_pow(r,a,n->exponent);break;
        case C_LIN:
            poly_constant(r,n->constant);
            for(slong j=0;j<n->nargs;j++)nmod_poly_scalar_addmul_nmod(r,w+c->nodes[n->refs[j]].slot,n->coeffs[j]);break;
        }
    }
    nmod_poly_set(f,w+c->nodes[c->outputs[0]].slot);nmod_poly_set(g,w+c->nodes[c->outputs[1]].slot);
}
static int specialize_gr(gr_poly_t f,gr_poly_t g,gr_srcptr point,gr_ctx_t ctx,const void *data)
{
    const dr_circuit *c=data;gr_poly_struct *w=flint_malloc(c->slots*sizeof(*w));int status=GR_SUCCESS;
    for(slong i=0;i<c->slots;i++)gr_poly_init(w+i,ctx);
    gr_poly_t tmp;gr_poly_init(tmp,ctx);
    for(slong i=0;i<c->count && status==GR_SUCCESS;i++) {
        const circuit_node *n=c->nodes+i;gr_poly_struct *r=w+n->slot,*a=w+c->nodes[n->a].slot,*b=w+c->nodes[n->b].slot;
        switch(n->op) {
        case C_X:case C_Y:
            status|=gr_poly_zero(r,ctx);
            if((n->op==C_Y)==c->eliminate)status|=gr_poly_set_coeff_ui(r,1,1,ctx);
            else status|=gr_poly_set_coeff_scalar(r,0,point,ctx);break;
        case C_CONST:status|=gr_poly_set_ui(r,n->constant,ctx);break;
        case C_ADD:status|=gr_poly_add(r,a,b,ctx);break;
        case C_SUB:status|=gr_poly_sub(r,a,b,ctx);break;
        case C_MUL:status|=gr_poly_mul(r,a,b,ctx);break;
        case C_POW:status|=gr_poly_pow_ui(r,a,n->exponent,ctx);break;
        case C_LIN:
            status|=gr_poly_set_ui(r,n->constant,ctx);
            for(slong j=0;j<n->nargs;j++){status|=gr_poly_mul_ui(tmp,w+c->nodes[n->refs[j]].slot,n->coeffs[j],ctx);status|=gr_poly_add(r,r,tmp,ctx);}break;
        }
    }
    status|=gr_poly_set(f,w+c->nodes[c->outputs[0]].slot,ctx);status|=gr_poly_set(g,w+c->nodes[c->outputs[1]].slot,ctx);
    gr_poly_clear(tmp,ctx);for(slong i=0;i<c->slots;i++)gr_poly_clear(w+i,ctx);flint_free(w);return status;
}
void dr_circuit_source(dr_bezout_source *source,const dr_circuit *c)
{
    *source=(dr_bezout_source){c,workspace_init,workspace_clear,specialize,specialize_gr};
}
static int scalar_check(const dr_circuit *c,ulong x,ulong y)
{
    ulong *w=flint_calloc(c->slots,sizeof(*w));
    for(slong i=0;i<c->count;i++) {
        const circuit_node *n=c->nodes+i;ulong *r=w+n->slot,a=w[c->nodes[n->a].slot],b=w[c->nodes[n->b].slot];
        switch(n->op) {
        case C_X:*r=x;break;case C_Y:*r=y;break;case C_CONST:*r=n->constant;break;
        case C_ADD:*r=nmod_add(a,b,c->mod);break;case C_SUB:*r=nmod_sub(a,b,c->mod);break;
        case C_MUL:*r=nmod_mul(a,b,c->mod);break;case C_POW:*r=nmod_pow_ui(a,n->exponent,c->mod);break;
        case C_LIN:*r=n->constant;for(slong j=0;j<n->nargs;j++)*r=nmod_add(*r,nmod_mul(n->coeffs[j],w[c->nodes[n->refs[j]].slot],c->mod),c->mod);break;
        }
    }
    int ok=!w[c->nodes[c->outputs[0]].slot] && !w[c->nodes[c->outputs[1]].slot];flint_free(w);return ok;
}

int dr_circuit_cli_main(int argc,char **argv)
{
    double start=dr_profile_now();const char *input=NULL,*output=NULL;int only=0,threads=1,verbose=1;
#ifdef _OPENMP
    threads=(omp_get_max_threads()+1)/2;
#endif
    for(int i=1;i<argc;i++) {
        if(!strcmp(argv[i],"--circuit") && i+1<argc && !input)input=argv[++i];
        else if(!strcmp(argv[i],"-o") && i+1<argc && !output)output=argv[++i];
        else if(!strcmp(argv[i],"--resultant-only"))only=1;
        else if(!strcmp(argv[i],"--time")){} /* timings always reported here */
        else if((!strcmp(argv[i],"--threads") || !strcmp(argv[i],"--bezout-initial-points") || !strcmp(argv[i],"-v")) && i+1<argc) {
            const char *key=argv[i];ulong v;
            if(!natural(argv[++i],&v) || v>(ulong)WORD_MAX)goto usage;
            if(!strcmp(key,"--threads")){if(!v || v>INT_MAX)goto usage;threads=(int)v;}
            else if(!strcmp(key,"-v")){if(v>2)goto usage;verbose=(int)v;}
            else dr_bezout_block_set_initial_points((slong)v);
        } else goto usage;
    }
    if(!input || (output && !strcmp(input,output)))goto usage;
    FILE *in=fopen(input,"r");if(!in){perror(input);return 1;}
    char error[512];dr_circuit *c=dr_circuit_read(in,error,sizeof(error));fclose(in);
    if(!c){fprintf(stderr,"Circuit error: %s\n",error);return 1;}
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    flint_set_num_threads(threads);
    if(verbose)printf("Circuit: prime=%lu nodes=%ld workspace_slots=%ld degree=%ld eliminate=%c threads=%d\n",c->mod.n,c->count,c->slots,c->degree,c->eliminate?'y':'x',threads);
    double parsed=dr_profile_now();
    dr_bezout_source source;dr_circuit_source(&source,c);
    nmod_poly_t resultant;nmod_poly_init(resultant,c->mod.n);
    flint_rand_t state;flint_rand_init(state);dr_bezout_block_stats stats;
    int ok=dr_nmod_resultant_bezout_source_try(resultant,c->degree,&source,state,&stats);
    flint_rand_clear(state);double eliminated=dr_profile_now();
    if(!ok) {
        fprintf(stderr,"Circuit resultant could not be reconstructed and verified (including field-point availability). No solution claim is made.\n");
        nmod_poly_clear(resultant);dr_circuit_clear(c);flint_cleanup_master();return 1;
    }
    FILE *out=output?fopen(output,"w"):NULL;
    if(output && !out){perror(output);nmod_poly_clear(resultant);dr_circuit_clear(c);flint_cleanup_master();return 1;}
    if(verbose)printf("Circuit resultant: degree=%ld block=%ld points=%ld attempts=%ld\nSampling+inverse %.6fs | interpolation %.6fs | determinant %.6fs | verification %.6fs (field degree %ld, %ld checks)\n",
        nmod_poly_degree(resultant),stats.block,stats.points,stats.attempts,stats.inverse,stats.interpolation,stats.determinant,stats.verification,stats.verification_extension_degree,stats.verification_checks);
    slong solutions=0;double roots_time=0,back_time=0;
    if(only) {
        FILE *dest=out?out:stdout;fprintf(dest,"Resultant:\n");nmod_poly_fprint_pretty(dest,resultant,c->eliminate?"x":"y");fputc('\n',dest);
    } else {
        nmod_roots_t roots;nmod_roots_init(roots);double t=dr_profile_now();
        our_nmod_poly_roots(roots,resultant,1);roots_time=dr_profile_now()-t;
        t=dr_profile_now();void *w=workspace_init(c);
        nmod_poly_t f,g,h;nmod_poly_init(f,c->mod.n);nmod_poly_init(g,c->mod.n);nmod_poly_init(h,c->mod.n);
        for(slong i=0;i<roots->num;i++) {
            specialize(f,g,roots->roots[i],w,c);nmod_poly_gcd(h,f,g);
            nmod_roots_t second;nmod_roots_init(second);our_nmod_poly_roots(second,h,1);
            for(slong j=0;j<second->num;j++) {
                ulong x=c->eliminate?roots->roots[i]:second->roots[j],y=c->eliminate?second->roots[j]:roots->roots[i];
                if(!scalar_check(c,x,y)){fprintf(stderr,"Circuit scalar verification failed\n");ok=0;continue;}
                solutions++;
                if(verbose)printf("[%ld] x = %lu, y = %lu [PASS]\n",solutions,x,y);
                if(out)fprintf(out,"[%ld] x = %lu, y = %lu [PASS]\n",solutions,x,y);
            }
            nmod_roots_clear(second);
        }
        nmod_poly_clear(f);nmod_poly_clear(g);nmod_poly_clear(h);workspace_clear(w,c);nmod_roots_clear(roots);
        back_time=dr_profile_now()-t;
        if(!solutions && ok){if(verbose)puts("No solutions over the finite field.");if(out)fputs("No solutions over the finite field.\n",out);}
    }
    if(out){int failed=ferror(out);if(fclose(out))failed=1;if(failed)ok=0;}
    nmod_poly_clear(resultant);dr_circuit_clear(c);flint_cleanup_master();
    if(verbose)printf("Circuit timing: parse=%.6f resultant=%.6f roots=%.6f back_substitution=%.6f total=%.6f\n",parsed-start,eliminated-parsed,roots_time,back_time,dr_profile_now()-start);
    return !ok;
usage:
    fprintf(stderr,"Usage: %s --circuit input.drc [--threads N] [--bezout-initial-points N] [--resultant-only] [--time] [-v 0|1|2] [-o output]\n",argv[0]);return 1;
}

#include "../../determinant/mq_quotient_impl.h"
/* Exercise tile boundaries, signs, zeros and the IFMA product-width guard
 * against FLINT's independent modular matrix multiplication. */
static void check_packed_tiles(void)
{
    const ulong moduli[]={37,65537,67108859,67108879,536870923};
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,9876,54321);
    for(size_t p=0;p<sizeof(moduli)/sizeof(*moduli);p++) {
        prime=moduli[p];nmod_init(&mod,prime);
        nmod_mat_t a,b,expected;nmod_mat_init(a,43,47,prime);
        nmod_mat_init(b,47,53,prime);nmod_mat_init(expected,43,53,prime);
        for(slong i=0;i<a->r;i++)for(slong j=0;j<a->c;j++)
            nmod_mat_entry(a,i,j)=((i+j)%7==0)?0:prime-1-n_randint(rng,FLINT_MIN(prime,17));
        for(slong i=0;i<b->r;i++)for(slong j=0;j<b->c;j++)
            nmod_mat_entry(b,i,j)=prime-1-n_randint(rng,FLINT_MIN(prime,17));
        nmod_mat_mul(expected,a,b);
        packed_matrix pa,pb;pack_matrix(&pa,a);pack_matrix(&pb,b);
        ulong *acc=flint_calloc(43*53,sizeof(ulong));
        for(int kernel=4;kernel<=5;kernel++)for(int neg=0;neg<2;neg++) {
            for(slong i=0;i<43*53;i++)acc[i]=3;
            recurrence_kernel=kernel;packed_mat_addmul(acc,&pa,&pb,neg);
            for(slong i=0;i<43;i++)for(slong j=0;j<53;j++) {
                ulong c=nmod_mat_entry(expected,i,j);if(neg)c=nmod_neg(c,mod);
                assert(nmod_set_ui(acc[i*53+j],mod)==nmod_add(c,3,mod));
            }
        }
        flint_free(acc);flint_free(pa.data);flint_free(pb.data);nmod_mat_clear(expected);
    }
    flint_rand_clear(rng);
}
static void check_skipped_pivot(void)
{
    assert(quotient_dimensions(4));prime=37;nmod_init(&mod,prime);
    ulong q[64]={0};
    /* Pure-x relations x0*x1, x0*x2, x1^2, x1*x2 leave column x0^2 free. */
    const slong ii[]={0,0,1,1},jj[]={1,2,1,2};
    for(slong f=0;f<N;f++)for(slong j=0;j<N;j++)
        q[(f*N+ii[j])*N+jj[j]]=(f==j?2:1);
    assert(leading_init(q));assert(frames[2].pivot[0]==1);
    for(slong sample=0;sample<3;sample++) {
        for(slong f=0;f<N;f++)for(slong i=0;i<N;i++)
            q[(f*N+i)*N+N-1]=(sample*11+f*7+i*3)%prime;
        grade *cached=flint_calloc(L+1,sizeof(grade));
        grade *full=flint_calloc(L+1,sizeof(grade));
        frames_enabled=1;assert(normal_forms(cached,q));
        frames_enabled=0;assert(normal_forms(full,q));
        for(slong d=0;d<=L;d++) {
            assert(!memcmp(cached[d].free,full[d].free,cached[d].h*sizeof(slong)));
            assert(nmod_mat_equal(cached[d].reduction,full[d].reduction));
        }
        clear_grades(cached);clear_grades(full);flint_free(cached);flint_free(full);
    }
    leading_clear();normal_time=map_time=dp_time=0;
}
int main(int argc,char **argv)
{
    check_packed_tiles();check_skipped_pivot();
    assert(argc>=2);FILE *in=fopen(argv[1],"r");assert(in);
    assert(fscanf(in,"%ld %lu",&N,&prime)==2);assert(quotient_dimensions(N));
    M=N/2;K=N-M;L=K;threads=argc>2?atol(argv[2]):1;int mode=argc>3?atoi(argv[3]):0;
    recurrence_kernel=argc>4?atoi(argv[4]):-1;
    if(recurrence_kernel==-1)recurrence_kernel=packed_admitted()?5:1;
    assert(recurrence_kernel>=0 && recurrence_kernel<=6);
    if(recurrence_kernel>=3)assert(packed_admitted());
    int interpolate=mode!=0;
    flint_set_num_threads(1);
#ifdef _OPENMP
    omp_set_num_threads(threads);
#endif
    nmod_init(&mod,prime);
    slong w=N+1;ulong *input=flint_calloc(N*w*w,sizeof(ulong));
    for(slong f=0;f<N;f++)for(slong i=0;i<w;i++)for(slong j=i;j<w;j++)assert(fscanf(in,"%lu",input+(f*w+i)*w+j)==1);
    nmod_poly_t reference,det;nmod_poly_init(reference,prime);nmod_poly_init(det,prime);
    slong refdegree;int has_reference=fscanf(in,"%ld",&refdegree)==1;
    if(has_reference)for(slong i=0;i<=refdegree;i++){ulong c;assert(fscanf(in,"%lu",&c)==1);nmod_poly_set_coeff_ui(reference,i,c);}
    fclose(in);slong h=choose(N,M),samples=interpolate?N+2:1;
    assert((ulong)(samples+2)<prime);ulong *xs=flint_malloc(samples*sizeof(ulong)),*ys=flint_malloc(h*h*samples*sizeof(ulong));
    uint64_t *basis=flint_malloc(h*sizeof(uint64_t)),*basis_y=flint_malloc(h*sizeof(uint64_t));
    nmod_poly_mat_t polynomial;nmod_poly_mat_init(polynomial,h,h,prime);
    double total=now(),interpolation_time=0,det_time=0,leading_time=0;
    if(mode==3) {
        ulong *q=flint_calloc(N*N*N,sizeof(ulong));specialize(q,input,0);
        flint_set_num_threads(argc>5?atol(argv[5]):1);
        double t=now();assert(leading_init(q));leading_time=now()-t;flint_set_num_threads(1);flint_free(q);frames_enabled=1;
    }
    slong start=1;
    if(mode==3) {
        leading_basis(basis,M);leading_basis(basis_y,K);start=0;
    } else {
    nmod_mat_t first;nmod_mat_init(first,h,h,prime);evaluate(first,input,0,basis,basis_y,1);xs[0]=0;
    for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)ys[(i*h+j)*samples]=nmod_mat_entry(first,i,j);
    nmod_mat_clear(first);
    }
    #pragma omp parallel for schedule(dynamic) num_threads(threads) if(mode>=2 && threads>1)
    for(slong sample=start;sample<samples;sample++) {
        nmod_mat_t core;nmod_mat_init(core,h,h,prime);evaluate(core,input,sample,basis,basis_y,0);xs[sample]=sample;
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)ys[(i*h+j)*samples+sample]=nmod_mat_entry(core,i,j);
        nmod_mat_clear(core);
    }
    int equal=-1;
    if(interpolate) {
        double t=now();for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
            nmod_poly_struct *entry=nmod_poly_mat_entry(polynomial,i,j);
            nmod_poly_interpolate_nmod_vec(entry,xs,ys+(i*h+j)*samples,samples);
            slong wi=M-(basis[i]>>(key_bits*(N-1))),wj=K-(basis_y[j]>>(key_bits*(N-1)));
            assert(nmod_poly_degree(entry)<=N+1-wi-wj);
        }
        interpolation_time=now()-t;t=now();assert(nmod_poly_mat_det_hnf(det,polynomial));det_time=now()-t;
        nmod_poly_make_monic(det,det);if(has_reference){equal=nmod_poly_equal(det,reference);assert(equal);}
    }
    double elapsed=now()-total,nf=normal_time,mp=map_time,dp=dp_time,validation=now();
    verify_recurrence=1;
    if(interpolate)for(slong sample=samples;sample<samples+2;sample++) {
        nmod_mat_t core;nmod_mat_init(core,h,h,prime);evaluate(core,input,sample,basis,basis_y,0);
        for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)assert(nmod_poly_evaluate_nmod(nmod_poly_mat_entry(polynomial,i,j),sample)==nmod_mat_entry(core,i,j));
        nmod_mat_clear(core);
    }
    if(frames_enabled) {
        ulong *q=flint_calloc(N*N*N,sizeof(ulong));specialize(q,input,samples+2);
        grade *cached=flint_calloc(L+1,sizeof(grade)),*full=flint_calloc(L+1,sizeof(grade));assert(normal_forms(cached,q));
        frames_enabled=0;assert(normal_forms(full,q));frames_enabled=1;
        for(slong d=0;d<=L;d++) {
            assert(!memcmp(cached[d].free,full[d].free,cached[d].h*sizeof(slong)));
            assert(nmod_mat_equal(cached[d].reduction,full[d].reduction));
        }
        clear_grades(cached);clear_grades(full);flint_free(cached);flint_free(full);flint_free(q);
    }
    validation=now()-validation;
    fprintf(stderr,"constant leading frames: %.6f seconds\n",leading_time);
    printf("{\"n\":%ld,\"prime\":%lu,\"threads\":%ld,\"mode\":%d,\"kernel\":%d,\"samples\":%ld,\"h\":%ld,\"normal_forms_work\":%.6f,\"multiplication_maps_work\":%.6f,\"quotient_dp_work\":%.6f,\"interpolation\":%.6f,\"determinant\":%.6f,\"compute_seconds\":%.6f,\"validation_seconds\":%.6f,\"degree\":%ld,\"reference_equal\":%d}\n",N,prime,threads,mode,recurrence_kernel,samples,h,nf,mp,dp,interpolation_time,det_time,elapsed,validation,nmod_poly_degree(det),equal);
    flint_free(input);flint_free(xs);flint_free(ys);flint_free(basis);flint_free(basis_y);
    nmod_poly_mat_clear(polynomial);nmod_poly_clear(reference);nmod_poly_clear(det);if(frames_enabled)leading_clear();flint_cleanup_master();return 0;
}

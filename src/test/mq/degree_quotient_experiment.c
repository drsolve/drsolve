#include "determinant/degree_quotient_impl.h"

static void write_input(const char *prefix) {
    char path[1024];snprintf(path,sizeof(path),"%s.dr",prefix);FILE *f=fopen(path,"w");assert(f);
    for(slong v=0;v<n-1;v++)fprintf(f,"%sx%ld",v?",":"",v);fprintf(f,"\n%lu\n",prime);
    for(slong eq=0;eq<n;eq++) {
        int first=1;
        for(slong t=0;t<input_count;t++)if(input[t].coeff[eq]) {
            fprintf(f,"%s%lu",first?"":"+",input[t].coeff[eq]);first=0;
            for(slong v=0;v<n;v++)if(exponent(input[t].key,v))fprintf(f,"*x%ld^%ld",v,exponent(input[t].key,v));
        }
        fprintf(f,"\n");
    }
    fclose(f);
}
int main(int argc,char **argv) {
    if(argc<5 || argc>8)fail("usage: N degrees seed output-prefix [prime] [memory-MiB] [grouped|reference]");
    if(argc==8) {
        if(!strcmp(argv[7],"reference"))grouped_kernel=0;
        else if(strcmp(argv[7],"grouped"))fail("kernel must be grouped or reference");
    }
    n=positive_integer(argv[1]);
    if(n<2 || n>8)fail("N must be between 2 and 8 for packed monomial keys");
    char *end;errno=0;ulong seed=strtoul(argv[3],&end,10);
    if(errno || !*argv[3] || *end || argv[3][0]=='-')fail("invalid seed");
    prime=argc>=6?(ulong)positive_integer(argv[5]):65537;
    if(!n_is_prime(prime))fail("coefficient modulus must be prime");
    init_dimensions(argv[2],argc>=7?positive_integer(argv[6]):512);
    slong samples=socle+3;
    if(prime<=(ulong)(samples+2))fail("prime too small for distinct interpolation and validation points");
    nmod_init(&mod,prime);flint_set_num_threads(1);
    init_spaces();input=flint_calloc(binomial(n+max_degree,max_degree),sizeof(*input));
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,seed,seed+1);input_enum(0,max_degree,0,rng);write_input(argv[4]);
    slong h=spaces[left_degree].h;assert(h==spaces[right_degree].h);
    ulong **q=flint_malloc(n*sizeof(*q));
    for(slong f=0;f<n;f++)q[f]=flint_malloc(spaces[degrees[f]].size*sizeof(ulong));
    ulong *xs=flint_malloc(samples*sizeof(ulong)),*ys=flint_malloc(h*h*samples*sizeof(ulong));
    double normal_seconds=0,reduced_seconds=0,direct_seconds=0,setup_seconds=0,update_seconds=0;
    slong products_total=0;
    nmod_mat_t core,direct;nmod_mat_init(core,h,h,prime);nmod_mat_init(direct,h,h,prime);
    for(slong s=0;s<samples+2;s++) {
        double t=seconds();specialize(q,s);if(!normal_forms(q,s==0))fail("non-generic Macaulay rank or changing quotient basis");double normal=seconds()-t;
        entry *entries=difference_matrix(q);
        t=seconds();bezout_core(core,entries,1);double reduced_time=seconds()-t;
        t=seconds();bezout_core(direct,entries,0);double direct_time=seconds()-t;
        assert(nmod_mat_equal(core,direct));
        if(s<samples) {
            setup_seconds+=operator_setup_seconds;update_seconds+=operator_update_seconds;products_total+=matrix_products;
            normal_seconds+=normal;reduced_seconds+=reduced_time;direct_seconds+=direct_time;xs[s]=s;
            for(slong i=0;i<h;i++)for(slong j=0;j<h;j++)ys[(i*h+j)*samples+s]=nmod_mat_entry(core,i,j);
        }
        for(slong i=0;i<n*n;i++)flint_free(entries[i].terms);flint_free(entries);
        if(s>=samples) {
            nmod_poly_t p;nmod_poly_init(p,prime);
            for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
                nmod_poly_interpolate_nmod_vec(p,xs,ys+(i*h+j)*samples,samples);
                if(nmod_poly_evaluate_nmod(p,s)!=nmod_mat_entry(core,i,j))fail("interpolation failed at an extra validation point");
            }
            nmod_poly_clear(p);
        }
    }
    double t=seconds();nmod_poly_mat_t matrix;nmod_poly_mat_init(matrix,h,h,prime);
    for(slong i=0;i<h;i++)for(slong j=0;j<h;j++) {
        nmod_poly_struct *entry=nmod_poly_mat_entry(matrix,i,j);
        nmod_poly_interpolate_nmod_vec(entry,xs,ys+(i*h+j)*samples,samples);
        slong wi=left_degree-exponent(spaces[left_degree].keys[spaces[left_degree].free[i]],n-1);
        slong wj=right_degree-exponent(spaces[right_degree].keys[spaces[right_degree].free[j]],n-1);
        if(nmod_poly_degree(entry)>socle+1-wi-wj)fail("interpolated core violates the weighted degree bound");
    }
    double interpolation=seconds()-t;t=seconds();nmod_poly_t det;nmod_poly_init(det,prime);
    if(!nmod_poly_mat_det_hnf(det,matrix))fail("HNF determinant failed");
    if(nmod_poly_degree(det)!=expected_degree)fail("resultant degree differs from the product of equation degrees");nmod_poly_make_monic(det,det);double determinant=seconds()-t;
    char path[1024];snprintf(path,sizeof(path),"%s.coeff",argv[4]);FILE *f=fopen(path,"w");assert(f);
    for(slong i=0;i<det->length;i++)fprintf(f,"%lu%c",det->coeffs[i],i+1==det->length?'\n':' ');fclose(f);
    printf("{\"kernel\":\"%s\",\"operator_setup\":%.9f,\"operator_updates\":%.9f,\"update_matrix_products\":%ld,\"degrees\":[",
        grouped_kernel?"grouped":"reference",setup_seconds,update_seconds,products_total);
    for(slong f=0;f<n;f++)printf("%s%ld",f?",":"",degrees[f]);
    printf("],\"n\":%ld,\"degree\":%ld,\"seed\":%lu,\"core\":%ld,\"full_rows\":%ld,\"full_cols\":%ld,\"samples\":%ld,\"normal\":%.9f,\"reduced\":%.9f,\"direct\":%.9f,\"interpolation\":%.9f,\"determinant\":%.9f,\"resultant_degree\":%ld}\n",n,max_degree,seed,h,spaces[left_degree].size,spaces[right_degree].size,samples,normal_seconds,reduced_seconds,direct_seconds,interpolation,determinant,nmod_poly_degree(det));
    nmod_poly_clear(det);nmod_poly_mat_clear(matrix);nmod_mat_clear(core);nmod_mat_clear(direct);
    for(slong f=0;f<n;f++)flint_free(q[f]);
    flint_free(q);flint_free(xs);flint_free(ys);
    flint_free(degrees);flint_free(prefix_degree);flint_free(hilbert_dims);
    for(slong d=0;d<=limit;d++){flint_free(spaces[d].keys);flint_free(spaces[d].free);nmod_mat_clear(spaces[d].reduce);}flint_free(spaces);
    for(slong i=0;i<input_count;i++)flint_free(input[i].coeff);flint_free(input);flint_rand_clear(rng);flint_cleanup_master();return 0;
}

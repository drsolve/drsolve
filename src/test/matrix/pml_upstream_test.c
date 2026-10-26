/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <flint/ulong_extras.h>
#include "nmod_mat_poly.h"
#include "nmod_poly_mat_utils.h"
#include "nmod_poly_mat_forms.h"
int g_dixon_verbose_level=0, g_dixon_debug_mode=0;

/* Original conversion loops, retained here as independent references. */
static void forward(nmod_mat_poly_t out,const nmod_poly_mat_t in,slong n)
{
    slong len=nmod_poly_mat_max_length(in);n=FLINT_MIN(n,len);
    nmod_mat_poly_fit_length(out,n);_nmod_mat_poly_set_length(out,n);
    for(slong k=0;k<n;k++)for(slong i=0;i<in->r;i++)for(slong j=0;j<in->c;j++)
        nmod_mat_poly_entry(out,k,i,j)=nmod_poly_get_coeff_ui(nmod_poly_mat_entry(in,i,j),k);
    if(n<len)_nmod_mat_poly_normalise(out);
}
static void backward(nmod_poly_mat_t out,const nmod_mat_poly_t in,slong n)
{
    n=FLINT_MIN(n,in->length);
    for(slong i=0;i<out->r;i++)for(slong j=0;j<out->c;j++)
        nmod_poly_fit_length(nmod_poly_mat_entry(out,i,j),n);
    for(slong k=0;k<n;k++)for(slong i=0;i<out->r;i++)for(slong j=0;j<out->c;j++)
        nmod_poly_mat_entry(out,i,j)->coeffs[k]=nmod_mat_poly_entry(in,k,i,j);
    for(slong i=0;i<out->r;i++)for(slong j=0;j<out->c;j++) {
        _nmod_poly_set_length(nmod_poly_mat_entry(out,i,j),n);
        _nmod_poly_normalise(nmod_poly_mat_entry(out,i,j));
    }
}
static void pivots(void)
{
    orientation_t orientations[]={ROW_LOWER,ROW_UPPER,COL_LOWER,COL_UPPER};
    nmod_poly_struct vec[3];slong shift[]={2,-1,5};
    for(int j=0;j<3;j++)nmod_poly_init(vec+j,7);
    for(int o=0;o<4;o++)for(int shifted=0;shifted<2;shifted++) {
        slong index=-99,degree=-99;
        _nmod_poly_vec_pivot_profile(&index,&degree,vec,shifted?shift:NULL,3,orientations[o]);
        assert(degree==-1);
        assert(index==((orientations[o]==ROW_LOWER || orientations[o]==COL_UPPER)?-1:3));
    }
    for(int j=0;j<3;j++)nmod_poly_clear(vec+j);
    nmod_poly_mat_t m;nmod_poly_mat_init(m,2,2,7);
    nmod_poly_one(nmod_poly_mat_entry(m,0,0));
    slong zero_shift[2]={0,0},indices[2];
    assert(nmod_poly_mat_weak_popov_iter(m,zero_shift,NULL,indices,NULL,ROW_UPPER)==1);
    assert(indices[0]==0);
    assert(nmod_poly_is_zero(nmod_poly_mat_entry(m,1,0)));
    assert(nmod_poly_is_zero(nmod_poly_mat_entry(m,1,1)));
    nmod_poly_mat_clear(m);

}
static void conversions(flint_rand_t rng)
{
    const slong sizes[]={0,1,7,8,9,16,33},lengths[]={0,1,7,8,9,31,32,65,257};
    for(size_t a=0;a<7;a++)for(size_t b=0;b<7;b++) {
        slong r=sizes[a],c=sizes[b];
        nmod_poly_mat_t input,out,ref;
        nmod_poly_mat_init(input,r,c,65537);nmod_poly_mat_init(out,r,c,65537);
        nmod_poly_mat_init(ref,r,c,65537);
        nmod_mat_poly_t mp,expected;
        nmod_mat_poly_init(mp,r,c,65537);nmod_mat_poly_init(expected,r,c,65537);
        for(size_t t=0;t<9;t++) {
            nmod_poly_mat_randtest(input,rng,lengths[t]);
            if(r&&c)nmod_poly_zero(nmod_poly_mat_entry(input,r-1,c-1));
            for(int truncate=0;truncate<3;truncate++) {
                slong n=truncate==0?lengths[t]+5:truncate==1?lengths[t]/2:0;
                nmod_mat_poly_set_trunc_from_poly_mat(mp,input,n);forward(expected,input,n);
                assert(mp->length==expected->length);
                for(slong k=0;k<mp->length;k++)assert(nmod_mat_equal(mp->coeffs+k,expected->coeffs+k));
                /* Exercise logical row order, including old FLINT's pointer swaps. */
                slong *perm=flint_malloc(FLINT_MAX(r,1)*sizeof(slong));
                for(slong i=0;i<r;i++)perm[i]=r-1-i;
                nmod_mat_poly_permute_rows(mp,perm,NULL);flint_free(perm);
                nmod_poly_mat_set_trunc_from_mat_poly(out,mp,n);backward(ref,mp,n);
                assert(nmod_poly_mat_equal(out,ref));
            }
        }
        nmod_mat_poly_clear(mp);nmod_mat_poly_clear(expected);
        nmod_poly_mat_clear(input);nmod_poly_mat_clear(out);nmod_poly_mat_clear(ref);
    }
}
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static void bench(flint_rand_t rng)
{
    const slong sizes[]={12,32,64},lengths[]={32,512,1024};
    for(int z=0;z<3;z++) {
        slong n=sizes[z],len=lengths[z];int repeats=z==2?20:100;
        nmod_poly_mat_t a,b;nmod_mat_poly_t mp;
        nmod_poly_mat_init(a,n,n,65537);nmod_poly_mat_init(b,n,n,65537);
        nmod_mat_poly_init(mp,n,n,65537);nmod_poly_mat_randtest(a,rng,len);
        forward(mp,a,len);backward(b,mp,len);
        for(int direction=0;direction<2;direction++) {
            double seconds[2];
            for(int method=0;method<2;method++) {
                double start=now();
                for(int k=0;k<repeats;k++) {
                    if(direction==0){if(method) nmod_mat_poly_set_trunc_from_poly_mat(mp,a,len);else forward(mp,a,len);}
                    else {if(method)nmod_poly_mat_set_trunc_from_mat_poly(b,mp,len);else backward(b,mp,len);}
                }
                seconds[method]=(now()-start)/repeats;
            }
            printf("%ldx%ld len=%ld %s old=%.6g new=%.6g speedup=%.2f\n",n,n,len,direction?"backward":"forward",seconds[0],seconds[1],seconds[0]/seconds[1]);
        }
        nmod_poly_mat_clear(a);nmod_poly_mat_clear(b);nmod_mat_poly_clear(mp);
    }
}
int main(int argc,char **argv)
{
    flint_rand_t rng;flint_rand_init(rng);flint_rand_set_seed(rng,813,987);
    if(argc==2 && strcmp(argv[1],"--bench")==0)bench(rng);
    else {pivots();conversions(rng);puts("PASS: zero shifted pivots and 1323 conversion cases");}
    flint_rand_clear(rng);flint_cleanup();return 0;
}

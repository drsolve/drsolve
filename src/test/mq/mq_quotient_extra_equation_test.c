/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Pointwise experiment for the Bezoutian backend, NOT a complete solver.
 * Quotient both tensor factors by additional equations; compare the new
 * recurrence with explicit projection of the original square core.
 * Full RREF deliberately permits parameter-dependent pivots. */
#include "../../determinant/mq_quotient_impl.h"

static void extra_specialize(ulong *q,const ulong *input,slong count,ulong value)
{
    slong w=N+1;
    for(slong f=0;f<count;f++) {
        for(slong i=0;i<N-1;i++)for(slong j=i;j<N-1;j++)q[(f*N+i)*N+j]=input[(f*w+i)*w+j];
        for(slong i=0;i<N-1;i++)q[(f*N+i)*N+N-1]=nmod_add(input[(f*w+i)*w+N],nmod_mul(value,input[(f*w+i)*w+N-1],mod),mod);
        ulong c=input[(f*w+N)*w+N];
        c=nmod_add(c,nmod_mul(value,input[(f*w+N-1)*w+N],mod),mod);
        q[(f*N+N-1)*N+N-1]=nmod_add(c,nmod_mul(nmod_mul(value,value,mod),input[(f*w+N-1)*w+N-1],mod),mod);
    }
}
static void extra_normal_forms(grade *g,const ulong *q,slong count)
{
    for(slong d=0;d<=L;d++) {
        grade_support(g+d,d);
        if(d<2) {
            for(slong i=0;i<g[d].s;i++){g[d].free[i]=i;nmod_mat_entry(g[d].reduction,i,i)=1;}
            continue;
        }
        nmod_mat_t relations;nmod_mat_init(relations,count*g[d-2].s,g[d].s,prime);
        for(slong f=0;f<count;f++)for(slong a=0;a<g[d-2].s;a++)
            for(slong i=0;i<N;i++)for(slong j=i;j<N;j++) {
                uint64_t key=g[d-2].keys[a]+(UINT64_C(1)<<(key_bits*i))+(UINT64_C(1)<<(key_bits*j));
                nmod_mat_entry(relations,f*g[d-2].s+a,index_key(g+d,key))=q[(f*N+i)*N+j];
            }
        slong rank=nmod_mat_rref(relations),r=0,nf=0;
        g[d].h=g[d].s-rank;
        g[d].free=flint_realloc(g[d].free,FLINT_MAX(1,g[d].h)*sizeof(slong));
        slong *pivot=flint_malloc(FLINT_MAX(1,rank)*sizeof(slong));
        nmod_mat_clear(g[d].reduction);nmod_mat_init(g[d].reduction,g[d].s,g[d].h,prime);
        for(slong j=0;j<g[d].s;j++) {
            if(r<rank && nmod_mat_entry(relations,r,j))pivot[r++]=j;
            else g[d].free[nf++]=j;
        }
        assert(r==rank && nf==g[d].h);
        for(slong j=0;j<g[d].h;j++) {
            nmod_mat_entry(g[d].reduction,g[d].free[j],j)=1;
            for(slong i=0;i<rank;i++)nmod_mat_entry(g[d].reduction,pivot[i],j)=nmod_neg(nmod_mat_entry(relations,i,g[d].free[j]),mod);
        }
        flint_free(pivot);nmod_mat_clear(relations);
    }
    for(slong d=0;d<L;d++) {
        nmod_mat_init(g[d].variables,N,g[d].h*g[d+1].h,prime);
        for(slong v=0;v<N;v++)for(slong i=0;i<g[d].h;i++) {
            slong c=index_key(g+d+1,g[d].keys[g[d].free[i]]+(UINT64_C(1)<<(key_bits*v)));
            for(slong j=0;j<g[d+1].h;j++)
                nmod_mat_entry(g[d].variables,v,i*g[d+1].h+j)=nmod_mat_entry(g[d+1].reduction,c,j);
        }
    }
}
static void compare_projection(const nmod_mat_t old,const nmod_mat_t core,const grade *g,const grade *base)
{
    nmod_mat_t left,right,tmp,expected;
    nmod_mat_init(left,g[M].h,base[M].h,prime);nmod_mat_init(right,base[K].h,g[K].h,prime);
    for(slong j=0;j<base[M].h;j++)for(slong i=0;i<g[M].h;i++)
        nmod_mat_entry(left,i,j)=nmod_mat_entry(g[M].reduction,base[M].free[j],i);
    for(slong i=0;i<base[K].h;i++)for(slong j=0;j<g[K].h;j++)
        nmod_mat_entry(right,i,j)=nmod_mat_entry(g[K].reduction,base[K].free[i],j);
    nmod_mat_init(tmp,g[M].h,base[K].h,prime);nmod_mat_init(expected,g[M].h,g[K].h,prime);
    nmod_mat_mul(tmp,left,old);nmod_mat_mul(expected,tmp,right);assert(nmod_mat_equal(expected,core));
    nmod_mat_clear(left);nmod_mat_clear(right);nmod_mat_clear(tmp);nmod_mat_clear(expected);
}
int main(int argc,char **argv)
{
    assert(argc==2 || argc==3);FILE *in=fopen(argv[1],"r");assert(in);
    slong count;ulong root;
    assert(fscanf(in,"%ld %ld %lu %lu",&N,&count,&prime,&root)==4);assert(quotient_dimensions(N));
    assert(count>N);nmod_init(&mod,prime);threads=1;flint_set_num_threads(1);
    recurrence_kernel=packed_admitted()?5:1;assert(recurrence_kernel==5);
    slong w=N+1;ulong *input=flint_calloc(count*w*w,sizeof(ulong)),*q=flint_calloc(count*N*N,sizeof(ulong));
    for(slong f=0;f<count;f++)for(slong i=0;i<w;i++)for(slong j=i;j<w;j++)assert(fscanf(in,"%lu",input+(f*w+i)*w+j)==1);
    fclose(in);extra_specialize(q,input,count,0);assert(leading_init(q));frames_enabled=1;
    const ulong points[]={0,1,root,argc==3?strtoul(argv[2],NULL,10):0};
    for(slong sample=0;sample<(argc==3?4:3);sample++) {
        extra_specialize(q,input,count,points[sample]);
        grade *base=flint_calloc(L+1,sizeof(grade)),*g=flint_calloc(L+1,sizeof(grade));
        double t=now();assert(normal_forms(base,q));double old_normal=now()-t;
        nmod_mat_t old,core;nmod_mat_init(old,base[M].h,base[K].h,prime);
        t=now();quotient_det(old,base,NULL,q);double old_dp=now()-t;
        t=now();extra_normal_forms(g,q,count);double new_normal=now()-t;
        nmod_mat_init(core,g[M].h,g[K].h,prime);
        t=now();quotient_det(core,g,NULL,q);double new_dp=now()-t;
        compare_projection(old,core,g,base);
        slong old_rank=nmod_mat_rank(old),new_rank=nmod_mat_rank(core);
        printf("{\"n\":%ld,\"equations\":%ld,\"parameter\":%lu,\"known_root\":%s,\"old_size\":%ld,\"rows\":%ld,\"cols\":%ld,\"old_rank\":%ld,\"rank\":%ld,\"old_normal\":%.6f,\"old_recurrence\":%.6f,\"new_normal\":%.6f,\"new_recurrence\":%.6f,\"projection_equal\":true}\n",
            N,count,points[sample],sample==2?"true":"false",base[M].h,g[M].h,g[K].h,old_rank,new_rank,old_normal,old_dp,new_normal,new_dp);
        fflush(stdout);nmod_mat_clear(old);nmod_mat_clear(core);clear_grades(base);clear_grades(g);flint_free(base);flint_free(g);
    }
    /* The old backend interpolates core entries with N+2 samples. Test
     * that rule explicitly after the extra quotient; do not assume that
     * parameter-dependent RREF denominators preserve it. */
    if(N==8) {
        slong samples=N+2,hx=0,hy=0,mismatches=0,basis_changes=0;
        ulong *xs=flint_malloc(samples*sizeof(ulong)),*values=NULL;
        uint64_t *bx=NULL,*by=NULL;
        for(slong sample=0;sample<=samples;sample++) {
            extra_specialize(q,input,count,sample);
            grade *g=flint_calloc(L+1,sizeof(grade));extra_normal_forms(g,q,count);
            if(sample==0) {
                hx=g[M].h;hy=g[K].h;values=flint_malloc(hx*hy*samples*sizeof(ulong));
                bx=flint_malloc(hx*sizeof(uint64_t));by=flint_malloc(hy*sizeof(uint64_t));
                for(slong i=0;i<hx;i++)bx[i]=g[M].keys[g[M].free[i]];
                for(slong j=0;j<hy;j++)by[j]=g[K].keys[g[K].free[j]];
            }
            assert(hx==g[M].h && hy==g[K].h);
            for(slong i=0;i<hx;i++)basis_changes+=bx[i]!=g[M].keys[g[M].free[i]];
            for(slong j=0;j<hy;j++)basis_changes+=by[j]!=g[K].keys[g[K].free[j]];
            nmod_mat_t core;nmod_mat_init(core,hx,hy,prime);quotient_det(core,g,NULL,q);
            if(sample<samples) {
                xs[sample]=sample;
                for(slong i=0;i<hx;i++)for(slong j=0;j<hy;j++)values[(i*hy+j)*samples+sample]=nmod_mat_entry(core,i,j);
            } else {
                nmod_poly_t interpolated;nmod_poly_init(interpolated,prime);
                for(slong i=0;i<hx;i++)for(slong j=0;j<hy;j++) {
                    nmod_poly_interpolate_nmod_vec(interpolated,xs,values+(i*hy+j)*samples,samples);
                    mismatches+=nmod_poly_evaluate_nmod(interpolated,sample)!=nmod_mat_entry(core,i,j);
                }
                nmod_poly_clear(interpolated);
            }
            nmod_mat_clear(core);clear_grades(g);flint_free(g);
        }
        printf("{\"n\":%ld,\"equations\":%ld,\"old_interpolation_samples\":%ld,\"extra_point_mismatches\":%ld,\"entries\":%ld,\"basis_changes\":%ld}\n",
            N,count,samples,mismatches,hx*hy,basis_changes);
        flint_free(xs);flint_free(values);flint_free(bx);flint_free(by);
    }
    leading_clear();flint_free(q);flint_free(input);flint_cleanup_master();return 0;
}

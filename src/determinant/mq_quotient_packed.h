/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Canonical 32-bit residues, 64-bit accumulation. Admission covers every
 * update in the recurrence, including signs represented as p-c. */
#if defined(__AVX2__) || defined(__AVX512F__)
#include <immintrin.h>
#endif
static int packed_admitted(void)
{
    return prime<=UINT32_MAX && (ulong)(2*N*choose(N,L))<=UWORD_MAX/(prime-1)/(prime-1);
}
static inline void packed_addmul(ulong *restrict out,const uint32_t *restrict in,
                                uint32_t scalar,slong count)
{
    slong j=0;
#if defined(__AVX512F__) && FLINT_BITS == 64
    __m512i c=_mm512_set1_epi64(scalar);
    for(;j+8<=count;j+=8) {
        __m512i v=_mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i *)(in+j)));
        __m512i sum=_mm512_add_epi64(_mm512_loadu_si512(out+j),_mm512_mul_epu32(c,v));
        _mm512_storeu_si512(out+j,sum);
    }
#elif defined(__AVX2__) && FLINT_BITS == 64
    __m256i c=_mm256_set1_epi64x(scalar);
    for(;j+4<=count;j+=4) {
        __m256i v=_mm256_cvtepu32_epi64(_mm_loadu_si128((const __m128i *)(in+j)));
        __m256i sum=_mm256_add_epi64(_mm256_loadu_si256((const __m256i *)(out+j)),_mm256_mul_epu32(c,v));
        _mm256_storeu_si256((__m256i *)(out+j),sum);
    }
#endif
    for(;j<count;j++)out[j]+=(ulong)scalar*in[j];
}
#if defined(__AVX512F__) && FLINT_BITS == 64
static inline __m512i packed_madd(__m512i sum,__m512i a,__m512i b)
{
#if defined(__AVX512IFMA__)
    /* IFMA adds the low 52 product bits. This equals the full product
     * only when both residues are below 2^26. Sum overflow is covered
     * separately by packed_admitted(). */
    if(prime<=(UWORD(1)<<26))return _mm512_madd52lo_epu64(sum,a,b);
#endif
    return _mm512_add_epi64(sum,_mm512_mul_epu32(a,b));
}
#endif
/* Keep a 4-by-16 output tile in registers through the inner product.
 * Reuse each right-hand vector across four rows, without packing a second
 * copy of the state. The same admission bound covers these unreduced sums. */
static void packed_mat_addmul(ulong *restrict out,const packed_matrix *left,
                             const packed_matrix *right,int neg)
{
    slong rows=left->r,inner=left->c,cols=right->c,i=0;
#if defined(__AVX512F__) && FLINT_BITS == 64
    if(recurrence_kernel>=5 && inner>=32 && cols>=32) {
        for(;i+4<=rows;i+=4) {
            slong j=0;
            for(;j+16<=cols;j+=16) {
                __m512i s00=_mm512_loadu_si512(out+(i+0)*cols+j+0);
                __m512i s01=_mm512_loadu_si512(out+(i+0)*cols+j+8);
                __m512i s10=_mm512_loadu_si512(out+(i+1)*cols+j+0);
                __m512i s11=_mm512_loadu_si512(out+(i+1)*cols+j+8);
                __m512i s20=_mm512_loadu_si512(out+(i+2)*cols+j+0);
                __m512i s21=_mm512_loadu_si512(out+(i+2)*cols+j+8);
                __m512i s30=_mm512_loadu_si512(out+(i+3)*cols+j+0);
                __m512i s31=_mm512_loadu_si512(out+(i+3)*cols+j+8);
                for(slong r=0;r<inner;r++) {
                    uint32_t c0=left->data[(i+0)*inner+r];
                    uint32_t c1=left->data[(i+1)*inner+r];
                    uint32_t c2=left->data[(i+2)*inner+r];
                    uint32_t c3=left->data[(i+3)*inner+r];
                    if(!(c0|c1|c2|c3))continue;
                    if(neg) {
                        if(c0)c0=(uint32_t)prime-c0;
                        if(c1)c1=(uint32_t)prime-c1;
                        if(c2)c2=(uint32_t)prime-c2;
                        if(c3)c3=(uint32_t)prime-c3;
                    }
                    __m512i b0=_mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i *)(right->data+r*cols+j+0)));
                    __m512i b1=_mm512_cvtepu32_epi64(_mm256_loadu_si256((const __m256i *)(right->data+r*cols+j+8)));
                    __m512i a0=_mm512_set1_epi64(c0);
                    s00=packed_madd(s00,a0,b0);
                    s01=packed_madd(s01,a0,b1);
                    __m512i a1=_mm512_set1_epi64(c1);
                    s10=packed_madd(s10,a1,b0);
                    s11=packed_madd(s11,a1,b1);
                    __m512i a2=_mm512_set1_epi64(c2);
                    s20=packed_madd(s20,a2,b0);
                    s21=packed_madd(s21,a2,b1);
                    __m512i a3=_mm512_set1_epi64(c3);
                    s30=packed_madd(s30,a3,b0);
                    s31=packed_madd(s31,a3,b1);
                }
                _mm512_storeu_si512(out+(i+0)*cols+j+0,s00);
                _mm512_storeu_si512(out+(i+0)*cols+j+8,s01);
                _mm512_storeu_si512(out+(i+1)*cols+j+0,s10);
                _mm512_storeu_si512(out+(i+1)*cols+j+8,s11);
                _mm512_storeu_si512(out+(i+2)*cols+j+0,s20);
                _mm512_storeu_si512(out+(i+2)*cols+j+8,s21);
                _mm512_storeu_si512(out+(i+3)*cols+j+0,s30);
                _mm512_storeu_si512(out+(i+3)*cols+j+8,s31);
            }
            for(slong u=0;u<4;u++)for(slong r=0;r<inner;r++) {
                uint32_t c=left->data[(i+u)*inner+r];if(!c)continue;
                if(neg)c=(uint32_t)prime-c;
                packed_addmul(out+(i+u)*cols+j,right->data+r*cols+j,c,cols-j);
            }
        }
    }
#endif
    for(;i<rows;i++)for(slong r=0;r<inner;r++) {
        uint32_t c=left->data[i*inner+r];if(!c)continue;
        if(neg)c=(uint32_t)prime-c;
        packed_addmul(out+i*cols,right->data+r*cols,c,cols);
    }
}
/* Experimental parent batching: each equation-column map is reused by
 * up to four destination minors before moving to the next map. */
static void packed_layer_batched(packed_matrix *state,const grade *g,
                                 const linear_map *op,slong k,slong workers)
{
    slong count=choose(N,k),*masks=flint_malloc(count*sizeof(slong)),n=0;
    for(slong mask=1;mask<(1L<<N);mask++)if(__builtin_popcountl(mask)==k)masks[n++]=mask;
    assert(n==count);
    #pragma omp parallel for schedule(dynamic) num_threads(workers) if(workers>1)
    for(slong block=0;block<(count+3)/4;block++) {
        slong size=FLINT_MIN(4,count-block*4);
        for(slong a=FLINT_MAX(0,k-K);a<=FLINT_MIN(M,k);a++) {
            slong b=k-a,rows=g[a].h,cols=g[b].h;
            ulong *acc[4];
            for(slong u=0;u<size;u++) {
                slong mask=masks[block*4+u];
                state[mask*(M+1)+a]=(packed_matrix){rows,cols,flint_malloc(rows*cols*sizeof(uint32_t))};
                acc[u]=flint_calloc(rows*cols,sizeof(ulong));
            }
            for(slong j=0;j<N;j++)for(slong u=0;u<size;u++) {
                slong mask=masks[block*4+u];if(!(mask&(1L<<j)))continue;
                slong child=mask^(1L<<j);
                int neg=(k-1+__builtin_popcountl(mask&((1L<<j)-1)))&1;
                if(a)packed_mat_addmul(acc[u],&op[((k-1)*N+j)*L+a-1].px,
                    state+child*(M+1)+a-1,neg);
                if(b)packed_mat_addmul(acc[u],state+child*(M+1)+a,
                    &op[((k-1)*N+j)*L+b-1].py,neg);
            }
            for(slong u=0;u<size;u++) {
                packed_matrix *dest=state+masks[block*4+u]*(M+1)+a;
                for(slong i=0;i<rows*cols;i++)dest->data[i]=(uint32_t)nmod_set_ui(acc[u][i],mod);
                flint_free(acc[u]);
            }
        }
    }
    flint_free(masks);
}
static void quotient_det_packed(nmod_mat_t result,const grade *g,const linear_map *full,const ulong *q)
{
    double start=now(),map_seconds=0;slong masks=1L<<N,workers=quotient_recurrence_workers();
    packed_matrix *state=flint_calloc(masks*(M+1),sizeof(*state));
    state[0]=(packed_matrix){1,1,flint_malloc(sizeof(uint32_t))};state[0].data[0]=1;
    for(slong k=1;k<=N;k++) {
        double t=now();
        linear_map *layer=recurrence_kernel>=4?make_linear_maps_for_row(g,q,k-1):NULL;
        const linear_map *op=layer?layer:full;
        map_seconds+=now()-t;
        if(recurrence_kernel==6)packed_layer_batched(state,g,op,k,workers);
        else {
            #pragma omp parallel for schedule(dynamic) num_threads(workers) if(workers>1)
            for(slong mask=1;mask<masks;mask++)if(__builtin_popcountl(mask)==k) {
                for(slong a=FLINT_MAX(0,k-K);a<=FLINT_MIN(M,k);a++) {
                    slong b=k-a,rows=g[a].h,cols=g[b].h,pos=0;
                    packed_matrix *dest=state+mask*(M+1)+a;
                    *dest=(packed_matrix){rows,cols,flint_malloc(rows*cols*sizeof(uint32_t))};
                    ulong *acc=flint_calloc(rows*cols,sizeof(ulong));
                    for(slong j=0;j<N;j++)if(mask&(1L<<j)) {
                        slong child=mask^(1L<<j);int neg=(k-1+pos++)&1;
                        if(a) {
                            const packed_matrix *src=state+child*(M+1)+a-1;
                            const packed_matrix *map=&op[((k-1)*N+j)*L+a-1].px;
                            packed_mat_addmul(acc,map,src,neg);
                        }
                        if(b) {
                            const packed_matrix *src=state+child*(M+1)+a;
                            const packed_matrix *map=&op[((k-1)*N+j)*L+b-1].py;
                            packed_mat_addmul(acc,src,map,neg);
                        }
                    }
                    for(slong i=0;i<rows*cols;i++)dest->data[i]=(uint32_t)nmod_set_ui(acc[i],mod);
                    flint_free(acc);
                }
            }
        }
        clear_linear_maps(layer);
        for(slong mask=0;mask<masks;mask++)if(__builtin_popcountl(mask)==k-1)
            for(slong a=0;a<=M;a++){flint_free(state[mask*(M+1)+a].data);state[mask*(M+1)+a].data=NULL;}
    }
    packed_matrix *out=state+(masks-1)*(M+1)+M;
    for(slong i=0;i<out->r;i++)for(slong j=0;j<out->c;j++)nmod_mat_entry(result,i,j)=out->data[i*out->c+j];
    flint_free(out->data);flint_free(state);
    double elapsed=now()-start-map_seconds;
    #pragma omp atomic update
    dp_time+=elapsed;
}

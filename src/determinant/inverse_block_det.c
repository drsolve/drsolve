/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "inverse_block_det.h"
#include <math.h>
#include <stdint.h>
#include <time.h>
#include <flint/nmod_mat.h>
#include <flint/nmod.h>
#include <flint/ulong_extras.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef HAVE_PML
#include "pml_det.h"
#endif

static double seconds(void)
{
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9*t.tv_nsec;
}

slong dr_inverse_block_choose(slong n, slong bound, int threads)
{
    if (n < 1 || bound < 1) return 0;
    /* Balance dense sampling (~n^3 D/m / workers) with reconstruction
     * (~D m^2). The coefficient is calibrated on the bundled PML backend,
     * not a complexity theorem. Use 8-column steps only when worthwhile. */
    double target = 0.18*(double)n/cbrt(FLINT_MAX(1, threads));
    slong m = target >= 8 ? (slong)(target/8 + 0.5)*8 : (slong)(target + 0.5);
    return FLINT_MAX(1, FLINT_MIN(n, FLINT_MIN(bound, m)));
}

slong dr_nmod_poly_mat_det_degree_bound(const nmod_poly_mat_t a)
{
    slong n=a->r, maxdeg=0;
    if (n!=a->c || n<0 || n==WORD_MAX || (size_t)n+1>SIZE_MAX/sizeof(slong)) return -1;
    if (!n) return 0;
    for (slong i=0;i<n;i++) for (slong j=0;j<n;j++)
        maxdeg=FLINT_MAX(maxdeg,nmod_poly_degree(nmod_poly_mat_entry(a,i,j)));
    if (!maxdeg) return 0;
    if (maxdeg>WORD_MAX/4/n) return -1;
    /* Hungarian shortest augmenting paths for costs -degree. Allowing a
     * nonexistent entry to have weight zero can only increase the bound. */
    slong *u=flint_calloc(n+1,sizeof(slong)), *v=flint_calloc(n+1,sizeof(slong));
    slong *p=flint_calloc(n+1,sizeof(slong)), *way=flint_calloc(n+1,sizeof(slong));
    slong *dist=flint_malloc((n+1)*sizeof(slong));
    unsigned char *used=flint_malloc(n+1);
    for (slong i=1;i<=n;i++) {
        p[0]=i; slong j0=0;
        for (slong j=0;j<=n;j++) {dist[j]=WORD_MAX/2; used[j]=0;}
        do {
            used[j0]=1; slong i0=p[j0], delta=WORD_MAX/2, j1=0;
            for (slong j=1;j<=n;j++) if (!used[j]) {
                slong d=FLINT_MAX(0,nmod_poly_degree(nmod_poly_mat_entry(a,i0-1,j-1)));
                slong cur=-d-u[i0]-v[j];
                if (cur<dist[j]) {dist[j]=cur;way[j]=j0;}
                if (dist[j]<delta) {delta=dist[j];j1=j;}
            }
            for (slong j=0;j<=n;j++) {
                if (used[j]) {u[p[j]]+=delta;v[j]-=delta;}
                else dist[j]-=delta;
            }
            j0=j1;
        } while (p[j0]);
        do {slong j1=way[j0];p[j0]=p[j1];j0=j1;} while (j0);
    }
    slong bound=0;
    for (slong j=1;j<=n;j++)
        bound+=FLINT_MAX(0,nmod_poly_degree(nmod_poly_mat_entry(a,p[j]-1,j-1)));
    flint_free(u);flint_free(v);flint_free(p);flint_free(way);flint_free(dist);flint_free(used);
    return bound;
}

#ifdef HAVE_PML
static ulong determinant_at(nmod_mat_t value, const nmod_poly_mat_t a, ulong point)
{
    for (slong i=0;i<a->r;i++) for (slong j=0;j<a->c;j++)
        nmod_mat_entry(value,i,j)=nmod_poly_evaluate_nmod(nmod_poly_mat_entry(a,i,j),point);
    return nmod_mat_det(value);
}

static int verify(nmod_poly_t candidate, const nmod_poly_mat_t a, slong bound,
                  int workers, flint_rand_t state, dr_inverse_block_stats *stats)
{
    if (nmod_poly_is_zero(candidate) || nmod_poly_degree(candidate)>bound) return 0;
    ulong p=a->modulus; nmod_t mod; nmod_init(&mod,p);
    nmod_mat_t value; nmod_mat_init(value,a->r,a->c,p);
    int scaled=0;
    for (int attempt=0;attempt<32 && !scaled;attempt++) {
        ulong point=n_randint(state,p), actual=nmod_poly_evaluate_nmod(candidate,point);
        if (!actual) continue;
        ulong expected=determinant_at(value,a,point);
        if (!expected) continue;
        nmod_poly_scalar_mul_nmod(candidate,candidate,nmod_div(expected,actual,mod));
        scaled=1;
    }
    nmod_mat_clear(value);
    if (!scaled) return 0;
    /* 66 bits per candidate leaves room for the two-attempt union bound. */
    slong checks=FLINT_MAX(8,(slong)ceil(66.0/log2((double)p/(double)bound)));
    ulong *points=flint_malloc(checks*sizeof(ulong));
    for (slong k=0;k<checks;k++) points[k]=n_randint(state,p);
    int failed=0;
#ifdef _OPENMP
    #pragma omp parallel num_threads(workers) reduction(|:failed) if(workers>1)
#endif
    {
        int previous=flint_set_num_workers(0);
        nmod_mat_t val;nmod_mat_init(val,a->r,a->c,p);
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (slong k=0;k<checks;k++)
            failed |= determinant_at(val,a,points[k]) != nmod_poly_evaluate_nmod(candidate,points[k]);
        nmod_mat_clear(val);flint_reset_num_workers(previous);
    }
    flint_free(points);stats->checks+=checks;return !failed;
}

static int reconstruct(nmod_poly_t candidate, const nmod_poly_mat_t a,
    slong m, slong count, ulong root, ulong dilation, int workers,
    flint_rand_t state, dr_inverse_block_stats *stats)
{
    slong n=a->r;ulong prime=a->modulus; nmod_t mod;nmod_init(&mod,prime);
    if ((size_t)n>SIZE_MAX/(size_t)n || (size_t)n*n>SIZE_MAX/(size_t)count/sizeof(ulong) ||
        (size_t)count>SIZE_MAX/sizeof(nmod_mat_struct) ||
        (size_t)m>SIZE_MAX/(size_t)m/4/sizeof(nmod_poly_struct)) return 0;
    double start=seconds();
    ulong *values=flint_malloc((size_t)n*n*count*sizeof(ulong));
#ifdef _OPENMP
    #pragma omp parallel num_threads(workers) if(workers>1)
#endif
    {
        int previous=flint_set_num_workers(0);
        nmod_geometric_progression_t geo;
        _nmod_geometric_progression_init_function(geo,root,count,mod,1);
        nmod_poly_t scaled;nmod_poly_init(scaled,prime);
#ifdef _OPENMP
        #pragma omp for collapse(2) schedule(static)
#endif
        for (slong i=0;i<n;i++) for (slong j=0;j<n;j++) {
            const nmod_poly_struct *entry=nmod_poly_mat_entry(a,i,j);
            ulong *dest=values+((size_t)i*n+j)*count;
            nmod_poly_fit_length(scaled,entry->length);ulong power=1;
            for (slong k=0;k<entry->length;k++) {
                scaled->coeffs[k]=nmod_mul(entry->coeffs[k],power,mod);
                power=nmod_mul(power,dilation,mod);
            }
            scaled->length=entry->length;
            if (entry->length) _nmod_poly_evaluate_geometric_nmod_vec_fast_precomp(
                dest,scaled->coeffs,scaled->length,geo,count,mod);
            else for (slong k=0;k<count;k++) dest[k]=0;
        }
        nmod_poly_clear(scaled);nmod_geometric_progression_clear(geo);
        flint_reset_num_workers(previous);
    }
    stats->evaluation+=seconds()-start;start=seconds();
    nmod_mat_t u,v;nmod_mat_init(u,m,n,prime);nmod_mat_init(v,n,m,prime);
    for (slong i=0;i<n;i++) for (slong j=0;j<m;j++) {
        nmod_mat_entry(u,j,i)=n_randint(state,prime);
        nmod_mat_entry(v,i,j)=n_randint(state,prime);
    }
    nmod_mat_struct *blocks=flint_malloc(count*sizeof(*blocks));int failed=0;
#ifdef _OPENMP
    #pragma omp parallel num_threads(workers) reduction(|:failed) if(workers>1)
#endif
    {
        int previous=flint_set_num_workers(0);
        nmod_mat_t value,x,h;nmod_mat_init(value,n,n,prime);
        nmod_mat_init(x,n,m,prime);nmod_mat_init(h,m,m,prime);
#ifdef _OPENMP
        #pragma omp for schedule(static)
#endif
        for (slong k=0;k<count;k++) {
            nmod_mat_init(blocks+k,2*m,m,prime);
            for (slong i=0;i<n;i++) for (slong j=0;j<n;j++)
                nmod_mat_entry(value,i,j)=values[((size_t)i*n+j)*count+k];
            if (!nmod_mat_solve(x,value,v)) {failed=1;continue;}
            nmod_mat_mul(h,u,x);
            for (slong i=0;i<m;i++) {
                for (slong j=0;j<m;j++) nmod_mat_entry(blocks+k,i,j)=nmod_mat_entry(h,i,j);
                nmod_mat_entry(blocks+k,m+i,i)=1;
            }
        }
        nmod_mat_clear(value);nmod_mat_clear(x);nmod_mat_clear(h);
        flint_reset_num_workers(previous);
    }
    nmod_mat_clear(u);nmod_mat_clear(v);flint_free(values);
    stats->inverse+=seconds()-start;int ok=0;
    if (!failed) {
        slong *shift=flint_calloc(2*m,sizeof(slong)), *order=flint_malloc(2*m*sizeof(slong));
        nmod_poly_mat_t basis,den;
        nmod_poly_mat_init(basis,2*m,2*m,prime);nmod_poly_mat_init(den,m,m,prime);
        start=seconds();nmod_poly_mat_pmintbasis_geometric(basis,shift,NULL,blocks,root,count);
        for (slong i=0;i<2*m;i++) order[i]=i;
        for (slong i=1;i<2*m;i++) {
            slong row=order[i],j=i;
            while (j>0 && shift[order[j-1]]>shift[row]) {order[j]=order[j-1];j--;}
            order[j]=row;
        }
        for (slong i=0;i<m;i++) for (slong j=0;j<m;j++)
            nmod_poly_set(nmod_poly_mat_entry(den,i,j),nmod_poly_mat_entry(basis,order[i],j));
        stats->interpolation+=seconds()-start;start=seconds();
        ok=nmod_poly_mat_det_hnf(candidate,den);
        if (ok) {
            ulong inv=nmod_inv(dilation,mod),power=1;
            for (slong i=0;i<candidate->length;i++) {
                candidate->coeffs[i]=nmod_mul(candidate->coeffs[i],power,mod);
                power=nmod_mul(power,inv,mod);
            }
        }
        stats->determinant+=seconds()-start;
        nmod_poly_mat_clear(basis);nmod_poly_mat_clear(den);flint_free(shift);flint_free(order);
    }
    for (slong k=0;k<count;k++) nmod_mat_clear(blocks+k);
    flint_free(blocks);return ok;
}
#endif

int dr_nmod_poly_mat_det_inverse_block_try(nmod_poly_t out,
    const nmod_poly_mat_t a, slong block, int threads,
    flint_rand_t state, dr_inverse_block_stats *stats)
{
    dr_inverse_block_stats empty={0};if (!stats) stats=&empty;*stats=empty;
#ifndef HAVE_PML
    (void)out;(void)a;(void)block;(void)threads;(void)state;return 0;
#else
    double start=seconds();slong n=a->r;ulong p=a->modulus;
    if (n<1 || n!=a->c || out->mod.n!=p || block<0 || block>n || !n_is_prime(p)) return 0;
    int workers=1;
#ifdef _OPENMP
    if (!omp_in_parallel()) workers=threads>0 ? threads : omp_get_max_threads();
#else
    (void)threads;
#endif
    slong bound=dr_nmod_poly_mat_det_degree_bound(a);
    stats->degree_bound=bound;
    if (bound<1 || (ulong)bound>(p-1)/2) return 0;
    slong m=block ? block : dr_inverse_block_choose(n,bound,workers);
    if (m<1 || m>WORD_MAX/2) return 0;
    stats->block=m;slong half=bound/m+(bound%m!=0);
    if (half>(WORD_MAX-1)/2) return 0;
    slong count=2*half+1;
    if ((ulong)count>(p-1)/2) return 0;
    ulong root=n_primitive_root_prime(p);
    stats->setup=seconds()-start;
    nmod_poly_t candidate;nmod_poly_init(candidate,p);int ok=0;
    for (int attempt=0;attempt<2;attempt++) {
        stats->attempts++;stats->points=count;
        ulong dilation=1+n_randint(state,p-1);
        nmod_poly_zero(candidate);
        if (reconstruct(candidate,a,m,count,root,dilation,workers,state,stats)) {
            start=seconds();ok=verify(candidate,a,bound,workers,state,stats);
            stats->verification+=seconds()-start;
            if (ok) {nmod_poly_swap(out,candidate);break;}
        }
        if (count>(WORD_MAX-1)/2 || (ulong)(2*count+1)>(p-1)/2) break;
        count=2*count+1;
    }
    nmod_poly_clear(candidate);return ok;
#endif
}

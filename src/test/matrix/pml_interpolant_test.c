/* Independent module certificate: membership + ordered weak Popov form +
 * determinant degree equal to the codimension of the interpolation module. */
#include <assert.h>
#include <stdio.h>
#include <flint/ulong_extras.h>
#include <flint/nmod.h>
#include "../../../pml_det/pml_det.h"
#include "nmod_poly_mat_multiply.h"

int g_dixon_verbose_level = 0;
int g_dixon_debug_mode = 0;

static void check(const nmod_poly_mat_t b, const slong *s, const slong *out,
                  const ulong *points, const nmod_mat_struct *e, slong count)
{
    slong degree = 0, ranks = 0, rows = b->r;
    nmod_mat_t v, product;
    nmod_mat_init(v, rows, rows, b->modulus);
    nmod_mat_init(product, rows, count ? e[0].c : 0, b->modulus);
    for (slong k = 0; k < count; ++k) {
        nmod_poly_mat_evaluate_nmod(v, b, points[k]);
        nmod_mat_mul(product, v, e+k);
        assert(nmod_mat_is_zero(product));
        ranks += nmod_mat_rank(e+k);
    }
    for (slong i = 0; i < rows; ++i) {
        slong pivot = -1, best = WORD_MIN;
        for (slong j = 0; j < rows; ++j) {
            slong d = nmod_poly_degree(nmod_poly_mat_entry(b,i,j));
            if (d >= 0 && d+s[j] >= best) { best=d+s[j]; pivot=j; }
        }
        assert(pivot == i && best == out[i]);
        degree += best-s[i];
    }
    assert(degree == ranks);
    nmod_mat_clear(v); nmod_mat_clear(product);
}

static void trial(ulong p, slong count, slong rows, slong cols, int pattern,
                  flint_rand_t rng)
{
    nmod_t mod; nmod_init(&mod,p);
    ulong root = p == 2 ? 1 : n_primitive_root_prime(p);
    ulong *pts = flint_malloc(FLINT_MAX(count,1)*sizeof(ulong));
    ulong *got = flint_malloc(FLINT_MAX(count,1)*sizeof(ulong));
    nmod_mat_struct *e = flint_malloc(FLINT_MAX(count,1)*sizeof(*e));
    nmod_mat_struct *saved = flint_malloc(FLINT_MAX(count,1)*sizeof(*e));
    ulong point = 1, ratio = nmod_mul(root,root,mod);
    for(slong k=0;k<count;k++) {
        pts[k]=point; point=nmod_mul(point,ratio,mod);
        nmod_mat_init(e+k,rows,cols,p); nmod_mat_randtest(e+k,rng);
        if(pattern==1) nmod_mat_zero(e+k);
        if(pattern==2 && cols>1)
            for(slong i=0;i<rows;i++) nmod_mat_entry(e+k,i,1)=nmod_mat_entry(e+k,i,0);
        nmod_mat_init_set(saved+k,e+k);
    }
    slong *s=flint_malloc(rows*sizeof(slong)), *out=flint_malloc(rows*sizeof(slong));
    for(slong i=0;i<rows;i++) s[i]=pattern==3 ? 0 : 3*i-4;
    for(int method=0;method<4;method++) {
        nmod_poly_mat_t b; nmod_poly_mat_init(b,rows,rows,p);
        for(slong i=0;i<rows;i++) out[i]=s[i];
        if(method==0) nmod_poly_mat_mintbasis(b,out,pts,e,count);
        if(method==1) nmod_poly_mat_pmintbasis(b,out,pts,e,count);
        if(method==2) nmod_poly_mat_pmintbasis_geometric(b,out,got,e,root,count);
        if(method==3) nmod_poly_mat_pmintbasis_geometric_auto(b,out,got,e,count);
        if(method>=2) for(slong k=0;k<count;k++) assert(got[k]==pts[k]);
        check(b,s,out,pts,e,count);
        for(slong k=0;k<count;k++) assert(nmod_mat_equal(e+k,saved+k));
        nmod_poly_mat_clear(b);
    }
    for(slong k=0;k<count;k++) {nmod_mat_clear(e+k);nmod_mat_clear(saved+k);}
    flint_free(s);flint_free(out);flint_free(pts);flint_free(got);flint_free(e);flint_free(saved);
}
static void large_ntt(flint_rand_t rng)
{
    nmod_poly_mat_t a,b,c;
    nmod_poly_t expected;
    nmod_poly_mat_init(a,1,1,998244353);
    nmod_poly_mat_init(b,1,1,998244353);
    nmod_poly_mat_init(c,1,1,998244353);
    nmod_poly_init(expected,998244353);
    nmod_poly_randtest(nmod_poly_mat_entry(a,0,0),rng,33000);
    nmod_poly_randtest(nmod_poly_mat_entry(b,0,0),rng,33000);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(a,0,0),32999,1);
    nmod_poly_set_coeff_ui(nmod_poly_mat_entry(b,0,0),32999,1);
    nmod_poly_mul(expected,nmod_poly_mat_entry(a,0,0),nmod_poly_mat_entry(b,0,0));
    assert(drsolve_nmod_poly_mat_mul_window_ntt(c,a,b,0,65999));
    assert(nmod_poly_equal(expected,nmod_poly_mat_entry(c,0,0)));
    nmod_poly_clear(expected);
    nmod_poly_mat_clear(a);nmod_poly_mat_clear(b);nmod_poly_mat_clear(c);
}

int main(void)
{
    flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng,781,101);
    const ulong primes[]={65537,998244353,1000000007,
#if FLINT_BITS == 64
        UWORD(18446744069414584321),
#endif
    };
    const slong counts[]={0,1,31,32,33,65,97};
    slong cases=0;
    for(size_t p=0;p<sizeof(primes)/sizeof(primes[0]);p++) for(int d=0;d<7;d++)
        for(int shape=0;shape<3;shape++) for(int pattern=0;pattern<4;pattern++) {
            trial(primes[p],counts[d],shape==0?1:4,shape==0?1:shape==1?2:5,pattern,rng);cases++;
        }
    trial(17,8,4,2,0,rng);trial(2,1,4,2,0,rng);trial(2,0,4,2,0,rng);
    large_ntt(rng);
    flint_rand_clear(rng);flint_cleanup();
    printf("PASS: %ld interpolation cases, four backends each\n",cases+3);
    return 0;
}

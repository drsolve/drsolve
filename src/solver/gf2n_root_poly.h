/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Private packed polynomial multiplication for repeated Frobenius reduction.
 * n is a power of two; output and inputs must not overlap. Scratch: 4*n. */
static void gf2n_root_mul(gf2128_t *out, const gf2128_t *a,
                         const gf2128_t *b, slong n, gf2128_t *scratch,
                         const gf2n_native_ctx *ctx)
{
    memset(out, 0, 2*n*sizeof(*out));
    if (n <= 16) {
        for (slong i=0; i<n; ++i) {
            if (!(a[i].low | a[i].high)) continue;
            for (slong j=0; j<n; ++j) {
                gf2128_t v=gf2n_native_mul_fast(a[i],b[j],ctx);
                out[i+j].low ^= v.low;
                out[i+j].high ^= v.high;
            }
        }
        return;
    }
    slong h=n/2;
    gf2128_t *as=scratch, *bs=as+h, *mid=bs+h, *work=mid+n;
    for (slong i=0; i<h; ++i) {
        as[i]=(gf2128_t){a[i].low^a[h+i].low,a[i].high^a[h+i].high};
        bs[i]=(gf2128_t){b[i].low^b[h+i].low,b[i].high^b[h+i].high};
    }
    gf2n_root_mul(out,a,b,h,work,ctx);
    gf2n_root_mul(out+n,a+h,b+h,h,work,ctx);
    gf2n_root_mul(mid,as,bs,h,work,ctx);
    for (slong i=0; i<n; ++i) {
        mid[i].low ^= out[i].low^out[n+i].low;
        mid[i].high ^= out[i].high^out[n+i].high;
    }
    for (slong i=0; i<n; ++i) {
        out[h+i].low ^= mid[i].low;
        out[h+i].high ^= mid[i].high;
    }
}

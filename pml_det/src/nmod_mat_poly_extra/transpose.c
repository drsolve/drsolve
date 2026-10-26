/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include "transpose.h"

#define _PML_LOADPAD(buf, s, l, k, W)                                        \
    do {                                                                     \
        const slong _av = (l) - (k);                                         \
        if (_av >= (W))                                                      \
            memcpy((buf), (s) + (k), (W) * sizeof(ulong));                   \
        else                                                                 \
        {                                                                    \
            memset((buf), 0, (W) * sizeof(ulong));                           \
            if (_av > 0)                                                     \
                memcpy((buf), (s) + (k), (size_t) _av * sizeof(ulong));      \
        }                                                                    \
    } while (0)

/* --- 8x8 scalar kernel (portable fallback) --- */
#define _PML_KERNEL_SCALAR(e, k)                                             \
    do {                                                                     \
        ulong _b[8][8];                                                      \
        for (slong _t = 0; _t < 8; _t++)                                     \
            _PML_LOADPAD(_b[_t], src[(e) + _t], slen[(e) + _t], (k), 8);     \
        for (slong _u = 0; _u < 8; _u++)                                     \
        {                                                                    \
            nn_ptr _d = dst[(k) + _u] + (e);                                 \
            for (slong _t = 0; _t < 8; _t++)                                 \
                _d[_t] = _b[_t][_u];                                         \
        }                                                                    \
    } while (0)

#define _PML_MK_TRANSPOSE(NAME, W, KERNEL, DST_MAJOR)                           \
static void NAME(nn_ptr * dst, slong len,                                    \
                 nn_srcptr * src, const slong * slen, slong n)               \
{                                                                            \
    const slong nW = n - (n % (W));                                          \
    const slong lW = len - (len % (W));                                      \
                                                                             \
    if (DST_MAJOR)                                                         \
        for (slong k = 0; k < lW; k += (W))                                  \
            for (slong e = 0; e < nW; e += (W))                              \
                KERNEL(e, k);                                                \
    else                                                                     \
        for (slong e = 0; e < nW; e += (W))                                  \
            for (slong k = 0; k < lW; k += (W))                              \
                KERNEL(e, k);                                                \
                                                                             \
    for (slong k = lW; k < len; k++)                                         \
    {                                                                        \
        nn_ptr d = dst[k];                                                   \
        for (slong e = 0; e < n; e++)                                        \
            d[e] = (k < slen[e]) ? src[e][k] : UWORD(0);                     \
    }                                                                        \
    for (slong e = nW; e < n; e++)                                           \
    {                                                                        \
        nn_srcptr s = src[e];                                                \
        const slong l = FLINT_MIN(slen[e], lW);                              \
        for (slong k = 0; k < l; k++)                                        \
            dst[k][e] = s[k];                                                \
        for (slong k = l; k < lW; k++)                                       \
            dst[k][e] = UWORD(0);                                            \
    }                                                                        \
}

_PML_MK_TRANSPOSE(_pml_tr_sca_dm, 8, _PML_KERNEL_SCALAR, 1)
_PML_MK_TRANSPOSE(_pml_tr_sca_sm, 8, _PML_KERNEL_SCALAR, 0)

void drsolve_pml_transpose(nn_ptr *dst, slong ndst, nn_srcptr *src,
                          const slong *slen, slong nsrc, int dst_major)
{
    if (dst_major) _pml_tr_sca_dm(dst,ndst,src,slen,nsrc);
    else _pml_tr_sca_sm(dst,ndst,src,slen,nsrc);
}

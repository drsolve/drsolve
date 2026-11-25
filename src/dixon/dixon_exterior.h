/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Finite-field exterior ranks used by the mixed-degree predictor. */
#ifndef DIXON_EXTERIOR_H
#define DIXON_EXTERIOR_H
#include <flint/nmod_mat.h>

/* Advance a sorted k-subset of 0..n-1 in lexicographic order. */
static inline int next_subset(slong *a, slong k, slong n)
{
    slong j = k;
    while (j && a[j - 1] == n - k + j - 1)
        --j;
    if (!j)
        return 0;
    ++a[j - 1];
    for (; j < k; ++j)
        a[j] = a[j - 1] + 1;
    return 1;
}

slong dixon_exterior_rank(const slong *d, slong n, slong r, slong s, int top);
#endif

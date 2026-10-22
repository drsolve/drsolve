/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Same PMBasis/HNF algorithms, compact fmpz storage and whole-operation
 * arithmetic dispatch. Polynomial products accumulate over Z, then reduce
 * once per output coefficient instead of once per scalar multiply/add. */
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mod_poly.h>
#include <flint/fmpz_poly_mat.h>

static void big_mul_window(unified_mat_poly_t out, const unified_mat_poly_t a,
                           const unified_mat_poly_t b, slong start, slong len);
#define UM_NAME(x) big_##x
#define UM_E fmpz
#define UM_BYTE 0
#define UM_INIT(x) fmpz_init(x)
#define UM_CLEAR(x) fmpz_clear(x)
#define UM_ZERO(x) fmpz_zero(x)
#define UM_ONE(x) fmpz_one(x)
#define UM_SET(x, y) fmpz_set(x, y)
#define UM_ISZERO(x) fmpz_is_zero(x)
#define UM_ADD(x, y, z) fmpz_mod_add(x, y, z, ctx->ctx.big_ctx)
#define UM_SUB(x, y, z) fmpz_mod_sub(x, y, z, ctx->ctx.big_ctx)
#define UM_NEG(x, y) fmpz_mod_neg(x, y, ctx->ctx.big_ctx)
#define UM_MUL(x, y, z) fmpz_mod_mul(x, y, z, ctx->ctx.big_ctx)
#define UM_INV(x, y) fmpz_mod_inv(x, y, ctx->ctx.big_ctx)
#define UM_CUSTOM_MUL_WINDOW 1
#define UM_SUBMUL(x, y, z)                                                                         \
    do {                                                                                           \
        fmpz_submul(x, y, z);                                                                      \
        fmpz_mod_set_fmpz(x, x, ctx->ctx.big_ctx);                                                 \
    } while (0)
#include "unified_mat_poly_template.h"
#undef UM_NAME
#undef UM_E
#undef UM_BYTE
#undef UM_INIT
#undef UM_CLEAR
#undef UM_ZERO
#undef UM_ONE
#undef UM_SET
#undef UM_ISZERO
#undef UM_ADD
#undef UM_SUB
#undef UM_NEG
#undef UM_MUL
#undef UM_INV
#undef UM_CUSTOM_MUL_WINDOW
#undef UM_SUBMUL

static void big_to_integer(fmpz_poly_mat_t out, const unified_mat_poly_t a, slong len)
{
    const fmpz *data = a->data;
    for (slong i = 0; i < a->r; i++)
        for (slong j = 0; j < a->c; j++)
            for (slong d = 0; d < FLINT_MIN(a->length, len); d++)
                fmpz_poly_set_coeff_fmpz(fmpz_poly_mat_entry(out, i, j), d,
                                         data + (d * a->r + i) * a->c + j);
}

static void big_mul_window(unified_mat_poly_t out, const unified_mat_poly_t a,
                           const unified_mat_poly_t b, slong start, slong len)
{
    slong total = a->length && b->length ? a->length + b->length - 1 : 0;
    len = FLINT_MIN(len, FLINT_MAX(0, total - start));
    unified_mat_poly_t t;
    unified_mat_poly_init(t, a->r, b->c, a->ctx);
    big_resize(t, len);
    if (len && a->r == 1 && a->c == 1 && b->c == 1) {
        /* Compact 1x1 storage is already a FLINT coefficient vector. Borrow
         * read-only views; output remains separate even for aliased calls. */
        fmpz_mod_poly_struct x = {.coeffs = a->data, .alloc = a->alloc, .length = a->length};
        fmpz_mod_poly_struct y = {.coeffs = b->data, .alloc = b->alloc, .length = b->length};
        fmpz_mod_poly_t z;
        fmpz_mod_poly_init(z, a->ctx->ctx.big_ctx);
        fmpz_mod_poly_mullow(z, &x, &y, start + len, a->ctx->ctx.big_ctx);
        for (slong d = 0; d < len; d++)
            fmpz_mod_poly_get_coeff_fmpz((fmpz *)t->data + d, z, start + d, a->ctx->ctx.big_ctx);
        fmpz_mod_poly_clear(z, a->ctx->ctx.big_ctx);
    } else if (len && a->c) {
        fmpz_poly_mat_t x, y, z;
        fmpz_poly_mat_init(x, a->r, a->c);
        fmpz_poly_mat_init(y, b->r, b->c);
        fmpz_poly_mat_init(z, a->r, b->c);
        big_to_integer(x, a, start + len);
        big_to_integer(y, b, start + len);
        fmpz_poly_mat_mul(z, x, y);
        fmpz *dest = t->data;
        fmpz_t v;
        fmpz_init(v);
        for (slong i = 0; i < t->r; i++)
            for (slong j = 0; j < t->c; j++)
                for (slong d = 0; d < len; d++) {
                    fmpz_poly_get_coeff_fmpz(v, fmpz_poly_mat_entry(z, i, j), start + d);
                    fmpz_mod_set_fmpz(dest + (d * t->r + i) * t->c + j, v, a->ctx->ctx.big_ctx);
                }
        fmpz_clear(v);
        fmpz_poly_mat_clear(z);
        fmpz_poly_mat_clear(y);
        fmpz_poly_mat_clear(x);
    }
    big_normalise(t);
    big_swap(out, t);
    big_clear(t);
}

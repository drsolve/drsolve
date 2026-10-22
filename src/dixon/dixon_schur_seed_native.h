/* SPDX-License-Identifier: GPL-2.0-or-later */
#define SEED_NAME(x) seed_n8_##x
#define SEED_E uint8_t
#define SEED_MUL(a,b) gf2n_native_mul_word(a,b,field)
#define SEED_INV(a) gf2n_native_inv_word(a,field)
#define SEED_ADD(a,b) ((a)^(b))
#define SEED_ISZERO(a) (!(a))
#define SEED_FROM_WORD(a) (a)
#define SEED_FROM_FQ(a) gf2n_native_from_fq(a).low
#include "dixon_schur_seed_native_template.h"

#define SEED_NAME(x) seed_n16_##x
#define SEED_E uint16_t
#define SEED_MUL(a,b) gf2n_native_mul_word(a,b,field)
#define SEED_INV(a) gf2n_native_inv_word(a,field)
#define SEED_ADD(a,b) ((a)^(b))
#define SEED_ISZERO(a) (!(a))
#define SEED_FROM_WORD(a) (a)
#define SEED_FROM_FQ(a) gf2n_native_from_fq(a).low
#include "dixon_schur_seed_native_template.h"

#define SEED_NAME(x) seed_n32_##x
#define SEED_E uint32_t
#define SEED_MUL(a,b) gf2n_native_mul_word(a,b,field)
#define SEED_INV(a) gf2n_native_inv_word(a,field)
#define SEED_ADD(a,b) ((a)^(b))
#define SEED_ISZERO(a) (!(a))
#define SEED_FROM_WORD(a) (a)
#define SEED_FROM_FQ(a) gf2n_native_from_fq(a).low
#include "dixon_schur_seed_native_template.h"

#define SEED_NAME(x) seed_n64_##x
#define SEED_E uint64_t
#define SEED_MUL(a,b) gf2n_native_mul_word(a,b,field)
#define SEED_INV(a) gf2n_native_inv_word(a,field)
#define SEED_ADD(a,b) ((a)^(b))
#define SEED_ISZERO(a) (!(a))
#define SEED_FROM_WORD(a) (a)
#define SEED_FROM_FQ(a) gf2n_native_from_fq(a).low
#include "dixon_schur_seed_native_template.h"

#define SEED_NAME(x) seed_n128_##x
#define SEED_E gf2128_t
#define SEED_MUL(a,b) gf2n_native_mul(a,b,field)
#define SEED_INV(a) gf2n_native_inv(a,field)
#define SEED_ADD(a,b) ((gf2128_t){(a).low^(b).low,(a).high^(b).high})
#define SEED_ISZERO(a) (!((a).low|(a).high))
#define SEED_FROM_WORD(a) ((gf2128_t){a,0})
#define SEED_FROM_FQ(a) gf2n_native_from_fq(a)
#include "dixon_schur_seed_native_template.h"


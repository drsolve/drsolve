/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Packed relation elimination specialized for native binary fields. */
#define MQGF_TYPE uint8_t
#define MQGF_NAME quotient_insert_gf24
#define MQGF_BYTE
#define MQGF_MEMBER gf24
#define MQGF_ZERO(a) ((a) == 0)
#define MQGF_ADD(a,b) ((a) ^ (b))
#define MQGF_PREPARE(a) const uint8_t *scalar_row = gf24_get_scalar_row(a)
#define MQGF_MUL(a,b) scalar_row[b]
#include "quotient_gf2n_template.h"

#define MQGF_TYPE uint8_t
#define MQGF_NAME quotient_insert_gf28
#define MQGF_BYTE
#define MQGF_MEMBER gf28
#define MQGF_ZERO(a) ((a) == 0)
#define MQGF_ADD(a,b) ((a) ^ (b))
#define MQGF_PREPARE(a) const uint8_t *scalar_row = gf28_get_scalar_row(a)
#define MQGF_MUL(a,b) scalar_row[b]
#include "quotient_gf2n_template.h"

#define MQGF_TYPE uint16_t
#define MQGF_NAME quotient_insert_gf216
#define MQGF_MEMBER gf216
#define MQGF_ZERO(a) ((a) == 0)
#define MQGF_ADD(a,b) ((a) ^ (b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf216_mul(a,b)
#include "quotient_gf2n_template.h"

#define MQGF_TYPE gf232_t
#define MQGF_NAME quotient_insert_gf232
#define MQGF_MEMBER gf232
#define MQGF_ZERO(a) gf232_is_zero(&(a))
#define MQGF_ADD(a,b) gf232_add(&(a), &(b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf232_mul(&(a), &(b))
#include "quotient_gf2n_template.h"

#define MQGF_TYPE gf264_t
#define MQGF_NAME quotient_insert_gf264
#define MQGF_MEMBER gf264
#define MQGF_ZERO(a) gf264_is_zero(&(a))
#define MQGF_ADD(a,b) gf264_add(&(a), &(b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf264_mul(&(a), &(b))
#include "quotient_gf2n_template.h"

#define MQGF_TYPE gf2128_t
#define MQGF_NAME quotient_insert_gf2128
#define MQGF_MEMBER gf2128
#define MQGF_ZERO(a) gf2128_is_zero(&(a))
#define MQGF_ADD(a,b) gf2128_add(&(a), &(b))
#define MQGF_PREPARE(a) ((void)0)
#define MQGF_MUL(a,b) gf2128_mul(&(a), &(b))
#include "quotient_gf2n_template.h"


#define MQGF_PREFIX quotient_insert
#define MQGF_CONTEXT e->arithmetic
#define MQGF_TEMPLATE "solver/quotient_gf2n_template.h"
#define MQGF_USE_BYTE
#include "gf2n_native_instances.h"


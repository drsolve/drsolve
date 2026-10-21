/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Share the certified seed/completion algorithm while preserving native
 * prime arithmetic and full extension coefficients. */
#define SCHUR_SEED_NAME dixon_schur_seed_prime
#define SEED_MAT nmod_mat_t
#define SEED_INIT(a,r,c) nmod_mat_init(a,r,c,prime)
#define SEED_CALL(fn,...) nmod_mat_##fn(__VA_ARGS__)
#define SEED_ENTRY(a,i,j) nmod_mat_entry(a,i,j)
#define SEED_ZERO(a,i,j) (!SEED_ENTRY(a,i,j))
#define SEED_STORE(a,i,j,v) (SEED_ENTRY(a,i,j)=nmod_poly_get_coeff_ui(v,0))
#define SEED_ACCUMULATE(a,i,j,v) (SEED_ENTRY(a,i,j)=nmod_add(SEED_ENTRY(a,i,j),nmod_poly_get_coeff_ui(v,0),(a)->mod))
#include "dixon_schur_seed_template.h"
#undef SCHUR_SEED_NAME
#undef SEED_MAT
#undef SEED_INIT
#undef SEED_CALL
#undef SEED_ENTRY
#undef SEED_ZERO
#undef SEED_STORE
#undef SEED_ACCUMULATE

#define SCHUR_SEED_EXTENSION
#define SCHUR_SEED_NAME dixon_schur_seed_extension
#define SEED_MAT fq_nmod_mat_t
#define SEED_INIT(a,r,c) fq_nmod_mat_init(a,r,c,ctx)
#define SEED_CALL(fn,...) fq_nmod_mat_##fn(__VA_ARGS__,ctx)
#define SEED_ENTRY(a,i,j) fq_nmod_mat_entry(a,i,j)
#define SEED_ZERO(a,i,j) fq_nmod_is_zero(SEED_ENTRY(a,i,j),ctx)
#define SEED_STORE(a,i,j,v) fq_nmod_set(SEED_ENTRY(a,i,j),v,ctx)
#define SEED_ACCUMULATE(a,i,j,v) fq_nmod_add(SEED_ENTRY(a,i,j),SEED_ENTRY(a,i,j),v,ctx)
#include "dixon_schur_seed_template.h"
#undef SCHUR_SEED_EXTENSION
#undef SCHUR_SEED_NAME
#undef SEED_MAT
#undef SEED_INIT
#undef SEED_CALL
#undef SEED_ENTRY
#undef SEED_ZERO
#undef SEED_STORE
#undef SEED_ACCUMULATE

#include "dixon_schur_seed_byte.h"

static int dixon_schur_seed_select(unified_mpoly_struct ***matrix,
    slong nr, slong nc, const monom_t *rm, const monom_t *cm, slong m,
    const long *degrees, const fq_nmod_ctx_t ctx,
    slong **out_rows, slong **out_cols, slong *out_size)
{
    if (fq_nmod_ctx_prime(ctx)==2 &&
        (fq_nmod_ctx_degree(ctx)==4 || fq_nmod_ctx_degree(ctx)==8)) {
        seed_byte_field byte;
        field_ctx_init(&byte.field,ctx);
        if (byte.field.field_id!=FIELD_ID_GF24 && byte.field.field_id!=FIELD_ID_GF28) {
            field_ctx_clear(&byte.field);
            return dixon_schur_seed_extension(matrix,nr,nc,rm,cm,m,degrees,ctx,
                out_rows,out_cols,out_size);
        }
        int small=byte.field.field_id==FIELD_ID_GF24;
        byte.row=small ? gf24_get_scalar_row : gf28_get_scalar_row;
        byte.inverse=small ? gf24_inv : gf28_inv;
        quotient_byte_tables(byte.tables,small ? 16 : 256,byte.row);
        int ok=dixon_schur_seed_byte(matrix,nr,nc,rm,cm,m,degrees,ctx,
            out_rows,out_cols,out_size,&byte);
        field_ctx_clear(&byte.field);
        return ok;
    }
    return fq_nmod_ctx_degree(ctx) == 1
        ? dixon_schur_seed_prime(matrix,nr,nc,rm,cm,m,degrees,ctx,out_rows,out_cols,out_size)
        : dixon_schur_seed_extension(matrix,nr,nc,rm,cm,m,degrees,ctx,out_rows,out_cols,out_size);
}

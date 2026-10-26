#include "transpose.h"
#include "nmod_mat_poly.h"
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

void nmod_mat_poly_init_set_from_nmod_mat(nmod_mat_poly_t matp,
                                          const nmod_mat_t cmat)
{
    nmod_mat_poly_init2_preinv(matp, cmat->r, cmat->c, cmat->mod.n, cmat->mod.ninv, 1);
    // if cmat is zero: do nothing more
    // otherwise set length to 1 and copy data
    if (! nmod_mat_is_zero(cmat))
    {
        nmod_mat_init_set(matp->coeffs + 0, cmat);
        matp->length = 1;
    }
}

void nmod_mat_poly_set_from_nmod_mat(nmod_mat_poly_t matp,
                                     const nmod_mat_t cmat)
{
    if (nmod_mat_is_zero(cmat))
        nmod_mat_poly_zero(matp);
    else
    {
        nmod_mat_poly_fit_length(matp, 1);
        _nmod_mat_poly_set_length(matp, 1);
        nmod_mat_set(matp->coeffs + 0, cmat);
    }
}

void nmod_mat_poly_set_trunc_from_poly_mat(nmod_mat_poly_t matp,
                                           const nmod_poly_mat_t pmat,
                                           slong order)
{
    const slong len = nmod_poly_mat_max_length(pmat);
    if (order > len)
        order = len;

    // allocate and init coefficients
    nmod_mat_poly_fit_length(matp, order);
    _nmod_mat_poly_set_length(matp, order);

    /* Upstream's blocked transpose, using logical rows for compatibility
     * with the bundled matrix storage and old FLINT row permutations. */
    if (order < 8 || matp->c < 8 || (double)order*matp->r*matp->c < 512) {
        for (slong k=0;k<order;k++)
            for (slong i=0;i<matp->r;i++)
                for (slong j=0;j<matp->c;j++)
                    nmod_mat_poly_entry(matp,k,i,j)=nmod_poly_get_coeff_ui(nmod_poly_mat_entry(pmat,i,j),k);
    } else {
        nn_srcptr *src=flint_malloc(matp->c*sizeof(nn_srcptr));
        nn_ptr *dst=flint_malloc(order*sizeof(nn_ptr));
        slong *lengths=flint_malloc(matp->c*sizeof(slong));
        for (slong i=0;i<matp->r;i++) {
            for (slong j=0;j<matp->c;j++) {
                const nmod_poly_struct *p=nmod_poly_mat_entry(pmat,i,j);
                src[j]=p->coeffs; lengths[j]=FLINT_MIN(p->length,order);
            }
            for (slong k=0;k<order;k++) dst[k]=nmod_mat_entry_ptr(matp->coeffs+k,i,0);
            drsolve_pml_transpose(dst,order,src,lengths,matp->c,0);
        }
        flint_free(src);flint_free(dst);flint_free(lengths);
    }

    // normalize (useless if order==len)
    if (order < len)
        _nmod_mat_poly_normalise(matp);
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

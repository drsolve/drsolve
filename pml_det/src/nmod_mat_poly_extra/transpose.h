/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

/* Portable subset of PML 22bc6b81, adapted to drsolve's nmod_mat headers.
 * The caller provides logical row pointers, including permuted/window rows. */
#ifndef DRSOLVE_PML_TRANSPOSE_H
#define DRSOLVE_PML_TRANSPOSE_H
#include <string.h>
#include <flint/flint.h>
void drsolve_pml_transpose(nn_ptr *dst, slong ndst, nn_srcptr *src,
                          const slong *slen, slong nsrc, int dst_major);
#endif

# drsolve bundled PML determinant subset

This directory contains the **minimal PML subset** that drsolve vendors for its
prime-field polynomial-matrix determinant path.

- Upstream project: <https://github.com/vneiger/pml>
- Upstream version copied from this repository snapshot: `0.0.5-dev`
- License: see `COPYING` and `COPYING_FLINT`

This is **not** an untouched copy of upstream PML. It is a drsolve-local,
trimmed, modified subset containing only the determinant-related code that we
actually use.

## Retained source files

Only the files needed by drsolve's `fq_poly_mat_det.c` integration are kept.
This is based on the full 13-source determinant chain, plus one additional
support file (`nmod_mat_poly_mem.c`) that is still required for a clean drsolve
build and link.

The retained source files under `pml_det/src/` are:

- `nmod_poly_mat_extra/nmod_poly_mat_det.c`
- `nmod_poly_mat_extra/nmod_poly_mat_weak_popov.c`
- `nmod_poly_mat_extra/nmod_poly_mat_utils.c`
- `nmod_poly_mat_extra/degrees_pivots_leadingmatrix.c`
- `nmod_poly_mat_extra/kernel.c`
- `nmod_poly_mat_extra/approximant_basis.c`
- `nmod_poly_mat_extra/middle_product.c`
- `nmod_mat_poly_extra/nmod_mat_poly_mbasis.c`
- `nmod_mat_poly_extra/nmod_mat_poly_mem.c`
- `nmod_mat_poly_extra/nmod_mat_poly_arith.c`
- `nmod_mat_poly_extra/nmod_mat_poly_set_from.c`
- `nmod_mat_poly_extra/nmod_mat_poly_shift.c`
- `nmod_mat_extra/left_nullspace.c`
- `nmod_extra/nmod_find_root.c`
- `nmod_vec_extra/nmod_vec_dot_product.c`

## Local interface

`pml_det.h` is a drsolve-local minimal public header exposing only the
determinant entry points used by drsolve:

- `nmod_poly_mat_det_iter`
- `nmod_poly_mat_det_hnf`
- `nmod_poly_mat_det_hnf_knowing_degree`

The internal headers under `pml_det/src/` are trimmed copies of the upstream
PML headers required to build the retained source files. For convenience, the
subset also keeps a reduced umbrella header:

- `pml_det/src/nmod_poly_mat_extra.h`

which re-exports only the determinant-related internal headers still present in
this subset.

## drsolve-specific modifications

The internal full left-nullspace helper is named
`pml_nmod_mat_left_nullspace` to avoid colliding with FLINT when linking static
libraries. The compact helper used by the determinant routines is unchanged.

The most important local work is in:

- `pml_det/src/nmod_poly_mat_extra/nmod_poly_mat_det.c`

### Main algorithmic addition: `nmod_poly_mat_det_hnf`

We added and use:

- `nmod_poly_mat_det_hnf`
- `nmod_poly_mat_det_hnf_knowing_degree`

These are our own FLINT / `nmod_poly_mat` determinant routines for polynomial
matrices over prime fields.

They are based on the original PML `ntl-extras`
`determinant_generic_knowing_degree` idea, but were **completed, adapted, and
integrated on the FLINT side** for drsolve. So this is not just a direct copy
of an upstream FLINT routine: it is our finished FLINT implementation of that
determinant recursion.

### Why this path exists

The older iterative weak-Popov / Mulders-Storjohann determinant path
`nmod_poly_mat_det_iter` is still kept as a fallback, but drsolve also wanted:

- a recursive block determinant routine,
- a path with explicit success/failure return value,
- a degree-checking variant,
- a cleaner fit for drsolve's prime-field determinant dispatch.

### High-level structure of `nmod_poly_mat_det_hnf`

The routine recursively splits a square matrix
\( P = [P_L \mid P_R] \) into left and right column blocks:

1. build the left block `P_L`,
2. compute a kernel basis of `P_L^T`,
3. convert that kernel basis into row form and normalize it into ordered weak
   Popov form,
4. identify pivot rows and complementary rows,
5. construct:
   - a leading block from complementary rows of `P_L`,
   - a Schur-like block from the kernel basis times `P_R`,
   - a pivot submatrix extracted from the kernel basis,
6. recurse on these smaller square polynomial matrices,
7. combine the recursive determinants through
   \[
   \det(P) = \pm \frac{\det(\text{leading}) \det(\text{schur})}
   {\det(\text{kernel-pivots})},
   \]
   with the sign determined from the induced block row ordering.

### What was implemented here

Compared with the original `ntl-extras` determinant idea, this drsolve / FLINT
version adds the practical machinery needed to make the recursion usable:

- row / column extraction helpers,
- transpose / row-form conversion for the kernel basis,
- pivot-row bookkeeping and complementary-row reconstruction,
- block-order permutation-sign handling,
- exact FLINT polynomial division with remainder check,
- a boolean recursive success/failure contract,
- the wrapper `nmod_poly_mat_det_hnf_knowing_degree`.

### Semantics

- `nmod_poly_mat_det_hnf(det, mat)`
  - returns `1` on success, `0` on failure,
  - on success, `det` is the determinant polynomial.

- `nmod_poly_mat_det_hnf_knowing_degree(det, mat, degree)`
  - runs the same recursion,
  - additionally checks that the output determinant degree is exactly
    `degree`,
  - returns `1` only if both the recursion succeeds and the degree matches.

### How drsolve uses it

Inside `src/fq_poly_mat_det.c`, drsolve first tries this HNF-style recursive
path on the prime-field `nmod_poly_mat` image of the matrix. If it fails,
drsolve falls back to `nmod_poly_mat_det_iter`.

## Other subset / integration changes

Besides the determinant work above, the bundled copy was also modified so it
fits drsolve cleanly:

- only the determinant dependency chain is kept,
- unused modules such as I/O, Dixon, Hermite, random/test/timing code, and
  unrelated extras were removed,
- internal headers were trimmed to match the reduced subset,
- `pml_det.h` was added as a tiny drsolve-local interface instead of exposing
  a broader upstream PML API,
- drsolve public headers no longer expose PML headers,
- build scripts use an explicit minimal source list instead of relying on the
  whole upstream PML tree,
- the old `third_party/pml` layout was replaced by this smaller root-level
  `pml_det/` layout.

There is also one drsolve-side integration change outside the determinant chain:

- `src/fq_sparse_interpolation.c` now uses the restored
  `nmod_vec_dot_product_unbalanced` implementation from a minimal vendored
  `nmod_vec_extra` subset, since the temporary local fallback was noticeably
  worse for sparse interpolation performance.

## Windows notes

Windows-specific changes are intentionally secondary:

- machine-vector feature macros are forced off for this bundled subset,
- a few macro/config guards were adjusted to avoid FLINT / toolchain clashes,
- the repository layout was simplified so the bundled MinGW files now live in
  `mingw/` instead of under `third_party/`.

### Prime-field Step 4 products and degree ordering

The PMBasis middle product now discards input coefficients above the requested
window before multiplying. It first tries a generic radix-2 NTT matrix
product: each input entry is transformed once, pointwise matrix products use
FLINT, and the output entries are inverse-transformed. The same backend handles
PMBasis's final full product. For a middle window `[start, hi)`, the transform
length also satisfies `deg(A_trunc B_trunc) < length + start`, so cyclic wrap
cannot contaminate a requested coefficient. This remains valid for shifted
approximants whose degree exceeds `start`.

Automatic NTT dispatch uses the same arithmetic and admission checks for all
prime fields admitting the necessary roots of unity. There are no
modulus-specific multiplication or dispatch branches. It checks primality,
root availability, dimensions and workspace before writing the output. Transforms are limited to 65536 points and approximately
256 MiB of working arrays; rejected attempts fall back to FLINT. No extension
field backend is changed. The existing geometric middle product now checks its
input-degree assumptions before entering its reversal/interpolation code.

With OpenMP enabled, NTT input transforms, pointwise matrix products, and
output transforms share a worker team controlled by `--threads`. Each worker
has private transform and matrix scratch. Small products stay serial, and
calls inside an existing OpenMP team do not create another worker team. The
worker count is also limited by transform length and the 256 MiB workspace
budget, including per-worker scratch. Barriers finish all input reads before
output writes, preserving in-place and shared-window multiplication. This
parallelizes the products inside HNF/PMBasis; the dependent PMBasis recursion
remains sequential.

M-Basis kernel-row updates also distribute independent coefficient blocks
across workers. Dispatch requires at least 262144 multiply-add terms per
worker, keeps small updates serial, and avoids nested OpenMP teams. Workers
share the existing gathered-pivot scratch budget instead of multiplying it
by the thread count. Pivot-row shifts and the approximation-order loop remain
sequential because later iterations depend on the updated residuals.

For the LNZ determinant backend, matrices of 9--192 rows are stably sorted by
increasing row degree before the first split, with the permutation sign restored
and the input preserved. Already ordered matrices incur no copy. Sorting is a
performance heuristic, not a rank assumption: failures retry the original
order, and existing determinant verification/fallback remains in place. Larger
matrices retain their order by default because a 300-row benchmark incurred
more ZLS recursion after sorting.

Ablation switches (read by the bundled PML implementation):

- `DRSOLVE_PML_MUL=auto` (default): generic NTT dispatch plus bounded-input fallback.
- `DRSOLVE_PML_MUL=flint`: bounded-input FLINT products, without NTT.
- `DRSOLVE_PML_MUL=ntt`: attempt NTT for all eligible shapes and prime fields.
- `DRSOLVE_PML_MUL=legacy`: original full products without input truncation.
- `DRSOLVE_PML_DET_ROW_ORDER=auto` (default): degree ordering up to 192 rows.
- `DRSOLVE_PML_DET_ROW_ORDER=degree`: request degree ordering at any size.
- `DRSOLVE_PML_DET_ROW_ORDER=legacy`: retain original row order.

`make test-prime-step4` compares coefficient windows (including in-place and
shared-storage outputs) against classical multiplication at one, two, and four
threads, tests unsupported transforms and geometric degree bounds, and compares
complete determinants against FFLU over six prime fields, including singular
and nonprimitive inputs.
M-Basis tests compare serial and parallel bases and shifted degrees, including
rank-deficient inputs, and verify that the basis times the input vanishes
modulo the requested power of the variable.
The CMake/CTest equivalent is `drsolve_prime_step4` when bundled PML is enabled.
For end-to-end result equality and alternating single-thread timing runs:

```sh
python3 src/test/pml_prime_step4_bench.py --degrees 12 16 20 --repeats 3
```

Earlier timing measurements used an experimental F_65537-specific butterfly
reduction. That specialization and its modulus-specific dispatch have been
removed; use the benchmark script above to measure the generic backend.

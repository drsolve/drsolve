# Direct univariate matrix construction for verified projections

Step 2 now bypasses the intermediate `fq_mvpoly_t ***full_matrix` when:

* Step 1 has verified the projected target block;
* there is one parameter;
* the caller requests a univariate polynomial matrix;
* the collected row and column support sizes agree.

The branch is automatic. Unverified projections, multivariate coefficient
outputs and other callers retain the generic extraction path.

`src/dixon/dixon_projected_matrix.h` scans the existing term-to-row/column
indices to compute row parameter valuations, then column valuations after
row extraction. A further scan computes the post-extraction maximum degrees
for each axis. Sorting uses the same degree/index comparator as the old path;
`DRSOLVE_PREDICT_REORDER=0` still preserves the original order. The final
`fq_nmod_poly_mat` is allocated once, then filled in parallel by source row.
Each entry's coefficient capacity is fitted to its maximum remaining degree
before insertion. Each output row belongs to exactly one thread.

This preserves the extracted parameter power, selected row/column indices,
and the labels/permutation parity passed to the existing MQ Schur Step 4
preparation. Step 3 reuses Step 1's certificate, as before.

The direct path allocates no generic coefficient matrix pointer grid, per-cell
`fq_mvpoly_t` objects, per-term generic monomial arrays, or cloned per-term
parameter exponent arrays. It retains term row/column maps and a row-grouping
index, plus O(matrix order) auxiliary arrays per thread. The final matrix still
uses `fq_nmod_poly_mat`, including for prime fields; no `nmod_poly_mat` conversion
or Step 4 algorithm change is included. The source Dixon polynomial is not
modified or consumed and is still freed by its caller at the original point.
Consequently this removes an overlapping representation, not all possible
memory bottlenecks or an assurance that the n=11 computation fits in 64 GiB.

Validation commands:

```sh
make test-dixon-projected-matrix
make build/dixon_mq_step4_test
./build/dixon_mq_step4_test
```

The direct-path test uses a compile-time-only switch to run the old verified
projection path as an oracle. Twenty comparisons cover random projected MQ
systems with 3–5 eliminated variables over F257, nonzero row and column
parameter content, F7^2 coefficients, ordering enabled/disabled, and 1/4
threads. They compare every matrix entry, the extracted power, row/column
indices, and all prepared Schur profile fields. Synthetic cases also verify
that source parameter exponents are unchanged. The existing Step 4 integration
test exercises both resultant APIs and Schur versus the full determinant.
No n=11 or other large-memory run is required by these tests.

Result: all 20 direct/generic comparisons and the existing Step 4 integration
test passed. The CLI and library rebuilt successfully. Large-case peak memory
and timing were not measured on this machine.

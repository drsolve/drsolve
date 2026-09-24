# Direct matrix construction for verified projections

The first section records the initial fq implementation. The native, consuming
prime-field path described below supersedes it inside the resultant pipeline.

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

## Native prime-field construction and early source release

Both resultant APIs now use `nmod_poly_mat` for verified single-parameter
prime-field projections when the univariate determinant path is selected.
The implementation keeps the same content and degree-order scans, then packs
each retained source term into `(column, degree, ulong coefficient)` records,
grouped by source row. During this pass it frees the original exponent arrays
and field coefficient allocations. It releases the source term array and both
term-to-axis maps before allocating the native matrix. The consumed polynomial
has NULL terms and zero length/capacity and remains safe to clear normally.
Monomial labels have already been copied into independently owned storage.
The two caller index arrays formerly sized by the full term count are no longer
allocated on this path; selected indices are sized by matrix order inside the
extractor and freed after preparing Schur metadata.

This is a private consuming path. The public extraction API still treats its
input as const, and the direct fq path for extension fields remains unchanged.
An explicit alternative Step 4 backend still retains its original selection.

Each native scalar is one machine word rather than an `fq_nmod` coefficient
object. Step 4 accepts native matrices without making an fq-to-nmod full copy.
For Schur, row/column permutations are applied in place, preserving their sign
in the determinant factor. On success the original large matrix is released
before the smaller core determinant runs; on rejection the permuted native
matrix goes directly to the determinant backend, with the same sign correction.
The existing auto/HNF/iter selection, HNF normalization check and iterative
fallback policy are preserved on PML builds. Non-PML builds use FLINT's native
polynomial-matrix determinant.

For the reported 97,708,413-term input, 24-byte packed records occupy about
2.18 GiB. They coexist briefly with source terms while those terms are being
consumed, and remain until matrix filling finishes. The source's 20-coordinate
exponent vectors alone previously held about 14.56 GiB of payload, in addition
to coefficients, term structures and allocation overhead. The change removes
those live source objects before final matrix allocation and also eliminates
the two term-count-sized caller index buffers (about 1.46 GiB in that case).
These are storage calculations, not measured RSS: allocator retention, native
matrix capacity, and Step 4 workspaces still affect the process peak. No claim
is made that the entire n=11 computation is guaranteed below 64 GiB.

Validation completed on small inputs:

* The original 20 generic/direct comparisons still pass, including F7^2.
* 16 consuming native cases compare every scalar, parameter content, row/column
  order and all Schur profile fields against the original generic path. They
  also check that consumed source polynomials can be cleared safely.
* 96 native determinant comparisons cover auto/HNF/iter, Schur and forced
  Schur rejection, with and without degree reordering, at 1/4 threads.
* The existing Step 4 integration test passes for both resultant APIs.
* A five-equation F257 CLI run with seed 1790242205 reached the early-release
  marker, completed Step 2/3, compressed 32 to 10 with Schur, and finished Step 4.

The library/CLI rebuilt without warnings. No large-memory n=11 run or peak-RSS
measurement was performed on this machine. With `-v 2`, the new path reports:

```
Released source Dixon terms and term maps; allocating native prime-field matrix...
```

## Step 2 scan and fill experiment (packing superseded below)

Native metadata scans now use contiguous logical chunks with private row/column
minima/maxima and row counts. Row counts are collected in the row-valuation
pass, removing a separate full term scan. Prefix sums assign disjoint packed
ranges for each (chunk,row); chunk ordering preserves the original term order
within each row, including duplicate last-write semantics. The extra metadata
storage is O(threads * matrix order), not O(threads * term count). Native inputs
with at least 65,536 terms use the requested OpenMP worker count; smaller and
non-native inputs retain one metadata worker.

Packing is parallel, but source destruction is serial. An initial experiment
with parallel destruction showed no stable total improvement; allocator lock
contention is a possible explanation, not a proven allocator profile. The
retained version separates these phases. The existing compact record buffer
can now become fully populated before source destruction, so although no new
term-sized array was added, unchanged peak RSS is not claimed. The original
source and term maps are still released before final matrix allocation.

Each matrix cell is fitted and zeroed once at its known maximum degree. Scalar
coefficients are then written directly to its buffer, followed by length
normalization. This avoids repeated per-scalar length/gap checks while retaining
canonical polynomials. Each matrix row remains exclusively owned by one worker.

Verbose level 2 now prints independent wall times for direct metadata, packing,
source release, matrix initialization, matrix filling and buffer cleanup.
Existing support collection and degree-bound timings remain available.

Validation passed: 21 generic/direct comparisons, 17 consuming native matrix
comparisons, 102 native determinant comparisons, and the existing Step 4
integration suite. An added 80,000-term duplicate fixture crosses the parallel
threshold and checks deterministic ordering, content and cleanup with 4 threads.

Five-repeat interleaved CLI measurements compare pinned old/new executable and
library pairs, using eight equations over F257, seed 1790242205:

| Threads | Previous Step 2 median | Current Step 2 median |
|---:|---:|---:|
| 1 | 0.094 s | 0.093 s |
| 4 | 0.090 s | 0.086 s |

These small differences are within observed run-to-run variability; they do
not establish a significant overall speedup. Current source-release medians
were 0.038/0.039 s, compared with matrix-fill medians of 0.006/0.004 s at 1/4
threads. All old/new CLI result files matched after excluding timing lines.
No n=11 measurement was performed. Larger-case speed and peak-memory behavior
remain unverified; the new phase logs are intended to identify the actual
large-case bottleneck before further changes to the Step 1/2 representation.

The driver `src/test/mq_step2_bench.py` accepts `--reference-dir` containing a
previous `drsolve` and `libdrsolve.so`, pins library loading via
`LD_LIBRARY_PATH`, and alternates versions. Local raw results are ignored at
`src/test/data/mq_layout/step2_speed_final.json`; the earlier parallel-free
experiment is in `step2_speed.json`. Separate runs had different machine load
and must not be combined as a paired comparison.


## Restore streaming consumption after n=11 regression

User measurements at n=11, F257, 16 threads show total Step 2 wall time rising
from 64.127 s to 91.987 s. In the experimental version, packing alone took
61.416 s and separate source release took 13.283 s; metadata took 0.864 s and
matrix fill took 4.802 s. The small n=8 measurements did not predict this
regression. These timings do not establish whether page faults, memory
bandwidth or allocator/locality effects were responsible.

Native packing now restores the original serial streaming traversal: copy a
term into its row's compact buffer, then immediately free that term's source
allocations. This removes the second traversal for destruction and avoids
requiring every source allocation to remain live until the compact buffer is
fully populated. Per-worker packing offsets are removed. Parallel metadata
scans and direct coefficient-buffer filling remain enabled. No full-sized
intermediate coefficient matrix is reintroduced.

Verbose output reports `Step 2 direct pack/release` as one combined phase,
without per-term clock overhead. Matrix initialization, filling and cleanup
retain separate timers. The large n=11 case must be measured on the user's
machine; no restored n=11 speed or peak RSS is claimed here.

After restoring streaming consumption, the dynamic library and CLI rebuilt
successfully. All 21 generic/direct comparisons, 17 consuming native matrix
comparisons, 102 native determinant checks and the Step 4 integration suite
passed. No large-memory performance test was run.


## Optional bounded row staging

A subsequent user n=11 run measured 73.866 s total Step 2, of which streaming
pack/release took 58.422 s and matrix filling took 2.172 s. Restoring streaming
helped versus 91.987 s, but did not recover the original 64.127 s measurement.
Run-to-run effects and the remaining code differences have not been isolated.

`DRSOLVE_STEP2_PACK_BUFFER=1` enables experimental row staging. Up to 32 compact
records per row accumulate in a small buffer before a contiguous copy to that
row's final packed range. This aims to reduce scattered writes into the large
packed buffer. Source allocations are still released immediately per term;
row order, full-width indices/degrees/coefficients and duplicate last-write
semantics are preserved. Partial batches flush at the end. Staging records
plus occupancy counters are bounded by 8 MiB; capacity shrinks as matrix order
increases, falling back to ordinary streaming when fewer than two records per
row fit. There is no new restriction on supported n or field size.

The default remains ordinary serial streaming. Both paths report staging
capacity and combined pack/release time. This experiment cannot eliminate the
per-term allocator calls, and its large-case performance is not yet measured.
Compare the same seed, field and thread count on the large-memory machine:

```sh
DRSOLVE_STEP2_PACK_BUFFER=0 ./drsolve -r '[2]*11' 257 --seed 1790242205 --threads 16 --time -v 2
DRSOLVE_STEP2_PACK_BUFFER=1 ./drsolve -r '[2]*11' 257 --seed 1790242205 --threads 16 --time -v 2
```

Validation: both disabled and enabled staging passed all 21 generic/direct,
17 consuming-native and 102 determinant comparisons. The repeated-term
fixture now varies coefficients and uses 80,004 terms so that row tails also
exercise partial batch flushing. The Step 4 integration suite passed with
staging enabled. Library/CLI builds completed without warnings. No n=11 test
was run locally.


## Compact Step 1-to-Step 2 handoff (current default)

The optional row-staging experiment was slower on the user's large case. The
new default removes the representation conversion that it attempted to speed
up. It applies to the existing prime-field, one-parameter, verified MQ minor-DP
path when Step 4 uses a native univariate matrix. Explicit alternative Step 1
backends retain their previous paths. `DRSOLVE_MQ_COMPACT=0` restores the legacy
streaming conversion; `DRSOLVE_STEP2_PACK_BUFFER` has no effect on compact data.

The filtered native determinant is canonical ORD_LEX, so terms with equal
original-variable support are already contiguous. Output assembly creates:

* one array of `(slong column, slong degree, ulong coefficient)` records;
* row offsets, and one exponent vector per distinct row/column support.

Support IDs preserve first-occurrence order, including the column order used
by the old extractor. Full exponent vectors are hashed; no additional packed
axis/word-size admission condition is introduced. The records occupy 24 bytes
per term on the 64-bit build (about 2.18 GiB for 97,708,413 terms), plus small
support/offset arrays. This describes the handoff payload, not total peak RSS.
Assembly temporarily coexists with the native DP result, and matrix filling
coexists with the compact records; neither phase creates a generic polynomial
or a second term-sized index/packing array.

Rank verification evaluates the records directly into the same candidate
matrix, at the same points and in the same target ordering. If a deficient
candidate requires local Schur repair, only then materialize its conventional
polynomial and continue the old repair path. Rejected candidates release all
compact storage before full-Dixon fallback. Successful compact candidates keep
an empty, safely clearable conventional placeholder in the caller.

Step 2 scans contiguous records for valuations and degree ordering, then fills
native polynomial rows in parallel. It preserves parameter content and Schur
permutation signs. No support recollection, term-map allocation, per-term
object destruction or row repacking occurs. The compact term buffer is freed
before Step 3, and support labels after profile preparation. Debug degree
reports use the support labels and cached parameter degree. Very small verbose
polynomial displays may materialize at most 100 terms for printing.

Validation:

* `make test-mq-compact`: primes 2, 3, 257, eliminated-variable counts 2–5,
  1/4 threads, rectangular projections, exact materialization, native matrix
  coefficients and Schur profiles, reorder enabled/disabled, positive content,
  repaired/rejected candidates and an identically singular full fallback.
* Existing projected-matrix tests: 21 generic/direct, 17 consuming native,
  102 native determinant checks; existing Step 4 integration passed.
* Complete n=7/8 CLI results match with compact disabled/enabled (F257,
  seed 1790242205, 4 threads). Three interleaved runs per mode:

| Equations | Legacy Step 2 median | Compact Step 2 median | Legacy Step 1+2 median | Compact Step 1+2 median |
|---:|---:|---:|---:|---:|
| 7 | 0.034 s | 0.005 s | 0.083 s | 0.045 s |
| 8 | 0.151 s | 0.022 s | 0.574 s | 0.330 s |

These are small local measurements, with millisecond-resolution CLI timers;
they do not predict n=11 speed or memory usage. No n=11 run was attempted.
Reproduce with `python3 src/test/mq_compact_bench.py --output /tmp/compact.json`.
At verbosity 2, `MQ compact output` reports assembly time and payload size;
`Using compact Step 1 rows directly` confirms the new Step 2 path.

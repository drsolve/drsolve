# MQ Step 1 shared-support experiment (2026-09-24)

The n=7/8 measurements support using a common monomial layout for all minors
in a DP layer. **The guarded production implementation is now default on**;
the original ablations remain available in the benchmark. The executable accepts n=4..8 and generates prime-field
MQ inputs with one retained parameter. It is not a general-purpose determinant
API, and it does not establish a new asymptotic complexity bound.

## Statistics before optimization

The baseline computes exactly the existing sparse suffix-minor DP. After each
layer, an observer constructs a structural support envelope and counts actual
occupancy. Observer work is excluded from the reported layer time. The total
time and RSS of a profiling run include observer work and **must not be used
as uninstrumented performance measurements**.

Let E_i be the union of entry supports in row i. Starting from S_0={0}, form
S_k = (S_{k-1}+E_{n-k}) intersected with the applicable downward coefficient
ideal. This ignores column exclusivity and determinant cancellations, so it
is a safe envelope for every minor of those rows. It uses input entry supports,
not a sampled minor's nonzero pattern. The original layer-one copy is not
filtered, which the observer preserves. Subsequent supports are filtered as
in the baseline. The per-layer `projected` field in the raw sparse profile
means **the closure filter was active in that layer**; the final run record
identifies whether projection was requested for the whole determinant.

Dense fixture: F65537, random seed (132,941), n equations in n variables with
n-1 eliminated variables. The existing `random_mq` fixture allows zero random
coefficients and forces the x_0^2 coefficient to one in each polynomial to
ensure elimination degree two even over F2.

For n=8, the default projected path gave:

| Linear layer k | Minors | Structural monomials per minor | Occupancy | Sparse allocated MiB | Dense coefficients MiB |
|---:|---:|---:|---:|---:|---:|
| 3 | 56 | 336 | 100% | 0.86 | 0.14 |
| 4 | 70 | 1,581 | 99.99819% | 5.14 | 0.84 |
| 5 | 56 | 6,805 | 99.99843% | 20.87 | 2.91 |
| 6 | 28 | 27,167 | 99.99882% | 48.50 | 5.80 |
| 7 | 8 | 98,954 | 99.99848% | 67.95 | 6.04 |

All linear layers in both n=7 and n=8 exceeded 99.99% occupancy. Occupancy is
sum(nonzero coefficients)/(number of minors times structural support size).
The actual union fills the entire envelope in these linear layers. Root
occupancy is lower because of structural determinant cancellations.

Sparse allocation counts FLINT coefficient and exponent capacities of the
live layer. Dense coefficient sizes exclude shared keys, hash tables, maps,
the previous layer and temporary outputs. These are not process-RSS estimates.
The shared profiles separately report explicitly accounted workspace arrays;
that counter excludes the packed output polynomial and FLINT scratch space.

Raw profiles, including the quadratic-first control and exact comparisons:
[profiles.json](data/mq_layout/profiles.json).

## Implementation and ablations

The compile-time macro `DRSOLVE_MQ_LAYOUT_TEST` adds private ablation hooks to the
existing layered DP. Normal builds contain none of the experiment or observer,
but do contain the guarded production shared-index implementation.
The shared implementation is in [mq_layout_experiment.h](mq_layout_experiment.h).

The benchmark's `shared` argument selects:

* **0:** Existing sparse DP.
* **1:** Shared support, uint32 transition maps, and coefficient arrays in the
  linear layers. Convert the final child minors to sparse polynomials and use
  the existing root multiplication and tree reduction.
* **2:** Also use shared indices for the root. Each worker computes a cofactor
  contribution in a coefficient array; reduce coefficients, then pack/sort
  the final polynomial once.
* **3:** Mode 2 plus delayed modular reduction when its bound fits one limb.
  All other cases use ordinary `nmod_mul`/`nmod_add`.
* **4:** The production backend, including input validation and workspace
  preflight. In n=7/8 F65537 audits the driver asserts that this backend was
  actually entered, preventing unnoticed fallback from passing the audit.

For delayed reduction each fixed monomial shift is injective. A destination
coefficient receives at most k times the number of shifts products in a
linear layer, or the number of shifts in one root cofactor. With canonical
residues, checking `contributions*(p-1)^2 <= ULONG_MAX` proves no overflow;
the code checks this with division to avoid overflowing the bound itself.
Each array is reduced before the next layer, and root cofactor sums use
modular addition. Negative contributions use their canonical field residues.
There is no division by the layer index, so small characteristic is supported.

Keys and maps are constructed afresh for every determinant call. There is no
cross-run warm cache and no assumption that coefficients in every allocated
slot will be nonzero. Map construction currently runs serially; minors within
a layer run in parallel. Structural support is substantially smaller than a
dense rectangular exponent grid.

## Uninstrumented measurements

Three repetitions, sequential and interleaved with reversed variant order on
alternate repetitions; n=7/8, F65537, seed (132,941). Host: Intel i7-11700,
Linux under a Microsoft hypervisor, GCC -O3 -march=native -flto=auto, OpenMP,
FLINT 3.7.0-dev. No other experiment, build, or regression was run concurrently
by this benchmark driver. System load was not controlled.

Times include candidate construction, the determinant, conversion to the
ordinary output representation, and candidate verification. They exclude the
common cancellation-matrix construction and divided-difference row operations,
and exclude Step 2 onward. Therefore these are not whole-solver timings or
complete CLI Step 1 timings. All timed projected candidates verified.

| n | Threads | Sparse (0), s | Shared linear (1), s | Shared all (2), s | Shared all + delayed reduction (3), s | Mode 0 / mode 3 |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 1 | 0.18525 | 0.13221 | 0.06461 | 0.05761 | 3.22x |
| 7 | 4 | 0.07283 | 0.05798 | 0.05774 | 0.05085 | 1.43x |
| 8 | 1 | 1.53201 | 1.02563 | 0.63672 | 0.57978 | 2.64x |
| 8 | 4 | 0.51216 | 0.38188 | 0.41859 | 0.37751 | 1.36x |

The gains are primarily from shared representation, not delayed reduction.
Root map construction introduces a serial cost: at n=8 with four threads,
mode 2 is slower than mode 1, and mode 3 has nearly the same time as mode 1.
Mode 3's memory use remains substantially lower.

| n | Threads | Sparse process peak MiB | Mode 3 process peak MiB |
|---:|---:|---:|---:|
| 7 | 1 | 33.45 | 25.07 |
| 7 | 4 | 43.62 | 24.79 |
| 8 | 1 | 191.67 | 109.45 |
| 8 | 4 | 250.58 | 109.12 |

RSS is `getrusage(RUSAGE_SELF).ru_maxrss`, sampled after the timed computation
and before exact-reference auditing. It includes input polynomials, matrix
construction, libraries, and output representation; it is not backend-only
workspace. Timing runs disable both profiling and auditing. Canonical output
fingerprints, term counts and verification outcomes agree across all variants
and thread counts. Fingerprints are a repeated-run check, not the correctness
proof; the audit suite compares complete polynomials exactly.

Raw repetitions and medians: [timings.json](data/mq_layout/timings.json).
These local measurements do not predict n=9/10, different fields, or arbitrary
MQ sparsity. In particular small fields can make array occupancy much lower.

## Where to put the unique quadratic row

The divided-difference matrix has one quadratic row and n-1 affine-linear
rows. In the existing layout the quadratic row is **matrix row zero**, but
suffix DP processes from the bottom, so it participates **last**.

The control cyclically moves row zero to the bottom, preserving the relative
order of the linear rows, and corrects the determinant by (-1)^(n-1). The
original safe-filter-layer certificate is invalid after reordering, so that
skip optimization is disabled in the control. The downward ideal itself
remains valid. All reordered results match the original ordering exactly.

Putting the quadratic row first in DP raises intermediate degree bounds from
k to k+1 and introduces all original x variables immediately. For k>=2 its
combination with the last divided-difference row already activates all x/y/t
variables; the normal suffix of k linear rows has at most n+k active variables.
The eventual cheaper linear root does not compensate for larger middle layers.

For n=8, normal versus quadratic-first structural support sizes are 27,167
versus 98,515 in layer six, and 98,954 versus 211,428 in layer seven. The root
envelope is the same (381,766), as expected for this input and target ideal.

A **separate** three-repeat, interleaved row-order experiment gave median
seconds below. Do not compare its absolute times with the preceding experiment:
system load changed substantially between runs.

| n | Backend | Quadratic last in DP | Quadratic first in DP | First / last |
|---:|---|---:|---:|---:|
| 7 | Sparse | 0.26563 | 0.59078 | 2.22x |
| 7 | Shared (3) | 0.10403 | 0.12235 | 1.18x |
| 8 | Sparse | 2.32704 | 5.78933 | 2.49x |
| 8 | Shared (3) | 0.87585 | 1.55388 | 1.77x |

Thus retain the existing matrix-first / DP-last position for this backend.
This is not a theorem that every possible row order and input behaves alike.
Raw paired measurements: [row_order.json](data/mq_layout/row_order.json).

## Correctness and reproduction

The full audit comprises 42 exact polynomial comparisons: all shared variants
at n=7/8, one/four threads, full and projected outputs, two seeds, reordered
rows, and n=5 over F2, F7, F101 and F18446744073709551557. The last prime forces
ordinary modular arithmetic instead of delayed reduction. Each audit compares
against the original sparse algorithm with the original row order, including
all parameter coefficients. All passed. Eight additional profiled runs also
passed exact comparison. Raw audits: [audits.json](data/mq_layout/audits.json).

`make test-mq-layout test-minor-dp test-mq-filter` also passed, covering the
new smoke target and the existing determinant scheduling, projected-coefficient,
integration, fallback, packing, and filter-certificate regressions.

```
make build/mq_layout_bench
# n threads shared projected quadratic-first profile audit seed prime
./build/mq_layout_bench 8 4 3 1 0 0 1 132 65537
make test-mq-layout
python3 src/test/mq_layout_bench.py profile --output /tmp/mq-layout-profiles.json
python3 src/test/mq_layout_bench.py timing --output /tmp/mq-layout-timings.json
python3 src/test/mq_layout_bench.py row-order --output /tmp/mq-layout-row-order.json
python3 src/test/mq_layout_bench.py audit --output /tmp/mq-layout-audits.json
```

## Production default

The retained implementation is
[`mq_shared_layout.h`](../determinant/mq_shared_layout.h), privately included by
the determinant backend. The library and CLI default to sharing, with explicit
`--no-mq-step1-shared` / `--mq-step1-shared` opt-out/re-enable flags. Existing
explicit Step 1 backends and the minor-cache entry limit retain precedence.
The coefficient-filter switch remains independent.

Preflight checks lex packing on a 64-bit build, at most three exponent words,
2n-1 variables, and total degree at most two in the first row and one in every
other row. Each layer's active-variable/degree simplex bounds the support
without assuming projection success or accidental coefficient zeros. The
workspace estimate includes both coefficient layers, transition maps, support
keys/hash capacities, key reallocation overlap, factors, and an output packing
allowance. A bound above 256 MiB rejects the shared path before large arrays
are allocated, and the existing sparse DP computes the result instead.
This is a conservative admission limit, not a whole-process RSS guarantee;
inputs, allocator overhead, and FLINT sorting scratch are outside that estimate.

All repacking operations execute outside assertions, including in release
builds. Empty support returns zero, and rejected calls preserve caller output.
The filter kernel regression checks n=8 admission, n=9 budget rejection,
degree-three rejection, output preservation, and zero input.

The production audit passed 24 exact comparisons against sharing-disabled
sparse DP: n=7/8, full/projected, one/four threads, two seeds; and n=5 over
F2, F7, F101, and a near-full-word prime. See
[production_audits.json](data/mq_layout/production_audits.json). These timings
are correctness-run timings, not a new paired performance measurement.

```
make test-mq-layout test-mq-filter test-minor-dp test-mq-shared-cli
python3 src/test/mq_layout_bench.py production --output /tmp/mq-shared-audits.json
./drsolve --no-mq-step1-shared --threads 4 -r --seed 12345 -n 8 '[2]*8' 65537
```

The next optimization targets are faster/parallel index construction and
reducing the root's per-cofactor coefficient buffers. Both need fresh paired
measurements; the historical tables above describe the experimental variants.

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
`data/mq_layout/profiles.json`.

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
* **4:** The original production native-key/serial backend, including input
  validation and workspace preflight, retained as the baseline for later
  ablations. In n=7/8 F65537 audits, modes >=4 assert that the shared backend
  was actually entered, preventing unnoticed fallback from passing the audit.

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

Raw repetitions and medians: `data/mq_layout/timings.json`.
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
Raw paired measurements: `data/mq_layout/row_order.json`.

## Correctness and reproduction

The full audit comprises 42 exact polynomial comparisons: all shared variants
at n=7/8, one/four threads, full and projected outputs, two seeds, reordered
rows, and n=5 over F2, F7, F101 and F18446744073709551557. The last prime forces
ordinary modular arithmetic instead of delayed reduction. Each audit compares
against the original sparse algorithm with the original row order, including
all parameter coefficients. All passed. Eight additional profiled runs also
passed exact comparison. Raw audits: `data/mq_layout/audits.json`.

`make test-mq-layout test-minor-dp test-mq-filter` also passed, covering the
new smoke target and the existing determinant scheduling, projected-coefficient,
integration, fallback, packing, and filter-certificate regressions.

```
make build/mq_layout_bench
# n threads shared projected quadratic-first profile audit seed prime
./build/mq_layout_bench 8 4 12 1 0 0 1 132 65537
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
`data/mq_layout/production_audits.json`. These timings
are correctness-run timings, not a new paired performance measurement.

```
make test-mq-layout test-mq-filter test-minor-dp test-mq-shared-cli
python3 src/test/mq_layout_bench.py production --output /tmp/mq-shared-audits.json
./drsolve --no-mq-step1-shared --threads 4 -r --seed 12345 -n 8 '[2]*8' 65537
```

The following iteration optimizes index construction. Reducing the root's
per-cofactor coefficient buffers remains a separate candidate.

## Compact keys and parallel map construction (historical)

The production policy measured in this section kept n<=7 on its previous
shared-array path (see the later default-policy update below).
For n>=8 it used a four-bit internal exponent representation when all fields
fit one limb and the admitted total degree is below 16. For n=8, the 15
variables fit in 60 bits, replacing two native FLINT exponent words with one.
Addition is safe because the degree bound rules out carries between fields.
Input and output still use ordinary FLINT packing; nibble spreading converts
internal fields back to eight-bit fields, including the x/y filter keys.
Larger layouts retain native exponent packing.

Maps with at least 1,000,000 shift/child pairs can be built in parallel. The
number of shards is the largest power of two not exceeding the thread count,
capped at eight. Hash ownership partitions the resulting monomials into
disjoint sets. Each shard scans the packed sums and inserts only its own keys
into a private hash table; it also owns the corresponding map entries. There
are no hash-table locks or cross-shard duplicate keys. Local hash buckets use
the hash bits above the ownership bits to avoid clustering. After construction,
concatenate the disjoint supports, remap local IDs, and discard shard tables.
Arithmetic uses the ordinary shared coefficient-array recurrence.

The conservative 256 MiB admission limit is retained with 64 KiB reserved for
small shard tables/rounding. Map construction finishes before current-layer
coefficient arrays are allocated. That unused coefficient allowance, together
with the existing support bound, covers the shard keys/hash tables and their
concatenation buffer; compact packing only decreases the native-width bound.
The limit remains a workspace admission estimate, not a whole-process RSS cap.

Three experiments were kept separately to avoid selecting results from
different runs as if they were paired:

* `data/mq_layout/maps_timings.json`: three-repeat ablation
  of native/compact keys and serial/sharded maps, using a 32,768-pair threshold.
* `data/mq_layout/ordered_timings.json`: adds a common
  radix-sorted monomial order and remaps transitions. This makes each fixed
  shift's writes monotone and reduces arithmetic time, but extra sorting/map
  remapping does not consistently improve total time. It remains test-only.
* `data/mq_layout/maps_final.json`: five-repeat check of the
  compact/sharded variant with the original low threshold. n=8 improved, but
  n=7 four-thread median time increased from 0.0742 to 0.0800 seconds. This
  motivated retaining the original n<=7 path and raising the shard threshold.

The recorded **five-repeat interleaved** comparison used benchmark mode 4 for
the previous production shared-index implementation and mode 12 for the policy
at measurement time. Field F65537, seed (132,941), same host and timing boundaries as above:

| n | Threads | Previous shared, s | Retained policy, s | Previous / retained | Previous RSS MiB | Retained RSS MiB |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 1 | 0.09339 | 0.09817 | Same algorithm | 25.21 | 25.28 |
| 7 | 4 | 0.08070 | 0.09554 | Same algorithm | 24.89 | 24.91 |
| 8 | 1 | 0.88817 | 0.67590 | 1.31x | 109.35 | 109.42 |
| 8 | 4 | 0.81259 | 0.50147 | 1.62x | 109.06 | 117.32 |

In that recorded comparison, n=7 took the same native-key serial construction
in both modes;
its timing differences illustrate the system-load noise, not an algorithmic
benefit or a reason to enable the new plan there. For n=8 with four threads,
median map-building time fell from 0.5023 to 0.1758 seconds. Arithmetic plus
packing rose from 0.1920 to 0.2135 seconds, while total time improved. The
extra roughly 8 MiB process peak is a tradeoff of the parallel shard buffers
and allocator lifetimes. Stage medians need not sum to the total median.

Raw final comparison: `data/mq_layout/maps_policy.json`.
These are local results for the tested input, not portable speed guarantees.

The variant audit passed 37 exact comparisons, including the rejected sorting
experiment: `data/mq_layout/maps_audits.json`. The retained
policy passed another 29 exact comparisons against sparse DP, including full
and projected n=7/8 outputs, two seeds, n=8 over F2/F7/a near-full-word prime,
and thread counts 3 and 16:
`data/mq_layout/maps_policy_audits.json`.
The kernel regression also compares every transition against serial construction
for a three-word layout, with/without filtering and 3/4/16 requested threads.

Benchmark modes 5/6/7 select compact, sharded, and combined construction;
9/11 add experimental ordering; **12 follows the actual production policy**.
Modes 4..11 are benchmark ablations, not public solver switches. Reproduction:

```
make build/mq_layout_bench
python3 src/test/mq_layout_bench.py maps-policy --repeats 5 --output /tmp/maps-policy.json
python3 src/test/mq_layout_bench.py maps-audit --output /tmp/maps-audit.json
python3 src/test/mq_layout_bench.py production --output /tmp/maps-policy-audit.json
make test-mq-layout test-mq-filter test-minor-dp test-mq-shared-cli
```

## Size-independent default policy (superseded)

At user request, production now enables compact keys and sharded map
construction for every admitted n, removing the n>=8 performance gate. Compact
keys still require all four-bit fields to fit in one word and total degree
below 16; otherwise native packing is used. Parallel construction still
requires at least 1,000,000 transition pairs and enough threads. The matrix
shape, native packing and 256 MiB workspace admission checks remain in force;
ineligible inputs fall back to sparse DP.

`--no-mq-step1-shared` disables the shared backend and its optimizations;
`--mq-step1-shared` re-enables it. Benchmark mode 12 follows this updated
production policy. No tests or benchmarks were rerun for this policy change,
as requested; the measurements and audits above refer to the earlier policy
and do not establish performance for larger n or the new n<=7 default.

## Current policy: native keys and shared indices

The small-layout four-bit encoding, nibble-spreading filter and output
conversion have been removed, including their benchmark implementations.
Shared indices and hash-sharded map construction remain enabled by default
for all admitted sizes, using native FLINT packing, including multiword keys.
The 1,000,000-pair parallel threshold and existing memory/packing admission
limits remain. This removes repeated per-minor indexing, not the exponential
subset count; it does not establish unrestricted-size support or a new
asymptotic bound for the whole determinant algorithm.

Current benchmark modes are 4 (serial), 6 (sharded), 8 (serial with experimental
ordering), 10 (sharded with experimental ordering), and 12 (production policy).
The retired compact flag has no effect. Historical timings above describe
removed or superseded variants and must not be used as current speed claims.
No tests or timings were rerun for this change, following the earlier request.

`src/test/data/mq_layout/` holds generated local measurements, not required test
fixtures, and is now ignored and removed from the Git index; local files are preserved.
Historical data paths above refer to local experiment artifacts, which are not
distributed with the repository. New runs should write locally or to `/tmp`
using the commands above.

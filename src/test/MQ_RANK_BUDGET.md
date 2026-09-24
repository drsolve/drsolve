# MQ direct-index feasibility and map-time budget

Direction 3 replaces support hashing with direct indices in a known support
envelope. The feasibility section below records the initial budget analysis. The later
implementation section records the integrated configurable backend and its tests.

The n convention and target candidate are those of `MQ_SUPPORT_COUNT.md`:
n is the determinant order, and targets are the initial production candidate
rectangle before numerical verification/repair. Feasibility uses the dense
row-support fixture, while timings use the existing random F65537 fixture,
seed (132,941), projected output, quadratic row last.

## Rankable envelope

At layer k, restrict each axis's existing target downward closure to variables
appearing in the first k processed rows. Call the resulting sets X_k and Y_k.
Use every (a,b,t) with a in X_k, b in Y_k, t>=0 and
|a|+|b|+t<=d_k, where d_k=k for linear layers and d_n=n+1.
This contains the actual row-union support: inactive variables cannot occur,
row-degree sums bound total degree, and the existing ideal bounds each axis.
Extra unreachable monomials remain zero under the ordinary recurrence.

A direct contiguous index exists without packed-word hash keys. Order both
axes lexicographically. For each residual degree B, precompute prefix sums
over Y_k of max(0,B-|b|+1). Precompute an x-block offset by summing those block
sizes for preceding x exponents. The index is then

```
x_block_offset[a] + y_prefix[d_k - degree(a)][b] + t
```

Axis lookup can use exponent-vector comparisons or a trie; shifts can also
use precomputed per-axis transitions. None of these requires all exponents
to fit one machine word. There are still storage and indexing limits, and
explicitly enumerated axis sets can grow exponentially. This is not a claim
of polynomial total complexity or already unrestricted-size implementation.

| n | Layer | Existing support | Envelope | Extra slots |
|---:|---:|---:|---:|---:|
| 7 | 1 | 9 | 9 | 0% |
| 7 | 2 | 53 | 54 | 1.89% |
| 7 | 3 | 260 | 273 | 5.00% |
| 7 | 4 | 1,156 | 1,259 | 8.91% |
| 7 | 5 | 4,725 | 5,211 | 10.29% |
| 7 | 6 | 17,747 | 17,747 | 0% |
| 7 | 7 | 71,645 | 71,645 | 0% |
| 8 | 1 | 10 | 10 | 0% |
| 8 | 2 | 64 | 65 | 1.56% |
| 8 | 3 | 336 | 350 | 4.17% |
| 8 | 4 | 1,581 | 1,700 | 7.53% |
| 8 | 5 | 6,805 | 7,537 | 10.76% |
| 8 | 6 | 27,167 | 30,112 | 10.84% |
| 8 | 7 | 98,954 | 98,954 | 0% |
| 8 | 8 | 381,766 | 381,766 | 0% |

Weighting slots by the number of coefficient arrays in each layer (binomial
for inner layers, n for the current cofactor-parallel root) gives cumulative
coefficient-slot increases of 1.84% and 2.59% for n=7 and n=8. These are neither
peak-RSS increases nor arithmetic-time estimates. Middle-layer arithmetic may
pay more than this cumulative slot ratio suggests.

## Current native-key map budget

Five-repeat interleaved measurements compare benchmark modes 4 and 12.
The table below reports mode 12, the current native-key production policy.
The map timer covers support/transition construction; it excludes row support
extraction, eligibility checks and output arithmetic. Fractions and ceilings
are computed per run and then medianed, so ratios of the displayed stage
medians need not equal the displayed fraction.

| n | Threads | Total median s | Map median s | Median map share | Zero-map ceiling |
|---:|---:|---:|---:|---:|---:|
| 7 | 1 | 0.10516 | 0.04839 | 48.30% | 1.93x |
| 7 | 4 | 0.10180 | 0.05920 | 58.15% | 2.39x |
| 8 | 1 | 0.92211 | 0.56093 | 59.08% | 2.44x |
| 8 | 4 | 0.63261 | 0.27023 | 45.51% | 1.84x |

The zero-map ceiling is 1/(1-f), assuming all measured map cost disappears and
all other costs remain unchanged. It is not a measured optimization speedup
or an upper bound on changes that also accelerate arithmetic. Runs show
noticeable machine-load variation, especially at n=7; these are local budget
estimates, not portable performance guarantees. No removal of the subset-DP
exponential factor is claimed.

For illustration, halving map cost at n=8 while leaving everything else
unchanged would yield about 1.42x at one thread and 1.29x at four threads.
Quartering it would yield about 1.80x and 1.52x. Actual direct indices have
construction/lookup costs and padding overhead, which must be measured.

Unlike the backward-pruning experiment's roughly 0.13% arithmetic opportunity,
this is a substantial enough budget to justify an experimental direct-index
backend. A useful first version would retain a transition map, replacing only
hash-based construction. Removing the map itself is a separate tradeoff:
recomputing indices for every minor could lose more than it saves.

## Reproduction and validation

```sh
make build/mq_layout_bench
python3 src/test/mq_layout_bench.py maps-policy --repeats 5 --output src/test/data/mq_layout/native_rank_budget.json
python3 src/test/mq_rank_budget.py --supports src/test/data/mq_layout/support_backward.jsonl --timings src/test/data/mq_layout/native_rank_budget.json --output src/test/data/mq_layout/rank_budget.json
```

Generate `support_backward.jsonl` first using `MQ_SUPPORT_COUNT.md` if absent.
The rank-budget script checks n=4 envelope counts by explicit enumeration and
checks that the rank formula is a bijection onto contiguous indices. Those
checks passed. Forty timing runs completed, with the benchmark driver's
canonical fingerprints matching between serial and production map policies.
No new coefficient algorithm was implemented or claimed validated. Generated
JSON files remain in the ignored local measurements directory.

## Integrated direct-index experiment

The solver accepts `--mq-step1-rank` / `--no-mq-step1-rank`, now default on
for eligible projected shared DP. The measurements below predate this default
change and used explicit benchmark modes.
The option operates inside shared projected DP; disabling sharing, disabling
projection, or selecting an explicit different determinant backend preserves
those choices. There is no n-specific speed gate. Existing shared-backend
packing, filter and workspace limits still apply; this does not lift their
size limits. Full, unprojected determinants continue using hash construction.

`mq_rank_layout.h` enumerates axis ideals using native multiword keys, sorts
them for binary lookup, and assigns a base offset to every admissible x/y
pair. The parameter exponent is an offset within that pair's block. It
precomputes per-axis shifted indices, then builds the existing full transition
map by combining those small indices. The coefficient recurrence is unchanged;
unreachable envelope entries start at zero. Native full keys remain available
for subsequent fallback layers and final FLINT packing.

Degree limits come from the actual maximum degrees of the preceding support
and current row shifts, not an assumption that every row is fully dense.
The parameter must already occur in the processed support. The ordinary
preflight simplex bound remains valid for the padded envelope. A second
256 MiB workspace check includes axis tables, pair offsets, source identifiers,
axis transitions, previous/new coefficients, map, keys, factors, and output
packing allowance. Failure leaves the destination untouched and falls back to
hash construction. This is a workspace estimate, not a process-RSS guarantee.

Five-repeat interleaved comparison, F65537, seed (132,941), same projected
fixture and timing boundaries as the feasibility experiment. Benchmark mode 13
forces direct indices, mode 14 forces the previous production hash policy;
mode 12 continues to follow production configuration. Fingerprints agree.

| n | Threads | Hash total s | Direct total s | Speedup | Hash map s | Direct map s |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 1 | 0.09357 | 0.07138 | 1.31x | 0.03996 | 0.00659 |
| 7 | 4 | 0.09690 | 0.05043 | 1.92x | 0.04749 | 0.00623 |
| 8 | 1 | 1.04091 | 0.45884 | 2.27x | 0.59814 | 0.04314 |
| 8 | 4 | 0.62998 | 0.34018 | 1.85x | 0.27671 | 0.03928 |

At n=8/4 threads, median process peak RSS changed from 117.58 MiB to 109.41 MiB;
at one thread both were about 109.58 MiB. These peaks include the complete
projection invocation, not just live DP arrays. n=7 timings show noticeable
load variation; numbers are local results, not guarantees for other inputs or
larger n. Stage medians do not sum to total medians. The method reduces index
construction overhead without removing the exponential subset-DP factor.

Validation completed:

* 17 exact sparse-DP comparisons: projected n=7/8 at two seeds and 1/4 threads,
  full-output fallback, projected n=5 over F2/F7/a near-full-word prime, and
  projected n=8 with 16 threads. Usage counters confirm direct layers ran in
  projected cases; full-output cases used fallback. Small-field cases include
  the existing candidate-repair flow.
* Private transition checks for native one-, two-, and three-word layouts,
  comparing every transition against serial hash construction; direct-to-hash
  and hash-to-direct transitions; rejection of sources outside the ideal;
  memory rejection preserving destination and map; missing-filter fallback.
* 21 CLI checks for exact output equality, activation, opt-out ordering,
  disabled sharing/filtering, explicit pencil precedence, n=7/8 and 1/4 threads.

Reproduce:

```sh
make test-mq-rank test-mq-rank-cli
python3 src/test/mq_rank_bench.py timing --repeats 5 --output src/test/data/mq_layout/rank_timings.json
```

Local raw audit and timing files are `src/test/data/mq_layout/rank_audits.json`
and `src/test/data/mq_layout/rank_timings.json`; both are ignored by Git.

After enabling the default, the CLI suite was expanded to 28 checks including
flag-free one- and four-thread runs, explicit enable/disable, and precedence.
All passed with exact output comparisons. Explicit hash benchmark modes remain
hash-only; mode 12 follows the updated production default.

## Removal of the fixed workspace cap

The default `DRSOLVE_MQ_SHARED_WORKSPACE_BYTES` is now `SIZE_MAX`, removing
both the shared-DP and direct-index 256 MiB admission gates. Index-width,
packing, input-shape and addressability checks remain; estimates accumulate
without size_t summation overflow. An explicit compile-time budget override
is still supported. This does not guarantee sufficient physical RAM. Earlier
measurements and validation above describe the former 256 MiB policy.
The benchmark driver now accepts n=9 and uses a sufficiently large exponent
buffer for its fingerprint, including the 17-variable native layout.

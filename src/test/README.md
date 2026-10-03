# Tests

All test cases and benchmarks live in this directory, grouped by module.
Random-polynomial generation APIs remain in their implementation modules.

With Make:

```sh
make test-minor-dp
make test-components
./build/component_tests --list
./build/component_tests rational-solver
```

With CMake (`BUILD_TESTING=ON`, the default):

```sh
cmake --build <build-dir> --target det_minor_dp_test component_tests
ctest --test-dir <build-dir> -R 'drsolve_(minor_dp|component_)' --output-on-failure
```

`det_minor_dp.c` contains automated determinant and scheduling regressions.
Its executable compiles both determinant backends with test-only observations
to check their actual parallel scheduling and bounded-memory fallback.
The `*_test.c` files extracted from implementation modules retain their
original diagnostic output; many are examples or benchmarks rather than
assertion-based tests. The component runner exposes them individually.
The sparse interpolation case and roots benchmark are opt-in because they
can be expensive. Existing CLI suites such as `--test 1` and `--test 5`
also perform large computations and are not part of the smoke tests above.

`dixon_test.c` and `polynomial_system_solver_test.c` remain linked into the
library to preserve existing CLI test entry points and exported functions.
The former also provides random-system generation used by the CLI's `-r`
mode. The existing `./drsolve --test <n>` interface is unchanged.

`make test-mq-layout` checks the default MQ shared-support DP backend against
the existing determinant backend. Its n=7/8 layer statistics, timing ablations,
row-order comparison, raw data, and reproduction commands are documented in
[MQ_LAYOUT_EXPERIMENT.md](MQ_LAYOUT_EXPERIMENT.md). Private profiling and
ablations are compiled only into `build/mq_layout_bench`; the production
backend is default on with `--no-mq-step1-shared` as an explicit opt-out.
`make test-mq-shared-cli` compares complete solver outputs and flag precedence.

`make test-prime-step4` checks the bundled PML prime-field NTT products,
coefficient windows, aliasing, fallback, geometric input bounds and degree-order
determinants. With CMake, build `pml_prime_step4_test` and run
`ctest --test-dir <build-dir> -R drsolve_prime_step4 --output-on-failure`.
`pml_prime_step4_bench.py` compares complete CLI resultants and reports Step 4
medians, alternating baseline/optimized runs. It accepts degrees, seeds, prime,
repeat count, optional ablations and a JSON output path. The backend switches
are documented in [the bundled PML README](../../pml_det/README.md).

`make test-poly-mat-interpolation` covers both automatic and forced batch
matrix evaluation. It compares individual entries against scalar FLINT
evaluation, including rectangular matrix windows, zero/constant entries,
partial final batches, small fields with repeated points, word-size primes,
input preservation and parallel readers. Full determinant tests compare
against FFLU, retain the existing small-field value-polynomial semantics,
and verify rejection of an incorrect degree bound. CMake provides
`drsolve_poly_mat_interpolation_auto`, `_scalar`, and `_batch`.

The interpolation backend packs entry coefficients once, groups them into up
to eight degree ranges, and multiplies batches of ordinary point powers by
those coefficient blocks using FLINT. Results are stored by point so that
assembling each numeric matrix reads contiguous data. This uses no special
modulus, roots of unity or assumptions about the number of equations. Each
batch has at most 128 points; coefficient packing is limited to 64 MiB and
explicit worker buffers to approximately 192 MiB in total. Tiny problems or
plans exceeding these limits fall back to scalar evaluation. These limits do
not include FLINT's internal multiplication workspace or the input matrix.

`DRSOLVE_INTERP_EVAL=auto` is the default. Use `scalar` to reproduce the old
point-by-point evaluator, or `batch` to bypass the small-shape threshold
(the memory limits still apply). At verbosity 2, the additional evaluation
line separates coefficient packing, matrix evaluation and numeric determinants.
Worker times are summed across workers, while the existing interpolation line
continues to report wall time. This dispatch applies to `--fq-det-method interp`;
it does not change the selection of a determinant backend.

For alternating scalar/batch timing runs with complete output comparison:

```sh
python3 src/test/poly_mat_interpolation_bench.py --degrees 12 16 \
    --primes 65537 1000003 998244353 --repeats 3 --json /tmp/interpolation.json
```

On the development machine, degree-16 three-equation systems, `--seed 12345`,
`--threads 1 --fq-det-method interp`, three alternating runs per variant gave:

| Prime | Scalar Step 4 median | Batched Step 4 median | Time reduction |
| --- | ---: | ---: | ---: |
| 65537 | 13.882 s | 6.698 s | 52% |
| 1000003 | 14.071 s | 6.172 s | 56% |
| 998244353 | 17.564 s | 9.348 s | 47% |

All 36 degree-12/16 benchmark runs matched the complete scalar resultants.
These timings compare the interpolation backend, not the default HNF backend;
they depend on the machine and input. No modulus-specific arithmetic was used.

## Partial Dixon relations and bounded degree closure

`make test-dixon-closure` builds a standalone FLINT experiment and runs 330
exact native regressions. This does not change the production solver or its
dispatch. The question is whether generating some Dixon relations reduces the
total cost of closing the relations from an overdetermined quadratic system.

The same incremental finite-field RREF closure is used for all six methods:

| Method | Initial relations |
| --- | --- |
| `seed` | All original equations |
| `hybrid3` | Original equations plus Dixon coefficients of degree at most 3 |
| `hybrid4` | Original equations plus Dixon coefficients of degree at most 4 |
| `macaulay` | All original-equation multiples through total degree D |
| `augment3` | The full Macaulay input plus the degree-3 Dixon coefficients |
| `augment4` | The full Macaulay input plus the degree-4 Dixon coefficients |

For the Dixon construction, the first n equations are used and the last of n
variables is retained. The adjacent divided-difference rows precede the original
row. Subset determinant DP discards terms that cannot reach the requested
auxiliary-degree band during construction; it does not first build the full
Dixon polynomial. The highest auxiliary-degree band is omitted because it only
gives constant linear combinations of the original equations. Both prototypes
then use all n variables as scalar monomial columns, rather than working over
a univariate polynomial coefficient ring.

Closure reduces incoming rows against existing pivots, echelonizes the residual
on free columns, and prolongs only previously unprocessed low-degree relations.
The full Macaulay variants mark the original multiples through degree D-1 as
already prolonged. There are no field equations, Gröbner basis calls, or planted
coordinates inside the algorithms. A unique point is returned only after all
coordinate linear relations have been obtained and the point passes substitution
in every input equation. For other cases, the output is a relation space, not a
claim to have enumerated all solutions.

There is a useful limitation on these specific truncations. Laplace expansion
along the original row gives `g_beta = sum_i a_(i,beta) f_i`, where
`deg(a_(i,beta)) <= n-1-|beta|`. Selecting `|beta| >= n+1-d` for d=3 or 4 thus
gives multiplier degree at most d-2. When D>=d, every selected relation already
belongs to the degree-D Macaulay row space. It cannot increase that initial
rank or change the saturated degree-D closure. It could still accelerate how
the space is generated; the experiment measures that possibility, including
the cost of constructing the extra relations.

Reproduce a Sage/native paired experiment:

```sh
make test-dixon-closure
DOT_SAGE=/tmp/drsolve-closure-sage OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  sage -python src/test/dixon_closure_experiment.py \
    --variables 6 --extra 1 --degree 5 --seeds 5 --repeats 3 \
    --methods seed,hybrid3,hybrid4,macaulay,augment3,augment4 \
    --export /tmp/closure-inputs --c-binary ./build/dixon_closure_bench \
    --output /tmp/closure-results.json
```

Add `--audit` for exact Macaulay row membership, equality of final Sage row
spaces, and comparison of truncated coefficients with the full determinant
for n<=5. Native and Sage runs also compare rank histories, recovered points,
canonical row-space fingerprints, and work counters. The native regression
suite covers primes 2, 3, 7, 257 and 65537, two planted points, inconsistency,
zero/duplicate equations, and exhaustive small-field solution checks. The
optional final argument to the native executable dumps the entire RREF matrix
for coefficient-by-coefficient cross-checking.

The recorded inputs, raw paired runs and commands are in
[`data/dixon_closure`](data/dixon_closure). Regenerate the independent native
timing pass and its summary without Sage:

```sh
python3 src/test/dixon_closure_native_bench.py
python3 src/test/dixon_closure_summary.py
```

The independent C pass shuffles method order, warms up each native process,
then records five samples of both wall and process CPU time. Summary values
take medians within each seed before taking a median across seeds. Timings
include basis setup, partial Dixon construction, row assembly, RREF, closure,
and single-point recovery/verification; input parsing and optional audits are
outside the timer. RREF cell counts and dense matrix multiply-add counts are
work proxies, not measurements of FLINT's actual field-operation count. The
Dixon term-product count is a naive convolution work estimate, not a claim
that FLINT uses naive polynomial multiplication.

On the recorded F257 dense planted systems, the independent single-thread C
pass gave the following wall-time medians (milliseconds, medians within each
seed followed by a median across seeds):

| Variables / equations | D | Seeds | `seed` | `hybrid3` | `hybrid4` | `macaulay` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 5 / 6 | 4 | 6 | 1.573 | 2.820 | 3.553 | 1.212 |
| 6 / 7 | 5 | 5 | 40.945 | 53.171 | 64.308 | 34.670 |
| 7 / 8 | 5 | 3 | 238.556 | 272.234 | 314.510 | 190.739 |
| 7 / 10 | 4 | 5 | 16.574 | 45.915 | 99.904 | 10.454 |

These truncations did not provide a stable total-time improvement. Process CPU
times support the same conclusion. Full-Macaulay augmentation also increased
the paired median cost. The three-seed double-point and shared-quadratic-part
tests likewise gave no aggregate advantage. These are small finite-field
experiments, not a universal claim about every choice of Dixon relations.

For example, the six-variable degree-5 `seed` and `hybrid3` runs both reach
rank 461 in 462 columns. Their rank histories are respectively
`7,49,175,441,455,461` and `35,151,411,441,455,461`. Both have six recorded
stages; the higher initial rank does not save a stage. Submitted rows increase
from 1225 to 1267, RREF cell volume from 633570 to 686346, and the dense
matrix-multiplication work proxy from 109983342 to 117367326, before charging
any partial Dixon construction.

At the next lower tested degree bounds, all variants stop at the same remaining
dimension: 20 for 5 variables at D=3, 35 for 6 variables at D=4, and 70 for
7 variables / 8 equations at D=4. No lower-degree success was obtained by
adding these coefficients. The double-point tests retain dimension 2 and do
not incorrectly report a unique point; this experiment does not time their
subsequent multi-root extraction.

[`summary.json`](data/dixon_closure/summary.json) includes paired per-seed speed
ratios, both timing sources, phase costs and work counts. The native-only pass
is [`native_paired.json`](data/dixon_closure/native_paired.json); the individual
`n*.json` experiment files retain the earlier Sage/C crosschecks. All 498
Sage/C comparisons match rank histories, work counters, recovered coordinates
and row-space fingerprints. The separate `audit_*.json` files additionally
record full RREF coefficient comparisons on representative 5-, 6- and
7-variable systems. The native 330-run regression suite also passes with
AddressSanitizer, LeakSanitizer and UndefinedBehaviorSanitizer enabled.

## Overdetermined route sweep

`dixon_over_routes.py` extends the fixed-band experiment with adaptive relation
spaces, certified quotient closure, shared cofactors, polynomial-module
compression, streaming core quotients, and block linearization. These benchmarks
are standalone experiments. Production dispatch was unchanged during the sweep;
the subsequent prime-field integration is described below.

```bash
make test-dixon-over
python3 src/test/dixon_over_suite.py
python3 src/test/dixon_over_native_bench.py
python3 src/test/dixon_over_oracle_native.py
DOT_SAGE=/tmp/drsolve-closure-sage OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  sage -python src/test/dixon_over_production.py
python3 src/test/dixon_over_summary.py
```

Run timing commands serially. The suite records exact commands, raw inputs,
input hashes, all successful and unsuccessful attempts, and per-stage results
under [`data/dixon_over_routes`](data/dixon_over_routes). It uses planted dense,
structured, shared-leading-part, two-point, and linearly redundant families over
F257. Inputs with different equation counts share the same prefix. Independent
GB audits establish the final quotient dimension and eliminant, so cases with
the same certified simple point have the same ideal despite different generating
sets. The redundant control adds only constant linear combinations of the first
n+1 equations. It changes neither that ideal nor any complete Macaulay row space.

Methods and scope:

| Name | What is timed |
| --- | --- |
| `scalar` | Original-only bounded closure from the preceding experiment |
| `quotient` | Incremental Macaulay degrees, closure in shrinking residual coordinates, and an algebra certificate |
| `hybrid4`, `hybrid5` | Low original relations, shared Dixon coefficient batches, and quotient/ordinary closure between batches |
| `hybrid_support4` | Initially add only the support of requested coefficients; a bounded-closure fallback can enlarge it |
| `shared` | Identical coefficient bands from several fixed-core subsystems, cached versus independently constructed |
| `module`, `module_units` | Original x-multiples over F257[t], with/without unit-pivot elimination before residual HNF |
| `module_dixon`, `module_unit_dixon`, `module_stream` | Add complete Dixon coefficient relations to a smaller original module, optionally in batches |
| `module_nf`, `module_nf_full` | Reduce inside cofactor DP using original-derived unit-leading polynomial rules, then certify the resulting module |
| `module_nf2_full` | The same with a larger original preprocessing space |
| `core_full`, `core_stream` | Build the square-core GB/normal forms, add extra-equation images all at once or in batches, then certify the quotient |
| `linear2` | Cancel a two-variable quadratic block, substitute, and compute both determinant branches using Singular GB |
| `oracle` | Offline selection of individual/batched coefficients and early/middle injection schedules, excluding generation and selection cost |

`module_nf_full` includes the full auxiliary-degree range (apart from the
constant input combinations), not merely the cubic/quartic band. Its internal
division rules generate a proven subideal J of the input ideal. Each division
changes a polynomial by an element of J, so arithmetic remains congruent modulo
J and its output coefficients remain valid ideal relations. This does **not**
assume that a linear residual projection is multiplicative, or that the rules
already form a GB. Reduced coefficients need not equal unreduced Dixon
coefficients; their ideal membership is independently audited.

The quotient certificate checks all multiplication-operator commutators, input
actions, and the cyclic basis generated by 1. It preserves multiple roots and
nonreduced algebras; it is not a nullity-stabilization heuristic. The exact
minimal polynomial of the last-coordinate operator is included. No separating
coordinate or radicality assumption is made. The `linear2` prototype certifies
a unique point only when its exceptional determinant branch is inconsistent;
other outcomes are reported as incomplete rather than silently discarding that
branch. This linearization experiment uses GB for the transformed equations,
not a new higher-degree Dixon implementation.

All preprocessing and certificates are charged. Independent GB audits and input
generation are outside timings. Sage runs warm the methods and shuffle order;
the main configurations use three seeds and three repetitions. Seven-variable
dense runs use one repetition; the full-degree seven-variable module supplement
uses two seeds. Deterministic failures at the specified cap and timeouts are
retained, without repeating their timings. `oracle` is only a finite diagnostic
search, not a globally optimal subset oracle or an implementable speedup.

`dixon_over_native.c` implements residual quotient closure and both shared
cofactor generation variants using native FLINT arithmetic. It reuses the
audited basis and closure kernel from `dixon_closure_bench.c`. The native pass
uses one warmup and five samples in each of three shuffled rounds. Fixed-degree
baselines receive an **untimed search for their smallest successful degree**;
the separate `adaptive` method charges its degree search. Native coefficient
equality and quotient/eliminant results are checked exactly. Native oracle
timings still exclude generation/selection; the generator's cache setup is
reported separately, not presented as an end-to-end oracle speedup.

For nonunique cases, `seed`/`macaulay` stop at a bounded relation space; the
native degree calibration only checks its nullity against the independent
quotient dimension. Those baseline timings do not include constructing a
complete multiplication representation or extracting multiple roots. Full
solve-time ratios are consequently restricted to the certified simple-point
cases. `quotient`/`adaptive` always require their complete algebra certificate.

The optional production comparison times **one square-core resultant** from
the existing CLI. It is a different output task from solving all overdetermined
equations and must not be used as a full-solver speed ratio. Both CLI wall time
and the reported internal step times are saved.

### Route sweep results (2026-10-03)

The sweep contains 60 distinct inputs, 1,641 Sage result records, and 3,780
measured native samples (756 per-process medians, excluding warmups). All 1,323
certified Sage results match the independent full-system GB oracle. Three
measured full-degree module attempts hit their 30-second limit; these are
retained as timeouts. Construction-only and offline-diagnostic records are not
counted as solver failures.

The clearest end-to-end gain is from **original equations alone**, using the
small quotient for closure. Adaptive C timings include the search over degrees
and the final algebra/point certificate. This native adaptive version rebuilds
the Macaulay space at each degree, charging every attempt.

| Variables | Equations | Final degree | Monomial columns | Adaptive C total, ms | Paired speedup |
| --- | --- | --- | --- | --- | --- |
| 5 | 6 → 8 | 4 → 3 | 126 → 56 | 0.761 → 0.139 | 5.55x |
| 6 | 7 → 9 | 5 → 4 | 462 → 210 | 18.997 → 1.928 | 10.25x |
| 7 | 8 → 10 | 5 → 4 | 792 → 330 | 73.345 → 5.279 | 13.78x |

Times are medians of per-seed medians; speedups are medians of paired seed
ratios. Both systems in each pair define the **same reduced point ideal**, as
checked by dimension-one certificates and identical coordinates. These gains
therefore do not merely reflect discarding some solutions of a square core.
In contrast, adding only constant linear combinations left the final degree
at 4 and the five-variable time around 0.76 ms.

For the same weakly overdetermined inputs and the same independently selected
best degree, residual quotient closure also improves the existing full-Macaulay
closure baseline: 1.194 → 0.674 ms (5/6), 29.918 → 15.134 ms (6/7), and
109.702 → 66.368 ms (7/8). These are comparisons to the supplied relation-closure
implementation, not to an optimal GB solver.

Other routes have narrower benefits:

* Native shared coefficient construction is 1.66–1.75x faster when producing
  16 fixed-core subsystem coefficient bands (6/21 and 5/20). With only two
  subsystems it is 1.4–2.3x slower. This is a construction-stage improvement,
  not a measured full-solver improvement.
* Full-degree unit-rule reduction **inside** Dixon construction solves both
  7-variable/8-equation module supplement cases in 8.52 and 11.39 seconds,
  versus 16.34 and 22.49 seconds for the degree-three original-multiplier
  module. It remains slower than original-only scalar/quotient closure in the
  same Sage environment. The seven-variable supplement has only two seeds.
* Streaming extra-equation images sometimes improves the core-quotient route
  (e.g. structured 6/9: 35.89 → 26.75 ms), but is not consistently faster.
  Core GB preparation is included, and solving the full system with Singular
  remains faster than either core-first route on these tested configurations.
* Two-variable block linearization typically raises the remaining polynomial
  degree from 2 to 6. Both determinant branches are computed; it does not beat
  the full-system GB baseline here, and the two-point collision cases do not
  obtain its specialized unique-point certificate.
* The offline coefficient search can reduce the closure work proxy by roughly
  5–16% in some configurations, but the C savings are small. For 5/6 it reduces
  fixed-degree closure from 1.474 to 1.395 ms with 24 precomputed relations;
  just building this generator's cofactor cache costs about 2.775 ms. For 6/9
  the corresponding saving is about 0.82 ms against a 29.5 ms cache cost.
  These cache comparisons concern this implementation, not a lower bound on
  all possible coefficient generators. The fixed-degree oracle also fails to
  beat the adaptive original-only baseline on the more overdetermined cases.

No tested mixed Dixon route improves on the best original-only closure
baseline after all construction costs are charged. This does not rule out a
different coefficient-selection or structured-generation algorithm, and the
small n=5–7 experiments do not establish asymptotic complexity bounds.

`make test-dixon-over` checks 50 zero-dimensional/inconsistent systems plus
positive-dimensional rejection fixtures: 135 native calls, 32 exact shared
coefficient comparisons, 100 polynomial-module checks, and 50 internally
reduced coefficient checks. Fixtures include projection collisions and
`(x^2,(t-x)^2)`, whose retained-coordinate minimal polynomial distinguishes
characteristic 2 from odd characteristic. The native checks also pass with
AddressSanitizer and UndefinedBehaviorSanitizer; leak detection is disabled
for the sandboxed sanitizer run. See `audit.json`, `asan_audit.json`,
`summary.json`, and `manifest.json` in the data directory.

## Production prime-field quotient solver (2026-10-04)

`src/solver/quotient_solver.c` ports the original-only adaptive native route
into the solver library. It accepts general polynomial degrees and dynamic
monomial supports; the standalone experiment's eight-variable/six-degree caps
are not inherited. The backend has no dependency on test sources.

The initial Macaulay reduction and shrinking quotient closure use nmod linear
algebra. When the low basis does not span the residual space, bounded
prolongation is attempted before increasing degree. A successful result checks
commutation, cyclicity and every original polynomial's action. Rational points
are then extracted through successive common eigenspaces of all multiplication
operators, so a separating last coordinate and a radical ideal are unnecessary.
Each extracted point is verified against every parsed input equation.

`--solver auto` chooses this backend for overdetermined systems over machine-word
prime fields. `--quotient`/`--solver quotient` also enables it on square systems;
`--solver dixon` requests the previous route. Other coefficient domains retain
the previous solvers. Quotient failure is reported as incomplete, with a nonzero
CLI exit status; it is never silently reclassified as inconsistency or positive
dimension. Random input preserves the declared variables, including variables
absent from all generated polynomials. Explicit resultant mode is separate.

```bash
make test-quotient
make test-quotient-sage
```

`quotient_cli_test.py` compares complete point sets against exhaustive searches,
checks automatic/explicit dispatch, file and random input, reproducible seeds,
resource/degree failures, unsupported fields, sparse missing-variable inputs,
and cases beyond the old experimental caps. `quotient_solver_test.py` uses
independent Singular GBs to audit quotient dimension, the exact elimination
polynomial (including nonreduced structure), and all rational points; it also
checks positive-dimensional rejection and a prime near the word-size limit.
The driver can be compiled with `-DNDEBUG` or sanitizers and passed through
`--binary` to test those builds against the same independent oracle.

Validation on 2026-10-04: 35 exhaustive CLI point-set comparisons plus dispatch
and error-path checks passed; 70 production-backend cases passed the independent
algebra/point audit (the near-word-limit linear fixture uses its known exact
algebra because Sage 9.5's Singular backend does not support that characteristic).
The same 70 cases passed with `-DNDEBUG` and AddressSanitizer/UndefinedBehaviorSanitizer;
leak detection was disabled for the sandbox sanitizer run. All 15 existing
`make check` regressions also passed.

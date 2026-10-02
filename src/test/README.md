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

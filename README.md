# DRSolve: Dixon Resultant & Polynomial System Solver
A C implementation for computing Dixon resultants and solving polynomial systems over finite fields and the rationals ℚ, based on the FLINT and PML libraries.

Website: <https://drsolve.github.io>

Author: Haohai Suo (<haohai.suo@mail.sdu.edu.cn>)

## Features
- Dixon resultant computation for variable elimination
- Polynomial system solver
- Finite fields:
  - Prime fields F_p (any size): Implemented with FLINT modular arithmetic, optionally accelerated by PML.
  - Extension fields F_{p^k}: Further optimized for binary fields F_{2^n} with n in {8, 16, 32, 64, 128}.
- Rational field ℚ: Rational reconstruction via multi-prime CRT. Set field_size = 0 to enable.
- Complexity analysis — estimates Dixon matrix size, Bezout degree bound, and operation count before computing

---

## Dependencies
- **FLINT** (recommended version: 3.6.0)  
  <https://github.com/flintlib/flint>

```bash
sudo apt install libgmp-dev libmpfr-dev make autoconf libtool-bin
git clone https://github.com/flintlib/flint.git && cd flint
./bootstrap.sh
./configure 
make
sudo make install
```
  
- **PML** (built in)  
  <https://github.com/vneiger/pml>

---

## Build
```bash
git clone https://github.com/drsolve/drsolve.git && cd drsolve
./configure
make
make check                         # optional
sudo make install                       # optional
```
For more options, run `./configure --help` or `make help`.

We also provide a Windows GUI, which can be built with CMake.
```bash
cmake -B build-win -DCMAKE_TOOLCHAIN_FILE="$(pwd)/cmake/toolchain-mingw64.cmake"
cmake --build build-win -j$(nproc)
```
---

## Usage
### BASIC USAGE

#### Elimination / resultant mode
```bash
./drsolve "polynomials" "eliminate_vars" field_size
```
Example:
```bash
./drsolve "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
```
- Default output file: `out/solution_YYYYMMDD_HHMMSS.dr`

#### Polynomial system solver
```bash
./drsolve "polynomials" field_size
```
Example:
```bash
./drsolve "x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1" 0
```
- Writes all solutions to `out/solution_YYYYMMDD_HHMMSS.dr`

### FILE FORMAT
#### File input/output
```bash
./drsolve input_file
./drsolve -f input_file -o output_file
```
#### Dixon resultant elimination (multiline)
``` 
Line 1 : variables to ELIMINATE (comma-separated)
Line 2 : field size (prime or p^k; use 0 for Q; generator defaults to 't')
Line 3+: polynomials (comma-separated, may span multiple lines, #eliminate = #equations - 1)
```
Example:
```bash
x0,x1
257
x0^3+x1^3+x2^3, x0*x1+x1*x2+x2*x1, x1*x2*x0+1
```
Run:
```bash
./drsolve examples/example.dr
./drsolve -f examples/example.dr -o out/my_result.dr
```
- If line 1 lists `n` variables for `n` equations, compatibility mode uses the first `n-1` variables

#### Polynomial solver mode (multiline)
```
Line 1 : field size
Line 2+: polynomials (comma-separated, may span multiple lines)
```
Example:
```bash
0
x^2+y^2+z^2-6, x+y+z-4, x*y*z-x-1
```
Run:
```bash
./drsolve examples/example_solve.dr
./drsolve -f examples/example_solve.dr -o out/my_solutions.dr
```

### OTHER MODES

#### Extension fields
```bash
./drsolve "x + y^2 + t, x*y + t*y + 1" "y" 2^8
```
The default settings use `t` as the extension field generator and FLINT's built-in field polynomial.
```bash
./drsolve "x^2 + t*y, x*y + t^2" "2^8: t^8 + t^4 + t^3 + t + 1"
```
- Example: AES polynomial for `GF(2^8)`
- In `Q` and prime fields, `t` is treated as an ordinary variable; only extension fields reserve it as the generator

#### Complexity analysis
Estimates the difficulty of a Dixon resultant computation without performing it.
Reports equation count, variable count, degree sequence, Dixon matrix size,
Bezout degree bound, and complexity in bits.

```bash
./drsolve -c "polynomials" "eliminate_vars" field_size
./drsolve -c -f input.dr
```
- Prints complexity information
- For n quadratics in n variables with one remaining parameter, `-v 2` also
  compares the legacy cached Laplace surrogate with total-degree and per-layer
  MQ bounds. Step 1 selects the smallest estimate including total-degree simplex
  interpolation, charging entry evaluation, numerical determinants, Newton
  transforms, and any extension-field arithmetic. Step 4 compares rank-sized
  methods with blocked Schur formation, the core determinant, and verification;
  its MQ entry-degree bound is `n+1`. The overall estimate is
  `max(Step 1/2, Step 4)`, with no separate Step 3 extraction charge.
  `-v 1` retains the original method selection. For example:
  `./drsolve -c -r '[2]*10' 257 -v 2`.
- Default output file: `out/comp_YYYYMMDD_HHMMSS.dr`
- Add `--omega <value>` or `-w <value>` to set the matrix-multiplication exponent

Example:
```bash
./drsolve -c "x^3+y^3+z^3, x^2*y+y^2*z+z^2*x, x+y+z-1" "x,y" 257
```

#### Random mode
Generate random polynomial systems for testing and benchmarking.

```bash
./drsolve -r       "[d1,d2,...,dn]" field_size
./drsolve -r       "[d]*n"          field_size
./drsolve -r -n 4 --density 0.5 "[d]*3" field_size
./drsolve -r -s    "[d1,...,dn]"    field_size
./drsolve -r -c    "[d]*n"         field_size
```
- Add `-n <num_vars>` to set the variable count
- Add `--density <ratio>` with `0 <= ratio <= 1`
- Add `--seed <num>` for reproducible random systems
- Mixed degree specifications such as `"[2]*5+[3]*6"` are supported

Examples:
```bash
./drsolve -r "[3,3,2]" 257
./drsolve -r "[3]*3" 0
./drsolve -r -n 4 --density 0.5 "[3]*3" 257
./drsolve -r --seed 12345 "[3]*3" 257
./drsolve -r "[2]*4+[3]*2" 257
./drsolve -r -s "[2]*3" 257
./drsolve -r --comp --omega 2.373 "[4]*4" 257
```

#### Dixon with ideal reduction
```bash
./drsolve --ideal "ideal_generators" "polynomials" "eliminate_vars" field_size
./drsolve --ideal -f input.dr -o output.dr
```
- `ideal_generators` is a comma-separated list of relations with `=`
- In file mode, lines after the first two lines containing `=` are treated as ideal generators

Example:
```bash
./drsolve --ideal "a2^3=2*a1+1, a3^3=a1*a2+3" "a1^2+a2^2+a3^2-10, a3^3-a1*a2-3" "a3" 257
```

#### Field-equation reduction mode
After each multiplication, reduces `x^q -> x` for every variable.

```bash
./drsolve --field-equation "polynomials" "eliminate_vars" field_size
./drsolve --field-equation -r "[d1,d2,...,dn]" field_size
```

Example:
```bash
./drsolve --field-equation "x0 + x0*x2, 1 + x1, x1 + x0*x1" "x0,x1" 2
./drsolve --field-equation -r "[3]*5" 3 --density 0.5
```

### OPTIONS

#### Method selection
```bash
./drsolve --method <num> <args>
./drsolve --step1 <num> --step4 <num> <args>
```
- Available methods: `0.Recursive`, `1.Kronecker+HNF`, `2.Interpolation`, `3.Sparse interpolation`, `4.Bareiss`, `5.Recursive Dixon construction`
- `--method` sets both Step 1 and Step 4 for backward compatibility

#### Experimental construction of the MQ Schur core

A standalone prototype requests coefficient panels from projected Step 1,
solves complementary degree blocks, and accumulates the exact Schur core
without constructing the full candidate matrix first:

```bash
make test-mq-direct-core
# n, threads, use-pencil (0=minor DP, 1=pencil), seed, prime
./build/mq_direct_core_test 8 4 0 12345 65537
```

The executable then constructs the ordinary candidate independently and checks
every core coefficient and the determinant multiplier. Small cases also compare
full determinants and exercise rejection without publishing a partial core.
This is an experiment, not a solver option: it still generates all complementary
coefficients and repeats projected determinant work across panels. Local n=8
measurements are slower than ordinary construction plus Schur compression.
It does not establish a smaller asymptotic Step 1 bound or an algebraic reduction
that bypasses complementary coefficients. Details and measurements:
[direct-core experiment](paper/rank/DIRECT_CORE_EXPERIMENT.md).

#### Experimental MQ Step 1 degree recurrence

```bash
./drsolve --mq-step1-pencil --threads 4 -f input.dr -v 2
```

`--mq-step1-pencil` enables a Faddeev-LeVerrier degree recurrence instead of
subset-minor DP. It normalizes the parameter coefficient matrix by constant
column operations to obtain a lower block `tI+L`, then constructs the determinant
and adjugate border terms using one polynomial matrix and per-thread column
scratch. Boundary vectors are formed from the current matrix.
Columns are overwritten only after their old entries have all been consumed;
the trace and diagonal correction follow the column-update barrier.
The constant determinant scale is preserved exactly. Parameter coefficients
are accumulated in separate degree buckets. With the default MQ prediction,
intermediate matrices are projected to the candidate's downward closure using
the shared packed coefficient filter. Final contributions are filtered
before bucket accumulation. Candidate verification and degree-aware Schur repair
remain enabled; missing repair strips are computed by projected minor DP.

This option is off by default; `--no-mq-step1-pencil` disables it. It currently
requires prime characteristic `p > n-1`, one parameter, and full row rank of the
linear rows' parameter coefficient matrix (`n` is the equation count). Other
explicit Step 1 backends take precedence. Unsupported inputs fall back to the
existing backend. The conservative eligibility estimate retains the old
two-matrix bound of 268,435,456 coefficient slots; this permits `n<=10`, while `n>=11` currently
falls back. This is an eligibility bound, not a process-memory limit.
`--no-mq-step1-filter` computes the complete pencil determinant.
If candidate repair fails, the complete determinant is computed as a fallback.
This remains an experimental backend; see the paired timings in the research
note before selecting it for performance.

`--threads` parallelizes independent columns and border products, using
at most `n-1` worker threads. `-v 2` reports normalization, recurrence, assembly,
and the largest sampled matrix-plus-column-scratch term count (excluding
vectors, multiplication temporaries, and retained allocation capacity).
When both experimental Step 1 options are supplied, the last enabling option
(`--mq-step1-pencil` or `--mq-step1-simplex`) wins. Run
`make test-mq-pencil test-mq-pencil-cli` for exactness and fallback tests.
See [the degree-DP research and measurements](paper/rank/STEP1_DP_STRUCTURE.md).

#### MQ Step 4 Schur compression
```bash
./drsolve --mq-step4-schur -f input.dr -v 2
```
Checked complementary-block compression is enabled by default for
prime-field MQ systems with one retained parameter. It applies to automatic
Step 4 and `--step4 1`; the smaller determinant uses the existing backend,
including `--fq-det-method`. Other explicit Step 4 methods are preserved.
`--no-mq-step4-schur` disables it explicitly; `--mq-step4-schur` re-enables it.

The compressor factors each constant degree-diagonal block once, solves
`E X = V` by block back substitution, then forms `A - U X`. It does not
copy or update the entire polynomial matrix for scalar pivot elimination.
Prime-field, single-parameter MQ Step 1 uses shared monomial indices and
coefficient arrays by default when the matrix shape, exponent packing and
index bounds allow it. The unique quadratic row remains last in the DP evaluation order. There
is no fixed workspace cap by default. Preflight retains integer/index and
addressability bounds for arrays, maps and support tables. Unsupported shapes,
packing or index bounds fall back to the existing sparse minor DP. Large inputs
can require substantial memory; no physical-memory availability check is made.
Builds may override `DRSOLVE_MQ_SHARED_WORKSPACE_BYTES` to impose a budget.

For all admitted sizes, shared indices use native FLINT exponent packing,
including multiword keys. Transition maps with at least 1,000,000 pairs can use
independent hash shards to build their indices in parallel. There is no n-based
performance gate; the shape, packing and workspace admission checks above still
apply. Shared indices avoid duplicating monomial indexing across minors, but do
not remove the exponential subset count of the DP.

Direct indices replace hash-based support construction by default for eligible
projected shared Step 1 DP: disable with `--no-mq-step1-rank`, re-enable with
`--mq-step1-rank`. They use separate x/y axis transitions and
parameter offsets with native multiword exponent keys, retaining the full
transition map. The shared backend's eligibility limits still apply; missing
projection filters, absent parameter support, or excess workspace fall back
to ordinary hash construction. `--no-mq-step1-shared` also disables this path.
See [the implementation measurements and checks](src/test/MQ_RANK_BUDGET.md).

For a verified Step 1 projection with one parameter, the prime-field resultant
pipeline constructs `nmod_poly_mat` directly. It packs the Dixon terms into
scalar records and releases the source terms and term maps before allocating
the matrix. Row/column content, degree ordering and Schur labels are preserved.
Step 4 stays in the native prime-field representation; successful Schur
compression releases the large matrix before taking the core determinant.
Public extraction calls and extension fields retain the non-consuming
`fq_nmod_poly_mat` path; other extraction cases keep the generic implementation.
With `-v 2`, direct Step 2 reports metadata, packing, source-release, matrix
initialization, filling and cleanup times separately.
See [direct Step 2 construction](src/test/MQ_STEP2_DIRECT.md).

Use `--no-mq-step1-shared` to select the previous sparse DP, or
`--mq-step1-shared` to re-enable sharing. This is independent of
`--no-mq-step1-filter`, and explicit Step 1 backends such as pencil, simplex,
interpolation and Bareiss keep precedence. Small characteristic is supported;
delayed modular reduction is used only when a proved one-limb bound is safe.
See [the n=7/8 measurements and production checks](src/test/MQ_LAYOUT_EXPERIMENT.md).
The sparse fallback retains its packed linear-product streams and coefficient
filtering.

Enable total-degree interpolation in Step 1 explicitly with:

```bash
./drsolve --mq-step1-simplex --mq-step4-schur --threads 4 -f input.dr -v 2
```

`--mq-step1-simplex` is off by default; `--no-mq-step1-simplex` disables it.
It uses `--threads` for matrix-entry evaluation, numerical determinants, and
independent interpolation fibers. Candidate verification and Schur repair reuse
the full interpolated polynomial, so failed selection does not recompute it.
`--no-mq-step1-filter` retains the full polynomial for the ordinary extraction
path. Other explicit Step 1 backends take precedence.

The current implementation requires a prime field with `p > n+1`, one retained
parameter, and MQ divided differences, where `n` is the equation count. It falls
back to the existing Step 1 backend for other inputs or when the dense simplex
workspace exceeds 16,777,216 points (covering through `n=9`). It does not yet
implement extension-field interpolation. Verbosity two reports each stage's
wall time and the fallback reason when applicable.

`make build/mq_simplex_bench` builds a full-coefficient differential benchmark;
for example `./build/mq_simplex_bench 8 65537 4` compares the production engine
with DP using four threads. The whole-minor merge experiments remain separate:
`make build/mq-sum-direct-cli` and `make build/mq-sum-products-cli`. Neither merge
kernel is enabled in the normal solver. See
[the experiment report](paper/rank/STEP1_EXPERIMENTS.md) for timings and commands.

The selected matrix's actual monomial labels are used after reordering or
repair. Degree checks and nonzero constant pivots certify each compression;
if the profile is unavailable or the complement is singular, the original
matrix goes to the normal determinant backend. This option also works with
`--no-mq-step1-filter`. It does not assume that every random MQ candidate has
an invertible complement.

#### Resultant construction
```bash
./drsolve --dixon <args>
./drsolve --macaulay <args>
./drsolve --subres <args>
```
- `--dixon`, `--macaulay`, and `--subres` are direct method selectors
- `--subres` is for exactly 2 polynomials and 1 elimination variable

#### Verbosity
```bash
./drsolve -v 0 <arguments>
./drsolve -v 2 <arguments>
./drsolve --time <args>

```
`-v 0` prints nothing but still writes the output file. `-v 1` is the default. `-v 2` restores the debug-level console output and timing. `-v 3` prints detailed profiling for Dixon construction. `-v 4` additionally prints the cancellation matrix, Dixon matrix, maximal-rank submatrix when each is <= 100 x 100. `--time` prints per-step timing

Example:
```bash
./drsolve -v 2 -f in.dr -o out.dr
```

#### Parallelism
```bash
./drsolve --threads <num> <args>
```
- Sets the number of threads for parallel computation

---

## SageMath Interface

`drsolve_sage_interface.sage` lets you call DRsolve directly from SageMath with Sage polynomial objects.

- Load the interface with `load("drsolve_sage_interface.sage")`, then set the binary path once with `set_drsolve_path("./drsolve")`.
- Main entry points:
  - `DixonRes(F, elim_vars, ...)` / `DixonResultant(...)`
  - `DixonSolve(F, ...)`
  - `DixonComplexity(F, elim_vars, ...)`
  - `DixonIdeal(F, ideal_gens, elim_vars, ...)`
- Common options include `field_size`, `verbosity`, `time`, `threads`, `debug`, `live_output`, `timeout`, and `output_dir`.
- Sage interface calls now keep autogenerated output files in `./out/` by default. Use `set_drsolve_output_dir("results")` once, `output_dir="results"` for one call, or `foutput="results/result.dr"` for an exact filename.
- `field_size` may be an integer prime, a string such as `"2^8"` or `"2^8: t^8+t^4+t^3+t+1"`, a Sage `GF(...)` object, or `0` for ℚ. If omitted, it is inferred from the Sage polynomial ring when possible.
- Resultants are returned as strings, so iterative elimination works naturally by feeding one `DixonRes(...)` output into the next call.
- For a fuller Sage reference with examples and options, see the top docstring in `drsolve_sage_interface.sage`.

---

## Development Notes
Parts of this project were developed with the assistance of AI-based coding tools. All AI-assisted contributions were reviewed and tested by the project author.

## License
DRSolve is distributed under the GNU General Public License version 2.0 (GPL-2.0-or-later). See the file COPYING.

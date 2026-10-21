# DRSolve: Dixon Resultant & Polynomial System Solver
A C implementation for computing Dixon resultants and solving polynomial systems over finite fields and the rationals ℚ, based on the FLINT and PML libraries.

Website: <https://drsolve.github.io>

Author: Haohai Suo (<haohai.suo@mail.sdu.edu.cn>)

## Features
- Dixon resultant computation for variable elimination
- Polynomial system solver
- Finite fields:
  - Prime fields F_p (any size): Implemented with FLINT modular arithmetic, optionally accelerated by PML.
  - Extension fields F_{p^k}: Further optimized for binary fields F_{2^n} with n in {4, 8, 16, 32, 64, 128}.
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

Add `--sift` to enable quotient-based removal of extraneous resultant factors.

Supported square MQ systems over prime fields (at least three variables, with
one variable retained) use Dixon quotient / Bezoutian construction by default.
Use `--no-quotient` to disable it in resultant and solving modes;
`--no-mq-quotient` remains an alias. Explicit method selections take priority
over automatic selection. `--quotient` explicitly enables the construction;
inputs rejected by its rank or degree checks fall back to Dixon.

Explicit `--quotient` also supports square systems over extension fields, for
quadratic and general positive degrees, with one variable retained. There is no
fixed equation-count or memory-budget cutoff; dimension and address-space
overflow checks reject only workspaces that cannot be represented.
MQ systems reuse the prime kernel's shared leading Macaulay frames and linear
multiplication-map recurrence. GF(16)/GF(256) additionally use compact byte minors,
SIMD field accumulation, and maps allocated one recurrence layer at a time.
General degrees use grouped recurrence with dynamically indexed monomials;
packed MQ keys that do not fit select this general representation.
Both paths use the native binary-field
backends for GF(2^4), GF(2^8), GF(2^16), GF(2^32), GF(2^64), and GF(2^128), and the
unified field interface for other extensions. `--threads` parallelizes samples.
The existing prime-field MQ fast kernel and automatic selection are unchanged.

```bash
./drsolve -r '[2]*8' '2^8' --quotient --seed 12345 --threads 4
./drsolve -r '[3]*6' '2^8' --quotient --seed 12345 --threads 4
./drsolve -r '[3]*3' '3^3' --quotient --seed 12345
make test-extension-quotient
```

For degrees `d_i`, the general-degree extension path needs at least `sum(d_i - 1) + 5`
distinct field elements; the MQ fast path needs `N + 4` for `N` equations. Both
include two validation points. They enumerate elements
in a polynomial basis rather than embedding integer constants in the prime
subfield. Insufficient field size, nongeneric ranks, changing quotient bases,
or failed interpolation/determinant degree checks cause a diagnostic fallback
to Dixon. This construction is distinct from the `--solver quotient` closure
solver and does not use its degree or memory budgets.

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
- Fewer variables than equations automatically selects full-system solving
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

Finite-field elimination of two variables from three equations in three
variables defaults to method 5. Explicit method options take precedence;
use `--dixon` to select ordinary Dixon construction.

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

#### Step 4 determinant options

Checked Schur compression is enabled by default for eligible single-parameter
Dixon matrices over prime and extension fields, including recursive Dixon.
Step 3 uses complement-first (Schur-first) selection in both prime and
extension fields; `--no-schur-seed` disables this selection while retaining
Step 4 compression. Extension completion samples the extension generator and
the generator plus one, retaining full field coefficients. GF(16)/GF(256)
selection uses compact byte matrices, SIMD row operations and parallel matrix
products; other extensions retain the general FLINT matrix backend.
If the predicted complement fails, exact degree-layer pivoting selects a
smaller usable complement. GF(2^n) uses the native arithmetic kernels; other
extensions use the general field interface. A failed degree or invertibility
certificate falls back to the original determinant backend. Disable with
`--no-mq-step4-schur`; run `make test-extension-schur` for regression checks.

```bash
./drsolve --fq-det-method hnf <args>
./drsolve --fq-det-method iter <args>
./drsolve --fq-det-method interp <args>
```

Example:
```bash
./drsolve -r '[12]*3' 65537 --fq-det-method interp --threads 16 --time
```
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

## Paper

This library implements the methods described in:

> **Efficient Polynomial System Solving via Dixon Resultants: Applications to AO Primitives**
- ePrint: <https://eprint.iacr.org/2026/1969>
- reproduction materials: <https://github.com/drsolve/dixon_asiacrypt>

## Development Notes

The standalone general-degree quotient experiment accepts equal or mixed
positive equation degrees over a prime field. It generates dense square
systems with one variable retained and compares both Bezoutian construction
paths against the production Dixon resultant, coefficient by coefficient.

```bash
make experiment-degree-quotient
# Choose cases: number of equations followed by one degree or a degree list.
python3 src/test/mq/degree_quotient_experiment.py --cases 3:4 3:5 3:2,3,4 --repeat 3
# One generated system: N, degrees, seed, output prefix, prime, memory MiB.
./build/degree_quotient_experiment 3 2,3,4 1 /tmp/quotient_mixed 65537 512
```

Results and per-run timings are written to `out/degree_quotient_experiment/`.
The experiment uses dynamically indexed monomials, subject to address-space
checks and its optional estimated storage limit (512 MiB by default). It requires generic Macaulay ranks and a
fixed quotient basis across interpolation points; unsupported cases are
reported as errors. This remains an experimental executable, separate from
the optimized MQ backend selected by `drsolve --quotient`.

Parts of this project were developed with the assistance of AI-based coding tools. All AI-assisted contributions were reviewed and tested by the project author.

## License
DRSolve is distributed under the GNU General Public License version 2.0 (GPL-2.0-or-later). See the file COPYING.

# Usage Guide

[Documentation index](README.md)

`-h` / `--help` shows the quick reference. `--help-full` lists all options by topic, with examples. This guide explains selection rules and important constraints.

## Solving systems

Omitting elimination variables selects solving; `-s` requests it explicitly. `--solver auto|dixon|quotient` both enables solving and selects a backend.

| Backend | Use |
| --- | --- |
| `auto` (default) | Selects a path based on the input |
| `dixon` | Dixon elimination and back substitution |
| `quotient` | Finite-field systems with at least as many equations as variables; eligible square MQ systems use Bezoutian construction, while other supported systems use quotient closure |

`--quotient-max-degree D` and `--quotient-memory MiB` limit quotient closure degree and memory. Both default to 0, meaning no imposed budget. Hitting a budget makes auto fall back to Dixon; explicit quotient mode reports an incomplete solve, not absence of solutions. Nonzero budgets or `--no-quotient` retain closure for square systems. These budgets do not limit the separate quotient resultant construction.

For underdetermined systems, `--specialize-vars` fixes selected variables. This solves a chosen slice, not a complete parametric solution set.

```bash
./drsolve --solver quotient "x^2-1, y-x" 257
./drsolve --solver auto --quotient-memory 256 "x^2-1, y-x, x*y-1" 257
./drsolve -s --specialize-vars "z=0,w=1" "x+z, y+w" 257
```

## Elimination and resultants

Dixon is the default. `--macaulay` selects Macaulay; `--subres` handles two polynomials with one elimination variable. The equivalent selector is `--resultant dixon|macaulay|subres`. Automatic dispatch may use a pair-specific fast path depending on the field and input structure; explicit method options affect dispatch.

`--bezout-fast` / `--no-bezout-fast` controls fast finite-field pair resultants, enabled automatically for eligible inputs. Compatibility aliases are `--bivar-ntt` / `--no-bivar-ntt`. NTT or multipoint interpolation is used only when its conditions hold; otherwise computation falls back to a general path.

`--quotient` (alias `--mq-quotient`) enables quotient/Bezoutian resultant construction; it is **different from `--solver quotient`**. Supported square MQ systems with at least three variables use it automatically. Explicit selection also supports general-degree systems with one parameter. Failed applicability, rank, or degree checks fall back to Dixon.

`--resultant-only` skips root finding. It applies to elimination, not system solving.

```bash
./drsolve --subres "x^2+y, x+y+1" "x" 257
./drsolve --resultant-only "x+y, x-y^2+2" "x" 0
./drsolve --quotient -r "[2]*3" 257 --seed 1234
```

## Methods and advanced optimizations

`--method N` selects construction/determinant methods; `--step1 N` and `--step4 N` override individual stages. Not every method is supported over every field or at every stage.

| N | Method | Notes |
| --- | --- | --- |
| 0 | Minor expansion | Layered dynamic programming and caching |
| 1 | HNF | Polynomial-matrix determinants; multivariate paths may use Kronecker substitution |
| 2 | Interpolation | Evaluation and reconstruction using degree bounds |
| 3 | Sparse interpolation | Requires a suitable backend and sparse result structure |
| 4 | Bareiss | Fraction-free elimination |
| 5 | Recursive Dixon | Changes Dixon construction, rather than naming an ordinary Step 4 determinant backend |
| 6 | Balanced Laplace expansion | Experimental |

Finite-field elimination of two variables from three equations in three variables defaults to method 5. `--fq-det-method auto|hnf|iter|interp` selects the univariate polynomial-matrix backend; interp is prime-field only.

Eligible inputs use projected MQ indices, shared indices, coefficient filtering, and checked Schur compression by default. For comparisons, disable them with `--no-mq-step1-rank`, `--no-mq-step1-shared`, `--no-mq-step1-filter`, `--no-mq-step4-schur`, or `--no-schur-seed`. Schur compression supports single-parameter prime/extension-field systems of general degree, not just MQ.

`--mq-step1-pencil`, `--mq-step1-simplex`, `--fast-ksy`, and `--rank-pred` are disabled by default and are experimental or require additional conditions. `--step3-verify-second` adds a second verification pass. Consult complete help before changing advanced defaults.

```bash
./drsolve --step1 0 --step4 1 "x+y+z, x*y+y*z+z*x, x*y*z+1" "x,y" 257
./drsolve --no-quotient --no-mq-step4-schur -r "[2]*3" 257 --seed 1234
```

## Roots and reconstruction

Over Q, `--solve-rational-only` retains exact rational solutions. No rational solution does not imply no real solution. For roots after elimination, `--approx-roots` enables approximate real-root finding. `--root-precision BITS` also enables it; the default is 128 bits, with a range of 32 to 1048576.

`--complex` supports the Q-only two-equation, two-variable solver and complex roots of univariate rational resultants; it is not a general multivariate complex solver.

`--rational-root-scan auto|off|force` controls rational-root scanning before approximate real-root finding. Auto skips scans above 100 candidates; force still has a hard limit of 1000. These limits concern scan candidates, not the number of polynomial roots.

`--max-primes N` sets the reconstruction prime budget, defaulting to 64 for rational reconstruction. Large prime fields use native large-integer modular arithmetic by default. `--large-prime-reconstruct` selects reconstruction instead and cannot be combined with `--subres`.

```bash
./drsolve --solve-rational-only "x^2-1, y-x" 0
./drsolve --approx-roots --root-precision 256 "x+y, x-y^2+2" "x" 0
./drsolve --complex "x^2+y^2-1, x-y" 0
```

## Reduction and extraneous factors

`--ideal "a^3=b+1,..."` reduces using substitution relations. `--field-equation` applies `x^q=x` during multiplication; `--field-equation-s` applies it only to the final output. These change the quotient structure of the computation; they are not semantically neutral speed switches or a general Gröbner-basis interface. Q and large prime fields do not support these reduction options.

`--sift` uses quotient algebra to remove extraneous elimination factors and is disabled by default. `--sift-max-degree` defaults to 4 and must be positive. Sift applies to ordinary elimination over finite fields with machine-word characteristic and retained parameters. It is not combined with solving, complexity, ideal, or field-equation modes, and does not guarantee a generator of the elimination ideal for arbitrary inputs.

```bash
./drsolve --sift "x^2+y, x+y+1" "x" 257
./drsolve --ideal "a^3=b+1" "a^3+b, a+b^2" "a" 257
./drsolve --field-equation "x*y+z, y*z+x, x*z+y" "x,y" 2
```

## Random input and complexity

`-r` accepts total-degree specifications such as `"[2,3,4]"`, `"[2]*3"`, and `"[2]*4+[3]*2"`. `-n` sets the variable count. With `--vardeg`, degrees become per-variable exponent bounds and `-m` sets the equation count. `--homogeneous` keeps terms of the specified total degree and is incompatible with vardeg.

`--density` ranges from 0 to 1, defaulting to 0.5 over F_2 and 1 otherwise. Overdetermined random systems normally contain a planted solution; `--pure-random` disables planting. Random n+1 equations in n variables default to scalar consistency screening; add `-s` for solving.

`--seed` fixes random input; save the actual polynomials when comparing versions or backends. `-c` / `--comp` estimates Dixon matrix sizes, degree bounds, and operation counts. `--omega` defaults to 2.81. Estimates are not runtime or memory guarantees, and this complexity mode does not support large prime fields.

```bash
./drsolve -r -s "[2]*2" 257 --seed 1234
./drsolve -r --vardeg -m 3 "[1,1]" 257 --seed 1234
./drsolve -c -r --omega 2.81 "[2]*3" 257
```

## Parallelism, memory, and diagnostics

`--threads N` sets a positive thread count for supported parallel stages in suitable builds. `--cache N` limits live minor-DP entries for method 0, defaulting to 1024; 0 disables DP. `--array-limit-k K` caps extension-field multiplication tables at 2^K entries, for K from 0 to 62. This is an entry count, not a byte count.

`-v 0..4` ranges from silent to matrix output, defaulting to 1. Level 2 adds details and timing, level 3 profiling, and level 4 bounded matrix dumps and recursive traces. `--time` reports stage timings; some stages distinguish CPU time, wall time, and thread counts. Interpret memory statistics according to their reported scope and platform.

```bash
./drsolve --threads 2 --time -r "[2]*3" 257 --seed 1234
./drsolve -v 2 "x^2-1, y-x" 257
```

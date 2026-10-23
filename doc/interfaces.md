# Interfaces and Troubleshooting

[Documentation index](README.md)

## SageMath

The wrapper is [drsolve_sage_interface.sage](../drsolve_sage_interface.sage). Build the CLI first, then load it in SageMath:

```python
load("/path/to/drsolve/drsolve_sage_interface.sage")
set_drsolve_path("/path/to/drsolve/drsolve")
R.<x,y> = PolynomialRing(GF(257))
F = [x^2-1, y-x]
solutions = DixonSolve(F)
resultant = DixonRes(F, [x])
```

| Function | Purpose |
| --- | --- |
| `DixonRes(F, elim_vars, ...)` | Eliminate variables and return the parsed resultant string |
| `DixonSolve(F, ...)` | Solve and parse the output |
| `DixonComplexity(F, elim_vars, ...)` | Analyze complexity |
| `DixonIdeal(F, ideal_gens, elim_vars, ...)` | Eliminate with substitution relations |
| `set_drsolve_path(path)` | Set the executable path |
| `set_drsolve_output_dir(path)` | Set the wrapper's default output directory |

Common arguments include `field_size`, `drsolve_path`, `finput`, `foutput`, `output_dir`, `timeout`, `verbosity`, `time`, and `threads`; consult each function's signature. The field is normally inferred from Sage polynomials; specify it explicitly for plain strings. Check generator representations when reusing extension-field results.

The return value of `DixonSolve` depends on the result status. Callers must handle no solutions, incomplete computation, or parsing failure rather than assuming a fixed-size solution table. The wrapper launches the CLI and parses output; it is not an in-process C binding.

The Sage interface does not expose every advanced CLI flag. Use the CLI directly for full control of newer quotient/Bezout options. Undeclared option names cannot be passed as arbitrary Python keyword arguments.

## C integration and source layout

The build can produce `libdrsolve`, but headers under `src/` are internal and do not promise a stable ABI. Pin the source version and match FLINT and build settings when integrating directly. Use the CLI or Sage wrapper when a process-level boundary is preferable.

| Path | Contents |
| --- | --- |
| [apps/](../apps) | CLI entry point and argument parsing |
| [src/dixon/](../src/dixon) | Dixon construction and complexity analysis |
| [src/determinant/](../src/determinant) | Determinants, HNF, Schur, quotient construction |
| [src/field/](../src/field) | Prime/extension fields and unified arithmetic |
| [src/solver/](../src/solver) | Roots, back substitution, quotient solving |
| [src/test/](../src/test) | Regression tests |
| [examples/](../examples) | Sample inputs |

## Common questions

**Why is an expected solution missing?**

Check the coefficient field, extension modulus, and solving mode. No solution over F_p does not imply none over an extension; rational-only mode excludes irrational real solutions. Check specialization, budgets, and random consistency-screening behavior. Incomplete, degenerate, and unsatisfiable are different outcomes.

**Why does elimination return extra roots or zero?**

Elimination supplies necessary conditions and can introduce extraneous factors. Zero may reflect input or matrix degeneracy. Use solving with back substitution and verification. Where applicable, try `--sift`, without treating factor filtering as a general completeness proof.

**Why do more threads provide little speedup?**

Confirm OpenMP support in the build. Parsing, some construction stages, and small computations may be serial or dominated by scheduling overhead. Compare wall time for identical input with `--time`, not accumulated CPU time.

**When should precision or reconstruction budgets increase?**

Nearly repeated roots and large coefficients can make numerical root finding or exact reconstruction harder. `--root-precision` controls the applicable approximate-root path; `--max-primes` controls reconstruction. They serve different purposes.

**What if FLINT or another dynamic library cannot be found?**

Check runtime library paths and ensure versions and architectures match the build. Inspect dependency installation and linking mode; preserve the DLL layout in Windows release packages.

**How should an issue be reported?**

Include `--version`, OS, compiler, FLINT version, the exact command, a minimal `.dr` input, and relevant `-v 2 --time` output. For random input, supply both the seed and actual polynomials. Distinguish incorrect output, incomplete computation, and performance problems.

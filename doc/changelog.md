# Changelog

[Documentation index](README.md) · [Usage guide](usage.md) · [Algorithms](algorithms.md)

## Contents

- [0.6.0](#v060)
- [0.5.x](#v05)
- [0.4.x](#v04)
- [0.3.x](#v03)
- [0.2.x](#v02)
- [0.1.x and 0.0.x](#early)

<a id="v060"></a>
## 0.6.0

- Optimize square MQ solving and quotient/Bezoutian construction, including binary-extension arithmetic.
- Improve Schur compression and extend it to single-parameter systems of general degree.
- Add HNF determinant computation over extension fields.
- Add native large-integer modular arithmetic for large prime fields.
- Optimize extension-field parsing, root finding, and bivariate resultants using NTT or multipoint interpolation where applicable.
- Add memory reporting and improve macOS/Windows build and test compatibility.
- Add `--help-full` and a structured manual under `doc/`, keeping `-h` / `--help` concise.

<a id="v05"></a>
## 0.5.x

| Version | Date | Main changes |
| --- | --- | --- |
| 0.5.5 | 2026-10-12 | Add `--sift`; accelerate overdetermined back substitution with quotient algebra; optimize binary-extension computations |
| 0.5.4 | 2026-10-08 | Improve extension-field MQ, especially binary extensions; support efficient overdetermined quotient solving |
| 0.5.3 | 2026-09-30 | Checked MQ Schur compression; Step 1 degree filtering; recursive construction by default for three-variable systems |
| 0.5.2 | 2026-09-22 | Optimize MQ Step 3 and prime-field square-system Step 4; add `--vardeg` |
| 0.5.1 | 2026-09-07 | Layered Step 1 dynamic programming, Step 2 monomial-index reuse, and fixes |
| 0.5.0 | 2026-08-15 | Sparse prime-field determinants and rank-prediction optimizations; integrate Windows GUI and examples; reorganize sources |

<a id="v04"></a>
## 0.4.x

| Version | Date | Main changes |
| --- | --- | --- |
| 0.4.5 | 2026-08-05 | Add `-v 4` and homogeneous random input; improve rational reconstruction and timing |
| 0.4.4 | 2026-07-27 | Rank heuristics, F_2 random density, and real-root fixes |
| 0.4.3 | 2026-07-22 | Random input over arbitrary-size primes; reconstruction prime budgets; parallel CRT fixes |
| 0.4.2 | 2026-07-17 | Parallel CRT over Q and large prime fields; FLINT 3.6 root support; extension multiplication-table limits |
| 0.4.1 | 2026-07-13 | Small-field rank and timing fixes; split CLI implementation |
| 0.4.0 | 2026-07-10 | Prime-field HNF and submatrix-rank improvements; revise three-variable defaults |

<a id="v03"></a>
## 0.3.x

| Version | Date | Main changes |
| --- | --- | --- |
| 0.3.5 | 2026-07-07 | Improve Step 3 memory use and rank heuristics; refine small-resultant output |
| 0.3.4 | 2026-07-03 | Extension-field, complexity-analysis, and memory-management fixes |
| 0.3.3 | 2026-06-18 | Add variable specialization; improve small binary fields, parsing, zero handling, and ideal reduction |
| 0.3.2 | 2026-06-11 | Parallel minor reuse; complex solving for two equations in two variables over Q; Windows GUI fixes |
| 0.3.1 | 2026-06-05 | M-basis residual-update optimization in HNF |
| 0.3.0 | 2026-05-19 | ASIACRYPT submission stage; stored-minor algorithms, arbitrary-size prime solving/random input, test and output improvements |

<a id="v02"></a>
## 0.2.x

| Version | Date | Main changes |
| --- | --- | --- |
| 0.2.5 | 2026-05-14 | Parallel Dixon matrix construction and Step 3; improve Sage integration |
| 0.2.4 | 2026-05-10 | Add Zhao–Fu recursive construction and Qin et al. complexity analysis; expand help |
| 0.2.3 | 2026-05-07 | Improve complexity analysis, random input, detailed output, and default output directory |
| 0.2.2 | 2026-05-01 | Add verbosity/output options and bundled PML subset; legacy input compatibility; rational/sparse fixes |
| 0.2.1 | 2026-04-29 | Large-prime roots, Macaulay resultants, and Step 1 complexity analysis |
| 0.2.0 | 2026-04-25 | Rename DixonRes to DRSolve; expand prime-field determinants, rational/real solving, and file formats |

<a id="early"></a>
## 0.1.x and 0.0.x

| Version | Main changes |
| --- | --- |
| 0.1.5 | Add derivative-based sparse interpolation |
| 0.1.4 | Add fast HNF determinants and PML |
| 0.1.3 | Extension-field evaluation for small fields, method/thread options, automatic solving entry point |
| 0.1.2 | Rational/real/complex root support and `.dr` files |
| 0.1.1 | Sage interface, macOS, and older FLINT compatibility |
| 0.1.0 | Early integration of solving and elimination; see the corresponding release tag |
| 0.0.5 | Experimental rational CRT |
| 0.0.4 | Field-equation reduction and `make check` |
| 0.0.3 | Random input, complexity analysis, and multivariate binary-extension support |
| 0.0.2 | Memory, build, and license fixes |
| 0.0.1 | Initial release, corresponding to the CRYPTO submission stage |

See the [release tags](https://github.com/drsolve/drsolve/tags) and website for the full historical records.

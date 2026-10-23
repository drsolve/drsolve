# DRSolve 0.6.0 Documentation

DRSolve is a C tool for polynomial elimination and system solving over the rationals, prime fields, and extension fields, built on FLINT and PML. This manual covers version 0.6.0, its main algorithms, and their practical scope.

## Contents

| Document | Topics |
| --- | --- |
| [Installation](installation.md) | Dependencies, Makefile builds, platforms, deployment |
| [Input and coefficient fields](input.md) | CLI syntax, `.dr` files, field representations, output |
| [Usage guide](usage.md) | Solving, elimination, random input, complexity, options |
| [Algorithms](algorithms.md) | Dixon stages, MQ optimizations, Schur compression, quotient algebras |
| [Interfaces and troubleshooting](interfaces.md) | SageMath, C integration, common questions |
| [References](references.md) | Papers, method attribution, dependencies |
| [Changelog](changelog.md) | Version 0.6.0 and earlier releases |

## Quick start

After building in the repository root:

```bash
# Solve over F_257
./drsolve "x^2-1, y-x" 257
# Eliminate x, retaining y
./drsolve "x^2+y, x+y+1" "x" 257
# Quick reference and complete help
./drsolve -h
./drsolve --help-full
```

Omitting elimination variables selects solving; supplying them selects elimination. Automatic method selection is suitable for most inputs. Results go to `out/` by default; use `-o FILE` to choose a path.

Roots of an elimination polynomial are candidates, not necessarily solutions of the original system. Solving also requires back substitution and verification. Finite-field solving searches the specified field, not all its algebraic extensions. Rational mode distinguishes exact rational solutions from approximate real solutions; see the [usage guide](usage.md).

[Website](https://drsolve.github.io/) · [Source](https://github.com/drsolve/drsolve) · [Research paper](https://eprint.iacr.org/2026/1969). Licensed under [GPL-2.0-or-later](../COPYING).

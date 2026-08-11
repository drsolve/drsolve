# DRSolve examples

These files use the .dr file format documented in the main README.
For solver-mode files, the first line is the field size: 0 means the
rationals and a positive value selects a finite field.

Examples:

- cyclic5_rational.dr and cyclic5_finite.dr: classical cyclic-5 benchmarks.
- eco6_rational.dr: a moderate ECO benchmark over the rationals.
- eco10_finite.dr: a larger ten-variable ECO stress case.
- kat6_finite.dr: a Katsura-style structured quadratic system.
- mq5_f2.dr: a small multivariate-quadratic system over F_2.
- dense3_degree4.dr: a dense three-variable, degree-four system over F_257.

Run a solver example with:

    ./drsolve examples/cyclic5_rational.dr

# References

[Documentation index](README.md) · [Algorithms](algorithms.md)

These references cover the project paper and methods explicitly cited in source comments or project descriptions. Attribution does not mean every optimization comes directly from the cited paper or that its authors endorse this implementation.

<a id="r1"></a>
## R1 · DRSolve and Dixon solving

Haohai Suo, Jiamin Cui. **Efficient Polynomial System Solving via Dixon Resultants: Applications to AO Primitives**. Cryptology ePrint Archive, 2026/1969.

[Paper](https://eprint.iacr.org/2026/1969) · [Reproduction repository](https://github.com/drsolve/dixon_asiacrypt)

The overall method and cryptanalytic applications. Later software changes are documented in the changelog.

<a id="r2"></a>
## R2 · Recursive Dixon construction

Shizhong Zhao, Hongguang Fu. **An extended fast algorithm for constructing the Dixon resultant matrix**. Science in China Series A: Mathematics, 48, 131–143, 2005.

[Publisher / DOI](https://doi.org/10.1360/04ys0166)

Algorithmic background for recursive Dixon construction.

<a id="r3"></a>
## R3 · Dixon construction complexity

Xiaolin Qin et al. **Complexity of constructing Dixon resultant matrix**. International Journal of Computer Mathematics, 2017.

[Publisher / DOI](https://doi.org/10.1080/00207160.2016.1276572)

Matrix-construction complexity analysis; theoretical estimates do not directly determine runtime.

<a id="r4"></a>
## R4 · Polynomial determinants and minor reuse

W. M. Gentleman, S. C. Johnson. **Analysis of Algorithms, A Case Study: Determinants of Matrices with Polynomial Entries**. ACM Transactions on Mathematical Software, 1976.

[Publisher / DOI](https://doi.org/10.1145/355694.355696)

Background for determinant computation using stored minors.

<a id="r5"></a>
## R5 · HNF and polynomial-matrix determinants

George Labahn, Vincent Neiger, Wei Zhou. **Fast, deterministic computation of the Hermite normal form and determinant of a polynomial matrix**. Journal of Complexity, 2017; arXiv:1607.04176.

[Author preprint](https://arxiv.org/abs/1607.04176)

Fast HNF/determinant algorithms, implemented through PML and field-specific adaptations.

<a id="r6"></a>
## R6 · Approximant bases

Pascal Giorgi, Claude-Pierre Jeannerod, Gilles Villard. **On the Complexity of Polynomial Matrix Computations**. ISSAC, 2003.

[Author-hosted paper](https://www.lirmm.fr/~giorgi/issac03.pdf)

Approximant-basis algorithms underlying PML and related polynomial-matrix operations.

<a id="r7"></a>
## R7 · Kernel bases

Wei Zhou, George Labahn, Arne Storjohann. **Computing minimal nullspace bases**. ISSAC, 2012.

[Publisher / DOI](https://doi.org/10.1145/2442829.2442881)

Algorithmic source for PML kernel-basis computation.

<a id="r8"></a>
## R8 · Sparse interpolation

Qiao-Long Huang. **Sparse polynomial interpolation based on derivatives**. Journal of Symbolic Computation, 114, 359–375, 2023.

[Publisher / DOI](https://doi.org/10.1016/j.jsc.2022.06.002)

The derivative-based sparse interpolation method referenced by the implementation.

## Software dependencies

- [FLINT](https://flintlib.org/) / [source](https://github.com/flintlib/flint): integer, rational, finite-field, polynomial, and numerical ball arithmetic.
- [PML](https://github.com/vneiger/pml): polynomial-matrix algorithms. The repository includes a required subset with its copyright and license notices.
- [GMP](https://gmplib.org/) and [MPFR](https://www.mpfr.org/): underlying multiprecision arithmetic used by FLINT.

For research using DRSolve, cite R1 and record the software version, coefficient field, input, and important options for reproducibility.

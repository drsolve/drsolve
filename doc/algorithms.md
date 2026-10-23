# Algorithms

[Documentation index](README.md)

This chapter outlines the computation rather than reproducing full proofs. The overall method and cryptanalytic applications are described in the [Suo–Cui paper](references.md#r1). Later implementation improvements are recorded in the [changelog](changelog.md); the paper should not be assumed to cover every version 0.6.0 path.

## Four stages of Dixon elimination

Let K be the coefficient field. To eliminate x₁,…,xₙ from n+1 polynomials, treat the remaining variables z as parameters. Introduce y₁,…,yₙ and form a matrix by successively replacing x variables with y variables. Its determinant contains the factor ∏ᵢ(xᵢ−yᵢ). After division by this factor, expansion in x and y monomials yields a coefficient matrix with entries polynomial in z.

| Stage | Computation | Main costs |
| --- | --- | --- |
| Step 1 | Construct the Dixon polynomial | Polynomial determinants, products, intermediate terms |
| Step 2 | Extract coefficients into a matrix | Monomial indexing, storage, sparse data handling |
| Step 3 | Select a suitable full-rank square submatrix | Rank computation, pivots, minor selection |
| Step 4 | Compute its polynomial determinant | Polynomial-matrix arithmetic or evaluation/interpolation |

A nonzero elimination polynomial generally supplies a necessary condition on projected common zeros. Degenerate specializations, degree drops, and extraneous factors require care; candidates must be checked against the original equations. An identically zero result can indicate degeneracy and does not establish that every parameter value is a solution.

General construction uses minor expansion with shared intermediate results. Recursive construction exploits smaller subproblems and block structure; see [Zhao–Fu](references.md#r2). Complexity analysis uses matrix-size and degree-bound models related to [Qin et al.](references.md#r3). Actual costs also depend on sparsity, rank, and the field.

Source: [Dixon](../src/dixon/unified_dixon.c), [recursive construction](../src/dixon/dixon_recursive.c), [complexity analysis](../src/dixon/dixon_complexity.c).

## Minor reuse and MQ structure

MQ denotes multivariate quadratic systems. Degree and monomial constraints can exclude impossible coefficients before construction, reducing both intermediate work and matrix size. Projected indices, shared indices, and coefficient filtering avoid repeatedly generating the same monomials and minors.

Method 0 reuses minors by layers, with a cache controlling retained intermediate results. Stored-minor approaches to polynomial determinants are discussed by [Gentleman–Johnson](references.md#r4). Cache layouts, filtering policies, and specialized field kernels are implementation choices in DRSolve.

Compact MQ paths require applicability and rank checks. Different degrees, degenerate inputs, or unsupported coefficient representations may use general construction; being quadratic alone does not activate every optimization.

## Polynomial-matrix determinants

With one parameter, matrix entries lie in K[z]. HNF-based computation uses polynomial-matrix transformations and recursive decomposition, supported by approximant and kernel bases. The main determinant/HNF reference is [Labahn–Neiger–Zhou](references.md#r5); relevant basis algorithms include [Giorgi–Jeannerod–Villard](references.md#r6) and [Zhou–Labahn–Storjohann](references.md#r7).

Version 0.6.0 extends HNF capabilities to extension fields, with field-specific arithmetic backends. The bundled PML subset supplies part of the prime-field matrix machinery; not every field uses an identical PML routine.

Evaluation/interpolation computes numerical determinants at selected points and reconstructs the polynomial using degree bounds. Sparse interpolation exploits the support of the result; see [Huang](references.md#r8). Small fields may require extension-field evaluation, while multiple parameters introduce encoding and degree-bound considerations. Sparse input does not guarantee a sparse determinant.

Source: [unified HNF](../src/determinant/unified_hnf.c), [sparse interpolation](../src/determinant/fq_sparse_interpolation.c), [PML subset](../pml_det).

## Schur compression

After row and column permutations, write

```text
M = [ A  B ]
    [ C  D ]
```

For invertible A, `det(M) = det(A) · det(D − C A⁻¹ B)`. If a large, inexpensive constant block A is available, costly parameter-polynomial arithmetic can be confined to a smaller Schur complement.

DRSolve checks block selection and compression and accounts for determinant factors and permutations. Version 0.6.0 improves complement-first minor selection, field-specific kernels, and degenerate-case handling, with applicability extending to single-parameter systems of general degree. A general path remains necessary when no suitable block is found. Compression reduces computation; it does not replace solution verification.

Source: [Schur coefficient kernels](../src/determinant/mq_coeff_schur.h), [recursive extension-field Schur](../src/dixon/dixon_recursive_extension_schur.h).

## Quotient algebras, Bezoutians, and sift

For I=⟨f₁,…,fₘ⟩, a finite-dimensional quotient algebra A=K[x₁,…,xₙ]/I represents variable multiplication by finite matrices. Relations satisfied by these operators constrain common zeros. Root finding, back substitution, and checks against the original equations recover solutions.

The closure path incrementally generates monomials and relations and checks for a usable finite-dimensional representation. Having at least as many equations as variables does not guarantee zero-dimensionality or completion within a finite budget. Eligible square MQ systems can instead use Bezoutian structure to obtain compact representations.

Three related features serve different purposes:

- `--solver quotient`: solve using a quotient-algebra path.
- `--quotient`: construct a resultant using a quotient/Bezoutian fast path.
- `--sift`: remove extraneous factors from an elimination result using quotient relations.

Version 0.6.0 improves square MQ construction, extension-field arithmetic, and back substitution on top of earlier quotient functionality. This is a conceptual account: closure checks, shared layouts, and filtering details are defined by the implementation and are not collectively attributed to a single external paper.

Source: [quotient solver](../src/solver/quotient_solver.c), [MQ construction](../src/determinant/mq_quotient.c), [general-degree construction](../src/determinant/degree_quotient.c).

## Fast bivariate resultants

For two polynomials eliminating one variable, the remaining variable can serve as an interpolation parameter. Compute univariate resultants at selected points and recover the parameter polynomial. Version 0.6.0 selects Dixon/Bezout fast paths according to input structure, using NTT where modulus and transform conditions permit, or multipoint evaluation/interpolation.

Correct degree bounds, usable evaluation points, and nondegeneracy conditions are required. NTT is not universally available over finite fields; unsupported cases fall back to general computation. Explicit resultant or determinant selections affect automatic dispatch.

Source: [pair-resultant fast path](../src/resultant/bezout_fast.h).

## Field arithmetic and roots

Machine-word prime fields use modular integer arithmetic; larger primes use FLINT large-integer modular arithmetic. Eligible binary extensions use native binary representations and specialized multiplication. Other extensions use the corresponding unified interfaces. The generator definition must stay consistent through parsing, computation, and output.

Some rational elimination paths combine modular results with CRT and reconstruction. Prime budgets and coefficient heights affect completion. Large prime fields compute natively by default, with reconstruction available as an option.

Solving combines recursive elimination, root finding, back substitution, and verification. Finite-field roots are sought in the specified field. Rational mode distinguishes exact rational roots from approximate real roots, with complex roots available on selected paths. Read numerical results together with precision and residuals; decimal output is not an exact algebraic representation.

"""Degree counts for the arbitrary-n Dixon rank investigation.

This program never constructs a Dixon matrix or runs drsolve. The five ranks
for n=5 below are supplied by the user, not measurements made by this script.
No production code imports this file.

Notation: n degree-d equations, m=n-1 eliminated variables, one hidden
parameter; sigma=n*(d-1)+1. Work with the top homogeneous forms f_i in
S=k[x_1,...,x_m]. The characteristic-zero generic Hilbert function is

    H_p = [(1-z^d)^n/(1-z)^(n-1)]_p^+.

The Hilbert-series assertion is known for n forms in n-1 variables. It does
not identify the ranks of the individual Dixon maps. Special coefficients
or positive characteristic need their own justification.

EXACT FACTORIZATION (valid over any field):
After divided-difference row operations the cancellation matrix has first
row (f_1(x),...,f_n(x)). Let K_i be its signed first-row cofactors. Then

    Delta(x,y) = sum_i f_i(x)*K_i(x,y),   deg(K_i)=m*(d-1).

For p+q=sigma, let C_p: S_q^* -> S_(p-d)^n send a y-coefficient functional
to the corresponding coefficient vector of (K_1,...,K_n). Let

    M_p: S_(p-d)^n -> S_p,  (u_i) -> sum_i f_i*u_i.

Using full monomial bases, with zero rows/columns retained, B_(p,q)=M_p*C_p.
Writing U_p=im(C_p) and Z_p=ker(M_p), linear algebra gives

    rank(B_(p,q)) = dim(U_p) - dim(U_p intersect Z_p).                 (1)

The original Hilbert model gives dimensions of ideals and syzygies, but not
U_p or its intersection with Z_p. These spaces depend on the SAME f_i;
assuming random independent relative position is not justified.

GENERIC SYZYGY COUNTS ON THE SMALL-DEGREE SIDE:
For 0<=p<=floor(sigma/2), the untruncated Hilbert coefficient is nonnegative
and agrees with H_p. Indeed, (1-z^d)^n/(1-z)^(n-1) is (1-z) times the
symmetric unimodal polynomial (1+z+...+z^(d-1))^n. Hence

    dim Z_p = n*dim S_(p-d) - dim S_p + H_p
            = sum_{j=2..floor(p/d)} (-1)^j C(n,j)*dim S_(p-j*d).       (2)

The j-th Koszul module can contribute in this range precisely when
    j*d <= floor(sigma/2), equivalently (n-2*j)*d >= n-1.
These are graded Koszul terms, not a count of independent minimal generators.

For n=4 the range ends at 2*d-2, so Z_p=0 throughout the small-degree half:
all boundary rank loss there is already present in C_p. For n=5 it ends
at floor((5*d-4)/2)<3*d, so dim Z_p=10*C(p-2*d+3,3), with C=0 when the
upper argument is too small. This first becomes nonzero at d=4.

Summing dim Z_min(p,q) over ordered boundary pairs gives the ambient count
implemented by syzygy_ambient_sum(). For n=5 it simplifies to

    5*d^2*(d^2-4)/96                    (d even),
    5*(d+3)*(d+1)*(d-1)*(d-3)/96        (d odd).

These counts are 0,0,10,20,60 for d=2,...,6. The observed TOTAL corrections
are 0,0,4,14,46. The ambient counts are NOT correction terms. In particular,
at n=5,d=5 the previously recorded boundary blocks (8,13) and (9,12) have
defects 1 and 2 despite Z_8=Z_9=0. Counting syzygies alone cannot suffice.
Nor does the naive n=4 Hilbert-tail substitution work: the ternary tail
[(1-z^5)^5/(1-z)^3]^+ starting at 12-5+1=8 is zero, but defect(9,12)=2.

RELATION TO THE FULL MATRIX:
With b0(p,q) the existing boundary-rank model, h=sum H_p, and T the Schur
remainder after selecting actual boundary pivots, the exact identity is

    rho(D) = rho_old - sum_(p+q=sigma)
                  [b0(p,q)-dim U_p+dim(U_p intersect Z_p)]
                    + (rank(T)-h).                                  (3)

The sum here uses every ordered block, not just the small-degree half.
The bracket in (3) is a block defect; its two pieces need not be separately
nonnegative. There is no generic-rank proof for U_p, its intersection, or
rank(T)=h in this file. In particular (3) is a structural identity, NOT a
closed correction formula in n,d. Five supplied total ranks cannot identify
all these blockwise quantities or justify extrapolation to d>=7.

Background for the known generic Hilbert function, not a Dixon rank theorem:
https://arxiv.org/html/2503.16155v4 (Introduction).
"""

from itertools import accumulate
from math import comb


# Provenance: values supplied by the user; do not relabel as new measurements.
SUPPLIED_FIVE_VARIABLE_RANKS = {2: 32, 3: 165, 4: 541, 5: 1315, 6: 2750}


def choose(a, b):
    return comb(a, b) if a >= b >= 0 else 0


def monomials(variables, degree):
    return choose(degree + variables - 1, variables - 1)


def hilbert(n, variables, d):
    """Formal positive truncation; not a generic Hilbert theorem for all inputs."""
    values = []
    for p in range(n * d + 1):
        value = sum((-1)**j * comb(n, j) * monomials(variables, p-j*d)
                    for j in range(min(n, p // d) + 1))
        if value <= 0:
            break
        values.append(value)
    return values


def rank_model(n, d):
    """The existing degree-only model, with no proposed correction applied."""
    support = [1]
    for k in range(1, n):
        prefix = list(accumulate(support))
        support = [prefix[min(p, len(prefix)-1)]
                   for p in range(k * (d-1) + 1)]
    h = hilbert(n, n-1, d)
    sigma = n * (d-1) + 1
    boundary = []
    for p, rp in enumerate(support):
        q = sigma - p
        if 0 <= q < len(support):
            small = min(p, q)
            hp = h[small] if small < len(h) else 0
            boundary.append((p, q, max(0, min(rp, support[q])-hp)))
    return sum(h) + sum(b for _, _, b in boundary), boundary


def small_side_syzygies(n, d, p):
    """Generic characteristic-zero dim ker M_p, only on the small-degree side."""
    if not 0 <= p <= (n*(d-1)+1)//2:
        raise ValueError('p must be in the small-degree half')
    return sum((-1)**j * comb(n, j) * monomials(n-1, p-j*d)
               for j in range(2, min(n, p // d) + 1))


def syzygy_ambient_sum(n, d):
    """Hockey-stick sum of small-side syzygy dimensions; NOT rank correction."""
    sigma = n * (d-1) + 1
    low = (sigma - 1) // 2
    paired = 2 * sum((-1)**j * comb(n, j) * choose(low-j*d+n-1, n-1)
                     for j in range(2, n+1))
    center = small_side_syzygies(n, d, sigma//2) if sigma % 2 == 0 else 0
    return paired + center


def main():
    # Check two independently written degree-count identities; no matrix ranks.
    for n in range(3, 13):
        for d in range(2, 13):
            sigma = n * (d-1) + 1
            h = hilbert(n, n-1, d)
            for p in range(sigma//2 + 1):
                hp = h[p] if p < len(h) else 0
                assert small_side_syzygies(n, d, p) == (
                    n*monomials(n-1, p-d) - monomials(n-1, p) + hp)
            direct = sum(small_side_syzygies(n, d, min(p, sigma-p))
                         for p in range(sigma+1))
            assert syzygy_ambient_sum(n, d) == direct

    print('n=5; actual ranks supplied by user; no Dixon computations')
    print('d  model  supplied  correction  ambient_syzygies (NOT correction)')
    for d, rank in SUPPLIED_FIVE_VARIABLE_RANKS.items():
        model, _ = rank_model(5, d)
        print(f'{d}  {model:5}  {rank:8}  {model-rank:10}  '
              f'{syzygy_ambient_sum(5, d):16}')
    print('n=5 small-side degrees with nonzero syzygies:')
    for d in (4, 5, 6):
        print(d, {p: small_side_syzygies(5, d, p)
                  for p in range(2*d, (5*d-4)//2 + 1)})


if __name__ == '__main__':
    main()

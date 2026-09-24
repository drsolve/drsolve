# MQ forward/backward support count

This experiment counts monomials and transitions only. It never computes a
minor coefficient, a determinant, or a numerical rank. The production DP and
its filtering policy are unchanged.

Reproduce from the repository root:

```sh
make build/mq_support_count
python3 src/test/mq_support_count.py --output src/test/data/mq_layout/support_backward.jsonl
```

Here n is the Step 1 determinant order; there are n-1 eliminated variables
and one parameter. The fixture includes every degree-at-most-two input
monomial with coefficient one over F65537. These coefficients are placeholders
for generic supports, not a generic numerical system: identical equations
would have zero determinant. The divided-difference row supports are extracted
before any determinant computation. Every column is checked to have the same
support as its row union. No observed coefficient cancellations are used.

Targets are the initial canonical candidate rectangle selected by the existing
production degree-profile/mirror selector (205 by 205 for n=7, 480 by 480 for
n=8). No rank verification or repair is performed. Results therefore describe
this candidate projection, not all possible repaired target sets or arbitrary
MQ inputs. The quadratic row is processed last.

Let A_k be the union of monomial exponent vectors in the row processed at
layer k. Let I be the existing downward closure of the x/y target rectangle;
parameter exponents are unrestricted except by reachability. Forward supports
are F_0={0}, F_k=(F_(k-1)+A_k) intersect I. Applying I even in the certified
safe early layers is equivalent to the production filter's skipped checks.
At the terminal layer, retain only exact target x/y exponents. Going backward,
retain alpha in F_(k-1) iff some shift in A_k reaches a retained state in F_k.

This is exact forward/backward reachability in the row-union graph, not a
loose total-degree or prefix-degree test. It deliberately ignores column
exclusivity and determinant cancellations. Consequently it measures the
strongest pruning obtainable from remaining row unions and these targets
alone, not an upper bound on every algebraic optimization.

| n | Layer k | Existing support | Backward-live support | Removed |
|---:|---:|---:|---:|---:|
| 7 | 1 | 9 | 9 | 0 |
| 7 | 2 | 53 | 53 | 0 |
| 7 | 3 | 260 | 260 | 0 |
| 7 | 4 | 1,156 | 1,156 | 0 |
| 7 | 5 | 4,725 | 4,725 | 0 |
| 7 | 6 | 17,747 | 17,739 | 8 |
| 7 | 7 | 71,645 | 71,425 | 220 |
| 8 | 1 | 10 | 10 | 0 |
| 8 | 2 | 64 | 64 | 0 |
| 8 | 3 | 336 | 336 | 0 |
| 8 | 4 | 1,581 | 1,581 | 0 |
| 8 | 5 | 6,805 | 6,805 | 0 |
| 8 | 6 | 27,167 | 27,167 | 0 |
| 8 | 7 | 98,954 | 98,901 | 53 |
| 8 | 8 | 381,766 | 380,600 | 1,166 |

Count each accepted (source monomial, row shift) edge, then weight layer k by
binomial(n,k)*k, including root weight n. Since the fixture has identical
support in every column, this counts scalar product contributions in the
shared DP model. It does not count modular reductions, allocation, index
construction, rejected-edge checks, root summation, or exploit zero minor
coefficients. Therefore this is a workload estimate, not a timing prediction.

| n | Existing weighted edges | Retained weighted edges | Reduction |
|---:|---:|---:|---:|
| 7 | 5,924,009 | 5,916,771 | 0.122181% |
| 8 | 51,378,384 | 51,311,648 | 0.129891% |

Even with free pruning, savings in this arithmetic model are only about
0.12–0.13%. A backward pass itself scans transitions and stores extra supports.
These results do not justify implementing this pruning in the production
backend for the tested generic candidate projections. Earlier, heavily
weighted layers do not shrink at all. Different, much smaller target sets
could behave differently; no asymptotic claim for all n follows from two sizes.

Validation: an independent Python tuple-set oracle at n=4 reconstructs the
forward supports and performs reverse exponent subtraction, comparing every
layer's support count, live count, accepted-edge count and live-edge count
against the C implementation. It passed. Both n=7/8 analyses completed with
assertions enabled. Raw model supports, targets, counts and summaries are
written to the ignored local data directory. No determinant timing or
coefficient-equality benchmark was run for this experiment.

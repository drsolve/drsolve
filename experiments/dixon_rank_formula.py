"""Experimental rank conjecture for four dense degree-d equations.

Notation: n equations, m=n-1 eliminated variables. No production code imports
this file. Numerical evidence here is not a proof of generic maximal rank.

The production model for n=4 simplifies to (5*d**3-2*d)/3. Let
    H2(a) = [(1-z**d)**4/(1-z)**2]_a^+.
For a>=d this is max(d-3*(a-d+1), 0). The candidate extra nullity of
boundary block (p,q), p+q=4*d-3, is
    delta(p,q) = sum(H2(a), a>=max(p,q)-d+1).
Summing the two halves yields
    Delta = 2*sum(i*(d-3*i), i=1..floor((d-1)/3))
          = k*(k+1)*(d-2*k-1), k=floor((d-1)/3).
Thus the candidate rank is (5*d**3-2*d)/3 - Delta.

PROOF GAPS (do not treat the calculation of the sum as closing them):
1. Identify the excess kernel of EACH Dixon boundary block with the proposed
   binary Hilbert-function tail, or give matching upper/lower rank bounds.
2. Prove that the remaining Schur block has rank sum(H_m), with H_m the
   Hilbert function of n generic degree-d forms in m variables.
These hold for the samples recorded here, but are not proved in general.
The special-family correction is not a formula for arbitrary coefficients;
sparse/special systems can have smaller rank at the same (n,d).
No formula for arbitrary n or mixed degrees is claimed. In particular the
naive extension replacing H2 by H_(n-2) is false: for n=5,d=5, block (9,12)
has model rank 135 and actual rank 133, while the proposed ternary tail
starts at 12-5+1=8 and is zero. This is a recorded counterexample.

An exact matrix identity, distinct from the conjecture: choose independent
rows/columns inside every degree-boundary block (p,sigma-p). Their combined
pivot matrix is triangular in degree blocks and invertible. If S is the
Schur remainder, rank(D)=sum(rank(B_p))+rank(S). The diagnostic records
rank(S) as rank(D)-sum(rank(B_p)); it does NOT assume it equals sum(H_m).

Generic rank evidence is subject to specialization: a modular rank provides
an attained rank, not a symbolic upper bound in characteristic zero.

Relevant background, not a proof of the Dixon conjecture:
https://arxiv.org/abs/2503.16155 (Hilbert series versus graded syzygies)
https://arxiv.org/abs/math/0007036 (Bezout maps and resultant complexes)
"""
from math import comb
from pathlib import Path
import argparse
import json
import subprocess


def four_equation_rank(degree):
    k = (degree - 1) // 3
    return (5 * degree**3 - 2 * degree) // 3 - k * (k + 1) * (degree - 2 * k - 1)


def binary_hilbert(degree, a):
    """Positive truncation, stopping at the first nonpositive coefficient."""
    for t in range(a + 1):
        coefficient = sum((-1)**j * comb(4, j) * (t-j*degree+1)
                          for j in range(min(4, t // degree) + 1))
        if coefficient <= 0:
            return 0
    return coefficient


def run_case(root, n, d, prime, seed):
    result = subprocess.run([str(root / 'build/dixon_rank_probe'), str(n), str(d),
                             str(prime), str(seed)], cwd=root, check=True,
                            capture_output=True, text=True, timeout=300)
    return json.loads(result.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--output', type=Path, default=Path('experiments/dixon_rank_samples.jsonl'))
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    subprocess.run(['make', '-f', 'experiments/Makefile', 'build/dixon_rank_probe'],
                   cwd=root, check=True)
    cases = [(4, d, 65537, 12345) for d in range(2, 13)]
    cases += [(4, d, p, s) for d in (4, 7, 10, 12)
              for p, s in ((65521, 77), (1000003, 99))]
    cases += [(3, d, 65537, 12345) for d in (2, 3, 4, 7)]
    cases += [(5, d, 65537, 12345) for d in (2, 3, 4, 5)]
    cases += [(6, d, 65537, 12345) for d in (2, 3)]
    discrepancies = []
    with args.output.open('w') as output:
        for n, d, prime, seed in cases:
            record = run_case(root, n, d, prime, seed)
            output.write(json.dumps(record, separators=(',', ':')) + '\n')
            output.flush()
            candidate = four_equation_rank(d) if n == 4 else None
            if n == 4:
                if record['rank'] != candidate:
                    discrepancies.append((n, d, prime, seed, 'total rank'))
                for block in record['blocks']:
                    start = max(block['p'], block['q']) - d + 1
                    defect = sum(binary_hilbert(d, a) for a in range(start, 2*d))
                    if block['rank'] != block['predicted'] - defect:
                        discrepancies.append((n, d, prime, seed, block['p'], block['q']))
            if record['residual_rank'] != record['h']:
                discrepancies.append((n, d, prime, seed, 'residual != H'))
            print(f'[{d}]*{n} p={prime} seed={seed}: '
                  f'model={record["predicted"]} actual={record["rank"]} '
                  f'candidate={candidate} residual={record["residual_rank"]} H={record["h"]}',
                  flush=True)
    print(f'{len(cases)} samples; discrepancies: {discrepancies}', flush=True)
    return bool(discrepancies)


if __name__ == '__main__':
    raise SystemExit(main())

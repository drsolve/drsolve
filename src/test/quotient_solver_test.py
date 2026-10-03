"""Sage/GB audit of the production C quotient and complete rational points."""
from sage.all import GF, PolynomialRing
from dixon_over_routes import solve_gb
from dixon_over_routes_bench import inputs
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--binary', type=Path, default=Path('build/quotient_solver_driver'))
    args = ap.parse_args()
    count = 0
    with tempfile.TemporaryDirectory(prefix='quotient-production-audit-') as tmp:
        source = Path(tmp) / 'case.txt'

        def check(fs, degree=6, expected=None, analytic=None):
            nonlocal count
            P = fs[0].parent()
            prime, n = P.base_ring().characteristic(), P.ngens()
            source.write_text(f'{prime} {n} {len(fs)}\n' + '\n'.join(map(str, fs)) + '\n')
            result = json.loads(subprocess.check_output(
                [str(args.binary.resolve()), str(source), str(degree), '512'], text=True))
            if analytic is None and P.ideal(fs).dimension() > 0:
                assert not result['certified'], result
                return
            assert result['certified'], (fs, result)
            reference = solve_gb(fs) if analytic is None else analytic
            assert result['dimension'] == reference['dimension'], (fs, result, reference)
            elim = sum(P.base_ring()(c)*P.gen(n-1)**i for i, c in enumerate(reference['eliminant']))
            assert P(result['eliminant']) == elim, (fs, result, elim)
            got = {tuple(a) for a in result['points']}
            assert len(got) == len(result['points']), result
            assert all(all(f(*a) == 0 for f in fs) for a in got)
            if expected is None:
                if prime**n <= 3000:
                    expected = {a for a in itertools.product(range(prime), repeat=n)
                                if all(f(*a) == 0 for f in fs)}
                else:
                    expected = {tuple(int(a[x]) for x in P.gens())
                                for a in P.ideal(fs).variety()}
            assert got == expected, (fs, got, expected)
            count += 1

        for prime in (2, 3, 7, 257):
            for n in (3, 4):
                for seed in (0, 1):
                    for family in ('dense', 'two'):
                        fs, _ = inputs(n, n+2, seed, family, prime)
                        check(fs)
            P = PolynomialRing(GF(prime), names=('x0', 'x1'))
            x, y = P.gens()
            for fs in ([x*x, y*y], [x*x, (y-x)**2], [x*x-x, y*y],
                       [x*x-x, y*y-y], [x*x, x*y, y*y],
                       [x*x, y*y, P(1)], [x*x, y*y, P(0), x*x],
                       [x**3-y, y**2-1], [x**2+1, y]):
                check(fs)
            check([x*x, P(0)], degree=4)
        P = PolynomialRing(GF(257), names=('x0', 'x1'))
        x, y = P.gens()
        check([x**7-1, y-x], degree=9)
        # More than eight variables and characteristic near the word limit.
        P = PolynomialRing(GF(18446744073709551557), names=tuple(f'x{i}' for i in range(9)))
        check([x-i for i, x in enumerate(P.gens())] + [P.gen(0)], degree=2,
              expected={tuple(range(9))},
              analytic={"dimension": 1, "eliminant": [-8, 1]})
    print(f'PASS: {count} exact quotient/minimal-polynomial/complete-point audits; positive-dimensional rejection')


if __name__ == '__main__':
    main()

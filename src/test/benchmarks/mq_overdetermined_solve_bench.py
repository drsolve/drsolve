#!/usr/bin/env python3
"""Paired planted-root benchmarks: identical square prefix plus extra quadrics.

This preserves a chosen root, not necessarily the complete algebraic zero set.
Every returned point is independently evaluated in every input polynomial.
"""
import argparse
import json
import os
from pathlib import Path
import random
import re
import subprocess
import time


def row_rank(rows, prime):
    a = [row[:] for row in rows]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [v * inv % prime for v in a[rank]]
        for i in range(rank + 1, len(a)):
            c = a[i][col]
            a[i] = [(v - c*w) % prime for v, w in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def generate(n, prime, seed):
    rng = random.Random(seed)
    root = [rng.randrange(1, prime) for _ in range(n)]
    monomials = [(i, j) for i in range(n) for j in range(i, n)] + [(i,) for i in range(n)] + [()]
    coefficients = []
    for _ in range(n + 2):
        row = [rng.randrange(prime) for _ in monomials[:-1]]
        value = sum(c * product(root, m) for c, m in zip(row, monomials))
        row.append(-value % prime)
        coefficients.append(row)
    assert row_rank(coefficients, prime) == n + 2
    return root, monomials, coefficients


def product(point, monomial):
    value = 1
    for i in monomial:
        value *= point[i]
    return value


def polynomial(row, monomials):
    return '+'.join(str(c) + ''.join('*x' + str(i) for i in m)
                    for c, m in zip(row, monomials) if c) or '0'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, nargs='+', default=[8])
    parser.add_argument('--seeds', type=int, nargs='+', default=[12345])
    parser.add_argument('--extra', type=int, nargs='+', choices=[0, 1, 2], default=[1, 0, 2])
    parser.add_argument('--methods', nargs='+', choices=['quotient', 'mq'], default=['quotient'])
    parser.add_argument('--threads', type=int, default=16)
    parser.add_argument('--timeout', type=int, default=45)
    parser.add_argument('--prime', type=int, default=65537)
    parser.add_argument('--redundant', action='store_true')
    parser.add_argument('--out', default='out/jouanolou/overdetermined_solve')
    args = parser.parse_args()
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, LD_LIBRARY_PATH='.:/usr/local/lib')
    records = []
    for n in args.n:
        for seed in args.seeds:
            root, monomials, rows = generate(n, args.prime, seed)
            for extra in args.extra:
                chosen = rows[:n + extra]
                if args.redundant:
                    assert extra == 1
                    chosen = rows[:n] + [[sum(row[j] for row in rows[:n]) % args.prime
                                         for j in range(len(monomials))]]
                for method in args.methods:
                    stem = f'n{n}_m{n+extra}_s{seed}_{method}' + ('_redundant' if args.redundant else '')
                    source, result = out / (stem + '.dr'), out / (stem + '_result.dr')
                    source.write_text(str(args.prime) + '\n' + ',\n'.join(polynomial(row, monomials) for row in chosen) + '\n')
                    (out / (stem + '_input.json')).write_text(json.dumps(dict(root=root, coefficients=chosen, monomials=monomials)))
                    rss = out / (stem + '.rss')
                    command = ['/usr/bin/time', '-f', '%M', '-o', str(rss), 'timeout', str(args.timeout),
                               './drsolve', '-s', '-f', str(source), '--threads', str(args.threads), '--time',
                               '--quotient-memory', '0', '-o', str(result)]
                    command += ['--solver', 'quotient'] if method == 'quotient' else ['--solver', 'auto', '--mq-quotient']
                    start = time.perf_counter()
                    p = subprocess.run(command, text=True, capture_output=True, env=env, timeout=args.timeout+10)
                    elapsed = time.perf_counter()-start
                    log = p.stdout+p.stderr; (out / (stem + '.log')).write_text(log)
                    record = dict(n=n, equations=len(chosen), seed=seed, method=method, redundant=args.redundant,
                                  status=p.returncode, process_seconds=elapsed, peak_rss_kib=int(rss.read_text().splitlines()[-1]))
                    if p.returncode == 0:
                        points = []
                        for line in log.splitlines():
                            if '[PASS]' not in line:
                                continue
                            coordinates = dict((int(i), int(v)) for i, v in re.findall(r'x(\d+) = (\d+)', line))
                            assert len(coordinates) == n, line
                            point = [coordinates[i] for i in range(n)]
                            assert all(sum(c*product(point, m) for c, m in zip(row, monomials)) % args.prime == 0 for row in chosen)
                            points.append(point)
                        assert root in points, (stem, root, points, log)
                        declared = re.search(r'Found (\d+) solution set', log)
                        assert declared and int(declared[1]) == len(points), log
                        record['points'] = points
                        summary = re.search(r'certified dimension (\d+), degree (\d+), (\d+) monomials', log)
                        if summary:
                            record.update(dimension=int(summary[1]), degree=int(summary[2]), monomials=int(summary[3]))
                        wall = re.search(r'Total -.*Wall time: ([0-9.]+)', log)
                        record['wall_seconds'] = float(wall[1])
                    elif p.returncode != 124:
                        raise RuntimeError(log)
                    records.append(record)
                    (out / 'results.json').write_text(json.dumps(records, indent=2))
                    print(json.dumps(record), flush=True)


if __name__ == '__main__':
    main()

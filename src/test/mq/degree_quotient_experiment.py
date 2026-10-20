#!/usr/bin/env python3
"""Reproducible general-degree quotient checks and single-thread timings.

Run `make experiment-degree-quotient` from the repository root.
The C prototype times normal forms, quotient recurrence, full Bezoutian block
construction followed by reduction, interpolation, and HNF separately.
The CLI comparison records its reported computation wall time (1 ms precision),
not process startup. Raw per-run records and median summaries are saved as JSON.
"""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess


def coefficients(text, variable, prime):
    expression = text.split('Resultant:\n', 1)[1].splitlines()[0].replace(' ', '')
    values = {}
    for term in expression.split('+'):
        if variable in term:
            coefficient, _, power = term.partition(variable)
            c = int(coefficient.rstrip('*') or '1')
            exponent = int(power[1:]) if power else 1
        else:
            c, exponent = int(term), 0
        values[exponent] = (values.get(exponent, 0) + c) % prime
    result = [values.get(i, 0) for i in range(max(values) + 1)]
    inverse = pow(result[-1], -1, prime)
    return [c * inverse % prime for c in result]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', default='build/degree_quotient_experiment')
    parser.add_argument('--drsolve', default='./drsolve')
    parser.add_argument('--output', default='out/degree_quotient_experiment')
    parser.add_argument('--repeat', type=int, default=3)
    parser.add_argument('--cases', nargs='+', default=[
        '3:2', '3:3', '2:4', '3:4', '3:5', '3:8', '4:4',
        '3:2,3,4', '3:1,3,5', '3:1,1,5', '4:2,3,4,5'],
        help='N:degree or N:d1,d2,...,dN')
    parser.add_argument('--seeds', nargs='+', type=int, default=[1, 12345, 54321])
    parser.add_argument('--prime', type=int, default=65537)
    parser.add_argument('--memory-mib', type=int, default=512)
    parser.add_argument('--timeout', type=float, default=60)
    parser.add_argument('--kernel', choices=['grouped', 'reference'], default='grouped')
    args = parser.parse_args()
    cases = []
    for spec in args.cases:
        n_text, degree_text = spec.split(':', 1)
        n = int(n_text)
        degrees = list(map(int, degree_text.split(',')))
        if len(degrees) == 1:
            degrees *= n
        if len(degrees) != n or n < 2 or min(degrees) < 1:
            parser.error('invalid case: ' + spec)
        cases.append((n, degrees))
    assert args.repeat > 0
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    binary, drsolve = str(Path(args.binary).resolve()), str(Path(args.drsolve).resolve())
    env = dict(os.environ, OMP_NUM_THREADS='1')
    # Rejections must happen before expensive allocations or polynomial work.
    invalid_prefix = str(output / 'rejected')
    for command, reason in [
        (['9', '2', '1', invalid_prefix], 'N must be between'),
        (['3', '2,3', '1', invalid_prefix], 'degree list'),
        (['3', '0', '1', invalid_prefix], 'positive integer'),
        (['3', '2', '1', invalid_prefix, '15'], 'must be prime'),
        (['3', '4', '1', invalid_prefix, '7'], 'prime too small'),
        (['8', '34', '1', invalid_prefix], 'packed exponent capacity'),
        (['4', '4', '1', invalid_prefix, '65537', '1'], 'exceeds budget'),
        (['3', '3', '1', invalid_prefix, '65537', '512', 'invalid'], 'kernel must'),
    ]:
        rejected = subprocess.run([binary, *command], capture_output=True,
                                  text=True, env=env, timeout=10)
        assert rejected.returncode == 2 and reason in rejected.stderr, rejected
    records = []
    for n, degrees in cases:
        degree_spec = ','.join(map(str, degrees))
        for seed in args.seeds:
            prefix = output / f'n{n}_d{degree_spec.replace(",", "-")}_s{seed}'
            for repeat in range(args.repeat):
                p = subprocess.run([binary, str(n), degree_spec, str(seed), str(prefix),
                                    str(args.prime), str(args.memory_mib), args.kernel],
                                   capture_output=True, text=True, env=env, timeout=args.timeout)
                if p.returncode:
                    raise RuntimeError(f'prototype failed: {prefix}\n{p.stdout}\n{p.stderr}')
                row = json.loads(p.stdout)
                reference = str(prefix) + '_reference.dr'
                p = subprocess.run([drsolve, '-f', str(prefix) + '.dr', '--no-quotient',
                                    '--resultant-only', '--threads', '1', '--time', '-v', '1',
                                    '-o', reference], capture_output=True, text=True, env=env, timeout=args.timeout)
                if p.returncode:
                    raise RuntimeError(f'Dixon failed: {prefix}\n{p.stdout}\n{p.stderr}')
                Path(str(prefix) + '_reference.log').write_text(p.stdout + p.stderr)
                expected = list(map(int, Path(str(prefix) + '.coeff').read_text().split()))
                actual = coefficients(Path(reference).read_text(), f'x{n-1}', args.prime)
                assert actual == expected, f'coefficient mismatch: {prefix}'
                row['repeat'] = repeat
                row['dixon_wall'] = float(re.search(r'Total -.*Wall time: ([0-9.]+)', p.stdout)[1])
                row['prototype_total'] = sum(row[k] for k in ('normal', 'reduced', 'interpolation', 'determinant'))
                row['full_block_total'] = sum(row[k] for k in ('normal', 'direct', 'interpolation', 'determinant'))
                records.append(row)
                (output / 'runs.json').write_text(json.dumps(records, indent=2) + '\n')
            print(f'PASS degrees={degrees}, seed={seed}: '
                  f'core={row["core"]}, resultant degree={row["resultant_degree"]}', flush=True)
    summaries = []
    for n, degrees in cases:
        selected = [r for r in records if r['degrees'] == degrees]
        summary = {k: selected[0][k] for k in ('n', 'degrees', 'core', 'full_rows', 'full_cols', 'resultant_degree')}
        summary['runs'] = len(selected)
        for key in ('normal', 'reduced', 'direct', 'interpolation', 'determinant',
                    'prototype_total', 'full_block_total', 'dixon_wall', 'operator_setup', 'operator_updates'):
            summary[key + '_median_ms'] = 1000 * statistics.median(r[key] for r in selected)
        summary['update_matrix_products_median'] = statistics.median(r['update_matrix_products'] for r in selected)
        summaries.append(summary)
        print(json.dumps(summary), flush=True)
    report = {
        'prime': args.prime, 'threads': 1, 'kernel': args.kernel,
        'notes': [
            'Prototype totals sum normal forms, core construction, interpolation and HNF; exclude input generation and validation.',
            'Direct means constructing the homogeneous Bezoutian middle block then reducing; it is not the production Dixon algorithm.',
            'Production Dixon times are CLI computation wall times rounded to 1 ms; process startup is excluded.',
            'The prototype has no shared leading frames, packed kernels or sample parallelism.',
            'Grouped kernel setup is included in the recurrence time; update_matrix_products counts update multiplications only, excluding operator construction.',
            'Every core is checked against full-block reduction, interpolation at two extra points, and the monic resultant against production Dixon.',
        ],
        'records': records, 'summary': summaries,
    }
    (output / 'results.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()

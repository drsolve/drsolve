#!/usr/bin/env python3
"""Compare complete resultants and single-thread Step 4 medians.

Example: python3 src/test/pml_prime_step4_bench.py --degrees 8 12 16 20
Legacy switches isolate these changes in the same binary and FLINT build.
"""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--degrees', type=int, nargs='+', default=[16])
parser.add_argument('--seeds', type=int, nargs='+', default=[12345])
parser.add_argument('--prime', type=int, default=65537)
parser.add_argument('--repeats', type=int, default=3)
parser.add_argument('--ablations', action='store_true')
parser.add_argument('--json', type=Path)
args = parser.parse_args()
if args.repeats < 1:
    parser.error('--repeats must be positive')
root = Path(__file__).resolve().parents[2]
variants = [('baseline', 'legacy', 'legacy'), ('optimized', 'auto', 'auto')]
if args.ablations:
    variants += [('multiplication', 'auto', 'legacy'), ('row_order', 'legacy', 'auto')]
records = []
with tempfile.TemporaryDirectory(prefix='drsolve-prime-step4-') as directory:
    for degree in args.degrees:
        for seed in args.seeds:
            reference = None
            for repeat in range(args.repeats):
                # Alternate order to reduce warm-up and frequency bias.
                for name, multiplication, row_order in variants[::1 if repeat % 2 == 0 else -1]:
                    output = Path(directory) / 'result.dr'
                    env = dict(os.environ, DRSOLVE_PML_MUL=multiplication,
                               DRSOLVE_PML_DET_ROW_ORDER=row_order,
                               DRSOLVE_PML_MIDDLE_PRODUCT='', DRSOLVE_PML_DET_ALGO='lnz')
                    command = [str(root / 'drsolve'), '-r', f'[{degree}]*3', str(args.prime),
                               '--seed', str(seed), '--threads', '1', '-v', '2', '-o', str(output)]
                    run = subprocess.run(command, cwd=root, env=env, check=True,
                                         capture_output=True, text=True)
                    result = output.read_text().split('Resultant:\n', 1)[1].split('\n', 1)[0]
                    if reference is None:
                        reference = result
                    if result != reference:
                        raise RuntimeError(f'Resultant mismatch: {degree=} {seed=} {name=}')
                    seconds = float(re.search(r'Step 4 time: ([0-9.]+)', run.stdout)[1])
                    records.append(dict(degree=degree, prime=args.prime, seed=seed,
                                        variant=name, repeat=repeat, step4=seconds))
                    print(f'd={degree} p={args.prime} seed={seed} {name} '
                          f'run={repeat+1}: {seconds:.3f}s; complete resultant matches', flush=True)
            for name, _, _ in variants:
                samples = [r['step4'] for r in records if r['degree'] == degree and
                           r['seed'] == seed and r['variant'] == name]
                print(f'MEDIAN d={degree} p={args.prime} seed={seed} {name}: '
                      f'{statistics.median(samples):.3f}s', flush=True)
if args.json:
    args.json.write_text(json.dumps(records, indent=2) + '\n')

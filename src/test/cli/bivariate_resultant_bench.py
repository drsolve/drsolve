"""Compare full monic resultants and timings on identical seeded inputs."""
import argparse
import json
import re
import statistics
import subprocess
import tempfile
import time
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', default='./drsolve')
parser.add_argument('--degree', type=int, default=80)
parser.add_argument('--threads', type=int, nargs='+', default=[1, 16])
parser.add_argument('--seeds', type=int, nargs='+', default=[1, 2, 3])
parser.add_argument('--timeout', type=float, default=180)
parser.add_argument('--output', default='/tmp/bivariate-resultant-bench.json')
args = parser.parse_args()
binary = str(Path(args.binary).resolve())
methods = {'dixon': ['--dixon', '--step4', '1'],
           'subres': ['--subres'], 'ntt': ['--bezout-fast']}
records = []
references = {}
with tempfile.TemporaryDirectory() as tmp:
    result_path = Path(tmp) / 'result.dr'
    for threads in args.threads:
        for seed in args.seeds:
            for name, flags in methods.items():
                cmd = [binary, '-r', f'[{args.degree}]*2', '65537', '--seed', str(seed),
                       '--threads', str(threads), '--resultant-only', '--time', '-v', '2',
                       '-o', str(result_path), *flags]
                start = time.perf_counter()
                proc = subprocess.run(cmd, capture_output=True, text=True,
                                      check=True, timeout=args.timeout)
                wall = time.perf_counter() - start
                result = re.sub(r'\s+', '', result_path.read_text().split('Resultant:\n', 1)[1])
                if seed in references:
                    assert result == references[seed], (threads, seed, name, 'coefficient mismatch')
                else:
                    references[seed] = result
                if name == 'ntt':
                    assert 'Dixon/Bezout fast time:' in proc.stdout and 'fallback' not in proc.stdout
                timing = [line for line in proc.stdout.splitlines()
                          if 'time:' in line.lower() or 'Evaluation ' in line or 'Memory:' in line]
                record = dict(method=name, threads=threads, seed=seed, wall=wall, timing=timing)
                records.append(record)
                Path(args.output).write_text(json.dumps(records, indent=2))
                print(f'{name:6s} threads={threads:2d} seed={seed} wall={wall:.6f}s coefficients=PASS', flush=True)
for threads in args.threads:
    for method in methods:
        times = [r['wall'] for r in records if r['threads'] == threads and r['method'] == method]
        print(f'{method:6s} threads={threads:2d} median={statistics.median(times):.6f}s')

#!/usr/bin/env python3
"""Interleaved compact/legacy Step 1+2 CLI comparison, with exact result checks."""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--repeats', type=int, default=3)
p.add_argument('--sizes', type=int, nargs='+', default=[7, 8])
p.add_argument('--threads', type=int, default=4)
a = p.parse_args()
records = []
with tempfile.TemporaryDirectory(prefix='mq-compact-') as temp:
    for n in a.sizes:
        expected = None
        for repeat in range(a.repeats):
            modes = [0, 1] if repeat % 2 == 0 else [1, 0]
            for compact in modes:
                env = dict(os.environ, DRSOLVE_MQ_COMPACT=str(compact), DRSOLVE_STEP2_PACK_BUFFER='0',
                           LD_LIBRARY_PATH=str(ROOT)+':'+os.environ.get('LD_LIBRARY_PATH', ''))
                output = Path(temp)/'result.dr'
                command = [str(ROOT/'drsolve'), '-r', f'[2]*{n}', '257', '--seed', '1790242205',
                           '--threads', str(a.threads), '--time', '-v', '2', '-o', str(output)]
                run = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, text=True,
                                     check=True, timeout=180)
                marker = 'Using compact Step 1 rows directly' if compact else 'Released source Dixon terms'
                assert marker in run.stdout, run.stdout
                result = '\n'.join(line for line in output.read_text().splitlines() if not line.startswith('Time: '))
                if expected is None: expected = result
                assert result == expected
                times = {}
                for step in (1, 2):
                    match = re.search(rf'Step {step} time: CPU time: ([\d.]+) seconds \| Wall time: ([\d.]+)', run.stdout)
                    assert match, run.stdout
                    times[f'step{step}'] = float(match[2])
                records.append(dict(n=n, compact=compact, repeat=repeat, threads=a.threads, **times))
                print(records[-1], flush=True)
                a.output.parent.mkdir(parents=True, exist_ok=True)
                a.output.write_text(json.dumps(records, indent=2)+'\n')
for n in a.sizes:
    for compact in (0, 1):
        subset = [r for r in records if r['n'] == n and r['compact'] == compact]
        print(n, compact, 'median Step 2:', statistics.median(r['step2'] for r in subset),
              'median Step 1+2:', statistics.median(r['step1']+r['step2'] for r in subset))

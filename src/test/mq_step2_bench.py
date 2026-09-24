#!/usr/bin/env python3
"""Small CLI Step 2 timing runs with a pinned executable/library pair."""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
p = argparse.ArgumentParser(__doc__)
p.add_argument('--binary', type=Path, default=ROOT/'drsolve')
p.add_argument('--library-dir', type=Path, default=ROOT)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--repeats', type=int, default=3)
p.add_argument('--reference-dir', type=Path)
a = p.parse_args()
records = []
with tempfile.TemporaryDirectory(prefix='mq-step2-') as tmp:
    for threads in (1, 4):
        for repeat in range(a.repeats):
            variants = [('current', a.binary.resolve(), a.library_dir.resolve())]
            if a.reference_dir:
                variants.insert(0, ('reference', a.reference_dir.resolve()/'drsolve', a.reference_dir.resolve()))
                if repeat % 2: variants.reverse()
            for label, binary, library in variants:
                env = dict(os.environ, LD_LIBRARY_PATH=str(library)+':'+os.environ.get('LD_LIBRARY_PATH',''))
                out = Path(tmp)/'result.dr'
                command = [str(binary), '-r', '[2]*8', '257', '--seed', '1790242205',
                           '--threads', str(threads), '--time', '-v', '2', '-o', str(out)]
                result = subprocess.run(command, cwd=ROOT, env=env, text=True, capture_output=True,
                                        check=True, timeout=180)
                assert 'Released source Dixon terms' in result.stdout
                match = re.search(r'Step 2 time: CPU time: ([\d.]+) seconds \| Wall time: ([\d.]+)', result.stdout)
                assert match, result.stdout
                canonical = '\n'.join(s for s in out.read_text().splitlines() if not s.startswith('Time: '))
                if records: assert canonical == records[0]['result']
                record = dict(label=label, threads=threads, repeat=repeat, cpu=float(match[1]), wall=float(match[2]),
                              phases=[s for s in result.stdout.splitlines() if 'Step 2 direct' in s], result=canonical)
                records.append(record)
                a.output.parent.mkdir(parents=True, exist_ok=True)
                a.output.write_text(json.dumps(records, indent=2)+'\n')
                print(label, threads, repeat, record['wall'], flush=True)

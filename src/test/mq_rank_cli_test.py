#!/usr/bin/env python3
"""Direct-index CLI activation, opt-out, precedence and exact output checks."""
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='mq-rank-cli-') as directory:
    cases = [(7, [], True), (8, [], True),
             (4, ['--no-mq-step1-shared'], False),
             (4, ['--no-mq-step1-filter'], False),
             (4, ['--mq-step1-pencil'], False),
             (4, ['--no-mq-step1-rank'], False),
             (4, ['--no-mq-step1-rank','--mq-step1-rank'], True)]
    for case, (n, extra, expected) in enumerate(cases):
        outputs = []
        for rank, threads in [(False, 1), (None, 1), (None, 4), (True, 1)]:
            path = Path(directory)/f'{case}-{rank}-{threads}.dr'
            flags = [] if rank is None else ['--mq-step1-rank' if rank else '--no-mq-step1-rank']
            command = [str(ROOT/'drsolve'), *flags,
                       *extra, '--threads', str(threads), '-v', '2', '-r', '--seed', '12345',
                       '-n', str(n), '-o', str(path), f'[2]*{n}', '65537']
            p = subprocess.run(command,cwd=ROOT,text=True,capture_output=True,check=True,timeout=180)
            active = 'MQ direct-index layer' in p.stdout
            want = ((rank is not False) or (extra and extra[-1]=='--mq-step1-rank')) and expected
            assert active == bool(want), (command,p.stdout)
            outputs.append('\n'.join(l for l in path.read_text().splitlines() if not l.startswith('Time: ')))
        assert all(s==outputs[0] for s in outputs), (n,extra)
print('MQ direct-index CLI: 28 default/activation/precedence/output checks PASS')

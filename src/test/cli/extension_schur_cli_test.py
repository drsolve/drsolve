#!/usr/bin/env python3
"""Compare complete Dixon resultants with and without extension Schur."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve()
env = dict(os.environ, LD_LIBRARY_PATH=str(binary.parent) + ':' + os.environ.get('LD_LIBRARY_PATH', ''))
with tempfile.TemporaryDirectory(prefix='extension-schur-') as work:
    output = Path(work) / 'result.dr'
    cases = [('[10]*3', '256', []), ('[10]*3', '257', []),
             ('[3]*3', '101^3', []), ('[2,3,4]', '3^3', []),
             ('[3]*4', '256', []), ('[2]*4', '256', ['--step1', '5']), ('[3]*3', '256', ['--step1', '0'])]
    cases += [('[3]*3', field, []) for field in
              ('2^4', '2^4:t^4+t^3+1', '2^16', '2^32', '2^64', '2^128', '2^5',
               '2^8:t^8+t^4+t^3+t+1', '2^9', '2^17', '2^65', '2^127',
               '2^64:t^64+t^4+t^3+t+1')]
    for degrees, field, extra in cases:
        results = []
        for enabled, threads, no_seed in ((False, 1, False), (True, 1, False),
                                          (True, 4, False), (True, 1, True)):
            command = [str(binary), '-r', degrees, field, '--seed', '12345',
                       '--no-quotient', '--threads', str(threads), '-v', '1',
                       '--mq-step4-schur' if enabled else '--no-mq-step4-schur',
                       '-o', str(output), *extra]
            if no_seed:
                command.append('--no-schur-seed')
            run = subprocess.run(command, env=env, capture_output=True, text=True, timeout=180)
            log = run.stdout + run.stderr
            assert run.returncode == 0, log
            assert 'Schur-first selection:' not in log, (command, log)
            results.append(output.read_text().split('Resultant:\n', 1)[1].splitlines()[0])
            if enabled:
                assert 'Schur' in log and ' -> ' in log, (command, log)
                assert 'compression ' not in log, (command, log)
                if field != '257':
                    assert 'Schur (extension field)' in log, log
            else:
                assert 'compression ' not in log, log
        assert len(set(results)) == 1, (degrees, field, extra)
        print('PASS', degrees, field, extra, flush=True)

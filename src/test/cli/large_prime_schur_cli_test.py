"""Large-prime Step 3 selection and exact Step 4 Schur compression."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
prime = '18446744073709551629'
with tempfile.TemporaryDirectory() as tmp:
    output = Path(tmp) / 'result.dr'
    for degrees in ('[3]*3', '[2]*4', '[2,3,4]'):
        results = []
        for flags in (['--no-mq-step4-schur'], [], ['--no-schur-seed'], ['--step4', '2']):
            run = subprocess.run([binary, '-r', degrees, prime, '--seed', '12345',
                                  '--resultant-only', '-v', '2', '--time', *flags,
                                  '-o', str(output)], capture_output=True, text=True,
                                 timeout=120, check=True)
            log = run.stdout + run.stderr
            assert re.findall(r'^Step ([1-4]):', log, re.M) == ['1', '2', '3', '4'], log
            assert re.findall(r'^Step ([1-4]) time:', log, re.M) == ['1', '2', '3', '4'], log
            assert ('Schur-first selection:' in log) == (flags in ([], ['--step4', '2'])), log
            compressed = re.search(r'Step 4 adaptive Schur: (\d+) -> (\d+)', log)
            if flags == ['--no-mq-step4-schur']:
                assert compressed is None, log
            else:
                assert compressed and int(compressed[1]) > int(compressed[2]), log
            results.append(output.read_text().split('Resultant:\n', 1)[1].strip())
        assert len(set(results)) == 1, degrees
        print('PASS large-prime Schur', degrees, flush=True)
print('Large-prime selection, Schur switches, determinant and timing checks passed.')

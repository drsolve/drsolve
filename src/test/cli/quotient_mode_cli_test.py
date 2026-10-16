"""Check construction-only --quotient and automatic -s MQ dispatch."""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve()
env = dict(os.environ, LD_LIBRARY_PATH=str(binary.parent) + ':/usr/local/lib')
with tempfile.TemporaryDirectory(prefix='quotient-mode-') as directory:
    output = Path(directory) / 'result.dr'

    def run(*options):
        p = subprocess.run([str(binary), '-r', '[2]*8', '257', '--seed', '12345',
                            '-v', '2', '-o', str(output), *options],
                           env=env, capture_output=True, text=True, timeout=120)
        log = p.stdout + p.stderr
        assert p.returncode == 0, log
        return log, output.read_text()

    log, result = run('--quotient')
    assert 'Experimental Macaulay quotient / Bezoutian' in log
    assert 'Resultant:\n' in result
    assert 'Roots in ' not in log + result
    assert 'Solution set ' not in result
    log, reference = run('--mq-quotient', '--resultant-only')
    extract = lambda text: text.split('Resultant:\n', 1)[1].splitlines()[0]
    assert extract(result) == extract(reference)

    log, result = run('-s')
    assert 'Solver: Macaulay quotient / Bezoutian' in log
    assert '[PASS]' in log and 'Found 2 solution set(s)' in log
    points = lambda text: set(re.findall(r'x0 = [^\n]*\[PASS\]', text))
    expected = points(log)
    for options in [('--quotient', '-s'), ('-s', '--quotient'),
                    ('--solver', 'dixon', '--mq-quotient')]:
        log, _ = run(*options)
        assert points(log) == expected
    log, _ = run('-s', '-n', '7')
    assert 'Solver: quotient algebra (8 equations, 7 variables)' in log
    assert '[PASS]' in log
print('Quotient mode CLI tests passed')

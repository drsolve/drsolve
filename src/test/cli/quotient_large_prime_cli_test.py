"""Large-prime quotient dispatch, recursive back substitution, and resource limits."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
prime = '18446744073709551629'
with tempfile.TemporaryDirectory() as tmp:
    def run(*args, success=True):
        proc = subprocess.run([binary, *args, '-o', str(Path(tmp) / 'result.dr')],
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True, timeout=120)
        assert (proc.returncode == 0) == success, proc.stdout
        return proc.stdout

    text = run('-r', '[2]*5', '-n', '4', '-s', prime, '--seed', '1')
    assert 'Solver: quotient algebra' in text and 'certified dimension 1,' in text
    planted = re.search(r'Planted solution: (.*)', text)[1]
    assert planted + ' [pass]' in text
    text = run('-r', '[2]*5', '-s', prime, '--seed', '1')
    assert 'Compute Dixon resultant for 5 variable(s)' in text
    assert 'Back substitution: quotient algebra' in text
    assert 'Compute Dixon resultant for 4 variable(s)' not in text
    assert '[pass]' in text and 'falling back' not in text
    text = run('-r', '[2]*5', '-n', '4', '-s', '257', '--seed', '1')
    assert 'Quotient algebra over F_257:' in text and '[PASS]' in text

    system = 'x^2+y^2-5,x*y-2'
    results = []
    for flags in ([], ['--solver', 'quotient'], ['--solver', 'dixon'],
                  ['--quotient-max-degree', '1'],
                  ['--large-prime-reconstruct', '--max-primes', '8']):
        text = run('-s', system, prime, *flags)
        if flags == ['--solver', 'quotient']:
            assert 'quotient algebra' in text
        else:
            assert 'quotient algebra' not in text
        results.append(sorted(re.findall(r'^  Solution \d+: (.*) \[pass\]$', text, re.M)))
    assert len(results[0]) == 4 and all(r == results[0] for r in results), results
    text = run('-s', 'x,y,x+y-1', prime)
    assert 'certified dimension 0,' in text and '0 rational point(s)' in text
    text = run('-s', 'x^2-x,y^2-y,(x^2-x)*y', prime,
               '--quotient-max-degree', '1', success=False)
    assert 'Quotient closure incomplete: degree limit reached' in text
    text = run('-s', 'x^20-1,y^20-1', prime, '--solver', 'quotient',
               '--quotient-memory', '1', success=False)
    assert 'Quotient closure incomplete: memory budget reached' in text
    text = run('-s', 'x-1,y-2,x+y-3', prime, '--silent')
    assert 'Progress:' not in text
print('Large-prime quotient CLI checks passed.')

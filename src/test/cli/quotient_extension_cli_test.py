"""Regression coverage for extension-field quotient dispatch and CLI limits."""
import os
import re
from pathlib import Path
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1]).resolve()
env = dict(os.environ, OMP_NUM_THREADS='1',
           LD_LIBRARY_PATH=str(binary.parent) + ':' + os.environ.get('LD_LIBRARY_PATH', ''))
with tempfile.TemporaryDirectory(prefix='drsolve-quotient-') as work:
    output = Path(work) / 'result.dr'

    def run(*args, success=True):
        result = subprocess.run([str(binary), *args, '-o', str(output)],
                                text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, env=env, timeout=60)
        assert (result.returncode == 0) == success, result.stdout
        return result.stdout

    text = run('-r', '[2]*8', '256', '-n', '7', '--solver', 'quotient', '--seed', '123')
    assert 'GF(2^8) lookup tables' in text and 'certified dimension 1,' in text
    assert '[PASS]' in text and '1 rational point(s)' in text
    for q, backend in [('9', 'General finite field'), ('1062961', 'General finite field')]:
        text = run('x-t,y-t-1', q, '--solver', 'quotient', '--quotient-max-degree', '2')
        assert backend in text and 'x = t, y = t+1 [PASS]' in text
    text = run('x-2,y-3', '257', '--solver', 'quotient')
    assert 'Quotient algebra over F_257:' in text and '[PASS]' in text
    for q in ['9', '256']:
        text = run('x*y,0', q, '--solver', 'quotient', '--quotient-max-degree', '3', success=False)
        assert 'Quotient closure incomplete: degree limit reached' in text
        text = run('x*y,0', q, '--solver', 'quotient', '--quotient-memory', '1', success=False)
        assert 'Quotient closure incomplete: memory budget reached' in text
    text = run('x-t,y-t-1', '9', '--solver', 'quotient', '--silent')
    assert 'Progress:' not in text and 'Field backend:' not in text
    text = run('-r', '[2]*7', '256', '-s', '--seed', '123')
    assert 'Solver: Dixon elimination (7 equations, 7 variables)' in text
    assert 'Back substitution: quotient algebra (7 equations, 6 variables)' in text
    assert '[PASS]' in text and 'falling back' not in text
    for q in ['9', '256']:
        system = 'x+y,(y-t)*(y-t-1)'
        text = run(system, q)
        assert 'quotient algebra' not in text and 'Found 2 solution set(s)' in text
        text = run(system, q, '--quotient-max-degree', '1')
        assert 'quotient algebra' not in text and 'Found 2 solution set(s)' in text
        text = run(system, q, '--solver', 'dixon')
        assert 'quotient algebra' not in text and 'Found 2 solution set(s)' in text

        system = '(a-t)*x,(a-t)*(x-1)'
        def resultant():
            return re.search(r'\nResultant:\n([^\n]+)', output.read_text())[1]

        text = run(system, 'x', q, '--no-sift', '--resultant-only')
        raw = resultant()
        text = run(system, 'x', q, '--sift')
        sifted = resultant()
        assert raw != sifted and 'Quotient sift: certified removed multiplier = 1' not in text
        for content in [text, output.read_text()]:
            assert 'Roots in a (degree 1)' in content and '(Multiplicity: 1)' in content
            assert '(Multiplicity: 2)' not in content
        text = run(system, 'x', q, '--sift', '--sift-max-degree', '1')
        assert 'certified removed multiplier = 1' in text
        assert 'Roots in a (degree 2)' in text and '(Multiplicity: 2)' in text
        text = run(system, 'x', q, '--sift', '--silent', '--resultant-only')
        assert text.strip() == resultant() == sifted
        # The coefficient generator is not a remaining parameter.
        text = run('x-t,x-t-1', 'x', q, '--sift', '--sift-max-degree', '4', success=False)
        assert 'sift limits require active sifting' in text
print('Quotient extension CLI tests passed')

"""Exercise sift activation, output consistency, limits, and paper fixtures."""
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1]).resolve()
root = Path(__file__).resolve().parents[3]
fixtures = root / 'examples' / 'sift'
env = dict(os.environ, OMP_NUM_THREADS='4')

with tempfile.TemporaryDirectory(prefix='drsolve-sift-') as work:
    output = Path(work) / 'result.dr'

    def run(name, *options, code=0):
        p = subprocess.run([str(binary), str(fixtures / (name + '.dr')),
                            *options, '-o', str(output)],
                           text=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, env=env, timeout=60)
        assert p.returncode == code, (name, options, p.returncode, p.stdout)
        saved = output.read_text() if code == 0 else ''
        result = re.search(r'\nResultant:\n([^\n]+)', saved)
        return p.stdout, saved, result[1] if result else None

    expected = 'b^3*a+b^2*a^2+65536*b^2*a+65536*b*a^2+a'
    stdout, saved, default = run('two_parameters', '--resultant-only')
    assert 'Quotient sift:' not in stdout + saved
    stdout, saved, sifted = run('two_parameters', '--sift', '--resultant-only')
    assert '(limit 4)' in saved
    assert sifted == expected
    assert 'certified removed multiplier = b^3*a^2' in saved
    _, _, raw = run('two_parameters', '--no-sift', '--resultant-only')
    assert raw == default
    assert raw.replace(' ', '') != sifted
    stdout, saved, budgeted = run('two_parameters', '--sift', '--sift-max-degree', '1', '--resultant-only')
    assert budgeted.replace(' ', '') == raw.replace(' ', '')
    assert 'Unproved factors retained' in saved
    stdout, _, silent = run('two_parameters', '--sift', '--silent', '--resultant-only')
    assert stdout.strip() == silent == expected

    expected = {
        'example34': 'a^3+59740*a^2+15124*a',
        'example35': 'z^2+32768*z+32768',
        'example36_b5_slice': 'T^6+27274*T^5+10988*T^4+44007*T^3+18108*T^2+25450*T+1190',
        'example36_a1_slice': 'T^6+1207*T^5+48563*T^4+2519*T^3+9806*T^2+20400*T+60521',
    }
    for name, wanted in expected.items():
        stdout, saved, _ = run(name, '--resultant-only')
        assert 'Quotient sift:' not in stdout + saved  # Single parameter is opt-in.
        _, saved, sifted = run(name, '--sift', '--resultant-only')
        assert sifted == wanted, (name, sifted)
        assert 'Quotient sift:' in saved

    stdout, saved, _ = run('example35', '--sift')
    for text in (stdout, saved):
        assert 'Roots in z (degree 2)' in text
        assert 'Find 2 roots:' in text
        assert '(Multiplicity: 3)' not in text
    stdout, saved, _ = run('example35', '--sift', '--silent')
    assert stdout.strip() == '' and 'Roots in z (degree 2)' in saved
    run('example34', '--sift-max-degree', '0', code=1)
    run('example34', '--sift-max-degree', code=1)
    run('example34', '--sift-max-degree', '3', code=1)  # Limits don't implicitly enable --sift.
    run('two_parameters', '--sift-max-degree', '3', code=1)
    run('example34', '--sift', '--solver', 'quotient', code=1)

    # Hold construction fixed while comparing Step 4 dispatch; the --method
    # flag also prevents explicit --step4 from changing the construction mode.
    # Disable inverse-block so the HNF selection under test is not intercepted.
    for field in ['101', '2^3']:
        native_results = []
        for extra in [[], ['--step4', '0'], ['--step4', '1']]:
            proc = subprocess.run([
                str(binary), 'x^2+y+z,x*y+y*z+1,x+y+z', 'x,y', field,
                '--method', '5', '--no-inverse-block', '--no-sift', '--resultant-only', '-v', '3',
                '-o', str(output), *extra], text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, env=env, timeout=60)
            assert proc.returncode == 0, proc.stdout
            step4 = proc.stdout.split('Step 4: Compute resultant', 1)[1]
            method = 'fast HNF' if extra == ['--step4', '1'] else 'minor expansion'
            assert 'Determinant method: ' + method in step4, proc.stdout
            native_results.append(re.search(r'\nResultant:\n([^\n]+)', output.read_text())[1])
        assert len(set(native_results)) == 1, (field, native_results)

    # The full ten-parameter example has an 8x8 minor and uses expansion. It
    # uses an explicit small degree bound to keep this regression inexpensive.
    stdout, saved, result = run('example36', '--sift', '--resultant-only', '--sift-max-degree', '3')
    assert result and 'Determinant method: minor expansion' in stdout.split('Step 4:', 1)[1]
    assert 'Quotient sift:' in saved and 'Unproved factors retained' in saved
    assert 'Sifted resultant: 282930 terms (not displayed)' in stdout
    assert 'workspace budget' not in stdout + saved
    assert 'characters (not displayed)' not in stdout

print('Quotient sift CLI and paper fixture tests passed')

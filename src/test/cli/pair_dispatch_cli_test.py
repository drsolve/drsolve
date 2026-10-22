"""Default pair elimination policy, solver recursion, and explicit overrides."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
big = '18446744073709551629'
system = 'x+y,y^2-y'
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / 'result.dr'
    def run(*args):
        p = subprocess.run([binary, *args, '-o', str(path)], capture_output=True,
                           text=True, check=True, timeout=120)
        return p.stdout + p.stderr, path.read_text()

    for field in ('2', '257', '18446744073709551557', '16', '256',
                  '2^8:t^8+t^4+t^3+t+1', '65536', '2^127', '9', big):
        want_subres = field in ('9', big)
        log, saved = run(system, 'x', field, '--resultant-only')
        assert ('Subresultant: computing resultant' in log) == want_subres, (field, log)
        assert ('Step 1:' in log) != want_subres, (field, log)
        _, explicit = run(system, 'x', field, '--subres', '--resultant-only')
        if field != '2':  # Existing subres parsing reduces y^2-y by the F_2 field equation.
            assert re.sub(r'\s+', '', saved.split('Resultant:\n')[1]) == re.sub(r'\s+', '', explicit.split('Resultant:\n')[1]), field
        log, _ = run('-s', system, field)
        if field == big:
            assert 'Compute subresultant for 2 variable(s)' in log, log
        else:
            expected = 'subresultant' if field == '9' else 'Dixon'
            assert 'Pair elimination: ' + expected in log, (field, log)
        assert 'Found 2 solution set(s)' in log, (field, log)

    log, _ = run('-s', 'x+y,(y-t)*(y-t-1)', '256')
    assert 'Pair elimination: Dixon' in log and 'Found 2 solution set(s)' in log

    for flags in (['--dixon'], ['--step4', '1'], ['--large-prime-reconstruct', '--max-primes', '8']):
        log, _ = run(system, 'x', big, '--resultant-only', *flags)
        assert 'Subresultant: computing resultant' not in log, log
    for flags in (['--solver', 'dixon'], ['--dixon'], ['--step4', '1'],
                  ['--large-prime-reconstruct', '--max-primes', '8']):
        log, _ = run('-s', system, big, *flags)
        assert 'Compute Dixon resultant for 2 variable(s)' in log, log
        assert 'Found 2 solution set(s)' in log, log
    log, _ = run('-s', system, big, '--solver', 'quotient')
    assert 'Solver: quotient algebra' in log and 'Compute subresultant' not in log
    # Large-prime pairs with multiple remaining parameters also default to subres.
    log, _ = run('x+y,x+z', 'x', big, '--resultant-only')
    assert 'Subresultant: computing resultant' in log
    log, _ = run('-r', '[30]*2', big, '--seed', '1', '--resultant-only')
    assert 'Subresultant: computing resultant' in log and '901 terms, total degree 900' in log
print('Pair default dispatch and explicit overrides passed.')

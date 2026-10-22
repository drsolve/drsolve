"""Direct FLINT subresultants over large prime fields, matching Dixon output."""
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
prime = '18446744073709551629'
with tempfile.TemporaryDirectory() as tmp:
    path = Path(tmp) / 'result.dr'
    def run(args, success=True):
        proc = subprocess.run([binary, *args, '-o', str(path)], text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=180)
        assert (proc.returncode == 0) == success, proc.stdout
        result = path.read_text().split('Resultant:\n', 1)[1].strip() if success else None
        return proc.stdout, result

    cases = [(['-r', '[3]*2', '--seed', '1', prime], None),
             (['x^2+y^2-5,x*y-2', 'x', prime], None),
             (['x+y,y+z', 'y', prime], None),
             (['x^2+1,x^2+2', 'x', prime], '1'),
             (['x*y,x*(y+1)', 'x', prime], '0'),
             (['y,x^2+y', 'x', prime], None),
             (['x+y,x-y', 'x', '170141183460469231731687303715884105727'], None)]
    for args, expected in cases:
        log, subres = run([*args, '--subres', '--resultant-only', '--time'])
        assert 'Subresultant: computing resultant' in log and 'FLINT fmpz_mod' in log
        assert 'Step 1:' not in log and 'Subresultant time:' in log
        _, dixon = run([*args, '--resultant', 'dixon', '--resultant-only'])
        assert subres == dixon, (args, subres, dixon)
        if expected is not None:
            assert subres == expected
    log, _ = run(['-r', '[30]*2', prime, '--resultant', 'subres', '--seed', '1', '--resultant-only'])
    assert '901 terms, total degree 900' in log and 'Step 1:' not in log
    log, result = run(['x^2+y^2-5,x*y-2', 'x', prime, '--subres', '--silent', '--resultant-only'])
    assert not log.strip() and result
    log, _ = run(['x+y,x-y', 'x', prime, '--subres'])
    assert 'Roots in y' in path.read_text()
    for args, message in [(['x+y,x-y,x-1', 'x', prime, '--subres'], 'exactly 2 polynomials'),
                          (['x+y,x-y', 'x', prime, '--subres', '--large-prime-reconstruct'], 'cannot combine')]:
        log, _ = run(args, success=False)
        assert message in log, log
print('Large-prime FLINT subresultant CLI checks passed.')

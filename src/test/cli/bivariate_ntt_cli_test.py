"""Default Dixon/Bezout fast resultant over prime and extension fields."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())

with tempfile.TemporaryDirectory() as tmp:
    output = Path(tmp) / 'result.dr'

    def run(args, flags):
        proc = subprocess.run([binary, *args, *flags, '--resultant-only',
                               '-o', str(output)], capture_output=True, text=True,
                              check=True, timeout=120)
        result = re.sub(r'\s+', '', output.read_text().split('Resultant:\n', 1)[1])
        return proc.stdout + proc.stderr, result

    cases = [(['-r', f'[{d}]*2', '65537', '--seed', str(seed)])
             for d in (2, 3, 5, 8) for seed in (1, 7)]
    cases += [
        ['x+y,x-y', 'x', '65537'],
        ['2*x^3+y^3+1,3*x^2+y+2', 'x', '65537'],
        ['x^2+y^5+1,x+y^3+2', 'x', '65537'],
        ['(x+y)*(x+1),(x+y)*(x+2)', 'x', '65537'],
        ['x^2+y^2,x^2+y^2+1', 'x', '65537'],
        ['x+y,x+y+1', 'x', '65537'],
        ['x+y,x-y', 'x', '18446744069414584321'],
        ['x^3+y^3+1,x^2+y+2', 'x', '101'],
        ['x^3+y^3+1,x^2+y+2', 'x', '18446744073709551629'],
        ['x^3+y^3+1,x^2+y+2', 'x', '170141183460469231731687303715884105727'],
        ['x^3+(t+1)*y^3+1,x^2+t*y+2', 'x', '65537^2'],
        ['x^3+(t+1)*y^3+1,x^2+t*y+2', 'x', '25'],
        ['x^3+(t+1)*y^3+1,x^2+t*y+1', 'x', '256'],
        ['x^3+(t+1)*y^3+1,x^2+t*y+1', 'x', '2^8:t^8+t^4+t^3+t+1'],
        ['x^3+(t+1)*y^3+1,x^2+t*y+1', 'x', '2^127'],
        ['x+y,x-y', 'x', '9'],
    ]
    cases += [['-r', '[8]*2', field, '--seed', '19']
              for field in ('256', '65537^2', '18446744073709551629')]
    for args in cases:
        _, reference = run(args, ['--subres'])
        for threads in ('1', '4'):
            log, actual = run(args, ['--threads', threads])
            assert 'Dixon/Bezout fast time:' in log and 'fallback' not in log, log
            assert actual == reference, (args, actual, reference)

    # Insufficient points, variable leading coefficients, zero input and
    # extra parameters retain the original matrix algorithm.
    for args in [
        ['x^3+y^3,x^3+y+1', 'x', '7'],
        ['x^4+y^4+1,x^4+y+2', 'x', '9'],
        ['y*x^2+x+1,x^2+y', 'x', '65537'],
        ['0,x+y', 'x', '65537'],
        ['x+y+z,x-y', 'x', '65537'],
    ]:
        _, reference = run(args, ['--dixon'])
        log, actual = run(args, ['--bivar-ntt'])
        assert 'Dixon/Bezout fast time:' not in log and 'Step 1:' in log, log
        assert actual == reference, (args, actual, reference)
    for flags in (['--no-bezout-fast'], ['--no-bivar-ntt'], ['--step4', '1']):
        log, reference = run(cases[0], flags)
        assert 'Dixon/Bezout fast time:' not in log and 'Step 1:' in log, log
        _, actual = run(cases[0], [])
        assert actual == reference
    for flag in ('--bezout-fast', '--bivar-ntt'):
        log, actual = run(cases[0], [flag])
        assert 'Dixon/Bezout fast time:' in log and actual == reference
    log, result = run(cases[0], ['--bezout-fast', '--silent'])
    assert re.sub(r'\s+', '', log) == result, log

print('Dixon/Bezout fast default, field, parallel and fallback checks passed.')

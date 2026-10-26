#!/usr/bin/env python3
"""Extension Bezoutian construction, compared independently with Dixon."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve()
env = dict(os.environ, OMP_NUM_THREADS='1',
           LD_LIBRARY_PATH=str(binary.parent) + ':' + os.environ.get('LD_LIBRARY_PATH', ''))
with tempfile.TemporaryDirectory(prefix='extension-bezoutian-') as work:
    output = Path(work) / 'result.dr'

    def run(arguments, flag=None, threads=1):
        p = subprocess.run([str(binary), *arguments, *([flag] if flag else []), '--threads', str(threads),
                            '-v', '3' if '-s' in arguments else '1', '-o', str(output)], env=env,
                           capture_output=True, text=True, timeout=120)
        log = p.stdout + p.stderr
        assert p.returncode == 0, log
        if '-s' in arguments:
            return '', log
        resultant = output.read_text().split('Resultant:\n', 1)[1].splitlines()[0]
        return resultant, log

    def compare(arguments, fallback=None):
        # Bivariate block and fast resultants bypass the quotient construction and its
        # fallback checks; compare the quotient and Dixon matrix paths here.
        arguments = [*arguments, '--no-bezout-block', '--no-bezout-fast']
        reference, _ = run(arguments, '--no-quotient')
        for threads in (1, 4):
            result, log = run(arguments, '--quotient', threads)
            assert result == reference, (arguments, threads, log)
            if fallback:
                assert 'fallback:' in log and fallback in log, log
            else:
                assert 'Extension field;' in log and 'fallback:' not in log, log
                assert 'normalize to monic.' in log, log
        print('PASS', ' '.join(arguments), flush=True)

    for field in ('2^4', '2^8', '2^8:t^8+t^4+t^3+t+1', '2^16', '2^32', '2^64', '2^128', '101^3',
                  '2^7', '2^9', '2^17', '2^31', '2^63', '2^65', '2^127',
                  '2^64:t^64+t^4+t^3+t+1'):
        compare(['-r', '[2]*3', field, '--seed', '12345'])
    for degrees, field in (('[3]*3', '2^8'), ('[2,3,4]', '2^8'),
                           ('[1,2,3]', '2^8'), ('[3]*3', '3^3'),
                           ('[2]*3', '2^5'), ('[2]*2', '3^2'),
                           ('[2]*4', '101^3'), ('[2]*4', '2^16')):
        compare(['-r', degrees, field, '--seed', '12345'])
    # Small MQ inputs must actually use the shared-frame kernel.
    for count in (3, 4):
        arguments = ['-r', f'[2]*{count}', '256', '--seed', '12345', '--time']
        reference, log = run(arguments, '--quotient', 4)
        assert 'shared leading-frame MQ kernel' in log and 'Compact byte minors' in log, log
        assert 'fallback:' not in log, log
        automatic, log = run(arguments, threads=4)
        assert automatic == reference
        assert 'shared leading-frame MQ kernel' in log and 'fallback:' not in log, log
    for options in (['--no-quotient'], ['--no-mq-quotient'], ['--step4', '1']):
        _, log = run(['-r', '[2]*3', '256', '--seed', '12345', *options])
        assert 'Experimental Dixon resultant via quotient-space Bezoutian' not in log, log
    for degrees, field in (('[3]*3', '256'), ('[2]*3', '4')):
        _, log = run(['-r', degrees, field, '--seed', '12345'])
        assert 'Experimental Dixon resultant via quotient-space Bezoutian' not in log, log
    # Coefficients explicitly outside the prime subfield, with an exact oracle.
    compare(['x-t,y-t-1', 'x', '2^8'])
    # Too few distinct field elements, including the two validation points.
    compare(['-r', '[2]*3', '4', '--seed', '12345'], 'field too small')
    # RREF admission and interpolation-degree failure both preserve Dixon fallback.
    compare(['x^3+y^3+z,x^3+y^3+x+1,x^3+y^3+y+t', 'x,y', '2^8'], 'rank')
    compare(['-r', '[2]*3', '9', '--seed', '12345'], 'leading Macaulay')
    # Square solve dispatch must retain all solutions, including back substitution.
    system = ['-r', '[2]*3', '2^8', '--seed', '12345', '-s']
    solutions = []
    for flag in ('--quotient', '--no-quotient'):
        _, log = run(system, flag)
        if flag == '--quotient':
            assert 'Extension field;' in log and 'fallback:' not in log, log
        solutions.append(output.read_text().split('=== Polynomial System Solutions ===', 1)[1]
                         .split('=== Solution Complete ===', 1)[0])
    assert solutions[0] == solutions[1]
print('Extension quotient CLI tests passed')

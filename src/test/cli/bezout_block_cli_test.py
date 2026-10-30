"""Default nmod Dixon/Bezout inverse-Hankel reconstruction and safe fallback."""
import re
import os
from math import comb
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / 'result.dr'
    def run(args, flags, env=None):
        proc = subprocess.run([binary, *args, *flags, '--threads', '1',
                               '--resultant-only', '-o', str(out)],
                              capture_output=True, text=True, check=True, timeout=120, env=env)
        result = re.sub(r'\s+', '', out.read_text().split('Resultant:\n', 1)[1])
        return proc.stdout + proc.stderr, result
    successes = 0
    for prime in ['65537', '998244353', '1000000007', '18446744069414584321']:
        for d in [2, 5, 8, 16]:
            args = ['-r', f'[{d}]*2', prime, '--seed', '19']
            _, ref = run(args, ['--no-bezout-block'])
            log, got = run(args, [])
            explicit_log, explicit = run(args, ['--bezout-block'])
            assert explicit == got and ('Dixon/Bezout block time:' in explicit_log) == ('Dixon/Bezout block time:' in log)
            assert got == ref, (args, log, got, ref)
            successes += 'Dixon/Bezout block time:' in log
    assert successes >= 12, successes
    unequal_successes = 0
    for prime in ['65537', '1000000007']:
        for a, b in [(5, 1), (9, 4), (4, 9), (25, 8), (35, 27)]:
            args = ['-r', f'[{a},{b}]', prime, '--seed', '19']
            _, ref = run(args, ['--no-bezout-block'])
            log, got = run(args, [])
            assert got == ref, (args, log)
            unequal_successes += 'C reduced spaced projection' in log
            _, got = run(args, ['--no-pml-fft'])
            assert got == ref, args
    assert unequal_successes >= 8, unequal_successes
    # Parser insertion order and elimination-variable permutation must not
    # make FLINT read a zero degree from an unnormalized native polynomial.
    degree = 27
    for variable in ('x', 'y'):
        terms = [f'{comb(degree, k)}*x^{k}*y^{degree-k}'
                 for k in range(degree+1)]
        for ordered in (terms, list(reversed(terms))):
            power = '+'.join(ordered)
            args = [f'{power}+x+1,{power}+y+2', variable, '2130706433']
            _, ref = run(args, ['--no-bezout-block'])
            log, got = run(args, [])
            assert got == ref, (args, log)
            match = re.search(r'C adaptive interpolation, block \d+, (\d+) points', log)
            assert match and int(match[1]) <= 2*degree+1, log
    # New backend can be explicitly disabled, then re-enabled (last flag wins).
    for prime in [998244353, 1000000007]:
        args = ['-r', '[32]*2', str(prime), '--seed', '19']
        _, ref = run(args, ['--bezout-block', '--no-pml-fft'])
        for flags in [['--bezout-block', '--pml-fft'],
                      ['--bezout-block', '--no-pml-fft', '--pml-fft']]:
            _, got = run(args, flags)
            assert got == ref, (prime, flags)

    # Unsupported field/shape, shared factors and too few interpolation points.
    for args in [
        ['x+y,x-y', 'x', '17'],
        ['(x+y)*(x+1),(x+y)*(x+2)', 'x', '65537'],
        ['y*x^2+x+1,x^2+y', 'x', '65537'],
        ['x^3+y^3,x^3+y+1', 'x', '7'],
        ['x^3+y^3+1,x^2+y+2', 'x', '65537'],
        ['3*x^7+2*x*y^6+y^7+1,5*x^3+x*y+2', 'x', '65537'],
        ['5*x^3+x*y+2,3*x^7+2*x*y^6+y^7+1', 'x', '65537'],
        ['x*y^6+y^7+1,5*x^3+x*y+2', 'x', '65537'],
        ['(x+y)*(x^4+y^4+1),3*(x+y)', 'x', '65537'],
        ['x^7+y^7+1,y*x+y^2+1', 'x', '65537'],
        ['x^2+y^2+1,x^2+x+y+2', 'x', '65537^2'],
        ['x^2+y^2+1,x^2+x+y+2', 'x', '18446744073709551629'],
    ]:
        _, ref = run(args, ['--no-bezout-block'])
        _, got = run(args, [])
        assert got == ref, (args, got, ref)
    args = ['-r', '[8]*2', '65537', '--seed', '19']
    _, ref = run(args, [])
    log, got = run(args, ['--bezout-block', '--no-bezout-block'])
    assert 'block-Hankel' not in log and got == ref, log
    log, got = run(args, ['--no-bezout-block', '--bezout-block'])
    assert 'Dixon/Bezout block time:' in log and got == ref, log
    # An underestimated initial budget must retry, never accept a wrong result.
    for initial in ('1', '1000', '-1', 'not-a-number'):
        env = dict(os.environ, DRSOLVE_BEZOUT_INITIAL_POINTS=initial,
                   DRSOLVE_PROFILE_BEZOUT='1')
        log, got = run(args, [], env)
        assert got == ref, (initial, log)
        if initial == '1':
            assert 'attempt 1: points=1' in log and 'status=accepted' in log, log
            assert re.search(r'attempt [2-9]:', log), log
        if initial in ('-1', 'not-a-number'):
            assert 'Invalid DRSOLVE_BEZOUT_INITIAL_POINTS' in log, log
    # CLI overrides the environment; zero explicitly restores automatic sizing.
    env = dict(os.environ, DRSOLVE_BEZOUT_INITIAL_POINTS='not-a-number',
               DRSOLVE_PROFILE_BEZOUT='1')
    for flags in (['--bezout-initial-points', '1'], ['--bezout-initial-points=1'],
                  ['--bezout-initial-points', '0']):
        log, got = run(args, flags, env)
        assert got == ref and 'Invalid DRSOLVE_BEZOUT_INITIAL_POINTS' not in log, log
        if flags[-1] != '0':
            assert 'attempt 1: points=1 ' in log, log
        else:
            assert 'attempt 1: points=17 ' in log, log
    for bad in ('-1', 'abc', '18446744073709551616', ''):
        proc = subprocess.run([binary, '--bezout-initial-points='+bad, *args],
                              capture_output=True, text=True, timeout=10)
        assert proc.returncode != 0 and 'nonnegative integer' in proc.stderr, proc
    log, got = run(args, ['--silent'])
    assert re.sub(r'\s+', '', log) == got == ref, log
print(f'Bezout block CLI: {successes} equal and {unequal_successes} unequal direct successes; exact outputs, fallback and silent mode PASS')

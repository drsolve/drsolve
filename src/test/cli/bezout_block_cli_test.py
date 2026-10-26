"""Default nmod Dixon/Bezout inverse-Hankel reconstruction and safe fallback."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / 'result.dr'
    def run(args, flags):
        proc = subprocess.run([binary, *args, *flags, '--threads', '1',
                               '--resultant-only', '-o', str(out)],
                              capture_output=True, text=True, check=True, timeout=120)
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
    log, got = run(args, ['--silent'])
    assert re.sub(r'\s+', '', log) == got == ref, log
print(f'Bezout block CLI: {successes} equal and {unequal_successes} unequal direct successes; exact outputs, fallback and silent mode PASS')

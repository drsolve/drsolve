"""Opt-in nmod Dixon/Bezout inverse-Hankel reconstruction and safe fallback."""
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
            _, ref = run(args, ['--bezout-fast'])
            log, got = run(args, ['--bezout-block'])
            assert got == ref, (args, log, got, ref)
            successes += 'Dixon/Bezout block time:' in log
    assert successes >= 12, successes
    # Unsupported field/shape, shared factors and too few interpolation points.
    for args in [
        ['x+y,x-y', 'x', '17'],
        ['(x+y)*(x+1),(x+y)*(x+2)', 'x', '65537'],
        ['y*x^2+x+1,x^2+y', 'x', '65537'],
        ['x^3+y^3,x^3+y+1', 'x', '7'],
        ['x^3+y^3+1,x^2+y+2', 'x', '65537'],
        ['x^2+y^2+1,x^2+x+y+2', 'x', '65537^2'],
        ['x^2+y^2+1,x^2+x+y+2', 'x', '18446744073709551629'],
    ]:
        _, ref = run(args, [])
        _, got = run(args, ['--bezout-block'])
        assert got == ref, (args, got, ref)
    args = ['-r', '[8]*2', '65537', '--seed', '19']
    _, ref = run(args, [])
    log, got = run(args, ['--bezout-block', '--no-bezout-block'])
    assert 'block-Hankel' not in log and got == ref, log
    log, got = run(args, ['--bezout-block', '--silent'])
    assert re.sub(r'\s+', '', log) == got == ref, log
print(f'Bezout block CLI: {successes} direct successes; exact outputs, fallback and silent mode PASS')

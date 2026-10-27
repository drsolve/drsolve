"""Dixon inverse-block dispatch, exact outputs, manual sizes and safe fallback."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / 'result.dr'

    def run(args, flags=(), threads=1):
        p = subprocess.run([binary, *args, *flags, '--threads', str(threads),
                            '--resultant-only', '-v', '2', '-o', str(out)],
                           capture_output=True, text=True, check=True, timeout=180)
        value = re.sub(r'\s+', '', out.read_text().split('Resultant:\n', 1)[1])
        return p.stdout + p.stderr, value

    compressed = lambda text: bool(re.search(r'(?:Step 4 (?:adaptive )?Schur:|Direct coefficient core:).*->', text))
    for degrees in ['[8]*3', '[4]*4', '[3]*5']:
        args = ['-r', degrees, '1073741827', '--seed', '1234']
        baseline_log, ref = run(args, ['--no-inverse-block'])
        assert compressed(baseline_log) and 'Dixon inverse-block time:' not in baseline_log, baseline_log
        log, got = run(args)
        assert got == ref
        assert ('Dixon inverse-block time:' in log) == (degrees != '[3]*5'), log
        assert compressed(log), log
        if degrees != '[3]*5':
            compression = re.search(r'(?:Step 4 (?:adaptive )?Schur:|Direct coefficient core:).*->', log)
            assert compression.start() < log.index('Dixon inverse-block:'), log
        for flags, enabled in [(['--step4-schur'], True),
                               (['--no-step4-schur'], False),
                               (['--no-step4-schur', '--step4-schur'], True),
                               (['--step4-schur', '--no-step4-schur'], False)]:
            log, got = run(args, flags)
            assert got == ref and compressed(log) == enabled, (flags, log)
        for flags, threads in [(['--inverse-block-size', '8'], 1),
                               (['--inverse-block'], 4),
                               (['--inverse-block-size=8'], 1),
                               (['--no-inverse-block', '--inverse-block'], 1),
                               (['--inverse-block', '--no-inverse-block'], 1)]:
            log, got = run(args, flags, threads)
            assert got == ref, (degrees, flags, log)
            assert ('Dixon inverse-block time:' in log) == (flags[-1] != '--no-inverse-block'), log
        log, got = run(args, ['--inverse-block-size', '99999'])
        assert got == ref and 'fallback' in log and compressed(log), log
        compression = re.search(r'(?:Step 4 (?:adaptive )?Schur:|Direct coefficient core:).*->', log)
        assert compression.start() < log.index('Dixon inverse-block:'), log
        assert 'Dixon inverse-block time:' not in log, log
        log, got = run(args, ['--inverse-block-size', '99999', '--no-step4-schur'])
        assert got == ref and 'fallback' in log and not compressed(log), log
        log, got = run(args, ['--inverse-block', '--inverse-block-auto'])
        assert got == ref
        assert ('Dixon inverse-block time:' in log) == (degrees != '[3]*5'), log

    args = ['-r', '[8]*3', '1073741827', '--seed', '1234']
    _, ref = run(args, ['--no-inverse-block'])
    log, got = run(args, ['--fq-det-method', 'hnf'])
    assert got == ref and 'Dixon inverse-block time:' not in log and compressed(log)
    log, got = run(args, ['--no-step4-schur', '--inverse-block'])
    assert got == ref and 'Dixon inverse-block time:' in log

    for args in [['-r', '[4]*3', '17', '--seed', '19'],
                 ['-r', '[10]*3', '2^8', '--seed', '19'],
                 ['-r', '[8]*2', '65537', '--seed', '19']]:
        _, ref = run(args, ['--no-inverse-block'])
        log, got = run(args, ['--inverse-block'])
        assert got == ref and 'Dixon inverse-block time:' not in log, log

    for bad in ['-1', 'abc', '1.5', '999999999999999999999999999999']:
        p = subprocess.run([binary, '-r', '[4]*3', '65537', '--inverse-block-size', bad],
                           capture_output=True, text=True)
        assert p.returncode != 0 and 'nonnegative integer' in p.stderr, p.stderr
print('inverse-block CLI: Schur-first defaults, fallback/overrides and exact outputs passed')

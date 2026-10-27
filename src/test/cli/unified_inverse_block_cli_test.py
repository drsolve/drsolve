"""Cross-field inverse-block defaults, [10]*3 resultants and safe fallback."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
with tempfile.TemporaryDirectory(prefix='unified-inverse-block-') as tmp:
    output = Path(tmp) / 'result.dr'

    def run(q, flags=(), degrees='[10]*3', threads=4):
        command = [binary, '-r', degrees, q, '--seed', '1234', '--threads', str(threads),
                   '--resultant-only', '-v', '2', '-o', str(output), *flags]
        p = subprocess.run(command, capture_output=True, text=True, timeout=900)
        log = p.stdout + p.stderr
        assert p.returncode == 0, (command, log)
        value = re.sub(r'\s+', '', output.read_text().split('Resultant:\n', 1)[1])
        return log, value

    for q in ('65537', '101^3', '2^16', '2^32', '18446744073709551629', '2^8'):
        baseline, reference = run(q, ['--no-inverse-block'])
        log, got = run(q)
        assert got == reference, q
        success = 'Dixon inverse-block time:' in log
        assert success == (q != '2^8'), (q, log)
        if not success:
            assert 'fallback' in log and 'Schur' in log, log
        else:
            compression = re.search(r'(?:Step 4 (?:adaptive )?Schur|Direct coefficient core:).*->', log)
            assert compression and compression.start() < log.index('Dixon inverse-block:'), log
        details = re.findall(r'Dixon inverse-block(?: time)?:[^\n]+', log)
        print('PASS [10]*3', q, '; '.join(details), flush=True)

    for q in ('101^3', '2^16', '18446744073709551629'):
        _, reference = run(q, ['--no-inverse-block'], '[4]*3')
        for flags in (['--inverse-block-size', '3'], ['--step4-schur'], ['--no-step4-schur'],
                      ['--inverse-block-size', '99999'],
                      ['--inverse-block-size', '99999', '--no-step4-schur'],
                      ['--fq-det-method', 'hnf']):
            log, got = run(q, flags, '[4]*3', 1)
            assert got == reference, (q, flags, log)
            if '99999' in flags or 'hnf' in flags:
                assert 'Dixon inverse-block time:' not in log, log
            else:
                assert 'Dixon inverse-block time:' in log, log
            if flags == ['--no-step4-schur']:
                assert not re.search(r'Step 4 (?:adaptive )?Schur.*->', log), log
            if flags == ['--step4-schur']:
                assert re.search(r'Step 4 (?:adaptive )?Schur.*->', log), log
        for degrees, expected in (('[2]*4', True), ('[2]*5', False)):
            _, reference = run(q, ['--no-inverse-block', '--no-quotient'], degrees)
            log, got = run(q, ['--no-quotient'], degrees=degrees)
            assert got == reference, (q, degrees)
            assert ('Dixon inverse-block time:' in log) == expected, (q, degrees, log)
            if expected:
                compression = re.search(r'Step 4 (?:adaptive )?Schur.*->', log)
                assert compression and compression.start() < log.index('Dixon inverse-block:'), log
        print('PASS switches/default equation counts', q, flush=True)

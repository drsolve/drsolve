"""Degree-only complexity reports must avoid unused Schur layer sums."""
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = Path(sys.argv[1]).resolve()
env = os.environ.copy()
libpath = "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
env[libpath] = str(binary.parent) + os.pathsep + env.get(libpath, "")
with tempfile.TemporaryDirectory() as tmp:
    def report(degrees, *flags, field="65537"):
        output = Path(tmp) / 'report.dr'
        result = subprocess.run(
            [str(binary), '-c', '-r', degrees, *([] if field is None else [field]), *flags, '-o', str(output)],
            cwd=tmp, env=env, text=True, capture_output=True, check=True, timeout=10)
        assert 'Generating random polynomial' not in result.stdout
        return result.stdout, output.read_text()

    for flags in ((), ('-v', '0'), ('--rank-pred',), ('-v', '3')):
        # The unused cubic Schur sums make this case prohibitively slow.
        log, saved = report('[4000]*3', *flags)
        assert re.search(r'(?:Overall|Total).*complexity', saved)
        assert 'degrees=[4000]*3' in saved
        if flags != ('-v', '0'):
            assert 'Degrees: [4000, 4000, 4000]' in log
            if flags == ('-v', '3'):
                assert '64000000000' in saved

    log, saved = report('[4000000]*3', '-v', '3')
    for text in (log, saved):
        assert 'matrix size: 20000000000000' in text
        assert 'Core size: 12000000000000;' in text
        assert 'Formation budget: soft-O(sigma*rho^omega)' in text
        assert 'Step 4 = max(formation budget' in text

    if sys.maxsize > 2**32:
        log, saved = report('[4000000000]*3', '-v', '3')
        for text in (log, saved):
            assert 'matrix size: 20000000000000000000' in text
            assert 'Core size: 12000000000000000000;' in text
            assert 'Delta=64000000000000000000000000000.' in text

    log, saved = report('[3,9,9,9]', '-v', '3')
    for text in (log, saved):
        assert 'Core size: 181;' in text
        assert 'Step 3 generic candidate selection/rank verification' in text
        assert 'Delta=2187' in text

    def total(text):
        return re.search(r'(?:Overall (?:Dixon )?|Total )complexity[^\n]*: ([0-9.]+)', text)[1]

    basic, saved = report('[400]*3')
    normal, _ = report('[400]*3', '-v', '2')
    detailed, _ = report('[400]*3', '-v', '3')
    debug, _ = report('[400]*3', '-v', '4')
    quiet, _ = report('[400]*3', '-v', '0')
    assert len({total(text) for text in (basic, normal, detailed, debug, quiet)}) == 1

    body = basic.split('=== Complexity Analysis ===\n', 1)[1].split('\nUse -v 2', 1)[0].strip()
    lines = body.splitlines()
    assert len(lines) == 6, body
    prefixes = ['Dixon matrix size: ', 'Rank prediction: ', 'Schur core size: ',
                'Step 1 complexity (log2): ', 'Step 4 complexity (log2): ', 'Total complexity (log2): ']
    assert all(line.startswith(prefix) for line, prefix in zip(lines, prefixes)), body
    assert lines[:3] == ['Dixon matrix size: 239800', 'Rank prediction: 200000', 'Schur core size: 120000']
    assert body == saved.rsplit('\n\n', 1)[-1].strip()
    for text in (basic, saved):
        assert 'Comparison' not in text and 'Space Complexity' not in text
        assert 'Formation budget:' not in text
    for decoration in ('DRSolve v', 'Random system spec:', 'Degrees: [400, 400, 400]',
                       'Mode: Complexity analysis  |  Field: F_65537', 'Omega',
                       'Report saved to:', 'Memory:', 'Time:'):
        assert decoration in basic, (decoration, basic)
    assert basic.index('Use -v 2 or higher for a more detailed report.') < basic.index('Report saved to:')
    assert '--- Comparison:' in normal and 'Space Complexity' not in normal
    assert 'Formation budget:' not in normal
    assert 'Space Complexity' in detailed and 'Formation budget:' in detailed
    assert 'Step 1 ordinary dense interpolation details:' not in detailed
    assert 'Step 1 ordinary dense interpolation details:' in debug
    assert 'Use -v 2 or higher' not in normal + detailed + debug
    alias, _ = report('[400]*3', '--debug')
    override, _ = report('[400]*3', '--debug', '-v', '2')
    assert 'Formation budget:' in alias and 'Formation budget:' not in override

    disabled, _ = report('[400]*3', '--no-rank-pred')
    assert 'Rank prediction: disabled' in disabled
    assert 'Schur core size: disabled' in disabled
    assert total(disabled) != total(basic)
    no_schur, _ = report('[400]*3', '--no-step4-schur')
    assert 'Rank prediction: 200000' in no_schur and 'Schur core size: disabled' in no_schur

    for flags in ((), ('-v', '0'), ('-v', '3')):
        generic, saved = report('[400]*3', *flags, field=None)
        for text in (generic, saved):
            assert 'sufficiently large' in text
            assert 'field-size logarithms' in text
            assert 'Field characteristic p: -1' not in text
            assert 'Field order q: 0' not in text
        if flags != ('-v', '0'):
            if flags == ():
                assert 'Rank prediction: 200000' in generic and 'Schur core size: 120000' in generic
            else:
                assert 'rank prediction (uniform) matrix size: 200000' in generic
                assert 'Core size: 120000;' in generic
                assert re.search(r'large-field soft-O, log2\): [0-9.]+', generic)
            assert 'unavailable without a field size' not in generic
    generic_mq, _ = report('[2]*6', '-v', '3', field=None)
    assert 'extension degree: 1.' in generic_mq
    assert 'Step 1 MQ simplex interpolation (log2):' in generic_mq
    small_mq, _ = report('[2]*6', '-v', '3', field='2')
    assert 'extension degree: 3.' in small_mq
    extension, _ = report('[2]*3', '-v', '2', field='2^8')
    assert 'Field characteristic p: 2' in extension and 'Field order q: 256' in extension
    rational, _ = report('[400]*3', field='0')
    assert 'Field: Q' in rational and 'sufficiently large' not in rational

    output = Path(tmp) / 'inline.dr'
    result = subprocess.run([str(binary), '-c', 'x^2+y^2+1,x*y+t,x+y+t^2',
                             'x,y', '-o', str(output)],
                            cwd=tmp, env=env, text=True, capture_output=True, check=True, timeout=10)
    assert 'sufficiently large' in result.stdout and 'Schur core size: 3' in result.stdout
    assert total(result.stdout) == total(output.read_text())
    missing = subprocess.run([str(binary), '-r', '[2]*3'], cwd=tmp, env=env,
                             text=True, capture_output=True, timeout=10)
    assert missing.returncode != 0, 'Random generation still requires an explicit field'

print('Complexity summaries, optional large-field model, rank controls and detailed reports passed')

"""Direct large-prime Dixon reports its four stages and opt-in stage timings."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
with tempfile.TemporaryDirectory() as tmp:
    results = []
    for flags in ([], ['--time'], ['--time', '--silent']):
        output = Path(tmp) / 'result.dr'
        run = subprocess.run(
            [binary, '-r', '--seed', '12345', '[3]*3', '18446744073709551629',
             '--resultant-only', *flags, '-o', str(output)],
            capture_output=True, text=True, check=True, timeout=120)
        stages = re.findall(r'^Step ([1-4]): .*$', run.stdout, re.M)
        times = re.findall(r'^Step ([1-4]) time: (.*)$', run.stdout, re.M)
        assert stages == ([] if '--silent' in flags else ['1', '2', '3', '4']), run.stdout
        assert [step for step, _ in times] == (['1', '2', '3', '4'] if flags == ['--time'] else []), run.stdout
        for _, timing in times:
            assert re.fullmatch(r'[0-9.]+ seconds \| Memory: [0-9.]+ MiB', timing), timing
        results.append(output.read_text().split('Resultant:\n', 1)[1].strip())
    assert results[0] == results[1] == results[2]
print('Large-prime stage output and timing checks passed.')

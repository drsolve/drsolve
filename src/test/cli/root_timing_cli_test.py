"""Root timing follows --time and does not run under --no-roots."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
with tempfile.TemporaryDirectory() as tmp:
    for field in ('65537', '2^32', '0', '18446744073709551629'):
        for flags in ([], ['--time'], ['--time', '--no-roots'], ['--time', '--silent']):
            run = subprocess.run(
                [binary, 'x+y,y^2+y', 'x', field, *flags,
                 '-o', str(Path(tmp) / 'result.dr')],
                capture_output=True, text=True, check=True)
            lines = re.findall(r'^Root finding time:.*$', run.stdout, re.M)
            assert len(lines) == (flags == ['--time']), (field, flags, run.stdout)
            if lines:
                assert re.search(r'CPU time: [0-9.]+ seconds \| Wall time: [0-9.]+ seconds.*Memory: [0-9.]+ MiB', lines[0]), lines
print('Root timing CLI checks passed.')

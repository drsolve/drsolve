"""A dynamic CLI must load its matching library, independent of cwd/LD_LIBRARY_PATH."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve()
if not sys.platform.startswith('linux'):
    print('Linux dynamic loader regression: skipped on this platform')
    sys.exit(0)

with tempfile.TemporaryDirectory() as tmp:
    directory = Path(tmp)
    # If the loader selects this conflicting library, startup must fail.
    (directory / 'libdrsolve.so').write_bytes(b'conflicting library\n')
    env = dict(os.environ)
    env['LD_LIBRARY_PATH'] = tmp + ':' + env.get('LD_LIBRARY_PATH', '')
    linked = subprocess.check_output(['ldd', str(binary)], cwd=tmp, env=env, text=True)
    assert str(binary.parent / 'libdrsolve.so') in linked, linked
    subprocess.run([str(binary), '--help'], cwd=tmp, env=env,
                   capture_output=True, check=True, timeout=30)
    output = directory / 'result.dr'
    subprocess.run([str(binary), 'x+y,x-y', 'x', '65537',
                    '--resultant-only', '-o', str(output)],
                   cwd=tmp, env=env, capture_output=True, check=True, timeout=30)
    result = ''.join(output.read_text().split('Resultant:\n', 1)[1].split())
    assert result == 'y', result  # CLI normalizes the resultant to monic form.
print('Dynamic library selection and CLI computation passed')

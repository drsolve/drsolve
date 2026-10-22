"""Check quick/full help dispatch without running a computation."""
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1]).resolve())
# Makefile builds locate libdrsolve relative to the working directory.
env = os.environ.copy()
libpath = "DYLD_LIBRARY_PATH" if sys.platform == "darwin" else "LD_LIBRARY_PATH"
env[libpath] = str(Path(binary).parent) + os.pathsep + env.get(libpath, "")
with tempfile.TemporaryDirectory() as tmp:
    def help_text(*args):
        result = subprocess.run([binary, *args], cwd=tmp, env=env, text=True,
                                capture_output=True, check=True, timeout=10)
        assert not result.stderr, result.stderr
        return result.stdout

    short = help_text('-h')
    assert short == help_text('--help') == help_text()
    assert '--help-full' in short and 'OTHER OPTIONS:' in short
    assert 'ADVANCED AND EXPERIMENTAL:' not in short
    full = help_text('--help-full')
    assert 'FILE FORMAT:' in full and 'ADVANCED AND EXPERIMENTAL:' in full
    assert '--root-precision' in full and '--large-prime-reconstruct' in full
    sections = re.split(r'(?m)^[A-Z][A-Z AND]+:\n', full)[2:]
    assert len(sections) == 9
    for section in sections:
        assert '  Examples:\n' in section, section
        assert len(section.split('  Examples:\n', 1)[1].strip().splitlines()) >= 2
    assert '\nEXAMPLES:\n' not in full
    assert help_text('-v', '1', '--help') == short
    assert help_text('-v', '1', '--help-full') == full
    source = (Path(__file__).resolve().parents[3] / 'apps' / 'drsolve.c').read_text()
    for option in re.findall(r'\{"([\w-]+)"', source):
        assert re.search(r'--' + re.escape(option) + r'\b', full), option
    assert not list(Path(tmp).iterdir()), 'Help must not create output files'
print('Quick and full help passed')

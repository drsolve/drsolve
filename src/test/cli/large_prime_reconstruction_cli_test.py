"""The opt-in legacy large-prime path agrees with direct elimination and solving."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
prime = '18446744073709551629'
with tempfile.TemporaryDirectory() as tmp:
    def run(args, name):
        path = Path(tmp) / (name + '.dr')
        result = subprocess.run([binary, *args, '-o', str(path)],
                                capture_output=True, text=True, timeout=120, check=True)
        return result.stdout, path.read_text()

    system = 'x^2+y^2-5,x*y-2'
    outputs = []
    solutions = []
    for reconstruct in (False, True):
        flags = ['--large-prime-reconstruct', '--max-primes', '8'] if reconstruct else []
        stdout, saved = run(['--resultant-only', system, 'x', prime, *flags],
                            'resultant_' + str(reconstruct))
        assert ('Reconstructing over Q' in stdout) == reconstruct, stdout
        outputs.append(re.sub(r'\s+', '', saved.split('Resultant:\n', 1)[1]))
        stdout, _ = run(['-v', '2', '-s', '--solver', 'dixon', system, prime, *flags],
                        'solve_' + str(reconstruct))
        assert ('Reconstructing over Q' in stdout) == reconstruct, stdout
        solutions.append(sorted(re.findall(r'^  Solution \d+: (.*) \[pass\]$', stdout, re.M)))
    assert outputs[0] == outputs[1], outputs
    assert len(solutions[0]) == 4 and solutions[0] == solutions[1], solutions
    for field in ('257', '0'):
        invalid = subprocess.run([binary, '--large-prime-reconstruct', 'x+y,x-y', 'x', field],
                                 capture_output=True, text=True)
        assert invalid.returncode != 0 and 'requires a prime exceeding' in invalid.stderr
print('Large-prime reconstruction option checks passed.')

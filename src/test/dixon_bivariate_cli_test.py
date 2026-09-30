#!/usr/bin/env python3
"""Compare complete resultants across construction, compression and fallback."""
import os
from pathlib import Path
import subprocess
import tempfile

os.chdir(Path(__file__).resolve().parents[2])
cases = [('[3]*3', '65537', 0, []), ('[5]*3', '65537', 42, []),
         ('[7]*3', '65537', 9, []), ('[3]*3', '3', 2, []),
         ('[3]*3', '2', 5, []), ('[3]*3', '2^3', 7, []),
         ('[4]*3', '65537', 17, ['--density', '0.3']),
         ('[3,4,5]', '65537', 11, [])]
with tempfile.TemporaryDirectory() as directory:
    for degrees, field, seed, extra in cases:
        results = []
        for method, compression, native in [(0, False, True), (0, True, True),
                                             (5, False, False), (5, False, True),
                                             (5, True, True)]:
            path = Path(directory) / 'result.dr'
            env = dict(os.environ, DRSOLVE_FAST_NATIVE='1' if native else '0')
            args = ['./drsolve', '-r', degrees, field, '--seed', str(seed),
                    '--threads', '2', '--resultant-only', '-o', str(path), *extra]
            if method == 5:
                args += ['--method', '5']
            args += ['--mq-step4-schur' if compression else '--no-mq-step4-schur']
            subprocess.run(args, env=env, capture_output=True, text=True, check=True)
            results.append(path.read_text().split('Resultant:\n', 1)[1].split('\n', 1)[0])
        assert results[0] == results[1], (degrees, field, seed, 'ordinary compression')
        assert results[2] == results[3] == results[4], (degrees, field, seed, 'recursive arithmetic/compression')
        # Sparse/mixed systems may select different minors in the two existing
        # constructors, with different extraneous factors. Preserve each path.
        if not extra and degrees != '[3,4,5]':
            assert results[0] == results[2], (degrees, field, seed, 'construction')
        print(degrees, field, seed, 'complete resultants match within each construction')
print('Bivariate construction/compression/native/fallback CLI PASS')

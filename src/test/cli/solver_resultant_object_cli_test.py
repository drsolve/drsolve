"""Elimination passes its resultant object to one root-finding stage."""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

binary = str(Path(sys.argv[1] if len(sys.argv) > 1 else './drsolve').resolve())
with tempfile.TemporaryDirectory() as tmp:
    output = Path(tmp) / 'result.dr'

    def solve(system, field, expected, flags=(), root_stages=1):
        proc = subprocess.run([binary, '--solve-verbose', '--time', system, field,
                               *flags, '-o', str(output)], capture_output=True,
                              text=True, check=True, timeout=120)
        log = proc.stdout + proc.stderr
        # Legacy resultant interfaces print their own root report before the
        # solver parses and solves the same resultant again.
        assert 'Roots in ' not in log, log
        assert len(re.findall(r'^Root finding time:', log, re.M)) == root_stages, log
        actual = set()
        for line in re.findall(r'^  \[\d+\] (.*?) \[PASS\]$', log, re.M):
            actual.add(tuple(sorted(part.strip() for part in line.split(','))))
        assert actual == {tuple(sorted(point)) for point in expected}, log
        assert '[FAIL]' not in log, log
        return log

    # Prime, Zech, native binary and general extension-field backends.
    for field, prime in [('257', 257), ('9', 3), ('256', 2), ('65537^2', 65537)]:
        solve('x+y,y^2-y', field, [('x = 0', 'y = 0'),
                                  (f'x = {prime-1}', 'y = 1')])
    solve('x+y,(y-1)^3', '257', [('x = 256', 'y = 1')])
    solve('x+y,y^2+1', '3', [])
    solve('x+y,x+y+1', '257', [])
    log = solve('x+y,2*x+2*y', '257', [], root_stages=0)
    assert 'dimension' in log.lower(), log

    # A zero first combination must not overwrite a later nonzero resultant;
    # subsequent combinations must not replace the selected object's roots.
    solve('x+y,2*x+2*y,y^2-y', '257',
          [('x = 0', 'y = 0'), ('x = 256', 'y = 1')], ['--solver', 'dixon'])
    solve('x+y,2*x+2*y,y^2-y,y-1', '257',
          [('x = 256', 'y = 1')], ['--solver', 'dixon'])
    solve('x+y+z,y-z,z^2-z', '257',
          [('x = 0', 'y = 0', 'z = 0'), ('x = 255', 'y = 1', 'z = 1')],
          ['--solver', 'dixon'], root_stages=3)

print('Solver object resultants, single root stage and solution sets passed')

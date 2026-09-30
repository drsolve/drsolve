#!/usr/bin/env python3
"""Compare complete resultants across construction, compression and fallback."""
import os
from pathlib import Path
import subprocess
import re
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
                                             (5, True, True), (None, True, True),
                                             (None, True, True), (5, False, True),
                                             (0, True, True)]:
            path = Path(directory) / 'result.dr'
            env = dict(os.environ, DRSOLVE_FAST_NATIVE='1' if native else '0')
            args = ['./drsolve', '-r', degrees, field, '--seed', str(seed),
                    '--threads', '2', '--resultant-only', '-o', str(path), *extra]
            if method == 5:
                args += ['--method', '5']
            elif method == 0:
                args += ['--dixon']
            args += ['--mq-step4-schur' if compression else '--no-mq-step4-schur']
            if len(results) >= 6 and field not in ('2', '3'):
                args += ['--fq-det-method', 'interp']
            run = subprocess.run(args, env=env, capture_output=True, text=True, check=True)
            if method is None:
                assert 'Recursive block construction builds the Dixon matrix directly.' in run.stdout
            results.append(path.read_text().split('Resultant:\n', 1)[1].split('\n', 1)[0])
        assert results[0] == results[1], (degrees, field, seed, 'ordinary compression')
        assert results[2] == results[3] == results[4], (degrees, field, seed, 'recursive arithmetic/compression')
        assert results[6] == results[7] == results[4], (degrees, field, seed, 'recursive interpolation')
        assert results[8] == results[0], (degrees, field, seed, 'ordinary interpolation')
        assert results[5] == results[4], (degrees, field, seed, 'automatic constructor')
        # Sparse/mixed systems may select different minors in the two existing
        # constructors, with different extraneous factors. Preserve each path.
        if not extra and degrees != '[3,4,5]':
            assert results[0] == results[2], (degrees, field, seed, 'construction')
        print(degrees, field, seed, 'complete resultants match within each construction')
    # Small-field interpolants retain base-field roots. CLI output is monic,
    # so coefficients/values can differ by a nonzero scalar after reduction.
    for prime in [2,3,5,13]:
        polynomials=[]
        for backend in ['hnf','interp']:
            out=Path(directory)/'small.dr'
            run=subprocess.run(['./drsolve','-r','[3]*3',str(prime),'--seed','2',
                                '--threads','2','--fq-det-method',backend,'-v','2',
                                '--resultant-only','-o',str(out)],capture_output=True,text=True,check=True)
            polynomials.append(out.read_text().split('Resultant:\n',1)[1].split('\n',1)[0])
            if backend=='interp':
                assert f'points={prime},' in run.stdout and 'field-value polynomial' in run.stdout
        for x in range(prime):
            zeros=[]
            for poly in polynomials:
                env={name:x for name in re.findall(r'\b[A-Za-z_]\w*',poly)}
                zeros.append(eval(poly.replace('^','**'),{'__builtins__':{}},env)%prime==0)
            assert zeros[0]==zeros[1],(prime,x,polynomials)
    for threads,level,expected in [(2,2,False),(4,1,False),(4,2,True)]:
        run=subprocess.run(['./drsolve','-r','[3]*3','65537','--seed','2',
                            '--threads',str(threads),'-v',str(level),'--resultant-only',
                            '-o',str(Path(directory)/'automatic.dr')],capture_output=True,text=True,check=True)
        assert ('Step 4 interpolation:' in run.stdout)==expected,(threads,level)
    # Explicit constructor/step options win; shape alone must not match a
    # three-equation system with zero or two retained variables.
    for extra, recursive in [([], True), (['--dixon'], False),
                             (['--method', '5'], True), (['--method', '0'], False),
                             (['--step1', '0'], False), (['--step4', '1'], False),
                             (['-n', '4'], False), (['-n', '2'], False)]:
        args = ['./drsolve', '-r', '[1]*3', '257', '--seed', '42',
                '--resultant-only', '-o', str(Path(directory)/'selection.dr'), *extra]
        run = subprocess.run(args, capture_output=True, text=True, check=True)
        assert ('Recursive block construction builds the Dixon matrix directly.' in run.stdout) == recursive, extra
    file = Path(directory)/'explicit.dr'
    file.write_text('# options: --method 0\nx,y\n257\nx+y+z+1, x+2*y+3*z+2, 2*x+4*y+z+3\n')
    run = subprocess.run(['./drsolve', '-f', str(file), '--resultant-only',
                          '-o', str(Path(directory)/'file.dr')], capture_output=True, text=True, check=True)
    assert 'Recursive block construction builds the Dixon matrix directly.' not in run.stdout
print('Bivariate construction/compression/native/fallback CLI PASS')

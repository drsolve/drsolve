"""Independent oracle: sage -python src/test/cli/quotient_sift_sage_test.py ./drsolve."""
from sage.all import GF, PolynomialRing, QQ
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

binary = Path(sys.argv[1]).resolve()
fixtures = Path(__file__).resolve().parents[3] / 'examples' / 'sift'
cases = ['two_parameters', 'example34', 'example35',
         'example36_b5_slice', 'example36_a1_slice']
with tempfile.TemporaryDirectory(prefix='sift-sage-') as work:
    for name in cases:
        source = fixtures / (name + '.dr')
        lines = [s.strip() for s in source.read_text().splitlines()
                 if s.strip() and not s.startswith('#')]
        elim = lines[0].split(',')
        parameters = ['a', 'b'] if name == 'two_parameters' else [
            'a' if name == 'example34' else 'z' if name == 'example35' else 'T']
        R = PolynomialRing(GF(int(lines[1])), names=elim + parameters, order='lex')
        I = R.ideal([R(s) for s in lines[2:]])
        E = I.elimination_ideal(list(R.gens()[:len(elim)]))
        oracle = list(E.groebner_basis())
        assert len(oracle) == 1
        outputs = []
        for flag in ['--no-sift', '--sift']:
            output = Path(work) / (name + flag + '.dr')
            proc = subprocess.run([str(binary), str(source), flag, '--resultant-only',
                                   '-o', str(output)], text=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  env=dict(os.environ, OMP_NUM_THREADS='4'), timeout=60)
            assert proc.returncode == 0, proc.stdout
            saved = output.read_text()
            f = R(re.search(r'\nResultant:\n([^\n]+)', saved)[1])
            outputs.append(f)
            if flag == '--sift':
                removed = R(re.search(r'certified removed multiplier = ([^\n]+)', saved)[1])
                assert outputs[0] == f * removed
                assert f in I
                assert f * f.lc()**(-1) == oracle[0] * oracle[0].lc()**(-1), (name, f, oracle)
        print(name + ': quotient-sift output equals Sage elimination generator', flush=True)

# The printed first specialization of Example 36 has a valid branch omitted
# by the paper's quartic: all four equations vanish at (-5,0,0,-90).
R = PolynomialRing(QQ, names=('x1', 'x2', 'x3', 'T'))
lines = [s.strip() for s in (fixtures / 'example36_b5_slice.dr').read_text().splitlines()
         if s.strip() and not s.startswith('#')]
assert all(R(s)(-5, 0, 0, -90) == 0 for s in lines[2:])
T = PolynomialRing(QQ, 'T').gen()
paper_g = 12397*T**4+502080*T**3-95020200*T**2+1960028928*T+9550921680
assert paper_g(-90) != 0
assert GF(65537)(paper_g(-90)) != 0
print('Example 36 paper quartic omits a verified original-system solution at T=-90')

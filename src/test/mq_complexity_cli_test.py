#!/usr/bin/env python3
"""Detailed estimates minimize complete stage costs; verbosity one stays legacy."""
import math
import re
import subprocess
from pathlib import Path

root = Path(__file__).resolve().parents[2]

def report(n, q, verbosity=2, degree=2):
    return subprocess.run([str(root/'drsolve'), '-c', '-v', str(verbosity), '-r',
                           '-n', str(n), f'[{degree}]*{n}', str(q)], cwd=root,
                          text=True, capture_output=True, check=True).stdout

def number(text, label):
    return float(re.search(re.escape(label)+r'.*?\(log2(?:: |\): )([\d.]+)', text)[1])

for n, q in [(4,65537),(10,65537),(32,65537),(6,2)]:
    text=report(n,q)
    step1=number(text,'Best Step 1 estimate:')
    step4=number(text,'Best Step 4 estimate:')
    simplex=number(text,'Step 1 MQ simplex interpolation')
    schur=number(text,'Step 4 MQ blocked Schur + core determinant')
    layered=number(text,'Step 1 MQ cached Laplace layered total-degree surrogate')
    assert step1<=min(simplex,layered)+1e-6
    assert step4<=schur+1e-6
    total=number(text,'Overall complexity =')
    assert abs(total-max(step1,step4))<1e-6
    if n==32: assert 'Best Step 1 estimate: MQ total-degree simplex interpolation' in text
    if n==10: assert 'Best Step 4 estimate: MQ blocked Schur + core determinant' in text
    if q==2: assert 'extension degree: 3.' in text
    old=report(n,q,1)
    assert 'MQ simplex interpolation' not in old and 'MQ blocked Schur' not in old
assert 'Step 1 MQ simplex interpolation' not in report(4,65537,degree=3)
print('MQ detailed minimum, complete Schur cost, extension cost, and verbosity checks passed')

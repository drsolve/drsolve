#!/usr/bin/env python3
"""Native exact regressions; no Sage dependency. Tiny fields are enumerated."""
import argparse
import itertools
import json
from pathlib import Path
import random
import subprocess
import tempfile


def exps(n,d):
    return [e for e in itertools.product(range(d+1),repeat=n) if sum(e)<=d]


def evaluate(f,a,p):
    value=0
    for e,c in f.items():
        term=c
        for x,k in zip(a,e): term=term*pow(x,k,p)%p
        value=(value+term)%p
    return value


def planted(n,p,seed,points=1):
    rng=random.Random(seed); a=tuple(rng.randrange(p) for _ in range(n)); fs=[]
    b=((a[0]+1)%p,)+a[1:]
    for _ in range(n+1):
        f={e:rng.randrange(p) for e in exps(n,2) if sum(e)}
        if points==2:
            e=(1,)+(0,)*(n-1)
            f[e]=(f[e]-(evaluate(f,b,p)-evaluate(f,a,p))*pow((b[0]-a[0])%p,-1,p))%p
        f[(0,)*n]=-evaluate(f,a,p)%p
        fs.append(f)
    return fs


def text_poly(f):
    terms=[]
    for e,c in f.items():
        if c: terms.append(str(c)+''.join('*x%d^%d'%(i,k) for i,k in enumerate(e) if k))
    return '+'.join(terms) or '0'


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--binary',default='./build/dixon_closure_bench'); args=ap.parse_args()
    cases=[]
    for p in (2,3,7,257,65537):
        for seed in range(3):
            for points in (1,2): cases.append((3,p,planted(3,p,seed,points)))
    cases += [(2,7,[{(1,0):1},{(0,1):1},{(0,0):1}]),
              (3,7,[{(1,0,0):1},{(0,1,0):1},{(0,0,1):1},{}]),
              (3,7,[{(2,0,0):1}]*4)]
    total=0
    with tempfile.TemporaryDirectory() as directory:
        path=Path(directory)/'input.txt'
        for index,(n,p,fs) in enumerate(cases):
            path.write_text('%d %d %d\n'%(p,n,len(fs))+'\n'.join(map(text_poly,fs))+'\n')
            expected=None
            if p<=7: expected={a for a in itertools.product(range(p),repeat=n) if all(evaluate(f,a,p)==0 for f in fs)}
            for degree in (3,4):
                results=[]
                methods=['seed','hybrid3','macaulay','augment3']
                if degree==4: methods+=['hybrid4','augment4']
                for method in methods:
                    run=subprocess.run([args.binary,str(path),method,str(degree),'1'],capture_output=True,text=True,check=True)
                    r=json.loads(run.stdout); results.append(r); total+=1
                    if expected is not None:
                        if r['recovered'] is not None: assert expected=={tuple(r['recovered'])},(index,method,r,expected)
                        if r['inconsistent']: assert not expected,(index,method,r,expected)
                    if r['recovered'] is not None: assert all(evaluate(f,r['recovered'],p)==0 for f in fs)
                assert len({r['rowspace_hash'] for r in results})==1,(index,degree,results)
    print('PASS:',total,'native closure runs; exact truncated/full Dixon equality, Macaulay membership, identical final row spaces, and small-field enumeration')


if __name__=='__main__': main()

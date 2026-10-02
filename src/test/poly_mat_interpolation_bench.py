#!/usr/bin/env python3
"""Compare scalar and degree-grouped batch evaluation in determinant interpolation."""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--degrees',type=int,nargs='+',default=[12,16])
parser.add_argument('--primes',type=int,nargs='+',default=[65537,1000003,998244353])
parser.add_argument('--seed',type=int,default=12345)
parser.add_argument('--threads',type=int,default=1)
parser.add_argument('--repeats',type=int,default=3)
parser.add_argument('--json',type=Path)
args=parser.parse_args()
if args.repeats<1 or args.threads<1: parser.error('repeats and threads must be positive')
root=Path(__file__).resolve().parents[2]
records=[]
with tempfile.TemporaryDirectory(prefix='drsolve-batch-eval-') as directory:
    for degree in args.degrees:
        for prime in args.primes:
            reference=None
            for repeat in range(args.repeats):
                for mode in (['scalar','auto'] if repeat%2==0 else ['auto','scalar']):
                    file=Path(directory)/'result.dr'
                    run=subprocess.run([str(root/'drsolve'),'-r',f'[{degree}]*3',str(prime),
                        '--seed',str(args.seed),'--threads',str(args.threads),'-v','2',
                        '--fq-det-method','interp','-o',str(file)],cwd=root,
                        env=dict(os.environ,DRSOLVE_INTERP_EVAL=mode),capture_output=True,text=True,check=True)
                    result=file.read_text().split('Resultant:\n',1)[1].split('\n',1)[0]
                    if reference is None: reference=result
                    if result!=reference: raise RuntimeError(f'Result mismatch: {degree=} {prime=} {mode=}')
                    seconds=float(re.search(r'Step 4 time: ([0-9.]+)',run.stdout)[1])
                    detail=re.search(r'matrix evaluation=([0-9.]+)s, determinants=([0-9.]+)s',run.stdout)
                    records.append(dict(degree=degree,prime=prime,seed=args.seed,threads=args.threads,
                        mode=mode,repeat=repeat,step4=seconds,evaluation_work=float(detail[1]),
                        determinant_work=float(detail[2])))
                    print(f'd={degree} p={prime} {mode} run={repeat+1}: {seconds:.3f}s; full result matches',flush=True)
            for mode in ['scalar','auto']:
                times=[r['step4'] for r in records if r['degree']==degree and r['prime']==prime and r['mode']==mode]
                print(f'MEDIAN d={degree} p={prime} {mode}: {statistics.median(times):.3f}s',flush=True)
if args.json: args.json.write_text(json.dumps(records,indent=2)+'\n')

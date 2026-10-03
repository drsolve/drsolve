#!/usr/bin/env python3
"""Native-only timing pass after Sage checks, in paired shuffled order."""
import argparse
import json
import os
from pathlib import Path
import random
import statistics
import subprocess


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('directory',nargs='?',type=Path,default=Path('src/test/data/dixon_closure'))
    ap.add_argument('--binary',default='./build/dixon_closure_bench')
    ap.add_argument('--rounds',type=int,default=3); ap.add_argument('--samples',type=int,default=5)
    args=ap.parse_args(); output=[]
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    timing_keys=('seconds','setup_seconds','closure_seconds','dixon_seconds',
                 'cpu_seconds','cpu_setup_seconds','cpu_closure_seconds','cpu_dixon_seconds')
    compare_keys=('rank_history','nullity','recovered','rowspace_hash','submitted_rows','rref_cells','dense_matmul_multiply_adds')
    for path in sorted(args.directory.glob('n[0-9]*.json')):
        if 'boundary' in path.stem: continue  # Correctness/degree-limit trials, not timing claims.
        records=json.loads(path.read_text()); seeds=sorted({r['seed'] for r in records})
        methods=sorted({r['method'] for r in records})
        print('START',path.stem,flush=True)
        for seed in seeds:
            for repeat in range(args.rounds):
                order=methods[:]; random.Random(100*seed+repeat+937).shuffle(order)
                for method in order:
                    reference=next(r for r in records if r['seed']==seed and r['method']==method)
                    command=[args.binary,str(args.directory/'inputs'/path.stem/('seed%d.txt'%seed)),
                             method,str(reference['degree']),'0',str(args.samples+1)]
                    runs=[json.loads(line) for line in subprocess.check_output(command,env=env,text=True).splitlines()][1:]
                    for run in runs:
                        for key in compare_keys: assert run[key]==reference[key],(path.stem,seed,method,key)
                    r=dict(runs[-1]); r.update(dataset=path.stem,seed=seed,repeat=repeat)
                    r['timing_samples']=[{k:run[k] for k in timing_keys} for run in runs]
                    for key in timing_keys: r[key]=statistics.median(run[key] for run in runs)
                    output.append(r)
            (args.directory/'native_paired.json').write_text(json.dumps(output,indent=2)+'\n')
        print('DONE',path.stem,flush=True)


if __name__=='__main__': main()

#!/usr/bin/env python3
"""Exact CLI comparisons for the simplex path, its cache, and fallbacks."""
from pathlib import Path
import subprocess
import tempfile
import re

root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='mq-simplex-cli-') as tmp:
    cases=[(4,65537,12345,[]),(5,65537,12345,[]),(6,65537,12345,[]),
           (4,65537,12345,['--no-mq-step1-filter']),
           (4,65537,12345,['--step1','1']),
           (4,2,12345,[]),(5,101,0,[]),(5,101,1,[]),(5,101,2,[])]
    for n,q,seed,extra in cases:
        outputs=[]
        for simplex,threads in [(False,1),(True,1),(True,4)]:
            out=Path(tmp)/f'{n}-{q}-{seed}-{int(simplex)}-{threads}.dr'
            cmd=[str(root/'drsolve'),'--mq-step4-schur','--threads',str(threads),'-v','2',
                 '-r','--seed',str(seed),'-n',str(n),'-o',str(out),*extra,f'[2]*{n}',str(q)]
            if simplex:cmd.insert(1,'--mq-step1-simplex')
            result=subprocess.run(cmd,cwd=root,text=True,capture_output=True,check=True)
            outputs.append('\n'.join(l for l in out.read_text().splitlines() if not l.startswith('Time: ')))
            if simplex:
                count=result.stdout.count('MQ Step 1 simplex: points=')
                if q>n+1 and '--step1' not in extra:
                    assert count==1,(n,q,extra,result.stdout)
                    assert f'threads={threads},' in result.stdout
                    # Candidate failure/repair must reuse the full polynomial.
                    assert 'recomputing full Dixon polynomial' not in result.stdout
                else:assert 'MQ Step 1 simplex: fallback' in result.stdout
        assert outputs[0]==outputs[1]==outputs[2],(n,q,seed,extra)
    result=subprocess.run([str(root/'drsolve'),'--mq-step1-simplex','--no-mq-step1-simplex',
        '-v','2','-r','--seed','12345','-n','4','[2]*4','65537'],cwd=root,text=True,capture_output=True,check=True)
    assert 'MQ Step 1 simplex:' not in result.stdout
print('MQ simplex CLI: full outputs, one/four threads, cache reuse, fallback, explicit opt-out PASS')

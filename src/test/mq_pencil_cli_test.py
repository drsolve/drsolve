#!/usr/bin/env python3
"""Exact CLI comparisons for the pencil path, its cache, and fallbacks."""
from pathlib import Path
import subprocess
import tempfile
import re

root=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='mq-pencil-cli-') as tmp:
    cases=[(4,65537,12345,[]),(5,65537,12345,[]),(6,65537,12345,[]),
           (4,65537,12345,['--no-mq-step1-filter']),
           (4,65537,12345,['--step1','1']),
           (4,3,12345,[]),(5,7,7,[]),(5,101,0,[]),(5,101,1,[]),(5,101,2,[])]
    for n,q,seed,extra in cases:
        outputs=[]
        for pencil,threads in [(False,1),(True,1),(True,4)]:
            out=Path(tmp)/f'{n}-{q}-{seed}-{int(pencil)}-{threads}.dr'
            cmd=[str(root/'drsolve'),'--mq-step4-schur','--threads',str(threads),'-v','2',
                 '-r','--seed',str(seed),'-n',str(n),'-o',str(out),*extra,f'[2]*{n}',str(q)]
            if pencil:cmd.insert(1,'--mq-step1-pencil')
            result=subprocess.run(cmd,cwd=root,text=True,capture_output=True,check=True)
            outputs.append('\n'.join(l for l in out.read_text().splitlines() if not l.startswith('Time: ')))
            if pencil:
                count=result.stdout.count('MQ Step 1 pencil: size=')
                if q>n-1 and '--step1' not in extra:
                    assert count==1,(n,q,extra,result.stdout)
                    assert f'threads={threads},' in result.stdout
                    if '--no-mq-step1-filter' not in extra:
                        assert '(closure projected)' in result.stdout
                    if (n,q,seed)==(5,7,7):
                        assert 'MQ Step 1 Schur repair verified' in result.stdout
                else:assert 'MQ Step 1 pencil: fallback' in result.stdout
        assert outputs[0]==outputs[1]==outputs[2],(n,q,seed,extra)
    result=subprocess.run([str(root/'drsolve'),'--mq-step1-pencil','--no-mq-step1-pencil',
        '-v','2','-r','--seed','12345','-n','4','[2]*4','65537'],cwd=root,text=True,capture_output=True,check=True)
    assert 'MQ Step 1 pencil:' not in result.stdout
# Last explicit experimental backend wins; neither is enabled by default.
for flags,active,inactive in [(['--mq-step1-simplex','--mq-step1-pencil'],'pencil','simplex'),
                              (['--mq-step1-pencil','--mq-step1-simplex'],'simplex','pencil')]:
    p=subprocess.run([str(root/'drsolve'),*flags,'-v','2','--threads','1','-r',
        '--seed','12345','-n','4','[2]*4','65537'],cwd=root,text=True,capture_output=True,check=True)
    assert f'MQ Step 1 {active}:' in p.stdout and f'MQ Step 1 {inactive}:' not in p.stdout
# An actual MQ input whose linear rows have no t coefficient must fall back.
polys='x^2+y^2+t^2+1,x^2+x*y+t^2+2,y^2+x+t^2+3'
logs=[]
for flags in ([],['--mq-step1-pencil']):
    p=subprocess.run([str(root/'drsolve'),*flags,'-v','2','--threads','1',polys,'x,y','101'],
        cwd=root,text=True,capture_output=True,check=True)
    logs.append(p.stdout)
assert 'parameter coefficient matrix is rank deficient' in logs[1]
a=[re.findall(r'^  Final Resultant = (.*)$',s,re.M) for s in logs]
assert a[0] and a[0]==a[1]
print('MQ pencil CLI: full outputs, one/four threads, closure projection, fallback, explicit opt-out PASS')

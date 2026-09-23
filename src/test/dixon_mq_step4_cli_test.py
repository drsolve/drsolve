#!/usr/bin/env python3
"""End-to-end opt-in/opt-out and repaired-candidate fallback regressions."""
from pathlib import Path
import os
os.chdir(Path(__file__).resolve().parents[2])
import subprocess,re
cases=[(4,65537,12345,[]),(5,65537,12345,[]),(6,65537,12345,[]),
       (5,65537,12345,['--no-mq-step1-filter']),
       (4,65537,12345,['--fq-det-method','iter']),
       (4,65537,12345,['--step4','1']),
       (4,65537,12345,['--step4','4']),
       *[(4,3,s,[]) for s in range(5)],*[(5,101,s,[]) for s in range(5)]]
for n,q,seed,extra in cases:
    base=['./drsolve','-v','2','--threads','1','-r','--seed',str(seed),'-n',str(n),*extra,f'[2]*{n}',str(q)]
    logs=[]
    for option in [['--no-mq-step4-schur'],[],['--mq-step4-schur']]:
        p=subprocess.run(base[:1]+option+base[1:],capture_output=True,text=True,check=True)
        logs.append(p.stdout)
    results=[re.findall(r'^  Final Resultant = (.*)$',s,re.M) for s in logs]
    assert results[0] and results[0]==results[1]==results[2],(n,q,seed,extra,results)
    status='; '.join(l.strip() for l in logs[1].splitlines() if 'MQ Step 4 Schur:' in l)
    assert 'MQ Step 4 Schur:' not in logs[0]
    assert ('MQ Step 4 Schur:' in logs[1]) == ('MQ Step 4 Schur:' in logs[2])
    if n==4 and q==65537 and not extra:
        assert 'compression ' in status, logs[1]
    repaired='Step 1 Schur repair verified' in logs[1]
    print(n,q,seed,extra,'repair',repaired,status or 'explicit backend retained')
# Explicit opt-out wins; no environment flag is involved.
p=subprocess.run(['./drsolve','--mq-step4-schur','--no-mq-step4-schur','-v','2','--threads','1','-r','--seed','12345','-n','4','[2]*4','65537'],capture_output=True,text=True,check=True)
assert 'MQ Step 4 Schur:' not in p.stdout
print('CLI default/on/off exact resultant comparisons and explicit opt-out PASS')

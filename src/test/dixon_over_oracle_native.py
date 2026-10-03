"""Measure optimistic selected rows in C; construction/search stay excluded.

Also measure this shared generator's cache setup cost. This is a diagnostic,
not an implementable end-to-end speedup and not a universal complexity bound.
"""
from pathlib import Path
import json,random,statistics,subprocess

root=Path(__file__).resolve().parents[2];dest=root/'src/test/data/dixon_over_routes'
binary=root/'build/dixon_over_native';records=[]
for label in ('oracle5','oracle6'):
    for case in json.loads((dest/(label+'.json')).read_text()):
        if case.get('timeout'):continue
        assert case['best']['inject_after']==0,'native diagnostic implements initial injection only'
        original=Path(case['input_file']);lines=original.read_text().splitlines()
        p,n,m=map(int,lines[0].split());extra=case['extra_polynomials']
        exported=dest/'inputs'/label/('selected_seed%d_m%d.txt'%(case['seed'],m))
        exported.write_text('%d %d %d\n'%(p,n,m+len(extra))+'\n'.join(lines[1:]+extra)+'\n')
        for repeat in range(3):
            variants=[('seed',original),('oracle_seed',exported),('shared',original)]
            random.Random(case['seed']*100+repeat).shuffle(variants)
            for method,path in variants:
                cmd=[str(binary),str(path),method,str(case['degree']),'1','6']
                runs=[json.loads(s) for s in subprocess.check_output(cmd,text=True).splitlines()][1:]
                r=dict(runs[-1])
                for k in list(r):
                    if k.endswith('seconds'):r[k]=statistics.median(x[k] for x in runs)
                if method!='shared':assert r['recovered']==case['baseline']['recovered']
                r.update(n=n,m=m,seed=case['seed'],repeat=repeat,selected=len(extra),diagnostic_only=True)
                records.append(r)
        (dest/'oracle_native.json').write_text(json.dumps(records,indent=2)+'\n')
        print('DONE oracle native',n,m,case['seed'],len(extra),flush=True)

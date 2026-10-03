"""Time one production square-core Dixon resultant on the same dense inputs.

This is deliberately a different output task: one subset resultant, before
extra-equation filtering or full-system solving. It is not a full-solver ratio.
"""
from sage.all import GF,PolynomialRing
from dixon_over_routes import solve_gb
from pathlib import Path
import json,os,re,subprocess,tempfile,time

root=Path(__file__).resolve().parents[2];dest=root/'src/test/data/dixon_over_routes';records=[]
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
with tempfile.TemporaryDirectory(prefix='dixon-over-production-') as directory:
    for n in (5,6,7):
        rows=json.loads((dest/('dense%d.json'%n)).read_text())
        for seed in range(3):
            ref=next(r for r in rows if r['seed']==seed and r['m']==n+1 and r['method']=='gb')
            inp=Path(ref['input_file']).read_text().splitlines();p=int(inp[0].split()[0]);fs=inp[1:1+n]
            R=PolynomialRing(GF(p),'x%d'%(n-1));a=R.base_ring()(ref['recovered'][-1])
            P=PolynomialRing(GF(p),names=tuple('x%d'%j for j in range(n)))
            core=solve_gb([P(f) for f in fs]);core_eliminant=R(core['eliminant'])
            cli_input=Path(directory)/'input.dr'
            cli_input.write_text(','.join('x%d'%j for j in range(n-1))+'\n'+str(p)+'\n'+'\n'.join(fs)+'\n')
            for repeat in range(4):
                out=Path(directory)/'result.dr'
                cmd=[str(root/'drsolve'),'--resultant-only','--threads','1','--time','-v','2','-o',str(out),'-f',str(cli_input)]
                start=time.perf_counter()
                try:run=subprocess.run(cmd,cwd=root,env=env,check=True,capture_output=True,text=True,timeout=30)
                except subprocess.TimeoutExpired:
                    records.append(dict(n=n,seed=seed,timeout=True,seconds=time.perf_counter()-start));break
                elapsed=time.perf_counter()-start
                text=out.read_text();poly=R(text.split('Resultant:\n',1)[1].split('\n',1)[0])
                assert poly and poly(a)==0 and poly % core_eliminant == 0
                steps={int(k):float(t) for k,t in re.findall(r'Step ([1-4]) time: ([0-9.]+)',run.stdout)}
                steps.update({int(k):float(t) for k,t in re.findall(r'Step ([1-4]) time: CPU time: [0-9.]+ seconds \| Wall time: ([0-9.]+)',run.stdout)})
                total=re.search(r'Total - CPU time: ([0-9.]+) seconds \| Wall time: ([0-9.]+)',run.stdout)
                if repeat:
                    records.append(dict(n=n,seed=seed,repeat=repeat,seconds=elapsed,step_seconds=steps,
                                        solver_cpu_seconds=float(total[1]),solver_wall_seconds=float(total[2]),
                                        resultant_degree=int(poly.degree()),planted_coordinate_is_root=True,
                                        core_eliminant_divides=True))
                if repeat==1:(dest/('production_n%d_seed%d.log'%(n,seed))).write_text(run.stdout)
            (dest/'production.json').write_text(json.dumps(records,indent=2)+'\n')
            print('DONE production',n,seed,flush=True)

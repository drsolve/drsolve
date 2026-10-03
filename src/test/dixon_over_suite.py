"""Run the Sage route sweep serially; concurrent timing jobs are prohibited."""
import json,os,subprocess,sys,time
from pathlib import Path

root=Path(__file__).resolve().parents[2]
dest=root/'src/test/data/dixon_over_routes';dest.mkdir(parents=True,exist_ok=True)
all_methods='scalar,quotient,hybrid4,hybrid_support4,module,module_units,module_dixon,module_unit_dixon,module_stream,module_nf,core_full,core_stream,gb,shared,linear2'
configs=[
    ('dense5',5,4,'dense',3,3,'1,3',all_methods,30),
    ('dense6',6,5,'dense',3,3,'1,3',all_methods+',hybrid5',30),
    ('dense7',7,5,'dense',3,1,'1,3','scalar,quotient,hybrid4,hybrid5,module,module_units,module_nf,core_full,core_stream,gb,shared,linear2',30),
    ('two5',5,4,'two',3,3,'1,3',all_methods,30),
    ('partial5',5,4,'partial',3,3,'1,3',all_methods,30),
    ('structured6',6,5,'structured',3,3,'1,3','scalar,quotient,hybrid4,module,module_units,module_nf,core_full,core_stream,gb,shared,linear2',30),
    ('structured7',7,5,'structured',3,3,'1,3','quotient,hybrid4,module_nf,core_full,core_stream,gb,shared',30),
    ('redundant5',5,4,'redundant',3,3,'1,3','scalar,quotient,module_nf,core_full,core_stream,gb,shared',30),
    ('square5',5,6,'dense',3,2,'0','quotient,core_full,core_stream,gb',30),
    ('shared5',5,4,'dense',3,3,'1,3,7,15','shared',30),
    ('shared6',6,4,'dense',3,3,'1,3,7,15','shared',30),
    ('oracle5',5,4,'dense',3,1,'1,3','oracle',60),
    ('oracle6',6,5,'dense',3,1,'1,3','oracle',90),
    ('nf_full5',5,4,'dense',3,2,'1,3','module_low,module_nf_full,module_nf2_full',30),
    ('nf_full6',6,5,'dense',3,1,'1,3','module_low,module_nf5,module_nf_full,module_nf2_full,module3',30),
    ('nf_full7',7,5,'dense',2,1,'1,3','module_nf_full,module_nf2_full,module3',30),
]
env=dict(os.environ,DOT_SAGE='/tmp/drsolve-closure-sage',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
manifest=[]
for label,n,d,fam,seeds,repeats,extras,methods,timeout in configs:
    if len(sys.argv)>1 and label not in sys.argv[1:]:continue
    cmd=['sage','-python','src/test/dixon_over_routes_bench.py','--variables',str(n),'--degree',str(d),
         '--family',fam,'--seeds',str(seeds),'--repeats',str(repeats),'--extras',extras,'--methods',methods,
         '--timeout',str(timeout),'--output',str(dest/(label+'.json'))]
    print('START',label,flush=True);start=time.monotonic()
    with (dest/(label+'.log')).open('w') as log:
        proc=subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT)
    record=dict(label=label,command=cmd,returncode=proc.returncode,seconds=time.monotonic()-start)
    manifest.append(record);(dest/'suite_runs.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('DONE',record,flush=True)
    if proc.returncode:raise SystemExit(proc.returncode)

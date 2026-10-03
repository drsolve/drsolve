"""Paired, shuffled benchmarks and independent exact audits for over routes."""
from sage.all import GF, PolynomialRing, matrix
from dixon_closure_cases import make_system
from dixon_over_routes import *
from pathlib import Path
import argparse, hashlib, json, random, signal, statistics, sys, traceback


def inputs(n,m,seed,family,prime=257):
    P=PolynomialRing(GF(prime),names=tuple('x%d'%i for i in range(n)),order='degrevlex')
    fs,points=make_system(P,seed,max(m,n+1),'partial' if family=='partial' else 'dense',2 if family=='two' else 1)
    if family=='structured':
        rng=random.Random(104729+seed);xs=P.gens();a=points[0]
        for i in range(n):
            f=xs[i]**2+sum(P.base_ring()(rng.randrange(prime))*x for x in xs)
            fs[i]=f-f(*a)
    elif family=='redundant':
        rng=random.Random(104729+seed);base=fs[:n+1]
        fs=base+[sum(P.base_ring()(rng.randrange(prime))*f for f in base) for _ in range(max(0,m-n-1))]
    return fs[:m],points


def timed_method(fs,method,degree,oracle_limit):
    if method=='scalar':return solve_scalar(fs,degree)
    if method=='quotient':return solve_quotient(fs,degree)[0]
    if method.startswith('hybrid'):
        rd=int(method[-1]);return solve_quotient(fs,degree,rd,support_only='support' in method)[0]
    if method=='gb':return solve_gb(fs)
    if method=='core_full':return solve_core_quotient(fs)
    if method=='core_stream':return solve_core_quotient(fs,True)
    if method=='shared':return shared_construction(fs,min(4,degree))
    if method=='oracle':return optimistic_oracle(fs,degree,oracle_limit)
    if method.startswith('linear'):return solve_linearized(fs,int(method[-1]))
    if method=='module':return solve_module(fs,multiplier_degree=2)
    if method=='module_units':return solve_module(fs,units=True,multiplier_degree=2)
    if method=='module_dixon':return solve_module(fs,min(4,degree))
    if method=='module_unit_dixon':return solve_module(fs,min(4,degree),units=True)
    if method=='module_stream':return solve_module(fs,min(4,degree),units=True,stream=True)
    if method=='module_nf':return solve_module_nf(fs,min(4,degree))
    if method=='module_nf5':return solve_module_nf(fs,5)
    if method=='module_nf_full':return solve_module_nf(fs,fs[0].parent().ngens()+1)
    if method=='module_nf2_full':return solve_module_nf(fs,fs[0].parent().ngens()+1,multiplier_degree=2)
    if method=='module_low':return solve_module(fs,units=True,multiplier_degree=1)
    if method=='module3':return solve_module(fs,units=True,multiplier_degree=3)
    raise ValueError(method)


def timeout_handler(sig,frame):
    raise TimeoutError('per-method time limit')


def audit_result(fs,r,reference):
    if '_audit_relations' in r:
        relations=r.pop('_audit_relations');I=fs[0].parent().ideal(fs)
        assert all(I.reduce(g)==0 for g in relations)
        r['nf_relation_membership_checked']=len(relations)
    if r.get('certified'):
        assert r['dimension']==reference['dimension'],(r,reference)
        if 'eliminant' in r:
            assert r['eliminant']==reference['eliminant'],(r,reference)
        if r.get('recovered') is not None:
            assert r['recovered']==reference['recovered']
            assert all(f(*r['recovered'])==0 for f in fs)
        r['independent_gb_equal']=True


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--variables',type=int,default=5);ap.add_argument('--extras',default='0,1,3')
    ap.add_argument('--family',default='dense');ap.add_argument('--prime',type=int,default=257)
    ap.add_argument('--degree',type=int,default=5);ap.add_argument('--seeds',type=int,default=3)
    ap.add_argument('--repeats',type=int,default=3);ap.add_argument('--timeout',type=int,default=30)
    ap.add_argument('--oracle-limit',type=int,default=24)
    ap.add_argument('--methods',default='scalar,quotient,hybrid4,hybrid_support4,module,module_units,module_dixon,module_unit_dixon,module_stream,core_full,core_stream,gb,shared,linear2')
    ap.add_argument('--output',type=Path,required=True);ap.add_argument('--no-warmup',action='store_true')
    args=ap.parse_args();methods=args.methods.split(',');extras=list(map(int,args.extras.split(',')))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    inputdir=args.output.parent/'inputs'/args.output.stem;inputdir.mkdir(parents=True,exist_ok=True)
    signal.signal(signal.SIGALRM,timeout_handler);results=[];disabled=set()
    if not args.no_warmup:
        fs,_=inputs(args.variables,args.variables+max(extras),999,args.family,args.prime)
        for method in methods:
            try:
                signal.alarm(args.timeout);timed_method(fs,method,args.degree,min(4,args.oracle_limit))
            except TimeoutError:pass
            finally:signal.alarm(0)
    for seed in range(args.seeds):
        for extra in extras:
            n=args.variables;m=n+extra;fs,points=inputs(n,m,seed,args.family,args.prime)
            assert all(all(f(*a)==0 for f in fs) for a in points)
            payload='%d %d %d\n'%(args.prime,n,m)+'\n'.join(map(str,fs))+'\n'
            path=inputdir/('seed%d_m%d.txt'%(seed,m));path.write_text(payload)
            reference=solve_gb(fs)
            for repeat in range(args.repeats):
                order=methods[:];random.Random(seed*1000+m*10+repeat).shuffle(order)
                for method in order:
                    if (seed,m,method) in disabled:continue
                    # Oracle is an offline diagnostic, repeated timings are
                    # not meaningful; its output contains the chosen subset.
                    if method=='oracle' and repeat:continue
                    start=time.perf_counter()
                    try:
                        signal.alarm(args.timeout);r=timed_method(fs,method,args.degree,args.oracle_limit)
                        signal.alarm(0);audit_result(fs,r,reference)
                    except TimeoutError:
                        r=dict(seconds=time.perf_counter()-start,certified=False,timeout=True)
                        disabled.add((seed,m,method))
                    finally:signal.alarm(0)
                    if not r.get('certified',False) and method not in ('shared','oracle'):
                        # Exact deterministic failure at this cap does not need
                        # repeated timing samples; retain the failed run.
                        disabled.add((seed,m,method))
                    r.update(n=n,m=m,prime=args.prime,family=args.family,seed=seed,repeat=repeat,
                             method=method,degree=args.degree,reference_dimension=reference['dimension'],
                             input_sha256=hashlib.sha256(payload.encode()).hexdigest(),input_file=str(path))
                    results.append(r);args.output.write_text(json.dumps(results,indent=2)+'\n')
                    print(json.dumps({k:r[k] for k in ('n','m','family','seed','repeat','method','seconds','certified','dimension','timeout') if k in r}),flush=True)
    for extra in extras:
        for method in methods:
            rr=[r for r in results if r['m']==args.variables+extra and r['method']==method]
            if rr:print('SUMMARY',extra,method,statistics.median(r['seconds'] for r in rr),
                        sum(r.get('certified',False) for r in rr),'/',len(rr),flush=True)


if __name__=='__main__':main()

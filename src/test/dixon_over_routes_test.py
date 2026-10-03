"""Independent ideal and native audits, including nonreduced/multiple roots."""
from sage.all import GF, PolynomialRing
from dixon_over_routes import *
from dixon_over_routes_bench import inputs, audit_result
from pathlib import Path
import argparse, json, subprocess, tempfile


def native(fs, binary, method, degree, directory):
    P=fs[0].parent();path=Path(directory)/'case.txt'
    path.write_text('%d %d %d\n'%(P.base_ring().characteristic(),P.ngens(),len(fs))+'\n'.join(map(str,fs))+'\n')
    return json.loads(subprocess.check_output([str(binary),str(path),method,str(degree),'1'],text=True))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--binary',type=Path,default=Path('build/dixon_over_native'))
    ap.add_argument('--output',type=Path);args=ap.parse_args();counts=dict(systems=0,native=0,shared=0,module=0,nf=0)
    with tempfile.TemporaryDirectory(prefix='dixon-over-audit-') as directory:
        for p in (2,3,7,257):
            for n in (3,4):
                for seed in (0,1):
                    for family in ('dense','two'):
                        fs,_=inputs(n,n+2,seed,family,p);I=fs[0].parent().ideal(fs)
                        if I.dimension()>0:continue
                        ref=solve_gb(fs);r,_=solve_quotient(fs,4);audit_result(fs,r,ref)
                        assert r['certified'],(p,n,seed,family,r)
                        for method in ('quotient','adaptive'):
                            nr=native(fs,args.binary,method,4,directory);audit_result(fs,nr,ref)
                            assert nr['certified'],(p,n,seed,family,nr)
                            counts['native']+=1
                        nr=native(fs,args.binary,'shared',4,directory);assert nr['audit'];counts['native']+=1
                        shared_construction(fs,4);counts['shared']+=1
                        for units in (False,True):
                            mr=solve_module(fs,units=units,multiplier_degree=2);audit_result(fs,mr,ref);counts['module']+=1
                        nr=solve_module_nf(fs,4);audit_result(fs,nr,ref);counts['nf']+=1
                        M,cols=module_matrix(fs)
                        for i,f in enumerate(fs):
                            reconstructed=sum(f.parent()({e+(k,):c}) for a,e in zip(M[i],cols) for k,c in enumerate(a.list()))
                            assert reconstructed==f
                        counts['systems']+=1
        for p in (2,7,257):
            P=PolynomialRing(GF(p),names=('x0','x1'));x,y=P.gens()
            fixtures=[[x*x,y*y],[x*x,(y-x)**2],[x*x-x,y*y],[x*x-x,y*y-y,y-x],
                      [x*x,y*y,P.one()],[x*x,y*y,P.zero(),x*x]]
            for fs in fixtures:
                ref=solve_gb(fs);r,_=solve_quotient(fs,4);audit_result(fs,r,ref);assert r['certified']
                for units in (False,True):
                    mr=solve_module(fs,units=units,multiplier_degree=2);audit_result(fs,mr,ref);assert mr['certified'];counts['module']+=1
                nr=solve_module_nf(fs,4);audit_result(fs,nr,ref);counts['nf']+=1
                for method in ('quotient','adaptive'):
                    nr=native(fs,args.binary,method,4,directory);audit_result(fs,nr,ref);assert nr['certified'];counts['native']+=1
                counts['systems']+=1
            fs=[x*x,P.zero()];r,_=solve_quotient(fs,4);assert not r['certified']
            nr=native(fs,args.binary,'adaptive',4,directory);assert not nr['certified'];counts['native']+=1
    print('PASS',json.dumps(counts),flush=True)
    if args.output:args.output.write_text(json.dumps(dict(passed=True,**counts),indent=2)+'\n')


if __name__=='__main__':main()

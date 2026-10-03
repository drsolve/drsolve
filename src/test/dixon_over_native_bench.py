"""Native paired timings, run after the serial Sage sweep has finished.

Fixed-degree baselines get an untimed degree search (a deliberately strong
baseline); adaptive includes every attempted degree. One warmup then five
samples per process, three shuffled rounds. Parsing/process startup are excluded
for every native method. Reference dimensions/eliminants come from Sage GB.
"""
from pathlib import Path
import argparse, hashlib, json, random, statistics, subprocess


def invoke(binary,path,method,degree,repeats=1,audit=False):
    return [json.loads(s) for s in subprocess.check_output(
        [str(binary),str(path),method,str(degree),str(int(audit)),str(repeats)],text=True).splitlines()]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,default=Path('src/test/data/dixon_over_routes'))
    ap.add_argument('--binary',type=Path,default=Path('build/dixon_over_native'))
    ap.add_argument('--rounds',type=int,default=3);ap.add_argument('--samples',type=int,default=5)
    ap.add_argument('--labels',default='dense5,dense6,dense7,two5,partial5,structured6,structured7,redundant5,square5,shared5,shared6')
    ap.add_argument('--output',default='native.json');args=ap.parse_args();inputs={}
    for label in args.labels.split(','):
        for r in json.loads((args.directory/(label+'.json')).read_text()):
            key=r['input_sha256']
            if key not in inputs:inputs[key]=dict(meta=r,reference=None,labels=[])
            if label not in inputs[key]['labels']:inputs[key]['labels'].append(label)
            if r.get('certified') and 'eliminant' in r:inputs[key]['reference']=r
    results=[]
    for key,case in inputs.items():
        meta=case['meta'];reference=case['reference'];path=Path(meta['input_file'])
        assert hashlib.sha256(path.read_bytes()).hexdigest()==key
        maxdegree=meta['degree'];methods=[('shared',min(4,maxdegree))]
        if reference and meta['m']<=meta['n']+3:
            methods += [('adaptive',maxdegree)]
            for method in ('quotient','seed','macaulay'):
                for d in range(2,maxdegree+1):
                    probe=invoke(args.binary,path,method,d)[0]
                    complete=probe.get('certified',False) if method=='quotient' else probe['nullity']==reference['dimension']
                    if complete:
                        methods.append((method,d));break
        for round_number in range(args.rounds):
            order=methods[:];random.Random(int(key[:8],16)+round_number).shuffle(order)
            for method,d in order:
                runs=invoke(args.binary,path,method,d,1+args.samples,audit=True)[1:]
                for r in runs:
                    if method in ('quotient','adaptive'):
                        assert r['certified'] and r['dimension']==reference['dimension'],(meta,r,reference)
                        assert r['eliminant']==reference['eliminant'] and r['recovered']==reference['recovered']
                    elif method in ('seed','macaulay'):
                        assert r['nullity']==reference['dimension']
                        if reference['dimension']==1:assert r['recovered']==reference['recovered']
                r=dict(runs[-1]);r['timing_samples']=[{k:v for k,v in rr.items() if k.endswith('seconds')} for rr in runs]
                for k in r['timing_samples'][0]:r[k]=statistics.median(rr[k] for rr in r['timing_samples'])
                r.update({k:meta[k] for k in ('n','m','prime','family','seed','input_sha256','reference_dimension')})
                r.update(round=round_number,degree_cap=d,degree=r.get('degree',d),labels=case['labels'],fixed_degree_selected_untimed=method in ('quotient','seed','macaulay'))
                results.append(r)
        (args.directory/args.output).write_text(json.dumps(results,indent=2)+'\n')
        print('DONE',meta['family'],meta['n'],meta['m'],meta['seed'],methods,flush=True)


if __name__=='__main__':main()

"""Summarize per-seed medians; do not pool repeated measurements as instances."""
from pathlib import Path
import collections,json,statistics,math

root=Path(__file__).resolve().parents[2]/'src/test/data/dixon_over_routes'
summary={}
for engine in ('sage','native'):
    rows=[]
    if engine=='sage':
        for path in root.glob('*.json'):
            if path.stem in ('audit','native','summary','suite_runs','manifest','production','asan_audit','oracle_native'):continue
            data=json.loads(path.read_text())
            if isinstance(data,list):
                rows += [dict(r,dataset=path.stem) for r in data if isinstance(r,dict) and 'method' in r]
    elif (root/'native.json').exists():rows=json.loads((root/'native.json').read_text())
    groups=collections.defaultdict(list)
    for r in rows:
        groups[(r.get('dataset','native'),r['family'],r['n'],r['m'],r['method'])].append(r)
    table=[]
    for key,rr in sorted(groups.items()):
        byseed=collections.defaultdict(list)
        for r in rr:byseed[r['seed']].append(r)
        times=[statistics.median(r['seconds'] for r in s) for s in byseed.values()]
        item=dict(dataset=key[0],family=key[1],n=key[2],m=key[3],method=key[4],seeds=len(byseed),
                  runs=len(rr),seconds=statistics.median(times),min_seed_seconds=min(times),max_seed_seconds=max(times),
                  certified=sum(r.get('certified',r.get('recovered') is not None) for r in rr),
                  timeouts=sum(r.get('timeout',False) for r in rr))
        for field in ('cpu_seconds','dixon_seconds','dixon_coefficients','core_setup_seconds','relations',
                      'unit_pivots','shared_term_products','repeated_term_products','repeated_seconds',
                      'free_generation_work_ratio'):
            vals=[statistics.median(r[field] for r in s if field in r) for s in byseed.values() if any(field in r for r in s)]
            if vals:item[field]=statistics.median(vals)
        if key[4]=='shared':
            item['shared_over_repeated_ratio']=statistics.median(
                statistics.median(r['seconds']/r['repeated_seconds'] for r in s) for s in byseed.values())
        table.append(item)
    summary[engine]=table
if (root/'native.json').exists():
    native=json.loads((root/'native.json').read_text());paired=[]
    for n in (5,6,7):
        times={};degrees={};points={}
        for m in (n+1,n+3):
            for seed in range(3):
                rr=[r for r in native if r['family']=='dense' and r['n']==n and r['m']==m and r['seed']==seed and r['method']=='adaptive']
                assert rr and all(r['certified'] and r['dimension']==1 for r in rr)
                times[m,seed]=statistics.median(r['seconds'] for r in rr)
                degrees[m,seed]=rr[0]['degree'];points[m,seed]=rr[0]['recovered']
                assert all(r['recovered']==points[m,seed] and r['degree']==degrees[m,seed] for r in rr)
        assert all(points[n+1,s]==points[n+3,s] for s in range(3))
        paired.append(dict(n=n,equations=[n+1,n+3],
                           seconds=[statistics.median(times[m,s] for s in range(3)) for m in (n+1,n+3)],
                           paired_speedup=statistics.median(times[n+1,s]/times[n+3,s] for s in range(3)),
                           degree=[degrees[n+1,0],degrees[n+3,0]],
                           monomial_columns=[math.comb(n+degrees[m,0],degrees[m,0]) for m in (n+1,n+3)],
                           same_simple_point_ideal=True))
    summary['native_equation_gain']=paired
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
for engine,table in summary.items():
    if engine not in ('sage','native'):continue
    print(engine.upper())
    for r in table:
        print(r['dataset'],r['family'],r['n'],r['m'],r['method'],round(1000*r['seconds'],3),
              '%d/%d'%(r['certified'],r['runs']),
              ('shared/repeated=%.3f'%r['shared_over_repeated_ratio']) if 'shared_over_repeated_ratio' in r else '')

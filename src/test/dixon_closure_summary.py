#!/usr/bin/env python3
"""Summarize paired seeds, retaining per-seed medians before aggregation."""
import argparse
import json
from pathlib import Path
import statistics


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('directory',nargs='?',type=Path,default=Path('src/test/data/dixon_closure'))
    args=ap.parse_args(); summaries=[]
    native_path=args.directory/'native_paired.json'
    native=json.loads(native_path.read_text()) if native_path.exists() else []
    for path in sorted(args.directory.glob('n[0-9]*.json')):
        records=json.loads(path.read_text()); methods=sorted({r['method'] for r in records})
        seeds=sorted({r['seed'] for r in records}); per_seed={}
        report={'dataset':path.stem,'n':records[0]['n'],'m':records[0]['m'],
                'degree':records[0]['degree'],'prime':records[0]['prime'],
                'family':records[0]['family'],'planted_points':records[0]['points'],
                'seeds':seeds,'paired_runs':len(records),'methods':{},'comparisons':{}}
        for method in methods:
            rr=[r for r in records if r['method']==method]
            ss={seed:[r for r in rr if r['seed']==seed] for seed in seeds}
            nr=[r for r in native if r['dataset']==path.stem and r['method']==method]
            native_groups={seed:[r for r in nr if r['seed']==seed] for seed in seeds}
            if not nr: native_groups={seed:[r['c'] for r in group] for seed,group in ss.items()}
            per_seed[method]={seed:statistics.median(r['seconds'] for r in group) for seed,group in native_groups.items()}
            stats={
                'c_seconds':statistics.median(per_seed[method].values()),
                'sage_seconds':statistics.median(statistics.median(r['seconds'] for r in group) for group in ss.values()),
                'native_timing_pass':'independent paired C' if nr else 'Sage/C crosscheck',
                'c_dixon_seconds':statistics.median(statistics.median(r['dixon_seconds'] for r in group) for group in native_groups.values()),
                'c_closure_seconds':statistics.median(statistics.median(r['closure_seconds'] for r in group) for group in native_groups.values()),
                'nullities':sorted({r['nullity'] for r in rr}),
                'unique_point_runs':sum(r['recovered'] is not None for r in rr),
                'per_seed_c_seconds':per_seed[method],
                'example_work':{k:rr[0][k] for k in ('rank_history','initial_rows','columns','submitted_rows','rref_cells','dense_matmul_multiply_adds')},
                'example_dixon_work':rr[0]['dixon'],
            }
            if nr:
                stats['c_cpu_seconds']=statistics.median(statistics.median(r['cpu_seconds'] for r in group) for group in native_groups.values())
            report['methods'][method]=stats
        for baseline,variant in [('seed','hybrid3'),('seed','hybrid4'),('macaulay','augment3'),('macaulay','augment4'),('macaulay','hybrid3'),('macaulay','hybrid4')]:
            if variant not in methods: continue
            ratios=[per_seed[baseline][s]/per_seed[variant][s] for s in seeds]
            report['comparisons'][baseline+'/'+variant]={
                'paired_speedup_median':statistics.median(ratios),'minimum':min(ratios),'maximum':max(ratios),
                'seeds_faster':sum(r>1 for r in ratios),'seeds_faster_by_5_percent':sum(r>1.05 for r in ratios)}
        summaries.append(report)
        print(path.stem,'seeds',len(seeds),'nullity',report['methods']['seed']['nullities'])
        print('  C milliseconds:',', '.join('%s %.3f'%(method,1000*report['methods'][method]['c_seconds']) for method in methods))
    (args.directory/'summary.json').write_text(json.dumps(summaries,indent=2)+'\n')


if __name__=='__main__': main()

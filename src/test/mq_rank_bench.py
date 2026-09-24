#!/usr/bin/env python3
"""Interleaved direct-index/hash timings and exact sparse-DP audits."""
import argparse
import json
from pathlib import Path
import statistics
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('phase', choices=['audit', 'timing'])
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    configs = []
    if args.phase == 'audit':
        for n in (7, 8):
            for threads in (1, 4):
                for seed in (132, 12345):
                    configs.append((n, threads, 13, 1, 0, 0, 1, seed, 65537))
            configs.append((n, 4, 13, 0, 0, 0, 1, 132, 65537))
        for prime in (2, 7, 18446744073709551557):
            for seed in (132, 12345):
                configs.append((5, 3, 13, 1, 0, 0, 1, seed, prime))
        configs.append((8, 16, 13, 1, 0, 0, 1, 132, 65537))
    else:
        for n in (7, 8):
            for threads in (1, 4):
                for repeat in range(args.repeats):
                    for mode in ((14, 13) if repeat % 2 == 0 else (13, 14)):
                        configs.append((n, threads, mode, 1, 0, 0, 0, 132, 65537))
    output = dict(phase=args.phase, runs=[])
    fingerprints = {}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for i, config in enumerate(configs):
        p = subprocess.run([str(ROOT/'build/mq_layout_bench'), *map(str,config)],
                           cwd=ROOT, text=True, capture_output=True, timeout=300)
        if p.returncode:
            raise RuntimeError(f'{config}: {p.stdout}\n{p.stderr}')
        records = [json.loads(line) for line in p.stdout.splitlines() if line.startswith('{')]
        r = next(r for r in records if r['kind']=='run')
        key = (r['n'], r['q'], r['seed'], r['projected'])
        fingerprint = (r['fingerprint'], r['terms'], r['verified'])
        assert fingerprints.setdefault(key, fingerprint) == fingerprint
        output['runs'].append(dict(config=config, records=records))
        args.output.write_text(json.dumps(output,indent=2)+'\n')
        print(f'{i+1}/{len(configs)} n={config[0]} threads={config[1]} mode={config[2]} {r["seconds"]:.5f}s',flush=True)
    if args.phase == 'timing':
        groups = {}
        for item in output['runs']:
            r = next(r for r in item['records'] if r['kind']=='run')
            groups.setdefault((r['n'], r['threads'], r['shared']), []).append(r)
        output['medians'] = [dict(n=k[0], threads=k[1], mode=k[2],
            **{field:statistics.median(r[field] for r in group)
               for field in ('seconds','plan_seconds','arithmetic_seconds','peak_rss_kib_before_audit')})
            for k,group in sorted(groups.items())]
        args.output.write_text(json.dumps(output,indent=2)+'\n')
        print(json.dumps(output['medians'],indent=2))


if __name__ == '__main__':
    main()

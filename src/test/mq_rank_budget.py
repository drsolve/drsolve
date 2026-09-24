#!/usr/bin/env python3
"""Count a rankable axis/degree envelope and report the map-time ceiling.

This is a feasibility analysis, not a replacement determinant backend.
"""
import argparse
from bisect import bisect_left
from collections import Counter
import itertools
import json
from pathlib import Path
import statistics


def closure(targets):
    return {e for t in targets for e in itertools.product(*(range(v+1) for v in t))}


def envelope(model, layers):
    n = model['n']
    m = n-1
    axes = [closure(model[name]) for name in ('rows', 'cols')]
    active = set()
    result = []
    old_slots = new_slots = 0
    for k, shifts in enumerate(model['shifts'], 1):
        active.update(i for e in shifts for i, v in enumerate(e) if v)
        degree = k + (k == n)
        allowed = [sorted(e for e in axis
                          if all(not v or i+offset in active for i, v in enumerate(e)))
                   for axis, offset in zip(axes, (0, m))]
        hist = [Counter(map(sum, axis)) for axis in allowed]
        count = sum(a*b*max(0, degree-i-j+1)
                    for i, a in hist[0].items() for j, b in hist[1].items())
        actual = layers[k]['support']
        assert count >= actual
        minors = layers[k]['weight']//k if k != n else n
        old_slots += minors*actual
        new_slots += minors*count
        result.append(dict(k=k, support=actual, envelope=count, ratio=count/actual))
        # Independent enumeration audits the count and a hash-free rank formula
        # on the small fixture. Use lex axes and per-budget y prefix sums.
        if n == 4:
            x, y = allowed
            prefix_y = {}
            for budget in range(degree+1):
                prefix = [0]
                for e in y:
                    prefix.append(prefix[-1] + max(0, budget-sum(e)+1))
                prefix_y[budget] = prefix
            offsets = [0]
            for e in x:
                budget = degree-sum(e)
                offsets.append(offsets[-1] + (prefix_y[budget][-1] if budget >= 0 else 0))
            assert offsets[-1] == count
            ranks = set()
            for ex, ey in itertools.product(x, y):
                for t in range(max(0, degree-sum(ex)-sum(ey)+1)):
                    rank = offsets[bisect_left(x, ex)] + prefix_y[degree-sum(ex)][bisect_left(y, ey)] + t
                    ranks.add(rank)
            assert ranks == set(range(count))
    return dict(n=n, layers=result, coefficient_slot_ratio=new_slots/old_slots)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--supports', type=Path, required=True)
    parser.add_argument('--timings', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    records = [json.loads(line) for line in args.supports.read_text().splitlines()]
    supports = []
    for model in (r for r in records if r['kind'] == 'model'):
        layers = {r['k']: r for r in records if r['kind'] == 'layer' and r['n'] == model['n']}
        supports.append(envelope(model, layers))
    timings = json.loads(args.timings.read_text())
    groups = {}
    for item in timings['runs']:
        run = next(r for r in item['records'] if r['kind'] == 'run')
        if run['shared'] == 12:
            groups.setdefault((run['n'], run['threads']), []).append(run)
    budget = []
    for (n, threads), runs in sorted(groups.items()):
        fractions = [r['plan_seconds']/r['seconds'] for r in runs]
        budget.append(dict(n=n, threads=threads, repeats=len(runs),
                           total_median=statistics.median(r['seconds'] for r in runs),
                           map_median=statistics.median(r['plan_seconds'] for r in runs),
                           map_fraction_median=statistics.median(fractions),
                           zero_map_ceiling_median=statistics.median(1/(1-f) for f in fractions)))
    output = dict(envelopes=supports, map_budget=budget,
                  audit='n=4 envelope enumeration and contiguous rank bijection passed')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()

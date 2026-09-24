#!/usr/bin/env python3
"""Support-only MQ analysis and independent small tuple-set oracle."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def run(n):
    text = subprocess.check_output([str(ROOT / 'build/mq_support_count'), str(n)],
                                   cwd=ROOT, text=True)
    return [json.loads(line) for line in text.splitlines()]


def oracle(records):
    model = next(r for r in records if r['kind'] == 'model')
    n = model['n']
    m = n - 1
    targets = [set(map(tuple, model[axis])) for axis in ('rows', 'cols')]
    closures = [{e for t in axis for e in itertools.product(*(range(v+1) for v in t))}
                for axis in targets]
    forward = [{(0,) * (2*n-1)}]
    shifts = [set(map(tuple, row)) for row in model['shifts']]
    for row in shifts:
        forward.append({z for a in forward[-1] for b in row
                        if (z := tuple(x+y for x, y in zip(a, b)))[:m] in closures[0]
                        and z[m:2*m] in closures[1]})
    live = {a for a in forward[-1] if a[:m] in targets[0] and a[m:2*m] in targets[1]}
    layers = {r['k']: r for r in records if r['kind'] == 'layer'}
    for k in range(n, 0, -1):
        assert len(forward[k]) == layers[k]['support']
        assert len(live) == layers[k]['live']
        # Reverse subtraction rather than the C forward-edge lookup.
        previous = set()
        kept_edges = 0
        for dest in live:
            for shift in shifts[k-1]:
                source = tuple(a-b for a, b in zip(dest, shift))
                if source in forward[k-1]:
                    previous.add(source)
                    kept_edges += 1
        assert kept_edges == layers[k]['live_edges']
        edges = sum(tuple(a+b for a, b in zip(src, shift)) in forward[k]
                    for src in forward[k-1] for shift in shifts[k-1])
        assert edges == layers[k]['edges']
        live = previous
    assert live == forward[0]


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    small = run(4)
    oracle(small)
    records = small + run(7) + run(8)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(''.join(json.dumps(r) + '\n' for r in records))
    print('n=4 independent tuple oracle: passed')
    for r in records:
        if r['kind'] == 'summary':
            print(f"n={r['n']}: weighted work reduction {(1-r['work_ratio'])*100:.6f}%")


if __name__ == '__main__':
    main()

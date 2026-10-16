#!/usr/bin/env python3
"""Compare the research backend with drsolve on exactly the same input.

Run from the repository root after make build/mq_quotient_bezout_test.
Example: python3 src/test/mq_quotient_bezout_bench.py --n 8 --threads 16
Validation samples are reported separately from experimental computation.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess


def convert(source, n, destination):
    text = Path(source).read_text()
    prime = int(re.search(r"Field: F_(\d+)", text)[1])
    polys = next(s[13:] for s in text.splitlines() if s.startswith("Polynomials: ")).split(",")
    assert len(polys) == n
    lines = [f"{n} {prime}"]
    for poly in polys:
        coefficients = {}
        for term in poly.strip().split(" + "):
            factors = term.split("*")
            c = int(factors.pop(0)) if factors[0].isdigit() else 1
            exps = [0] * n
            for factor in factors:
                m = re.fullmatch(r"x(\d+)(?:\^(\d+))?", factor)
                assert m, factor
                exps[int(m[1])] += int(m[2] or 1)
            assert sum(exps) <= 2
            coefficients[tuple(exps)] = c
        values = []
        for i in range(n + 1):
            for j in range(i, n + 1):
                exps = [0] * n
                if i < n:
                    exps[i] += 1
                if j < n:
                    exps[j] += 1
                values.append(coefficients.get(tuple(exps), 0))
        lines.append(" ".join(map(str, values)))
    result = text.split("Resultant:\n")[1].split("\n")[0]
    coefficients = {}
    for term in result.split(" + "):
        m = re.fullmatch(rf"(?:(\d+)\*)?x{n-1}(?:\^(\d+))?", term)
        if m:
            coefficients[int(m[2] or 1)] = int(m[1] or 1)
        else:
            coefficients[0] = int(term)
    degree = max(coefficients)
    lines += [str(degree), " ".join(str(coefficients.get(i, 0)) for i in range(degree + 1))]
    Path(destination).write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=8)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--prime", type=int, default=65537)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--mode", type=int, choices=(1, 2, 3), default=3)
    parser.add_argument("--kernel", type=int, choices=(-1, 0, 1, 2, 3, 4, 5, 6), default=-1)
    args = parser.parse_args()
    out = Path("out/jouanolou")
    out.mkdir(exist_ok=True, parents=True)
    env = dict(os.environ, LD_LIBRARY_PATH=".:/usr/local/lib")
    records = []
    for trial in range(args.repeat):
        stem = f"n{args.n}_p{args.prime}_t{args.threads}_s{args.seed}_m{args.mode}_k{args.kernel}_{trial}"
        source, numeric = out / (stem + ".dr"), out / (stem + ".txt")
        p = subprocess.run(["./drsolve", "-r", f"[2]*{args.n}", str(args.prime),
                            "--seed", str(args.seed), "--threads", str(args.threads),
                            "--time", "-o", str(source)], env=env, capture_output=True, text=True, check=True)
        (out / (stem + "_baseline.log")).write_text(p.stdout + p.stderr)
        total = next(line for line in p.stdout.splitlines() if line.startswith("Total -"))
        baseline = float(re.search(r"Wall time: ([0-9.]+)", total)[1])
        convert(source, args.n, numeric)
        p = subprocess.run(["./build/mq_quotient_bezout_test", str(numeric),
                            str(args.threads), str(args.mode), str(args.kernel)], env=env, capture_output=True, text=True, check=True)
        (out / (stem + "_experiment.log")).write_text(p.stdout + p.stderr)
        record = json.loads(p.stdout)
        assert record["reference_equal"] == 1
        record.update(seed=args.seed, trial=trial, baseline_seconds=baseline)
        records.append(record)
        (out / f"results_n{args.n}_p{args.prime}_t{args.threads}_s{args.seed}_m{args.mode}_k{args.kernel}.json").write_text(json.dumps(records, indent=2))
        print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact parallel quotient comparisons and optional reproducible timings."""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--seeds", nargs="+", type=int, default=[1, 12345])
    parser.add_argument("--threads", nargs="+", type=int, default=[1, 2, 4, 8])
    parser.add_argument("--cases", nargs="+", default=["[3]*5", "[8]*3", "[2,3,4]"])
    args = parser.parse_args()
    assert args.repeat > 0 and min(args.threads) > 0
    out = Path("out/degree_quotient_parallel")
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, LD_LIBRARY_PATH=".:/usr/local/lib", OMP_DYNAMIC="FALSE")
    records = []
    for case in args.cases:
        for seed in args.seeds:
            reference = None
            for repeat in range(args.repeat):
                # Both algorithms receive the same thread count. No concurrent
                # benchmark processes; compare reported compute wall times.
                for threads in args.threads:
                    for quotient in (False, True):
                        stem = f"{case}_{seed}_{repeat}_{threads}_{quotient}"
                        path = out / (stem + ".dr")
                        p = subprocess.run(
                            ["./drsolve", "-r", case, "65537", "--seed", str(seed),
                             "--quotient" if quotient else "--no-quotient",
                             "--threads", str(threads), "--time", "-v", "1", "-o", str(path)],
                            env=env, capture_output=True, text=True, timeout=120, check=True)
                        log = p.stdout + p.stderr
                        (out / (stem + ".log")).write_text(log)
                        result = path.read_text().split("Resultant:\n", 1)[1].splitlines()[0]
                        if reference is None:
                            reference = result
                        assert result == reference, stem
                        row = dict(case=case, seed=seed, repeat=repeat, threads=threads,
                                   quotient=quotient,
                                   wall=float(re.search(r"Total -.*Wall time: ([0-9.]+)", log)[1]))
                        if quotient:
                            assert "fallback:" not in log, log
                            row["workers"] = int(re.search(r"sample workers: ([0-9]+)", log)[1])
                            assert 1 <= row["workers"] <= threads
                            row["sampling"] = float(re.search(r"Step 2 time: ([0-9.]+)", log)[1])
                        records.append(row)
                        (out / "results.json").write_text(json.dumps(records, indent=2) + "\n")
            print(f"PASS {case} seed={seed}", flush=True)
    summary = []
    for case in args.cases:
        for threads in args.threads:
            row = dict(case=case, threads=threads)
            for quotient, name in [(False, "dixon_ms"), (True, "quotient_ms")]:
                chosen = [r for r in records if r["case"] == case and
                          r["threads"] == threads and r["quotient"] == quotient]
                row[name] = 1000 * statistics.median(r["wall"] for r in chosen)
                if quotient:
                    row["workers"] = sorted(set(r["workers"] for r in chosen))
            summary.append(row)
            print(json.dumps(row), flush=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()

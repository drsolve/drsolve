#!/usr/bin/env python3
"""Compare two --mq-quotient executables, including evaluation time and RSS.

Requires a baseline executable saved before rebuilding, e.g.
python3 src/test/benchmarks/mq_quotient_eval_bench.py --baseline /tmp/drsolve-before
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--n", type=int, nargs="+", default=[9, 10])
    parser.add_argument("--repeat", type=int, default=2)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--out", default="out/jouanolou/packed_comparison")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, LD_LIBRARY_PATH=".:/usr/local/lib")
    records = []
    for n in args.n:
        reference = None
        for trial in range(args.repeat):
            for label, executable in (("before", args.baseline), ("after", "./drsolve")):
                stem = f"n{n}_t{args.threads}_{trial}_{label}"
                output, rss = out / (stem + ".dr"), out / (stem + ".rss")
                cmd = ["/usr/bin/time", "-f", "%M", "-o", str(rss), "timeout", "180",
                       executable, "-r", f"[2]*{n}", "65537", "--seed", "12345",
                       "--threads", str(args.threads), "--time", "--mq-quotient", "-o", str(output)]
                p = subprocess.run(cmd, env=env, capture_output=True, text=True, check=True, timeout=190)
                log = p.stdout + p.stderr
                (out / (stem + ".log")).write_text(log)
                assert "  Core:" in log and "MQ quotient fallback:" not in log
                result = output.read_text().split("Resultant:\n", 1)[1].splitlines()[0]
                if reference is None:
                    reference = result
                assert result == reference
                total = next(line for line in log.splitlines() if line.startswith("Total -"))
                record = dict(n=n, threads=args.threads, trial=trial, version=label,
                              leading=float(re.search(r"Shared decomposition ([0-9.]+)s", log)[1]),
                              evaluation=float(re.search(r"evaluations ([0-9.]+)s", log)[1]),
                              total=float(re.search(r"Wall time: ([0-9.]+)", total)[1]),
                              peak_rss_kib=int(rss.read_text().strip()))
                records.append(record)
                (out / "results.json").write_text(json.dumps(records, indent=2))
                print(json.dumps(record), flush=True)
    print("All resultants match exactly.")


if __name__ == "__main__":
    main()

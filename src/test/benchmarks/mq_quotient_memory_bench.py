#!/usr/bin/env python3
"""Compare complete CLI peak RSS (GNU time), runtime and resultants.

Example: python3 src/test/benchmarks/mq_quotient_memory_bench.py --n 8 9 --threads 1 16
Each measurement runs in a fresh process. RSS includes the runtime and libraries.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", nargs="+", type=int, default=[8, 9])
    parser.add_argument("--threads", nargs="+", type=int, default=[1, 16])
    parser.add_argument("--timeout", type=int, default=180)
    args = parser.parse_args()
    out = Path("out/jouanolou/memory_dynamic")
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, LD_LIBRARY_PATH=".:/usr/local/lib")
    records = []
    for n in args.n:
        for threads in args.threads:
            results = []
            for enabled in (True, False):
                stem = f"n{n}_t{threads}_q{int(enabled)}"
                source, rss = out / (stem + ".dr"), out / (stem + ".rss")
                cmd = ["/usr/bin/time", "-f", "%M", "-o", str(rss),
                       "timeout", str(args.timeout), "./drsolve", "-r", f"[2]*{n}",
                       "65537", "--seed", "12345", "--threads", str(threads),
                       "--time", "-o", str(source)]
                if enabled:
                    cmd.append("--mq-quotient")
                p = subprocess.run(cmd, env=env, capture_output=True, text=True,
                                   check=True, timeout=args.timeout + 10)
                log = p.stdout + p.stderr
                (out / (stem + ".log")).write_text(log)
                if enabled:
                    assert "  Core:" in log and "MQ quotient fallback:" not in log, log
                results.append(source.read_text().split("Resultant:\n", 1)[1].splitlines()[0])
                total = next(line for line in log.splitlines() if line.startswith("Total -"))
                record = dict(n=n, threads=threads, quotient=enabled,
                              seconds=float(re.search(r"Wall time: ([0-9.]+)", total)[1]),
                              peak_rss_kib=int(rss.read_text().strip()))
                records.append(record)
                path = out / ("results_n" + "_".join(map(str, args.n)) +
                              "_t" + "_".join(map(str, args.threads)) + ".json")
                path.write_text(json.dumps(records, indent=2))
                print(json.dumps(record), flush=True)
            assert results[0] == results[1], (n, threads)
    print("All resultants match exactly.")


if __name__ == "__main__":
    main()

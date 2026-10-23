#!/usr/bin/env python3
"""Exact CLI comparisons and fallback checks for --mq-quotient.

Small deterministic systems exercise construction and fallback paths.
"""
import argparse
import json
import os
from pathlib import Path
import re
import statistics
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--binary", default="./drsolve")
    args = parser.parse_args()
    out = Path("out/jouanolou/cli")
    out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, LD_LIBRARY_PATH=".:/usr/local/lib")
    records = []

    def check(name, arguments, accepted, repeat=1):
        for trial in range(repeat):
            results, times = [], []
            for enabled in (False, True):
                stem = f"{name}_{trial}_{int(enabled)}"
                path = out / (stem + ".dr")
                cmd = [args.binary, *arguments, "--time", "--threads", "4", "-o", str(path)]
                if enabled:
                    cmd.append("--mq-quotient")
                else:
                    cmd.append("--no-quotient")
                p = subprocess.run(cmd, env=env, capture_output=True, text=True, check=True, timeout=120)
                log = p.stdout + p.stderr
                (out / (stem + ".log")).write_text(log)
                if enabled:
                    assert ("; normalize to monic." in log) == accepted, log
                    assert ("quotient fallback:" in log) != accepted, log
                else:
                    assert "Experimental Dixon resultant via quotient-space Bezoutian" not in log
                results.append(path.read_text().split("Resultant:\n", 1)[1].splitlines()[0])
                total = next(line for line in log.splitlines() if line.startswith("Total -"))
                times.append(float(re.search(r"Wall time: ([0-9.]+)", total)[1]))
            assert results[0] == results[1], name
            record = dict(case=name, trial=trial, accepted=accepted,
                          dixon_seconds=times[0], quotient_seconds=times[1])
            records.append(record)
            (out / "results.json").write_text(json.dumps(records, indent=2))
            print(json.dumps(record), flush=True)

    for n, prime in ((3, 65537), (4, 536870923), (4, 37)):
        check(f"n{n}_p{prime}", ["-r", f"[2]*{n}", str(prime), "--seed", "12345"],
              True, args.repeat)
    check("seed54321", ["-r", "[2]*4", "65537", "--seed", "54321"], True)
    check("recursive", ["-r", "[2]*4", "65537", "--seed", "12345", "--method", "5"], True)
    check("field_equation_final", ["-r", "[2]*3", "37", "--seed", "12345", "--field-equation-s"], False)
    check("cubic", ["-r", "[3]*3", "65537", "--seed", "12345"], True)
    check("small_field", ["-r", "[2]*3", "3", "--seed", "12345"], False)
    check("extension", ["-r", "[2]*3", "2^8", "--seed", "12345"], True)
    check("singular_leading", ["x^2+y^2+z,x^2+y^2+x+1,x^2+y^2+y+2", "x,y", "65537"], False)
    check("degree_drop", ["x^2+y,x^2+y+1", "x", "65537"], False)
    for name in ("n3_p65537", "n4_p37"):
        subset = [r for r in records if r["case"] == name]
        print(name, {key: statistics.median(r[key] for r in subset)
                     for key in ("dixon_seconds", "quotient_seconds")})


if __name__ == "__main__":
    main()

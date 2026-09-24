#!/usr/bin/env python3
"""Sequential MQ layout measurements; never time competing runs concurrently."""
import argparse
import json
import pathlib
import platform
import statistics
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("profile", "timing", "row-order", "audit", "production"))
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    configs = []
    if args.phase == "production":
        for n in (7, 8):
            for threads in (1, 4):
                for projected in (0, 1):
                    for seed in (132, 12345):
                        configs.append((n, threads, 4, projected, 0, 0, 1, seed, 65537))
        for prime in (2, 7, 101, 18446744073709551557):
            for projected in (0, 1):
                configs.append((5, 1, 4, projected, 0, 0, 1, 12345, prime))
    elif args.phase == "profile":
        for n in (7, 8):
            for rotate in (0, 1):
                for mode in (0, 3):
                    configs.append((n, 1, mode, 1, rotate, 1, 1, 132, 65537))
    elif args.phase == "timing":
        for n in (7, 8):
            for threads in (1, 4):
                for repeat in range(args.repeats):
                    # Alternate ordering to reduce consistent first-run bias.
                    for mode in ((0, 1, 2, 3) if repeat % 2 == 0 else (3, 2, 1, 0)):
                        configs.append((n, threads, mode, 1, 0, 0, 0, 132, 65537))
            for repeat in range(args.repeats):
                for mode in ((0, 3) if repeat % 2 == 0 else (3, 0)):
                    configs.append((n, 1, mode, 1, 1, 0, 0, 132, 65537))
    elif args.phase == "row-order":
        for n in (7, 8):
            for repeat in range(args.repeats):
                pairs = [(0, 0), (0, 1), (3, 0), (3, 1)]
                if repeat % 2: pairs.reverse()
                for mode, rotate in pairs:
                    configs.append((n, 1, mode, 1, rotate, 0, 0, 132, 65537))
    else:
        for n in (7, 8):
            for mode in (1, 2):
                for threads in (1, 4):
                    configs.append((n, threads, mode, 1, 0, 0, 1, 132, 65537))
            for seed in (132, 12345):
                for projected in (0, 1):
                    for threads in (1, 4):
                        configs.append((n, threads, 3, projected, 0, 0, 1, seed, 65537))
            configs.append((n, 1, 3, 0, 1, 0, 1, 132, 65537))
        for prime in (2, 7, 101, 18446744073709551557):
            for projected in (0, 1):
                for rotate in (0, 1):
                    configs.append((5, 1, 3, projected, rotate, 0, 1, 12345, prime))
    output = {"phase": args.phase, "platform": platform.platform(), "runs": []}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fingerprints = {}
    for i, config in enumerate(configs):
        command = [str(ROOT / "build/mq_layout_bench"), *map(str, config)]
        proc = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, timeout=300)
        if proc.returncode:
            raise RuntimeError(f"{command}\n{proc.stdout}\n{proc.stderr}")
        records = [json.loads(line) for line in proc.stdout.splitlines() if line.startswith("{")]
        run = next(record for record in records if record["kind"] == "run")
        key = (run["n"], run["q"], run["seed"], run["projected"])
        fingerprint = (run["fingerprint"], run["terms"], run["verified"])
        if key in fingerprints:
            assert fingerprints[key] == fingerprint, (key, run)
        fingerprints[key] = fingerprint
        output["runs"].append({"command": command[1:], "records": records})
        args.output.write_text(json.dumps(output, indent=2) + "\n")
        print(f"{i+1}/{len(configs)} n={config[0]} threads={config[1]} mode={config[2]} "
              f"rotate={config[4]} {run['seconds']:.4f}s", flush=True)
    if args.phase in ("timing", "row-order"):
        groups = {}
        for item in output["runs"]:
            run = item["records"][-1]
            key = (run["n"], run["threads"], run["shared"], run["quadratic_first"])
            groups.setdefault(key, []).append(run)
        output["medians"] = [dict(n=key[0], threads=key[1], shared=key[2], quadratic_first=key[3],
                                  seconds=statistics.median(r["seconds"] for r in runs),
                                  peak_rss_mib=statistics.median(r["peak_rss_kib_before_audit"] for r in runs)/1024)
                             for key, runs in sorted(groups.items())]
        args.output.write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()

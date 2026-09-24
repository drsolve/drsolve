#!/usr/bin/env python3
"""Production shared-index default, opt-out, precedence and exact CLI outputs."""
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
CASES = [(7, 65537, 12345, []), (8, 65537, 12345, []),
         (5, 7, 7, []), (4, 2, 12345, []),
         (4, 65537, 12345, ["--no-mq-step1-filter"]),
         (4, 65537, 12345, ["--step1", "1"]),
         (4, 65537, 12345, ["--mq-step1-pencil"]),
         (4, 65537, 12345, ["--mq-step1-simplex"])]

with tempfile.TemporaryDirectory(prefix="mq-shared-cli-") as directory:
    for case, (n, q, seed, extra) in enumerate(CASES):
        outputs = []
        for off, threads in [(True, 1), (False, 1), (False, 4)]:
            path = Path(directory) / f"{case}-{off}-{threads}.dr"
            command = [str(ROOT / "drsolve"), "--threads", str(threads), "-v", "2",
                       "-r", "--seed", str(seed), "-n", str(n), "-o", str(path),
                       *extra, f"[2]*{n}", str(q)]
            if off:
                command.insert(1, "--no-mq-step1-shared")
            result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True,
                                    check=True, timeout=180)
            used = "MQ shared-index DP:" in result.stdout
            explicit = any(flag in extra for flag in ("--step1", "--mq-step1-pencil", "--mq-step1-simplex"))
            if off or explicit:
                assert not used, (command, result.stdout)
            elif q == 65537:
                assert used, (command, result.stdout)
            outputs.append("\n".join(line for line in path.read_text().splitlines()
                                     if not line.startswith("Time: ")))
        assert all(value == outputs[0] for value in outputs), (n, q, extra)

for flags, expected in [(["--no-mq-step1-shared", "--mq-step1-shared"], True),
                         (["--mq-step1-shared", "--no-mq-step1-shared"], False)]:
    result = subprocess.run([str(ROOT / "drsolve"), *flags, "--threads", "1", "-v", "2",
                             "-r", "--seed", "12345", "-n", "4", "[2]*4", "65537"],
                            cwd=ROOT, text=True, capture_output=True, check=True, timeout=60)
    assert ("MQ shared-index DP:" in result.stdout) == expected
print("MQ shared CLI: default, opt-out, precedence, n=7/8 and exact full outputs PASS")

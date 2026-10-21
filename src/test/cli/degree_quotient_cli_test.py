#!/usr/bin/env python3
"""Compare explicitly enabled grouped quotient resultants with Dixon."""
import os
from pathlib import Path
import subprocess

output = Path("out/degree_quotient_cli")
output.mkdir(parents=True, exist_ok=True)
env = dict(os.environ, LD_LIBRARY_PATH=".:/usr/local/lib", OMP_NUM_THREADS="1")


def check(name, arguments, grouped, fallback=False):
    results = []
    for enabled in (False, True):
        path = output / f"{name}_{enabled}.dr"
        p = subprocess.run(
            ["./drsolve", *arguments, "--quotient" if enabled else "--no-quotient",
             "--threads", "1", "-v", "1", "-o", str(path)],
            capture_output=True, text=True, env=env, timeout=120, check=True)
        log = p.stdout + p.stderr
        (output / f"{name}_{enabled}.log").write_text(log)
        if enabled:
            assert ("grouped quotient recurrence" in log) == grouped, log
            if fallback:
                assert "fallback:" in log, log
            elif grouped:
                assert "fallback:" not in log, log
                assert "; normalize to monic." in log, log
        results.append(path.read_text().split("Resultant:\n", 1)[1].splitlines()[0])
    assert results[0] == results[1], name
    print("PASS", name, flush=True)


for degrees in ("[3]*3", "[3]*5", "[8]*3", "[2,3,4]", "[1,3,5]", "[20]*2"):
    for seed in (1, 12345):
        check(f"{degrees}_{seed}", ["-r", degrees, "65537", "--seed", str(seed)], True)
check("mq", ["-r", "[2]*4", "65537", "--seed", "12345"], False)
check("small_prime", ["-r", "[3]*3", "7", "--seed", "12345"], False)
check("extension", ["-r", "[3]*3", "2^8", "--seed", "12345"], True)
# Rejected genericity should safely fall back even after the core path starts.
check("singular", ["x^3+y^3+z,x^3+y^3+x+1,x^3+y^3+y+2", "x,y", "65537"], True, True)

# Non-MQ automatic selection stays on the established Dixon path.
p = subprocess.run(["./drsolve", "-r", "[3]*3", "65537", "--seed", "12345",
                    "--threads", "1", "-v", "1"], env=env, capture_output=True,
                   text=True, timeout=30, check=True)
assert "grouped quotient recurrence" not in p.stdout
print("PASS non-MQ automatic selection unchanged")

solutions = []
for enabled in (False, True):
    path = output / f"solve_{enabled}.dr"
    p = subprocess.run(
        ["./drsolve", "-r", "[3]*3", "257", "--seed", "12345", "-s",
         "--quotient" if enabled else "--no-quotient", "--threads", "4", "-v", "3",
         "-o", str(path)],
        env=env, capture_output=True, text=True, timeout=30, check=True)
    assert ("grouped quotient recurrence" in p.stdout) == enabled
    solutions.append(path.read_text().split("=== Polynomial System Solutions ===", 1)[1]
                     .split("=== Solution Complete ===", 1)[0])
assert solutions[0] == solutions[1]
print("PASS solving uses grouped first resultant and preserves solutions")

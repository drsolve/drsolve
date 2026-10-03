"""Prime-field quotient solver CLI regression, with exhaustive point oracles."""
import itertools
import os
from pathlib import Path
import random
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
BIN = Path(os.environ.get("DRSOLVE_BIN", ROOT / "drsolve")).resolve()


def run(args, directory, ok=True, content=None):
    output = directory / "result.dr"
    if output.exists():
        output.unlink()
    if content is not None:
        source = directory / "input.dr"
        source.write_text(content)
        args = [*args, "-f", str(source)]
    proc = subprocess.run([str(BIN), "--threads", "1", *map(str, args),
                           "-o", str(output)], text=True, capture_output=True,
                          timeout=40, cwd=directory)
    text = proc.stdout + proc.stderr
    assert (proc.returncode == 0) == ok, (args, proc.returncode, text)
    saved = output.read_text() if output.exists() else ""
    return text, saved


def points(saved, names):
    result = []
    for block in re.split(r"Solution set \d+:\n", saved)[1:]:
        values = dict(re.findall(r"^  (\w+) = (\d+)\s*$", block, re.M))
        result.append(tuple(int(values[name]) for name in names))
    assert len(result) == len(set(result)), result
    return set(result)


def oracle(polys, names, prime):
    code = [compile(f.strip().replace("^", "**"), "fixture", "eval") for f in polys]
    return {a for a in itertools.product(range(prime), repeat=len(names))
            if all(eval(f, {"__builtins__": {}}, dict(zip(names, a))) % prime == 0
                   for f in code)}


def main():
    count = 0
    with tempfile.TemporaryDirectory(prefix="drsolve-quotient-cli-") as tmp:
        directory = Path(tmp)
        fixtures = [
            (["x^2-1", "y^2-1"], 7),
            (["x^2-x", "y^2-y", "x*y"], 2),
            (["x^2-x", "y^2-y", "0", "x^2-x"], 7),
            (["x^2", "(y-x)^2"], 3),
            (["x^2", "(y-x)^2"], 2),
            (["x^2", "x*y", "y^2"], 7),
            (["x^2+1", "y"], 3),  # algebraic solutions, no rational points
            (["x^2-1", "y", "1"], 7),
            (["x^3-x", "y^3-y"], 7),
            (["x^3-y", "y^2-1", "x^6-1"], 7),
            (["x+y-1", "x-y"], 7),
        ]
        for polys, p in fixtures:
            text, saved = run(["--quotient", ",".join(polys), p], directory)
            assert "Solver: quotient algebra" in text
            assert "certified dimension" in saved
            assert points(saved, ("x", "y")) == oracle(polys, ("x", "y"), p)
            count += 1

        # Both routing modes and declared random variables matter independently
        # of the number of variables that occur in a generated polynomial.
        text, saved = run(["x^2-1,y-x,x*y-1", 257], directory)
        assert "Solver: quotient algebra" in text and len(points(saved, ("x", "y"))) == 2
        text, _ = run(["x^2-1,y-x", 257], directory)
        assert "Solver: Dixon elimination" in text and "quotient algebra" not in text
        text, _ = run(["--solver", "dixon", "x^2-1,y-x,x*y-1", 257], directory)
        assert "Solver: Dixon elimination" in text and "Solver: quotient" not in text
        text, saved = run([], directory, content="7\nx^2-1\ny-x\nx*y-1\n")
        assert "Solver: quotient algebra" in text and len(points(saved, ("x", "y"))) == 2
        text, _ = run(["--quotient", "x+y", 7], directory, ok=False)
        assert "at least as many equations" in text
        text, _ = run(["--quotient", "x*y,x^2*y", 7, "--quotient-max-degree", 3], directory, ok=False)
        assert "incomplete" in text and "No solutions over" not in text
        text, _ = run(["--quotient", "x^4-1,y^4-1", 7, "--quotient-max-degree", 2], directory, ok=False)
        assert "degree limit reached" in text
        text, _ = run(["--quotient", "x^4-1,y^4-1,z^4-1", 7, "--quotient-memory", 1], directory, ok=False)
        assert "memory budget reached" in text
        for field in (0, "2^3", "18446744073709551629"):
            text, _ = run(["--quotient", "x^2-1,y-x", field], directory, ok=False)
            assert "machine-word prime fields" in text
        for option in ("--quotient-max-degree", "--quotient-memory"):
            run([option, 0, "--quotient", "x^2-1,y-x", 7], directory, ok=False)
        run(["--solver", "invalid", "x^2-1,y-x", 7], directory, ok=False)
        text, saved = run(["--quotient", "--silent", "x^2-1,y-x", 7], directory)
        assert not text and len(points(saved, ("x", "y"))) == 2

        for args in (["-r", "-n", 2, "[2]*4", 7, "--seed", 91],
                     ["-r", "--vardeg", "-m", 4, "[1,1]", 7, "--seed", 91]):
            text, saved = run(args, directory)
            assert "Solver: quotient algebra (4 equations, 2 variables)" in text
            polynomials = re.search(r"^Polynomials: (.*)$", saved, re.M)[1].split(",")
            assert points(saved, ("x0", "x1")) == oracle(polynomials, ("x0", "x1"), 7)
            _, again = run(args, directory)
            assert re.search(r"^Polynomials: (.*)$", saved, re.M)[1] == re.search(r"^Polynomials: (.*)$", again, re.M)[1]
        text, _ = run(["-r", "-n", 2, "[2]*3", 7, "--density", 0, "--homogeneous",
                       "--quotient-max-degree", 3, "--seed", 8], directory, ok=False)
        assert "3 equations, 2 variables" in text and "incomplete" in text
        run(["-r", "-n", 2, "[2]*4", 7, "--resultant-only"], directory, ok=False)
        run(["-r", "-n", 2, "[2]*4", 7, "--comp"], directory, ok=False)
        # Old eight-variable/six-degree experiment caps must not leak into production.
        names = tuple(f"x{i}" for i in range(9))
        fs = [f"{x}-{i}" for i, x in enumerate(names)] + ["x0+x1-1"]
        _, saved = run([",".join(fs), 257], directory)
        assert points(saved, names) == {tuple(range(9))}
        _, saved = run(["--quotient", "x^9-1,x^3-1", 7], directory)
        assert points(saved, ("x",)) == {(1,), (2,), (4,)}

        # Independent brute-force oracles, including characteristic two,
        # inconsistent systems, arbitrary fibers and redundant generators.
        for p in (2, 3, 7):
            for seed in range(8):
                rng = random.Random(100*p + seed)
                fs = [f"x^2+{rng.randrange(p)}*x+{rng.randrange(p)}",
                      f"y^2+{rng.randrange(p)}*y+{rng.randrange(p)}"]
                fs += ["+".join(f"{rng.randrange(p)}*{term}" for term in
                               ("x^2", "x*y", "y^2", "x", "y", "1"))]
                _, saved = run([",".join(fs), p], directory)
                assert points(saved, ("x", "y")) == oracle(fs, ("x", "y"), p)
                count += 1
        print(f"PASS: {count} exhaustive point comparisons; routing, limits, file/random modes and CLI errors")


if __name__ == "__main__":
    main()

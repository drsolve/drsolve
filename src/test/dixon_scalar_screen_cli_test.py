"""Exercise screening indicator output and reproducible, unplanted random generation."""
import pathlib
import re
import subprocess
import sys
import tempfile


binary = pathlib.Path(sys.argv[1]).resolve()


def run(directory, name, args, code=None):
    path = directory / (name + ".dr")
    proc = subprocess.run(
        [str(binary), *args, "-o", str(path)],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60,
    )
    if code is not None:
        assert proc.returncode == code, (name, proc.returncode, proc.stdout)
    else:
        assert proc.returncode in (0, 2), (name, proc.returncode, proc.stdout)
    return proc.stdout, path.read_text() if path.exists() else ""


def polynomials(saved):
    return next(line for line in saved.splitlines() if line.startswith("Polynomials:"))


def nonconstant_terms(saved):
    # The fixtures are linear, so each term is either a constant or c*x_i.
    coefficients = []
    for p in polynomials(saved).translate(str.maketrans("", "", " ()")).split(","):
        terms = {}
        for coefficient, variable in re.findall(r"([+-]?\d*)\*?(x\d+)", p):
            if coefficient in ("", "+", "-"):
                coefficient += "1"
            terms[variable] = int(coefficient) % 257
        coefficients.append(terms)
    return coefficients


with tempfile.TemporaryDirectory(prefix="drsolve-scalar-cli-") as work:
    directory = pathlib.Path(work)
    for method in ([], ["--method", "5"]):
        stdout, saved = run(directory, "no_roots", method + ["x,y,3", "x,y", "7"], 0)
        assert "Step 4: Dixon decision result" in stdout and "Dixon screening value = 1" in stdout
        assert "Dixon scalar screen: no_common_zero" in saved
        stdout, saved = run(directory, "common_root", method + ["x^2-1,y-x,x*y-1", "x,y", "7"], 2)
        assert "Step 4: Dixon decision result" in stdout and "Dixon screening value = 0" in stdout
        assert "Dixon scalar screen: inconclusive" in saved
        stdout, saved = run(directory, "zero", method + ["x,y,x+y", "x,y", "7"], 2)
        assert "Step 4:" in stdout and "Dixon scalar screen: inconclusive" in saved

    stdout, saved = run(directory, "silent", ["-v", "0", "--resultant-only", "x,y,3", "x,y", "7"], 0)
    assert stdout.strip() == "Dixon scalar screen: no_common_zero"

    args = ["-r", "[1]*3", "-n", "2", "257", "--seed", "123"]
    stdout, planted = run(directory, "planted", args, 2)
    assert "Planted solution:" in stdout and "Planted solution:" in planted
    stdout, pure = run(directory, "pure", args + ["--pure-random"])
    assert "Planted solution:" not in stdout and "Planted solution:" not in pure
    assert "pure random (no planted solution)" in stdout
    assert "pure random (no planted solution)" in pure
    _, repeated = run(directory, "repeat", args + ["--pure-random"])
    assert polynomials(pure) == polynomials(repeated)
    assert polynomials(pure) != polynomials(planted)
    assert nonconstant_terms(pure) == nonconstant_terms(planted)

    for extra in (["--quotient"], ["-s"], ["--vardeg", "-m", "3"]):
        args = ["-r", "[1]*2", "257", "--seed", "123", "--pure-random", *extra]
        stdout, saved = run(directory, "other_mode", args)
        assert "Planted solution:" not in saved
        assert "pure random (no planted solution)" in saved

    stdout, _ = run(directory, "invalid", ["--pure-random", "x,y,1", "x,y", "257"], 1)
    assert "may only be used together with --random" in stdout

print("Screening indicator and pure-random CLI tests passed")

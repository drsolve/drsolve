# Input and Coefficient Fields

[Documentation index](README.md)

## Command-line syntax

```text
drsolve "polynomials" "elimination_variables" FIELD [options]
drsolve "polynomials" FIELD [options]
drsolve [-f] FILE [-o OUTPUT] [options]
```

Separate polynomials and elimination variables with commas, and quote each list. Use `*` for multiplication and `^` for nonnegative integer powers, as in `x^2+2*x*y-1`. Rational coefficients can be written as `1/2`. A typical elimination task uses n+1 equations to eliminate n variables; the remaining variables are parameters.

```bash
./drsolve "x^2+y, x+y+1" "x" 257
./drsolve "x^2-1, y-x" 0
```

## Coefficient fields

| FIELD | Meaning | Notes |
| --- | --- | --- |
| `0` | Rationals Q | Exact rational arithmetic and modular reconstruction; approximate real roots are available |
| `257` | Prime field F_257 | Integer coefficients are reduced modulo the prime |
| An arbitrarily large prime p | Prime field F_p | Primes beyond the machine-word range use large-integer modular arithmetic |
| `2^8`, `3^2` | Extension field F_(p^k) | The prime characteristic p must fit the machine-word backend; this is not the integer ring modulo p^k |
| `2^8:t^8+t^4+t^3+t+1` | Extension with an explicit defining polynomial | Fixes the representation using the specified irreducible polynomial |

In extension fields, `t` denotes the field generator. Over Q or a prime field, it can be an ordinary variable. When exchanging extension-field elements with another program, use the same defining polynomial and generator representation. Specifying only the field size does not ensure identical representations across programs.

```bash
./drsolve "x+t, y+t^2" "2^8:t^8+t^4+t^3+t+1"
./drsolve "x^2-1, y-x" 18446744073709551629
```

For `--specialize-vars`, integer values in extension fields encode generator coefficients: in binary extensions, `2=t`, `3=1+t`, and `4=t^2`. This differs from ordinary integer constants in polynomial expressions.

## `.dr` files

A solver file starts with FIELD, followed by polynomials:

```text
257
x^2-1,
y-x
```

An elimination file starts with elimination variables, then FIELD:

```text
x
257
x^2+y,
x+y+1
```

Polynomials may span lines or be comma-separated. A digit at the start of the first data line selects solver format. For compatibility, n equations accompanied by n elimination variables retain only the first n−1 listed variables.

An optional comment such as `# options: --method 5 --fq-det-method iter` supplies supported file options; not every CLI option can be embedded. With `--ideal`, lines containing `=` after the headers can supply substitution relations.

```bash
./drsolve -f examples/example.dr -o out/elimination.dr
./drsolve -s -f examples/example_solve.dr -o out/solutions.dr
```

## Output

Without `-o`, results normally go to `out/` under the working directory. Inline input uses `solution_TIMESTAMP.dr`, or `comp_TIMESTAMP.dr` for complexity analysis. File input derives a name with a `_solution` or `_comp` suffix from the input filename. If automatic directory creation fails, the program may fall back to the working directory. Use distinct explicit paths for batch jobs to avoid overwrites.

Output includes the input, coefficient field, and results. Read solution status, complete solution sets, and verification information together; candidate sets are not automatically verified complete solutions. Approximate and exact results are identified separately.

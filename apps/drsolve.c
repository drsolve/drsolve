/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <flint/flint.h>

#include "drsolve_cli.h"

#ifndef PROGRAM_VERSION
#define PROGRAM_VERSION "unknown"
#endif
#define DEFAULT_OUTPUT_DIR "out"

/* =========================================================================
 * Print usage
 * ========================================================================= */

void drsolve_cli_print_version(void)
{
    printf("===============================================\n");
    printf("DRSolve v%s\n", PROGRAM_VERSION);
    printf("FLINT version: %s (Recommended: 3.6.0)\n", FLINT_VERSION);
#ifdef HAVE_PML
    printf("PML support: ENABLED\n");
#else
    printf("PML support: DISABLED\n");
#endif
    printf("===============================================\n");
}

void drsolve_cli_print_usage(const char *prog_name)
{
    printf("USAGE:\n");
    printf("  %s \"polynomials\" \"eliminate_vars\" field_size\n", prog_name);
    printf("  %s \"polynomials\" field_size\n", prog_name);
    printf("  %s input_file -o output_file\n", prog_name);
    printf("FILE FORMAT:\n");
    printf("  Dixon resultant elimination:\n");
    printf("    Line 1 : variables TO ELIMINATE (comma-separated)\n");
    printf("    Line 2 : field size (prime or p^k, 0 means rational)\n");
    printf("    Line 3+: polynomials (comma-separated, may span multiple lines)\n");
    printf("  Solver mode:\n");
    printf("    Line 1 : field size (prime or p^k, 0 means rational)\n");
    printf("    Line 2+: polynomials (comma-separated, may span multiple lines)\n");
    printf("OPTIONS:\n");
    printf("  -r \"[d1,d2,...,dn]\" random polynomial generation\n");
    printf("  -s  solving mode (auto-enables when no vars given)\n");
    printf("  -c  complexity analysis mode\n");
    printf("EXAMPLES:\n");
    printf("  Dixon resultant elimination:\n");
    printf("    %s \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    %s \"x^2+y^2+z^2-1, x^2+y^2-2*z^2, x+y+z\" \"x,y\" 0\n", prog_name);
    printf("  Polynomial system solving:\n");
    printf("    %s \"x^3+y^2+z-8, x+y+z-6, x*y*z-6\" 0\n", prog_name);
    printf("    %s \"x^2 + t*y, x*y + t^2\" \"2^8: t^8 + t^4 + t^3 + t + 1\"\n", prog_name);
    printf("  Random input:\n");
    printf("    %s -r \"[3]*4\" 257\n", prog_name);
    printf("    %s -r -s \"[2,3,4]\" 2^8 --seed 1234\n", prog_name);
    printf("  Complexity analysis:\n");
    printf("    %s -c \"x^2+y^2+1, x*y+z, x+y+z^2\" \"x,y\" 257\n", prog_name);
    printf("    %s -c -r \"[10]*10\" 257\n", prog_name);
    printf("  File input:\n");
    printf("    %s examples/example.dr\n", prog_name);
    printf("    %s examples/example_solve.dr -o solution.dr\n", prog_name);
    printf("OTHER OPTIONS:\n");
    printf("  --method <n>      Determinant method selection (0:Expansion, 1:HNF, 2:Interpolation, 3:Sparse, 4:Bareiss, 5:Fdixon)\n");
    printf("  --step1, --step4  Override method <n> for specific algorithm steps\n");
    printf("  --quotient        Use Dixon quotient / Bezoutian construction (automatic for supported MQ)\n");
    printf("  --sift            Remove extraneous resultant factors using quotient algebra (finite fields)\n");
    printf("  --threads <num>   Set number of threads for parallel computation\n");
    printf("  --dixon           Use Dixon resultant (default)\n");
    printf("  --macaulay        Use Macaulay resultant\n");
    printf("  --subres          Use Subresultant (2 polys)\n");
    printf("  --field-equation  After each multiplication, reduces x^q -> x for every variable\n");
    printf("  --ideal <args>    After each multiplication, reduces using the given substitution\n");
    printf("  --complex         Output complex solutions (2x2 solver or complex roots over Q)\n");
    printf("  --max-primes <n>  Maximum primes for rational reconstruction (default: 64)\n");
    printf("  --test <n>        Run built-in tests (1: Dixon matrix size, 2: Bezout bound, 3: solver correctness, 4: performance)\n");
    printf("  --time            Print per-step timing information\n");
    printf("  -v, --verbose <n> Verbosity level (0:silent, 1:default, 2:detailed, 3:debug, 4:matrix)\n");
    printf("  -h, --help        Show this help\n");
    printf("  --help-full       Show all options, defaults, file format and examples\n");
    printf("  -V, --version     Print version and build information\n");
}

void drsolve_cli_print_full_usage(const char *prog_name)
{
    printf("USAGE:\n"
           "  %s \"polynomials\" \"eliminate_vars\" FIELD [options]\n"
           "  %s \"polynomials\" FIELD [options]\n"
           "  %s [-f] FILE [-o OUTPUT] [options]\n"
           "  %s -r \"degree_spec\" FIELD [options]\n\n",
           prog_name, prog_name, prog_name, prog_name);
    puts("INPUT AND OUTPUT:");
    puts("  Polynomials are comma-separated; quote polynomial and variable lists.");
    puts("  FIELD: 0 for Q, a prime p for F_p, or p^k for an extension field.");
    puts("  An explicit modulus is accepted, e.g. \"2^8: t^8+t^4+t^3+t+1\".");
    puts("  In extensions, t is the field generator; in Q/F_p it is a variable.");
    puts("  -f FILE                       Read input (or give FILE as the sole operand).");
    puts("  -o FILE                       Override the output filename.");
    puts("  Inline output: out/solution_YYYYMMDD_HHMMSS.dr (comp_... for complexity).");
    puts("  File output: input_file_solution.dr (input_file_comp.dr for complexity).");
    puts("  Examples:");
    printf("    %s \"x^2-1, y-x\" 257 -o solutions.dr\n", prog_name);
    printf("    %s \"x^2+t*y, x*y+t^2\" \"2^8: t^8+t^4+t^3+t+1\"\n", prog_name);
    puts("");
    puts("FILE FORMAT:");
    puts("  Solver:      line 1 FIELD; remaining lines polynomials.");
    puts("  Elimination: line 1 elimination variables; line 2 FIELD; then polynomials.");
    puts("  Polynomials may span lines or be comma-separated. A digit on the first");
    puts("  data line selects solver format. For n equations with n listed elimination");
    puts("  variables, compatibility mode retains only the first n-1 variables.");
    puts("  Optional comment: # options: --method 5 --fq-det-method iter");
    puts("  With --ideal, lines containing '=' after the headers supply relations.");
    puts("  Examples:");
    printf("    %s -f examples/example.dr -o elimination.dr\n", prog_name);
    printf("    %s -s -f examples/example_solve.dr -o solutions.dr\n", prog_name);
    puts("");
    puts("SOLVING:");
    puts("  -s, --solve                    Solve the system; automatic without elim vars.");
    puts("  --solver auto|dixon|quotient    Select solver and enable solving (default auto).");
    puts("                                Auto uses Dixon for 2 equations in 2 variables, including back substitution.");
    puts("  --solve-rational-only          Solve over Q, keeping exact rational solutions.");
    puts("  --specialize-vars \"z=0,w=1\"    Fix variables before solving. In extensions,");
    puts("                                integers encode generator coefficients:");
    puts("                                in F_2^k, 2=t, 3=1+t, 4=t^2.");
    puts("  --quotient-max-degree D        Quotient closure degree budget (0 = unlimited).");
    puts("  --quotient-memory MiB          Quotient memory budget (0 = unlimited).");
    puts("  Both budgets default to 0. Reaching a limit makes auto fall back to Dixon;");
    puts("  explicit --solver quotient reports an incomplete solve.");
    puts("  Quotient solving requires at least as many equations as variables. Supported");
    puts("  square MQ systems use Bezoutian construction; overdetermined systems use");
    puts("  quotient closure. --no-quotient or nonzero budgets retain closure for square");
    puts("  systems. Random n+1 equations in n variables default to scalar consistency");
    puts("  screening; use -s to request solving.");
    puts("  --sift / --no-sift             Enable/disable extraneous-factor removal using");
    puts("                                quotient algebra over finite fields (default off).");
    puts("  --sift-max-degree D            Positive sift degree limit (default 4).");
    puts("  Examples:");
    printf("    %s --solver quotient \"x^2-1, y-x\" 257\n", prog_name);
    printf("    %s --solve-rational-only \"x^2-1, y-x\" 0\n", prog_name);
    printf("    %s -s --specialize-vars \"z=0,w=1\" \"x+z, y+w\" 257\n", prog_name);
    puts("");
    puts("RESULTANTS AND DETERMINANTS:");
    puts("  --dixon                       Dixon resultant (default).");
    puts("  --pml-fft / --no-pml-fft      Enable supported FFT matrix backend (default) / fall back.");
    puts("  --inverse-block / --no-inverse-block / --inverse-block-auto");
    puts("                                Dixon inverse projection (nmod); auto for 3/4 equations");
    puts("                                with one remaining variable. Failure uses the existing backend.");
    puts("  --inverse-block-size N        Set block size and enable; 0 selects automatically.");
    puts("  --bezout-block / --no-bezout-block");
    puts("                                nmod Dixon/Bezout block-Hankel path (default on).");
    puts("                                Randomized verification; fallback on failure.");
    puts("  --bezout-fast / --no-bezout-fast");
    puts("                                Dixon/Bezout fast pair resultant (default auto):");
    puts("                                NTT or multipoint interpolation over finite fields.");
    puts("                                Explicit resultant/determinant methods override auto.");
    puts("                                --bivar-ntt / --no-bivar-ntt are compatibility aliases.");
    puts("  --macaulay                    Macaulay resultant.");
    puts("  --subres                      Subresultant: 2 polynomials, 1 elimination var.");
    puts("  --resultant NAME              dixon|macaulay|subres (alias --resultant-method).");
    puts("  --method N                    Set Step 1/4 methods (5 selects recursive Dixon).");
    puts("  --step1 N / --step4 N          Override the corresponding step.");
    puts("    0 minor expansion; 1 HNF; 2 interpolation; 3 sparse interpolation;");
    puts("    4 Bareiss; 5 recursive Dixon construction; 6 balanced split Laplace");
    puts("    (experimental). Applicability depends on the step and coefficient field.");
    puts("  Finite-field elimination of 2 variables from 3 equations in 3 variables");
    puts("  defaults to method 5; explicit method options take precedence.");
    puts("  --fq-det-method MODE          auto|hnf|iter|interp (default auto): univariate");
    puts("                                determinant backend; interp is prime-field only.");
    puts("  --quotient                    Enable quotient/Bezoutian resultant construction;");
    puts("                                unlike --solver quotient, does not select solving.");
    puts("                                Alias: --mq-quotient.");
    puts("  --no-quotient                  Disable it (alias --no-mq-quotient).");
    puts("  Automatic for supported prime/extension-field square MQ systems (3+ vars).");
    puts("  Explicit --quotient also supports general degrees with one parameter;");
    puts("  rejected inputs or rank/degree checks fall back to Dixon. Samples can use");
    puts("  --threads. Resultant construction has no quotient memory budget.");
    puts("  Examples:");
    printf("    %s --subres \"x^2+y, x+y+1\" \"x\" 257\n", prog_name);
    printf("    %s --step1 0 --step4 1 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    printf("    %s --quotient -r \"[3]*3\" 65537 --seed 1234\n", prog_name);
    puts("");
    puts("ROOTS AND RECONSTRUCTION:");
    puts("  --resultant-only               Output the resultant without root finding");
    puts("                                (alias --no-roots).");
    puts("  --approx-roots                Use approximate real-root finding over Q.");
    puts("  --root-precision BITS         Precision, default 128; range 32..1048576.");
    puts("                                Also enables --approx-roots.");
    puts("  --complex                     Complex solutions for the Q-only 2x2 solver;");
    puts("                                complex roots for univariate Q resultants.");
    puts("  --rational-root-scan MODE     auto|off|force: exhaustive rational-root scan");
    puts("                                before approximate real-root finding.");
    puts("  --no-rational-root-scan       Alias for --rational-root-scan off.");
    puts("  --force-rational-root-scan    Alias for --rational-root-scan force.");
    printf("  Default auto skips scans above %d candidates; force has a hard cap of %d.\n",
           FMPQ_ROOT_SEARCH_AUTO_MAX_CANDIDATES,
           FMPQ_ROOT_SEARCH_HARD_MAX_CANDIDATES);
    puts("  --max-primes N                Reconstruction prime budget, 1..1048576");
    puts("                                (default 64 for rational reconstruction).");
    puts("  --large-prime-reconstruct     Opt in to reconstruction for primes exceeding");
    puts("                                the nmod limit; incompatible with --subres.");
    puts("  Examples:");
    printf("    %s --resultant-only \"x+y, x-y^2+2\" \"x\" 0\n", prog_name);
    printf("    %s --approx-roots --root-precision 256 \"x+y, x-y^2+2\" \"x\" 0\n", prog_name);
    printf("    %s --complex \"x^2+y^2-1, x-y\" 0\n", prog_name);
    puts("");
    puts("RANDOM INPUT:");
    puts("  -r, --random                  Generate polynomials from a degree specification:");
    puts("                                \"[2,3,4]\", \"[3]*4\", or \"[2]*4+[3]*2\".");
    puts("  -n, --num-vars N              Total variables (alias --nvars).");
    puts("  --vardeg                      Interpret degrees as per-variable exponent bounds;");
    puts("                                default: n equations in n variables.");
    puts("  -m, --num-equations M         Equation count with --vardeg (alias --m).");
    puts("  --density R                   Fraction of monomials, 0..1; default 0.5 over F2,");
    puts("                                otherwise 1.");
    puts("  --homogeneous                 Use only terms of the requested total degree");
    puts("                                (alias --hom); incompatible with --vardeg.");
    puts("  --seed N                      Reproducible random input.");
    puts("  --pure-random                 Do not plant a guaranteed solution.");
    puts("  Fewer variables than equations auto-enables solving. Overdetermined random");
    puts("  input normally contains a planted solution. Resultant/complexity mode needs");
    puts("  at least (#equations - 1) variables. Combine -r with -s or -c as needed.");
    puts("  Examples:");
    printf("    %s -r -s \"[2,3,4]\" 257 --seed 1234\n", prog_name);
    printf("    %s -r --homogeneous --density 0.5 \"[3]*3\" 257\n", prog_name);
    printf("    %s -r --vardeg -m 3 \"[2,3]\" 257 --seed 1234\n", prog_name);
    puts("");
    puts("COMPLEXITY AND REDUCTION:");
    puts("  -c, --comp                    Complexity analysis (alias --complexity).");
    printf("  -w, --omega VALUE             Matrix exponent (default %.4g).\n", DIXON_OMEGA);
    puts("  --rank-pred                   Experimental rank prediction (default off;");
    puts("                                alias --rank-prediction).");
    puts("  --ideal \"a^3=b+1,...\"         Reduce by substitution relations after products.");
    puts("                                With file input, read relations from the file.");
    puts("  --field-equation              Reduce x^q -> x after each multiplication.");
    puts("  --field-equation-s            Apply x^q -> x only to the final output.");
    puts("  Examples:");
    printf("    %s -c -r --omega 2.81 \"[2]*3\" 257\n", prog_name);
    printf("    %s --field-equation \"x*y+z, y*z+x, x*z+y\" \"x,y\" 2\n", prog_name);
    printf("    %s --ideal \"a^3=b+1\" \"a^3+b, a+b^2\" \"a\" 257\n", prog_name);
    puts("");
    puts("PERFORMANCE AND DIAGNOSTICS:");
    puts("  --threads N                   Positive parallel thread count; otherwise use");
    puts("                                the runtime default.");
    puts("  --cache N                     Live minor-DP entry limit for method 0");
    puts("                                (default 1024; 0 disables DP).");
    puts("  --array-limit-k K             Cap extension-field multiplication tables at");
    puts("                                2^K entries, K in 0..62.");
    puts("  -v, --verbose N               Output level (default 1):");
    puts("                                0 silent; 1 normal; 2 detailed with timing;");
    puts("                                3 profiling; 4 matrix dumps (up to 100 x 100)");
    puts("                                and recursive fast-Dixon traces.");
    puts("  --silent / --debug            Aliases for -v 0 / -v 2.");
    puts("  --solve-verbose               Enable solving with -v 2.");
    puts("  --time                        Per-step timing; interpolation shows CPU/wall time");
    puts("                                and threads.");
    puts("  --test N                      Built-in checks: 0 help; 1 matrix size;");
    puts("                                2 Bezout bound; 3 solver correctness;");
    puts("                                4 performance; 5 XHash.");
    puts("  --test-solver                 Built-in solver correctness check.");
    puts("  -h, --help                    Show the quick reference.");
    puts("  --help-full                   Show this complete help.");
    puts("  -V, --version                 Show version and build information.");
    puts("  Examples:");
    printf("    %s --threads 4 --time -r \"[2]*3\" 257 --seed 1234\n", prog_name);
    printf("    %s -v 4 \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
    puts("");
    puts("ADVANCED AND EXPERIMENTAL:");
    puts("  The following optimizations apply only to eligible inputs.");
    puts("  --mq-step1-rank / --no-mq-step1-rank");
    puts("                                Projected MQ direct indices (default on).");
    puts("  --mq-step1-shared / --no-mq-step1-shared");
    puts("                                Shared-index Step 1 DP (default on).");
    puts("  --no-mq-step1-filter           Disable Step 1 coefficient filtering (default on)");
    puts("                                for single-parameter prime/extension-field MQ.");
    puts("  --mq-step1-pencil / --no-mq-step1-pencil");
    puts("                                Experimental degree recurrence (default off).");
    puts("  --mq-step1-simplex / --no-mq-step1-simplex");
    puts("                                Experimental total-degree interpolation (off).");
    puts("  --step4-schur / --no-step4-schur");
    puts("                                Enable/disable checked Step 4 Schur compression.");
    puts("                                Default: direct inverse-block first; Schur on fallback");
    puts("                                or when inverse-block is off. Any degree.");
    puts("  --no-schur-seed                Disable complement-first minor selection");
    puts("                                (default on), retaining Step 4 compression.");
    puts("  --fast-ksy / --no-fast-ksy      KSY precondition check for method 5 (default off;");
    puts("                                enable alias --ksy-precondition).");
    puts("  --fast-ksy-col I               KSY constant-column index (default 0).");
    puts("  --step3-verify-second / --no-step3-verify-second");
    puts("                                Second Step 3 verification pass (default off).");
    puts("  Examples:");
    printf("    %s --no-quotient --no-step4-schur -r \"[2]*3\" 257 --seed 1234\n", prog_name);
    printf("    %s --method 5 --fast-ksy --step3-verify-second \"x+y+z, x*y+y*z+z*x, x*y*z+1\" \"x,y\" 257\n", prog_name);
}

static const char *display_prog_name(const char *argv0)
{
    const char *env_name = getenv("DIXON_DISPLAY_NAME");
    if (env_name && env_name[0] != '\0') return env_name;
    return argv0;
}

static int validate_cli_options(int argc, char *argv[])
{
    enum { OPT_FLAG = 1000 };
    static const struct option options[] = {
        {"silent", no_argument, NULL, OPT_FLAG},
        {"solve-verbose", no_argument, NULL, OPT_FLAG},
        {"solve-rational-only", no_argument, NULL, OPT_FLAG},
        {"complex", no_argument, NULL, OPT_FLAG},
        {"resultant-only", no_argument, NULL, OPT_FLAG},
        {"no-roots", no_argument, NULL, OPT_FLAG},
        {"approx-roots", no_argument, NULL, OPT_FLAG},
        {"root-precision", required_argument, NULL, OPT_FLAG},
        {"solve", no_argument, NULL, 's'},
        {"sift", no_argument, NULL, OPT_FLAG},
        {"no-sift", no_argument, NULL, OPT_FLAG},
        {"sift-max-degree", required_argument, NULL, OPT_FLAG},
        {"quotient", no_argument, NULL, OPT_FLAG},
        {"solver", required_argument, NULL, OPT_FLAG},
        {"quotient-max-degree", required_argument, NULL, OPT_FLAG},
        {"quotient-memory", required_argument, NULL, OPT_FLAG},
        {"comp", no_argument, NULL, 'c'},
        {"complexity", no_argument, NULL, 'c'},
        {"rank-pred", no_argument, NULL, OPT_FLAG},
        {"rank-prediction", no_argument, NULL, OPT_FLAG},
        {"random", no_argument, NULL, 'r'},
        {"vardeg", no_argument, NULL, OPT_FLAG},
        {"homogeneous", no_argument, NULL, OPT_FLAG},
        {"hom", no_argument, NULL, OPT_FLAG},
        {"ideal", no_argument, NULL, OPT_FLAG},
        {"field-equation", no_argument, NULL, OPT_FLAG},
        {"field-equation-s", no_argument, NULL, OPT_FLAG},
        {"time", no_argument, NULL, OPT_FLAG},
        {"no-rational-root-scan", no_argument, NULL, OPT_FLAG},
        {"force-rational-root-scan", no_argument, NULL, OPT_FLAG},
        {"rational-root-scan", required_argument, NULL, OPT_FLAG},
        {"debug", no_argument, NULL, OPT_FLAG},
        {"verbose", required_argument, NULL, 'v'},
        {"dixon", no_argument, NULL, OPT_FLAG},
        {"macaulay", no_argument, NULL, OPT_FLAG},
        {"subres", no_argument, NULL, OPT_FLAG},
        {"resultant", required_argument, NULL, OPT_FLAG},
        {"resultant-method", required_argument, NULL, OPT_FLAG},
        {"omega", required_argument, NULL, 'w'},
        {"method", required_argument, NULL, OPT_FLAG},
        {"fq-det-method", required_argument, NULL, OPT_FLAG},
        {"fast-ksy", no_argument, NULL, OPT_FLAG},
        {"ksy-precondition", no_argument, NULL, OPT_FLAG},
        {"fast-ksy-col", required_argument, NULL, OPT_FLAG},
        {"no-fast-ksy", no_argument, NULL, OPT_FLAG},
        {"step3-verify-second", no_argument, NULL, OPT_FLAG},
        {"no-step3-verify-second", no_argument, NULL, OPT_FLAG},
        {"mq-quotient", no_argument, NULL, OPT_FLAG},
        {"no-quotient", no_argument, NULL, OPT_FLAG},
        {"no-mq-quotient", no_argument, NULL, OPT_FLAG},
        {"mq-step1-rank", no_argument, NULL, OPT_FLAG},
        {"no-mq-step1-rank", no_argument, NULL, OPT_FLAG},
        {"mq-step1-shared", no_argument, NULL, OPT_FLAG},
        {"no-mq-step1-shared", no_argument, NULL, OPT_FLAG},
        {"mq-step1-pencil", no_argument, NULL, OPT_FLAG},
        {"no-mq-step1-pencil", no_argument, NULL, OPT_FLAG},
        {"mq-step1-simplex", no_argument, NULL, OPT_FLAG},
        {"no-mq-step1-simplex", no_argument, NULL, OPT_FLAG},
        {"no-mq-step1-filter", no_argument, NULL, OPT_FLAG},
        {"step4-schur", no_argument, NULL, OPT_FLAG},
        {"no-step4-schur", no_argument, NULL, OPT_FLAG},
        /* Compatibility aliases; new help uses the degree-independent names. */
        {"mq-step4-schur", no_argument, NULL, OPT_FLAG},
        {"no-mq-step4-schur", no_argument, NULL, OPT_FLAG},
        {"pml-fft", no_argument, NULL, OPT_FLAG},
        {"no-pml-fft", no_argument, NULL, OPT_FLAG},
        {"inverse-block", no_argument, NULL, OPT_FLAG},
        {"no-inverse-block", no_argument, NULL, OPT_FLAG},
        {"inverse-block-auto", no_argument, NULL, OPT_FLAG},
        {"inverse-block-size", required_argument, NULL, OPT_FLAG},
        {"bezout-block", no_argument, NULL, OPT_FLAG},
        {"no-bezout-block", no_argument, NULL, OPT_FLAG},
        {"bezout-fast", no_argument, NULL, OPT_FLAG},
        {"no-bezout-fast", no_argument, NULL, OPT_FLAG},
        {"no-bivar-ntt", no_argument, NULL, OPT_FLAG},
        {"bivar-ntt", no_argument, NULL, OPT_FLAG},
        {"no-schur-seed", no_argument, NULL, OPT_FLAG},
        {"step1", required_argument, NULL, OPT_FLAG},
        {"step4", required_argument, NULL, OPT_FLAG},
        {"threads", required_argument, NULL, OPT_FLAG},
        {"max-primes", required_argument, NULL, OPT_FLAG},
        {"large-prime-reconstruct", no_argument, NULL, OPT_FLAG},
        {"cache", required_argument, NULL, OPT_FLAG},
        {"array-limit-k", required_argument, NULL, OPT_FLAG},
        {"nvars", required_argument, NULL, 'n'},
        {"num-vars", required_argument, NULL, 'n'},
        {"m", required_argument, NULL, OPT_FLAG},
        {"num-equations", required_argument, NULL, OPT_FLAG},
        {"density", required_argument, NULL, OPT_FLAG},
        {"seed", required_argument, NULL, OPT_FLAG},
        {"pure-random", no_argument, NULL, OPT_FLAG},
        {"specialize-vars", required_argument, NULL, OPT_FLAG},
        {"help", no_argument, NULL, 'h'},
        {"help-full", no_argument, NULL, OPT_FLAG},
        {"version", no_argument, NULL, 'V'},
        {"test", no_argument, NULL, OPT_FLAG},
        {"test-solver", no_argument, NULL, OPT_FLAG},
        {NULL, 0, NULL, 0}
    };
    char *args[argc + 1];
    int option_index = 0;
    int parsed;

    memcpy(args, argv, (size_t) argc * sizeof(*args));
    args[argc] = NULL;
    optind = 1;
    opterr = 0;

    while ((parsed = getopt_long(argc, args, "-scrv:w:n:m:f:o:hV", options,
                                  &option_index)) != -1) {
        if (parsed == '?') {
            const char *bad = (optind > 0 && optind <= argc) ? args[optind - 1] : "";
            fprintf(stderr, "Error: Unknown option or missing argument '%s'. Use --help-full to list supported options.\n",
                    bad);
            optind = 1;
            return 0;
        }
    }

    optind = 1;
    return 1;
}

int main(int argc, char *argv[])
{
    const char *prog_name = display_prog_name(argv[0]);

    if (argc == 1) {
        drsolve_cli_print_version();
        drsolve_cli_print_usage(prog_name);
        return 0;
    }

    if (argc == 2 &&
        (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0 ||
         strcmp(argv[1], "--help-full") == 0)) {
        drsolve_cli_print_version();
        if (strcmp(argv[1], "--help-full") == 0)
            drsolve_cli_print_full_usage(prog_name);
        else
            drsolve_cli_print_usage(prog_name);
        return 0;
    }

    if (!validate_cli_options(argc, argv)) {
        return 1;
    }

    if (argc == 2 &&
        (strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-V") == 0)) {
        drsolve_cli_print_version();
        return 0;
    }

    return drsolve_cli_main(argc, argv, prog_name);
}

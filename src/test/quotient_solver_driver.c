/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Read the standalone benchmark format; exercise the production backend.
 * JSON output is used by the independent Sage ideal audit. */
#include "quotient_solver.h"

int main(int argc, char **argv)
{
    if (argc != 4) return 2;
    FILE *fp = fopen(argv[1], "r");
    if (!fp) return 2;
    ulong p;
    slong n, m;
    if (fscanf(fp, "%lu %ld %ld\n", &p, &n, &m) != 3 || n < 1 || m < 1) return 2;
    fmpz_t prime;
    fmpz_init_set_ui(prime, p);
    fq_nmod_ctx_t ctx;
    fq_nmod_ctx_init(ctx, prime, 1, "a");
    fmpz_clear(prime);
    variable_info_t *vars = calloc(n, sizeof(*vars));
    char **fs = calloc(m, sizeof(char *));
    for (slong i = 0; i < n; i++) {
        vars[i].name = malloc(32);
        snprintf(vars[i].name, 32, "x%ld", i);
    }
    for (slong i = 0; i < m; i++) {
        size_t cap = 0;
        if (getline(fs+i, &cap, fp) < 0) return 2;
        fs[i][strcspn(fs[i], "\r\n")] = '\0';
    }
    fclose(fp);
    polynomial_solutions_t sols;
    polynomial_solutions_init(&sols, n, ctx);
    int ok = solve_by_quotient_closure(fs, m, vars, n, &sols, atol(argv[2]), atol(argv[3]));
    printf("{\"certified\":%s", ok ? "true" : "false");
    if (ok) {
        const char *s = strstr(sols.elimination_summary, "certified dimension ");
        slong dim = atol(s + strlen("certified dimension "));
        const char *mu = strstr(sols.resultant_steps[0], "): ") + 3;
        printf(",\"dimension\":%ld,\"eliminant\":\"%s\",\"points\":[", dim, mu);
        for (slong i = 0; i < sols.num_solution_sets; i++) {
            printf("%s[", i ? "," : "");
            for (slong j = 0; j < n; j++)
                printf("%s%lu", j ? "," : "", nmod_poly_get_coeff_ui(sols.solution_sets[i][j][0], 0));
            putchar(']');
        }
        putchar(']');
    }
    puts("}");
    polynomial_solutions_clear(&sols);
    for (slong i = 0; i < n; i++) free(vars[i].name);
    for (slong i = 0; i < m; i++) free(fs[i]);
    free(vars); free(fs);
    fq_nmod_ctx_clear(ctx);
    flint_cleanup_master();
    return 0;
}

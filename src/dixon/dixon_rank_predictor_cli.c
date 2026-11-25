/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "dixon_rank_predictor.h"
#include <errno.h>
#include <stdlib.h>
#include <string.h>

static slong positive(const char *s)
{
    char *end;
    errno = 0;
    long v = strtol(s, &end, 10);
    if (errno || *end || v <= 0 || v >= WORD_MAX) {
        fprintf(stderr, "Invalid positive integer: %s\n", s);
        exit(2);
    }
    return v;
}

int main(int argc, char **argv)
{
    int uniform = argc > 1 && !strcmp(argv[1], "--uniform");
    if (argc < 3 || (uniform && argc != 4)) {
        fprintf(stderr, "Usage: %s d1 d2 ... | --uniform n d\n", argv[0]);
        return 2;
    }
    slong n = uniform ? positive(argv[2]) : argc - 1;
    if (n < 2) {
        fputs("Need at least two degrees\n", stderr);
        return 2;
    }
    long *degrees = flint_calloc(n, sizeof(*degrees));
    for (slong i = 0; i < n; ++i)
        degrees[i] = positive(argv[uniform ? 3 : i + 1]);

    fmpz_t rank;
    fmpz_init(rank);
    int ok = dixon_generic_rank_from_degrees(rank, degrees, n, stdout);
    if (!ok)
        fputs("Degree index exceeds machine range\n", stderr);
    fmpz_clear(rank);
    flint_free(degrees);
    flint_cleanup();
    return ok ? 0 : 2;
}

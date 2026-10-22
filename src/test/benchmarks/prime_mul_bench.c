/* SPDX-License-Identifier: GPL-2.0-or-later */
/* cc -O3 -march=native -Isrc -Ipml_det -Ipml_det/src -I/usr/local/include
 * src/test/benchmarks/prime_mul_bench.c -L. -ldrsolve -L/usr/local/lib
 * -lflint -lgmp -lm -o /tmp/prime_mul_bench
 * LD_LIBRARY_PATH=.:/usr/local/lib OMP_NUM_THREADS=1 taskset -c 2 /tmp/prime_mul_bench
 * CSV medians of 7 warmed samples, >= 20 ms each; setup/conversion excluded.
 */
#include "unified_mat_poly.h"
#include <flint/fmpz_mod_poly.h>
#include <flint/nmod_poly.h>
#include <assert.h>
#include <time.h>

#define BATCH 1024
static volatile ulong sink;
static double now(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}
static int compare(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}
typedef void (*kernel)(void *, ulong);
static void measure(const char *p, const char *data, const char *op, slong length, slong batch,
                    kernel fn, void *ctx)
{
    ulong repeats = 1;
    double elapsed, samples[7];
    do {
        double start = now();
        fn(ctx, repeats);
        elapsed = now() - start;
        if (elapsed < 0.02)
            repeats *= 2;
    } while (elapsed < 0.02);
    for (int i = 0; i < 7; i++) {
        double start = now();
        fn(ctx, repeats);
        samples[i] = (now() - start) * 1e9 / (repeats * (double)batch);
    }
    qsort(samples, 7, sizeof(double), compare);
    printf("%s,%s,%s,%ld,%.3f,%.3f,%.3f,%lu\n", p, data, op, length, samples[3], samples[0],
           samples[6], repeats);
    fflush(stdout);
}
typedef struct {
    nmod_t mod;
    ulong a[BATCH], b[BATCH], c[BATCH];
    fmpz x[BATCH], y[BATCH], z[BATCH];
    fmpz_mod_ctx_t ctx;
} scalars;
static void scalar_nmod(void *arg, ulong repeats)
{
    scalars *s = arg;
    for (ulong r = 0; r < repeats; r++) {
        for (int j = 0; j < BATCH; j++)
            s->c[j] = nmod_mul(s->a[j], s->b[j], s->mod);
        __asm__ volatile("" ::: "memory");
    }
    sink = s->c[17];
}
static void scalar_fmpz(void *arg, ulong repeats)
{
    scalars *s = arg;
    for (ulong r = 0; r < repeats; r++) {
        for (int j = 0; j < BATCH; j++)
            fmpz_mod_mul(s->z + j, s->x + j, s->y + j, s->ctx);
        __asm__ volatile("" ::: "memory");
    }
    sink = fmpz_get_ui(s->z + 17);
}
typedef struct {
    nmod_poly_t a, b, c;
    fmpz_mod_poly_t x, y, z;
    fmpz_mod_ctx_t ctx;
    unified_mat_poly_t ua, ub, uc;
} polys;
static void poly_nmod(void *arg, ulong repeats)
{
    polys *s = arg;
    for (ulong r = 0; r < repeats; r++)
        nmod_poly_mul(s->c, s->a, s->b);
    sink = nmod_poly_get_coeff_ui(s->c, 0);
}
static void poly_fmpz(void *arg, ulong repeats)
{
    polys *s = arg;
    for (ulong r = 0; r < repeats; r++)
        fmpz_mod_poly_mul(s->z, s->x, s->y, s->ctx);
    sink = fmpz_get_ui(s->z->coeffs);
}
static void poly_unified(void *arg, ulong repeats)
{
    polys *s = arg;
    for (ulong r = 0; r < repeats; r++)
        unified_mat_poly_mul(s->uc, s->ua, s->ub);
    sink = s->uc->length;
}
static void run(const char *prime, int small, flint_rand_t rng)
{
    fmpz_t p, v;
    fmpz_init(p);
    fmpz_init(v);
    fmpz_set_str(p, prime, 10);
    int limb = fmpz_abs_fits_ui(p);
    const char *data = small ? "8bit" : "uniform";
    scalars s;
    fmpz_mod_ctx_init(s.ctx, p);
    if (limb)
        nmod_init(&s.mod, fmpz_get_ui(p));
    for (int j = 0; j < BATCH; j++) {
        fmpz_init(s.x + j);
        fmpz_init(s.y + j);
        fmpz_init(s.z + j);
        fmpz_randm(s.x + j, rng, p);
        fmpz_randm(s.y + j, rng, p);
        if (small) {
            fmpz_set_ui(s.x + j, 1 + fmpz_fdiv_ui(s.x + j, 255));
            fmpz_set_ui(s.y + j, 1 + fmpz_fdiv_ui(s.y + j, 255));
        }
        s.a[j] = fmpz_get_ui(s.x + j);
        s.b[j] = fmpz_get_ui(s.y + j);
    }
    scalar_fmpz(&s, 1);
    if (limb) {
        scalar_nmod(&s, 1);
        for (int j = 0; j < BATCH; j++)
            assert(fmpz_equal_ui(s.z + j, s.c[j]));
        measure(prime, data, "nmod_mul", 1, BATCH, scalar_nmod, &s);
    }
    measure(prime, data, "fmpz_mod_mul", 1, BATCH, scalar_fmpz, &s);
    for (int j = 0; j < BATCH; j++) {
        fmpz_clear(s.x + j);
        fmpz_clear(s.y + j);
        fmpz_clear(s.z + j);
    }
    fmpz_mod_ctx_clear(s.ctx);

    const slong lengths[] = {4, 16, 64, 256, 1024};
    field_ctx_t field;
    field_ctx_init_prime(&field, p);
    for (size_t k = 0; k < sizeof(lengths) / sizeof(*lengths); k++) {
        slong n = lengths[k];
        polys q;
        fmpz_mod_ctx_init(q.ctx, p);
        nmod_poly_init(q.a, limb ? fmpz_get_ui(p) : 257);
        nmod_poly_init(q.b, limb ? fmpz_get_ui(p) : 257);
        nmod_poly_init(q.c, limb ? fmpz_get_ui(p) : 257);
        fmpz_mod_poly_init(q.x, q.ctx);
        fmpz_mod_poly_init(q.y, q.ctx);
        fmpz_mod_poly_init(q.z, q.ctx);
        unified_poly_mat_t u[2];
        field_elem_u coeff;
        field_init_elem(&coeff, field.field_id, field_ctx_data(&field));
        for (int side = 0; side < 2; side++) {
            unified_poly_mat_init(u[side], 1, 1, &field);
            for (slong j = 0; j < n; j++) {
                fmpz_randm(v, rng, p);
                if (small)
                    fmpz_set_ui(v, 1 + fmpz_fdiv_ui(v, 255));
                if (j == n - 1 && fmpz_is_zero(v))
                    fmpz_one(v);
                fmpz_mod_poly_set_coeff_fmpz(side ? q.y : q.x, j, v, q.ctx);
                if (limb) {
                    nmod_poly_set_coeff_ui(side ? q.b : q.a, j, fmpz_get_ui(v));
                    coeff.nmod = fmpz_get_ui(v);
                } else
                    fmpz_set(&coeff.big, v);
                unified_poly_set_coeff(unified_poly_mat_entry(u[side], 0, 0), j, &coeff);
            }
        }
        unified_mat_poly_init(q.ua, 1, 1, &field);
        unified_mat_poly_init(q.ub, 1, 1, &field);
        unified_mat_poly_init(q.uc, 1, 1, &field);
        unified_mat_poly_from_unified(q.ua, u[0]);
        unified_mat_poly_from_unified(q.ub, u[1]);
        poly_fmpz(&q, 1);
        poly_unified(&q, 1);
        unified_poly_mat_t check;
        unified_poly_mat_init(check, 1, 1, &field);
        unified_mat_poly_to_unified(check, q.uc);
        for (slong j = 0; j < 2 * n - 1; j++) {
            unified_poly_get_coeff(&coeff, unified_poly_mat_entry(check, 0, 0), j);
            fmpz_mod_poly_get_coeff_fmpz(v, q.z, j, q.ctx);
            assert(limb ? fmpz_equal_ui(v, coeff.nmod) : fmpz_equal(v, &coeff.big));
        }
        if (limb) {
            poly_nmod(&q, 1);
            for (slong j = 0; j < 2 * n - 1; j++) {
                fmpz_mod_poly_get_coeff_fmpz(v, q.z, j, q.ctx);
                assert(fmpz_equal_ui(v, nmod_poly_get_coeff_ui(q.c, j)));
            }
            measure(prime, data, "nmod_poly_mul", n, 1, poly_nmod, &q);
        }
        measure(prime, data, "fmpz_mod_poly_mul", n, 1, poly_fmpz, &q);
        measure(prime, data, "unified_mat_poly_mul_1x1", n, 1, poly_unified, &q);
        unified_poly_mat_clear(check);
        unified_poly_mat_clear(u[0]);
        unified_poly_mat_clear(u[1]);
        unified_mat_poly_clear(q.ua);
        unified_mat_poly_clear(q.ub);
        unified_mat_poly_clear(q.uc);
        field_clear_elem(&coeff, field.field_id, field_ctx_data(&field));
        nmod_poly_clear(q.a);
        nmod_poly_clear(q.b);
        nmod_poly_clear(q.c);
        fmpz_mod_poly_clear(q.x, q.ctx);
        fmpz_mod_poly_clear(q.y, q.ctx);
        fmpz_mod_poly_clear(q.z, q.ctx);
        fmpz_mod_ctx_clear(q.ctx);
    }
    field_ctx_clear(&field);
    fmpz_clear(v);
    fmpz_clear(p);
}
int main(void)
{
    flint_set_num_threads(1);
    flint_rand_t rng;
    flint_rand_init(rng);
    flint_rand_set_seed(rng, 12345, 67890);
    puts("prime,data,operation,length,median_ns,min_ns,max_ns,repeats");
    const char *primes[] = {"257", "2305843009213693951", "18446744073709551557",
                            "18446744073709551629"};
    for (int small = 0; small < 2; small++)
        for (int i = 0; i < 4; i++)
            run(primes[i], small, rng);
    flint_rand_clear(rng);
    return 0;
}

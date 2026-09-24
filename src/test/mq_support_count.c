/* SPDX-License-Identifier: GPL-2.0-or-later */
#define main mq_support_fixture_main
#include "dixon_mq_filter_test.c"
#undef main
int mq_support_count(fq_mvpoly_t **matrix, slong n, const slong *rows,
                     const slong *cols, slong targets);

/* Coefficients are placeholders only: identical equations deliberately give
 * no useful determinant. Only their generic MQ term supports are analyzed. */
static fq_mvpoly_t *dense_support_fixture(slong n, const fq_nmod_ctx_t ctx, flint_rand_t state)
{
    fq_mvpoly_t *p = flint_malloc((size_t) (n + 1) * sizeof(*p));
    slong *e = flint_calloc((size_t) n + 1, sizeof(slong));
    fq_nmod_t c;
    fq_nmod_init(c, ctx);
    for (slong k = 0; k <= n; k++) {
        fq_mvpoly_init(p + k, n, 1, ctx);
        for (slong i = -1; i <= n; i++) {
            for (slong j = i; j <= n; j++) {
                if (i >= 0) e[i]++;
                if (j >= 0) e[j]++;
                ulong coeff = 1; /* Placeholder: all generic input monomials present. */
                /* Ensure elimination degree two even over F_2. */
                if (i == 0 && j == 0) coeff = 1;
                fq_nmod_set_ui(c, coeff, ctx);
                if (coeff) fq_mvpoly_add_term_fast(p + k, e, e + n, c);
                if (i >= 0) e[i]--;
                if (j >= 0) e[j]--;
            }
        }
    }
    fq_nmod_clear(c, ctx); flint_free(e);
    return p;
}

static int count_candidate(fq_mvpoly_t **matrix, slong nvars)
{
    long *degrees=flint_malloc((nvars+1)*sizeof(long));
    for(slong i=0;i<=nvars;i++) degrees[i]=2;
    slong *R = NULL, *H = NULL, rlen = 0, hlen = 0, sigma = 0, rank = 0;
    int model = dixon_rank_profile_from_degrees(&R, &rlen, &H, &hlen,
                                  &sigma, &rank, degrees, nvars + 1, nvars);

    slong count = 0;
    for (slong i = 0; model && i < rlen; i++) {
        /* Bound support construction before allocating/enumerating. */
        if (R[i] > (1L << 18) - count) model = 0;
        else count += R[i];
    }
    flint_free(R);
    if (!model || rank <= 0 || rank >= count) { flint_free(H); flint_free(degrees); return 0; }
    slong *exps = flint_malloc((size_t) count * nvars * sizeof(slong));
    slong *reverse = flint_malloc((size_t) count * nvars * sizeof(slong));
    slong *tmp = flint_calloc((size_t) nvars, sizeof(slong));
    slong actual = 0;
    if (!dixon_mq_support(exps, count, &actual, tmp, nvars, 0, nvars) || actual != count) {
        flint_free(tmp); flint_free(reverse); flint_free(exps); flint_free(H);
        flint_free(degrees); return 0;
    }
    monom_t *rm = NULL, *cm = NULL;
    slong nr = 0, nc = 0, rcap = 0, ccap = 0, rhs = 16, chs = 16;
    hash_entry_t **ri = flint_calloc(16, sizeof(*ri));
    hash_entry_t **ci = flint_calloc(16, sizeof(*ci));
    /* Canonical MQ staircase: squarefree monomials first, then increasing
     * repeated-variable degree, breaking ties by descending lex order.
     * Unlike the legacy first-occurrence order this exists before Step 1.
     * Its candidate must pass verification; no almost-revlex theorem is
     * assumed for random inputs or small characteristic. */
    for (slong excess = 0; excess <= nvars; excess++) {
        for (slong j = 0; j < count; j++) {
            slong e = 0;
            for (slong v = 0; v < nvars; v++) e += FLINT_MAX(0, exps[j * nvars + v] - 1);
            if (e != excess) continue;
            for (slong v = 0; v < nvars; v++)
                reverse[j * nvars + v] = exps[j * nvars + nvars - 1 - v];
            dixon_intern_monom(&rm, &nr, &rcap, &ri, &rhs, reverse + j * nvars, nvars);
            dixon_intern_monom(&cm, &nc, &ccap, &ci, &chs, exps + j * nvars, nvars);
        }
    }
    slong *rows = NULL, *cols = NULL, size = 0;
    int ok = dixon_build_predicted_mirror_indices(&rows, &cols, &size,
                    rm, nr, cm, nc, ri, rhs, ci, chs, nvars, H, hlen, sigma, 2, rank);
    flint_free(H);
    slong *target_rows = NULL, *target_cols = NULL;

    assert(ok && size==rank);
    target_rows=flint_malloc(rank*nvars*sizeof(slong));
    target_cols=flint_malloc(rank*nvars*sizeof(slong));
    for(slong i=0;i<rank;i++) {
        memcpy(target_rows+i*nvars,rm[rows[i]].exp,nvars*sizeof(slong));
        memcpy(target_cols+i*nvars,cm[cols[i]].exp,nvars*sizeof(slong));
    }
    int counted=mq_support_count(matrix,nvars+1,target_rows,target_cols,rank);
    flint_free(target_rows); flint_free(target_cols);
    flint_free(rows); flint_free(cols);
    free_monom_index(ri,rhs); free_monom_index(ci,chs);
    flint_free(rm); flint_free(cm);
    flint_free(tmp); flint_free(reverse); flint_free(exps); flint_free(degrees);
    return counted;
}
int main(int argc,char **argv)
{
    if(argc!=2) return 2;
    slong n=atol(argv[1]); if(n<4 || n>8) return 2;
    omp_set_num_threads(1); flint_set_num_threads(1); g_dixon_verbose_level=0;
    flint_rand_t rng; flint_rand_init(rng); flint_rand_set_seed(rng,132,941);
    fq_nmod_ctx_t fq; fq_nmod_ctx_init_ui(fq,65537,1,"a");
    fq_mvpoly_t *p=dense_support_fixture(n-1,fq,rng), **m, **a;
    build_fq_cancellation_matrix_mvpoly(&m,p,n-1,1);
    perform_fq_matrix_row_operations_mvpoly(&a,&m,n-1,1);
    int ok=count_candidate(a,n-1);
    clear_input(p,m,a,n-1); fq_nmod_ctx_clear(fq); flint_rand_clear(rng);
    flint_cleanup_master(); return !ok;
}

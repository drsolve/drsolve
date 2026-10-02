/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "dr_mpoly.h"
#include "zech_mpoly_access.h"
#include <pthread.h>

/* Share immutable ring contexts across the matrix. A ring cannot outlive the
 * last polynomial using it; the caller's field context must outlive its polys. */
typedef struct dr_mpoly_ring {
    const fq_nmod_ctx_struct *field;
    slong dimensions, references;
    field_ctx_t backend;
    unified_mpoly_ctx_t context;
    struct dr_mpoly_ring *next;
} dr_mpoly_ring;
static dr_mpoly_ring *rings;
static pthread_mutex_t ring_lock = PTHREAD_MUTEX_INITIALIZER;

static dr_mpoly_ring *ring_acquire(const fq_nmod_ctx_t field, slong dimensions)
{
    pthread_mutex_lock(&ring_lock);
    dr_mpoly_ring *r;
    for (r = rings; r; r = r->next)
        if (r->field == field && r->dimensions == dimensions)
            break;
    if (!r) {
        r = flint_calloc(1, sizeof(*r));
        r->field = field;
        r->dimensions = dimensions;
        /* Retain unified's binary-field and small-extension backend selection. */
        field_ctx_init_enhanced(&r->backend, field, 1L << 20);
        r->context = unified_mpoly_ctx_init(dimensions, ORD_LEX, &r->backend);
        r->next = rings;
        rings = r;
    }
    r->references++;
    pthread_mutex_unlock(&ring_lock);
    return r;
}

static void ring_release(dr_mpoly_ring *r)
{
    pthread_mutex_lock(&ring_lock);
    if (--r->references == 0) {
        dr_mpoly_ring **link = &rings;
        while (*link != r)
            link = &(*link)->next;
        *link = r->next;
        unified_mpoly_ctx_clear(r->context);
        field_ctx_clear(&r->backend);
        flint_free(r);
    }
    pthread_mutex_unlock(&ring_lock);
}

int g_field_equation_reduction;
int g_field_equation_final_only;
void dr_mpoly_set_field_equation_reduction(int enable) { g_field_equation_reduction = enable; }
void dr_mpoly_set_field_equation_final_only(int enable) { g_field_equation_final_only = enable; }

void dr_mpoly_clear(unified_mpoly_struct *p)
{
    if (p->ring) {
        if (p->field_id == FIELD_ID_NMOD)
            nmod_mpoly_clear(GET_NMOD_POLY(p), GET_NMOD_CTX(p->ctx_ptr));
        else if (p->field_id == FIELD_ID_FQ_ZECH)
            fq_zech_mpoly_clear(GET_ZECH_POLY(p), GET_ZECH_CTX(p->ctx_ptr));
        else
            fq_nmod_mpoly_clear(GET_FQ_POLY(p), GET_FQ_CTX(p->ctx_ptr));
        ring_release(p->ring);
    }
    memset(p, 0, sizeof(*p));
}

void dr_mpoly_init(unified_mpoly_struct *p, slong nvars, slong npars, const fq_nmod_ctx_t field)
{
    if (p->ring && p->ctx == field && p->nvars == nvars && p->npars == npars) {
        unified_mpoly_zero(p);
        p->canonical = 1;
        return;
    }
    dr_mpoly_ring *r = ring_acquire(field, nvars + npars);
    dr_mpoly_clear(p);
    p->ring = r;
    p->ctx_ptr = r->context;
    p->field_id = r->backend.field_id;
    p->nvars = nvars;
    p->npars = npars;
    p->ctx = field;
    p->canonical = 1;
    if (p->field_id == FIELD_ID_NMOD)
        nmod_mpoly_init(GET_NMOD_POLY(p), GET_NMOD_CTX(p->ctx_ptr));
    else if (p->field_id == FIELD_ID_FQ_ZECH)
        fq_zech_mpoly_init(GET_ZECH_POLY(p), GET_ZECH_CTX(p->ctx_ptr));
    else
        fq_nmod_mpoly_init(GET_FQ_POLY(p), GET_FQ_CTX(p->ctx_ptr));
}

void dr_mpoly_move(unified_mpoly_struct *out, unified_mpoly_struct *in)
{
    if (out == in)
        return;
    dr_mpoly_clear(out);
    *out = *in;
    memset(in, 0, sizeof(*in));
}

void dr_mpoly_fit_length(unified_mpoly_struct *p, slong length)
{
    if (p->field_id == FIELD_ID_NMOD)
        nmod_mpoly_fit_length(GET_NMOD_POLY(p), length, GET_NMOD_CTX(p->ctx_ptr));
    else if (p->field_id == FIELD_ID_FQ_ZECH)
        fq_zech_mpoly_fit_length(GET_ZECH_POLY(p), length, GET_ZECH_CTX(p->ctx_ptr));
    else
        fq_nmod_mpoly_fit_length(GET_FQ_POLY(p), length, GET_FQ_CTX(p->ctx_ptr));
}

void dr_mpoly_read_term(dr_mpoly_term_view *v, slong *exp, ulong *coefficient,
                        const unified_mpoly_struct *p, slong i)
{
    assert(i >= 0 && i < dr_mpoly_length(p));
    slong d;
    if (p->field_id == FIELD_ID_NMOD) {
        nmod_mpoly_get_term_exp_si(exp, &p->data.nmod_poly, i, GET_NMOD_CTX(p->ctx_ptr));
        v->coeff->coeffs = p->data.nmod_poly.coeffs + i;
        d = 1;
    } else if (p->field_id == FIELD_ID_FQ_ZECH) {
        fq_zech_mpoly_ctx_struct *ctx = GET_ZECH_CTX(p->ctx_ptr);
        dr_zech_get_exp(exp, GET_ZECH_POLY(p), i, ctx);
        d = fq_nmod_ctx_degree(p->ctx);
        v->coeff->coeffs = coefficient;
        v->coeff->alloc = d;
        v->coeff->length = 0;
        v->coeff->mod = p->ctx->modulus->mod;
        fq_zech_get_fq_nmod(v->coeff, p->data.zech_poly.coeffs + i, ctx->fqctx);
        d = v->coeff->length;
    } else {
        fq_nmod_mpoly_get_term_exp_si(exp, &p->data.fq_poly, i, GET_FQ_CTX(p->ctx_ptr));
        d = fq_nmod_ctx_degree(p->ctx);
        v->coeff->coeffs = p->data.fq_poly.coeffs + i * d;
    }
    v->coeff->alloc = d;
    v->coeff->length = d;
    v->coeff->mod = p->ctx->modulus->mod;
    while (v->coeff->length && !v->coeff->coeffs[v->coeff->length - 1])
        v->coeff->length--;
    v->var_exp = p->nvars ? exp : NULL;
    v->par_exp = p->npars ? exp + p->nvars : NULL;
}

void dr_mpoly_normalize(unified_mpoly_struct *p)
{
    if (!p->ctx_ptr || p->canonical)
        return;
    if (p->field_id == FIELD_ID_NMOD) {
        nmod_mpoly_sort_terms(GET_NMOD_POLY(p), GET_NMOD_CTX(p->ctx_ptr));
        nmod_mpoly_combine_like_terms(GET_NMOD_POLY(p), GET_NMOD_CTX(p->ctx_ptr));
    } else if (p->field_id == FIELD_ID_FQ_ZECH) {
        fq_zech_mpoly_sort_terms(GET_ZECH_POLY(p), GET_ZECH_CTX(p->ctx_ptr));
        fq_zech_mpoly_combine_like_terms(GET_ZECH_POLY(p), GET_ZECH_CTX(p->ctx_ptr));
    } else {
        fq_nmod_mpoly_sort_terms(GET_FQ_POLY(p), GET_FQ_CTX(p->ctx_ptr));
        fq_nmod_mpoly_combine_like_terms(GET_FQ_POLY(p), GET_FQ_CTX(p->ctx_ptr));
    }
    p->canonical = 1;
}

void dr_mpoly_add_term_fast(unified_mpoly_struct *p, const slong *vars, const slong *pars,
                            const fq_nmod_t coefficient)
{
    if (fq_nmod_is_zero(coefficient, p->ctx))
        return;
    ulong exp[FLINT_MAX(1, p->nvars + p->npars)];
    for (slong i = 0; i < p->nvars; i++)
        exp[i] = vars ? vars[i] : 0;
    for (slong i = 0; i < p->npars; i++)
        exp[p->nvars + i] = pars ? pars[i] : 0;
    if (p->field_id == FIELD_ID_NMOD)
        nmod_mpoly_push_term_ui_ui(GET_NMOD_POLY(p), nmod_poly_get_coeff_ui(coefficient, 0), exp,
                                   GET_NMOD_CTX(p->ctx_ptr));
    else if (p->field_id == FIELD_ID_FQ_ZECH) {
        fq_zech_t c;
        fq_zech_mpoly_ctx_struct *ctx = GET_ZECH_CTX(p->ctx_ptr);
        fq_zech_set_fq_nmod(c, coefficient, ctx->fqctx);
        dr_zech_push(GET_ZECH_POLY(p), c, exp, ctx);
    } else
        fq_nmod_mpoly_push_term_fq_nmod_ui(GET_FQ_POLY(p), coefficient, exp,
                                           GET_FQ_CTX(p->ctx_ptr));
    p->canonical = 0;
}

void dr_mpoly_add_term(unified_mpoly_struct *p, const slong *vars, const slong *pars,
                       const fq_nmod_t coefficient)
{
    if (fq_nmod_is_zero(coefficient, p->ctx))
        return;
    dr_mpoly_normalize(p);
    ulong exp[FLINT_MAX(1, p->nvars + p->npars)];
    for (slong i = 0; i < p->nvars; i++)
        exp[i] = vars ? vars[i] : 0;
    for (slong i = 0; i < p->npars; i++)
        exp[p->nvars + i] = pars ? pars[i] : 0;
    if (p->field_id == FIELD_ID_NMOD) {
        nmod_mpoly_struct *a = GET_NMOD_POLY(p);
        ulong c = nmod_mpoly_get_coeff_ui_ui(a, exp, GET_NMOD_CTX(p->ctx_ptr));
        c = nmod_add(c, nmod_poly_get_coeff_ui(coefficient, 0), p->ctx->modulus->mod);
        nmod_mpoly_set_coeff_ui_ui(a, c, exp, GET_NMOD_CTX(p->ctx_ptr));
    } else if (p->field_id == FIELD_ID_FQ_ZECH) {
        fq_zech_t c, value;
        fq_zech_mpoly_ctx_struct *ctx = GET_ZECH_CTX(p->ctx_ptr);
        fq_zech_set_fq_nmod(value, coefficient, ctx->fqctx);
        dr_zech_get_coeff(c, GET_ZECH_POLY(p), exp, ctx);
        fq_zech_add(c, c, value, ctx->fqctx);
        dr_zech_set_coeff(GET_ZECH_POLY(p), c, exp, ctx);
    } else {
        fq_nmod_t c;
        fq_nmod_init(c, p->ctx);
        fq_nmod_mpoly_get_coeff_fq_nmod_ui(c, GET_FQ_POLY(p), exp, GET_FQ_CTX(p->ctx_ptr));
        fq_nmod_add(c, c, coefficient, p->ctx);
        fq_nmod_mpoly_set_coeff_fq_nmod_ui(GET_FQ_POLY(p), c, exp, GET_FQ_CTX(p->ctx_ptr));
        fq_nmod_clear(c, p->ctx);
    }
}

void dr_mpoly_set_term_exp(unified_mpoly_struct *p, slong i, const slong *exp)
{
    if (p->field_id == FIELD_ID_NMOD)
        nmod_mpoly_set_term_exp_ui(GET_NMOD_POLY(p), i, (const ulong *)exp,
                                   GET_NMOD_CTX(p->ctx_ptr));
    else if (p->field_id == FIELD_ID_FQ_ZECH)
        dr_zech_set_exp(GET_ZECH_POLY(p), i, (const ulong *)exp, GET_ZECH_CTX(p->ctx_ptr));
    else
        fq_nmod_mpoly_set_term_exp_ui(GET_FQ_POLY(p), i, (const ulong *)exp,
                                      GET_FQ_CTX(p->ctx_ptr));
    p->canonical = 0;
}

/* Reuse output storage; FLINT arithmetic permits output/input aliasing. */
static void prepare_output(unified_mpoly_struct *out, const unified_mpoly_struct *in)
{
    if (out->ctx_ptr != in->ctx_ptr || out->nvars != in->nvars || out->npars != in->npars)
        dr_mpoly_init(out, in->nvars, in->npars, in->ctx);
}

void dr_mpoly_copy(unified_mpoly_struct *out, const unified_mpoly_struct *in)
{
    if (out == in)
        return;
    prepare_output(out, in);
    unified_mpoly_set(out, (unified_mpoly_t)in);
    out->canonical = in->canonical;
}

static void arithmetic(unified_mpoly_struct *out, const unified_mpoly_struct *a,
                       const unified_mpoly_struct *b, char op)
{
    assert(a->nvars == b->nvars && a->npars == b->npars && a->ctx == b->ctx);
    dr_mpoly_normalize((unified_mpoly_t)a);
    dr_mpoly_normalize((unified_mpoly_t)b);
    prepare_output(out, a);
    if (a->field_id == FIELD_ID_NMOD) {
        nmod_mpoly_ctx_struct *ctx = GET_NMOD_CTX(a->ctx_ptr);
        if (op == '+')
            nmod_mpoly_add(GET_NMOD_POLY(out), GET_NMOD_POLY(a), GET_NMOD_POLY(b), ctx);
        else if (op == '-')
            nmod_mpoly_sub(GET_NMOD_POLY(out), GET_NMOD_POLY(a), GET_NMOD_POLY(b), ctx);
        else
            nmod_mpoly_mul(GET_NMOD_POLY(out), GET_NMOD_POLY(a), GET_NMOD_POLY(b), ctx);
        out->canonical = 1;
        if (op == '*' && g_field_equation_reduction)
            dr_mpoly_reduce_field_equation(out);
    } else {
        if (op == '+')
            unified_mpoly_add(out, (unified_mpoly_t)a, (unified_mpoly_t)b);
        else if (op == '-')
            unified_mpoly_sub(out, (unified_mpoly_t)a, (unified_mpoly_t)b);
        else if (!unified_mpoly_mul(out, (unified_mpoly_t)a, (unified_mpoly_t)b))
            flint_throw(FLINT_ERROR, "multivariate multiplication failed");
        out->canonical = 1;
        if (op == '*' && out->field_id == FIELD_ID_FQ_ZECH && g_field_equation_reduction)
            dr_mpoly_reduce_field_equation(out);
    }
}
void dr_mpoly_add(unified_mpoly_struct *o, const unified_mpoly_struct *a,
                  const unified_mpoly_struct *b)
{
    arithmetic(o, a, b, '+');
}
void dr_mpoly_sub(unified_mpoly_struct *o, const unified_mpoly_struct *a,
                  const unified_mpoly_struct *b)
{
    arithmetic(o, a, b, '-');
}
void dr_mpoly_mul(unified_mpoly_struct *o, const unified_mpoly_struct *a,
                  const unified_mpoly_struct *b)
{
    arithmetic(o, a, b, '*');
}

void dr_mpoly_neg(unified_mpoly_struct *out, const unified_mpoly_struct *in)
{
    prepare_output(out, in);
    unified_mpoly_neg(out, (unified_mpoly_t)in);
    out->canonical = in->canonical;
}

void dr_mpoly_scalar_mul(unified_mpoly_struct *out, const unified_mpoly_struct *in,
                         const fq_nmod_t c)
{
    prepare_output(out, in);
    if (in->field_id == FIELD_ID_NMOD)
        nmod_mpoly_scalar_mul_ui(GET_NMOD_POLY(out), GET_NMOD_POLY(in),
                                 nmod_poly_get_coeff_ui(c, 0), GET_NMOD_CTX(in->ctx_ptr));
    else if (in->field_id == FIELD_ID_FQ_ZECH) {
        fq_zech_t scalar;
        fq_zech_mpoly_ctx_struct *ctx = GET_ZECH_CTX(in->ctx_ptr);
        fq_zech_set_fq_nmod(scalar, c, ctx->fqctx);
        fq_zech_mpoly_scalar_mul_fq_zech(GET_ZECH_POLY(out), GET_ZECH_POLY(in), scalar, ctx);
    } else
        fq_nmod_mpoly_scalar_mul_fq_nmod(GET_FQ_POLY(out), GET_FQ_POLY(in), c,
                                         GET_FQ_CTX(in->ctx_ptr));
    out->canonical = in->canonical;
}

void dr_mpoly_pow(unified_mpoly_struct *out, const unified_mpoly_struct *in, slong power)
{
    assert(power >= 0);
    unified_mpoly_struct a = {0}, r = {0};
    dr_mpoly_copy(&a, in);
    dr_mpoly_init(&r, in->nvars, in->npars, in->ctx);
    unified_mpoly_one(&r);
    while (power) {
        if (power & 1)
            dr_mpoly_mul(&r, &r, &a);
        power >>= 1;
        if (power)
            dr_mpoly_mul(&a, &a, &a);
    }
    dr_mpoly_clear(&a);
    dr_mpoly_move(out, &r);
}

void dr_mpoly_reduce_field_equation(unified_mpoly_struct *p)
{
    fmpz_t order;
    fmpz_init(order);
    fq_nmod_ctx_order(order, p->ctx);
    if (!fmpz_fits_si(order)) {
        fmpz_clear(order);
        return;
    }
    slong q = fmpz_get_si(order);
    fmpz_clear(order);
    for (slong i = 0; i < dr_mpoly_length(p); i++) {
        DR_MPOLY_TERM(t, p, i);
        for (slong j = 0; j < p->nvars + p->npars; j++)
            if (t_exponents[j] >= q)
                t_exponents[j] = 1 + (t_exponents[j] - 1) % (q - 1);
        dr_mpoly_set_term_exp(p, i, t_exponents);
    }
    dr_mpoly_normalize(p);
}

void dr_mpoly_divide_difference(unified_mpoly_struct *out, const unified_mpoly_struct *in,
                                slong variable, slong nvars, slong npars)
{
    unified_mpoly_struct divisor = {0}, q = {0};
    dr_mpoly_init(&divisor, nvars, npars, in->ctx);
    dr_mpoly_init(&q, nvars, npars, in->ctx);
    unified_mpoly_gen(&divisor, variable, divisor.ctx_ptr);
    unified_mpoly_gen(&q, variable + nvars / 2, q.ctx_ptr);
    unified_mpoly_sub(&divisor, &divisor, &q);
    dr_mpoly_normalize((unified_mpoly_t)in);
    if (!unified_mpoly_divides(&q, (unified_mpoly_t)in, &divisor))
        flint_throw(FLINT_ERROR, "Dixon difference division is not exact");
    dr_mpoly_clear(&divisor);
    dr_mpoly_move(out, &q);
}

void dr_mpoly_to_nmod_mpoly(nmod_mpoly_t out, const unified_mpoly_struct *in, nmod_mpoly_ctx_t ctx)
{
    assert(in->field_id == FIELD_ID_NMOD);
    nmod_mpoly_set(out, GET_NMOD_POLY(in), ctx);
    if (!in->canonical) {
        nmod_mpoly_sort_terms(out, ctx);
        nmod_mpoly_combine_like_terms(out, ctx);
    }
}
void nmod_mpoly_to_dr_mpoly(unified_mpoly_struct *out, const nmod_mpoly_t in, slong nv, slong np,
                            const nmod_mpoly_ctx_t ctx, const fq_nmod_ctx_t field)
{
    dr_mpoly_init(out, nv, np, field);
    nmod_mpoly_set(GET_NMOD_POLY(out), in, ctx);
}
void dr_mpoly_take_nmod(unified_mpoly_struct *out, nmod_mpoly_t in, slong nv, slong np,
                        const nmod_mpoly_ctx_t ctx, const fq_nmod_ctx_t field)
{
    dr_mpoly_init(out, nv, np, field);
    nmod_mpoly_swap(GET_NMOD_POLY(out), in, ctx);
}

void dr_mpoly_to_fq_nmod_mpoly(fq_nmod_mpoly_t out, const unified_mpoly_struct *in,
                               fq_nmod_mpoly_ctx_t ctx)
{
    if (in->field_id != FIELD_ID_NMOD && in->field_id != FIELD_ID_FQ_ZECH)
        fq_nmod_mpoly_set(out, GET_FQ_POLY(in), ctx);
    else {
        fq_nmod_mpoly_zero(out, ctx);
        fq_nmod_mpoly_fit_length(out, dr_mpoly_length(in), ctx);
        for (slong i = 0; i < dr_mpoly_length(in); i++) {
            DR_MPOLY_TERM(t, in, i);
            fq_nmod_mpoly_push_term_fq_nmod_ui(out, t.coeff, (ulong *)t_exponents, ctx);
        }
    }
    if (!in->canonical) {
        fq_nmod_mpoly_sort_terms(out, ctx);
        fq_nmod_mpoly_combine_like_terms(out, ctx);
    }
}
void fq_nmod_mpoly_to_dr_mpoly(unified_mpoly_struct *out, const fq_nmod_mpoly_t in, slong nv,
                               slong np, fq_nmod_mpoly_ctx_t ctx, const fq_nmod_ctx_t field)
{
    dr_mpoly_init(out, nv, np, field);
    if (out->field_id != FIELD_ID_NMOD && out->field_id != FIELD_ID_FQ_ZECH)
        fq_nmod_mpoly_set(GET_FQ_POLY(out), in, ctx);
    else {
        ulong exp[FLINT_MAX(1, nv + np)];
        fq_nmod_t c;
        fq_nmod_init(c, field);
        dr_mpoly_fit_length(out, in->length);
        for (slong i = 0; i < in->length; i++) {
            fq_nmod_mpoly_get_term_exp_ui(exp, in, i, ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(c, in, i, ctx);
            dr_mpoly_add_term_fast(out, (slong *)exp, (slong *)exp + nv, c);
        }
        fq_nmod_clear(c, field);
        out->canonical = 1;
    }
}

int compare_fq_degrees(const void *a, const void *b)
{
    fq_index_degree_pair *pa = (fq_index_degree_pair *)a;
    fq_index_degree_pair *pb = (fq_index_degree_pair *)b;

    /* Sort by degree (ascending) */
    if (pa->degree < pb->degree)
        return -1;
    if (pa->degree > pb->degree)
        return 1;

    /* Same degree, sort by index for stability */
    if (pa->index < pb->index)
        return -1;
    if (pa->index > pb->index)
        return 1;
    return 0;
}

void fq_nmod_print_with_parentheses(const fq_nmod_t a, const fq_nmod_ctx_t ctx,
                                    const char *gen_name)
{
    if (fq_nmod_is_zero(a, ctx)) {
        printf("0");
        return;
    }

    slong degree = fq_nmod_ctx_degree(ctx);

    if (degree == 1) {
        // Prime field - no parentheses needed
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, a, ctx);

        if (nmod_poly_degree(poly) >= 0) {
            printf("%lu", nmod_poly_get_coeff_ui(poly, 0));
        } else {
            printf("0");
        }
        nmod_poly_clear(poly);
    } else {
        // Extension field - ALWAYS use parentheses for multi-term expressions
        nmod_poly_t poly;
        nmod_poly_init(poly, fq_nmod_ctx_prime(ctx));
        fq_nmod_get_nmod_poly(poly, a, ctx);

        slong deg = nmod_poly_degree(poly);

        // Count non-zero terms to decide on parentheses
        int term_count = 0;
        for (slong i = deg; i >= 0; i--) {
            if (nmod_poly_get_coeff_ui(poly, i) != 0) {
                term_count++;
            }
        }

        // Use parentheses for multi-term expressions or when explicitly requested
        int use_parens = (term_count > 1);

        if (use_parens)
            printf("(");

        int first_term = 1;
        for (slong i = deg; i >= 0; i--) {
            mp_limb_t coeff = nmod_poly_get_coeff_ui(poly, i);
            if (coeff != 0) {
                if (!first_term) {
                    printf(" + ");
                }
                first_term = 0;

                if (i == 0) {
                    printf("%lu", coeff);
                } else if (i == 1) {
                    if (coeff == 1) {
                        printf("%s", gen_name ? gen_name : "t");
                    } else {
                        printf("%lu*%s", coeff, gen_name ? gen_name : "t");
                    }
                } else {
                    if (coeff == 1) {
                        printf("%s^%ld", gen_name ? gen_name : "t", i);
                    } else {
                        printf("%lu*%s^%ld", coeff, gen_name ? gen_name : "t", i);
                    }
                }
            }
        }

        if (first_term) {
            printf("0");
        }

        if (use_parens)
            printf(")");

        nmod_poly_clear(poly);
    }
}

void dr_mpoly_print_with_names(const unified_mpoly_struct *poly, const char *poly_name,
                               char **var_names, char **par_names, const char *gen_name,
                               int expanded_format)
{

    // Print polynomial name if provided
    if (poly_name && strlen(poly_name) > 0) {
        printf("%s = ", poly_name);
    }

    if (dr_mpoly_length(poly) == 0) {
        printf("0\n");
        return;
    }

    // Default names if not provided
    char default_var_names[] = {'x', 'y', 'z', 'w', 'v', 'u'};
    char default_par_names[] = {'a', 'b', 'c', 'd'};

    for (slong i = 0; i < dr_mpoly_length(poly); i++) {
        DR_MPOLY_TERM(term_15, poly, i);

        if (i > 0)
            printf(" + ");

        // Check if we have variables or parameters for this term
        int has_vars_or_pars = 0;

        // Check variables
        if (poly->nvars > 0 && term_15.var_exp) {
            for (slong j = 0; j < poly->nvars; j++) {
                if (term_15.var_exp[j] > 0) {
                    has_vars_or_pars = 1;
                    break;
                }
            }
        }

        // Check parameters
        if (!has_vars_or_pars && poly->npars > 0 && term_15.par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                if (term_15.par_exp[j] > 0) {
                    has_vars_or_pars = 1;
                    break;
                }
            }
        }

        // Handle coefficient
        fq_nmod_t one;
        fq_nmod_init(one, poly->ctx);
        fq_nmod_one(one, poly->ctx);

        if (fq_nmod_is_one(term_15.coeff, poly->ctx) && has_vars_or_pars) {
            // Coefficient is 1 and we have variables/parameters - don't print coefficient
        } else {
            // Print coefficient with proper parentheses
            fq_nmod_print_with_parentheses(term_15.coeff, poly->ctx, gen_name);

            if (has_vars_or_pars) {
                printf("*");
            }
        }

        fq_nmod_clear(one, poly->ctx);

        // Track whether we've printed anything for this term yet
        int term_printed = 0;

        // Print variables
        // Fix for dr_mpoly_print_with_names function
        // Replace the variable printing sections:

        // Print variables
        if (expanded_format && poly->nvars > 0 && poly->nvars % 2 == 0) {
            // Expanded format with dual variables (for Dixon polynomials)
            slong actual_nvars = poly->nvars / 2;

            // Regular variables
            for (slong j = 0; j < actual_nvars; j++) {
                if (term_15.var_exp && term_15.var_exp[j] > 0) {
                    // FIXED: Only add * if something was already printed for this term
                    if (term_printed) {
                        printf("*");
                    }

                    // Use provided name or default
                    if (var_names && var_names[j]) {
                        printf("%s", var_names[j]);
                    } else if (j < 6) {
                        printf("%c", default_var_names[j]);
                    } else {
                        printf("x_%ld", j);
                    }

                    if (term_15.var_exp[j] > 1) {
                        printf("^%ld", term_15.var_exp[j]);
                    }
                    term_printed = 1;
                }
            }

            // Dual variables with tilde
            for (slong j = actual_nvars; j < poly->nvars; j++) {
                if (term_15.var_exp && term_15.var_exp[j] > 0) {
                    // FIXED: Only add * if something was already printed for this term
                    if (term_printed) {
                        printf("*");
                    }
                    slong orig_idx = j - actual_nvars;

                    if (var_names && var_names[orig_idx]) {
                        printf("~%s", var_names[orig_idx]);
                    } else if (orig_idx < 6) {
                        printf("~%c", default_var_names[orig_idx]);
                    } else {
                        printf("~x_%ld", orig_idx);
                    }

                    if (term_15.var_exp[j] > 1) {
                        printf("^%ld", term_15.var_exp[j]);
                    }
                    term_printed = 1;
                }
            }
        } else {
            // Normal variable format
            for (slong j = 0; j < poly->nvars; j++) {
                if (term_15.var_exp && term_15.var_exp[j] > 0) {
                    // FIXED: Only add * if something was already printed for this term
                    if (term_printed) {
                        printf("*");
                    }

                    if (var_names && var_names[j]) {
                        printf("%s", var_names[j]);
                    } else if (j < 6) {
                        printf("%c", default_var_names[j]);
                    } else {
                        printf("x_%ld", j);
                    }

                    if (term_15.var_exp[j] > 1) {
                        printf("^%ld", term_15.var_exp[j]);
                    }
                    term_printed = 1;
                }
            }
        }

        // Print parameters with actual names
        for (slong j = 0; j < poly->npars; j++) {
            if (term_15.par_exp && term_15.par_exp[j] > 0) {
                // FIXED: Only add * if something was already printed for this term
                if (term_printed) {
                    printf("*");
                }

                // Use provided parameter name or default
                if (par_names && par_names[j]) {
                    printf("%s", par_names[j]);
                } else if (j < 4) {
                    printf("%c", default_par_names[j]);
                } else {
                    printf("p_%ld", j);
                }

                if (term_15.par_exp[j] > 1) {
                    printf("^%ld", term_15.par_exp[j]);
                }
                term_printed = 1;
            }
        }
    }
    printf("\n");
}

void dr_mpoly_print(const unified_mpoly_struct *poly, const char *name)
{
    dr_mpoly_print_with_names(poly, name, NULL, NULL, NULL, 0);
}

void dr_mpoly_print_expanded(const unified_mpoly_struct *poly, const char *name, int use_dual)
{
    dr_mpoly_print_with_names(poly, name, NULL, NULL, NULL, use_dual);
}

void dr_mpoly_matrix_print(unified_mpoly_struct **matrix, slong nrows, slong ncols,
                           const char *matrix_name, int show_details)
{
    printf("\n=== %s Matrix (%ld x %ld) ===\n", matrix_name, nrows, ncols);

    if (show_details) {
        for (slong i = 0; i < nrows; i++) {
            for (slong j = 0; j < ncols; j++) {
                printf("M[%ld][%ld]: ", i, j);
                if (dr_mpoly_length(&(matrix[i][j])) == 0) {
                    printf("0\n");
                } else if (dr_mpoly_length(&(matrix[i][j])) <= 10) {
                    dr_mpoly_print_expanded(&matrix[i][j], "", 1);
                } else {
                    printf("%ld terms", dr_mpoly_length(&(matrix[i][j])));

                    slong max_var_deg = 0, max_par_deg = 0;
                    for (slong t = 0; t < dr_mpoly_length(&(matrix[i][j])); t++) {
                        DR_MPOLY_TERM(term_21, &(matrix[i][j]), t);

                        if (term_21.var_exp) {
                            for (slong k = 0; k < matrix[i][j].nvars; k++) {
                                if (term_21.var_exp[k] > max_var_deg) {
                                    max_var_deg = term_21.var_exp[k];
                                }
                            }
                        }
                        if (term_21.par_exp && matrix[i][j].npars > 0) {
                            for (slong k = 0; k < matrix[i][j].npars; k++) {
                                if (term_21.par_exp[k] > max_par_deg) {
                                    max_par_deg = term_21.par_exp[k];
                                }
                            }
                        }
                    }
                    printf(" (max var deg: %ld, max par deg: %ld)\n", max_var_deg, max_par_deg);
                }
            }
        }
    } else {
        printf("Matrix term counts:\n");
        for (slong i = 0; i < nrows; i++) {
            printf("  Row %ld: ", i);
            for (slong j = 0; j < ncols; j++) {
                printf("%3ld", dr_mpoly_length(&(matrix[i][j])));
                if (j < ncols - 1)
                    printf(" ");
            }
            printf("\n");
        }
    }
    printf("=== End %s Matrix ===\n\n", matrix_name);
}

void evaluate_dr_mpoly_at_params(fq_nmod_t result, const unified_mpoly_struct *poly,
                                 const fq_nmod_t *param_vals)
{
    if (poly->field_id == FIELD_ID_NMOD) {
        ulong values[FLINT_MAX(1, poly->nvars + poly->npars)];
        for (slong i = 0; i < poly->nvars; i++)
            values[i] = 1;
        for (slong i = 0; i < poly->npars; i++)
            values[poly->nvars + i] = nmod_poly_get_coeff_ui(param_vals[i], 0);
        ulong value =
            nmod_mpoly_evaluate_all_ui(GET_NMOD_POLY(poly), values, GET_NMOD_CTX(poly->ctx_ptr));
        fq_nmod_set_ui(result, value, poly->ctx);
        return;
    }
    const field_ctx_t *unified_ctx = poly->ctx_ptr->field_ctx;

    // Get context pointer
    void *ctx_ptr = NULL;
    switch (unified_ctx->field_id) {
    case FIELD_ID_NMOD:
        ctx_ptr = (void *)&unified_ctx->ctx.nmod_ctx;
        break;
    case FIELD_ID_FQ_ZECH:
        ctx_ptr = (void *)unified_ctx->ctx.zech_ctx;
        break;
    default:
        ctx_ptr = (void *)unified_ctx->ctx.fq_ctx;
        break;
    }

    // Initialize unified result to zero
    field_elem_u unified_result;
    field_init_elem(&unified_result, unified_ctx->field_id, ctx_ptr);
    field_set_zero(&unified_result, unified_ctx->field_id, ctx_ptr);

    // Convert parameters to unified format
    field_elem_u *unified_params = NULL;
    if (poly->npars > 0) {
        unified_params = (field_elem_u *)malloc(poly->npars * sizeof(field_elem_u));
        for (slong i = 0; i < poly->npars; i++) {
            field_init_elem(&unified_params[i], unified_ctx->field_id, ctx_ptr);
            fq_nmod_to_field_elem(&unified_params[i], param_vals[i], unified_ctx);
        }
    }

    field_elem_u term_val;
    field_init_elem(&term_val, unified_ctx->field_id, ctx_ptr);
    // Evaluate parameter powers with the existing field backend.
    for (slong i = 0; i < dr_mpoly_length(poly); i++) {
        DR_MPOLY_TERM(term_32, poly, i);

        // Convert coefficient to unified format
        fq_nmod_to_field_elem(&term_val, term_32.coeff, unified_ctx);

        // Multiply by parameter powers
        if (poly->npars > 0 && term_32.par_exp) {
            for (slong j = 0; j < poly->npars; j++) {
                for (slong k = 0; k < term_32.par_exp[j]; k++) {
                    field_mul(&term_val, &term_val, &unified_params[j], unified_ctx->field_id,
                              ctx_ptr);
                }
            }
        }

        // Add to result
        field_add(&unified_result, &unified_result, &term_val, unified_ctx->field_id, ctx_ptr);
    }

    field_clear_elem(&term_val, unified_ctx->field_id, ctx_ptr);
    // Convert result back to fq_nmod
    field_elem_to_fq_nmod(result, &unified_result, unified_ctx);

    // Cleanup
    field_clear_elem(&unified_result, unified_ctx->field_id, ctx_ptr);
    if (unified_params) {
        for (slong i = 0; i < poly->npars; i++) {
            field_clear_elem(&unified_params[i], unified_ctx->field_id, ctx_ptr);
        }
        free(unified_params);
    }
}

void dr_mpoly_make_monic(unified_mpoly_struct *poly)
{
    dr_mpoly_normalize(poly);
    if (dr_mpoly_length(poly) == 0) {
        return; // Empty polynomial, nothing to do
    }

    // Find the leading term (highest total degree in parameters)
    slong leading_idx = 0;
    slong max_deg = -1;

    for (slong i = 0; i < dr_mpoly_length(poly); i++) {
        DR_MPOLY_TERM(term_33, poly, i);

        slong total_deg = 0;

        // Calculate total degree in parameters
        if (term_33.par_exp && poly->npars > 0) {
            for (slong j = 0; j < poly->npars; j++) {
                total_deg += term_33.par_exp[j];
            }
        }

        // Update leading term if this has higher degree
        if (total_deg > max_deg) {
            max_deg = total_deg;
            leading_idx = i;
        }
    }

    // Check if already monic
    DR_MPOLY_TERM(term_34, poly, leading_idx);
    if (fq_nmod_is_one(term_34.coeff, poly->ctx)) {
        return; // Already monic
    }

    // Check if leading coefficient is zero (shouldn't happen for valid polynomials)
    if (fq_nmod_is_zero(term_34.coeff, poly->ctx)) {
        printf("Warning: leading coefficient is zero in dr_mpoly_make_monic\n");
        return;
    }

    // Compute inverse of leading coefficient
    fq_nmod_t inv_leading;
    fq_nmod_init(inv_leading, poly->ctx);
    fq_nmod_inv(inv_leading, term_34.coeff, poly->ctx);

    // Multiply entire polynomial by inverse of leading coefficient
    dr_mpoly_scalar_mul(poly, poly, inv_leading);

    // Cleanup
    fq_nmod_clear(inv_leading, poly->ctx);
}

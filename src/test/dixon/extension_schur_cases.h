/* Exact determinant checks at several extension-field points, including
 * nontrivial permutations, deficient layers, empty cores and bad certificates. */
static void check_extension_schur(void)
{
    const ulong primes[] = {2,2,2,2,2,2,2,3};
    const slong degrees[] = {4,8,16,32,64,128,5,3};
    flint_rand_t state; flint_rand_init(state);
    fmpz_t prime; fmpz_init(prime);
    for (slong field = 0; field < 8; field++) {
        fmpz_set_ui(prime,primes[field]);
        fq_nmod_ctx_t ctx; fq_nmod_ctx_init(ctx,prime,degrees[field],"t");
        fq_nmod_t value, x, original, reduced, factor;
        fq_nmod_init(value,ctx); fq_nmod_init(x,ctx);
        fq_nmod_init(original,ctx); fq_nmod_init(reduced,ctx); fq_nmod_init(factor,ctx);
        for (int trial = 0; trial < 12; trial++) {
            slong n = 6, weights[] = {0,0,0,1,1,2};
            slong rows[] = {2,0,5,1,4,3}, cols[] = {1,4,2,0,5,3};
            slong rd[6], cd[6];
            fq_nmod_poly_mat_t a, core;
            fq_nmod_poly_mat_init(a,n,n,ctx); fq_nmod_poly_mat_init(core,1,1,ctx);
            for (slong i = 0; i < n; i++) for (slong j = 0; j < n; j++) {
                slong bound = trial == 0 ? 0 : 2-weights[i]-weights[j];
                for (slong k = 0; k <= bound; k++) {
                    fq_nmod_randtest(value,state,ctx);
                    fq_nmod_poly_set_coeff(fq_nmod_poly_mat_entry(a,i,j),k,value,ctx);
                }
            }
            /* Singular predicted complement, but a usable pivot elsewhere. */
            if (trial == 1) {
                fq_nmod_poly_zero(fq_nmod_poly_mat_entry(a,5,2),ctx);
                fq_nmod_one(value,ctx);
                fq_nmod_poly_set_coeff(fq_nmod_poly_mat_entry(a,5,0),0,value,ctx);
            }
            for (slong i = 0; i < n; i++) {
                rd[i] = trial == 0 ? 0 : weights[rows[i]];
                cd[i] = trial == 0 ? 0 : weights[cols[i]];
            }
            dixon_mq_step4_profile p = {.size=n, .h=2, .sigma=trial == 0 ? 0 : 2,
                .rows=rows, .cols=cols, .rd=rd, .cd=cd};
            if (trial >= 6) {
                for (slong i = 0; i < n; i++) rd[i] = cd[i] = weights[i];
                p.size = 0; p.adaptive_size = n; p.adaptive_sigma = 2;
                p.adaptive_rd = rd; p.adaptive_cd = cd;
            }
            int ok = dixon_extension_schur_adaptive(core,factor,a,&p,ctx);
            assert(ok);
            assert(core->r < n);
            fq_nmod_mat_t full, small;
            fq_nmod_mat_init(full,n,n,ctx); fq_nmod_mat_init(small,core->r,core->r,ctx);
            for (int point = 0; point < 5; point++) {
                fq_nmod_randtest(x,state,ctx);
                for (slong i = 0; i < n; i++) for (slong j = 0; j < n; j++)
                    fq_nmod_poly_evaluate_fq_nmod(fq_nmod_mat_entry(full,i,j),fq_nmod_poly_mat_entry(a,i,j),x,ctx);
                for (slong i = 0; i < core->r; i++) for (slong j = 0; j < core->r; j++)
                    fq_nmod_poly_evaluate_fq_nmod(fq_nmod_mat_entry(small,i,j),fq_nmod_poly_mat_entry(core,i,j),x,ctx);
                fq_nmod_mat_det(original,full,ctx); fq_nmod_mat_det(reduced,small,ctx);
                fq_nmod_mul(reduced,reduced,factor,ctx);
                assert(fq_nmod_equal(original,reduced,ctx));
            }
            /* Certificate rejection must not change the published factor. */
            fq_nmod_set(value,factor,ctx);
            fq_nmod_one(x,ctx);
            fq_nmod_poly_set_coeff(fq_nmod_poly_mat_entry(a,0,0),3,x,ctx);
            assert(!dixon_extension_schur_adaptive(core,factor,a,&p,ctx));
            assert(fq_nmod_equal(value,factor,ctx));
            fq_nmod_mat_clear(full,ctx); fq_nmod_mat_clear(small,ctx);
            fq_nmod_poly_mat_clear(core,ctx); fq_nmod_poly_mat_clear(a,ctx);
        }
        fq_nmod_clear(factor,ctx); fq_nmod_clear(reduced,ctx); fq_nmod_clear(original,ctx);
        fq_nmod_clear(x,ctx); fq_nmod_clear(value,ctx); fq_nmod_ctx_clear(ctx);
    }
    fmpz_clear(prime); flint_rand_clear(state);
}

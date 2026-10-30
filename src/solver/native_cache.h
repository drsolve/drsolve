/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Private per-solve storage. No global cache: recursive solves, variable
 * permutations and different fields cannot reuse one another's inputs. */
#ifndef DRSOLVE_SOLVER_NATIVE_CACHE_H
#define DRSOLVE_SOLVER_NATIVE_CACHE_H

typedef struct {
    char *expression;
    nmod_mpoly_t poly;
    int usable;
} solver_native_entry;

struct solver_native_cache {
    nmod_mpoly_ctx_t ctx;
    slong nvars, count, capacity;
    char **names;
    solver_native_entry **entries;
};

static void solver_native_cache_clear(struct solver_native_cache *cache)
{
    if(!cache)return;
    for(slong i=0;i<cache->count;i++) {
        solver_native_entry *entry=cache->entries[i];
        nmod_mpoly_clear(entry->poly,cache->ctx);
        free(entry->expression);flint_free(entry);
    }
    for(slong i=0;i<cache->nvars;i++)free(cache->names[i]);
    flint_free(cache->names);flint_free(cache->entries);
    nmod_mpoly_ctx_clear(cache->ctx);flint_free(cache);
}

static nmod_mpoly_struct *solver_native_get(polynomial_solutions_t *sols,
    variable_info_t *vars,slong n,const char *expression)
{
    if(!sols || n<=0 || fq_nmod_ctx_degree(sols->ctx)!=1 ||
       g_field_equation_reduction)return NULL;
    struct solver_native_cache *cache=sols->native_cache;
    int matches=cache && cache->nvars==n;
    for(slong i=0;matches && i<n;i++)matches=strcmp(cache->names[i],vars[i].name)==0;
    if(!matches) {
        solver_native_cache_clear(cache);
        cache=flint_calloc(1,sizeof(*cache));sols->native_cache=cache;
        cache->nvars=n;cache->names=flint_malloc(n*sizeof(*cache->names));
        for(slong i=0;i<n;i++)cache->names[i]=strdup(vars[i].name);
        nmod_mpoly_ctx_init(cache->ctx,n,ORD_LEX,fq_nmod_ctx_prime(sols->ctx));
    }
    for(slong i=0;i<cache->count;i++)
        if(strcmp(cache->entries[i]->expression,expression)==0)
            return cache->entries[i]->usable?cache->entries[i]->poly:NULL;

    if(cache->count==cache->capacity) {
        cache->capacity=cache->capacity?2*cache->capacity:4;
        cache->entries=flint_realloc(cache->entries,cache->capacity*sizeof(*cache->entries));
    }
    solver_native_entry *entry=flint_calloc(1,sizeof(*entry));
    cache->entries[cache->count++]=entry;entry->expression=strdup(expression);
    nmod_mpoly_init(entry->poly,cache->ctx);
    parser_state_t state={0};
    state.nvars=n;state.var_names=flint_malloc(n*sizeof(*state.var_names));
    for(slong i=0;i<n;i++)state.var_names[i]=strdup(vars[i].name);
    state.max_pars=16;state.par_names=malloc(state.max_pars*sizeof(char *));
    state.ctx=sols->ctx;state.generator_name=get_generator_name_for_solver(sols->ctx);
    state.input=expression;state.len=strlen(expression);
    fq_nmod_init(state.current.value,sols->ctx);
    unified_mpoly_struct parsed={0};dr_mpoly_init(&parsed,n,0,sols->ctx);
    next_token(&state);parse_polynomial_expression(&state,&parsed);
    if(parsed.field_id==FIELD_ID_NMOD && parsed.nvars==n && !parsed.npars &&
       !state.npars && state.current.type==TOK_EOF) {
        dr_mpoly_normalize(&parsed);
        nmod_mpoly_set(entry->poly,GET_NMOD_POLY(&parsed),cache->ctx);
        entry->usable=1;
    }
    dr_mpoly_clear(&parsed);
    for(slong i=0;i<n;i++)free(state.var_names[i]);
    for(slong i=0;i<state.npars;i++)free(state.par_names[i]);
    flint_free(state.var_names);free(state.par_names);free(state.generator_name);
    free(state.current.str);fq_nmod_clear(state.current.value,sols->ctx);
    return entry->usable?entry->poly:NULL;
}

static char *solver_native_substitute(polynomial_solutions_t *sols,
    variable_info_t *vars,slong n,const char *expression,const fq_nmod_t value)
{
    nmod_mpoly_struct *poly=solver_native_get(sols,vars,n,expression);
    if(!poly)return NULL;
    struct solver_native_cache *cache=sols->native_cache;
    nmod_mpoly_t reduced;nmod_mpoly_init(reduced,cache->ctx);
    nmod_mpoly_evaluate_one_ui(reduced,poly,0,nmod_poly_get_coeff_ui(value,0),cache->ctx);
    char *native=nmod_mpoly_get_str_pretty(reduced,(const char **)cache->names,cache->ctx);
    char *answer=strdup(native);flint_free(native);
    nmod_mpoly_clear(reduced,cache->ctx);return answer;
}

static int solver_native_evaluate(fq_nmod_t result,polynomial_solutions_t *sols,
    variable_info_t *vars,slong n,const char *expression,fq_nmod_t **values)
{
    nmod_mpoly_struct *poly=solver_native_get(sols,vars,n,expression);
    if(!poly)return 0;
    ulong *point=flint_malloc(n*sizeof(*point));
    for(slong i=0;i<n;i++)point[i]=nmod_poly_get_coeff_ui(values[i][0],0);
    fq_nmod_set_ui(result,nmod_mpoly_evaluate_all_ui(poly,point,sols->native_cache->ctx),sols->ctx);
    flint_free(point);return 1;
}

/* Reorder the cached pair to [elimination variable, parameter] without a
 * polynomial -> string -> polynomial round trip. Respect CLI method choice. */
static int solver_native_pair(unified_mpoly_struct *result,
    polynomial_solutions_t *sols,variable_info_t *vars,
    const char *a,const char *b,const char *eliminate)
{
    slong v=strcmp(vars[0].name,eliminate)==0?0:
            strcmp(vars[1].name,eliminate)==0?1:-1;
    if(v<0)return 0;
    nmod_mpoly_struct *pa=solver_native_get(sols,vars,2,a);
    nmod_mpoly_struct *pb=solver_native_get(sols,vars,2,b);
    if(!pa || !pb)return 0;
    unified_mpoly_struct pair[2]={{0},{0}};
    slong permutation[2]={v,1-v};
    for(int i=0;i<2;i++) {
        dr_mpoly_init(pair+i,1,1,sols->ctx);
        nmod_mpoly_compose_nmod_mpoly_gen(GET_NMOD_POLY(pair+i),i?pb:pa,
            permutation,sols->native_cache->ctx,GET_NMOD_CTX(pair[i].ctx_ptr));
        pair[i].canonical=1;
    }
    char *vn[]={vars[v].name},*pn[]={vars[1-v].name};
    char *generator=get_generator_name_for_solver(sols->ctx);
    if(g_resultant_method==RESULTANT_METHOD_MACAULAY)
        fq_macaulay_resultant_with_names(result,pair,1,1,vn,pn,generator);
    else if(g_resultant_method==RESULTANT_METHOD_DIXON_RECURSIVE)
        fq_dixon_fast_resultant_with_names(result,pair,1,1,vn,pn,generator);
    else fq_dixon_resultant_with_names(result,pair,1,1,vn,pn,generator);
    free(generator);dr_mpoly_clear(pair);dr_mpoly_clear(pair+1);
    return result->ctx_ptr!=NULL;
}

#endif

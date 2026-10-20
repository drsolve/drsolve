/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Group sum c_uv M_u S M_v^T by the side with fewer distinct monomials.
 * Aggregate the other side once; operators and transposes are reused across
 * all subset minors in one row layer. Pure-side groups skip the identity.
 */
typedef struct {
    slong source_degree;
    int left_identity,right_identity;
    nmod_mat_t left,right;
} grouped_term;
typedef struct grouped_operator {
    const entry *source;
    slong a,b,count;
    grouped_term *terms;
    struct grouped_operator *next;
} grouped_operator;

static grouped_operator *grouped_build(grouped_operator **cache,map_node **maps,
                                      const entry *e,slong a,slong b)
{
    for(grouped_operator *p=*cache;p;p=p->next)
        if(p->source==e && p->a==a && p->b==b)return p;
    grouped_operator *op=flint_calloc(1,sizeof(*op));
    op->source=e;op->a=a;op->b=b;op->next=*cache;*cache=op;
    op->terms=flint_calloc(FLINT_MAX(1,e->count),sizeof(grouped_term));
    uint64_t *xkeys=flint_malloc(FLINT_MAX(1,e->count)*sizeof(uint64_t));
    uint64_t *ykeys=flint_malloc(FLINT_MAX(1,e->count)*sizeof(uint64_t));
    for(slong ax=0;ax<=a;ax++) {
        slong by=-1,nx=0,ny=0;
        for(slong t=0;t<e->count;t++) {
            pair_term p=e->terms[t];
            if(key_degree(p.x)!=ax || key_degree(p.y)>b)continue;
            by=key_degree(p.y);
            slong i=0;while(i<nx && xkeys[i]!=p.x)i++;if(i==nx)xkeys[nx++]=p.x;
            i=0;while(i<ny && ykeys[i]!=p.y)i++;if(i==ny)ykeys[ny++]=p.y;
        }
        if(by<0)continue;
        slong ca=a-ax,cb=b-by;
        int group_x=nx<ny; /* ties prefer grouping by right monomial */
        slong count=group_x?nx:ny;
        for(slong i=0;i<count;i++) {
            uint64_t fixed=group_x?xkeys[i]:ykeys[i];
            grouped_term *g=op->terms+op->count++;g->source_degree=ca;
            g->left_identity=group_x && fixed==0;
            g->right_identity=!group_x && fixed==0;
            if(!g->left_identity)nmod_mat_init(g->left,spaces[a].h,spaces[ca].h,prime);
            if(!g->right_identity)nmod_mat_init(g->right,spaces[cb].h,spaces[b].h,prime);
            if(group_x && !g->left_identity)
                nmod_mat_set(g->left,multiplication(maps,ca,fixed));
            if(!group_x && !g->right_identity)
                nmod_mat_transpose(g->right,multiplication(maps,cb,fixed));
            /* Assemble the sparse polynomial action before reducing once.
             * C * normal_forms replaces a sum of separately reduced maps. */
            slong from=group_x?cb:ca,to=group_x?b:a;
            slong term_count=0;
            for(slong t=0;t<e->count;t++) {
                pair_term p=e->terms[t];
                if(key_degree(p.x)==ax && key_degree(p.y)==by &&
                   (group_x?p.x:p.y)==fixed)term_count++;
            }
            /* Short polynomials are cheaper as a linear combination of
             * cached maps; dense action/RREF multiplication suits long ones. */
            if(term_count<=16) {
                for(slong t=0;t<e->count;t++) {
                    pair_term p=e->terms[t];
                    if(key_degree(p.x)!=ax || key_degree(p.y)!=by ||
                       (group_x?p.x:p.y)!=fixed)continue;
                    const nmod_mat_struct *m=multiplication(maps,from,group_x?p.y:p.x);
                    if(group_x) {
                        for(slong r=0;r<m->r;r++)for(slong c=0;c<m->c;c++)
                            nmod_mat_entry(g->right,c,r)=nmod_add(nmod_mat_entry(g->right,c,r),
                                nmod_mul(p.c,nmod_mat_entry(m,r,c),mod),mod);
                    } else {
                        for(slong r=0;r<m->r;r++)for(slong c=0;c<m->c;c++)
                            nmod_mat_entry(g->left,r,c)=nmod_add(nmod_mat_entry(g->left,r,c),
                                nmod_mul(p.c,nmod_mat_entry(m,r,c),mod),mod);
                    }
                }
                continue;
            }
            nmod_mat_t coefficients,action;
            nmod_mat_init(coefficients,spaces[from].h,spaces[to].size,prime);
            nmod_mat_init(action,spaces[from].h,spaces[to].h,prime);
            for(slong t=0;t<e->count;t++) {
                pair_term p=e->terms[t];
                if(key_degree(p.x)!=ax || key_degree(p.y)!=by ||
                   (group_x?p.x:p.y)!=fixed)continue;
                const map_node *m=monomial_rows(maps,from,group_x?p.y:p.x);
                for(slong r=0;r<spaces[from].h;r++) {
                    slong c=m->rows[r];
                    nmod_mat_entry(coefficients,r,c)=nmod_add(nmod_mat_entry(coefficients,r,c),p.c,mod);
                }
            }
            nmod_mat_mul(action,coefficients,spaces[to].reduce);
            if(group_x)nmod_mat_set(g->right,action);else nmod_mat_transpose(g->left,action);
            nmod_mat_clear(coefficients);nmod_mat_clear(action);
        }
    }
    flint_free(xkeys);flint_free(ykeys);return op;
}
static void grouped_clear(grouped_operator *cache)
{
    while(cache) {
        grouped_operator *next=cache->next;
        for(slong i=0;i<cache->count;i++) {
            grouped_term *g=cache->terms+i;
            if(!g->left_identity)nmod_mat_clear(g->left);
            if(!g->right_identity)nmod_mat_clear(g->right);
        }
        flint_free(cache->terms);flint_free(cache);cache=next;
    }
}
/* Workspace matrices are reused throughout a layer for each needed shape. */
typedef struct grouped_workspace {
    slong rows,cols;
    nmod_mat_t matrix;
    struct grouped_workspace *next;
} grouped_workspace;
static nmod_mat_struct *grouped_scratch(grouped_workspace **cache,slong rows,slong cols)
{
    for(grouped_workspace *p=*cache;p;p=p->next)
        if(p->rows==rows && p->cols==cols)return p->matrix;
    grouped_workspace *p=flint_malloc(sizeof(*p));p->rows=rows;p->cols=cols;
    nmod_mat_init(p->matrix,rows,cols,prime);p->next=*cache;*cache=p;return p->matrix;
}
static void grouped_scratch_clear(grouped_workspace *cache)
{
    while(cache){grouped_workspace *next=cache->next;nmod_mat_clear(cache->matrix);flint_free(cache);cache=next;}
}
static void grouped_apply(nmod_mat_t dst,const grouped_operator *op,
                          const nmod_mat_struct *child,int negative,
                          grouped_workspace **temps,grouped_workspace **products)
{
    nmod_mat_struct *product=grouped_scratch(products,dst->r,dst->c);
    for(slong i=0;i<op->count;i++) {
        const grouped_term *g=op->terms+i;
        const nmod_mat_struct *src=child+g->source_degree;
        if(!src->r)continue;
        matrix_products+=(g->left_identity || g->right_identity)?1:2;
        if(g->left_identity)nmod_mat_mul(product,src,g->right);
        else if(g->right_identity)nmod_mat_mul(product,g->left,src);
        else {
            /* Choose the cheaper parenthesization for rectangular factors. */
            long double left_cost=(long double)g->left->r*src->r*src->c+
                                  (long double)g->left->r*src->c*g->right->c;
            long double right_cost=(long double)src->r*src->c*g->right->c+
                                   (long double)g->left->r*src->r*g->right->c;
            if(left_cost<=right_cost) {
                nmod_mat_struct *tmp=grouped_scratch(temps,g->left->r,src->c);
                nmod_mat_mul(tmp,g->left,src);nmod_mat_mul(product,tmp,g->right);
            } else {
                nmod_mat_struct *tmp=grouped_scratch(temps,src->r,g->right->c);
                nmod_mat_mul(tmp,src,g->right);nmod_mat_mul(product,g->left,tmp);
            }
        }
        if(negative)nmod_mat_sub(dst,dst,product);else nmod_mat_add(dst,dst,product);
    }
}

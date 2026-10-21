/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Monomials are ranks in the enumerated graded supports, not bit-packed
 * exponents. The immutable exponent table is shared by all sample workers. */
typedef slong dq_key;
static slong *key_offsets,*key_exponents,*key_degrees;
static slong exponent(dq_key key,slong v) {return key_exponents[key*n+v];}
static slong key_degree(dq_key key) {return key_degrees[key];}
static dq_key key_rank(const slong *e)
{
    slong d=0;for(slong v=0;v<n;v++)d+=e[v];
    assert(d<=limit);
    slong z=e[n-1],remaining=d-z;
    slong rank=binomial(n-1+d,d)-binomial(n-1+d-z,d-z);
    for(slong v=0;v<n-2;v++) {
        slong smaller=remaining-e[v];
        if(smaller)rank+=binomial(n-2-v+smaller-1,smaller-1);
        remaining-=e[v];
    }
    return key_offsets[d]+rank;
}
static dq_key key_set(dq_key key,slong v,slong value)
{
    slong e[n];memcpy(e,key_exponents+key*n,n*sizeof(slong));e[v]=value;
    return key_rank(e);
}
static dq_key key_add(dq_key a,dq_key b)
{
    slong e[n];for(slong v=0;v<n;v++)e[v]=exponent(a,v)+exponent(b,v);
    return key_rank(e);
}
static void key_tables_init(void)
{
    key_offsets=flint_calloc(limit+2,sizeof(slong));
    for(slong d=0;d<=limit;d++)key_offsets[d+1]=key_offsets[d]+binomial(n+d-1,d);
    slong count=key_offsets[limit+1];
    key_exponents=flint_calloc(count*n,sizeof(slong));
    key_degrees=flint_malloc(count*sizeof(slong));
}
static void key_tables_clear(void)
{
    flint_free(key_offsets);flint_free(key_exponents);flint_free(key_degrees);
    key_offsets=key_exponents=key_degrees=NULL;
}

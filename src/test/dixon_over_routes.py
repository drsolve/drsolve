"""Exact Sage experiments for overdetermined Dixon/quotient routes.

This is an experimental harness, not a production solver. Input generation and
independent ideal audits are outside timings; preparation, relation generation,
projection, and the algebra/point certificate are inside timings. No solver
receives planted coordinates. See dixon_over_routes_bench.py for paired inputs.
"""
from sage.all import (GF, PolynomialRing, TermOrder, matrix, vector, identity_matrix, prod)
from dixon_closure_experiment import Closure, exponents, original_multiples, partial_dixon
import itertools
import time
import copy


def monomial(P, e):
    return P({tuple(e): P.base_ring().one()})


def echelon(M):
    E = M.echelon_form()
    return E[:E.rank(), :]


def quotient_projection(E):
    """Rows v map to v*W; ker(W) is exactly rowspace(E)."""
    F = E.base_ring(); nc = E.ncols(); piv = list(E.pivots())
    ps = set(piv); free = [j for j in range(nc) if j not in ps]
    W = matrix(F, nc, len(free))
    for k, j in enumerate(free):
        W[j, k] = 1
    if piv and free:
        neg = -E.matrix_from_columns(free)
        for i, j in enumerate(piv):
            W[j, :] = neg[i, :]
    return W, free


class Space:
    """Incremental exact relation space with growing monomial support."""
    def __init__(self, P):
        self.P = P; self.F = P.base_ring(); self.n = P.ngens()
        self.ee = []; self.ix = {}; self.E = matrix(self.F, 0, 0)
        self.rows = 0; self.cells = 0; self.matmul = 0; self.history = []
        self.original_degree = 2

    def extend(self, support):
        new = sorted(set(self.ee).union(support), key=lambda e: (sum(e), e), reverse=True)
        if new == self.ee:
            return
        ix = {e: j for j, e in enumerate(new)}
        E = matrix(self.F, self.E.nrows(), len(new))
        for j, e in enumerate(self.ee):
            E[:, ix[e]] = self.E[:, j]
        self.ee = new; self.ix = ix; self.E = E

    def coefficients(self, polys):
        M = matrix(self.F, len(polys), len(self.ee))
        for i, f in enumerate(polys):
            for e, c in f.dict().items():
                M[i, self.ix[tuple(e)]] = c
        return M

    def insert(self, polys):
        polys = [f for f in polys if f]
        self.extend(tuple(e) for f in polys for e in f.dict())
        return self.insert_matrix(self.coefficients(polys))

    def insert_matrix(self, M):
        self.rows += M.nrows()
        if not M.nrows():
            return 0
        E = self.E; old = E.nrows(); nc = E.ncols()
        if not old:
            self.cells += M.nrows()*nc; self.E = echelon(M)
            return self.E.nrows()
        piv = list(E.pivots()); ps = set(piv); free = [j for j in range(nc) if j not in ps]
        R = M.matrix_from_columns(free)-M.matrix_from_columns(piv)*E.matrix_from_columns(free)
        self.matmul += M.nrows()*old*len(free)
        self.cells += R.nrows()*R.ncols(); S = echelon(R)
        if not S.nrows():
            return 0
        newp = [free[j] for j in S.pivots()]
        L = matrix(self.F, S.nrows(), nc)
        for j, k in enumerate(free):
            L[:, k] = S[:, j]
        E = E-E.matrix_from_columns(newp)*L
        self.matmul += old*S.nrows()*nc
        joined = E.stack(L); pp = piv+newp
        self.E = joined.matrix_from_rows(sorted(range(len(pp)), key=pp.__getitem__))
        return S.nrows()

    def stats(self):
        return dict(columns=len(self.ee), submitted_rows=int(self.rows),
                    rref_cells=int(self.cells), dense_matmul_multiply_adds=int(self.matmul))


def row_action(f, row, operators):
    result = vector(row.base_ring(), len(row))
    for e, c in f.dict().items():
        v = row
        for i, a in enumerate(e):
            for _ in range(a):
                v = v*operators[i]
        result += c*v
    return result


def algebra_certificate(fs, basis, operators, one):
    """Relations used to obtain this representation must already lie in I.

    Commutation + cyclic basis + input vanishing supplies the reverse ideal
    inclusion. Works with multiple roots, nilpotents, and projection collisions.
    """
    F = fs[0].base_ring(); dim = len(one)
    for i, A in enumerate(operators):
        for B in operators[:i]:
            if A*B != B*A:
                return None
    for j, b in enumerate(basis):
        e = vector(F, dim); e[j] = 1
        if row_action(b, one, operators) != e:
            return None
    if any(row_action(f, one, operators) for f in fs):
        return None
    result = dict(certified=True, dimension=dim, recovered=None)
    if dim == 1:
        a = [T[0, 0] for T in operators]
        assert all(f(*a) == 0 for f in fs)
        result['recovered'] = list(map(int, a))
    # Exact elimination polynomial, including multiplicities; no separation
    # assumption for the last coordinate. Matrix orientation is immaterial.
    result['eliminant'] = list(map(int, operators[-1].minimal_polynomial().list()))
    return result


def small_quotient(space, fs):
    """Close in residual coordinates; return a certified quotient if possible.

    All generated defects are ideal relations: x_i(m-N(m)) and reductions of
    x_i*b. Commutators and input actions have the same ideal-membership proof.
    Ambient relations are not rebuilt after each small-space compression.
    """
    P = space.P; F = space.F; ee = space.ee; ix = space.ix; n = space.n
    W, free = quotient_projection(space.E); dims = [W.ncols()]
    if not W.ncols():
        return dict(certified=True, dimension=0, recovered=None, eliminant=[1], dimensions=dims)
    zero = (0,)*n
    if zero in ix and not W[ix[zero]]:
        return dict(certified=True, dimension=0, recovered=None, eliminant=[1], dimensions=[0])
    low = [j for j, e in enumerate(ee)
           if all(tuple(a+(k == i) for k, a in enumerate(e)) in ix for i in range(n))]
    low.sort(key=lambda j: (sum(ee[j]), ee[j]))
    if not low:
        return dict(certified=False, reason='missing border', dimensions=dims)
    shifts = [[ix[tuple(a+(k == i) for k, a in enumerate(ee[j]))] for j in low] for i in range(n)]
    for iteration in range(len(free)+2):
        dim = W.ncols()
        local = list(W.matrix_from_rows(low).transpose().pivots())
        if len(local) != dim:
            return dict(certified=False, reason='low basis incomplete', dimensions=dims)
        basis_idx = [low[j] for j in local]
        W = W*W.matrix_from_rows(basis_idx).inverse()
        basis = [monomial(P, ee[j]) for j in basis_idx]
        operators = [W.matrix_from_rows([shift[j] for j in local]) for shift in shifts]
        WL = W.matrix_from_rows(low)
        C = matrix(F, 0, dim)
        for shift, T in zip(shifts, operators):
            C = C.stack(W.matrix_from_rows(shift)-WL*T)
            space.matmul += len(low)*dim*dim
        space.cells += C.nrows()*dim
        C = echelon(C)
        if not C.nrows():
            one = W[ix[zero]]
            cert = algebra_certificate(fs, basis, operators, one)
            if cert:
                cert['dimensions'] = dims
                return cert
            # A noncommuting candidate is not a stopping certificate. Its
            # defects still give valid ideal relations, and must be closed.
            for i, A in enumerate(operators):
                for B in operators[:i]:
                    C = C.stack(A*B-B*A)
            C = C.stack(matrix(F, [row_action(f, one, operators) for f in fs]))
            space.cells += C.nrows()*dim; C = echelon(C)
            if not C.nrows():
                return dict(certified=False, reason='cyclic certificate failed', dimensions=dims)
        projection, keep = quotient_projection(C)
        space.matmul += W.nrows()*dim*len(keep)
        W = W*projection; dims.append(W.ncols())
        if not W.ncols() or not W[ix[zero]]:
            return dict(certified=True, dimension=0, recovered=None, eliminant=[1], dimensions=dims)
    raise AssertionError('quotient compression did not decrease dimension')


def finish_space(space, fs):
    result=small_quotient(space,fs)
    if result['certified']:return result
    # Degree falls can prevent the current low-degree monomials from spanning
    # the cokernel even for a finite affine ideal. Saturate the bounded space
    # before increasing the degree; use the same audited closure as baseline.
    d=max(map(sum,space.ee));C=Closure(space.P,d)
    space.extend(C.ee);C.E=space.E
    C.processed=C.rref(C.coefficients(original_multiples(fs,space.original_degree-1)))
    oldrank=C.E.nrows();C.run(fs,[])
    space.E=C.E;space.cells+=C.rref_cells;space.matmul+=C.matmul_work;space.rows+=C.submitted_rows
    if C.E.nrows()>oldrank:
        result=small_quotient(space,fs);result['bounded_prolongation']=True
    return result


class SharedDixon:
    """Cache cofactor DP for n-1 fixed equations; stream complete coefficients.

    DP subsets index rows, columns are fixed input equations. All n cofactors
    are obtained from the same DP. Auxiliary-degree pruning is exact and only
    uses multiplication-compatible monomial truncation/reachability bounds.
    """
    def __init__(self, core, degree, min_degree=3, rules=()):
        start = time.perf_counter()
        self.P = P = core[0].parent(); self.n = n = P.ngens(); k = n-1
        assert len(core) == k
        self.lo = max(0, n+1-degree); self.hi = n+1-min_degree
        self.Q = Q = PolynomialRing(P.base_ring(), names=P.variable_names()+tuple('y%d'%i for i in range(k)),
                                    order=TermOrder('deglex',k)+TermOrder('deglex',n))
        self.z = z = Q.gens(); self.embed = P.hom(z[:n], Q)
        self.rules = [self.embed(f) for f in rules]
        self.reductions = 0
        columns = [self.column(f) for f in core]
        dp = {0: Q.one()}; self.products = 0; self.peak = 0
        for mask in range(1, (1 << n)-1):
            rows = [r for r in range(n) if mask >> r & 1]; col = len(rows)-1
            g = Q.zero()
            for j, r in enumerate(rows):
                a, b = columns[col][r], dp[mask ^ (1 << r)]
                self.products += len(a.dict())*len(b.dict())
                g += (-1)**(j+col)*a*b
            if self.rules and g:
                g=g.reduce(self.rules);self.reductions+=1
            remaining = n-len(rows)
            g = Q({e:c for e,c in g.dict().items()
                   if self.lo-remaining <= sum(e[n:]) <= self.hi})
            dp[mask] = g; self.peak = max(self.peak, len(g.dict()))
        self.cofactors = []
        for r in range(n):
            self.cofactors.append(self.group((-1)**(r+n-1)*dp[((1 << n)-1) ^ (1 << r)]))
        self.seconds = time.perf_counter()-start
        self.generated = 0; self.output_terms = 0

    def column(self, f):
        original = self.embed(f); previous = original; column = []
        for i in range(self.n-1):
            current = original.subs({self.z[j]:self.z[self.n+j] for j in range(i+1)})
            column.append((previous-current)//(self.z[i]-self.z[self.n+i]))
            previous = current
        return column+[original]

    def group(self, f):
        groups = {}
        for e, c in f.dict().items():
            groups.setdefault(tuple(e[self.n:]), {})[tuple(e[:self.n])] = c
        return {b:self.P(d) for b,d in groups.items()}

    def prepare(self, f):
        entries = [self.group(a) for a in self.column(f)]
        possible = set()
        for A, B in zip(entries, self.cofactors):
            for a in A:
                for b in B:
                    e = tuple(x+y for x,y in zip(a,b))
                    if self.lo <= sum(e) <= self.hi:
                        possible.add(e)
        return entries, sorted(possible, key=lambda b: (-sum(b), b))

    def coefficient(self, prepared, beta):
        entries, _ = prepared; g = self.P.zero()
        for A, B in zip(entries, self.cofactors):
            for a, f in A.items():
                b = tuple(x-y for x,y in zip(beta,a))
                h = B.get(b)
                if h is not None:
                    self.products += len(f.dict())*len(h.dict()); g += f*h
        self.generated += 1; self.output_terms += len(g.dict())
        return g


def solve_quotient(fs, max_degree=6, dixon_degree=0, batch=8, support_only=False):
    start = time.perf_counter(); P = fs[0].parent(); n = P.ngens()
    S = Space(P); stages = []; cache = None; generated = 0; dixon_seconds = 0.
    result = dict(certified=False)
    for d in range(2, max_degree+1):
        # Only new original multiples; prior echelon information is reused.
        rows = [monomial(P,e)*f for f in fs if f
                for e in exponents(n, d-int(f.total_degree()))
                if sum(e)+f.total_degree() == d or d == 2]
        S.extend(exponents(n,d)); S.insert(rows)
        S.original_degree=d
        result = finish_space(S,fs)
        stages.append(dict(degree=d, source='original', rank=int(S.E.nrows()), **result))
        if result['certified']:
            break
        if dixon_degree and d == 3:
            dt = time.perf_counter(); cache = SharedDixon(fs[:n-1],dixon_degree,min_degree=d+1)
            prepared = [cache.prepare(f) for f in fs[n-1:]]
            dixon_seconds += time.perf_counter()-dt
            if not support_only:
                S.extend(exponents(n,dixon_degree))
            # Interleave equations at each coefficient. Neither planted roots
            # nor a GB nor an offline rank oracle is used to order this stream.
            betas = sorted(set(b for _,bb in prepared for b in bb), key=lambda b:(-sum(b),b))
            pending = []
            for beta in betas:
                for pre in prepared:
                    dt = time.perf_counter(); g = cache.coefficient(pre,beta)
                    dixon_seconds += time.perf_counter()-dt; generated += 1
                    if g:
                        pending.append(g)
                    if len(pending) >= batch:
                        S.insert(pending); pending = []; result = finish_space(S,fs)
                        stages.append(dict(degree=dixon_degree, source='dixon', generated=generated,
                                           rank=int(S.E.nrows()), **result))
                        if result['certified']:
                            break
                if result['certified']:
                    break
            if pending and not result['certified']:
                S.insert(pending); result = finish_space(S,fs)
                stages.append(dict(degree=dixon_degree,source='dixon',generated=generated,**result))
            if result['certified']:
                break
    return dict(seconds=time.perf_counter()-start, **result, **S.stats(), stages=stages,
                dixon_seconds=dixon_seconds, dixon_coefficients=generated,
                dixon_term_products=cache.products if cache else 0,
                support_only=support_only), S


def solve_scalar(fs, degree, extra=()):
    start = time.perf_counter(); C = Closure(fs[0].parent(),degree)
    result = C.run(fs,list(fs)+list(extra)); result['seconds'] = time.perf_counter()-start
    result['certified'] = result['recovered'] is not None or result['nullity'] == 0
    result['dimension'] = 1 if result['recovered'] is not None else (0 if result['nullity'] == 0 else None)
    return result


def scalar_inject(fs,degree,extra,after=0):
    """Same bounded closure, optionally inject after several completed rounds."""
    start=time.perf_counter();C=Closure(fs[0].parent(),degree);C.insert(C.coefficients(fs))
    for _ in range(after):
        if C.E.nrows()>=C.nc-1:break
        low=C.E.matrix_from_rows([i for i,p in enumerate(C.E.pivots()) if sum(C.ee[p])<degree])
        delta=C.rref(C.residual(low,C.processed));C.processed=low
        if not delta.nrows():break
        rows=matrix(C.F,[[row[j] if j>=0 else C.F.zero() for j in perm] for row in delta.rows() for perm in C.inverse])
        rank=C.E.nrows();C.insert(rows)
        if C.E.nrows()==rank:break
    C.insert(C.coefficients(extra));r=C.run(fs,[])
    r.update(seconds=time.perf_counter()-start,inject_after=after,
             certified=r['recovered'] is not None or r['nullity']==0)
    return r


def solve_gb(fs):
    start = time.perf_counter(); I = fs[0].parent().ideal(fs)
    G = I.groebner_basis()
    if any(g == 1 for g in G):
        return dict(seconds=time.perf_counter()-start, certified=True, dimension=0, eliminant=[1], recovered=None)
    B = list(I.normal_basis()); T, one = multiplication_model(I,B,fs[0].parent())
    cert = algebra_certificate(fs,B,T,one)
    assert cert
    return dict(seconds=time.perf_counter()-start, **cert)


def multiplication_model(I, B, P):
    F = P.base_ring(); idx = {tuple(b.exponents()[0]):i for i,b in enumerate(B)}
    def nf(f):
        v = vector(F,len(B))
        for e,c in I.reduce(f).dict().items():
            v[idx[tuple(e)]] = c
        return v
    T = [matrix(F,[nf(x*b) for b in B]) for x in P.gens()]
    return T, nf(P.one())


def solve_core_quotient(fs, stream=False):
    """Villard-style quotient interface, with core GB/NF setup fully charged."""
    start = time.perf_counter(); P = fs[0].parent(); n = P.ngens(); F = P.base_ring()
    I = P.ideal(fs[:n]); I.groebner_basis(); B = sorted(I.normal_basis(),key=lambda b:(b.total_degree(),str(b)))
    idx = {tuple(b.exponents()[0]):i for i,b in enumerate(B)}; D = len(B)
    setup = time.perf_counter()-start; E = matrix(F,0,D); processed = 0; recovered = None
    def nf(f):
        v = vector(F,D)
        for e,c in I.reduce(f).dict().items():
            v[idx[tuple(e)]] = c
        return v
    one = nf(P.one()); batch = []; relation_start = time.perf_counter()
    for b in B:
        for g in fs[n:]:
            batch.append(nf(b*g)); processed += 1
            if stream and len(batch) >= 8:
                E = echelon(E.stack(matrix(F,batch))); batch = []
                if E.nrows() == D:
                    cert = dict(certified=True,dimension=0,recovered=None,eliminant=[1])
                    return dict(seconds=time.perf_counter()-start,core_setup_seconds=setup,
                                core_dimension=D,relations=processed,stream=stream,**cert)
                if E.nrows() == D-1:
                    W,_ = quotient_projection(E); a = one*W
                    if a[0]:
                        point = [(nf(x)*W)[0]/a[0] for x in P.gens()]
                        # Verifying the entire basis evaluation character makes
                        # W's kernel the exact ideal of this point in A_core.
                        if all(f(*point)==0 for f in fs) and all((nf(b)*W)[0]==a[0]*b(*point) for b in B):
                            recovered = list(map(int,point)); break
        if recovered is not None:
            break
    if recovered is not None:
        cert = dict(certified=True,dimension=1,recovered=recovered,
                    eliminant=[int(-F(recovered[-1])),1])
    else:
        if batch:
            E = echelon(E.stack(matrix(F,batch)))
        W,free = quotient_projection(E)
        if not free:
            cert = dict(certified=True,dimension=0,recovered=None,eliminant=[1])
        else:
            # NF(x*b) only for surviving basis elements; no full coordinate
            # matrices are necessary after the extra-ideal span is complete.
            Bsmall = [B[j] for j in free]
            T = [matrix(F,[nf(x*b)*W for b in Bsmall]) for x in P.gens()]
            cert = algebra_certificate(fs,Bsmall,T,one*W)
            assert cert
    return dict(seconds=time.perf_counter()-start, core_setup_seconds=setup,
                relations_seconds=time.perf_counter()-relation_start,core_dimension=D,
                relations=processed,full_relations=D*(len(fs)-n),stream=stream,**cert)


def module_matrix(polys):
    P = polys[0].parent(); F = P.base_ring(); n = P.ngens()
    R = PolynomialRing(F,'t')
    cols = sorted({tuple(e[:-1]) for f in polys for e in f.dict()}.union({(0,)*(n-1)}),
                  key=lambda e:(sum(e),e), reverse=True)
    ix = {e:j for j,e in enumerate(cols)}; M = matrix(R,len(polys),len(cols))
    for i,f in enumerate(polys):
        entries = {}
        for e,c in f.dict().items():
            e=tuple(e)
            entries.setdefault(ix[tuple(e[:-1])],{})[e[-1]] = c
        for j,d in entries.items():
            M[i,j] = R(d)
    return M,cols


def module_reduce(M, use_units):
    """Unimodular constant-pivot elimination followed by residual row HNF."""
    R = M.base_ring(); A = copy.copy(M); piv = []; rank = 0; products = 0
    if use_units:
        for col in range(A.ncols()):
            candidate = next((i for i in range(rank,A.nrows()) if A[i,col] and A[i,col].degree()==0),None)
            if candidate is None:
                continue
            A.swap_rows(rank,candidate); A.rescale_row(rank,~A[rank,col][0])
            for i in range(A.nrows()):
                if i != rank and A[i,col]:
                    products += A.ncols(); A.add_multiple_of_row(i,rank,-A[i,col])
            piv.append(col); rank += 1
    free = [j for j in range(A.ncols()) if j not in set(piv)]
    residual = A.matrix_from_rows_and_columns(range(rank,A.nrows()),free)
    H = residual.hermite_form(include_zero_rows=False) if free else matrix(R,0,0)
    rules = A[:rank,:]
    return H,rules,piv,free,dict(unit_pivots=rank,residual_rows=residual.nrows(),
                              residual_columns=len(free),unit_poly_products=products)


def module_certificate(fs, cols, H, rules, piv, free):
    P = fs[0].parent(); F = P.base_ring(); n = P.ngens(); R = H.base_ring()
    if H.nrows() != len(free):
        return dict(certified=False, reason='module has free generators')
    if any(not H[j,j] for j in range(len(free))):
        return dict(certified=False,reason='module not triangular')
    degrees = [int(H[j,j].degree()) for j in range(len(free))]
    slots = [(j,k) for j,d in enumerate(degrees) for k in range(d)]
    B = [monomial(P,cols[free[j]]+(k,)) for j,k in slots]
    ix = {e:j for j,e in enumerate(cols)}
    def nf(f):
        v = vector(R,len(cols))
        for e,c in f.dict().items():
            e=tuple(e)
            if tuple(e[:-1]) not in ix:
                return None
            v[ix[tuple(e[:-1])]] += c*R.gen()**e[-1]
        for i,j in enumerate(piv):
            v -= v[j]*rules[i]
        w = vector(R,[v[j] for j in free])
        for j in range(len(free)):
            q,_ = w[j].quo_rem(H[j,j]); w -= q*H[j]
        return vector(F,[w[j][k] for j,k in slots])
    one = nf(P.one())
    if not slots:
        return dict(certified=True,dimension=0,recovered=None,eliminant=[1])
    rows = [[nf(x*b) for b in B] for x in P.gens()]
    if any(v is None for rr in rows for v in rr):
        return dict(certified=False,reason='module missing border',dimension=len(B))
    T = [matrix(F,rr) for rr in rows]
    cert = algebra_certificate(fs,B,T,one)
    return cert or dict(certified=False,reason='module compatibility incomplete',dimension=len(B))


def solve_module(fs, dixon_degree=0, units=False, stream=False, multiplier_degree=1):
    start = time.perf_counter(); P = fs[0].parent(); n = P.ngens()
    polys = [monomial(P,e+(0,))*f for f in fs for e in exponents(n-1,multiplier_degree)]
    construction = 0.; stages = []; cache = None; generated = 0
    if dixon_degree:
        dt = time.perf_counter(); cache = SharedDixon(fs[:n-1],dixon_degree)
        prepared = [cache.prepare(f) for f in fs[n-1:]]
        construction += time.perf_counter()-dt
        betas = sorted(set(b for _,bb in prepared for b in bb),key=lambda b:(-sum(b),b))
        batches = ((pre,b) for b in betas for pre in prepared)
    else:
        batches = iter(())
    def finish():
        M,cols = module_matrix(polys); H,U,piv,free,stats = module_reduce(M,units)
        cert = module_certificate(fs,cols,H,U,piv,free)
        stages.append(dict(rows=M.nrows(),columns=M.ncols(),**stats,**cert))
        return cert
    cert = finish() if stream or not dixon_degree else dict(certified=False)
    if not cert['certified']:
        for pre,b in batches:
            dt=time.perf_counter();g=cache.coefficient(pre,b);construction+=time.perf_counter()-dt
            generated += 1
            if g: polys.append(g)
            if stream and generated % 16 == 0:
                cert=finish()
                if cert['certified']:break
        if not cert['certified']:
            cert=finish()
    return dict(seconds=time.perf_counter()-start,dixon_seconds=construction,
                dixon_coefficients=generated,stages=stages,**cert)


def solve_module_nf(fs, degree=4, multiplier_degree=1):
    """Reduce the cofactor DP by a proven subideal, not a linear projection.

    Each rewrite rule has a unit leading x coefficient and belongs to I. Every
    polynomial division changes its argument by an element of their ideal J.
    Thus DP arithmetic is congruent modulo J, whether or not these rules form a
    GB. Output coefficients are still in I; they need not equal raw Dixon
    coefficients. No claim is made that an arbitrary residual map is a ring map.
    """
    start=time.perf_counter();P=fs[0].parent();n=P.ngens()
    polys=[monomial(P,e+(0,))*f for f in fs for e in exponents(n-1,multiplier_degree)]
    M,cols=module_matrix(polys);H,U,piv,free,stats=module_reduce(M,True)
    cert=module_certificate(fs,cols,H,U,piv,free);stages=[dict(**stats,**cert)]
    if cert['certified']:
        return dict(seconds=time.perf_counter()-start,dixon_coefficients=0,dixon_seconds=0.,stages=stages,**cert)
    rules=[]
    for row in U.rows():
        lead=next((j for j,c in enumerate(row) if c),None)
        if lead is not None and row[lead].degree()==0:
            rules.append(sum(P({e+(k,):c}) for p,e in zip(row,cols) for k,c in enumerate(p.list())))
    dt=time.perf_counter();cache=SharedDixon(fs[:n-1],degree,rules=rules)
    construction=time.perf_counter()-dt;generated=0;audit_polys=[]
    for g in fs[n-1:]:
        dt=time.perf_counter();pre=cache.prepare(g);construction+=time.perf_counter()-dt
        for beta in pre[1]:
            dt=time.perf_counter();f=cache.coefficient(pre,beta)
            # Same subideal reduction after the final determinant expansion.
            if rules and f:f=P(cache.embed(f).reduce(cache.rules))
            construction+=time.perf_counter()-dt;generated+=1
            if f:polys.append(f);audit_polys.append(f)
            if generated%8==0:
                M,cols=module_matrix(polys);H,U,piv,free,stats=module_reduce(M,True)
                cert=module_certificate(fs,cols,H,U,piv,free);stages.append(dict(**stats,**cert))
                if cert['certified']:break
        if cert['certified']:break
    if not cert['certified']:
        M,cols=module_matrix(polys);H,U,piv,free,stats=module_reduce(M,True)
        cert=module_certificate(fs,cols,H,U,piv,free);stages.append(dict(**stats,**cert))
    return dict(seconds=time.perf_counter()-start,dixon_seconds=construction,
                dixon_coefficients=generated,nf_rules=len(rules),nf_reductions=cache.reductions,
                dixon_term_products=cache.products,stages=stages,
                _audit_relations=audit_polys,**cert)


def shared_construction(fs, degree=4):
    n=fs[0].parent().ngens(); start=time.perf_counter(); independent=[]; products=0
    for f in fs[n-1:]:
        gs,stats=partial_dixon(fs[:n-1]+[f],degree)
        independent.append(gs); products+=stats['candidate_term_products']
    repeated_seconds=time.perf_counter()-start
    start=time.perf_counter();C=SharedDixon(fs[:n-1],degree);shared=[]
    for f in fs[n-1:]:
        pre=C.prepare(f);gg={b:C.coefficient(pre,b) for b in pre[1]}
        shared.append({b:g for b,g in gg.items() if g})
    shared_seconds=time.perf_counter()-start
    assert independent==shared
    return dict(seconds=shared_seconds,repeated_seconds=repeated_seconds,
                shared_seconds=shared_seconds,repeated_term_products=int(products),
                shared_term_products=int(C.products),coefficients=C.generated,
                subsystems=len(shared),exact_coefficients_equal=True)


def optimistic_oracle(fs, degree=4, limit=32):
    """Diagnostic only: candidate generation/search are deliberately free.

    Exhaustively try individual candidates up to limit, plus several batches.
    Choose by an operation proxy (not by lucky measured timing). This is NOT an
    implementable speedup or a proof of the globally optimal coefficient set.
    """
    start=time.perf_counter();n=fs[0].parent().ngens();C=SharedDixon(fs[:n-1],degree)
    candidates=[]
    for f in fs[n-1:]:
        pre=C.prepare(f)
        for b in pre[1]:
            g=C.coefficient(pre,b)
            if g:candidates.append(g)
    # Only independently useful relations relative to the original equations.
    S=Space(fs[0].parent());S.extend(exponents(n,degree));S.insert(fs)
    unique=[]
    for g in candidates:
        if S.insert([g]):unique.append(g)
    if len(unique)>limit:
        unique=[unique[i*(len(unique)-1)//(limit-1)] for i in range(limit)]
    selections=[[]]+[[i] for i in range(len(unique))]
    for size in (2,4,8,16,len(unique)):
        if size<=len(unique):
            selections += [list(range(size)),list(range(len(unique)-size,len(unique)))]
    results=[]
    schedules=[(0,s) for s in selections]
    for after in (1,2):
        schedules += [(after,s) for s in selections if len(s)!=1 or s[0]%4==0]
    for after,indices in schedules:
        r=scalar_inject(fs,degree,[unique[i] for i in indices],after)
        if results:
            assert r['nullity']==results[0][2]['nullity']
            assert r['recovered']==results[0][2]['recovered']
        score=r['dense_matmul_multiply_adds']+r['rref_cells']
        results.append((score,indices,r))
    score,indices,best=min(results,key=lambda x:(x[0],len(x[1])))
    baseline=results[0][2]
    return dict(seconds=time.perf_counter()-start,diagnostic_only=True,
                candidates=len(candidates),independent_candidates=len(unique),trials=len(results),
                selected=indices,best=best,baseline=baseline,
                free_generation_work_ratio=score/(baseline['dense_matmul_multiply_adds']+baseline['rref_cells']),
                extra_polynomials=[str(unique[i]) for i in indices])


def solve_linearized(fs, block=2):
    """Exact variable-block elimination, charging both determinant branches.

    The remaining systems use Singular GB, so this measures the preprocessing
    route rather than claiming a new Dixon implementation for higher degrees.
    Only return a complete unique-point certificate if det(A)=0 is inconsistent.
    """
    start=time.perf_counter();P=fs[0].parent();n=P.ngens();F=P.base_ring();xs=P.gens()
    pure=[xs[i]*xs[j] for i in range(block) for j in range(i,block)]
    Q=matrix(F,[[f.monomial_coefficient(u) for u in pure] for f in fs])
    K=Q.left_kernel().basis_matrix();linear=[sum(c*f for c,f in zip(row,fs)) for row in K.rows()]
    if len(linear)<block:
        return dict(seconds=time.perf_counter()-start,certified=False,reason='not enough linear combinations')
    A=matrix(P,[[f.derivative(xs[j]) for j in range(block)] for f in linear])
    selected=list(A.transpose().pivots())[:block]
    if len(selected)<block:
        return dict(seconds=time.perf_counter()-start,certified=False,reason='linear block rank deficient')
    A=A.matrix_from_rows(selected);det=A.det()
    if not det:
        return dict(seconds=time.perf_counter()-start,certified=False,reason='singular linear block')
    b=vector(P,[linear[i].subs({x:0 for x in xs[:block]}) for i in selected])
    numerators=list(-A.adjugate()*b)
    R=PolynomialRing(F,names=tuple('v%d'%i for i in range(n-block))+('s',))
    vs=R.gens();embed=P.hom([R.zero()]*block+list(vs[:-1]),R)
    ds=embed(det);ns=[embed(a) for a in numerators]
    transformed=[]
    for f in fs:
        h=R.zero()
        for e,c in f.dict().items():
            du=sum(e[:block]);h+=c*ds**(2-du)*prod(a**k for a,k in zip(ns,e[:block]))*prod(a**k for a,k in zip(vs,e[block:]))
        if h:transformed.append(h)
    maxdegree=max(int(f.total_degree()) for f in transformed)
    generic=R.ideal(transformed+[vs[-1]*ds-1]);G=generic.groebner_basis()
    exceptional=P.ideal(fs+[det]);H=exceptional.groebner_basis()
    empty_exceptional=any(h==1 for h in H)
    if any(g==1 for g in G):
        gd=0;point=None
    else:
        B=list(generic.normal_basis());gd=len(B);point=None
        if gd==1:
            vv=[F(generic.reduce(v)) for v in vs]
            dd=ds(*vv)
            point=[a(*vv)/dd for a in ns]+vv[:-1]
            assert all(f(*point)==0 for f in fs)
    certified=empty_exceptional and gd<=1
    return dict(seconds=time.perf_counter()-start,certified=certified,
                dimension=gd if empty_exceptional else None,
                recovered=list(map(int,point)) if certified and point else None,
                pure_quadratic_rank=int(Q.rank()),linear_combinations=K.nrows(),block=block,
                remaining_variables=n-block,transformed_degree=maxdegree,
                determinant_degree=int(det.total_degree()),exceptional_empty=empty_exceptional,
                generic_dimension=gd,exceptional_basis_size=len(H))

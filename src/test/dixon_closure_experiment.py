"""Compare original-only relation closure with partial Dixon seeding.

Run with Sage's Python. Same exact scalar incremental RREF closure for every
method. All construction costs are timed. No planted coordinates or GB in the
solver. Audit uses a separately constructed full Dixon polynomial / Macaulay
row space, outside timing. The C benchmark consumes the same polynomial files.
"""
from sage.all import GF, PolynomialRing, matrix, prod
from pathlib import Path
from dixon_closure_cases import make_system
import argparse, json, random, statistics, subprocess, tempfile, time


def exponents(n, degree):
    if n == 0:
        yield ()
    else:
        for i in range(degree+1):
            for tail in exponents(n-1, degree-i):
                yield (i,)+tail


def partial_dixon(fs, relation_degree, full=False):
    """Adjacent differences first, original row last; prune during subset DP.

    For n quadratic inputs, a coefficient of auxiliary degree b has a
    representation sum a_i f_i with deg(a_i) <= n-1-b. Select degrees at
    most relation_degree and omit b=n-1 (constant combinations of inputs).
    """
    start = time.perf_counter()
    P = fs[0].parent(); n = P.ngens(); k = n-1
    Q = PolynomialRing(P.base_ring(), names=P.variable_names()+tuple('y%d'%i for i in range(k)), order='degrevlex')
    z = Q.gens(); embed = P.hom(z[:n], Q)
    original = [embed(f) for f in fs[:n]]; previous = original
    rows = []
    for i in range(k):
        current = [f.subs({z[j]:z[n+j] for j in range(i+1)}) for f in original]
        rows.append([(a-b)//(z[i]-z[n+i]) for a,b in zip(previous,current)])
        previous = current
    rows.append(original)
    lo, hi = (0,k) if full else (max(0,n+1-relation_degree), k-1)
    dp = {0:Q.one()}; products = 0; peak = 0
    for mask in range(1, 1<<n):
        r = mask.bit_count()-1; cols = [j for j in range(n) if mask>>j&1]
        g = Q.zero()
        for j,c in enumerate(cols):
            a,b = rows[r][c],dp[mask^(1<<c)]
            products += len(a.dict())*len(b.dict())
            g += (-1)**(r+j)*a*b
        remaining = max(0,k-r-1)
        if not full:
            g = Q({e:c for e,c in g.dict().items() if lo-remaining <= sum(e[n:]) <= hi})
        dp[mask] = g; peak = max(peak,len(g.dict()))
    groups = {}
    for e,c in dp[(1<<n)-1].dict().items():
        groups.setdefault(tuple(e[n:]),{})[tuple(e[:n])] = c
    result = {b:P(d) for b,d in sorted(groups.items())}
    return result, dict(seconds=time.perf_counter()-start, rows=len(result),
                        terms=sum(len(g.dict()) for g in result.values()),
                        candidate_term_products=products, peak_polynomial_terms=peak)


class Closure:
    def __init__(self, P, degree):
        self.P = P; self.F = P.base_ring(); self.n = P.ngens(); self.degree = degree
        self.ee = sorted(exponents(self.n,degree), key=lambda e:(sum(e),e), reverse=True)
        self.ix = {e:i for i,e in enumerate(self.ee)}; self.nc = len(self.ee)
        self.E = matrix(self.F,0,self.nc)
        self.processed = matrix(self.F,0,self.nc)
        self.history = []; self.rref_cells = 0; self.matmul_work = 0
        self.rref_shapes = []; self.submitted_rows = 0
        self.inverse = [[self.ix.get(tuple(a-(j==v) for j,a in enumerate(e)),-1) if e[v] else -1 for e in self.ee] for v in range(self.n)]

    def coefficients(self, polys):
        M = matrix(self.F,len(polys),self.nc)
        for i,f in enumerate(polys):
            for e,c in f.dict().items(): M[i,self.ix[tuple(e)]] = c
        return M

    def rref(self, M):
        self.rref_cells += M.nrows()*M.ncols()
        self.rref_shapes.append([int(M.nrows()),int(M.ncols())])
        E = M.echelon_form()
        return E[:E.rank(),:]

    def residual(self, M, E):
        if E.nrows() == 0: return M
        self.matmul_work += M.nrows()*E.nrows()*E.ncols()
        return M-M.matrix_from_columns(list(E.pivots()))*E

    def insert(self, M):
        self.submitted_rows += M.nrows()
        if not self.E.nrows():
            self.E = self.rref(M); return
        piv = list(self.E.pivots()); ps = set(piv)
        free = [j for j in range(self.nc) if j not in ps]
        self.matmul_work += M.nrows()*len(piv)*len(free)
        residual = M.matrix_from_columns(free)-M.matrix_from_columns(piv)*self.E.matrix_from_columns(free)
        S = self.rref(residual)
        if not S.nrows(): return
        newp = [free[j] for j in S.pivots()]
        lift = matrix(self.F,S.nrows(),self.nc)
        for i,row in enumerate(S.rows()):
            for j,c in zip(free,row): lift[i,j] = c
        self.matmul_work += self.E.nrows()*len(newp)*self.nc
        updated = self.E-self.E.matrix_from_columns(newp)*lift
        joined = updated.stack(lift); allp = piv+newp
        self.E = joined.matrix_from_rows(sorted(range(len(allp)),key=allp.__getitem__))

    def run(self, fs, initial):
        self.insert(self.coefficients(initial))
        candidate = None
        while True:
            self.history.append(int(self.E.nrows()))
            # A complete linear certificate, or 1 in the ideal, permits early exit.
            if self.E.nrows() == self.nc: break
            if self.E.nrows() == self.nc-1 and self.E.pivots()[-1] == self.nc-2:
                candidate = [-self.E[self.E.pivots().index(self.ix[tuple(int(i==j) for i in range(self.n))]),self.nc-1] for j in range(self.n)]
                assert all(f(*candidate)==0 for f in fs)
                break
            low = self.E.matrix_from_rows([i for i,p in enumerate(self.E.pivots()) if sum(self.ee[p]) < self.degree])
            delta = self.rref(self.residual(low,self.processed))
            self.processed = low
            if not delta.nrows(): break
            extra = matrix(self.F,[[row[j] if j>=0 else self.F.zero() for j in perm]
                                  for row in delta.rows() for perm in self.inverse])
            rank = self.E.nrows(); self.insert(extra)
            if self.E.nrows() == rank: break
        return dict(rank_history=self.history, nullity=int(self.nc-self.E.nrows()),
                    recovered=[int(c) for c in candidate] if candidate is not None else None,
                    columns=self.nc, submitted_rows=int(self.submitted_rows),
                    rref_cells=int(self.rref_cells), rref_shapes=self.rref_shapes,
                    dense_matmul_multiply_adds=int(self.matmul_work))


def original_multiples(fs, degree):
    P = fs[0].parent(); xs = P.gens()
    return [prod(x**a for x,a in zip(xs,e))*f for f in fs
            for e in exponents(P.ngens(),degree-int(f.total_degree())) if f]


def run(fs, degree, method):
    start = time.perf_counter(); construction = {}; partial = {}
    initial = fs[:]
    if method.startswith('hybrid') or method.startswith('augment'):
        rd = int(method[-1]); assert rd <= degree
        partial,construction = partial_dixon(fs,rd)
    if method in ('macaulay',) or method.startswith('augment'):
        initial = original_multiples(fs,degree)
    initial += list(partial.values())
    cl = Closure(fs[0].parent(),degree)
    if method == 'macaulay' or method.startswith('augment'):
        # These rows have already been prolonged in the degree-D input.
        # Avoid artificially weakening the full-Macaulay baseline by doing
        # all those multiplications a second time.
        cl.processed = cl.rref(cl.coefficients(original_multiples(fs,degree-1)))
    built = time.perf_counter()
    result = cl.run(fs,initial); done = time.perf_counter()
    result.update(method=method,degree=degree,seconds=done-start,setup_seconds=built-start,
                  closure_seconds=done-built,dixon=construction,initial_rows=len(initial))
    return result,cl,partial


def rowspace_hash(M):
    h = 14695981039346656037
    for row in M.rows():
        for c in row: h = ((h ^ int(c))*1099511628211) & ((1<<64)-1)
    return '%016x'%h


def audit(fs, degree, partial, cl):
    # Check exact membership in the ordinary, pre-closure Macaulay row space.
    M = cl.coefficients(original_multiples(fs,degree)); E = M.echelon_form(); E = E[:E.rank(),:]
    G = cl.coefficients(list(partial.values()))
    residual = G-G.matrix_from_columns(list(E.pivots()))*E
    assert residual.is_zero(), 'Dixon degree bound / coefficient assembly failure'
    return dict(macaulay_rank=int(E.nrows()), added_rank=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--variables',type=int,default=5); ap.add_argument('--extra',type=int,default=1)
    ap.add_argument('--degree',type=int,default=4); ap.add_argument('--seeds',type=int,default=3)
    ap.add_argument('--repeats',type=int,default=1); ap.add_argument('--points',type=int,default=1)
    ap.add_argument('--prime',type=int,default=257); ap.add_argument('--family',default='dense')
    ap.add_argument('--methods',default='seed,hybrid3,hybrid4,macaulay,augment3')
    ap.add_argument('--audit',action='store_true'); ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--export',type=Path)
    ap.add_argument('--c-binary',type=Path)
    ap.add_argument('--native-repeats',type=int,default=5)
    ap.add_argument('--matrix-audit',action='store_true',help='untimed full native RREF comparison and full Dixon comparison through n=7')
    args = ap.parse_args(); n = args.variables
    P = PolynomialRing(GF(args.prime),names=tuple('x%d'%i for i in range(n)),order='degrevlex')
    methods = args.methods.split(','); records = []
    if args.c_binary and not args.export: ap.error('--c-binary requires --export')
    if args.matrix_audit and not args.c_binary: ap.error('--matrix-audit requires --c-binary')
    warm,_ = make_system(P,999,n+args.extra,args.family,args.points)
    for method in methods: run(warm,args.degree,method)
    for seed in range(args.seeds):
        fs,points = make_system(P,seed,n+args.extra,args.family,args.points)
        if args.export:
            args.export.mkdir(parents=True,exist_ok=True)
            # Native benchmark format: prime n m, then m parseable polynomials.
            (args.export/('seed%d.txt'%seed)).write_text('%d %d %d\n'%(args.prime,n,len(fs))+'\n'.join(str(f) for f in fs)+'\n')
        for repeat in range(args.repeats):
            order = methods[:]; random.Random(seed*100+repeat).shuffle(order)
            final_spaces = []
            for method in order:
                result,cl,partial = run(fs,args.degree,method)
                assert all(all(f(*a)==0 for f in fs) for a in points)
                if result['recovered'] is not None:
                    assert len(points)==1 and result['recovered']==list(map(int,points[0]))
                result['rowspace_hash'] = rowspace_hash(cl.E)
                if args.audit:
                    result['audit'] = audit(fs,args.degree,partial,cl)
                    final_spaces.append(cl.E)
                    if partial and (n<=5 or (args.matrix_audit and n<=7)):
                        full,_ = partial_dixon(fs,n+1,full=True)
                        assert all(full[b]==g for b,g in partial.items())
                result.update(seed=seed,repeat=repeat,n=n,m=len(fs),prime=args.prime,points=args.points,family=args.family)
                records.append(result)
                print(json.dumps(result),flush=True)
                if args.c_binary:
                    # Sequential runs avoid CPU contention. Warm the process
                    # once, then take medians of the measured invocations.
                    cmd = [str(args.c_binary),str(args.export/('seed%d.txt'%seed)),method,str(args.degree),str(int(args.audit)),str(1+args.native_repeats)]
                    native_runs = [json.loads(line) for line in subprocess.check_output(cmd,text=True).splitlines()][1:]
                    for native in native_runs:
                        for key in ('rank_history','nullity','recovered','rowspace_hash','submitted_rows','rref_cells','dense_matmul_multiply_adds'):
                            assert native[key]==result[key], (method,seed,key,native[key],result[key])
                    native = dict(native_runs[-1])
                    if partial:
                        for ckey,skey in [('dixon_terms','terms'),('dixon_rows','rows'),('candidate_term_products','candidate_term_products')]:
                            assert native[ckey]==result['dixon'][skey], (ckey,native[ckey],result['dixon'][skey])
                    timing_keys = ('seconds','setup_seconds','closure_seconds','dixon_seconds','cpu_seconds','cpu_setup_seconds','cpu_closure_seconds','cpu_dixon_seconds')
                    native['timing_samples'] = [{k:r[k] for k in timing_keys if k in r} for r in native_runs]
                    for key in native['timing_samples'][0]: native[key] = statistics.median(r[key] for r in native_runs)
                    native['repeat'] = 'median after one warmup'
                    result['c'] = native
                    if args.matrix_audit:
                        with tempfile.TemporaryDirectory() as directory:
                            dump=Path(directory)/'rref.txt'
                            subprocess.check_output(cmd[:-1]+['1',str(dump)],text=True)
                            data=list(map(int,dump.read_text().split()))
                            assert matrix(P.base_ring(),data[0],data[1],data[2:])==cl.E
                        result['native_matrix_coefficients_equal']=True
            if args.audit: assert all(E==final_spaces[0] for E in final_spaces)
            args.output.write_text(json.dumps(records,indent=2)+'\n')
    for method in methods:
        rr = [r for r in records if r['method']==method]
        print(method,'median_seconds',statistics.median(r['seconds'] for r in rr),flush=True)


if __name__ == '__main__': main()

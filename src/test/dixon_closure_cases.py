"""Reproducible planted inputs for the partial Dixon/closure comparison."""
import random


def make_system(P, seed, count, family="dense", points=1):
    F, xs = P.base_ring(), P.gens()
    rng = random.Random(seed)
    elements = list(F) if F.cardinality() < 1024 else None
    def coefficient():
        return elements[rng.randrange(len(elements))] if elements else F(rng.randrange(F.cardinality()))
    planted = [tuple(coefficient() for _ in xs)]
    if points == 2:
        b = list(planted[0]); b[0] += 1
        planted.append(tuple(b))
    quadratics = [x*y for i, x in enumerate(xs) for y in xs[i:]]
    common = sum(coefficient()*u for u in quadratics)
    fs = []
    for i in range(count):
        q = common if family == "shared" or (family == "partial" and i < 2) else sum(coefficient()*u for u in quadratics)
        f = P(q + sum(coefficient()*x for x in xs))
        if points == 2:
            a, b = planted
            f -= (f(*b)-f(*a))/(b[0]-a[0])*xs[0]
        f -= f(*planted[0])
        fs.append(f)
    return fs, planted

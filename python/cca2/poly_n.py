"""poly_n -- exact multivariate integer polynomials (the N-D coefficient arrays of ratQ.m).

MATLAB stores a polynomial in nv variables as an nv-dimensional array A with
A(i1+1, ..., inv+1) the coefficient of x1^i1 ... xnv^inv, and keeps a second BigInteger
representation for values past 2^53. The port has ONE sparse representation:

    a dict  {(i1, ..., inv): coefficient}   with nonzero int coefficients; {} is zero.

toND / fromND convert to and from MATLAB's array (size vector + column-major values), which is
what the goldens and traces carry.
"""
from math import gcd

from .errors import error


def mono(e, c=1):
    return {tuple(e): c} if c != 0 else {}


def const(c, nv):
    return {(0,) * nv: c} if c != 0 else {}


def add(a, b):
    out = dict(a)
    for e, c in b.items():
        v = out.get(e, 0) + c
        if v == 0:
            out.pop(e, None)
        else:
            out[e] = v
    return out


def neg(a):
    return {e: -c for e, c in a.items()}


def sub(a, b):
    return add(a, neg(b))


def scale(a, k):
    if k == 0:
        return {}
    return {e: c * k for e, c in a.items()}


def mul(a, b):
    out = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = tuple(x + y for x, y in zip(ea, eb))
            out[e] = out.get(e, 0) + ca * cb
    return {e: c for e, c in out.items() if c != 0}


def content(a):
    g = 0
    for c in a.values():
        g = gcd(g, abs(c))
        if g == 1:
            return 1
    return g


def contentFree(a):
    """Divide out the (positive) gcd of the coefficients; the SIGN is left alone."""
    g = content(a)
    if g > 1:
        return {e: c // g for e, c in a.items()}
    return dict(a)


def degrees(a, nv):
    """Max exponent per variable (MATLAB: size(A) - 1)."""
    if not a:
        return [0] * nv
    return [max(e[d] for e in a) for d in range(nv)]


def maxCoeff(a):
    return max((abs(c) for c in a.values()), default=0)


def leadTerm(a):
    """(exponent, coefficient) of the lex-greatest term, first variable most significant."""
    if not a:
        return None, 0
    e = max(a)
    return e, a[e]


def divExact(P, D):
    """P / D when the division is exact; raises ratQ:notDivisible otherwise."""
    eD, cD = leadTerm(D)
    if eD is None:
        error('ratQ:divideByZero', 'the divisor is the zero polynomial.')
    R = dict(P)
    Q = {}
    guard = 0
    while R:
        eR, cR = leadTerm(R)
        guard += 1
        if guard > 100000:
            error('ratQ:divideNoConverge', 'exact division did not terminate.')
        eQ = tuple(x - y for x, y in zip(eR, eD))
        if any(x < 0 for x in eQ) or cR % cD != 0:
            error('ratQ:notDivisible',
                  'the division is not exact: the remainder\'s leading term is not a multiple '
                  'of the divisor\'s. Nothing is truncated or rounded here.')
        T = {eQ: cR // cD}
        Q = add(Q, T)
        R = sub(R, mul(T, D))
    return Q


def slices(a, d):
    """{k: polynomial in the other variables} -- the coefficients of variable d (0-based)."""
    out = {}
    for e, c in a.items():
        k = e[d]
        rest = e[:d] + e[d + 1:]
        out.setdefault(k, {})[rest] = c
    return out


def fromSlices(sl, d):
    out = {}
    for k, p in sl.items():
        for e, c in p.items():
            out[e[:d] + (k,) + e[d:]] = c
    return out


def evalAt(a, x):
    """Exact value at the point x (ints or Fractions)."""
    v = 0
    for e, c in a.items():
        t = c
        for xi, ei in zip(x, e):
            if ei:
                t = t * xi ** ei
        v += t
    return v


def toND(a, nv):
    """(size, column-major values) of MATLAB's array; at least 2 dimensions, as MATLAB reports."""
    deg = degrees(a, nv)
    sz = [d + 1 for d in deg]
    while len(sz) > 2 and sz[-1] == 1:
        sz.pop()
    if len(sz) < 2:
        sz.append(1)
    n = 1
    for s in sz:
        n *= s
    v = [0] * n
    for e, c in a.items():
        idx = 0
        m = 1
        for d in range(len(sz)):
            idx += e[d] * m
            m *= sz[d]
        v[idx] = c
    return sz, v


def fromND(sz, v, nv):
    out = {}
    sz = list(sz) + [1] * (nv - len(sz))
    for i, c in enumerate(v):
        if c == 0:
            continue
        e = []
        r = i
        for d in range(len(sz)):
            e.append(r % sz[d])
            r //= sz[d]
        out[tuple(e[:nv])] = c
    return out


def fromMatrix(M, nv=2):
    """A MATLAB 2-D coefficient array (list of rows) as a polynomial in nv >= 2 variables."""
    out = {}
    for i, row in enumerate(M):
        for j, c in enumerate(row):
            if c != 0:
                out[(i, j) + (0,) * (nv - 2)] = c
    return out

"""isolated_root_sign -- the exact sign of N(alpha) where alpha is the root of R isolated by
iv = (lo, hi] (port of isolatedRootSign.m).

Unlike rat_q.signAt it does not require R to be squarefree-reduced first and it meets RATIONAL
roots exactly; polynomials are integer lists, highest degree first.
"""
from fractions import Fraction

from . import rat_q
from .errors import error
from .rat_q import sign, polyval


def _strip(P):
    P = list(P)
    while P and P[0] == 0:
        P = P[1:]
    return P


def _prem(a, b):
    """MATLAB ratQ.bigPrem: a * |lc(b)|^(da-db+1) mod b (no content removal)."""
    da = len(a) - 1
    db = len(b) - 1
    if da < db:
        return list(a)
    lc = b[0]
    f = abs(lc) ** (da - db + 1)
    r = [c * f for c in a]
    while len(r) - 1 >= db:
        q = r[0] // lc
        for j in range(len(b)):
            r[j] -= q * b[j]
        r = r[1:]
        while r and r[0] == 0:
            r = r[1:]
        if not r:
            return r
    return r


def _bgcd(a, b):
    a = rat_q.polyPrimitive(a)
    b = rat_q.polyPrimitive(b)
    while True:
        r = _strip(_prem(a, b))
        if not r:
            return b
        if len(r) == 1:
            return [1]
        a = b
        b = rat_q.polyPrimitive(r)


def isolatedRootSign(N, R, iv):
    N = _strip(N)
    if not N:
        return 0
    Rb = _strip(R)
    lo = Fraction(iv[0])
    hi = Fraction(iv[1])
    if sign(polyval(Rb, hi)) == 0:
        return sign(polyval(N, hi))                 # the root IS hi, a rational
    while sign(polyval(Rb, lo)) == 0:               # lo is ANOTHER root: move it up
        mid = (lo + hi) / 2
        if sign(polyval(Rb, mid)) == 0:
            return sign(polyval(N, mid))
        if rat_q.countRootsIn(Rb, mid, hi) == 1:
            lo = mid
        else:
            hi = mid
    r = _strip(_prem(N, Rb))
    if not r:
        return 0
    if len(r) == 1:
        return sign(r[0])
    g = _bgcd(Rb, r)
    if len(g) > 1 and rat_q.countRootsIn(g, lo, hi) > 0:
        return 0
    for _ in range(400):
        if rat_q.countRootsIn(r, lo, hi) == 0 and sign(polyval(r, lo)) != 0:
            return sign(polyval(r, (lo + hi) / 2))
        mid = (lo + hi) / 2
        if sign(polyval(Rb, mid)) == 0:
            return sign(polyval(r, mid))            # a RATIONAL root, met exactly
        if rat_q.countRootsIn(Rb, lo, mid) == 1:
            hi = mid
        else:
            lo = mid
    error('isolatedRootSign:noConverge',
          '400 bisections did not separate the root from the zeros of the reduced polynomial; '
          'with the zero test above that is a defect.')

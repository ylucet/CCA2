"""dual_point -- the EXACT name of an isolated dual point: a common point of two conics
(port of dualPoint.m and dualPointSign.m).

After a shear s1 = u + k v, s2 = v the conics are quadratics in v with coefficients in u; their
resultant R(u) has the point's u as a root, isolated by a dyadic interval iv, and v is then the
rational function -S0(u)/S1(u). The name is {R, iv, S1, S0, k}; `s` is its double cache.
dualPointSign gives the exact sign of an integer polynomial Q(s1, s2) at the point.
"""
import math
from fractions import Fraction

from . import rat_q
from .errors import error
from .isolated_root_sign import isolatedRootSign
from .rat_q import sign, polyval


def _padd(a, b):
    n = max(len(a), len(b))
    a = [0] * (n - len(a)) + list(a)
    b = [0] * (n - len(b)) + list(b)
    return [x + y for x, y in zip(a, b)]


def _scal(a, k):
    return [x * k for x in a]


def _pmul(a, b):
    return rat_q.polyMul1(a, b)


def _inV(C, k):
    a, b, c, d, e, f = C
    return a * k * k + b * k + c, [2 * a * k + b, d * k + e], [a, d, f]


def dualPoint(c1, c2, sApprox):
    A = rat_q.conic(c1)
    B = rat_q.conic(c2)
    s = [float(sApprox[0]), float(sApprox[1])]
    if not all(math.isfinite(v) for v in s):
        error('dualPoint:approx', 'sApprox must be a finite 2-vector; got %r.' % (sApprox,))
    if all(x == 0 for x in A[0:3]) and all(x == 0 for x in B[0:3]):
        error('dualPoint:twoLines',
              'both curves are lines, so their crossing is RATIONAL and needs no algebraic '
              'name: %s and %s.' % (A, B))
    for k in range(13):
        A2, A1, A0 = _inV(A, k)
        B2, B1, B0 = _inV(B, k)
        if A2 == 0 and B2 == 0:
            continue
        X = _padd(_scal(B0, A2), _scal(A0, -B2))
        Y = _padd(_scal(B1, A2), _scal(A1, -B2))
        Z = _padd(_pmul(A1, B0), _scal(_pmul(B1, A0), -1))
        R = rat_q.polyTrim(_padd(_pmul(X, X), _scal(_pmul(Y, Z), -1)))
        if not R:
            error('dualPoint:commonComponent',
                  'the conics %s and %s share a component, so they meet in a curve and not in '
                  'a point.' % (A, B))
        R = rat_q.polyPrimitive(rat_q.polySquarefree(R))
        S1 = rat_q.polyTrim(_padd(_scal(A1, B2), _scal(B1, -A2)))
        S0 = rat_q.polyTrim(_padd(_scal(A0, B2), _scal(B0, -A2)))
        if not S1:
            continue
        u = s[0] - k * s[1]
        iv = _isolateNear(R, u)
        if iv is None:
            error('dualPoint:noMatch',
                  'no real root of the eliminant of %s and %s lies near u = %.17g; the '
                  'approximation is not a common point of the two conics.' % (A, B, u))
        if isolatedRootSign(S1, R, iv) == 0:
            continue                               # two roots on one vertical line
        u0 = _narrow(R, iv)
        v0 = -polyval([float(c) for c in S0], u0) / polyval([float(c) for c in S1], u0)
        s0 = [u0 + k * v0, v0]
        if math.hypot(s0[0] - s[0], s0[1] - s[1]) > 1e-6 * max(1.0, math.hypot(s[0], s[1])):
            error('dualPoint:noMatch',
                  'no common root of %s and %s lies near %r; the nearest is %r.' % (A, B, s, s0))
        return {'c1': A, 'c2': B, 'k': k, 'R': R, 'iv': iv, 'S1': S1, 'S0': S0, 's': s0}
    error('dualPoint:noShear',
          'no shear k = 0..12 gives the conics %s and %s a second coordinate determined at the '
          'root near %r.' % (A, B, s))


def _isolateNear(R, u):
    sc = max(1.0, abs(u))
    for e in range(-48, 5):
        d = sc * 2.0 ** e
        lo = Fraction(u - d)
        hi = Fraction(u + d)
        if polyval(R, lo) == 0 or polyval(R, hi) == 0:
            continue
        n = rat_q.countRootsIn(R, lo, hi)
        if n == 1:
            return (lo, hi)
        if n > 1:
            break
    return None


def _narrow(R, iv):
    lo = float(iv[0])
    hi = float(iv[1])
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid <= lo or mid >= hi:
            break
        if polyval(R, Fraction(mid)) == 0:
            return mid                             # a RATIONAL root, met exactly
        if rat_q.countRootsIn(R, Fraction(lo), Fraction(mid)) == 1:
            hi = mid
        else:
            lo = mid
    return (lo + hi) / 2


def dualPointSign(dp, Q):
    """The exact sign of Q(s1, s2) = sum Q[i][j] s1^i s2^j at the dual point dp."""
    Q = [list(r) for r in Q]
    rat_q.chk(Q, 'polynomial')
    if not any(x != 0 for r in Q for x in r):
        return 0
    k = dp['k']
    ni = len(Q)
    nj = len(Q[0])
    m = (ni - 1) + (nj - 1)
    Quv = [[0] * (m + 1) for _ in range(ni)]       # Quv[a][b]: coefficient of u^a v^b
    for i in range(ni):
        for j in range(nj):
            if Q[i][j] == 0:
                continue
            for t in range(i + 1):
                Quv[i - t][t + j] += Q[i][j] * math.comb(i, t) * k ** t
    while len(Quv[0]) > 1 and not any(r[-1] != 0 for r in Quv):
        Quv = [r[:-1] for r in Quv]
    m = len(Quv[0]) - 1

    def ppow(A, n):
        P = [1]
        for _ in range(n):
            P = _pmul(P, A)
        return P
    negS0 = [-x for x in dp['S0']]
    S1 = dp['S1']
    N = []
    for j in range(m + 1):
        Qj = rat_q.polyTrim([Quv[a][j] for a in range(ni - 1, -1, -1)])   # highest degree first
        if not Qj:
            continue
        N = _padd(N, _pmul(_pmul(Qj, ppow(negS0, j)), ppow(S1, m - j)))
    sN = isolatedRootSign(N, dp['R'], dp['iv'])
    if sN == 0 or m % 2 == 0:
        return sN
    return sN * isolatedRootSign(S1, dp['R'], dp['iv'])

"""conic_meet -- the real intersection points of two integer conics (port of conicMeet.m).

The conics are exact; the POINTS are doubles on purpose: a vertex of f* is generically irrational
(degree <= 4), so its coordinates are a rebuildable cache, never stored data. What is exact here is
the elimination -- the resultant quartic is built in integers -- and the canonical ORDER of the
returned points (sorted rows), which is what makes "root r of conics i and j" a NAME.
"""
import math

import numpy as np

from . import rat_q
from .errors import error, CCA2Error

EPS = 2.0 ** -52


_memo = {}


def clearMemo():
    _memo.clear()


def conicMeet(A, B, opts=None):
    """Returns (P, info): P a list of [x, y] doubles in canonical (sorted) order.

    A pure function of (A, B, opts), and the arrangement code asks for the same pair many times
    (conjQ on the benchmark: 79189 calls, 4358 distinct pairs), so results are remembered per
    ordered pair (and rat_q.EMULATE_2P53, which changes what the elimination may compute) and
    handed out as deep copies. Refusals are not remembered (they re-raise)."""
    opts = dict(opts or {})
    opts.setdefault('tol', 1e-9)
    opts.setdefault('polish', 8)
    # the key is the RAW pair when both are plain-int sequences (so validation is skipped on a hit:
    # a pair of ints that was accepted once is accepted again); anything else is validated first
    raw = all(type(x) is int for x in A) and all(type(x) is int for x in B)
    if not raw:
        A = rat_q.conic(A)
        B = rat_q.conic(B)
    try:
        key = (tuple(A), tuple(B), tuple(sorted(opts.items())), rat_q.EMULATE_2P53)
        hit = _memo.get(key)
    except TypeError:                               # an unhashable option: no memo
        return _conicMeet(rat_q.conic(A), rat_q.conic(B), opts)
    if hit is None:
        hit = _conicMeet(rat_q.conic(A), rat_q.conic(B), opts)
        if len(_memo) > 50000:
            _memo.clear()
        _memo[key] = hit
    return _copy(hit)


def _copy(hit):
    """Copy of a remembered (P, info): P is a list of [x, y] floats and info's values are scalars
    or flat lists of scalars, so a two-level copy is a deep copy."""
    P, info = hit
    return [p[:] for p in P], {k: (v[:] if type(v) is list else v) for k, v in info.items()}


def _conicMeet(A, B, opts):
    info = {'quartic': [], 'shear': 0, 'resid': [], 'sep': math.inf, 'degenerate': False}
    if A == B:
        info['degenerate'] = True
        return [], info
    if all(x == 0 for x in A[0:5]) or all(x == 0 for x in B[0:5]):
        bad = A if all(x == 0 for x in A[0:5]) else B
        error('conicMeet:constant',
              'one argument has no x or y in it, so it is a CONSTANT and not a curve: %s. Drop '
              'such a row before building the arrangement -- it is vacuous when it holds and '
              'empties the cell when it does not.' % (bad,))
    if all(x == 0 for x in A[0:3]) or all(x == 0 for x in B[0:3]):
        if all(x == 0 for x in A[0:3]):
            L, C = A, B
        else:
            L, C = B, A
        P, q = _lineMeetsConic(L, C, opts)
        info['quartic'] = [0, 0] + q
    else:
        P, info = _generalMeet(A, B, info, opts)
    if P:
        P = sorted([[_round12(p[0]), _round12(p[1])] for p in P])
        P, info['resid'] = _polishAndScore(P, A, B, opts)
        if len(P) > 1:
            info['sep'] = min(math.hypot(P[i][0] - P[j][0], P[i][1] - P[j][1])
                              for i in range(len(P) - 1) for j in range(i + 1, len(P)))
    return P, info


def _round12(x):
    """MATLAB round(x, 12): to 12 decimals, halves away from zero. Only an ORDERING key."""
    if not math.isfinite(x):
        return x
    y = abs(x) * 1e12
    if y >= 2.0 ** 52:
        return x
    r = math.floor(y + 0.5) / 1e12
    return r if x >= 0 else -r


def _lineMeetsConic(L, C, opts):
    d, e, f = L[3], L[4], L[5]
    if abs(d) >= abs(e):
        P0 = [-f, 0]
        scale = d
    else:
        P0 = [0, -f]
        scale = e
    u = [-e, d]
    q2 = C[0] * u[0] ** 2 + C[1] * u[0] * u[1] + C[2] * u[1] ** 2
    q1 = (2 * C[0] * P0[0] * u[0] + C[1] * (P0[0] * u[1] + P0[1] * u[0]) + 2 * C[2] * P0[1] * u[1]
          + scale * (C[3] * u[0] + C[4] * u[1]))
    q0 = (C[0] * P0[0] ** 2 + C[1] * P0[0] * P0[1] + C[2] * P0[1] ** 2
          + scale * (C[3] * P0[0] + C[4] * P0[1]) + scale ** 2 * C[5])
    q = rat_q.wall([q2, q1, q0], None, 'line-conic quadratic')
    t = _realRootsExact(q)
    P = [[P0[0] / scale + ti * u[0], P0[1] / scale + ti * u[1]] for ti in t]
    P = [p for p in P if math.isfinite(p[0]) and math.isfinite(p[1])]
    if len(P) == 2 and math.hypot(P[0][0] - P[1][0], P[0][1] - P[1][1]) <= \
            opts['tol'] * max(1.0, math.hypot(P[0][0], P[0][1])):
        P = P[0:1]                                  # a tangency is ONE point, not two
    return P, q


def _generalMeet(A, B, info, opts):
    k = _shearFor(A, B)
    info['shear'] = k
    As = _shearConic(A, k)
    Bs = _shearConic(B, k)
    q = _resultantY(As, Bs)
    info['quartic'] = q
    if all(x == 0 for x in q):
        info['degenerate'] = True
        return [], info
    xs = _realRootsExact(q)
    P = []
    for x in xs:
        ys = _realRoots([As[2], As[1] * x + As[4], As[0] * x * x + As[3] * x + As[5]])
        for y in ys:
            pt = _newtonOnPair(As, Bs, [x, y], opts['polish'])
            if abs(_evalConicAt(Bs, pt)) <= opts['tol'] * _conicScale(Bs, pt) and \
               abs(_evalConicAt(As, pt)) <= opts['tol'] * _conicScale(As, pt):
                P.append(pt)
    if P:
        P = [[p[0] + k * p[1], p[1]] for p in P]   # undo the shear x = x' + k*y'
        P = _dedupe(P, opts['tol'])
    return P, info


def _shearFor(A, B):
    for k in range(5):
        if (A[0] * k * k + A[1] * k + A[2]) != 0 and (B[0] * k * k + B[1] * k + B[2]) != 0:
            return k
    error('conicMeet:noShear',
          'no integer shear in 0..4 gives both conics a y^2 term, which is impossible unless one '
          'of them has a=b=c=0 and should have taken the line branch.')


def _shearConic(C, k):
    a, b, c, d, e, f = C
    return rat_q.wall([a, 2 * a * k + b, a * k * k + b * k + c, d, d * k + e, f],
                      [abs(a), abs(2 * a * k) + abs(b), abs(a) * k * k + abs(b * k) + abs(c),
                       abs(d), abs(d * k) + abs(e), abs(f)], 'sheared conic')


def _resultantY(A, B):
    """Resultant in y of the two conics, as a quartic in x (5 integer coefficients)."""
    a1, b1, c1 = [A[2]], [A[1], A[4]], [A[0], A[3], A[5]]
    a2, b2, c2 = [B[2]], [B[1], B[4]], [B[0], B[3], B[5]]
    z = [0]
    M = [[a1, b1, c1, z],
         [z, a1, b1, c1],
         [a2, b2, c2, z],
         [z, a2, b2, c2]]
    q = _trimPoly(_polyDet(M))
    return rat_q.wall([0] * (5 - len(q)) + q, None, 'resultant quartic')


def _polyDet(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    d = [0]
    for j in range(n):
        if all(x == 0 for x in M[0][j]):
            continue
        sub = [row[:j] + row[j + 1:] for row in M[1:]]
        term = rat_q.wall(_conv(M[0][j], _polyDet(sub)), None, 'polynomial product')
        if j % 2 == 1:
            term = [-x for x in term]
        d = rat_q.wall(rat_q.polyAdd1(d, term), None, 'polynomial sum')
    return d


def _conv(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def _trimPoly(p):
    i = next((k for k, x in enumerate(p) if x != 0), None)
    if i is None:
        return [0]
    return list(p[i:])


def _realRoots(p):
    """Distinct real roots of a real polynomial (highest first), sorted. Doubles."""
    p = _trimPoly(p)
    if len(p) <= 1:
        return []
    m = max(abs(float(x)) for x in p)
    pf = [float(x) / m for x in p]                 # scale only -- the ROOTS are unchanged
    rr = np.roots(pf)
    re = [float(z.real) for z in rr if abs(z.imag) <= 1e-7 * max(1.0, abs(z.real))]
    re.sort()
    if len(re) > 1:
        keep = [re[0]]
        for i in range(1, len(re)):
            if re[i] - re[i - 1] > 1e-10 * max(1.0, abs(re[i])):
                keep.append(re[i])
        re = keep
    return re


def _realRootsExact(p):
    return _realRoots(_squarefreePart(p))


def _squarefreePart(p):
    p = _trimPoly(p)
    if len(p) <= 2:
        return p
    try:
        g = _polyGcdInt(p, rat_q.polyder(p))
    except CCA2Error as ME:                         # only reachable under EMULATE_2P53
        if ME.identifier != 'ratQ:overflow':
            raise
        return p
    if len(g) <= 1:
        return p
    s = rat_q.polyDivExact1(p, g)
    if not s:
        return p
    return s


def _polyGcdInt(p, q):
    p = _primitivePart(_trimPoly(p))
    q = _primitivePart(_trimPoly(q))
    while len(q) > 1 or (len(q) == 1 and q[0] != 0):
        r = _pseudoRem(p, q)
        p = q
        q = _primitivePart(_trimPoly(r))
        if len(q) == 1 and q[0] == 0:
            break
    return _primitivePart(p)


def _pseudoRem(p, q):
    r = _trimPoly(p)
    q = _trimPoly(q)
    lq = q[0]
    while len(r) >= len(q) and any(x != 0 for x in r):
        shift = len(r) - len(q)
        qq = q + [0] * shift
        rat_q.wall([abs(lq) * abs(a) + abs(r[0]) * abs(b) for a, b in zip(r, qq)], None,
                   'pseudo-remainder')
        r = [lq * a - r[0] * b for a, b in zip(r, qq)]
        r = _trimPoly(r[1:])                        # the leading term cancels by construction
        r = _primitivePart(r)
        if len(r) == 1 and r[0] == 0:
            break
    return r


def _primitivePart(p):
    p = _trimPoly(p)
    g = 0
    for x in p:
        if x != 0:
            g = math.gcd(g, abs(x))
    if g > 1:
        p = [x // g for x in p]
    if p[0] < 0:
        p = [-x for x in p]
    return p


def _rcond2(J):
    a, b, c, d = J[0][0], J[0][1], J[1][0], J[1][1]
    det = a * d - b * c
    if det == 0 or not math.isfinite(det):
        return 0.0
    n1 = max(abs(a) + abs(c), abs(b) + abs(d))
    ninv = max(abs(d) + abs(c), abs(b) + abs(a)) / abs(det)
    if n1 == 0 or ninv == 0:
        return 0.0
    return 1.0 / (n1 * ninv)


def _solve2(J, r):
    """J \\ r for a 2x2 system, by elimination with partial pivoting (as LAPACK does)."""
    a, b, c, d = J[0][0], J[0][1], J[1][0], J[1][1]
    r0, r1 = r
    if abs(c) > abs(a):
        a, b, c, d = c, d, a, b
        r0, r1 = r1, r0
    if a == 0:
        return [math.nan, math.nan]
    m = c / a
    d2 = d - m * b
    y2 = r1 - m * r0
    if d2 == 0:
        return [math.nan, math.nan]
    x2 = y2 / d2
    x1 = (r0 - b * x2) / a
    return [x1, x2]


def _newton(A, B, p, steps, rc):
    for _ in range(steps):
        J = [_gradConic(A, p), _gradConic(B, p)]
        r = [_evalConicAt(A, p), _evalConicAt(B, p)]
        if not all(math.isfinite(v) for row in J for v in row) or _rcond2(J) < rc:
            break
        step = _solve2(J, r)
        if not all(math.isfinite(v) for v in step):
            break
        p = [p[0] - step[0], p[1] - step[1]]
        if math.hypot(step[0], step[1]) <= EPS * max(1.0, math.hypot(p[0], p[1])):
            break
    return p


def _newtonOnPair(A, B, p, steps):
    return _newton(A, B, p, steps, 1e-14)


def _evalConicAt(C, p):
    return (C[0] * p[0] ** 2 + C[1] * p[0] * p[1] + C[2] * p[1] ** 2
            + C[3] * p[0] + C[4] * p[1] + C[5])


def _conicScale(C, p):
    return max(1.0, max(abs(x) for x in C) * max(1.0, max(abs(p[0]), abs(p[1])) ** 2))


def _dedupe(P, tol):
    keep = [True] * len(P)
    for i in range(1, len(P)):
        for j in range(i):
            if keep[j] and math.hypot(P[i][0] - P[j][0], P[i][1] - P[j][1]) <= \
                    tol * max(1.0, math.hypot(P[j][0], P[j][1])):
                keep[i] = False
                break
    return [p for p, k in zip(P, keep) if k]


def _polishAndScore(P, A, B, opts):
    out = []
    resid = []
    for p in P:
        p = _newton(A, B, list(p), opts['polish'], 1e-12)
        out.append(p)
        resid.append(max(abs(_evalConicAt(A, p)) / _conicScale(A, p),
                         abs(_evalConicAt(B, p)) / _conicScale(B, p)))
    return out, resid


def _gradConic(C, p):
    return [2 * C[0] * p[0] + C[1] * p[1] + C[3], C[1] * p[0] + 2 * C[2] * p[1] + C[4]]

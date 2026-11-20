"""cell_has_interior -- does a cell cut out by conic sign conditions have interior?
(port of cellHasInterior.m)

rows is a list of [index (1-based) into EcQ, side]. Lines are decided exactly by
rat_q.feasible2; a curved cell is decided by looking for a witness point (pairwise conic
intersections, then a seeded random probe) -- the same screen MATLAB uses -- and, when that finds
nothing, by an exact and complete CAD (_cadPoint; the completeness argument is its docstring).
"""
from fractions import Fraction
from functools import cmp_to_key
import math

import numpy as np

from . import rat_q
from .conic_meet import conicMeet
from .errors import CCA2Error


def cellHasInterior(EcQ, rows):
    if not rows:
        return True
    isLine = [all(x == 0 for x in EcQ[r[0] - 1][0:3]) for r in rows]
    if all(isLine):
        P = [[r[1] * x for x in EcQ[r[0] - 1][3:6]] for r in rows]
        return rat_q.feasible2(P, True)
    if any(isLine):
        P = [[r[1] * x for x in EcQ[r[0] - 1][3:6]] for r, l in zip(rows, isLine) if l]
        if not rat_q.feasible2(P, True):
            return False
    for r, l in zip(rows, isLine):
        if r[1] != 0 and not l and _negSemidefinite([r[1] * x for x in EcQ[r[0] - 1][0:6]]):
            return False
    idx = sorted(set(r[0] for r in rows))
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            try:
                P, info = conicMeet(EcQ[idx[a] - 1], EcQ[idx[b] - 1])
            except CCA2Error as ME:                 # only reachable under rat_q.EMULATE_2P53:
                if ME.identifier != 'ratQ:overflow':    # MATLAB keeps a cell it cannot decide
                    raise
                return True
            if info['degenerate'] or not P:
                continue
            for p in P:
                if _satisfiesAll(EcQ, rows, p):
                    return True
    if _probeInterior(EcQ, rows):
        return True
    return _cadInterior(EcQ, rows)


def _cadInterior(EcQ, rows):
    return _cadPoint(EcQ, rows) is not None


def _cadPoint(EcQ, rows):
    """EXACT and COMPLETE decision of "does the open cell {side_i * q_i > 0} have a point?",
    run when the witness search finds nothing.  (Side-0 rows -- equalities -- have no open
    cell; the function answers False for them, as the witness stage alone decides those.)

    COMPLETENESS ARGUMENT (a cylindrical algebraic decomposition over Q, no floating point).
    1. Shear x -> x + k y (k a small integer, chosen so that every non-line conic has a nonzero
       y^2 coefficient).  A shear is a linear bijection, so the cell is non-empty iff the
       sheared cell is.  Afterwards every row is of one of four kinds: a *full conic*
       (a2 != 0), a *slanted line* (y-degree 1), a *vertical line* a3 x + a5 (y-degree 0), or
       a constant.
    2. The CRITICAL ABSCISSAE X are the real roots (a finite set, as every polynomial below is
       non-zero) of
         (a) disc_y(q) of every full conic q  (roots appear/merge: also the crossing of the two
             lines of a degenerate conic, which is its singular point);
         (b) a3 x + a5 of every vertical line  (a vertical line is invisible to every resultant,
             which was defect (a) of finding #2);
         (c) Res_y(q_i, q_j) of every pair with a full conic and a non-zero resultant, and the
             linear resultant of two slanted lines;
         (d) for a pair of non-proportional full conics with Res_y == 0 (they share a line L:
             q_i = L M_i, q_j = L M_j; the gcd of two rational conics is rational, since a
             shared irrational or complex line would bring its conjugate and make them
             proportional): the abscissa of M_i ^ M_j.  The other pair crossings L ^ M are
             singular points of q_i or q_j, covered by (a).  This is defect (b) of #2: the
             resultant is identically zero, so the old code dropped the pair.
       A slanted line sharing a factor with a conic is a factor of it, and adds nothing beyond
       (a).
    3. On each open interval I of R \\ X the number of distinct real roots of each q_i(x, .),
       and the order of all the roots of all the q_i(x, .) together, is constant: two root
       branches can only meet where a discriminant (a), a leading coefficient (constant after
       the shear, or the vertical lines (b)), or a resultant (c)/(d) vanishes.  Hence the sign
       of every row on the k-th open sector of the vertical line {x} x R (between the k-th and
       (k+1)-th root; below the first; above the last) is the same for all x in I.
    4. If the cell S is non-empty it is open, so its projection to x contains an open
       interval and therefore a point x0 outside X, and S meets the line {x0} x R in a
       non-empty open set, which contains a point off every root, i.e. a point of some sector.
       By 3 the same sector of any other column x_s of the same interval I is in S.  So testing
       one rational point per sector of one rational column per interval I (and beyond both
       ends) decides the cell.  The rational columns and sector points are taken from exact
       Sturm isolation (rat_q.isolateRoots); the sign tests are exact Fractions.
    _cadPoint returns that rational point (a certificate) or None (a proof of emptiness).  The random probe box of the witness stage is not needed for either."""
    if any(r[1] == 0 for r in rows):
        return None
    idx = sorted(set(r[0] for r in rows))
    cs = [[int(t) for t in EcQ[i - 1][0:6]] for i in idx]
    sides = {}
    for r in rows:
        sides.setdefault(r[0], set()).add(r[1])
    if any(len(v) > 1 for v in sides.values()):
        return None                                     # q > 0 and q < 0 on the same conic
    sg = [next(iter(sides[i])) for i in idx]
    k = 0
    while any((c[0] != 0 or c[1] != 0 or c[2] != 0) and c[0] * k * k + c[1] * k + c[2] == 0
              for c in cs):
        k += 1
    cs = [_shear(c, k) for c in cs]
    X, special = _criticalPolys(cs)
    for x0 in _samples(X):
        yp = _yPolys(cs, sg, x0, special)
        if yp is None:
            continue
        for y0 in _samples(yp):
            if all(s * _evalExact(c, x0, y0) > 0 for c, s in zip(cs, sg)):
                return (x0 + k * y0, y0)                # a rational point of the cell (unsheared)
    return None


def _shear(c, k):
    a, b, cc, d, e, f = c
    return [a, 2 * a * k + b, a * k * k + b * k + cc, d, d * k + e, f]


def _evalExact(c, x, y):
    return c[0] * x * x + c[1] * x * y + c[2] * y * y + c[3] * x + c[4] * y + c[5]


def _pmul(a, b):
    return rat_q.polyMul1(a, b)


def _padd(a, b, sb=1):
    return rat_q.polyTrim(rat_q.polyAdd1(a, [sb * x for x in b]))


def _criticalPolys(cs):
    """The integer polynomials in x whose real roots are the critical abscissae (step 2), and the
    set of row indices whose y-polynomials can share a root at EVERY column (a repeated factor, or
    a common one with another row): at any other column the y-polynomials are squarefree and
    pairwise coprime, because a resultant that is not identically zero has no root there."""
    polys = []
    special = set()
    full = [i for i, c in enumerate(cs) if c[2] != 0]
    slant = [i for i, c in enumerate(cs) if c[2] == 0 and c[0] == 0 and c[1] == 0 and c[4] != 0]
    for i, (a0, a1, a2, a3, a4, a5) in enumerate(cs):
        if a2 != 0:
            d = [a1 * a1 - 4 * a2 * a0, 2 * a1 * a4 - 4 * a2 * a3, a4 * a4 - 4 * a2 * a5]
            if not rat_q.polyTrim(d):
                special.add(i)
            polys.append(d)
        elif a4 == 0 and a2 == 0 and a1 == 0 and a0 == 0:
            polys.append([a3, a5])                      # vertical line (or a constant)
    for ii, i in enumerate(full):
        p = cs[i]
        al, be, ga = [p[2]], [p[1], p[4]], [p[0], p[3], p[5]]
        for j in full[ii + 1:]:
            q = cs[j]
            al2, be2, ga2 = [q[2]], [q[1], q[4]], [q[0], q[3], q[5]]
            u = _padd(_pmul(al2, be), _pmul(al, be2), -1)      # alpha' beta  - alpha beta'
            w = _padd(_pmul(al2, ga), _pmul(al, ga2), -1)      # alpha' gamma - alpha gamma'
            v = _padd(_pmul(be, ga2), _pmul(be2, ga), -1)      # beta gamma'  - beta' gamma
            res = _padd(_pmul(w, w), _pmul(u, v))              # = (a g' - a' g)^2 - (a b' - a' b)(b g' - b' g)
            if res:
                polys.append(res)
            else:
                special.update((i, j))
                if u or w:
                    polys.extend(_sharedLineCrossing(p, q, u, w))
        for j in slant:
            q = cs[j]
            B, C = [q[4]], [q[3], q[5]]
            res = _padd(_padd(_pmul(al, _pmul(C, C)), _pmul(be, _pmul(B, C)), -1),
                        _pmul(ga, _pmul(B, B)))
            if res:
                polys.append(res)
            else:
                special.update((i, j))
    for ii, i in enumerate(slant):
        p = cs[i]
        for j in slant[ii + 1:]:
            q = cs[j]
            res = _padd(_pmul([p[4]], [q[3], q[5]]), _pmul([q[4]], [p[3], p[5]]), -1)
            if res:
                polys.append(res)
            else:
                special.update((i, j))
    return polys, special


def _sharedLineCrossing(p, q, u, w):
    """p = L M_p and q = L M_q share the line L (resultant == 0, p, q not proportional).
    T = alpha' p - alpha q = u(x) y + w(x) = L * N with N vertical, so L = T / gcd(u, w).
    Returns [den*x - num] for the abscissa of M_p ^ M_q (nothing if they are parallel)."""
    g = rat_q.polyGcd(u, w) if w else rat_q.polyPrimitive(u)
    if not u:
        error('cellHasInterior:internal', 'shared line factor with no y term: %s %s' % (p, q))
    lu, lw = _fdiv(u, g), _fdiv(w, g)
    # L = ell y + pp x + rr
    ell = lu[-1]
    pp = lw[-2] if len(lw) >= 2 else Fraction(0)
    rr = lw[-1] if lw else Fraction(0)
    if len(lu) != 1 or len(lw) > 2 or ell == 0:
        error('cellHasInterior:internal', 'shared factor of two conics is not a line: %s %s' % (p, q))
    M = []
    for c in (p, q):
        a0, a1, a2, a3, a4, a5 = [Fraction(t) for t in c]
        m = a2 / ell
        kk = (a1 - pp * m) / ell
        hh = (a4 - rr * m) / ell
        M.append((m, kk, hh))
    (m1, k1, h1), (m2, k2, h2) = M
    det = m1 * k2 - m2 * k1
    if det == 0:
        return []
    xp = (m2 * h1 - m1 * h2) / det
    return [[xp.denominator, -xp.numerator]]


def _fdiv(p, g):
    """p / g as a list of Fractions (exact division of integer polynomials over Q); the
    quotient is empty when p is zero."""
    p = [Fraction(x) for x in rat_q.polyTrim(p)]
    g = [Fraction(x) for x in rat_q.polyTrim(g)]
    if not p:
        return []
    out = []
    for i in range(len(p) - len(g) + 1):
        c = p[i] / g[0]
        out.append(c)
        for j, gj in enumerate(g):
            p[i + j] -= c * gj
    if len(p) < len(g) or any(x != 0 for x in p):
        error('cellHasInterior:internal', 'inexact polynomial division in a shared line factor')
    return out


def _yPolys(cs, sg, x0, special):
    """The y-polynomials q(x0, .) of one column, scaled to integers, as (poly, special) pairs for
    _samples -- or None when the column holds no point of the cell.  Exact shortcuts: a row that is
    constant in y with the wrong sign, or a quadratic without real roots with the wrong sign,
    kills the column; a quadratic without real roots with the right sign holds everywhere, so its
    (absent) roots are not sample boundaries.  The rows in `special` need the coprime-base
    treatment in _samples; the others are squarefree and coprime to everything at this
    (non-critical) column."""
    out = []
    for i, (a0, a1, a2, a3, a4, a5) in enumerate(cs):
        co = [Fraction(a2), a1 * x0 + a4, a0 * x0 * x0 + a3 * x0 + a5]
        den = 1
        for t in co:
            den = den * t.denominator // math.gcd(den, t.denominator)
        co = rat_q.polyPrimitive([int(t * den) for t in co])
        if not co:
            return None
        if len(co) == 1:
            if sg[i] * co[0] <= 0:
                return None
        elif len(co) == 3 and co[1] * co[1] - 4 * co[0] * co[2] < 0:
            if sg[i] * co[0] < 0:
                return None
        else:
            out.append((co, i in special))
    return out


def _coprimeBase(polys):
    """Pairwise coprime squarefree primitive polynomials with the same set of real roots as the
    product of `polys` (so no two of them share a root)."""
    work = [rat_q.polyPrimitive(rat_q.polySquarefree(p)) for p in polys if len(rat_q.polyTrim(p)) > 1]
    base = []
    while work:
        p = rat_q.polyPrimitive(rat_q.polySquarefree(work.pop()))
        if len(p) <= 1:
            continue
        for i, b in enumerate(base):
            g = rat_q.polyGcd(p, b)
            if len(g) > 1:
                base.pop(i)
                work += [rat_q.polyDivExact1(b, g), rat_q.polyDivExact1(p, g), g]
                break
        else:
            base.append(p)
    return base


class _Root:
    """The unique real root of `poly` (squarefree) in the half-open interval (lo, hi]."""

    def __init__(self, poly, lo, hi, chain):
        self.poly, self.lo, self.hi, self.chain = poly, lo, hi, chain

    def refine(self):
        mid = (self.lo + self.hi) / 2
        if rat_q.sturmV(self.chain, mid) - rat_q.sturmV(self.chain, self.hi) == 1:
            self.lo = mid
        else:
            self.hi = mid


def _before(r, s):
    """r < s for roots of coprime polynomials (so r != s): refine until the intervals separate."""
    while True:
        if r.hi <= s.lo:
            return True
        if s.hi <= r.lo:
            return False
        (r if r.hi - r.lo >= s.hi - s.lo else s).refine()


def _samples(polys):
    """Rational sample points: one in every open interval between consecutive distinct real
    roots of the polynomials in `polys` (integer polynomials), one below and one above all
    roots (0 when there is none).  Exact: coprime base, Sturm isolation of each factor, pairwise
    separation by bisection."""
    if polys and isinstance(polys[0], tuple):
        polys = _coprimeBase([p for p, sp in polys if sp]) + [p for p, sp in polys if not sp]
    else:
        polys = _coprimeBase(polys)
    roots = []
    for p in polys:
        chain = rat_q.sturmChain(p)
        for lo, hi in rat_q.isolateRoots(p):
            roots.append(_Root(p, lo, hi, chain))
    if not roots:
        return [Fraction(0)]
    roots.sort(key=cmp_to_key(lambda r, s: -1 if _before(r, s) else 1))
    out = [roots[0].lo]                                 # r_1 > lo_1 and r_1 is the least root
    for r, s in zip(roots, roots[1:]):
        while r.hi >= s.lo:                             # r < s: bisect until a gap opens
            (r if r.hi - r.lo >= s.hi - s.lo else s).refine()
        out.append((r.hi + s.lo) / 2)                   # r <= hi < mid < lo < s
    out.append(roots[-1].hi + 1)                        # every root is <= its hi
    return out


def _negSemidefinite(c):
    M = [[-2 * c[0], -c[1], -c[3]], [-c[1], -2 * c[2], -c[4]], [-c[3], -c[4], -2 * c[5]]]
    flat = [x for r in M for x in r]
    if max(abs(x) for x in flat) > 2 ** 16 or not any(flat):
        return False
    if any(M[i][i] < 0 for i in range(3)):
        return False
    for p in ((0, 1), (0, 2), (1, 2)):
        if M[p[0]][p[0]] * M[p[1]][p[1]] - M[p[0]][p[1]] ** 2 < 0:
            return False
    d = (M[0][0] * M[1][1] * M[2][2] + 2 * M[0][1] * M[1][2] * M[0][2]
         - M[0][0] * M[1][2] ** 2 - M[1][1] * M[0][2] ** 2 - M[2][2] * M[0][1] ** 2)
    return d >= 0


def satisfiesAllDouble(EcQ, rows, x):
    """Do all the cell's sign conditions hold at the double point x, to a 1e-9 relative band?"""
    for r in rows:
        c = EcQ[r[0] - 1]
        v = (c[0] * x[0] ** 2 + c[1] * x[0] * x[1] + c[2] * x[1] ** 2
             + c[3] * x[0] + c[4] * x[1] + c[5])
        scale = max(1.0, max(abs(t) for t in c) * max(1.0, max(abs(x[0]), abs(x[1])) ** 2))
        if r[1] == 0:
            if abs(v) > 1e-9 * scale:
                return False
        elif r[1] * v < -1e-9 * scale:
            return False
    return True


_satisfiesAll = satisfiesAllDouble


def matlabRand(seed, n, m):
    """MATLAB ``rng(seed); rand(n, m)`` -- the same Mersenne Twister stream, column-major."""
    rs = np.random.RandomState(seed)
    u = rs.random_sample(n * m)
    return [[float(u[j * n + i]) for j in range(m)] for i in range(n)]


_PROBE = None


def _probeInterior(EcQ, rows):
    global _PROBE
    idx = sorted(set(r[0] for r in rows))
    pts = []
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            try:
                P, _ = conicMeet(EcQ[idx[a] - 1], EcQ[idx[b] - 1])
            except CCA2Error as ME:
                if ME.identifier != 'ratQ:overflow':
                    raise
                continue
            pts.extend(P)
    if not pts:
        ctr = [0.0, 0.0]
        r = 4.0
    else:
        ctr = [sum(p[0] for p in pts) / len(pts), sum(p[1] for p in pts) / len(pts)]
        r = max(2.0, 2 * max(max(abs(p[0] - ctr[0]), abs(p[1] - ctr[1])) for p in pts))
    if _PROBE is None:
        _PROBE = np.array(matlabRand(20260905, 400, 2), dtype=float)
    X0 = ctr[0] + (_PROBE[:, 0] - 0.5) * 2 * r
    X1 = ctr[1] + (_PROBE[:, 1] - 0.5) * 2 * r
    return bool(_satisfiesAllMany(EcQ, rows, X0, X1).any())


def _satisfiesAllMany(EcQ, rows, X0, X1):
    """satisfiesAllDouble at every point (X0[i], X1[i]) at once: the same double operations in the
    same order, elementwise, so each verdict is the one the scalar loop gives."""
    ok = np.ones(X0.shape, dtype=bool)
    m = np.maximum(1.0, np.maximum(np.abs(X0), np.abs(X1))) ** 2
    for r in rows:
        c = EcQ[r[0] - 1]
        f = [float(t) for t in c[0:6]]
        v = f[0] * X0 ** 2 + f[1] * X0 * X1 + f[2] * X1 ** 2 + f[3] * X0 + f[4] * X1 + f[5]
        scale = np.maximum(1.0, float(max(abs(t) for t in c)) * m)
        if r[1] == 0:
            ok &= ~(np.abs(v) > 1e-9 * scale)
        else:
            ok &= ~(r[1] * v < -1e-9 * scale)
        if not ok.any():
            break
    return ok

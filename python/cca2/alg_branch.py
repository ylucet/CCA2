"""alg_branch -- evaluate a NAMED real root of an implicit face polynomial (port of algBranch.m).

P(x, y, z) is an exact polynomial (a poly_n dict in three variables); the value at (x, y) is the
b-th real root in z of the fibre -- positive b counts up from the smallest real root, negative
b down from the largest (-1 is the top). Doubles: this is the EVALUATION of an AlgCrv face, the
one place a tolerance is allowed. Returns z (list, nan where the branch does not exist) or, with
wantInfo, (z, info) with info['allReal'] the sorted real roots per point.

Fidelity note: MATLAB evaluates the fibre in doubles when every coefficient fits a double and in
BigDecimal when the polynomial is stored in BigInteger; the port switches on the same condition
(max |coefficient| >= 2^53) so that the degenerate-leader tolerance test is the same one.
"""
import math
from fractions import Fraction

from .errors import error
from .numeric import roots

REALMIN = 2.2250738585072014e-308


def _layers(P):
    """[(k, [(i, j, c), ...])] by z-power, and the sizes nx, ny, nz."""
    nx = max(e[0] for e in P) + 1
    ny = max(e[1] for e in P) + 1
    nz = max(e[2] for e in P) + 1
    lay = [[] for _ in range(nz)]
    for e, c in P.items():
        lay[e[2]].append((e[0], e[1], c))
    return lay, nx, ny, nz


def _fibreDouble(lay, x, y, degTol):
    c = []
    isZero = []
    ca = []
    for terms in lay:
        v = 0.0
        a = 0.0
        for i, j, co in terms:
            t = float(co) * x ** i * y ** j
            v += t
            a += abs(t)
        c.append(v)
        ca.append(a)
        isZero.append(abs(v) <= degTol * a)
    return c, isZero, ca


def _fibreExact(lay, x, y):
    X = Fraction(x)
    Y = Fraction(y)
    c = []
    isZero = []
    for terms in lay:
        v = Fraction(0)
        for i, j, co in terms:
            v += co * X ** i * Y ** j
        c.append(float(v))
        isZero.append(v == 0)
    return c, isZero


def _absFibre(lay, x, y):
    return [sum(abs(float(co)) * abs(x) ** i * abs(y) ** j for i, j, co in terms)
            for terms in lay]


def _exactQuad(lay, x, y):
    X = Fraction(x)
    Y = Fraction(y)
    acc = []
    for k in range(3):
        v = Fraction(0)
        if k < len(lay):
            for i, j, co in lay[k]:
                v += co * X ** i * Y ** j
        acc.append(v)
    Dx = acc[1] * acc[1] - 4 * acc[2] * acc[0]
    c = [float(a) for a in acc]
    D = float(Dx)
    if D == 0 and Dx != 0:
        D = REALMIN if Dx > 0 else -REALMIN
    return c, D


def _doubleRootNoise(cz, v):
    """Is the complex root v a REAL DOUBLE root split by rounding? Rounding the fibre's
    coefficients by eps moves a double root t by about sqrt(eps * sum|c_k t^k| / |p''(t)/2|), in
    either direction: along the real axis (two close reals) or off it (a conjugate pair), as the
    LAPACK in use happens to round. A pair whose imaginary part is within that size is counted as
    the double root it is, so the branch numbering does not depend on the library."""
    t = v.real
    cabs = sum(abs(c) * abs(t) ** k for k, c in enumerate(cz))
    p2 = sum(k * (k - 1) / 2 * c * t ** (k - 2) for k, c in enumerate(cz) if k >= 2)
    if p2 == 0 or cabs == 0:
        return False
    return abs(v.imag) <= 8.0 * math.sqrt(2.220446049250313e-16 * cabs / abs(p2))


def _sgnz(v):
    """sign(v + (v == 0)): +1 at zero."""
    return 1.0 if v >= 0 else -1.0


def algBranch(P, X, b, opts=None, wantInfo=False):
    opts = dict(opts or {})
    imagTol = opts.get('imagTol', 1e-10)
    degTol = opts.get('degTol', 1e-12)
    if not P:
        error('algBranch:polynomial', 'P must be a non-empty polynomial.')
    if not isinstance(b, int) or isinstance(b, bool) or b == 0:
        error('algBranch:branch',
              'the branch index must be a nonzero integer -- positive counts up from the '
              'smallest real root, negative down from the largest; got %r.' % (b,))
    lay, nx, ny, nz = _layers(P)
    if nz < 2:
        error('algBranch:noZ',
              'P has degree 0 in z, so {p = 0} does not define a value: it is a condition on '
              '(x,y) alone. A face stored this way is a producer defect, not an evaluation one.')
    big = max(abs(c) for c in P.values()) >= 2 ** 53

    if not wantInfo and nz <= 3:
        return [_fibreLow(lay, float(p[0]), float(p[1]), b, imagTol, degTol, nz, big) for p in X]

    z = []
    allReal = []
    nreal = []
    for p in X:
        x = float(p[0])
        y = float(p[1])
        if big:
            c, isZero = _fibreExact(lay, x, y)
        else:
            c, isZero, _ = _fibreDouble(lay, x, y, degTol)
        zt = math.nan
        rr = []
        scale = max(abs(v) for v in c)
        keep = max((k for k in range(nz) if not isZero[k]), default=-1)
        cz = c[:keep + 1]
        if scale != 0 and len(cz) >= 2:
            r = roots(cz[::-1])                     # `roots` wants highest power first
            rr = sorted(v.real for v in r
                        if abs(v.imag) <= imagTol * max(1.0, abs(v)) or _doubleRootNoise(cz, v))
            if len(cz) == 3:
                ca = _absFibre(lay[:3], x, y)
                dq = cz[1] ** 2 - 4 * cz[2] * cz[0]
                T = ca[1] ** 2 + 4 * ca[2] * ca[0]
                if abs(dq) <= 1e-12 * T:
                    cx, dx = _exactQuad(lay, x, y)
                    if dx > 0:
                        qq = -(cx[1] + _sgnz(cx[1]) * math.sqrt(dx)) / 2
                        rr = sorted([qq / cx[2], cx[0] / qq])
                    else:
                        rr = [-cz[1] / (2 * cz[2])] * 2
            idx = b if b > 0 else len(rr) + 1 + b   # -1 is the largest, -2 the next below
            if 1 <= idx <= len(rr):
                zt = rr[idx - 1]
        z.append(zt)
        allReal.append(rr)
        nreal.append(len(rr))
    if wantInfo:
        return z, {'nreal': nreal, 'allReal': allReal}
    return z


def _fibreLow(lay, x, y, b, imagTol, degTol, nz, big):
    """The fibre is at most quadratic in z: closed form, one point."""
    ca = _absFibre(lay, x, y) + [0.0] * (3 - nz)
    if big:
        c, isZero = _fibreExact(lay, x, y)
    else:
        c, isZero, _ = _fibreDouble(lay, x, y, degTol)
    c = c + [0.0] * (3 - nz)
    isZero = isZero + [True] * (3 - nz)
    lead = 3 if not isZero[2] else (2 if not isZero[1] else 0)
    if lead == 2:
        if max(abs(v) for v in c) > 0 and b in (1, -1):
            return -c[0] / c[1]
        return math.nan
    if lead != 3:
        return math.nan
    a, bb, c0 = c[2], c[1], c[0]
    dq = bb * bb - 4 * a * c0
    T = ca[1] ** 2 + 4 * ca[2] * ca[0]
    win = abs(dq) <= 1e-12 * T
    if win:
        cx, dx = _exactQuad(lay, x, y)
        if dx > 0:
            a, bb, c0 = cx[2], cx[1], cx[0]
            dq = dx
        else:
            dq = 0.0
    mean2 = -bb / (2 * a)
    lo = hi = math.nan
    if win and dq <= 0:
        lo = hi = mean2
    elif dq > 0:
        qq = -(bb + _sgnz(bb) * math.sqrt(dq)) / 2
        r1 = qq / a
        r2 = c0 / qq if qq != 0 else math.nan
        lo = min(r1, r2)
        hi = max(r1, r2)
    elif dq < 0:                                    # a complex pair: real only within imagTol
        im = math.sqrt(-dq) / (2 * abs(a))
        if im <= imagTol * max(1.0, math.hypot(mean2, im)):
            lo = hi = mean2
    if b in (1, -2):
        return lo
    if b in (2, -1):
        return hi
    return math.nan

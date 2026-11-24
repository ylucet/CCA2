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
    """The fibre's coefficients at the double point (x, y), each computed EXACTLY and then
    rounded once: x = a/dx and y = c/dy exactly (doubles are dyadic), so with
    D = dx^(nx-1) dy^(ny-1) every layer's value times D is the integer
    sum co a^i dx^(nx-1-i) c^j dy^(ny-1-j). float(num / D) is Python's correctly rounded integer
    division -- the same double float(Fraction) gives -- and the value is zero iff num is."""
    a, dx = float(x).as_integer_ratio()
    cc, dy = float(y).as_integer_ratio()
    nx = 1 + max((i for terms in lay for i, _, _ in terms), default=0)
    ny = 1 + max((j for terms in lay for _, j, _ in terms), default=0)
    A = [a ** i * dx ** (nx - 1 - i) for i in range(nx)]
    C = [cc ** j * dy ** (ny - 1 - j) for j in range(ny)]
    D = dx ** (nx - 1) * dy ** (ny - 1)
    c = []
    isZero = []
    for terms in lay:
        v = 0
        for i, j, co in terms:
            v += co * A[i] * C[j]
        c.append(v / D)
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


_PROBES = ((1, 0), (0, 1), (1, 1), (1, -1))


def _nearDoubles(lay, x, y, cz, r, rr, imagTol):
    """Decide EXACTLY the near-double roots the doubles cannot (#108 residual, D-153).

    A NEAR-DOUBLE is a pair of roots within the double-root noise of each other: two real roots
    _doubleRootNoise could not tell apart, or a complex pair it accepted as one real double root.
    The doubles cannot say whether such a pair is a value of the face, and the numbering above
    every one of them depends on it.

    A face is the b-th real root of its fibre, and f** is closed and convex: its value at a point
    is the limit from the points around it. So a pair that is real only on a CURVE -- a complex
    pair touching the real axis there and leaving it on both sides -- is not the value of any
    branch, and must not be counted. Such a pair appears where a bridge's dual point runs off to
    infinity, on a boundary line of dom f**. For a continuous f it coincides with the true value
    there; for a discontinuous one it does not (D-153: the convex piece's own value 10.2578 on its
    hull edge next to the jump vertex, f** = 9.5177), and the noise rule also accepted it in a
    ~1e-5 strip inside the domain.

    The test: the fibre is evaluated EXACTLY (doubles are dyadic) at x +- h e for four directions
    e and a tiny dyadic h, and its distinct real roots near the pair are counted by Sturm. A pair
    with no real root on EITHER side in at least two directions is dropped. A crossing of two real
    branches, a fold (real on one side) and a double root that is real along a curve all keep it:
    at most one of the four directions can run along that curve."""
    cands = []
    for v in r:
        if v.imag > 0 and abs(v.imag) > imagTol * max(1.0, abs(v)) and _doubleRootNoise(cz, v):
            cands.append(v.real)
    for a, b in zip(rr[:-1], rr[1:]):
        if b > a and _doubleRootNoise(cz, complex((a + b) / 2, (b - a) / 2)):
            cands.append((a + b) / 2)
    if not cands:
        return rr
    out = list(rr)
    for z0 in cands:
        d = sorted(abs(complex(v) - z0) for v in r)[2:]       # the pair itself is nearest
        rad = d[0] / 2 if d and d[0] > 0 else 1.0 + abs(z0)
        if _isolatedPair(lay, x, y, z0, rad):
            for _ in range(2):
                k = min(range(len(out)), key=lambda i: abs(out[i] - z0), default=None)
                if k is not None and abs(out[k] - z0) <= rad:
                    out.pop(k)
    return _snapMultipleRoots(lay, x, y, cz, out)


def _snapMultipleRoots(lay, x, y, cz, rr):
    """A GENUINE multiple root of the fibre, to full precision (#108 residual: the jump vertices).

    Rounding returns a k-fold real root as k roots up to ~eps^(1/k) apart, and at a jump vertex
    f** is such a root: D-153's vertex (-3,-2), f** = 7.5 a double root, came back as 7.5000017
    (2.3e-7 relative), its tangency point (0,3), a triple root 12.5, as 12.500022, and the kink's
    jump vertex (1,2) as 1.0000000256. The multiple roots are the roots of gcd(p, p'), where they
    are SIMPLE: the fibre at the exact dyadic point is split by multiplicity (Yun), the real roots
    of each factor of multiplicity k >= 2 are isolated by Sturm and bisected exactly to a double,
    and the k entries nearest each one, within the rounding of a k-fold root
    8 (eps sum|c_j m^j| / |p^(k)(m)/k!|)^(1/k), take its value."""
    from .rat_q import isolateRoots, polyPrimitive
    if not rr:
        return rr
    p = polyPrimitive(_exactFibreAt(lay, Fraction(x), Fraction(y)))
    if len(p) <= 2:
        return rr
    out = list(rr)
    for k, a in _yun(p):
        if k < 2:
            continue
        for lo, hi in isolateRoots(a):
            m = _bisectRoot(a, lo, hi)
            w = _multipleRootNoise(cz, m, k)
            keep = [v for v in out if abs(v - m) > w]
            if len(keep) < len(out):            # rounding left something there: k copies of m
                out = keep + [m] * k
    return sorted(out)


def _multipleRootNoise(cz, m, k):
    """How far rounding the fibre's coefficients (z^0 first) by eps moves a k-fold root m:
    8 (eps sum|c_j m^j| / |p^(k)(m)/k!|)^(1/k) -- _doubleRootNoise's bound for any k."""
    cabs = sum(abs(c) * abs(m) ** j for j, c in enumerate(cz))
    pk = sum(math.comb(j, k) * c * m ** (j - k) for j, c in enumerate(cz) if j >= k)
    if pk == 0 or cabs == 0:
        return 0.0
    return 8.0 * (2.220446049250313e-16 * cabs / abs(pk)) ** (1.0 / k)


def _yun(p):
    """Yun's square-free decomposition of the integer polynomial p (highest first): the list of
    (k, a_k) with p = c * prod a_k^k, each a_k square-free, primitive and of positive degree."""
    from .rat_q import polyder, polyGcd, polyDivExact1, polyPrimitive
    g = polyGcd(p, polyder(p))
    if len(g) <= 1:
        return [(1, p)]
    out = []
    b = polyDivExact1(polyPrimitive(p), g)         # every distinct root once
    k = 1
    rest = g
    while len(b) > 1:
        a = polyGcd(b, rest)                        # the roots of multiplicity > k
        f = polyDivExact1(b, a) if len(a) > 1 else b
        if len(f) > 1:
            out.append((k, polyPrimitive(f)))       # multiplicity exactly k
        if len(a) <= 1:
            break
        rest = polyDivExact1(rest, a)
        b = a
        k += 1
    return out


def _bisectRoot(a, lo, hi):
    """The simple root of the square-free a in (lo, hi], bisected EXACTLY until the interval is
    below 2^-60 of its size (or 400 halvings); the double nearest the midpoint."""
    from .rat_q import countRootsIn
    lo = Fraction(lo)
    hi = Fraction(hi)
    for _ in range(400):
        if hi - lo <= abs(hi + lo) * Fraction(1, 2 ** 60):
            break
        mid = (lo + hi) / 2
        if countRootsIn(a, lo, mid) == 1:
            hi = mid
        else:
            lo = mid
    return float((lo + hi) / 2)


def _exactFibreAt(lay, X, Y):
    """The fibre at the exact rational point (X, Y), as an integer polynomial in z, highest
    power first (scaled by a positive common denominator: the same roots)."""
    c = []
    for terms in lay:
        v = Fraction(0)
        for i, j, co in terms:
            v += co * X ** i * Y ** j
        c.append(v)
    den = 1
    for v in c:
        den = den * v.denominator // math.gcd(den, v.denominator)
    return [int(v * den) for v in reversed(c)]


def _isolatedPair(lay, x, y, z0, rad):
    from .rat_q import countRootsIn
    X = Fraction(x)
    Y = Fraction(y)
    h = Fraction(2) ** (math.frexp(max(1.0, abs(x), abs(y)))[1] - 40)   # ~1e-12 relative
    lo = Fraction(z0 - rad)
    hi = Fraction(z0 + rad)
    both = 0
    for e in _PROBES:
        if all(countRootsIn(_exactFibreAt(lay, X + sg * h * e[0], Y + sg * h * e[1]), lo, hi) == 0
               for sg in (1, -1)):
            both += 1
    return both >= 2


def fibreVanishes(P, x):
    """Is p(x, y, .) identically zero at the point -- a BASE POINT, where every z is a root and a
    row comparing this value is undefined (#108 residual: D-153's corner (0,0))? Every coefficient
    zero by algBranch's own zero test (exact for a big P, within degTol of its terms otherwise)."""
    lay, nx, ny, nz = _layers(P)
    if max(abs(c) for c in P.values()) >= 2 ** 53:
        _, isZero = _fibreExact(lay, float(x[0]), float(x[1]))
    else:
        _, isZero, _ = _fibreDouble(lay, float(x[0]), float(x[1]), 1e-12)
    return all(isZero)


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
        layF = [[(i, j, float(co)) for i, j, co in terms] for terms in lay]
        return [_fibreLow(lay, float(p[0]), float(p[1]), b, imagTol, degTol, nz, big, layF, nx, ny)
                for p in X]

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
            rr = _nearDoubles(lay, x, y, cz, r, rr, imagTol)
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


def _fibreLow(lay, x, y, b, imagTol, degTol, nz, big, layF, nx, ny):
    """The fibre is at most quadratic in z: closed form, one point.

    The powers x**i, y**j (and of |x|, |y|) are computed once per point and each term is the same
    product, in the same order, as _fibreDouble / _absFibre form it: bit-identical, fewer pow calls."""
    xp = [x ** i for i in range(nx)]
    yp = [y ** j for j in range(ny)]
    ax = abs(x)
    ay = abs(y)
    axp = [ax ** i for i in range(nx)]
    ayp = [ay ** j for j in range(ny)]
    ca = [sum(abs(co) * axp[i] * ayp[j] for i, j, co in terms) for terms in layF] + [0.0] * (3 - nz)
    if big:
        c, isZero = _fibreExact(lay, x, y)
    else:
        c = []
        isZero = []
        for terms in layF:
            v = 0.0
            a = 0.0
            for i, j, co in terms:
                t = co * xp[i] * yp[j]
                v += t
                a += abs(t)
            c.append(v)
            isZero.append(abs(v) <= degTol * a)
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

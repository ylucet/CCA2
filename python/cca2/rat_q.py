"""rat_q -- EXACT rational arithmetic on coefficient vectors (port of ratQ.m).

Ported from ratQ.m. The MATLAB header holds the full rationale; the load-bearing points are:

  * A FACE FUNCTION has a scale: integer numerator vector over ONE positive integer denominator,
    gcd removed -- ``canon``.
  * A CONIC has NO scale: a PRIMITIVE integer vector whose first nonzero entry is positive --
    ``conic``. Two conics are the same curve exactly when the two vectors are equal.
  * Nothing here is approximate. ``fromDouble`` refuses what it cannot represent exactly.

PORT DIVERGENCE (deliberate, PORT_NOTES.md "int64/2^53 wall"): MATLAB stores integers in doubles
and raises ``ratQ:overflow`` at 2^53, with a second BigInteger path behind most routines. Python
``int`` is unbounded, so there is ONE path, ``chk``/``chkB`` no longer bound magnitudes, and
``ratQ:overflow`` is never raised by arithmetic. Interval endpoints are ``Fraction`` (dyadic), so
bisection is exact at any depth.

Conventions: a univariate polynomial is a list of ints, HIGHEST degree first (MATLAB's); the empty
list is the zero polynomial. A conic is [a b c d e f] for a x^2 + b xy + c y^2 + d x + e y + f.
A face row is 10 wide in the weighted cubic basis: slots 5..10 (1-based) hold
(c5 x^2 + 2 c6 xy + c7 y^2)/2 + c8 x + c9 y + c10.
"""
from fractions import Fraction
from math import gcd, floor, isfinite, inf
import math

from .errors import error, CCA2Error

LIMIT = 2 ** 53

# EMULATE_2P53: off by default. Turned on (by the parity tests only) it makes the routines whose
# MATLAB twins have NO BigInteger fallback raise ``ratQ:overflow`` where MATLAB does, so that a
# difference from a MATLAB golden can be PROVED to come from the 2^53 wall rather than assumed.
EMULATE_2P53 = False


def wall(v, bound=None, what='value'):
    """Raise ratQ:overflow where MATLAB's chk/chkB would -- only under EMULATE_2P53."""
    if EMULATE_2P53:
        m = max((abs(x) for x in _flat(v)), default=0)
        if bound is not None:
            m = max(m, max((abs(x) for x in _flat(bound)), default=0))
        if m >= LIMIT:
            error('ratQ:overflow', '%s reached 2^53 (emulated MATLAB wall).' % what)
    return v


# --------------------------------------------------------------------------------------------
# small helpers (no MATLAB counterpart: they spell MATLAB's matrix expressions on lists)
# --------------------------------------------------------------------------------------------
def sign(x):
    return (x > 0) - (x < 0)


def lcm(a, b):
    a = abs(a)
    b = abs(b)
    if a == 0 or b == 0:
        return 0
    return a // gcd(a, b) * b


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def mv(M, v):
    """M * v for a list-of-rows matrix and a vector."""
    return [dot(row, v) for row in M]


def quad(u, M, v):
    """u' * M * v."""
    return dot(u, mv(M, v))


def _isint(v):
    return isinstance(v, int) and not isinstance(v, bool)


def _flat(v):
    if isinstance(v, (list, tuple)):
        for x in v:
            yield from _flat(x)
    else:
        yield v


def toInt(v):
    """An exact integer from an int, an integral Fraction or an integral float; else None."""
    if _isint(v):
        return v
    if isinstance(v, Fraction):
        return v.numerator if v.denominator == 1 else None
    if isinstance(v, float):
        if isfinite(v) and v == int(v):
            return int(v)
        return None
    try:  # numpy scalars
        f = float(v)
        if isfinite(f) and f == int(f):
            return int(f)
    except (TypeError, ValueError):
        pass
    return None


def chk(v, what="value"):
    """Assert every entry is an exact integer. Returns v, so it can wrap an expression.

    MATLAB also bounds the magnitude by 2^53; the port has no such wall (see the module header).
    """
    for x in _flat(v):
        if _isint(x):
            continue
        if isinstance(x, float) and not isfinite(x):
            error('ratQ:notFinite', '%s must be real and finite.' % what)
        error('ratQ:notInteger', '%s must be integral; got a non-integer entry.' % what)
    return v


def chkB(v, bound=None, what="value"):
    """MATLAB bounds the terms that produced v; with unbounded ints only integrality remains."""
    return chk(v, what)


# --------------------------------------------------------------------------------------------
# canonical forms
# --------------------------------------------------------------------------------------------
def canon(n, d):
    """Canonical form of the rational vector n/d: gcd removed, denominator positive.

    n is a flat list (or a list of rows, returned in the same shape); d a nonzero int.
    """
    nested = bool(n) and isinstance(n[0], (list, tuple))
    flat = [x for row in n for x in row] if nested else list(n)
    chk(flat, 'numerator')
    chk(d, 'denominator')
    if d == 0:
        error('ratQ:badDenominator', 'the denominator must be a nonzero scalar.')
    if d < 0:
        flat = [-x for x in flat]
        d = -d
    g = d
    for x in flat:
        if x != 0:
            g = gcd(g, abs(x))
        if g == 1:
            break
    if g > 1:
        flat = [x // g for x in flat]
        d = d // g
    if nested:
        w = len(n[0])
        return [flat[i * w:(i + 1) * w] for i in range(len(n))], d
    return flat, d


def conic(c):
    """Canonical form of a conic: primitive, first nonzero entry positive."""
    c = list(c)
    chk(c, 'conic')
    nz = next((i for i, x in enumerate(c) if x != 0), None)
    if nz is None:
        error('ratQ:zeroConic',
              'an all-zero conic names no curve. In a QuaCon every edge carries its own '
              'equation, and a straight edge is [0 0 0 d e f], not all-zero.')
    g = 0
    for x in c:
        if x != 0:
            g = gcd(g, abs(x))
    c = [x // g for x in c]
    if c[nz] < 0:
        c = [-x for x in c]
    return c


def sameConic(c1, c2):
    return conic(c1) == conic(c2)


def add(n1, d1, n2, d2):
    """(n1/d1) + (n2/d2), canonically."""
    l = lcm(d1, d2)
    n = [a * (l // d1) + b * (l // d2) for a, b in zip(n1, n2)]
    return canon(n, l)


def sub(n1, d1, n2, d2):
    return add(n1, d1, [-x for x in n2], d2)


def scale(n1, d1, kn, kd=1):
    return canon([x * kn for x in n1], d1 * kd)


def eqRat(n1, d1, n2, d2):
    a, p = canon(n1, d1)
    b, q = canon(n2, d2)
    return a == b and p == q


def diffConic(n1, d1, n2, d2):
    """The conic on which two QUADRATIC faces agree, canonically."""
    if any(x != 0 for x in n1[0:4]) or any(x != 0 for x in n2[0:4]):
        error('ratQ:cubicFace',
              'diffConic is for QUADRATIC faces; a cubic numerator has no conic level set.')
    dn, _ = sub(n1, d1, n2, d2)
    c = [dn[4], 2 * dn[5], dn[6], 2 * dn[7], 2 * dn[8], 2 * dn[9]]
    return conic(c)


def discriminant(c):
    chk(c, 'conic')
    return c[1] ** 2 - 4 * c[0] * c[2]


def conicKind(c):
    if all(x == 0 for x in c[0:3]):
        return 'line'
    d = discriminant(c)
    if d == 0:
        return 'parabola'
    if d > 0:
        return 'hyperbolic'
    return 'elliptic'


def signQ(n, d):
    chk(n, 'numerator')
    chk(d, 'denominator')
    if d == 0:
        error('ratQ:badDenominator', 'the denominator must be a nonzero scalar.')
    if isinstance(n, (list, tuple)):
        return [sign(x) * sign(d) for x in n]
    return sign(n) * sign(d)


def div(n1, d1, n2, d2):
    if isinstance(n2, (list, tuple)):
        error('ratQ:badDivisor',
              'the divisor must be a scalar rational; dividing by a vector is not defined.')
    if n2 == 0:
        error('ratQ:divideByZero',
              'division by the exact zero: a defect in the caller, not a case to regularise.')
    return scale(n1, d1, d2, n2)


def combineDen(Ns, ds):
    """Rows Ns[i]/ds[i] rewritten over ONE common denominator. Returns (N, d)."""
    chk(Ns, 'numerators')
    chk(ds, 'denominators')
    ds = list(ds)
    if len(Ns) != len(ds):
        error('ratQ:sizeMismatch',
              'combineDen needs one denominator per row: %d rows, %d denominators.'
              % (len(Ns), len(ds)))
    if any(x == 0 for x in ds):
        error('ratQ:badDenominator', 'every denominator must be nonzero.')
    d = 1
    for x in ds:
        d = lcm(d, abs(x))
    # MATLAB: Ns .* (d ./ ds) -- a negative denominator flips the row, as there
    N = [[x * (d // abs(di)) * (1 if di > 0 else -1) for x in row] for row, di in zip(Ns, ds)]
    return N, d


def detExact(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    if n == 2:
        return A[0][0] * A[1][1] - A[0][1] * A[1][0]
    if n == 3:
        return (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]))
    error('ratQ:badSystem', 'detExact handles n <= 3; got %d.' % n)


detBig = detExact


def solve2(A, b):
    """Solve A x = b exactly (2x2 or 3x3) by Cramer. Returns (x, xd): x a list, x/xd the solution."""
    chk(A, 'system matrix')
    chk(b, 'right-hand side')
    n = len(A)
    if n not in (2, 3) or any(len(r) != n for r in A) or len(b) != n:
        error('ratQ:badSystem', 'solve2 takes a 2x2 or 3x3 matrix with a matching n x 1 column.')
    D = detExact(A)
    if D == 0:
        error('ratQ:singular',
              'the system is exactly singular (determinant 0), so it has no unique solution. '
              'The caller must handle the degenerate configuration -- parallel edges, a '
              'collapsed cell -- rather than perturbing it.')
    x = []
    for j in range(n):
        Aj = [list(r) for r in A]
        for i in range(n):
            Aj[i][j] = b[i]
        x.append(detExact(Aj))
    return canon(x, D)


def hessQ(n, dn):
    if any(x != 0 for x in n[0:4]):
        error('ratQ:cubicFace',
              'hessQ is for QUADRATIC faces: a cubic Hessian varies with the point, so there '
              'is no constant matrix to return.')
    return canon([[n[4], n[5]], [n[5], n[6]]], dn)


def gradQ(n, dn):
    if any(x != 0 for x in n[0:4]):
        error('ratQ:cubicFace',
              'gradQ is for QUADRATIC faces: a cubic gradient is quadratic, not affine, so it '
              'does not fit the 2x3 map.')
    return canon([[n[4], n[5], n[7]], [n[5], n[6], n[8]]], dn)


# --------------------------------------------------------------------------------------------
# linear feasibility in the plane (Fourier-Motzkin), exact
# --------------------------------------------------------------------------------------------
def feasible2(P, strict=False):
    """Is { s : a*s1 + b*s2 + c >= 0 for every row [a b c] } nonempty (strict: has interior)?

    Exact and complete: eliminate s2 by Fourier-Motzkin, then compare the 1-D bounds as fractions.
    """
    chk(P, 'half-plane')
    if not P:
        return True
    Q = []                                       # rows [a b] meaning a*s1 + b >= 0
    lo = [r for r in P if r[1] > 0]
    up = [r for r in P if r[1] < 0]
    for r in P:
        if r[1] == 0:
            Q.append((r[0], r[2]))
    for ri in lo:
        for rj in up:
            Q.append((-rj[1] * ri[0] + ri[1] * rj[0], -rj[1] * ri[2] + ri[1] * rj[2]))
    haveLo = haveHi = False
    loN = loD = hiN = hiD = None
    for a, b in Q:
        if a == 0:
            if (strict and b <= 0) or (not strict and b < 0):
                return False
        elif a > 0:
            if not haveLo or (-b * loD - loN * a) > 0:
                loN, loD, haveLo = -b, a, True
        else:
            nN, nD = b, -a
            if not haveHi or (nN * hiD - hiN * nD) < 0:
                hiN, hiD, haveHi = nN, nD, True
    if haveLo and haveHi:
        c = loN * hiD - hiN * loD
        if (strict and c >= 0) or (not strict and c > 0):
            return False
    return True


feasible2Double = feasible2
feasible2Big = feasible2


# --------------------------------------------------------------------------------------------
# definiteness
# --------------------------------------------------------------------------------------------
def _conicMatrix(c):
    return [[2 * c[0], c[1], c[3]],
            [c[1], 2 * c[2], c[4]],
            [c[3], c[4], 2 * c[5]]]


def isPSD2(H):
    return H[0][0] >= 0 and H[1][1] >= 0 and detExact(H) >= 0


def isPSD3(M):
    if any(M[i][i] < 0 for i in range(3)):
        return False
    for p in ((0, 1), (0, 2), (1, 2)):
        if detExact([[M[p[0]][p[0]], M[p[0]][p[1]]], [M[p[1]][p[0]], M[p[1]][p[1]]]]) < 0:
            return False
    return detExact(M) >= 0


def conicSign(c):
    """+1 if the conic's quadratic form is >= 0 on the whole plane, -1 if <= 0, else 0."""
    chk(c, 'conic')
    M = _conicMatrix(c)
    if isPSD3(M):
        return 1
    if isPSD3([[-x for x in r] for r in M]):
        return -1
    return 0


def _leadingPositive(M):
    return (M[0][0] > 0
            and detExact([[M[0][0], M[0][1]], [M[1][0], M[1][1]]]) > 0
            and detExact(M) > 0)


def conicDefinite(c):
    chk(c, 'conic')
    M = _conicMatrix(c)
    return _leadingPositive(M) or _leadingPositive([[-x for x in r] for r in M])


# --------------------------------------------------------------------------------------------
# univariate integer polynomials, highest degree first
# --------------------------------------------------------------------------------------------
def polyTrim(p):
    p = list(p)
    k = next((i for i, x in enumerate(p) if x != 0), None)
    if k is None:
        return []
    return p[k:]


def polyPrimitive(p):
    p = polyTrim(p)
    if not p:
        return p
    g = 0
    for x in p:
        if x != 0:
            g = gcd(g, abs(x))
    if g > 1:
        p = [x // g for x in p]
    return p


def polyder(p):
    """MATLAB polyder: the derivative, [0] for a constant (never empty)."""
    n = len(p) - 1
    if n <= 0:
        return [0]
    return [c * (n - i) for i, c in enumerate(p[:-1])]


def polyval(p, x):
    """Horner; exact when x is an int or a Fraction."""
    v = 0
    for c in p:
        v = v * x + c
    return v


def polyMul1(a, b):
    """Product of two univariate polynomials (MATLAB conv)."""
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def polyAdd1(a, b):
    n = max(len(a), len(b))
    a = [0] * (n - len(a)) + list(a)
    b = [0] * (n - len(b)) + list(b)
    return [x + y for x, y in zip(a, b)]


def polyDivExact1(p, g):
    """p / g for integer polynomials when the division is exact; None when it is not."""
    p = polyTrim(p)
    g = polyTrim(g)
    if not g:
        return None
    if len(p) < len(g):
        return [] if not p else None
    r = list(p)
    q = []
    for i in range(len(p) - len(g) + 1):
        c, rem = divmod(r[i], g[0])
        if rem != 0:
            return None
        q.append(c)
        if c != 0:
            for j, gj in enumerate(g):
                r[i + j] -= c * gj
    if any(x != 0 for x in r):
        return None
    return q


def polyPremPrimitive(a, b):
    """The primitive pseudo-remainder used by polyGcd (content removed at every step)."""
    a = polyPrimitive(a)
    b = polyPrimitive(b)
    if not b:
        error('ratQ:divideByZero', 'pseudo-remainder by the zero polynomial.')
    db = len(b) - 1
    r = a
    guard = 0
    while True:
        r = polyTrim(r)
        if not r or len(r) - 1 < db:
            return r
        guard += 1
        if guard > 10000:
            error('ratQ:premNoConverge', 'the reduced pseudo-remainder did not terminate.')
        lcb = b[0]
        lcr = r[0]
        r = [x * lcb for x in r]
        for i in range(len(b)):
            r[i] -= lcr * b[i]
        r[0] = 0
        r = polyPrimitive(polyTrim(r))


def polyPseudoRem(a, b):
    """prem(a, b) = lc(b)^(deg a - deg b + 1) * a  mod  b, then made primitive (sign kept)."""
    a = polyPrimitive(a)
    b = polyPrimitive(b)
    if not b:
        error('ratQ:divideByZero', 'pseudo-remainder by the zero polynomial.')
    da = len(a) - 1
    db = len(b) - 1
    if da < db:
        return a
    f = b[0] ** (da - db + 1)
    r = [x * f for x in a]
    while True:
        r = polyTrim(r)
        if not r or len(r) - 1 < db:
            break
        c, rem = divmod(r[0], b[0])
        if rem != 0:
            error('ratQ:pseudoRemNotIntegral',
                  'the pseudo-remainder step did not divide evenly, which the '
                  'lc(b)^(deg a - deg b + 1) scaling guarantees it must.')
        for k in range(len(b)):
            r[k] -= c * b[k]
        r[0] = 0
    return polyPrimitive(polyTrim(r))


polyPseudoRemDouble = polyPseudoRem
polyPseudoRemBig = polyPseudoRem


def polyGcd(a, b):
    """gcd of two integer polynomials, primitive; [1] when coprime."""
    a = polyPrimitive(a)
    b = polyPrimitive(b)
    if not a:
        return b
    if not b:
        return a
    while True:
        r = polyTrim(polyPremPrimitive(a, b))
        if not r:
            return polyPrimitive(b)
        a = b
        b = polyPrimitive(r)
        if len(b) <= 1:
            return [1]


def polySquarefree(p):
    """The squarefree part of p (same roots, each simple)."""
    p = polyTrim(p)
    if len(p) <= 2:
        return p
    g = polyGcd(p, polyder(p))
    if len(g) <= 1:
        return p
    q = polyDivExact1(p, g)
    if q is None:
        return p
    return polyTrim(q)


def sturmChain(p):
    """A Sturm chain of the squarefree part of p.

    PORT NOTE: MATLAB scales the pseudo-remainder by lc(b)^(gap+1) WITH its sign, so a negative
    leading coefficient with an even degree gap flips the remainder's sign and breaks the Sturm
    property (its BigInteger twin uses |lc|). The sign is restored here; see PORT_NOTES.md.
    """
    p = polyTrim(p)
    p = polySquarefree(p)
    chain = [p]
    if len(p) <= 1:
        return chain
    chain.append(polyTrim(polyder(p)))
    while len(chain[-1]) > 1:
        a = polyPrimitive(chain[-2])
        b = polyPrimitive(chain[-1])
        r = polyTrim(polyPseudoRem(a, b))
        if not r:
            break
        if b[0] < 0 and (len(a) - len(b) + 1) % 2 == 1:
            r = [-x for x in r]
        chain.append(polyPrimitive([-x for x in r]))
    return chain


def _signAtInf(c, x):
    """Sign of the polynomial c (highest first) at x = +inf or -inf."""
    if not c:
        return 0
    s = sign(c[0])
    if x < 0 and (len(c) - 1) % 2 == 1:
        s = -s
    return s


def sturmV(chain, x):
    """Number of sign variations of the chain at x (zeros skipped). x: int, Fraction or +-inf."""
    v = 0
    last = 0
    for c in chain:
        if isinstance(x, float) and not isfinite(x):
            s = _signAtInf(c, x)
        else:
            s = sign(polyval(c, x))
        if s == 0:
            continue
        if last != 0 and s != last:
            v += 1
        last = s
    return v


sturmVBig = sturmV
sturmVBigInf = sturmV


def toFrac(x):
    """An interval endpoint as an exact Fraction (floats are dyadic, so this is exact)."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, float) and not isfinite(x):
        return x
    return Fraction(x)


def countRootsIn(p, a, b):
    """Number of distinct real roots of p in the half-open interval (a, b]. Exact."""
    chain = sturmChain(p)
    if len(chain[0]) <= 1:
        return 0
    return sturmV(chain, toFrac(a)) - sturmV(chain, toFrac(b))


countRootsInHalfOpenBig = countRootsIn
countRootsInBig = countRootsIn


def _bisect(chain, lo, hi, n, depth):
    if n == 0:
        return []
    if n == 1:
        return [(lo, hi)]
    if depth > 200:
        error('ratQ:isolateNoConverge',
              'bisection did not separate %d roots in 200 levels. For a SQUAREFREE integer '
              'polynomial the roots are distinct and separated by a positive distance, so this '
              'is a defect rather than a precision limit.' % n)
    mid = (lo + hi) / 2
    vL = sturmV(chain, lo)
    vM = sturmV(chain, mid)
    vH = sturmV(chain, hi)
    return _bisect(chain, lo, mid, vL - vM, depth + 1) + _bisect(chain, mid, hi, vM - vH, depth + 1)


def isolateRoots(p):
    """Isolating intervals (lo, hi] for the real roots of p, sorted. Endpoints are Fractions."""
    p = polySquarefree(polyTrim(p))
    if len(p) <= 1:
        return []
    # MATLAB: M = ceil(1 + max(abs(p(2:end)))/abs(p(1))) -- the Cauchy bound
    M = 1 + -(-max(abs(x) for x in p[1:]) // abs(p[0]))
    chain = sturmChain(p)
    lo = Fraction(-M)
    hi = Fraction(M)
    I = _bisect(chain, lo, hi, sturmV(chain, lo) - sturmV(chain, hi), 0)
    return sorted(I)


isolateRootsBig = isolateRoots


def signAt(q, p, iv):
    """The sign of q at the root of p isolated by iv = (lo, hi]. Exact."""
    p = polySquarefree(polyTrim(p))
    q = polyTrim(q)
    if not q:
        return 0
    lo = toFrac(iv[0])
    hi = toFrac(iv[1])
    g = polyGcd(p, q)
    if len(g) > 1 and countRootsIn(g, lo, hi) > 0:
        return 0                               # alpha is a common root: q(alpha) = 0 exactly
    chainP = sturmChain(p)
    for _ in range(400):
        if countRootsIn(q, lo, hi) == 0:
            v = polyval(q, (lo + hi) / 2)
            if v != 0:
                return sign(v)
            v = polyval(q, hi)
            if v != 0:
                return sign(v)
            return sign(polyval(q, lo))
        mid = (lo + hi) / 2
        if sturmV(chainP, lo) - sturmV(chainP, mid) == 1:
            hi = mid
        else:
            lo = mid
    error('ratQ:signAtNoConverge',
          'bisection did not separate the root from q\'s zeros. The root and a zero of q are '
          'then closer than a dyadic rational of that depth, which for exact integer input '
          'means they are equal -- so the gcd test above should have caught it, and this is a '
          'defect rather than a precision limit.')


signAtBig = signAt


# --------------------------------------------------------------------------------------------
# reading the caller's doubles exactly
# --------------------------------------------------------------------------------------------
def _roundHalfAway(x):
    return int(math.floor(abs(x) + 0.5)) * (1 if x >= 0 else -1)


def denomOf(x, maxDen=10 ** 6):
    """The denominator of the rational a double was MEANT to be (continued fractions)."""
    if isinstance(x, Fraction):
        return x.denominator
    if _isint(x):
        return 1
    x = float(x)
    if not isfinite(x):
        error('ratQ:notFinite', 'value must be real and finite.')
    if x == round(x):
        return 1
    ax = abs(x)
    p0, q0, p1, q1 = 1.0, 0.0, float(floor(ax)), 1.0
    r = ax - p1
    for _ in range(64):
        if r == 0:
            break
        a = float(floor(1 / r))
        r = 1 / r - a
        p2 = a * p1 + p0
        q2 = a * q1 + q0
        if q2 > maxDen:
            break
        p0, q0, p1, q1 = p1, q1, p2, q2
        if p1 / q1 == ax:
            return int(q1)
    error('ratQ:notExact',
          '%.17g is not a rational with denominator <= %d to within a rounding. It is most '
          'likely an irrational produced upstream (a VERTEX, not a coefficient) -- and by '
          'CONJ_FIELD_PROOF.md Theorem 1 no COEFFICIENT of f* can be one, so this is a defect '
          'in the caller and not a case to round away.' % (x, maxDen))


def fromDouble(x, maxDen=10 ** 6):
    """The exact rational vector a caller's numbers stand for: returns (n, d), canonical.

    Accepts ints and Fractions (taken as they are) and floats (reconstructed as MATLAB does, and
    refused with ``ratQ:notExact`` when no denominator <= maxDen reproduces them).
    """
    scalar = not isinstance(x, (list, tuple))
    xs = [x] if scalar else list(x)
    d = 1
    for v in xs:
        d = lcm(d, denomOf(v, maxDen))
        if d > LIMIT:
            error('ratQ:overflow',
                  'the common denominator exceeded 2^53 while converting %.17g.' % float(v))
    n = []
    for v in xs:
        if isinstance(v, Fraction):
            n.append(int(v * d))
        elif _isint(v):
            n.append(v * d)
        else:
            v = float(v)
            ni = _roundHalfAway(v * d)
            if ni / d != v:
                error('ratQ:notExact',
                      '%.17g does not round-trip through the common denominator %d. The vector '
                      'mixes denominators whose lcm is too coarse for one of its entries; '
                      'convert the entries separately, or raise maxDen.' % (v, d))
            n.append(ni)
    n, d = canon(n, d)
    return n, d


# --------------------------------------------------------------------------------------------
# evaluation (doubles on purpose -- see QuaCon.eval)
# --------------------------------------------------------------------------------------------
_CW = (1 / 6, 1 / 2, 1 / 2, 1 / 6, 1 / 2, 1.0, 1 / 2, 1.0, 1.0, 1.0)


def evalPoly(c, X):
    """Value of the weighted cubic row c (6 or 10 wide) at each point of X (list of [x y])."""
    c = [0] * (10 - len(c)) + [float(v) for v in c]
    out = []
    for pt in X:
        x = float(pt[0])
        y = float(pt[1])
        mon = (x ** 3, x * x * y, x * y * y, y ** 3, x * x, x * y, y * y, x, y, 1.0)
        out.append(sum(m * w * ci for m, w, ci in zip(mon, _CW, c)))
    return out


def evalConic(c, X):
    return [c[0] * p[0] ** 2 + c[1] * p[0] * p[1] + c[2] * p[1] ** 2
            + c[3] * p[0] + c[4] * p[1] + c[5] for p in X]


def evalFace(n, d, X):
    return [v / d for v in evalPoly(n, X)]


def evalFaceExact(n, d, x):
    """Exact value of the face n/d at the rational point x = (x1, x2) (Fractions or ints)."""
    x1 = Fraction(x[0])
    x2 = Fraction(x[1])
    v = (Fraction(n[0], 6) * x1 ** 3 + Fraction(n[1], 2) * x1 ** 2 * x2
         + Fraction(n[2], 2) * x1 * x2 ** 2 + Fraction(n[3], 6) * x2 ** 3
         + Fraction(n[4], 2) * x1 ** 2 + n[5] * x1 * x2 + Fraction(n[6], 2) * x2 ** 2
         + n[7] * x1 + n[8] * x2 + n[9])
    return v / d


# --------------------------------------------------------------------------------------------
# the sign of a quadratic on a polyhedron
# --------------------------------------------------------------------------------------------
def recCone2(P):
    """Extreme rays of the recession cone of {P >= 0}; ok False when it is not a pointed cone
    generated by at most two rays (declined)."""
    N = [[r[0], r[1]] for r in P if r[0] != 0 or r[1] != 0]
    if not N:
        return [], False
    C = []
    for n in N:
        C.append([n[1], -n[0]])
        C.append([-n[1], n[0]])
    C = [c for c in C if all(n[0] * c[0] + n[1] * c[1] >= 0 for n in N)]
    U = []
    for c in C:
        new = True
        for u in U:
            if u[0] * c[1] - u[1] * c[0] == 0 and u[0] * c[0] + u[1] * c[1] > 0:
                new = False
                break
        if new:
            U.append(c)
    if len(U) == 2 and U[0][0] * U[1][1] - U[0][1] * U[1][0] == 0:
        return [], False
    if len(U) > 2:
        return [], False
    return U, True


def coneFormMinSign(H, u1, u2):
    A = quad(u1, H, u1)
    B = quad(u1, H, u2)
    C = quad(u2, H, u2)
    a2 = A - 2 * B + C
    a1 = 2 * (B - C)
    s = min(sign(A), sign(C))
    if a2 > 0 and a1 < 0 and -a1 < 2 * a2:
        s = min(s, sign(4 * a2 * C - a1 ** 2))
    return s


def flatRaySign(H, g0, P, u):
    a = mv(H, u)
    b = dot(g0, u)
    if all(x == 0 for x in a):
        return -1 if b < 0 else 0
    if feasible2([list(r) for r in P] + [[-a[0], -a[1], -b]], True):
        return -1
    return 0


def quadSignAtSolve(c, P, A, b):
    s = 1
    if detExact(A) == 0:
        return s
    xn, xd = solve2(A, b)
    if any(r[0] * xn[0] + r[1] * xn[1] + r[2] * xd < 0 for r in P):
        return s
    return sign(c[0] * xn[0] ** 2 + c[1] * xn[0] * xn[1] + c[2] * xn[1] ** 2
                + c[3] * xn[0] * xd + c[4] * xn[1] * xd + c[5] * xd ** 2)


def quadMinSign(c, P):
    """Sign of the minimum of the quadratic c over the polyhedron {P >= 0}: -1, 0, +1, or None
    (MATLAB NaN) when the configuration is declined."""
    H = [[2 * c[0], c[1]], [c[1], 2 * c[2]]]
    g0 = [c[3], c[4]]
    rec, ok = recCone2(P)
    if not ok:
        return None
    if len(rec) == 1:
        a = quad(rec[0], H, rec[0])
        if a < 0:
            return -1
        if a == 0:
            r = flatRaySign(H, g0, P, rec[0])
            if r < 0:
                return r
    elif len(rec) == 2:
        u1, u2 = rec
        t = coneFormMinSign(H, u1, u2)
        if t < 0:
            return -1
        if t == 0:
            f1 = quad(u1, H, u1) == 0
            f2 = quad(u2, H, u2) == 0
            if not f1 and not f2:
                return None
            if f1:
                r = flatRaySign(H, g0, P, u1)
                if r < 0:
                    return r
            if f2:
                r = flatRaySign(H, g0, P, u2)
                if r < 0:
                    return r
    elif len(rec) > 2:
        return None
    s = 1
    k = len(P)
    s = min(s, quadSignAtSolve(c, P, H, [-g0[0], -g0[1]]))
    for i in range(k):
        n = [P[i][0], P[i][1]]
        if n[0] == 0 and n[1] == 0:
            continue
        u = [n[1], -n[0]]
        uH = [u[0] * H[0][0] + u[1] * H[1][0], u[0] * H[0][1] + u[1] * H[1][1]]
        s = min(s, quadSignAtSolve(c, P, [uH, n], [-dot(u, g0), -P[i][2]]))
    for i in range(k - 1):
        for j in range(i + 1, k):
            s = min(s, quadSignAtSolve(c, P, [[P[i][0], P[i][1]], [P[j][0], P[j][1]]],
                                       [-P[i][2], -P[j][2]]))
    return s


def quadSignOnPoly(c, P):
    """(s, strict): s = +1 when q >= 0 on the polyhedron {P >= 0}, -1 when q <= 0, 0 undecided."""
    chk(c, 'conic')
    chk(P, 'half-plane')
    if len(c) != 6 or any(len(r) != 3 for r in P):
        error('ratQ:badQuadSignInput',
              'quadSignOnPoly takes a 1 x 6 conic and k x 3 half-plane rows.')
    mn = quadMinSign(c, P)
    if mn is not None and mn >= 0:
        return 1, mn > 0
    mx = quadMinSign([-x for x in c], P)
    if mx is not None and mx >= 0:
        return -1, mx > 0
    if any(x != 0 for x in c):
        s = conicSign(c)
        return s, (s != 0) and conicDefinite(c)
    return 0, False


# --------------------------------------------------------------------------------------------
# a degenerate conic as a NAMED pair of lines
# --------------------------------------------------------------------------------------------
def parallelPairData(c):
    a, b, cc, d, e, f = c
    if a != 0:
        u, v = 2 * a, b
        if 2 * a * e - d * b != 0:
            error('ratQ:notALinePair',
                  'the linear part is not a multiple of the repeated factor, so this conic is '
                  'a PARABOLA rather than two parallel lines.')
        quad_ = [1, 2 * d, 4 * a * f]
    elif cc != 0:
        u, v = 0, 1
        if d != 0:
            error('ratQ:notALinePair',
                  'the linear part is not a multiple of the repeated factor, so this conic is '
                  'a PARABOLA rather than two parallel lines.')
        quad_ = [cc, e, f]
    else:
        error('ratQ:notALinePair',
              'the quadratic part vanishes, so this conic is a single line, not a pair.')
    quad_ = polyTrim(quad_)
    I = isolateRoots(quad_)
    if len(I) != 2:
        error('ratQ:notALinePair',
              'the level quadratic has %d real roots, not two distinct lines.' % len(I))
    return {'kind': 'parallel', 'apex': [0, 0, 1], 'quad': quad_, 'iv': I,
            'vertical': False, 'form': [u, v]}


def pairLineData(c):
    """The two lines of a degenerate conic, NAMED: branch b is the line whose slope (or level)
    is the b-th root of info['quad'], isolated by info['iv'][b-1]."""
    c = conic(c)
    a, b, cc, d, e, f = c
    disc = b * b - 4 * a * cc
    if disc < 0:
        error('ratQ:notALinePair',
              'this conic\'s quadratic part is DEFINITE, so its zero set is a point or empty '
              '-- not two lines, and there is no branch to name.')
    if disc == 0:
        return parallelPairData(c)
    M3 = [[2 * a, b, d], [b, 2 * cc, e], [d, e, 2 * f]]
    if detExact(M3) != 0:
        error('ratQ:notDegenerate',
              'this conic is NOT degenerate, so it is one curve and not two lines. A branch '
              'index names nothing on it.')
    M2 = [[2 * a, b], [b, 2 * cc]]
    if detExact(M2) == 0:
        error('ratQ:parallelPair',
              'the two lines are PARALLEL, which the quadratic part being indefinite excludes '
              '-- so this cannot be reached and reaching it is a defect here.')
    xn, xd = solve2(M2, [-d, -e])
    if xd < 0:
        xn = [-x for x in xn]
        xd = -xd
    quad_ = polyTrim([cc, b, a])
    vertical = (cc == 0)
    if vertical:
        if b == 0:
            error('ratQ:notALinePair',
                  'both lines are vertical, so they are parallel and cannot cross.')
        t = Fraction(-a, b)
        iv = [(t - 1, t + 1), (Fraction(0), Fraction(0))]
    else:
        I = isolateRoots(quad_)
        if len(I) != 2:
            error('ratQ:notALinePair',
                  'the slope quadratic has %d real roots, not two, although the discriminant '
                  'is positive.' % len(I))
        iv = I
    return {'kind': 'crossing', 'apex': [xn[0], xn[1], xd], 'quad': quad_, 'iv': iv,
            'vertical': vertical}


def pairLineSign(c, b, xn, xd):
    """Exact sign of branch b of the line pair c at the rational point xn/xd."""
    info = pairLineData(c)
    if b not in (1, 2):
        error('ratQ:branch', 'the branch must be 1 or 2; got %r.' % (b,))
    if info['kind'] == 'parallel':
        tv = info['form'][0] * xn[0] + info['form'][1] * xn[1]
        return signAt(polyTrim([-xd, tv]), info['quad'], info['iv'][b - 1])
    p = info['apex']
    ux = xn[0] * p[2] - p[0] * xd
    uy = xn[1] * p[2] - p[1] * xd
    if info['vertical'] and b == 2:
        return sign(ux)
    if info['vertical']:
        td = info['quad'][0]
        tn = -info['quad'][1]
        return sign(uy * td - tn * ux) * sign(td)
    q = polyTrim([-ux, uy])
    if not q:
        return 0
    return signAt(q, info['quad'], info['iv'][b - 1])


_pairLineMemo = {}


def _sortedRealRoots2(quad_):
    """sort(roots(quad)) for a quadratic with two real roots, in doubles (an eval-only cache)."""
    a, b, c = [float(v) for v in quad_]
    disc = b * b - 4 * a * c
    sq = math.sqrt(max(disc, 0.0))
    q = -0.5 * (b + math.copysign(sq, b)) if b != 0 else 0.5 * sq
    if b != 0:
        r1 = q / a
        r2 = c / q if q != 0 else r1
    else:
        r1 = sq / (2 * a)
        r2 = -r1
    return sorted([r1, r2])


def pairLineEval(c, b, X):
    """Double-precision value of branch b of the line pair c at the points X (eval's cache)."""
    key = tuple(c[0:6])
    info = _pairLineMemo.get(key)
    if info is None:
        info = pairLineData(c)
        if len(_pairLineMemo) > 20000:
            _pairLineMemo.clear()
        _pairLineMemo[key] = info
    if info['kind'] == 'parallel':
        r = _sortedRealRoots2(info['quad'])
        return [info['form'][0] * p[0] + info['form'][1] * p[1] - r[b - 1] for p in X]
    p0 = (info['apex'][0] / info['apex'][2], info['apex'][1] / info['apex'][2])
    if info['vertical'] and b == 2:
        return [p[0] - p0[0] for p in X]
    if info['vertical']:
        t = -info['quad'][1] / info['quad'][0]
    else:
        t = _sortedRealRoots2(info['quad'])[b - 1]
    return [(p[1] - p0[1]) - t * (p[0] - p0[0]) for p in X]


# --------------------------------------------------------------------------------------------
# Q(t): elements [n1 n2 D] = (n1 + n2 t)/D with t a root of Q = [al be ga]
# --------------------------------------------------------------------------------------------
def qsRed(u):
    if len(u) != 3:
        error('ratQ:qsShape', 'a Q(t) element is [n1 n2 D]; got length %d.' % len(u))
    if u[2] == 0:
        error('ratQ:qsZeroDen', 'a Q(t) element cannot have a zero denominator.')
    u = list(u)
    if u[2] < 0:
        u = [-x for x in u]
    g = gcd(gcd(abs(u[0]), abs(u[1])), abs(u[2]))
    if g > 0:
        u = [x // g for x in u]
    return chk(u, 'Q(t) element')


def qsAdd(u, v):
    return qsRed([u[0] * v[2] + v[0] * u[2], u[1] * v[2] + v[1] * u[2], u[2] * v[2]])


def qsMul(u, v, Q):
    if len(Q) == 2:
        if u[1] != 0 or v[1] != 0:
            error('ratQ:qsField',
                  'a linear Q means t does not occur, but an operand has a t term.')
        return qsRed([u[0] * v[0], 0, u[2] * v[2]])
    al, be, ga = Q
    return qsRed([al * u[0] * v[0] - ga * u[1] * v[1],
                  al * (u[0] * v[1] + u[1] * v[0]) - be * u[1] * v[1],
                  al * u[2] * v[2]])


def qsSign(u, Q, iv):
    u = qsRed(u)
    if u[1] == 0:
        return sign(u[0])
    if len(Q) == 2:
        error('ratQ:qsField', 'a linear Q means t does not occur, but u has a t term.')
    return signAt(polyTrim([u[1], u[0]]), Q, iv)


def qs2Red(u):
    if len(u) != 5:
        error('ratQ:qs2Shape', 'a Q(t1,t2) element is [n1 n2 n3 n4 D]; got length %d.' % len(u))
    if u[4] == 0:
        error('ratQ:qs2ZeroDen', 'a Q(t1,t2) element cannot have a zero denominator.')
    u = list(u)
    if u[4] < 0:
        u = [-x for x in u]
    g = 0
    for x in u:
        g = gcd(g, abs(x))
    if g > 0:
        u = [x // g for x in u]
    return chk(u, 'Q(t1,t2) element')


def qs2Add(u, v):
    return qs2Red([u[0] * v[4] + v[0] * u[4], u[1] * v[4] + v[1] * u[4],
                   u[2] * v[4] + v[2] * u[4], u[3] * v[4] + v[3] * u[4], u[4] * v[4]])


def qs2Mul(u, v, Q1, Q2):
    al1, be1, ga1 = Q1
    al2, be2, ga2 = Q2
    a = u[0:4]
    b = v[0:4]
    c00 = a[0] * b[0]
    c10 = a[0] * b[1] + a[1] * b[0]
    c01 = a[0] * b[2] + a[2] * b[0]
    c11 = a[0] * b[3] + a[3] * b[0] + a[1] * b[2] + a[2] * b[1]
    c20 = a[1] * b[1]
    c02 = a[2] * b[2]
    c21 = a[1] * b[3] + a[3] * b[1]
    c12 = a[2] * b[3] + a[3] * b[2]
    c22 = a[3] * b[3]
    n1 = al1 * al2 * c00 - al2 * ga1 * c20 - al1 * ga2 * c02 + ga1 * ga2 * c22
    n2 = al1 * al2 * c10 - al2 * be1 * c20 - al1 * ga2 * c12 + be1 * ga2 * c22
    n3 = al1 * al2 * c01 - al1 * be2 * c02 - al2 * ga1 * c21 + ga1 * be2 * c22
    n4 = al1 * al2 * c11 - al2 * be1 * c21 - al1 * be2 * c12 + be1 * be2 * c22
    return qs2Red([n1, n2, n3, n4, al1 * al2 * u[4] * v[4]])


def qs2Interval(u, i1, i2):
    c = [float(x) for x in u[0:4]]
    xs = [float(i1[0]), float(i1[1])]
    ys = [float(i2[0]), float(i2[1])]
    v = [c[0] + c[1] * a + c[2] * b + c[3] * a * b for a in xs for b in ys]
    eps = 2.0 ** -52
    tol = 8 * eps * max(1.0, max(abs(x) for x in c)) * max(1.0, max(abs(x) for x in xs + ys)) ** 2
    if min(v) > tol:
        return 1, True
    if max(v) < -tol:
        return -1, True
    return 0, False


def qs2Bisect(Q, lo, hi):
    lo = Fraction(lo)
    hi = Fraction(hi)
    mid = (lo + hi) / 2
    sm = sign(polyval(Q, mid))
    sl = sign(polyval(Q, lo))
    if sm == 0:
        return mid, mid
    if sm == sl:
        return mid, hi
    return lo, mid


def qs2Sign(u, Q1, iv1, Q2, iv2):
    """(s, ok): the sign of an element of Q(t1, t2); ok False when it could not be decided."""
    u = qs2Red(u)
    d1 = Q1[1] ** 2 - 4 * Q1[0] * Q1[2]
    d2 = Q2[1] ** 2 - 4 * Q2[0] * Q2[2]
    if d1 > 0 and d2 > 0:
        r = math.isqrt(d1 * d2)
        if r * r == d1 * d2:
            return 0, False                    # one field in disguise: the zero test fails
    if all(x == 0 for x in u[1:4]):
        return sign(u[0]), True
    if all(x == 0 for x in u[0:4]):
        return 0, True
    lo1, hi1 = Fraction(iv1[0]), Fraction(iv1[1])
    lo2, hi2 = Fraction(iv2[0]), Fraction(iv2[1])
    for _ in range(60):
        s, decided = qs2Interval(u, (lo1, hi1), (lo2, hi2))
        if decided:
            return s, True
        lo1, hi1 = qs2Bisect(Q1, lo1, hi1)
        lo2, hi2 = qs2Bisect(Q2, lo2, hi2)
    return 0, False


def toRational(x):
    f = Fraction(x)
    return f.numerator, f.denominator

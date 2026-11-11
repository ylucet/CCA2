"""alg_row_fails -- a CERTIFICATE that a set of algebraic rows cannot all hold anywhere in the
part of a box that lies in a region (port of algRowFails.m).

A row is "side * (q - z) >= 0" with q a rational quadratic and z a named root of an implicit
face polynomial FP (or "z_k <= z_w" between two such roots, a row with field 'W'). For a fibre
of degree 1 or 2 in z the row's failure is decided by the SIGNS of a few bivariate polynomials
(the leading coefficient, the value and slope of the fibre at q, its discriminant); those signs
are certified on a box by Bernstein enclosures (poly_sign_on), and the box is subdivided until
every sub-box is either certified or misses the region. Square-free parts (Yun) are taken first,
since a squared factor makes a sign unprovable on the curve where it vanishes; polynomials that
vanish at a corner of the region are handled by BLOWING UP the corner.

True is a proof that the rows hold nowhere; False means "not proved" (never "holds").

PORT NOTES. One exact path on unbounded ints replaces MATLAB's double / BigInteger / symbolic
ladder (the Symbolic-Toolbox `factor` fallback is never needed: Yun's decomposition cannot
overflow here). The guards MATLAB places on its DOUBLE blow-up (bound < 2^52) are kept, so the
port proves what MATLAB proves. The rows against the TOP root of a face of degree >= 3 in z
(topRowPolys / topFaceRow) are ported too: the poly2 kink's bridge faces are cubic in z, and
without them the port could not prove what MATLAB proves there (334 faces against 317).
"""
import math

from . import _poly
from . import biv
from . import rat_q
from .errors import CCA2Error
from .numeric import roots
from .poly_sign_on import polySignOn
from .rat_q import sign

EPS = 2.0 ** -52
REALMIN = 2.2250738585072014e-308
_ALWAYS = (lambda _b: True)


# ---- polynomial values ---------------------------------------------------------------------------
def _num(M):
    """A row polynomial with integer (or dyadic) coefficients."""
    return {'t': 'num', 'M': M}


def _slices(FP):
    """The z-coefficients of FP(x, y, z) as full nx x ny matrices, lowest power first."""
    nx = max(e[0] for e in FP) + 1
    ny = max(e[1] for e in FP) + 1
    nz = max(e[2] for e in FP) + 1
    S = [biv.zeros(nx, ny) for _ in range(nz)]
    for e, c in FP.items():
        S[e[2]][e[0]][e[1]] = c
    return S


def _polyKey(FP):
    return tuple(sorted(FP.items()))


def _matKey(M):
    return tuple(tuple(r) for r in M)


# ---- the entry point -----------------------------------------------------------------------------
def algRowFails(a, box, inBox, depth=6, factors=None, corners=None):
    factors = [_asIntMatrix(f) for f in (factors or [])]
    corners = corners or []
    rows = []
    for at in a:
        r = _rowPolys(at, factors)
        if r is not None:
            rows.append(r)
    if not rows:
        return False

    tf = _boxRec(rows, list(box), inBox, depth, [])
    if tf or not corners:
        return tf

    # ---- corners where some sign polynomial VANISHES: blow them up --------------------------
    W = max(box[1] - box[0], box[3] - box[2])
    sq = []
    for p in corners:
        p = [float(p[0]), float(p[1])]
        q = [_round(p[0] * 256) / 256, _round(p[1] * 256) / 256]
        if all(abs(qi - pi) <= 1e-9 * max(1.0, abs(pi)) for qi, pi in zip(q, p)):
            p = q
        if p[0] < box[0] or p[0] > box[1] or p[1] < box[2] or p[1] > box[3]:
            continue
        if not _anyVanishes(rows, p):
            continue                                # nothing to blow up
        for k in range(2, 9):
            e = W / 2 ** k
            if _cornerOK(rows, p, e, factors):
                if not _isDyadic(p):
                    e = e * (1 - 1e-9)
                sq.append([p[0], p[1], e])
                break
    if not sq:
        return False
    tf = _boxRec(rows, list(box), inBox, depth, sq)
    if not tf:
        tf = _boxRec(rows, list(box), inBox, depth + 2, sq)
    return tf


def _round(x):
    return math.floor(abs(x) + 0.5) * (1 if x >= 0 else -1)


def _asIntMatrix(g):
    """A factor (the region's row polynomial) as given; integrality is checked where used."""
    return [list(r) for r in g]


# ---- the sign polynomials of one row ---------------------------------------------------------------
_memoSF = {}
_memoFace = {}
_memoSurd = {}
_memoSurdFace = {}
_memoFactorBU = {}


def _rowPolys(a, factors):
    if a.get('W'):
        return _faceRowPolys(a, factors)
    FP = a['FP']
    num = list(a['num'])
    if len(num) != 10 or any(x != 0 for x in num[0:4]) or not (a['den'] > 0):
        return None
    Q = [[2 * num[9], 2 * num[8], num[6]],
         [2 * num[7], 2 * num[5], 0],
         [num[4], 0, 0]]
    Dn = 2 * a['den']
    want = -a['side']                               # the sign q - z must have for the row to fail
    nz = max(e[2] for e in FP) + 1
    FB = a['FB']
    if nz - 1 >= 3:
        # THE TOP ROOT OF A POLYNOMIAL OF DEGREE n >= 3 in z (MATLAB 2026-11-04). Sufficient
        # conditions, by Taylor at q (G_n = lead in z): q < z_top when sign(G_n) G(q) < 0 (a root
        # above q); q > z_top when sign(G_n) G^(j)(q) > 0 for every j < n (no root at or above q).
        if FB != -1:
            return None
        key = ('top', want, _polyKey(FP), tuple(num), a['den'])
        P = _memoSF.get(key)
        if P is None:
            P = _topRowPolys(FP, Q, Dn, nz - 1, want)
            if len(_memoSF) > 5000:
                _memoSF.clear()
            _memoSF[key] = P
        if not P:
            return None
        P = [_num(_divideOut(p, factors)) for p in P]

        def test(s, want=want):
            return s[0] != 0 and all(s[0] * v == want for v in s[1:])
        return {'P': P, 'test': test}
    if nz - 1 == 2:
        hi = FB in (-1, 2)
        if not hi and FB not in (1, -2):
            return None

        def test(s, hi=hi, want=want):
            return _quadOK(s, hi, want)
    elif nz - 1 == 1:
        if FB not in (1, -1):
            return None

        def test(s, want=want):
            return s[0] != 0 and s[0] * s[1] == want
    else:
        return None
    key = (_polyKey(FP), tuple(num), a['den'])
    P = _memoSF.get(key)
    if P is None:
        P = [_squareFreeSign(p) for p in _signPolys(_slices(FP), Q, Dn, nz - 1)]
        if len(_memoSF) > 5000:
            _memoSF.clear()
        _memoSF[key] = P
    P = [_num(_divideOut(p, factors)) for p in P]
    return {'P': P, 'test': test}


def _signPolys(S, Q, Dn, deg):
    c2 = biv.conv2
    pa = biv.padAdd
    if deg == 2:
        A, B, C0 = S[2], S[1], S[0]
        return [A,                                                       # a
                pa(biv.scal(c2(A, Q), 2), biv.scal(B, Dn)),              # Dn (2aq + b)
                pa(pa(c2(A, c2(Q, Q)), biv.scal(c2(B, Q), Dn)), biv.scal(C0, Dn * Dn)),  # Dn^2 F(q)
                pa(c2(B, B), biv.scal(c2(A, C0), -4))]                   # D
    A, C0 = S[1], S[0]
    return [A, pa(c2(A, Q), biv.scal(C0, Dn))]                           # b, Dn (b q + c)


def _topRowPolys(FP, Q, Dn, n, want):
    """[G_n, H_0] (want -1) or [G_n, H_0, ..., H_(n-1)] (want +1), H_j = Dn^(n-j) G^(j)(Q / Dn) by
    Horner, exact, content and squared factors divided out; [] on failure (topRowPolys)."""
    try:
        G = _slices(FP)
        nj = n - 1 if want == 1 else 0
        P = [biv.contentFree(G[n])]
        for j in range(nj + 1):
            def c(i, j=j):
                return math.factorial(i) // math.factorial(i - j)
            R = biv.scal(G[n], c(n))
            for i in range(n - 1, j - 1, -1):
                R = biv.padAdd(biv.conv2(R, Q), biv.scal(G[i], c(i) * Dn ** (n - i)))
            P.append(_squareFreeSign(biv.contentFree(biv.trimN(R))))
        return P
    except CCA2Error:
        return []


# ---- the box recursion ---------------------------------------------------------------------------
def _boxRec(rows, b, inBox, depth, sq):
    for c in sq:
        if b[0] >= c[0] - c[2] and b[1] <= c[0] + c[2] and b[2] >= c[1] - c[2] \
                and b[3] <= c[1] + c[2]:
            return True
    for row in rows:
        s = [_polySign(P, b) for P in row['P']]
        if row['test'](s):
            return True
    if not inBox(b):
        return True                                 # misses the region: nothing to prove here
    if depth == 0:
        return False
    xm = (b[0] + b[1]) / 2
    ym = (b[2] + b[3]) / 2
    return (_boxRec(rows, [b[0], xm, b[2], ym], inBox, depth - 1, sq)
            and _boxRec(rows, [xm, b[1], b[2], ym], inBox, depth - 1, sq)
            and _boxRec(rows, [b[0], xm, ym, b[3]], inBox, depth - 1, sq)
            and _boxRec(rows, [xm, b[1], ym, b[3]], inBox, depth - 1, sq))


def _signM(M, b):
    v = polySignOn(M, b, _ALWAYS, 0)
    return 0 if v is None else v


def _polySign(P, b):
    if P['t'] == 'num':
        return _signM(P['M'], b)
    lo = _signM(P['Lo'], b)
    if lo != 0 and _signM(P['Hi'], b) == lo:
        return lo
    if P['E'] is None:
        return 0
    if P.get('Us') is not None:                     # a surd face row's: reduced, for their signs
        sU = _signM(P['Us'], b)
        sV = _signM(P['Vs'], b)
    else:
        sU = _signM(P['U'], b)
        sV = _signM(P['V'], b)
    sE = _signM(P['E'], b)
    return _sqSign(sU, sV, sE)


def _evalTerms(M, p):
    """(sum, sum of absolute values) of the monomials of M at the point p (real or complex)."""
    v = 0.0
    m = 0.0
    for i, r in enumerate(M):
        xi = p[0] ** i
        for j, c in enumerate(r):
            if c == 0:
                continue
            t = c * xi * p[1] ** j
            v += t
            m += abs(t)
    return v, m


def _vanishesAt(P, p):
    if P['t'] != 'num':
        return _vanishesAt(_num(P['U']), p) and _vanishesAt(_num(P['V']), p)
    v, m = _evalTerms(P['M'], p)
    return abs(v) <= 1e-9 * max(m, REALMIN)


def _anyVanishes(rows, p):
    return any(_vanishesAt(P, p) for row in rows for P in row['P'])


# ---- corners --------------------------------------------------------------------------------------
def _cornerOK(rows, p, e, factors):
    """Does every row fail on the square of half-width e at the corner p? Proved in the four
    blown-up charts (x - p1 = +-r, y - p2 = +-r t  and its transpose)."""
    if not (e > 0):
        return False
    for major in (1, 2):
        for sg in (1, -1):
            rowsC = []
            for row in rows:
                if row.get('corner') is False:
                    continue
                Ps = []
                ok = True
                for P in row['P']:
                    G = _blowUp(P, p, sg, major)
                    if G is None:
                        ok = False
                        break
                    Ps.append(G)
                if ok:
                    rowsC.append({'P': Ps, 'test': row['test']})
            if not rowsC:
                return False
            facC = []
            for g in factors:
                fk = (_matKey(g), p[0], p[1], sg, major)
                if fk not in _memoFactorBU:
                    if len(_memoFactorBU) > 20000:
                        _memoFactorBU.clear()
                    gn = _num(g)
                    _memoFactorBU[fk] = _blowUpRaw(gn, p, sg, major) if _vanishesAt(gn, p) else None
                G = _memoFactorBU[fk]
                if G is not None:
                    facC.append(G)

            def inC(cb, facC=facC):
                return not any(_polySign(G, cb) == -1 for G in facC)
            if not _boxRec(rowsC, [0.0, e, -1.0, 1.0], inC, 6, []):
                return False
    return True


def _isDyadic(p):
    return all(v * 256 == _round(v * 256) for v in p)


def _toIntegers(M):
    """A dyadic matrix scaled to integers (sign kept), or None."""
    M = [list(r) for r in M]
    for _ in range(30):
        if all(x == int(x) for r in M for x in r):
            break
        M = [[2 * x for x in r] for r in M]
    if any(x != int(x) for r in M for x in r):
        return None
    return [[int(x) for x in r] for r in M]


def _blowUp(P, p, sg, major):
    """The polynomial in the chart at corner p, divided by the power of r it vanishes to; None
    when it cannot be formed exactly (MATLAB's double-precision guard). Memoised on the row
    polynomial itself: the same blow-up is asked for at every trial size of the corner square."""
    memo = P.setdefault('_bu', {})
    key = (p[0], p[1], sg, major)
    if key not in memo:
        memo[key] = _blowUpRaw(P, p, sg, major)
    return memo[key]


def _blowUpRaw(P, p, sg, major):
    if P['t'] != 'num':
        return _blowUpPair(P, p, sg, major)
    M = _toIntegers(P['M'])
    if M is None:
        return None
    if not _isDyadic(p):
        return _blowUpSurd(M, p, sg, major)
    sc = 1
    while any(v * sc != _round(v * sc) for v in p) and sc < 256:
        sc *= 2
    n1, n2 = biv.shape(M)
    deg = n1 + n2 - 2
    bound = (sum(abs(x) for r in M for x in r) * (abs(p[0]) + 1) ** (n1 - 1)
             * (abs(p[1]) + 1) ** (n2 - 1) * 2.0 ** deg * float(sc) ** deg)
    if bound >= 2.0 ** 52:
        return None
    from fractions import Fraction
    px = Fraction(p[0])
    py = Fraction(p[1])
    Ps = biv.transpose(_shiftTo(biv.transpose(_shiftTo(M, px)), py))   # coefficient of u^i v^j
    G = [[Fraction(0)] * max(n1, n2) for _ in range(n1 + n2 - 1)]
    for i in range(n1):
        for j in range(n2):
            c = Ps[i][j] * sg ** (i + j)
            if major == 1:
                G[i + j][j] += c
            else:
                G[i + j][i] += c
    k = next((i for i, r in enumerate(G) if any(x != 0 for x in r)), None)
    if k is None:
        return _num([[0]])
    return _num([[float(x) for x in r] for r in G[k:]])


def _shiftTo(P, x0):
    n1 = len(P)
    S = [[0] * len(P[0]) for _ in range(n1)]
    for i in range(n1):
        for k in range(i + 1):
            f = math.comb(i, k) * x0 ** (i - k)
            for j, v in enumerate(P[i]):
                if v:
                    S[k][j] += f * v
    return S


def _chartZ(Z, p, sg, major):
    """Exact blow-up of an integer matrix at a dyadic p (scaled by sc^deg to stay integral)."""
    C = [list(r) for r in Z]
    sc = 1
    while any(v * sc != _round(v * sc) for v in p) and sc < 256:
        sc *= 2
    a = [int(_round(v * sc)) for v in p]
    n1, n2 = biv.shape(C)
    deg = n1 + n2 - 2
    for i in range(n1):
        for j in range(n2):
            e = deg - i - j
            if e > 0:
                C[i][j] *= sc ** e
    for j in range(n2):                             # shift in x, column by column
        for i in range(n1 - 1):
            for k in range(n1 - 2, i - 1, -1):
                C[k][j] += a[0] * C[k + 1][j]
    for i in range(n1):                             # shift in y, row by row
        for j in range(n2 - 1):
            for k in range(n2 - 2, j - 1, -1):
                C[i][k] += a[1] * C[i][k + 1]
    G = biv.zeros(n1 + n2 - 1, max(n1, n2))
    for k in range(n1):
        for l in range(n2):
            c = C[k][l]
            if c == 0:
                continue
            if k + l > 0:
                c *= sc ** (k + l)
            if sg ** (k + l) < 0:
                c = -c
            if major == 1:
                G[k + l][l] += c
            else:
                G[k + l][k] += c
    return G


def _trimS(Z):
    for i in range(len(Z)):
        if any(Z[i]):
            return biv.contentFree(Z[i:])
    return [[0]]


def _blowUpPair(P, p, sg, major):
    if P.get('EZ') is None or not _isDyadic(p):
        return None
    try:
        GU = _chartZ(P['UZ'], p, sg, major)
        GV = _chartZ(P['VZ'], p, sg, major)
        n1 = max(len(GU), len(GV))
        n2 = max(len(GU[0]), len(GV[0]))
        GU = biv.pad(GU, n1, n2)
        GV = biv.pad(GV, n1, n2)
        k = next((i for i in range(n1) if any(GU[i]) or any(GV[i])), None)
        if k is None:
            return _num([[0]])
        GE = _trimS(_chartZ(P['EZ'], p, sg, major))
        GUs = _trimS(_chartZ(P['UsZ'], p, sg, major))
        GVs = _trimS(_chartZ(P['VsZ'], p, sg, major))
        return _makeSurd(GU[k:], GV[k:], P['k'], GE, GUs, GVs)
    except CCA2Error:
        return None


# ---- surds: U + sqrt(k) V ---------------------------------------------------------------------------
def _fmat(M):
    return [[float(x) for x in r] for r in M]


def _makeSurd(UZ, VZ, k, EZ, UsZ, VsZ):
    E = biv.contentFree(EZ)
    if biv.isZero(E):
        E = None                                    # U = -sqrt(k) V is impossible unless both vanish
    rd = math.sqrt(k)
    n1 = max(len(UZ), len(VZ))
    n2 = max(len(UZ[0]), len(VZ[0]))
    U = _fmat(biv.pad(UZ, n1, n2))
    V = _fmat(biv.pad(VZ, n1, n2))
    lo = rd * (1 - 4 * EPS)
    hi = rd * (1 + 4 * EPS)
    return {'t': 'surd', 'U': UZ, 'V': VZ, 'E': E, 'Us': biv.contentFree(UsZ),
            'Vs': biv.contentFree(VsZ),
            'Lo': [[u + lo * v for u, v in zip(ru, rv)] for ru, rv in zip(U, V)],
            'Hi': [[u + hi * v for u, v in zip(ru, rv)] for ru, rv in zip(U, V)],
            'UZ': UZ, 'VZ': VZ, 'EZ': EZ, 'UsZ': UsZ, 'VsZ': VsZ, 'k': k}


def _surdPoly(U, V, k):
    n1 = max(len(U), len(V))
    n2 = max(len(U[0]), len(V[0]))
    U = biv.pad(U, n1, n2)
    V = biv.pad(V, n1, n2)
    g = 0
    for r in U + V:
        for x in r:
            g = math.gcd(g, abs(x))
    if g > 1:
        U = [[x // g for x in r] for r in U]
        V = [[x // g for x in r] for r in V]
    E = biv.padAdd(biv.conv2(U, U), biv.scal(biv.conv2(V, V), -k))
    return _makeSurd(U, V, k, _oddPartZ(E), _oddPartZ(U), _oddPartZ(V))


def _oddPartZ(Z):
    """The product of the factors of ODD multiplicity, with Z's sign (the same sign as Z
    wherever Z is nonzero); Z itself when that cannot be established."""
    D = biv.trimN(Z)
    if biv.isConst(D):
        return [[sign(D[0][0])]]
    try:
        R = biv.contentFree(biv.trimN(_oddPart(D, 1, 0)))
        sg = _relSign(D, R)
        if sg != 0:
            return biv.scal(R, sg)
    except CCA2Error:
        pass
    return Z


_memoSurdOf = {}


def _surdOf(p):
    """p = ((a1 + b1 sqrt(d))/q, (a2 + b2 sqrt(d))/q) recognised from doubles, or None.
    A pure function of the corner, searched by brute force (~5e5 trials when it fails), and the
    same corners come back box after box: remembered."""
    key = (float(p[0]), float(p[1]))
    if key not in _memoSurdOf:
        if len(_memoSurdOf) > 20000:
            _memoSurdOf.clear()
        _memoSurdOf[key] = _surdOfSearch(p)
    r = _memoSurdOf[key]
    return None if r is None else list(r)


def _surdOfSearch(p):
    D = [2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37,
         38, 39, 41, 42, 43, 46, 47]
    bs = list(range(-64, 65))
    for q in range(1, 65):
        for d in D:
            sd = math.sqrt(d)
            ab = []
            ok = True
            for c in range(2):
                best = None
                for b in bs:
                    v = q * p[c] - b * sd
                    a = _round(v)
                    if abs(v - a) < 1e-11 * q * max(1.0, abs(p[c])):
                        if best is None or abs(b) < abs(best[1]):
                            best = (int(a), b)
                if best is None:
                    ok = False
                    break
                ab.append(best)
            if ok and any(x[1] != 0 for x in ab):
                return [ab[0][0], ab[0][1], ab[1][0], ab[1][1], q, d]
    return None


def _surdMul(x, y, d):
    return [x[0] * y[0] + d * x[1] * y[1], x[0] * y[1] + x[1] * y[0]]


def _surdPow(A, n, d):
    M = [[1, 0]]
    for _ in range(n):
        M.append(_surdMul(M[-1], A, d))
    return M


def _blowUpSurd(P, p, sg, major):
    """Blow-up at a corner with QUADRATIC-SURD coordinates: the result is U + sqrt(d) V."""
    r = _surdOf(p)
    if r is None:
        return None
    A = r[0:2]
    B = r[2:4]
    q = r[4]
    d = r[5]
    n1, n2 = biv.shape(P)
    deg = n1 + n2 - 2
    NA = abs(A[0]) + math.sqrt(d) * abs(A[1])
    NB = abs(B[0]) + math.sqrt(d) * abs(B[1])
    bound = 0.0
    for i in range(n1):
        for j in range(n2):
            bound += abs(P[i][j]) * (NA + q) ** i * (NB + q) ** j * float(q) ** (deg - i - j)
    if bound >= 2.0 ** 52:
        return None
    pA = _surdPow(A, n1 - 1, d)
    pB = _surdPow(B, n2 - 1, d)
    Su = biv.zeros(n1, n2)
    Sv = biv.zeros(n1, n2)                          # coefficient of u^k v^l
    for k in range(n1):
        for l in range(n2):
            acc = [0, 0]
            for i in range(k, n1):
                for j in range(l, n2):
                    if P[i][j] == 0:
                        continue
                    m = _surdMul(pA[i - k], pB[j - l], d)
                    f = P[i][j] * math.comb(i, k) * math.comb(j, l) * q ** (k + l + deg - i - j)
                    acc[0] += f * m[0]
                    acc[1] += f * m[1]
            Su[k][l] = acc[0]
            Sv[k][l] = acc[1]
    U = biv.zeros(n1 + n2 - 1, max(n1, n2))
    V = biv.zeros(n1 + n2 - 1, max(n1, n2))
    for k in range(n1):
        for l in range(n2):
            f = sg ** (k + l)
            col = l if major == 1 else k
            U[k + l][col] += f * Su[k][l]
            V[k + l][col] += f * Sv[k][l]
    k0 = next((i for i in range(len(U)) if any(U[i]) or any(V[i])), None)
    if k0 is None:
        return _num([[0]])
    U = U[k0:]
    V = V[k0:]
    E = None                                        # U^2 - d V^2, when it is exact in doubles
    aU = [[abs(x) for x in r] for r in U]
    aV = [[abs(x) for x in r] for r in V]
    bnd = biv.padAdd(biv.conv2(aU, aU), biv.scal(biv.conv2(aV, aV), d))
    if max(x for r in bnd for x in r) < 2 ** 53:
        E = biv.padAdd(biv.conv2(U, U), biv.scal(biv.conv2(V, V), -d))
    rd = math.sqrt(d)
    lo = rd * (1 - 4 * EPS)
    hi = rd * (1 + 4 * EPS)
    Uf = _fmat(U)
    Vf = _fmat(V)
    return {'t': 'surd', 'U': U, 'V': V, 'E': E, 'Us': None, 'Vs': None,
            'Lo': [[u + lo * v for u, v in zip(ru, rv)] for ru, rv in zip(Uf, Vf)],
            'Hi': [[u + hi * v for u, v in zip(ru, rv)] for ru, rv in zip(Uf, Vf)],
            'EZ': None, 'k': d}


# ---- rows comparing two faces ----------------------------------------------------------------------
def _nzOf(FP):
    return max(e[2] for e in FP) + 1


def _faceRowPolys(a, factors):
    Fk = a['FP']
    Fw = a['W']['FP']
    nk = _nzOf(Fk)
    nw = _nzOf(Fw)
    FBk = a['FB']
    FBw = a['W']['FB']
    if nw >= 4:
        if nk not in (2, 3) or FBw != -1:
            return None
        hiK = False
        if nk == 3:
            hiK = FBk in (-1, 2)
            if not hiK and FBk not in (1, -2):
                return None
        elif FBk not in (1, -1):
            return None
        return _topFaceRow(Fk, Fw, nk, hiK, factors)
    if nk == 3 and nw == 2:
        mode = 'ql'
    elif nk == 2 and nw == 3:
        mode = 'lq'
    elif nk == 3 and nw == 3:
        mode = 'qq'
    else:
        return None
    hiK = hiW = False
    if nk == 3:
        hiK = FBk in (-1, 2)
        if not hiK and FBk not in (1, -2):
            return None
    elif FBk not in (1, -1):
        return None
    if nw == 3:
        hiW = FBw in (-1, 2)
        if not hiW and FBw not in (1, -2):
            return None
    elif FBw not in (1, -1):
        return None

    r = _surdFaceRow(Fk, Fw, nk, nw, hiK, hiW, factors)
    if r is not None:
        return r

    key = (mode, _polyKey(Fk), _polyKey(Fw))
    P = _memoFace.get(key)
    if P is None:
        try:
            P = [_squareFreeSign(p) for p in _faceRaw(_slices(Fk), _slices(Fw), -1, mode)]
        except CCA2Error:
            return None
        if len(_memoFace) > 2000:
            _memoFace.clear()
        _memoFace[key] = P
    P = [_num(_divideOut(p, factors)) for p in P]
    if mode == 'qq':
        def test(s):
            return _faceOK(s, hiK, hiW)
    elif mode == 'ql':
        def test(s):                                # z_w >= z_k
            return s[1] != 0 and _quadOK([s[0], s[1] * s[2], s[3], s[4]], hiK, 1)
    else:
        def test(s):                                # z_k <= z_w
            return s[1] != 0 and _quadOK([s[0], s[1] * s[2], s[3], s[4]], hiW, -1)
    return {'P': P, 'test': test}


def _faceRaw(Fk, Fw, sg, mode):
    c2 = biv.conv2
    pa = biv.padAdd
    sc = biv.scal
    if mode != 'qq':
        Q, L = (Fk, Fw) if mode == 'ql' else (Fw, Fk)
        A, B, C = Q[2], Q[1], Q[0]
        Bl, Cl = L[1], L[0]
        G = pa(c2(B, Bl), sc(c2(A, Cl), 2 * sg))
        H = pa(pa(c2(A, c2(Cl, Cl)), sc(c2(B, c2(Cl, Bl)), sg)), c2(C, c2(Bl, Bl)))
        D = pa(c2(B, B), sc(c2(A, C), 4 * sg))
        return [A, Bl, G, H, D]
    Ak, Bk, Ck = Fk[2], Fk[1], Fk[0]
    Aw, Bw, Cw = Fw[2], Fw[1], Fw[0]
    Dw = pa(c2(Bw, Bw), sc(c2(Aw, Cw), 4 * sg))
    Dk = pa(c2(Bk, Bk), sc(c2(Ak, Ck), 4 * sg))
    P1 = pa(c2(Aw, Bk), sc(c2(Ak, Bw), sg))
    E1 = pa(c2(P1, P1), sc(c2(c2(Ak, Ak), Dw), sg))
    P2 = pa(sc(c2(P1, Bw), sg), sc(c2(Aw, pa(c2(Aw, Ck), sc(c2(Ak, Cw), sg))), 2))
    E2 = pa(c2(P2, P2), sc(c2(c2(P1, P1), Dw), sg))
    return [Aw, Ak, Dw, Dk, P1, E1, P2, E2]


def _surdFaceRow(Fk, Fw, nk, nw, hiK, hiW, factors):
    """Both faces written as (a0 + eta sqrt(c) S)/den with a POLYNOMIAL S (the discriminant is a
    constant times a square): their difference is then U + sqrt(kappa) V, a surd polynomial."""
    key = (_polyKey(Fk), _polyKey(Fw), hiK, hiW)
    if key in _memoSurd:
        P = _memoSurd[key]
    else:
        if len(_memoSurd) > 2000:
            _memoSurd.clear()
        _memoSurd[key] = None
        P = None
        try:
            fk = _surdFace(Fk, nk, hiK)
            fw = _surdFace(Fw, nw, hiW)
            if fk is not None and fw is not None and not (fk['lin'] and fw['lin']):
                kap = _squareFreeInt(fw['c'] if not fw['lin'] else fk['c'])
                if kap.bit_length() <= 52:
                    muK = _fieldMu(fk, kap)
                    muW = _fieldMu(fw, kap)
                    if muK is not None and muW is not None:
                        k = kap
                        c2 = biv.conv2
                        U = biv.scal(biv.padAdd(c2(fk['den'], fw['a0']),
                                                biv.scal(c2(fw['den'], fk['a0']), -1)), k)
                        Wk = c2(fw['den'], fk['S'])     # den_w S_k
                        Ww = c2(fk['den'], fw['S'])     # den_k S_w
                        T = [None] * 4
                        for tk in (1, -1):
                            for tw in (1, -1):
                                V = biv.padAdd(biv.scal(Ww, muW * tw * fw['eta']),
                                               biv.scal(Wk, -muK * tk * fk['eta']))
                                T[(2 if tk < 0 else 0) + (1 if tw < 0 else 0)] = _surdPoly(U, V, k)
                        cf = biv.contentFree
                        P = [cf(fk['den']), cf(fw['den']), cf(fk['S']), cf(fw['S'])] + T
                        _memoSurd[key] = P
        except CCA2Error:
            P = None
    if P is None:
        return None
    Pn = [_num(_divideOut(p, factors)) for p in P[0:4]] + \
         [_divideOutSurd(T, factors) for T in P[4:8]]
    return {'P': Pn, 'test': _surdOK}


_memoTop = {}


def _topFaceRow(Fk, Fw, nk, hiK, factors):
    """A FACE ROW AGAINST A TOP ROOT (topFaceRow): z_k > z_w with z_w the LARGEST real root of
    G = Fw, of degree n >= 3 in z, and z_k a face over one quadratic field (surdFace) or rational.
    sign(G_n) G(z_k) < 0 puts a real root of G above z_k, so z_k < z_w and the row fails. With
    z_k = N / M, N = kappa a0 + eta tau mu S sqrt(kappa) = X + Y sqrt(kappa), M = kappa den,
          M^n G(N / M) = sum_i G_i N^i M^(n-i) = U + sqrt(kappa) V,
    by Horner over Z[sqrt kappa]; one (U, V) per tau_k. Sufficient, not equivalent."""
    key = (_polyKey(Fk), _polyKey(Fw), hiK)
    if key in _memoTop:
        P = _memoTop[key]
    else:
        if len(_memoTop) > 2000:
            _memoTop.clear()
        _memoTop[key] = None
        P = None
        try:
            fk = _surdFace(Fk, nk, hiK)
            if fk is not None:
                kap = 1 if fk['lin'] else _squareFreeInt(fk['c'])
                mu = _fieldMu(fk, kap)
                if kap.bit_length() <= 52 and mu is not None:
                    k = kap
                    G = _slices(Fw)
                    n = len(G) - 1
                    c2 = biv.conv2
                    X = biv.scal(fk['a0'], k)
                    M = biv.scal(fk['den'], k)
                    T = [None, None]
                    for tk in (1, -1):
                        Y = biv.scal(fk['S'], mu * tk * fk['eta'])
                        U = G[n]
                        V = [[0]]
                        Mp = [[1]]
                        for i in range(n - 1, -1, -1):
                            Mp = c2(Mp, M)
                            Un = biv.padAdd(c2(U, X), biv.scal(c2(V, Y), k))
                            V = biv.padAdd(c2(U, Y), c2(V, X))
                            U = biv.padAdd(Un, c2(G[i], Mp))
                        T[0 if tk > 0 else 1] = _surdPoly(biv.trimN(U), biv.trimN(V), k)
                    cf = biv.contentFree
                    P = [cf(fk['den']), cf(fk['S']), cf(G[n])] + T + [n % 2]
                    _memoTop[key] = P
        except CCA2Error:
            P = None
    if P is None:
        return None
    odd = P[-1]
    Pn = [_num(_divideOut(p, factors)) for p in P[0:3]] + \
         [_divideOutSurd(T, factors) for T in P[3:5]]
    # NOT blown up at corners (as MATLAB): the Horner images are of high degree with big
    # coefficients. Leaving a row out of a corner chart is sound.
    return {'P': Pn, 'test': (lambda s, odd=odd: _topOK(s, odd)), 'corner': False}


def _topOK(s, odd):
    """s = signs of [den_k, S_k, G_n, T(+), T(-)]: z_k < z_w (the top root) on the box."""
    if any(v == 0 for v in s[0:3]):
        return False
    t = s[3 + (1 if s[0] * s[1] < 0 else 0)]
    dn = s[0] if odd else 1
    return t * s[2] * dn == -1


def _surdOK(s):
    if any(v == 0 for v in s[0:4]):
        return False
    tk = s[0] * s[2]
    tw = s[1] * s[3]
    t = s[4 + (2 if tk < 0 else 0) + (1 if tw < 0 else 0)]
    return t != 0 and t * s[0] * s[1] == 1


def _surdFace(F, n, hi):
    key = (_polyKey(F), hi)
    if key in _memoSurdFace:
        return _memoSurdFace[key]
    f = None
    S = _slices(F)
    if n == 2:
        f = {'lin': True, 'den': S[1], 'a0': biv.scal(S[0], -1), 'S': [[1]], 'c': 0, 'eta': 1}
    else:
        A, B, C = S[2], S[1], S[0]
        D = biv.padAdd(biv.conv2(B, B), biv.scal(biv.conv2(A, C), -4))
        g = biv.content(D)
        if g > 0:
            Dp = [[x // g for x in r] for r in D]
            Sq = _zSqrt(Dp)
            if Sq is not None:
                f = {'lin': False, 'den': biv.scal(A, 2), 'a0': biv.scal(B, -1), 'S': Sq, 'c': g,
                     'eta': 2 * int(hi) - 1}
    if len(_memoSurdFace) > 2000:
        _memoSurdFace.clear()
    _memoSurdFace[key] = f
    return f


def _fieldMu(f, kap):
    if f['lin']:
        return 0
    n = f['c'] * kap
    s = math.isqrt(n)
    return s if s * s == n else None


_PRIMES = None


def _squareFreeInt(c):
    """The square-free part of the positive integer c (trial division below 5000, then a square
    test on what is left)."""
    global _PRIMES
    if _PRIMES is None:
        sieve = [True] * 5001
        _PRIMES = []
        for i in range(2, 5001):
            if sieve[i]:
                _PRIMES.append(i)
                for j in range(i * i, 5001, i):
                    sieve[j] = False
    k = 1
    r = c
    for p in _PRIMES:
        e = 0
        while r % p == 0:
            r //= p
            e += 1
        if e % 2 == 1:
            k *= p
        if r == 1:
            return k
    s = math.isqrt(r)
    if s * s != r:
        k *= r
    return k


def _zSqrt(Dp):
    """A polynomial S with S^2 = Dp exactly, or None."""
    try:
        P = biv.trimN(Dp)
        if biv.isConst(P):
            v = P[0][0]
            if v > 0 and math.isqrt(v) ** 2 == v:
                return [[math.isqrt(v)]]
            return None
        L = _yunParts(P, 1, 0)
        Q = [[1]]
        for fac, mult in L:
            if mult % 2 == 1:
                return None
            for _ in range(mult // 2):
                Q = biv.trimN(biv.conv2(Q, fac))
        Q = biv.contentFree(biv.trimN(Q))
        if biv.equal(biv.trimN(biv.conv2(Q, Q)), P):
            return Q
    except CCA2Error:
        pass
    return None


def _yunParts(S, d, lvl):
    """Yun's square-free decomposition: [(factor, multiplicity)]."""
    L = []
    S = biv.trimN(S)
    if biv.isConst(S) or lvl > 2:
        return L
    Sd = biv.diff(S, d)
    if biv.isZero(Sd):
        if d == 1:
            L = _yunParts(S, 2, lvl + 1)
        return L
    C = biv.trimN(biv.polyGcdBiv(S, Sd))
    W = biv.trimN(biv.polyDivExact(S, C))
    i = 1
    while not biv.isConst(W):
        if i > 60:
            raise CCA2Error('algRowFails:yun', 'square-free decomposition did not end')
        Y = biv.trimN(biv.polyGcdBiv(W, C))
        Z = biv.trimN(biv.polyDivExact(W, Y))
        if not biv.isConst(Z):
            L.append((Z, i))
        C = biv.trimN(biv.polyDivExact(C, Y))
        W = Y
        i += 1
    if d == 1 and not biv.isConst(C):
        L = L + _yunParts(C, 2, lvl + 1)
    return L


# ---- the sign logic --------------------------------------------------------------------------------
def _faceOK(s, hiK, hiW):
    if s[2] == -1 or s[3] == -1:
        return True                                 # a branch exists nowhere here
    if s[0] == 0 or s[1] == 0:
        return False
    e = s[0] if hiW else -s[0]                      # z_w = (-B_w + e sqrt D_w) / (2 A_w)
    u = s[0] * s[1] * _sqSign(s[4], e * s[1], s[5])   # sign(z_w - m_k)
    v = s[1] * _sqSign(s[6], e * s[4], s[7])          # sign(F_k(z_w) / A_k)
    if hiK:
        return u == 1 and v == 1                    # z_w >= the upper root
    return u == 1 or v == -1                        # z_w >= the lower root


def _sqSign(sP, sQ, sE):
    """The sign of P + Q sqrt(.) from the signs of P, Q and E = P^2 - (.) Q^2."""
    if sP != 0 and sP == sQ:
        return sP
    if sE == 1 and sP != 0:
        return sP
    if sE == -1 and sQ != 0:
        return sQ
    return 0


def _quadOK(s, hi, want):
    if s[0] == 0:
        return False
    u = s[1] * s[0]                                 # sign(q - m)
    v = s[2] * s[0]                                 # sign(F(q) / a): -1 is "between the roots"
    if hi and want == -1:
        return u == -1 or v == -1
    if hi:
        return u == 1 and v == 1
    if want == -1:
        return u == -1 and v == 1
    return u == 1 or v == -1


# ---- square-free parts -----------------------------------------------------------------------------
def _squareFreeSign(R):
    """A polynomial with the SAME SIGN as R wherever R is nonzero, without squared factors."""
    S = R
    try:
        S0 = biv.contentFree(biv.trimN(R))
        if biv.isConst(S0):
            return S0
        O = biv.contentFree(biv.trimN(_oddPart(S0, 1, 0)))
        sg = _relSign(S0, O)
        if sg != 0:
            return biv.scal(O, sg)
    except CCA2Error:
        pass
    try:
        S = biv.contentFree(biv.trimN(R))
        for _ in range(4):
            G = [[1]]
            for d in (1, 2):
                Sd = biv.diff(S, d)
                if not biv.isZero(Sd):
                    G = biv.trimN(biv.polyGcdBiv(S, Sd))
                if not biv.isConst(G):
                    break
            if biv.isConst(G):
                return S
            S = biv.contentFree(biv.trimN(biv.polyDivExact(S, biv.conv2(G, G))))
    except CCA2Error:
        pass
    return S


def _oddPart(S, d, lvl):
    """The product of the square-free factors of ODD multiplicity (Yun); through FLINT when
    python-flint is installed (cca2/_poly.py), which gives the same polynomial up to sign."""
    if lvl == 0 and d == 1 and _poly.available():
        O = _poly.oddPart(biv.trimN(S))
        if O is not None:
            return O
    O = [[1]]
    S = biv.trimN(S)
    if biv.isConst(S) or lvl > 2:
        return O
    Sd = biv.diff(S, d)
    if biv.isZero(Sd):
        if d == 1:
            O = _oddPart(S, 2, lvl + 1)
        return O
    C = biv.trimN(biv.polyGcdBiv(S, Sd))
    W = biv.trimN(biv.polyDivExact(S, C))
    i = 1
    while not biv.isConst(W):
        if i > 60:
            raise CCA2Error('algRowFails:yun', 'square-free decomposition did not end')
        Y = biv.trimN(biv.polyGcdBiv(W, C))
        Z = biv.trimN(biv.polyDivExact(W, Y))
        if i % 2 == 1:
            O = biv.trimN(biv.conv2(O, Z))
        C = biv.trimN(biv.polyDivExact(C, Y))
        W = Y
        i += 1
    if d == 1 and not biv.isConst(C):
        O = biv.trimN(biv.conv2(O, _oddPart(C, 2, lvl + 1)))
    return O


_RELPTS = [[0.6180339887, 0.3247179572], [-0.7182818285, 0.4142135624],
           [0.2360679775, -0.8660254038], [1.7320508076, 1.4142135624],
           [-1.2599210499, -0.5772156649], [3.1415926536, -2.7182818285]]


def _relSign(S, O):
    """sign(S) * sign(O) at a generic point (0 when no test point separates them from zero)."""
    for p in _RELPTS:
        vs, ms = _evalTerms(S, p)
        vo, mo = _evalTerms(O, p)
        if abs(vs) > 1e-6 * max(ms, REALMIN) and abs(vo) > 1e-6 * max(mo, REALMIN):
            return sign(vs) * sign(vo)
    return 0


# ---- dividing out the region's own boundary --------------------------------------------------------
def _factorInts(g):
    g = _toIntegers(g)
    if g is None or biv.isConst(g):
        return None
    g = biv.contentFree(biv.trimN(g))
    if biv.isConst(g):
        return None
    return g


def _divideOut(S, factors):
    """S with every factor of the region's boundary polynomials divided out (their sign is
    known on the region, and they vanish exactly where a proof is needed)."""
    for g0 in factors:
        g = _factorInts(g0)
        if g is None:
            continue
        try:
            if not _mayDivide(S, g):
                continue
            for _ in range(4):
                S = biv.trimN(biv.polyDivExact(S, g))
        except CCA2Error:
            pass
    return S


def _divideOutSurd(T, factors):
    U = T['UZ']
    V = T['VZ']
    ch = False
    for g0 in factors:
        g = _factorInts(g0)
        if g is None:
            continue
        try:
            for _ in range(4):
                if not _mayDivide(U, g) or (not biv.isZero(V) and not _mayDivide(V, g)):
                    break
                U2 = biv.trimN(biv.polyDivExact(U, g))
                V2 = biv.trimN(biv.polyDivExact(V, g)) if not biv.isZero(V) else [[0]]
                U, V = U2, V2
                ch = True
        except CCA2Error:
            pass
    E = T['EZ']
    Us = T['UsZ']
    Vs = T['VsZ']
    E2 = _divideOut(E, factors)
    ch = ch or not biv.equal(E2, E)
    E = E2
    U3 = _divideOut(Us, factors)
    ch = ch or not biv.equal(U3, Us)
    Us = U3
    if not biv.isZero(Vs):
        V3 = _divideOut(Vs, factors)
        ch = ch or not biv.equal(V3, Vs)
        Vs = V3
    if not ch:
        return T
    n1 = max(len(U), len(V))
    n2 = max(len(U[0]), len(V[0]))
    return _makeSurd(biv.pad(U, n1, n2), biv.pad(V, n1, n2), T['k'], E, Us, Vs)


def _mayDivide(S, g):
    """A cheap double-precision screen: does S vanish on the curve g = 0? (False only when it
    clearly does not, so an exact division is not attempted.)"""
    if not any(x != 0 for i, r in enumerate(S) for j, x in enumerate(r) if (i, j) != (0, 0)):
        return S[0][0] == 0                         # a constant: only 0 has a non-constant divisor
    if len(g[0]) < 2:                               # g has no y: solve for x instead
        g = biv.transpose(g)
        S = biv.transpose(S)
    if len(g[0]) < 2:
        return True
    for x in (0.6180339887, -1.3247179572, 2.7182818285):
        cy = [sum(g[i][j] * x ** i for i in range(len(g))) for j in range(len(g[0]))]
        while cy and cy[-1] == 0:
            cy.pop()
        if len(cy) < 2:
            continue
        for y in roots(cy[::-1]):
            v, m = _evalTerms(S, [x, y])
            if abs(v) > 1e-6 * max(m, REALMIN):
                return False
    return True

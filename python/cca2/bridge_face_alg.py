"""bridge_face_alg -- the bridge face of the convex envelope for one cell pair whose tie is
CURVED, as an IMPLICIT polynomial (port of bridgeFaceAlg.m).

For dual faces g1, g2 (over the common denominator d) the pair's contribution is
    z(x) = sup over the tie {g1 = g2} of <s, x> - g(s),
attained where the family M_mu = (1+mu) N1 - mu N2, mu in [-1, 0], is stationary. Eliminating mu
(the discriminant of Psi in mu) gives a polynomial FP(x, y, z) with integer coefficients that
vanishes on the face; FB = -1 names the TOP real root in z. The face also carries its region as
cells of conic rows, the pair's data (the eval-time predicate reads it) and a VERIFIED rational
witness. Returns (face, info).

One path on unbounded ints (MATLAB: a double path and a BigInteger path, PORT_NOTES.md 2).
"""
import math

from . import poly_n as pn
from . import rat_q
from .bridge_arc_cut import bridgeArcCut
from .bridge_phi import bridgePhi
from .bridge_witness import bridgeWitness
from .errors import error, CCA2Error
from .numeric import polyval, roots
from .rat_q import sign


def bridgeFaceAlg(n1, d1, n2, d2, P=None):
    P = [list(r) for r in (P or [])]
    if any(len(r) != 3 for r in P):
        error('bridgeFaceAlg:pairShape', 'P must be k x 3 integer half-planes.')
    n1 = _checkFace(n1, 'n1')
    n2 = _checkFace(n2, 'n2')
    d1 = _checkDen(d1, 'd1')
    d2 = _checkDen(d2, 'd2')

    d = rat_q.lcm(d1, d2)
    k1 = d // d1
    k2 = d // d2
    N1 = [[n1[4] * k1, n1[5] * k1], [n1[5] * k1, n1[6] * k1]]
    m1 = [n1[7] * k1, n1[8] * k1]
    q1 = n1[9] * k1
    N2 = [[n2[4] * k2, n2[5] * k2], [n2[5] * k2, n2[6] * k2]]
    m2 = [n2[7] * k2, n2[8] * k2]
    q2 = n2[9] * k2

    def zero(N):
        return all(x == 0 for r in N for x in r)
    r1 = rat_q.detExact(N1) == 0 and not zero(N1)
    r2 = rat_q.detExact(N2) == 0 and not zero(N2)
    if zero(N1) and r2:
        from .bridge_face_rat import bridgeFaceRat
        return bridgeFaceRat(n1, d1, n2, d2, P)
    if zero(N2) and r1:
        from .bridge_face_rat import bridgeFaceRat
        return bridgeFaceRat(n2, d2, n1, d1, P)

    Dl = [[N1[i][j] - N2[i][j] for j in range(2)] for i in range(2)]
    dl = [m1[0] - m2[0], m1[1] - m2[1]]
    kp = q1 - q2                                    # the tie's own data: Delta, delta, kappa
    tie = rat_q.diffConic(n1, d1, n2, d2)           # the tie conic, canonically spelled

    FP, ipInfo = implicitPoly(d, N1, m1, q1, N2, m2, q2)
    C1 = imageConic(tie, N1, m1, d)
    C2 = imageConic(tie, N2, m2, d)

    W, wInfo = bridgeWitness(d, N1, m1, N2, m2, Dl, dl, kp, C1, C2, P)
    FB = None
    lastErr = None
    for wrow in W:
        xn = wrow[0:2]
        xd = wrow[2]
        try:
            s1 = _sideOf(C1, xn, xd)
            s2 = _sideOf(C2, xn, xd)
            FB, bInfo = branchIndex(FP, d, N1, m1, q1, N2, m2, q2, Dl, dl, kp, xn, xd)
            break
        except CCA2Error as ME:
            if ME.identifier not in ('ratQ:overflow', 'bridgeFaceAlg:branchUndecided',
                                     'ratQ:sturmEndpointIsRoot'):
                raise
            lastErr = ME
            FB = None
    if FB is None:
        error('bridgeFaceAlg:branchUndecided',
              'none of the %d verified witnesses could certify the branch index; the last one '
              'reported "%s".' % (len(W), lastErr.message if lastErr else ''))

    base = []
    if C1:
        base.append(C1 + [s1, 0])
    if C2:
        base.append(C2 + [s2, 0])
    cutInfo = {'kind': 'none', 'nCut': 0, 'nCells': 1}
    cells = [base]
    if P:
        try:
            blocks, cutInfo = bridgeArcCut(tie, P, N1, m1, N2, m2, d, [xn[0], xn[1], xd],
                                           C1, C2, kp)
            cells = [base + b for b in blocks]
        except CCA2Error as ME:
            if ME.identifier != 'bridgeArcCut:conicTie':
                raise
            cutInfo = {'kind': 'predicate', 'nCut': 0, 'nCells': 1}

    pair = {'d': d, 'N1': N1, 'm1': m1, 'N2': N2, 'm2': m2, 'Dl': Dl, 'dl': dl, 'kp': kp, 'P': P}
    face = {'kind': 'algebraic', 'FP': FP, 'FB': FB, 'cells': cells, 'pair': pair,
            'wit': [xn[0], xn[1], xd]}
    info = {'tie': tie, 'degrees': ipInfo['degrees'], 'maxCoeff': ipInfo['maxCoeff'],
            'ceiling': ipInfo['maxCoeff'] / 2.0 ** 53, 'muDegree': ipInfo['muDegree'],
            'zContentRemoved': ipInfo['zContentRemoved'],
            'witnessRoots': wInfo['nRoots'], 'witnessTries': wInfo['tries'],
            'value': bInfo['value'], 'nRealZ': bInfo['nRealZ'], 'bracket': bInfo['bracket'],
            'above': bInfo['above'], 'cut': cutInfo['kind'], 'nCut': cutInfo['nCut'],
            'nCells': len(cells)}
    return face, info


# ------------------------------------------------------------------------------------------------
def _muMul(a, b):
    c = [{} for _ in range(len(a) + len(b) - 1)]
    for i, x in enumerate(a):
        if not x:
            continue
        for j, y in enumerate(b):
            if not y:
                continue
            c[i + j] = pn.add(c[i + j], pn.mul(x, y))
    return c


def _muAdd(a, b):
    n = max(len(a), len(b))
    c = [{} for _ in range(n)]
    for i in range(n):
        if i < len(a):
            c[i] = pn.add(c[i], a[i])
        if i < len(b):
            c[i] = pn.add(c[i], b[i])
    return c


def _muScale(a, k):
    return [pn.scale(x, k) for x in a]


def implicitPoly(d, N1, m1, q1, N2, m2, q2):
    """FP(x, y, z): the discriminant in mu of Psi, content and z-content removed."""
    X = pn.mono((1, 0, 0))
    Y = pn.mono((0, 1, 0))
    Z = pn.mono((0, 0, 1))
    U = pn.mono((0, 0, 0))

    def num(v):
        return [pn.scale(U, x) for x in v]
    Dl = [[N1[i][j] - N2[i][j] for j in range(2)] for i in range(2)]
    e11 = [N1[0][0], Dl[0][0]]
    e12 = [N1[0][1], Dl[0][1]]
    e22 = [N1[1][1], Dl[1][1]]

    def conv(a, b):
        out = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return out
    Dt = [x - y for x, y in zip(conv(e11, e22), conv(e12, e12))]   # det M_mu, degree 2 in mu
    a11 = e22
    a12 = [-x for x in e12]
    a22 = e11                                                       # adj(M_mu), symmetric
    c1 = [m1[0], m1[0] - m2[0]]
    c2 = [m1[1], m1[1] - m2[1]]
    kk = [q1, q1 - q2]

    W1 = [pn.add(pn.scale(X, d), pn.scale(U, -c1[0])), pn.scale(U, -c1[1])]   # d x - c_mu
    W2 = [pn.add(pn.scale(Y, d), pn.scale(U, -c2[0])), pn.scale(U, -c2[1])]
    term1 = _muScale(_muMul([pn.add(pn.scale(Z, d), pn.scale(U, kk[0])), pn.scale(U, kk[1])],
                            num(Dt)), 2)
    term2 = _muAdd(_muAdd(_muMul(num(a11), _muMul(W1, W1)),
                          _muScale(_muMul(num(a12), _muMul(W1, W2)), 2)),
                   _muMul(num(a22), _muMul(W2, W2)))
    Psi = _muAdd(term1, _muScale(term2, -1))
    last = 0
    for i, x in enumerate(Psi):
        if x:
            last = i + 1
    Psi = Psi[:max(last, 1)]

    deg = len(Psi) - 1
    M = pn.mul
    if deg == 3:
        A, B, C, E = Psi[3], Psi[2], Psi[1], Psi[0]
        p = pn.add(pn.scale(M(A, M(B, M(C, E))), 18), pn.scale(M(B, M(B, M(B, E))), -4))
        p = pn.add(p, M(M(B, B), M(C, C)))
        p = pn.add(p, pn.scale(M(A, M(C, M(C, C))), -4))
        p = pn.add(p, pn.scale(M(M(A, A), M(E, E)), -27))
    elif deg == 2:
        B, C, E = Psi[2], Psi[1], Psi[0]
        p = pn.add(M(C, C), pn.scale(M(B, E), -4))
    else:
        error('bridgeFaceAlg:straightTie',
              'Psi has degree %d in the multiplier, which happens when the tie is STRAIGHT -- a '
              'double line, a line pair, or no curve at all. And then there is no separate '
              'bridge face to build: a straight tie is an ordinary FACET of each side polygon, '
              'so the sup over that side is the facet or vertex candidate that '
              'envelopeCellRegions already produces once the tie is handed to it as a facet.'
              % deg)

    p = pn.contentFree(p)
    if pn.degrees(p, 3)[2] < 1:
        error('bridgeFaceAlg:noZ',
              'the eliminated polynomial has degree 0 in z, so it states a condition on x alone '
              'and names no value. That is a producer defect, not an input one.')

    # ---- the z-CONTENT: a factor in (x, y) alone common to every z-coefficient ------------------
    Cz = pn.slices(p, 2)
    nz = max(Cz)
    g = pn.contentFree(Cz[nz])
    removed = False
    isUnit = len(g) == 1 and (0, 0) in g and abs(g[(0, 0)]) == 1
    if not isUnit:
        Q = {}
        ok = True
        for k, c in Cz.items():
            try:
                Q[k] = pn.divExact(c, g)
            except CCA2Error as ME:
                if not ME.identifier.startswith('ratQ:'):
                    raise
                ok = False
                break
        if ok:
            p = pn.contentFree(pn.fromSlices(Q, 2))
            removed = True

    info = {'degrees': pn.degrees(p, 3), 'maxCoeff': pn.maxCoeff(p), 'muDegree': deg,
            'zContentRemoved': removed}
    return p, info


def imageConic(tie, Nk, mk, d):
    """The image of the tie conic under s -> grad g_k(s): one boundary of the face's region.
    [] when N_k is singular (the image is not a conic)."""
    if rat_q.detExact(Nk) == 0:
        return []
    adjN = [[Nk[1][1], -Nk[0][1]], [-Nk[1][0], Nk[0][0]]]
    Mn = [[d * x for x in r] for r in adjN]
    tn = [-(adjN[0][0] * mk[0] + adjN[0][1] * mk[1]), -(adjN[1][0] * mk[0] + adjN[1][1] * mk[1])]
    md = rat_q.detExact(Nk)
    q = tie
    n2, _ = rat_q.substAffine([0, 0, 0, 0, 2 * q[0], q[1], 2 * q[2], q[3], q[4], q[5]], 1,
                              Mn, tn, md)
    return rat_q.conic([n2[4], 2 * n2[5], n2[6], 2 * n2[7], 2 * n2[8], 2 * n2[9]])


def fibre(FP, xn, xd):
    """FP restricted to the vertical line over x = xn/xd: integer coefficients in z, highest
    first, content removed."""
    deg = max((e[0] + e[1] for e in FP), default=0)
    nz = max((e[2] for e in FP), default=0)
    acc = [0] * (nz + 1)
    for e, c in FP.items():
        acc[e[2]] += c * xn[0] ** e[0] * xn[1] ** e[1] * xd ** (deg - e[0] - e[1])
    Pz = acc[::-1]
    g = 0
    for v in Pz:
        g = math.gcd(g, abs(v))
    if g > 1:
        Pz = [v // g for v in Pz]
    return Pz, deg


def _dyadic(x, mode):
    from fractions import Fraction
    s = 2 ** 20
    v = math.floor(x * s) if mode == 'floor' else math.ceil(x * s)
    return Fraction(v, s)


def branchIndex(FP, d, N1, m1, q1, N2, m2, q2, Dl, dl, kp, xn, xd):
    """Certify that the face's value at the witness is the TOP real root of FP there.
    Returns (FB, info) with FB = -1."""
    PhiA, Dt, _, _ = bridgePhi(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd)
    Phi = rat_q.polyTrim(PhiA[::-1])
    x = [xn[0] / xd, xn[1] / xd]

    mus = [z.real for z in roots(Phi) if abs(z.imag) <= 1e-10 * max(1.0, abs(z))]
    mus = [m for m in mus if m >= -1 - 1e-12 and m <= 1e-12]
    if mus:
        dScale = max(abs(v) for v in Dt)
        DtR = [float(v) for v in Dt[::-1]]
        mus = [m for m in mus if abs(polyval(DtR, m)) > 1e-9 * max(1.0, dScale)]
    if not mus:
        error('bridgeFaceAlg:noStationary',
              'the witness was VERIFIED to be in the region, so Phi has a root in [-1,0], yet '
              'none survived the double-precision root finder. The exact test and the numeric '
              'one disagree, which is a defect here and not an input.')
    zs = [_psiValue(m, d, N1, m1, q1, N2, m2, q2, x) for m in mus]
    zs = [z for z in zs if math.isfinite(z)]
    if not zs:
        error('bridgeFaceAlg:noStationary',
              'every stationary multiplier at the witness gave a non-finite value.')
    zStar = max(zs)

    Pz, _ = fibre(FP, xn, xd)
    Pz = rat_q.polyTrim(Pz)
    if len(Pz) <= 1:
        error('bridgeFaceAlg:emptyFibre',
              'FP degenerates to a constant at the witness, so it names no value there.')
    rz = [z.real for z in roots(Pz) if abs(z.imag) <= 1e-10 * max(1.0, abs(z))]
    rr = sorted(rz)
    scale = max(1.0, abs(zStar))
    gap = min((abs(r - zStar) for r in rr), default=math.inf)
    if gap > 1e-6 * scale:
        error('bridgeFaceAlg:valueNotARoot',
              'the value %.12g computed from the multiplier is not a root of FP at the witness '
              '(the nearest root is %.3e away). The elimination and the closed form disagree, '
              'which is a defect in this routine and not in its input.' % (zStar, gap))

    w = max(1e-7, 1e-7 * scale)
    lo = _dyadic(zStar - w, 'floor')
    hi = _dyadic(zStar + w, 'ceil')
    M = 1 + -(-max(abs(v) for v in Pz[1:]) // abs(Pz[0])) + 1      # Cauchy: every real root inside
    inBracket = rat_q.countRootsIn(Pz, lo, hi)
    above = rat_q.countRootsIn(Pz, hi, M)
    if inBracket != 1 or above != 0:
        error('bridgeFaceAlg:branchUndecided',
              'the top branch could not be CERTIFIED at the witness: %d root(s) in the bracket '
              '(needs 1) and %d above it (needs 0). The double-precision value and the exact '
              'count disagree, so no branch index is claimed.' % (inBracket, above))
    info = {'value': zStar, 'nRealZ': len(rr), 'bracket': [float(lo), float(hi)], 'above': above}
    return -1, info


def _psiValue(mu, d, N1, m1, q1, N2, m2, q2, x):
    M = [[(1 + mu) * N1[i][j] - mu * N2[i][j] for j in range(2)] for i in range(2)]
    c = [(1 + mu) * m1[i] - mu * m2[i] for i in range(2)]
    k = (1 + mu) * q1 - mu * q2
    dt = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    if dt == 0:
        return math.inf
    adjM = [[M[1][1], -M[0][1]], [-M[1][0], M[0][0]]]
    w = [d * x[0] - c[0], d * x[1] - c[1]]
    wAw = w[0] * (adjM[0][0] * w[0] + adjM[0][1] * w[1]) + w[1] * (adjM[1][0] * w[0] + adjM[1][1] * w[1])
    return (wAw / (2 * dt) - k) / d


def _sideOf(C, xn, xd):
    if not C:
        return 0
    v = (C[0] * xn[0] ** 2 + C[1] * xn[0] * xn[1] + C[2] * xn[1] ** 2
         + C[3] * xn[0] * xd + C[4] * xn[1] * xd + C[5] * xd ** 2)
    return sign(v)


def _checkFace(n, name):
    n = list(n)
    if len(n) != 10:
        error('bridgeFaceAlg:shape', '%s must be a 1x10 face numerator.' % name)
    rat_q.chk(n, name)
    if any(x != 0 for x in n[0:4]):
        error('bridgeFaceAlg:cubicFace',
              '%s has a CUBIC part. A conjugate face is quadratic -- the sup of an affine '
              'function minus a quadratic is one -- so a cubic here means the caller passed '
              'something other than a conjugate cell.' % name)
    return n


def _checkDen(d, name):
    if not isinstance(d, int) or isinstance(d, bool) or d <= 0:
        error('bridgeFaceAlg:denominator', '%s must be a positive integer; got %r.' % (name, d))
    return d

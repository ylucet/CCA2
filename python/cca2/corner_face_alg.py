"""corner_face_alg -- the CORNER faces of one side of a cell pair (port of cornerFaceAlg.m).

Where a facet of the cell meets the tie, the sup over that side is attained AT that dual corner
s for every x in the wedge grad g(s) + cone(facet normal, +-grad tie): an affine face
<s, x> - g(s). The corner is generally a quadratic irrational (a root sigma of the tie restricted
to the facet line), so the face is stored implicitly -- FP(x, y, z) of degree 2 in z with FB
naming the root -- together with its wedge rows. A RATIONAL corner gets a degree-1 FP directly.
Returns a list of face dicts (kind 'corner').
"""
import math
from math import gcd

from . import poly_n as pn
from . import rat_q
from .alg_branch import algBranch
from .corner_wedge_rows import cornerWedgeRows
from .errors import error
from .rat_q import sign, dot, mv, quad


def cornerFaceAlg(n, d, P, tie, side):
    n = list(n)
    if len(n) != 10 or any(x != 0 for x in n[0:4]):
        error('cornerFaceAlg:shape', 'n must be a 1x10 QUADRATIC face numerator.')
    if side not in (-1, 1):
        error('cornerFaceAlg:side', 'side must be +1 or -1; got %r.' % (side,))
    tie = list(tie)
    N = [[n[4], n[5]], [n[5], n[6]]]
    m = [n[7], n[8]]
    qq = n[9]
    faces = []
    for i, row in enumerate(P):
        a = [row[0], row[1]]
        c = row[2]
        na = a[0] * a[0] + a[1] * a[1]
        if na == 0:
            continue
        w = [a[1], -a[0]]
        s0n = [-c * a[0], -c * a[1]]
        al, be, ga = _restrictTie(tie, s0n, w, na)
        if al == 0:
            if be == 0:
                continue
            faces += _rationalCorner(n, d, N, m, qq, P, a, s0n, w, na, be, ga, tie, side, i + 1)
            continue
        g = gcd(gcd(al, be), ga)
        al, be, ga = al // g, be // g, ga // g
        D = be * be - 4 * al * ga
        if D <= 0:
            continue                                # the facet misses the tie, or touches it
        E0, E1, K = _valueForm(N, m, qq, d, s0n, w, na, al, be, ga)
        FP = _implicitValue(E0, E1, K, al, be, D)
        for root in (1, -1):                        # sigma = (-be + root*sqrt(D)) / (2 al)
            sg = (-be + root * math.sqrt(D)) / (2 * al)
            s = [(s0n[0] + sg * w[0]) / na, (s0n[1] + sg * w[1]) / na]
            lim = -1e-9 * max(1.0, abs(s[0]), abs(s[1]))
            if not all(r[0] * s[0] + r[1] * s[1] + r[2] >= lim for r in P):
                continue                            # the corner is not in the cell at all
            gq = [2 * tie[0] * s[0] + tie[1] * s[1] + tie[3],
                  tie[1] * s[0] + 2 * tie[2] * s[1] + tie[4]]
            M = [[a[0], side * gq[0]], [a[1], side * gq[1]]]
            detM = M[0][0] * M[1][1] - M[0][1] * M[1][0]
            if abs(detM) < 1e-12 * max(1.0, max(abs(x) for r in M for x in r)):
                continue                            # the two active constraints are parallel
            p = [(N[0][0] * s[0] + N[0][1] * s[1] + m[0]) / d,
                 (N[1][0] * s[0] + N[1][1] * s[1] + m[1]) / d]   # grad g at the corner
            FB, ok = _branchOf(FP, s, p, d, N, m, qq)
            if not ok:
                continue
            face = {'kind': 'corner', 'FP': FP, 'FB': FB, 's': s,
                    'sigma': {'quad': [al, be, ga], 'root': root, 'facet': i + 1, 's0n': s0n,
                              'w': w, 'na': na},
                    'M': M, 'p': p, 'src': {'n': n, 'd': d, 'tie': tie, 'side': side}, 'con': []}
            face['con'] = cornerWedgeRows(face, n, d, tie, side)[0]
            faces.append(face)
    return faces


def _rationalCorner(n, d, N, m, qq, P, a, s0n, w, na, be, ga, tie, side, i):
    sn = [s0n[0] * be - ga * w[0], s0n[1] * be - ga * w[1]]
    sd = na * be
    if sd < 0:
        sn = [-sn[0], -sn[1]]
        sd = -sd
    gs = gcd(gcd(sn[0], sn[1]), sd)                 # LOWEST TERMS before anything is squared
    sn = [sn[0] // gs, sn[1] // gs]
    sd = sd // gs
    if any(r[0] * sn[0] + r[1] * sn[1] + r[2] * sd < 0 for r in P):
        return []                                   # exact: the corner is not in the cell

    # the corner's affine face, <s, x> - g(s), as a polynomial of degree 1 in z
    FP = {}
    for e, c in (((0, 0, 1), 2 * d * sd * sd), ((1, 0, 0), -2 * d * sd * sn[0]),
                 ((0, 1, 0), -2 * d * sd * sn[1]),
                 ((0, 0, 0), quad(sn, N, sn) + 2 * sd * dot(m, sn) + 2 * sd * sd * qq)):
        if c != 0:
            FP[e] = c
    FP = pn.contentFree(FP)

    T = [[2 * tie[0], tie[1], tie[3]], [tie[1], 2 * tie[2], tie[4]]]
    gqn = mv(T, [sn[0], sn[1], sd])
    Mn = [[a[0], side * gqn[0]], [a[1], side * gqn[1]]]
    dM = Mn[0][0] * Mn[1][1] - Mn[0][1] * Mn[1][0]
    if dM == 0:
        return []                                   # the two constraints are parallel
    pnum = [N[0][0] * sn[0] + N[0][1] * sn[1] + m[0] * sd,
            N[1][0] * sn[0] + N[1][1] * sn[1] + m[1] * sd]
    pd = d * sd                                     # grad g(s) = pnum/pd
    adjM = [[Mn[1][1], -Mn[0][1]], [-Mn[1][0], Mn[0][0]]]
    con = []
    sdm = sign(dM)
    for r in range(2):
        R = [-sdm * pd * adjM[r][0], -sdm * pd * adjM[r][1],
             sdm * (adjM[r][0] * pnum[0] + adjM[r][1] * pnum[1])]
        first = next(x for x in R if x != 0)
        con.append(rat_q.conic([0, 0, 0] + R) + [sign(first), 0])
    s = [sn[0] / sd, sn[1] / sd]
    p = [pnum[0] / pd, pnum[1] / pd]
    M = [[a[0], side * gqn[0] / sd], [a[1], side * gqn[1] / sd]]
    return [{'kind': 'corner', 'FP': FP, 'FB': -1, 's': s,
             'sigma': {'quad': [0, be, ga], 'root': 0, 'facet': i, 's0n': s0n, 'w': w, 'na': na},
             'M': M, 'p': p, 'src': {'n': n, 'd': d, 'tie': tie, 'side': side}, 'con': con}]


def _restrictTie(tie, s0n, w, na):
    """The tie along the facet line s = (s0n + sigma*w)/na: al*sigma^2 + be*sigma + ga."""
    A, B, C, Dc, E, F = tie
    al = A * w[0] ** 2 + B * w[0] * w[1] + C * w[1] ** 2
    be = (2 * A * s0n[0] * w[0] + B * (s0n[0] * w[1] + s0n[1] * w[0]) + 2 * C * s0n[1] * w[1]
          + na * (Dc * w[0] + E * w[1]))
    ga = (A * s0n[0] ** 2 + B * s0n[0] * s0n[1] + C * s0n[1] ** 2
          + na * (Dc * s0n[0] + E * s0n[1]) + F * na * na)
    return al, be, ga


def _valueForm(N, m, qq, d, s0n, w, na, al, be, ga):
    T0 = (al * quad(s0n, N, s0n) - ga * quad(w, N, w) + 2 * na * al * dot(m, s0n)
          + 2 * al * na * na * qq)
    T1 = 2 * al * quad(s0n, N, w) - be * quad(w, N, w) + 2 * na * al * dot(m, w)
    E0 = [2 * al * na * d * s0n[0], 2 * al * na * d * s0n[1], -T0]     # [gx gy const]
    E1 = [2 * al * na * d * w[0], 2 * al * na * d * w[1], -T1]
    K = 2 * al * na * na * d
    return E0, E1, K


def _affPlusZ(e, kz):
    out = {}
    for ex, c in (((1, 0, 0), e[0]), ((0, 1, 0), e[1]), ((0, 0, 0), e[2]), ((0, 0, 1), kz)):
        if c != 0:
            out[ex] = c
    return out


def _implicitValue(E0, E1, K, al, be, D):
    L = _affPlusZ([-2 * al * x + be * y for x, y in zip(E0, E1)], 2 * al * K)
    E = _affPlusZ(E1, 0)
    FP = pn.sub(pn.mul(L, L), pn.scale(pn.mul(E, E), D))
    return pn.contentFree(FP)


def _branchOf(FP, s, p, d, N, m, qq):
    z = s[0] * p[0] + s[1] * p[1] - ((s[0] * (N[0][0] * s[0] + N[0][1] * s[1])
                                      + s[1] * (N[1][0] * s[0] + N[1][1] * s[1])) / 2
                                     + m[0] * s[0] + m[1] * s[1] + qq) / d
    _, info = algBranch(FP, [p], -1, wantInfo=True)
    rr = info['allReal'][0]
    if not rr:
        return 0, False
    gaps = [abs(r - z) for r in rr]
    gap = min(gaps)
    j = gaps.index(gap)
    if gap > 1e-6 * max(1.0, abs(z)):
        return 0, False
    return (j + 1) - len(rr) - 1, True              # negative: counted from the TOP

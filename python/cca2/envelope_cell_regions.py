"""envelope_cell_regions -- the faces of  psi(x) = sup_{s in P} <s,x> - phi(s)  for ONE convex
quadratic phi(s) = (s'Ns/2 + m's + q)/dd on a polyhedral cell P = {P*[s;1] >= 0}
(port of envelopeCellRegions.m).

One face per place the maximiser can sit: the interior ('stationary'), each facet of P along
which phi curves ('facet'), each vertex of P ('vertex'). Each face carries its rows R on x
(the region where that maximiser is the right one), its function fN/fD, and the maximiser map
S: s*(x) = (Mn*x*... + tn)/d. Integer arithmetic throughout.
"""
from . import rat_q
from .errors import error
from .rat_q import sign


def envelopeCellRegions(N, m, q, dd, P):
    if len(N) != 2 or len(N[0]) != 2 or N[0][1] != N[1][0]:
        error('envelopeCellRegions:shape', 'N must be a symmetric 2x2.')
    if dd <= 0:
        error('envelopeCellRegions:denominator', 'dd must be positive; got %s.' % dd)
    m = list(m)
    faces = []
    k = len(P)
    if k == 0:
        return faces

    detN = N[0][0] * N[1][1] - N[0][1] * N[1][0]
    K = [[N[1][1], -N[0][1]], [-N[0][1], N[0][0]]]            # adj(N), so N*K = detN*I
    if detN > 0:
        R = []
        for r in P:
            A = [r[0] * K[0][0] + r[1] * K[1][0], r[0] * K[0][1] + r[1] * K[1][1]]
            R.append([dd * A[0], dd * A[1], -(A[0] * m[0] + A[1] * m[1]) + r[2] * detN])
        if rat_q.feasible2(R, True):
            Km = rat_q.mv(K, m)
            fN = [0] * 10
            fN[4] = 2 * dd ** 2 * K[0][0]
            fN[5] = 2 * dd ** 2 * K[0][1]
            fN[6] = 2 * dd ** 2 * K[1][1]
            fN[7] = -2 * dd * Km[0]
            fN[8] = -2 * dd * Km[1]
            fN[9] = rat_q.dot(m, Km) - 2 * q * detN
            S = _mkMap([[dd * x for x in row] for row in K], [-Km[0], -Km[1]], detN)
            faces.append(_mkFace('stationary', 0, 0, R, fN, 2 * detN * dd, S))

    for i in range(k):
        a = [P[i][0], P[i][1]]
        c = P[i][2]
        na = a[0] * a[0] + a[1] * a[1]
        if na == 0:
            continue
        w = [a[1], -a[0]]                                     # along the facet
        DND = rat_q.quad(w, N, w)
        if DND <= 0:
            continue                                          # phi affine along the facet
        P0 = [-c * a[0], -c * a[1]]                           # P0/na is on the line
        alpha = [na * dd * w[0], na * dd * w[1]]
        beta = -na * rat_q.dot(w, m) - rat_q.quad(w, N, P0)
        R = []
        dead = False
        for j in range(k):
            if j == i:
                continue
            aj = [P[j][0], P[j][1]]
            cj = P[j][2]
            Aj = rat_q.dot(aj, w)
            Bj = rat_q.dot(aj, P0) + cj * na
            if Aj == 0:
                if Bj < 0:
                    dead = True
                    break
                continue
            R.append([Aj * alpha[0], Aj * alpha[1], Aj * beta + Bj * DND])
        if dead:
            continue
        NaV = rat_q.mv(N, a)
        lhsX = [DND * na * dd * a[0], DND * na * dd * a[1]]
        lhsC = -DND * na * rat_q.dot(a, m) - DND * rat_q.dot(a, rat_q.mv(N, P0))
        wN = rat_q.dot(w, NaV)
        R.append([-(lhsX[0] - wN * alpha[0]), -(lhsX[1] - wN * alpha[1]), -(lhsC - wN * beta)])
        if not rat_q.feasible2(R, True):
            continue
        fN = [0] * 10
        fN[4] = 2 * alpha[0] ** 2
        fN[5] = 2 * alpha[0] * alpha[1]
        fN[6] = 2 * alpha[1] ** 2
        fN[7] = 2 * na * DND * dd * P0[0] + 2 * beta * alpha[0]
        fN[8] = 2 * na * DND * dd * P0[1] + 2 * beta * alpha[1]
        fN[9] = (-2 * na * DND * rat_q.dot(P0, m) - DND * rat_q.quad(P0, N, P0)
                 - 2 * na ** 2 * DND * q + beta ** 2)
        S = _mkMap([[w[0] * alpha[0], w[0] * alpha[1]], [w[1] * alpha[0], w[1] * alpha[1]]],
                   [P0[0] * DND + beta * w[0], P0[1] * DND + beta * w[1]], na * DND)
        faces.append(_mkFace('facet', i + 1, 0, R, fN, 2 * na ** 2 * DND * dd, S))

    for i in range(k):
        for j in range(i + 1, k):
            ai = [P[i][0], P[i][1]]
            aj = [P[j][0], P[j][1]]
            D2 = ai[0] * aj[1] - ai[1] * aj[0]
            if D2 == 0:
                continue                                      # parallel rows name no vertex
            rhs = [-P[i][2], -P[j][2]]
            # adj of the row matrix [ai'; aj'] times rhs
            V0 = [aj[1] * rhs[0] - ai[1] * rhs[1], -aj[0] * rhs[0] + ai[0] * rhs[1]]
            vd = D2
            if vd < 0:
                V0 = [-V0[0], -V0[1]]
                vd = -vd
            if not all(r[0] * V0[0] + r[1] * V0[1] + r[2] * vd >= 0 for r in P):
                continue                                      # not a vertex OF P
            # adj of the COLUMN matrix [ai, aj]
            adjA2 = [[aj[1], -aj[0]], [-ai[1], ai[0]]]
            sg = sign(D2)
            MC = [-vd * m[0] - (N[0][0] * V0[0] + N[0][1] * V0[1]),
                  -vd * m[1] - (N[1][0] * V0[0] + N[1][1] * V0[1])]
            R = []
            for row in adjA2:
                R.append([-sg * row[0] * vd * dd, -sg * row[1] * vd * dd,
                          -sg * (row[0] * MC[0] + row[1] * MC[1])])
            if not rat_q.feasible2(R, True):
                continue
            fN = [0] * 10
            fN[7] = 2 * vd * dd * V0[0]
            fN[8] = 2 * vd * dd * V0[1]
            fN[9] = -2 * vd * rat_q.dot(V0, m) - rat_q.quad(V0, N, V0) - 2 * q * vd ** 2
            S = _mkMap([[0, 0], [0, 0]], V0, vd)
            faces.append(_mkFace('vertex', i + 1, j + 1, R, fN, 2 * vd ** 2 * dd, S))
    return faces


def _mkFace(kind, i, j, R, fN, fD, S):
    if fD < 0:
        fN = [-x for x in fN]
        fD = -fD
    fN, fD = rat_q.canon(fN, fD)
    return {'who': {'kind': kind, 'facet': i, 'facet2': j}, 'R': R, 'fN': fN, 'fD': fD, 'S': S}


def _mkMap(Mn, tn, d):
    if d < 0:
        Mn = [[-x for x in r] for r in Mn]
        tn = [-x for x in tn]
        d = -d
    return {'Mn': Mn, 'tn': list(tn), 'd': d}

"""cell_sample_point -- a strictly interior point of a cell with curved rows, by sampling the
polygon of its straight rows (port of cellSamplePoint.m).

C is a table of [conic(6), branch]; rows a list of [index (1-based), side]; H the cell's
straight rows as (normalised) half-planes, R the half-width of the bounding box. Returns
(pt, ok, P): the sample of largest slack, whether one was found, and every sample inside the
cell (capped). Doubles -- a witness search, nothing is decided from it exactly.
"""
import math

import numpy as np

from . import rat_q
from .conic_meet import _solve2
from .numeric import linspace, roundHalfAway


def cellSamplePoint(C, rows, H, R, dense=False):
    pt = [0.0, 0.0]
    ok = False
    P = []
    V = _polygonVertices(H, R)
    if len(V) < 3:
        return pt, ok, P
    c0 = [sum(v[0] for v in V) / len(V), sum(v[1] for v in V) / len(V)]
    V = sorted(V, key=lambda v: math.atan2(v[1] - c0[1], v[0] - c0[0]))
    for n in (12, 48, 192):
        # MATLAB: [I, J] = meshgrid(0:n, 0:n); keep = I + J < n & I > 0 & J > 0 (column-major)
        bi, bj = _barycentric(n)
        qx = []
        qy = []
        for e in range(len(V)):
            a = V[e]
            b = V[(e + 1) % len(V)]
            # the same expression, elementwise: c0 + bi*(a - c0) + bj*(b - c0)
            qx.append(c0[0] + bi * (a[0] - c0[0]) + bj * (b[0] - c0[0]))
            qy.append(c0[1] + bi * (a[1] - c0[1]) + bj * (b[1] - c0[1]))
        Q = np.column_stack((np.concatenate(qx), np.concatenate(qy)))
        s = _slackArr(C, rows, Q[:, 0], Q[:, 1])
        pos = s > 0
        anyPos = bool(pos.any())
        if dense and n < 192 and (n < 48 or not anyPos):
            if anyPos and not ok:
                pt = Q[int(np.argmax(s))].tolist()
                ok = True
            P += Q[pos].tolist()
            continue
        if anyPos or ok:
            if not ok:
                pt = Q[int(np.argmax(s))].tolist()
                ok = True
            P += Q[pos].tolist()
            cap = 4000 if dense else 400
            if len(P) > cap:
                idx = [roundHalfAway(t) for t in linspace(1, len(P), cap)]
                P = [P[i - 1] for i in idx]         # enough for a check
            return pt, ok, P
    return pt, ok, P


_baryMemo = {}


def _barycentric(n):
    """(i/n, j/n) for i, j = 1..n-1 with i + j < n, i outer and j inner (as the list version)."""
    B = _baryMemo.get(n)
    if B is None:
        ij = [(i, j) for i in range(n + 1) for j in range(n + 1) if i + j < n and i > 0 and j > 0]
        B = (np.array([i / n for i, _ in ij]), np.array([j / n for _, j in ij]))
        _baryMemo[n] = B
    return B


def _slack(C, rows, X):
    """_slackArr on a list of points, as a list (kept for callers that hold lists)."""
    if len(X) == 0:
        return []
    A = np.asarray(X, dtype=float)
    return _slackArr(C, rows, A[:, 0], A[:, 1]).tolist()


def _slackArr(C, rows, x, y):
    """min over the rows of side * value / scale, at every point (vectorised; the same
    arithmetic, in the same order, as MATLAB's).

    A value within its own rounding is not a margin: a row counts as satisfied only when its
    value clears 32 eps times the sum of its terms' magnitudes; below that it contributes a slack
    of 0. Otherwise a sample lying exactly ON a line can "strictly" satisfy two rows that
    contradict each other there, and a zero-width cell is kept (cellSamplePoint.m, 2026-11-11)."""
    xx = x * x                                      # the shared products, formed once (the same
    xy = x * y                                      # doubles each row formed for itself)
    yy = y * y
    nx = np.maximum(1.0, np.sqrt(xx + yy))
    nx2 = None
    s = np.full(len(x), np.inf)
    eps32 = 32 * 2.0 ** -52
    for r in rows:
        c = C[r[0] - 1][0:6]
        br = C[r[0] - 1][6]
        sd = r[1]
        if br != 0:
            v, m = rat_q.pairLineEvalMArr(c, br, x, y)
            sc = nx
        else:
            v = c[0] * xx + c[1] * xy + c[2] * yy + c[3] * x + c[4] * y + c[5]
            m = (abs(c[0]) * xx + np.abs(c[1] * xy) + abs(c[2]) * yy
                 + np.abs(c[3] * x) + np.abs(c[4] * y) + abs(c[5]))
            if nx2 is None:
                nx2 = nx * nx
            sc = max(abs(t) for t in c) * nx2
        w = sd * v
        w = np.where((w > 0) & (w <= eps32 * m), 0.0, w)
        s = np.minimum(s, w / sc)
    return s


def _polygonVertices(H, R):
    A = [list(h) for h in H] + [[1, 0, R], [-1, 0, R], [0, 1, R], [0, -1, R]]
    V = []
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            M = [[A[i][0], A[i][1]], [A[j][0], A[j][1]]]
            if abs(M[0][0] * M[1][1] - M[0][1] * M[1][0]) < 1e-12:
                continue
            x = _solve2(M, [-A[i][2], -A[j][2]])
            lim = -1e-9 * max(1.0, math.sqrt(x[0] * x[0] + x[1] * x[1]))
            if all(a[0] * x[0] + a[1] * x[1] + a[2] >= lim for a in A):
                V.append(x)
    if not V:
        return V
    return _uniquetolRows(V, 1e-12)


def _uniquetolRows(V, tol):
    """MATLAB uniquetol(V, tol, 'ByRows', true, 'DataScale', 1): rows within tol of one another
    are merged; the result is sorted, each group represented by its first (smallest) row."""
    S = sorted(V)
    out = []
    for v in S:
        if any(abs(v[0] - u[0]) <= tol and abs(v[1] - u[1]) <= tol for u in out):
            continue
        out.append(v)
    return out


def cellWitness(C, rows):
    """A strictly interior point of the cell by a SECTOR search, or (pt, False).

    cellSamplePoint grids the polygon of the straight rows clipped to a box of half-width R, so a
    cell that is small next to that polygon (a unit-wide cell under rows that leave the polygon
    open to |x| = 1e3) falls between the grid points, and a live cell was dropped: co f = +inf
    at points of dom f (review 2026-11-19 #9c). This search does not depend on any box. It is a
    cylindrical decomposition of the rows' curves (a branch row as its own line): the critical
    abscissae are the roots of every curve's discriminant in y, of its leading coefficient in y,
    and of every pairwise resultant in y; an open cell meets the open sector between two
    consecutive critical abscissae and two consecutive curves above one of them, so the middle
    of every sector (and past both ends) is tested. Doubles, with the slack rule of
    cellSamplePoint (a value inside its own rounding is not a margin): a witness search, the
    point is checked, never trusted."""
    if not rows:
        return [0.0, 0.0], True
    P = []
    for r in rows:
        c = C[r[0] - 1][0:6]
        br = C[r[0] - 1][6]
        if br != 0:
            v, _ = rat_q.pairLineEvalM(c, br, [[0, 0], [1, 0], [0, 1]])
            P.append([0.0, 0.0, 0.0, v[1] - v[0], v[2] - v[0], v[0]])
        else:
            P.append([float(t) for t in c])
    P = [list(p) for p in sorted({tuple(p) for p in P})]       # deterministic order
    # each curve as a quadratic in y with polynomial coefficients in x (highest power first)
    Y = [([p[2]], [p[1], p[4]], [p[0], p[3], p[5]]) for p in P]
    xs = []
    for A, B, Cc in Y:
        if A[0] != 0:
            xs += _realRootsF(np.polysub(np.polymul(B, B), np.polymul([4 * A[0]], Cc)))
        else:
            xs += _realRootsF(B)
    for i in range(len(Y)):
        for j in range(i + 1, len(Y)):
            xs += _realRootsF(_resultantY(Y[i], Y[j]))
    xs = sorted(set(x for x in xs if math.isfinite(x)))
    if not xs:
        xs = [0.0]
    span = max(1.0, xs[-1] - xs[0])
    samp = [xs[0] - span, xs[-1] + span] + [(xs[i] + xs[i + 1]) / 2 for i in range(len(xs) - 1)]
    best = -math.inf
    pt = [0.0, 0.0]
    for x in samp:
        ys = []
        for A, B, Cc in Y:
            a = A[0]
            b = float(np.polyval(B, x))
            c = float(np.polyval(Cc, x))
            if a != 0:
                d = b * b - 4 * a * c
                if d >= 0:
                    sd = math.sqrt(d)
                    ys += [(-b - sd) / (2 * a), (-b + sd) / (2 * a)]
            elif b != 0:
                ys.append(-c / b)
        ys = sorted(set(y for y in ys if math.isfinite(y)))
        if not ys:
            ys = [0.0]
        sp = max(1.0, ys[-1] - ys[0])
        cand = [ys[0] - sp, ys[-1] + sp] + [(ys[i] + ys[i + 1]) / 2 for i in range(len(ys) - 1)]
        cy = np.array(cand)
        s = _slackArr(C, rows, np.full(cy.shape, x), cy)
        k = int(np.argmax(s))
        if s[k] > best:
            best = float(s[k])
            pt = [float(x), float(cy[k])]
    return pt, best > 0


def _realRootsF(p):
    p = list(np.trim_zeros(np.asarray(p, dtype=float), 'f'))
    if len(p) < 2:
        return []
    return [z.real for z in np.roots(p) if abs(z.imag) <= 1e-9 * max(1.0, abs(z))]


def _resultantY(p, q):
    """The resultant in y of two curves A y^2 + B y + C (A constant, B, C polynomials in x)."""
    A, B, Cc = p
    D, E, F = q
    a, d = A[0], D[0]
    if a == 0 and d == 0:                       # two lines in y: B y + C = 0, E y + F = 0
        return np.polysub(np.polymul(B, F), np.polymul(Cc, E))
    if a == 0:                                  # y = -C/B into the other: D C^2 - E B C + F B^2
        return np.polyadd(np.polysub(np.polymul([d], np.polymul(Cc, Cc)),
                                     np.polymul(E, np.polymul(B, Cc))),
                          np.polymul(F, np.polymul(B, B)))
    if d == 0:
        return _resultantY(q, p)
    # Sylvester, degree 2 x 2: (AF - CD)^2 - (AE - BD)(BF - CE)
    t1 = np.polysub(np.polymul([a], F), np.polymul(Cc, [d]))
    t2 = np.polysub(np.polymul([a], E), np.polymul(B, [d]))
    t3 = np.polysub(np.polymul(B, F), np.polymul(Cc, E))
    return np.polysub(np.polymul(t1, t1), np.polymul(t2, t3))

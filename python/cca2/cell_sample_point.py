"""cell_sample_point -- a strictly interior point of a cell with curved rows, by sampling the
polygon of its straight rows (port of cellSamplePoint.m).

C is a table of [conic(6), branch]; rows a list of [index (1-based), side]; H the cell's
straight rows as (normalised) half-planes, R the half-width of the bounding box. Returns
(pt, ok, P): the sample of largest slack, whether one was found, and every sample inside the
cell (capped). Doubles -- a witness search, nothing is decided from it exactly.
"""
import math

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
        B = [(i / n, j / n) for i in range(n + 1) for j in range(n + 1)
             if i + j < n and i > 0 and j > 0]
        Q = []
        for e in range(len(V)):
            a = V[e]
            b = V[(e + 1) % len(V)]
            for bi, bj in B:
                Q.append([c0[0] + bi * (a[0] - c0[0]) + bj * (b[0] - c0[0]),
                          c0[1] + bi * (a[1] - c0[1]) + bj * (b[1] - c0[1])])
        s = _slack(C, rows, Q)
        anyPos = any(v > 0 for v in s)
        if dense and n < 192 and (n < 48 or not anyPos):
            if anyPos and not ok:
                i = s.index(max(s))
                pt = Q[i]
                ok = True
            P += [q for q, v in zip(Q, s) if v > 0]
            continue
        if anyPos or ok:
            if not ok:
                i = s.index(max(s))
                pt = Q[i]
                ok = True
            P += [q for q, v in zip(Q, s) if v > 0]
            cap = 4000 if dense else 400
            if len(P) > cap:
                idx = [roundHalfAway(t) for t in linspace(1, len(P), cap)]
                P = [P[i - 1] for i in idx]         # enough for a check
            return pt, ok, P
    return pt, ok, P


def _slack(C, rows, X):
    s = [math.inf] * len(X)
    nx = [max(1.0, math.sqrt(x[0] * x[0] + x[1] * x[1])) for x in X]
    for r in rows:
        c = C[r[0] - 1][0:6]
        br = C[r[0] - 1][6]
        sd = r[1]
        if br != 0:
            v = rat_q.pairLineEval(c, br, X)
            sc = nx
        else:
            v = [c[0] * x[0] ** 2 + c[1] * x[0] * x[1] + c[2] * x[1] ** 2
                 + c[3] * x[0] + c[4] * x[1] + c[5] for x in X]
            m = max(abs(t) for t in c)
            sc = [m * t * t for t in nx]
        s = [min(a, sd * b / d) for a, b, d in zip(s, v, sc)]
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

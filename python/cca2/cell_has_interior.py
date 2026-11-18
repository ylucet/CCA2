"""cell_has_interior -- does a cell cut out by conic sign conditions have interior?
(port of cellHasInterior.m)

rows is a list of [index (1-based) into EcQ, side]. Lines are decided exactly by
rat_q.feasible2; a curved cell is decided by looking for a witness point (pairwise conic
intersections, then a seeded random probe) -- the same screen MATLAB uses.
"""
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
    """Complete sample search (a 2-D cylindrical decomposition), run when the witness search finds
    nothing: the cell is open, so if non-empty it meets a sector of the decomposition by the critical
    abscissae (conic discriminant roots, vertical asymptotes, pairwise intersections); a point in the
    middle of every sector (and past both ends) is tested. The random probe box is only a few units
    wide and misses a cell far from the origin when the conics do not meet (DECISIONS 2026-11-17).
    Skipped when a row is an equality (side 0), which has no interior sector."""
    if any(r[1] == 0 for r in rows):
        return False
    idx = sorted(set(r[0] for r in rows))
    cs = [[float(t) for t in EcQ[i - 1][0:6]] for i in idx]
    xs = []
    for c in cs:
        a0, a1, a2, a3, a4, a5 = c
        if a2 != 0:   # disc(x) = (a1 x + a4)^2 - 4 a2 (a0 x^2 + a3 x + a5)
            q = [a1 * a1 - 4 * a2 * a0, 2 * a1 * a4 - 4 * a2 * a3, a4 * a4 - 4 * a2 * a5]
        elif a1 != 0:
            xs.append(-a4 / a1)
            continue
        else:
            continue
        if q[0] != 0:
            d = q[1] ** 2 - 4 * q[0] * q[2]
            if d >= 0:
                sd = d ** 0.5
                xs += [(-q[1] - sd) / (2 * q[0]), (-q[1] + sd) / (2 * q[0])]
        elif q[1] != 0:
            xs.append(-q[2] / q[1])
    for a in range(len(idx)):
        for b in range(a + 1, len(idx)):
            try:
                P, _ = conicMeet(EcQ[idx[a] - 1], EcQ[idx[b] - 1])
            except CCA2Error as ME:
                if ME.identifier != 'ratQ:overflow':
                    raise
                continue
            xs += [float(p[0]) for p in P]
    xs = sorted(set(x for x in xs if np.isfinite(x)))
    if not xs:
        xs = [0.0]
    samp = [xs[0] - 1.0, xs[-1] + 1.0] + [(xs[i] + xs[i + 1]) / 2 for i in range(len(xs) - 1)]
    for x in samp:
        ys = []
        for a0, a1, a2, a3, a4, a5 in cs:
            B = a1 * x + a4
            C = a0 * x * x + a3 * x + a5
            if a2 != 0:
                d = B * B - 4 * a2 * C
                if d >= 0:
                    sd = d ** 0.5
                    ys += [(-B - sd) / (2 * a2), (-B + sd) / (2 * a2)]
            elif B != 0:
                ys.append(-C / B)
        ys = sorted(set(y for y in ys if np.isfinite(y)))
        if not ys:
            ys = [0.0]
        cand = [ys[0] - 1.0, ys[-1] + 1.0] + [(ys[i] + ys[i + 1]) / 2 for i in range(len(ys) - 1)]
        X1 = np.array(cand)
        X0 = np.full(X1.shape, x)
        if _satisfiesAllMany(EcQ, rows, X0, X1).any():
            return True
    return False


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

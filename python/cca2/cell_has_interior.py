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
    return _probeInterior(EcQ, rows)


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
        _PROBE = matlabRand(20260905, 400, 2)
    for u in _PROBE:
        x = [ctr[0] + (u[0] - 0.5) * 2 * r, ctr[1] + (u[1] - 0.5) * 2 * r]
        if _satisfiesAll(EcQ, rows, x):
            return True
    return False

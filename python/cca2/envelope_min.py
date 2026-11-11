"""EnvelopeMin -- the convex envelope as a MIN over atom pairs, each pair's envelope a MAX of
candidates, plus the contact TRIANGLES (port of EnvelopeMin.m).

It is the AlgCon that carries no subdivision: `biconj` returns it where the subdivision is
declined, and it can only be evaluated. eval returns the list of values.
"""
import math

import numpy as np

from .envelope_eval import envelopeEval
from .qua_con import AlgCon


class EnvelopeMin(AlgCon):
    def __init__(self, pairs, tris, atoms):
        self.pairs = list(pairs)      # envelopeAssemble results
        self.tris = list(tris)        # dicts: T (3x2 contact triangle), s (dual point), c (g(s)),
        #                               optional G (generators of the normal cone of dom f* at s)
        self.atoms = list(atoms)      # the atoms, for inspection
        self.nf = None

    def isSubdivided(self):
        return False

    def eval(self, X):
        n = len(X)
        val = [math.inf] * n
        for pr in self.pairs:
            v, _ = envelopeEval(pr, X)
            val = [min(a, b) for a, b in zip(val, v)]
        for t in self.tris:
            T = t['T']
            s = t['s']
            G = t.get('G') or []
            for i in range(n):
                x = [float(X[i][0]), float(X[i][1])]
                inside = _inTrianglePlusCone(x, T, G) if G else _inTriangle(x, T)
                if inside:
                    val[i] = min(val[i], x[0] * s[0] + x[1] * s[1] - t['c'])
        return val

    def __repr__(self):
        return 'EnvelopeMin(%d pairs, %d triangles)' % (len(self.pairs), len(self.tris))


def _inTriangle(x, T):
    c = []
    for i in range(3):
        a = T[i]
        b = T[(i + 1) % 3]
        c.append((b[0] - a[0]) * (x[1] - a[1]) - (b[1] - a[1]) * (x[0] - a[0]))
    tol = 1e-9 * max(1.0, max(abs(v) for r in T for v in r)) ** 2
    return all(v >= -tol for v in c) or all(v <= tol for v in c)


def _inTrianglePlusCone(x, T, G):
    """Is x in T + cone(G)? G is a list of generators [g1, g2] (erratum Case E)."""
    k = len(G)
    tol = 1e-9 * max(1.0, max(abs(v) for r in T for v in r)) ** 2
    A = []
    b = []
    for e in range(3):
        a = T[e]
        c = T[(e + 1) % 3]
        o = T[(e + 2) % 3]
        nrm = [-(c[1] - a[1]), c[0] - a[0]]          # left normal of edge a -> c
        sg = np.sign(nrm[0] * (o[0] - a[0]) + nrm[1] * (o[1] - a[1]))
        if sg == 0:
            sg = 1
        A.append([sg * (nrm[0] * g[0] + nrm[1] * g[1]) for g in G])
        b.append(sg * (nrm[0] * (x[0] - a[0]) + nrm[1] * (x[1] - a[1])))
    for j in range(k):
        A.append([-1.0 if i == j else 0.0 for i in range(k)])
        b.append(0.0)
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    cand = []
    if k == 1:
        cand = [np.array([bi / ai]) for ai, bi in zip(A[:, 0], b) if ai != 0]
        cand.append(np.array([0.0]))
    else:
        for p in range(A.shape[0]):
            for q in range(p + 1, A.shape[0]):
                M = A[[p, q], :]
                if M.shape[1] != 2 or abs(np.linalg.det(M)) < 1e-14:
                    continue
                cand.append(np.linalg.solve(M, b[[p, q]]))
    for c in cand:
        if np.all(A @ c <= b + tol):
            return True
    return False

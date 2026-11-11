"""domain_is_convex_d -- is the DOMAIN of a polyhedral subdivision convex? (port of
domainIsConvexD.m). Doubles with a tolerance: it feeds only QuaPol.dom['isConvex'].
"""
import math

from .errors import error


def domainIsConvexD(V, E, F, tol=None):
    if tol is None:
        tol = math.sqrt(2.0 ** -52)
    if not F or not E or not V:
        return True
    if all(x == 0 for r in F for x in r):
        error('PLQ:domainIsConvexD:thinSubdivision',
              'every edge has +infinity on both sides, so the domain has dimension < 2 and this '
              'test says nothing about it.')
    Vf = [[float(x) for x in r] for r in V]
    dirs = []
    for e in E:
        if e[2] == 0:
            u = [Vf[e[1] - 1][0] - Vf[e[0] - 1][0], Vf[e[1] - 1][1] - Vf[e[0] - 1][1]]
            n = math.hypot(u[0], u[1])
            if n > 0:
                dirs.append([u[0] / n, u[1] / n])
    for j in range(len(F)):
        iL, iR = F[j]
        if (iL != 0 and iR != 0) or (iL == 0 and iR == 0):
            continue
        v1 = Vf[E[j][0] - 1]
        v2 = Vf[E[j][1] - 1]
        a = [-(v2[1] - v1[1]), v2[0] - v1[0]]       # the LEFT normal of the edge direction
        na = math.hypot(a[0], a[1])
        if na == 0:
            continue
        c = v1[0] * v2[1] - v2[0] * v1[1]
        s = -1.0 if iL == 0 else 1.0                # the finite side is the RIGHT one
        a = [s * a[0] / na, s * a[1] / na]
        c = s * c / na
        for p in Vf:
            if a[0] * p[0] + a[1] * p[1] + c < -tol * max(1.0, math.hypot(p[0], p[1])):
                return False
        for d in dirs:
            if a[0] * d[0] + a[1] * d[1] < -tol:
                return False
    return True

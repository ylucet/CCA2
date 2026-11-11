"""is_convex_q -- is an exact QuaPol convex? Decided exactly (port of isConvexQ.m).

Three conditions: every face's Hessian is positive semidefinite (Sylvester on integers); the
domain is convex; and across every interior edge f is continuous with a gradient jump pointing
the right way at both ends (and not decreasing along a ray). Returns (tf, why).
"""
from . import rat_q
from .errors import error
from .qua_pol import QuaPol


def isConvexQ(obj):
    if not isinstance(obj, QuaPol):
        error('PLQ:isConvexQ:input',
              'isConvexQ takes a QuaPol (quadratic on polyhedral); got a %s.'
              % type(obj).__name__)
    obj.assertOperable()
    obj.assertExact()
    nf = len(obj.fN)
    for k in range(nf):
        H, _ = rat_q.hessQ(obj.fN[k], obj.fD[k])
        if not rat_q.isPSD2(H):
            return False, 'face %d has a Hessian that is not positive semidefinite' % (k + 1)
    if obj.nv == 0:
        return True, ''
    if obj.dom['dim'] < 2:
        error('PLQ:isConvexQ:thinDomain',
              'isConvexQ is for a 2-dimensional domain; this one has dimension %d. biconjQ '
              'routes thin domains to thinEnvelope before asking.' % obj.dom['dim'])
    ok, jBad = _domainIsConvex(obj)
    if not ok:
        return False, ('the domain is not convex: it leaves the half-plane of boundary edge %d'
                       % jBad)
    for j in range(len(obj.F)):
        iL, iR = obj.F[j]
        if iL == 0 or iR == 0 or iL == iR:
            continue
        dn, dd = rat_q.sub(obj.fN[iL - 1], obj.fD[iL - 1], obj.fN[iR - 1], obj.fD[iR - 1])
        P1, D1 = _vertexOf(obj, obj.E[j][0])
        P2, D2 = _vertexOf(obj, obj.E[j][1])
        a = [-(P2[1] * D1 - P1[1] * D2), P2[0] * D1 - P1[0] * D2]
        if a[0] == 0 and a[1] == 0:
            error('PLQ:isConvexQ:degenerateEdge',
                  'edge %d has both endpoints at the same vertex, so it has no normal.' % (j + 1))
        Pm = [P1[0] * D2 + P2[0] * D1, P1[1] * D2 + P2[1] * D1]
        Dm = 2 * D1 * D2
        if _evalNum(dn, P1, D1) != 0 or _evalNum(dn, P2, D2) != 0 or _evalNum(dn, Pm, Dm) != 0:
            return False, 'f is discontinuous across edge %d' % (j + 1)
        G, _ = rat_q.gradQ(dn, dd)                # grad d (x,y) = G*[x;y;1] / dd, with dd > 0
        s1 = _deltaNum(G, a, P1, D1)
        s2 = _deltaNum(G, a, P2, D2)
        if s1 < 0 or s2 < 0:
            return False, ('the subdifferential is empty on edge %d: the gradient jump points '
                           'the wrong way at an endpoint' % (j + 1))
        if obj.E[j][2] == 0:
            if s1 * D2 > s2 * D1:
                return False, ('the subdifferential empties out along the unbounded edge %d: '
                               'delta decreases' % (j + 1))
    return True, ''


def _vertexOf(obj, i):
    P = obj.VN[i - 1]
    D = obj.VD[i - 1]
    if D <= 0:
        error('PLQ:isConvexQ:badVertex', 'vertex %d has a nonpositive denominator.' % i)
    return P, D


def _evalNum(n, P, D):
    return (n[4] * P[0] ** 2 + 2 * n[5] * P[0] * P[1] + n[6] * P[1] ** 2
            + 2 * n[7] * P[0] * D + 2 * n[8] * P[1] * D + 2 * n[9] * D ** 2)


def _deltaNum(G, a, P, D):
    g = rat_q.mv(G, [P[0], P[1], D])
    return g[0] * a[0] + g[1] * a[1]


def _domainIsConvex(obj):
    rays = []
    for e in obj.E:
        if e[2] == 0:
            Q1, C1 = _vertexOf(obj, e[0])
            Q2, C2 = _vertexOf(obj, e[1])
            rays.append([Q2[0] * C1 - Q1[0] * C2, Q2[1] * C1 - Q1[1] * C2])
    for j in range(len(obj.F)):
        iL, iR = obj.F[j]
        if iL != 0 and iR != 0:
            continue                               # an interior edge is not a domain constraint
        if iL == 0 and iR == 0:
            error('PLQ:isConvexQ:danglingEdge',
                  'edge %d has +infinity on both sides, so the subdivision is not a subdivision.'
                  % (j + 1))
        P1, D1 = _vertexOf(obj, obj.E[j][0])
        P2, D2 = _vertexOf(obj, obj.E[j][1])
        a = [-(P2[1] * D1 - P1[1] * D2), P2[0] * D1 - P1[0] * D2]
        c = P1[0] * P2[1] - P2[0] * P1[1]
        sgn = -1 if iL == 0 else 1                 # the finite side is the RIGHT one
        for i in range(1, len(obj.VN) + 1):
            Q, C = _vertexOf(obj, i)
            if sgn * (a[0] * Q[0] + a[1] * Q[1] + c * C) < 0:
                return False, j + 1
        for r in rays:
            if sgn * (a[0] * r[0] + a[1] * r[1]) < 0:
                return False, j + 1
    return True, 0

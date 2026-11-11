"""triangle_face -- the contact TRIANGLE face of the envelope, exactly (port of triangleFace.m
and triangleFaceSign.m).

Three dual forms tie at one isolated dual point s0; their contact points grad g_i(s0) span a
triangle T on which co f is the affine function <s0, x> - g(s0). s0 is rational when both tie
conics are lines ('rational': corners Tn/Td exact), otherwise it is NAMED by dualPoint
('algebraic'), and every sign needed about it is decided exactly by dualPointSign.
"""
from . import biv
from . import rat_q
from .dual_point import dualPoint, dualPointSign
from .errors import error
from .rat_q import sign


def triangleFace(tri):
    if not isinstance(tri, dict) or not all(k in tri for k in ('Gn', 'Gd', 's')):
        error('triangleFace:shape',
              'expects a biconjAtoms(f, \'triangles\') entry with .Gn .Gd .s.')
    if tri.get('G'):
        error('triangleFace:unbounded',
              'the face holds on T + cone(G) (erratum Case E); only the bounded triangle is '
              'carried so far.')
    Gn = [list(r) for r in tri['Gn']]
    Gd = list(tri['Gd'])
    if len(Gn) != 3 or any(len(r) != 6 for r in Gn) or len(Gd) != 3:
        error('triangleFace:noExactForms',
              'the forms have no exact conic, so s0 cannot be named.')
    t12 = rat_q.conic([a * Gd[1] - b * Gd[0] for a, b in zip(Gn[0], Gn[1])])
    t13 = rat_q.conic([a * Gd[2] - b * Gd[0] for a, b in zip(Gn[0], Gn[2])])
    tf = {'kind': '', 's0n': None, 's0d': None, 'Tn': None, 'Td': None, 'dp': None, 'Gn': Gn,
          'Gd': Gd, 's': None, 'c': None, 'T': None, 'orient': 0}
    if all(x == 0 for x in t12[0:3]) and all(x == 0 for x in t13[0:3]):
        xn, xd = rat_q.solve2([t12[3:5], t13[3:5]], [-t12[5], -t13[5]])
        if xd < 0:
            xn = [-v for v in xn]
            xd = -xd
        tf['kind'] = 'rational'
        tf['s0n'] = xn
        tf['s0d'] = xd
        tf['s'] = [xn[0] / xd, xn[1] / xd]
        tf['Tn'] = []
        tf['Td'] = []
        for i in range(3):
            g = Gn[i]
            H = [[2 * g[0], g[1]], [g[1], 2 * g[2]]]
            tn = [H[0][0] * xn[0] + H[0][1] * xn[1] + g[3] * xd,
                  H[1][0] * xn[0] + H[1][1] * xn[1] + g[4] * xd]
            tn, td = rat_q.canon(tn, xd * Gd[i])
            tf['Tn'].append(tn)
            tf['Td'].append(td)
    else:
        tf['kind'] = 'algebraic'
        tf['dp'] = dualPoint(t12, t13, tri['s'])
        tf['s'] = tf['dp']['s']
    s = tf['s']
    G = [v / Gd[0] for v in Gn[0]]
    tf['c'] = G[0] * s[0] ** 2 + G[1] * s[0] * s[1] + G[2] * s[1] ** 2 + G[3] * s[0] + G[4] * s[1] + G[5]
    tf['T'] = []
    for i in range(3):
        G = [v / Gd[i] for v in Gn[i]]
        tf['T'].append([2 * G[0] * s[0] + G[1] * s[1] + G[3], G[1] * s[0] + 2 * G[2] * s[1] + G[4]])
    tf['orient'] = triangleFaceSign(tf, areaPoly(Gn, Gd))
    if tf['orient'] == 0:
        error('triangleFace:flat',
              'the contact triangle has zero area: its corners are collinear.')
    return tf


def contactPoly(g):
    """grad of the form g = [a b c d e f] as two coefficient matrices in (s1, s2)."""
    return [[[g[3], g[1]], [2 * g[0], 0]], [[g[4], 2 * g[2]], [g[1], 0]]]


def areaPoly(Gn, Gd):
    """Twice the signed area of the contact triangle, as a polynomial in s (cleared of the
    denominators)."""
    M = [contactPoly(g) for g in Gn]
    A = [[0]]
    for i, j in ((0, 1), (1, 2), (2, 0)):
        o = 3 - i - j
        cr = biv.padAdd(biv.conv2(M[i][0], M[j][1]), biv.scal(biv.conv2(M[i][1], M[j][0]), -1))
        A = biv.padAdd(A, biv.scal(cr, Gd[o]))
    return A


def triangleFaceSign(tf, P):
    """The exact sign of the polynomial P(s1, s2) at the triangle's dual point."""
    if tf['kind'] == 'algebraic':
        return dualPointSign(tf['dp'], P)
    deg = (len(P) - 1) + (len(P[0]) - 1)
    acc = 0
    for i, r in enumerate(P):
        for j, c in enumerate(r):
            if c == 0:
                continue
            acc += c * tf['s0n'][0] ** i * tf['s0n'][1] ** j * tf['s0d'] ** (deg - i - j)
    return sign(acc)

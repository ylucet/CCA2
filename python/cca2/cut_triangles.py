"""cut_triangles -- cut the contact triangles out of an algebraic subdivision and add each one as
its own face (port of cutTriangles.m).

A RATIONAL triangle has rational edges: each face it meets is cut by those three lines and the
triangle becomes an ordinary affine face. An ALGEBRAIC triangle has irrational corners: its
three sides enter as rows of table 3 (AlgCrv evaluates them from the double cache), bounded by a
dyadic box that provably contains the triangle; the triangle is a face of kind 'triangle'.
Returns (subdivision, cuts); with mode 'lines' only the cut lines are returned.
"""
import math
from fractions import Fraction

import numpy as np

from . import rat_q
from .cell_interior import cellInterior
from .conic_meet import _solve2
from .errors import error
from .rat_q import sign
from .triangle_face import triangleFaceSign


def cutTriangles(subdivision, tfs, mode=None):
    if mode == 'lines':
        cuts = []
        for tf in tfs:
            cuts += _triangleLines(tf)
        return [], cuts
    subdivision = [dict(f) for f in subdivision]
    for f in subdivision:
        f['tri'] = []
        f['tf'] = None
        f['tfs'] = tfs
    triF = []
    for t, tf in enumerate(tfs, start=1):
        keep = [True] * len(subdivision)
        out = []
        L = _triangleLines(tf)
        keys = [rat_q.conic([0, 0, 0] + list(l)) for l in L]
        flip = [sign(next(x for x in l if x != 0)) for l in L]   # canonical conic vs the line
        for k, f in enumerate(subdivision):
            sd = []
            for i in range(len(L)):
                v = _sideOf(f['con'], keys[i], flip[i], L[i])
                sd.append(v)
                if v < 0:
                    break                           # outside one line: T misses it
            if any(v < 0 for v in sd):
                continue                            # outside T, or outside its box
            if tf['kind'] == 'rational':
                keep[k] = False
                if all(v > 0 for v in sd):
                    continue                        # inside T: the plane owns it
                for e in range(3):
                    if sd[e] > 0:
                        continue
                    p = dict(f)
                    add = [keys[i] + [flip[i] if i < e else -flip[e], 0] for i in range(e + 1)]
                    p['con'] = [list(r) for r in f['con']] + add
                    out.append(p)
                continue
            keep[k] = False
            for e in range(1, 4):
                p = dict(f)
                p['tri'] = [list(r) for r in f['tri']] + [[t, i, 1] for i in range(1, e)] + [[t, e, -1]]
                out.append(p)
        triF.append(_triangleItself(tf, t, L, keys, flip, tfs))
        subdivision = [f for f, kp in zip(subdivision, keep) if kp] + out
    return subdivision + triF, []


def _sideOf(con, key, flip, line):
    """+1 / -1 when the whole face lies on that side of the line, 0 when the line crosses it."""
    for r in con:
        if r[0:6] == key and (len(r) < 8 or r[7] == 0):
            return r[6] * flip
    V = _polygonOf(con)
    if V:
        v = [p[0] * line[0] + p[1] * line[1] + line[2] for p in V]
        m = 1e-9 * max(1.0, math.sqrt(sum(float(x) ** 2 for x in line))) * \
            max(1.0, max(abs(c) for p in V for c in p))
        if all(x > m for x in v):
            return 1
        if all(x < -m for x in v):
            return -1
    kept = [r for r in con if not all(x == 0 for x in r[0:5])]
    C = [list(r[0:6]) + [r[7] if len(r) >= 8 else 0] for r in kept]
    rows = [[i + 1, r[6]] for i, r in enumerate(kept)]
    C.append(list(key) + [0])
    j = len(C)
    if not cellInterior(C, rows + [[j, -flip]])[0]:
        return 1
    if not cellInterior(C, rows + [[j, flip]])[0]:
        return -1
    return 0


def _polygonOf(con):
    H = [[r[6] * x for x in r[3:6]] for r in con
         if all(x == 0 for x in r[0:3]) and (len(r) < 8 or r[7] == 0) and (r[3] != 0 or r[4] != 0)]
    H = [[float(x) for x in h] for h in H] + [[1, 0, 1e6], [-1, 0, 1e6], [0, 1, 1e6], [0, -1, 1e6]]
    mh = max(abs(x) for h in H for x in h)
    V = []
    for i in range(len(H)):
        for j in range(i + 1, len(H)):
            M = [[H[i][0], H[i][1]], [H[j][0], H[j][1]]]
            nM = math.sqrt(sum(x * x for r in M for x in r))
            if abs(M[0][0] * M[1][1] - M[0][1] * M[1][0]) < 1e-12 * max(1.0, nM) ** 2:
                continue
            x = _solve2(M, [-H[i][2], -H[j][2]])
            lim = -1e-9 * max(1.0, math.hypot(x[0], x[1])) * max(1.0, mh)
            if all(h[0] * x[0] + h[1] * x[1] + h[2] >= lim for h in H):
                V.append(x)
    return V


def _triangleLines(tf):
    """The triangle's three edges as integer rows (inside is >= 0), or, for an algebraic
    triangle, the four sides of a dyadic box that provably contains it."""
    if tf['kind'] != 'rational':
        return _boxAround(tf)
    L = []
    for e in range(3):
        i = e
        j = (e + 1) % 3
        o = 3 - i - j
        d = rat_q.lcm(tf['Td'][i], tf['Td'][j])
        a = [x * (d // tf['Td'][i]) for x in tf['Tn'][i]]
        b = [x * (d // tf['Td'][j]) for x in tf['Tn'][j]]      # corners times d
        n = [-(b[1] - a[1]), b[0] - a[0]]                      # normal of the edge, times d
        r = [n[0] * d, n[1] * d, -(n[0] * a[0] + n[1] * a[1])]
        g = math.gcd(math.gcd(abs(r[0]), abs(r[1])), abs(r[2]))
        r = [x // g for x in r]
        if r[0] * tf['Tn'][o][0] + r[1] * tf['Tn'][o][1] + r[2] * tf['Td'][o] < 0:
            r = [-x for x in r]
        L.append(r)
    return L


def _boxAround(tf):
    T = tf['T']
    mn = [min(p[0] for p in T), min(p[1] for p in T)]
    mx = [max(p[0] for p in T), max(p[1] for p in T)]
    w = max(1.0, max(mx[0] - mn[0], mx[1] - mn[1]))
    for e in range(-8, 1, 2):
        dl = 2.0 ** (math.floor(math.log2(w)) + e)   # a power of two: the ends stay SHORT dyadics
        lo = [math.floor((mn[k] - dl) / dl) * dl for k in range(2)]
        hi = [math.ceil((mx[k] + dl) / dl) * dl for k in range(2)]
        if _allInside(tf, lo, hi):
            L = []
            for k in range(2):
                fl = Fraction(lo[k])
                fh = Fraction(hi[k])
                r1 = [0, 0, -fl.numerator]
                r1[k] = fl.denominator                 # ld x_k - ln >= 0
                r2 = [0, 0, fh.numerator]
                r2[k] = -fh.denominator                # hn - hd x_k >= 0
                L += [r1, r2]
            return L
    error('cutTriangles:box',
          'no dyadic box up to the triangle\'s own size provably contains it.')


def _allInside(tf, lo, hi):
    for i in range(3):
        g = tf['Gn'][i]
        M = [[[g[3], g[1]], [2 * g[0], 0]], [[g[4], 2 * g[2]], [g[1], 0]]]   # grad, over Gd(i)
        for k in range(2):
            fl = Fraction(lo[k])
            fh = Fraction(hi[k])
            A = [[fl.denominator * x for x in r] for r in M[k]]
            A[0][0] -= fl.numerator * tf['Gd'][i]
            B = [[-fh.denominator * x for x in r] for r in M[k]]
            B[0][0] += fh.numerator * tf['Gd'][i]
            if triangleFaceSign(tf, A) < 0 or triangleFaceSign(tf, B) < 0:
                return False
    return True


def _triangleItself(tf, t, L, keys, flip, tfs):
    f = {'FP': None, 'FB': 0, 'alg': [], 'tri': [], 'tf': None, 'tfs': tfs, 'ii': []}
    f['con'] = [keys[i] + [flip[i], 0] for i in range(len(keys))]
    if tf['kind'] == 'rational':
        num, den = _planeQ(tf)
        f.update(kind='quad', num=num, den=den)
        return f
    f.update(kind='triangle', num=[0] * 10, den=1, tf=tf)       # the box, then the sides
    f['tri'] = [[t, 1, 1], [t, 2, 1], [t, 3, 1]]
    assert len(L) == 4
    return f


def _planeQ(tf):
    sn = tf['s0n']
    sd = tf['s0d']
    g = tf['Gn'][0]
    gd = tf['Gd'][0]
    c = (g[0] * sn[0] ** 2 + g[1] * sn[0] * sn[1] + g[2] * sn[1] ** 2
         + (g[3] * sn[0] + g[4] * sn[1]) * sd + g[5] * sd ** 2)
    num = [0] * 10
    num[7] = sn[0] * sd * gd
    num[8] = sn[1] * sd * gd
    num[9] = -c
    return rat_q.canon(num, sd * sd * gd)

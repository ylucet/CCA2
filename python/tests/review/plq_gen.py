"""Input generator for the 2026-11-19 review's differential test (MATLAB vs Python vs Rust vs the
reviewer's exact oracle). Every coordinate and coefficient is an integer, so the three engines
read exactly the same rational input (MATLAB and the Rust golden reader take doubles).

    python3 plq_gen.py OUT.json [--seed N] [--nrand N]

Families
  E-*  exhaustive: Hessian entries in {-1,0,1}^3 (27 quadratics) x 6 domains
       (plane, triangle, square, L-shape, quadrant, segment)
  R1   one random quadratic on a random convex triangle/quadrilateral
  R2   one random quadratic on a random star-shaped (nonconvex) polygon
  R3   a convex quadrilateral cut by a diagonal, two pieces, CONTINUOUS (q2 = q1 + l m)
  R4   a convex pentagon fan-triangulated, three pieces, continuous
  R5   one quadratic on a pointed cone (PSD, or strictly bounded on the cone)
  R6   two continuous pieces on two adjacent cones (a pointed fan)
  D    two DIScontinuous pieces (acceptance test: what do the engines do?)
"""
import json
import random
import sys
from fractions import Fraction as Fr

import plq_oracle as O


def F(v):
    return Fr(v)


def tri_ccw(P):
    s = sum(P[i][0] * P[(i + 1) % len(P)][1] - P[(i + 1) % len(P)][0] * P[i][1] for i in range(len(P)))
    return P if s > 0 else P[::-1]


def q_from(H11, H12, H22, L1=0, L2=0, c=0):
    return tuple(F(v) for v in (H11, H12, H22, L1, L2, c))


def q_add_lm(q, l, m):
    """q + l(x) m(x), l = (a, b, g), m = (d, e, z) affine forms a x + b y + g."""
    a, b, g = l
    d, e, z = m
    H11, H12, H22, L1, L2, c = q
    return (H11 + 2 * a * d, H12 + a * e + b * d, H22 + 2 * b * e,
            L1 + a * z + g * d, L2 + b * z + g * e, c + g * z)


def line_through(p, r):
    """Affine form vanishing on the line through p and r."""
    a = p[1] - r[1]
    b = r[0] - p[0]
    return (a, b, -(a * p[0] + b * p[1]))


def row(q):
    H11, H12, H22, L1, L2, c = q
    return [0, 0, 0, 0, H11, H12, H22, L1, L2, c]


# ---------------------------------------------------------------------------------------------
# encoding: reviewer's pieces -> the engines' (V, E, f, F)
# ---------------------------------------------------------------------------------------------
def encode(pieces):
    kinds = {p['kind'] for p in pieces}
    if kinds == {'plane'}:
        return {'V': [], 'E': [], 'f': [row(pieces[0]['q'])], 'F': []}
    if kinds == {'point'}:
        return {'V': [list(pieces[0]['P'][0])], 'E': [], 'f': [row(pieces[0]['q'])], 'F': []}
    if kinds == {'segment'}:
        a, b = pieces[0]['P']
        return {'V': [list(a), list(b)], 'E': [[1, 2, 1]], 'f': [row(pieces[0]['q'])], 'F': [[0, 0]]}
    V, vid = [], {}

    def v(p):
        p = tuple(p)
        if p not in vid:
            vid[p] = len(V) + 1
            V.append(list(p))
        return vid[p]
    E, Fm = [], []
    if kinds == {'poly'}:
        edges = {}
        for k, p in enumerate(pieces, start=1):
            P = p['P']
            for i in range(len(P)):
                a, b = v(P[i]), v(P[(i + 1) % len(P)])
                if (b, a) in edges:
                    Fm[edges[(b, a)]][1] = k          # neighbour sees it reversed: right side
                else:
                    edges[(a, b)] = len(E)
                    E.append([a, b, 1])
                    Fm.append([k, 0])                  # CCW polygon: face on the LEFT
        return {'V': V, 'E': E, 'f': [row(p['q']) for p in pieces], 'F': Fm}
    if kinds == {'cone'}:
        # a pointed fan: pieces k = 1..m between rays r_{k-1} and r_k (CCW)
        apex = pieces[0]['apex']
        a = v(apex)
        rays = [pieces[0]['d1']] + [p['d2'] for p in pieces]
        for i, d in enumerate(rays):
            b = v((apex[0] + d[0], apex[1] + d[1]))
            E.append([a, b, 0])
            left = i + 1 if i < len(pieces) else 0
            right = i if i > 0 else 0
            Fm.append([left, right])
        return {'V': V, 'E': E, 'f': [row(p['q']) for p in pieces], 'F': Fm}
    raise ValueError(kinds)


# ---------------------------------------------------------------------------------------------
# families
# ---------------------------------------------------------------------------------------------
TRI = [(0, 0), (1, 0), (0, 1)]
SQ = [(0, 0), (1, 0), (1, 1), (0, 1)]
LSH = [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]


def P(pts):
    return [tuple(F(c) for c in p) for p in pts]


def cone_ok(q, d1, d2):
    """Exclude cones where the recession analysis would hinge on alpha == 0 with H d != 0."""
    piece = {'kind': 'cone', 'apex': (F(0), F(0)), 'd1': d1, 'd2': d2, 'q': q}
    H11, H12, H22 = q[0], q[1], q[2]
    psd = H11 >= 0 and H22 >= 0 and H11 * H22 - H12 * H12 >= 0
    if psd:
        return True
    A = [O.dHd(q, ((1 - w) * d1[0] + w * d2[0], (1 - w) * d1[1] + w * d2[1]))
         for w in [Fr(i, 64) for i in range(65)]]
    del piece
    return min(A) > 0 or min(A) < 0 and all(a != 0 for a in A)


def exhaustive():
    out = []
    doms = ['plane', 'tri', 'square', 'L', 'quadrant', 'segment']
    for H11 in (-1, 0, 1):
        for H12 in (-1, 0, 1):
            for H22 in (-1, 0, 1):
                q = q_from(H11, H12, H22)
                for d in doms:
                    if d == 'plane':
                        pcs = [{'kind': 'plane', 'q': q}]
                    elif d == 'tri':
                        pcs = [{'kind': 'poly', 'P': P(TRI), 'q': q}]
                    elif d == 'square':
                        pcs = [{'kind': 'poly', 'P': P(SQ), 'q': q}]
                    elif d == 'L':
                        pcs = [{'kind': 'poly', 'P': P(LSH), 'q': q}]
                    elif d == 'quadrant':
                        d1, d2 = (F(1), F(0)), (F(0), F(1))
                        if not cone_ok(q, d1, d2):
                            continue
                        pcs = [{'kind': 'cone', 'apex': (F(0), F(0)), 'd1': d1, 'd2': d2, 'q': q}]
                    else:
                        pcs = [{'kind': 'segment', 'P': P([(0, 0), (2, 1)]), 'q': q}]
                    out.append(('E-%s-H%d%d%d' % (d, H11, H12, H22), pcs))
    return out


def rq(rng, lo=-2, hi=2):
    return q_from(*(rng.randint(lo, hi) for _ in range(6)))


def rconvex(rng, n):
    while True:
        pts = [(rng.randint(-3, 3), rng.randint(-3, 3)) for _ in range(n)]
        hull = convex_hull(pts)
        if len(hull) == n:
            return tri_ccw(hull)


def convex_hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3:
        return pts

    def half(ps):
        h = []
        for p in ps:
            while len(h) >= 2 and O.cross(h[-2], h[-1], p) <= 0:
                h.pop()
            h.append(p)
        return h
    lo = half(pts)
    up = half(pts[::-1])
    return lo[:-1] + up[:-1]


def rstar(rng, n):
    import math
    while True:
        angs = sorted(rng.uniform(0, 2 * math.pi) for _ in range(n))
        pts = []
        for a in angs:
            r = rng.choice([1, 2, 3])
            pts.append((round(r * math.cos(a) * 2), round(r * math.sin(a) * 2)))
        if len(set(pts)) < n:
            continue
        if not simple(pts):
            continue
        pts = tri_ccw(pts)
        if any(O.cross(pts[i - 1], pts[i], pts[(i + 1) % n]) == 0 for i in range(n)):
            continue
        if all(O.cross(pts[i - 1], pts[i], pts[(i + 1) % n]) > 0 for i in range(n)):
            continue                                   # convex: not what this family is for
        return pts


def simple(pts):
    n = len(pts)

    def inter(a, b, c, d):
        d1, d2 = O.cross(a, b, c), O.cross(a, b, d)
        d3, d4 = O.cross(c, d, a), O.cross(c, d, b)
        if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)) and d1 * d2 < 0 and d3 * d4 < 0:
            return True
        return False
    for i in range(n):
        for j in range(i + 1, n):
            if abs(i - j) in (1, n - 1):
                continue
            if inter(pts[i], pts[(i + 1) % n], pts[j], pts[(j + 1) % n]):
                return False
    return True


def rlin(rng):
    return tuple(F(rng.randint(-2, 2)) for _ in range(3))


def random_family(rng, nrand):
    out = []
    for i in range(nrand):
        fam = ['R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'D'][i % 7]
        if fam == 'R1':
            Pp = rconvex(rng, rng.choice([3, 4]))
            pcs = [{'kind': 'poly', 'P': P(Pp), 'q': rq(rng)}]
        elif fam == 'R2':
            Pp = rstar(rng, rng.choice([5, 6]))
            pcs = [{'kind': 'poly', 'P': P(Pp), 'q': rq(rng)}]
        elif fam in ('R3', 'D'):
            Pp = P(rconvex(rng, 4))
            A = [Pp[0], Pp[1], Pp[2]]
            B = [Pp[0], Pp[2], Pp[3]]
            q1 = rq(rng)
            if fam == 'R3':
                m = rlin(rng)
                if m == (0, 0, 0):
                    m = (F(1), F(0), F(0))
                q2 = q_add_lm(q1, line_through(Pp[0], Pp[2]), m)
            else:
                q2 = rq(rng)
                if all(O.qval(q1, x) == O.qval(q2, x) for x in (Pp[0], Pp[2], ((Pp[0][0] + Pp[2][0]) / 2, (Pp[0][1] + Pp[2][1]) / 2))):
                    q2 = q_add_lm(q2, (F(0), F(0), F(1)), (F(0), F(0), F(1)))
            pcs = [{'kind': 'poly', 'P': A, 'q': q1}, {'kind': 'poly', 'P': B, 'q': q2}]
        elif fam == 'R4':
            Pp = P(rconvex(rng, 5))
            q = rq(rng)
            pcs = []
            for t in range(3):
                tri = [Pp[0], Pp[t + 1], Pp[t + 2]]
                if t > 0:
                    m = rlin(rng)
                    q = q_add_lm(q, line_through(Pp[0], Pp[t + 1]), m)
                pcs.append({'kind': 'poly', 'P': tri, 'q': q})
        elif fam == 'R5':
            while True:
                d1 = (F(rng.randint(-2, 2)), F(rng.randint(-2, 2)))
                d2 = (F(rng.randint(-2, 2)), F(rng.randint(-2, 2)))
                cr = d1[0] * d2[1] - d1[1] * d2[0]
                if cr <= 0:
                    continue
                q = rq(rng)
                if cone_ok(q, d1, d2):
                    break
            apex = (F(rng.randint(-1, 1)), F(rng.randint(-1, 1)))
            pcs = [{'kind': 'cone', 'apex': apex, 'd1': d1, 'd2': d2, 'q': q}]
        else:  # R6
            while True:
                ds = [(F(rng.randint(-2, 2)), F(rng.randint(-2, 2))) for _ in range(3)]
                ok = all(ds[j][0] * ds[j + 1][1] - ds[j][1] * ds[j + 1][0] > 0 for j in range(2))
                tot = ds[0][0] * ds[2][1] - ds[0][1] * ds[2][0]
                if not ok or tot <= 0:
                    continue
                q1 = q_from(*(rng.randint(0, 2) for _ in range(3)), rng.randint(-2, 2), rng.randint(-2, 2), 0)
                m = rlin(rng)
                q2 = q_add_lm(q1, line_through((F(0), F(0)), ds[1]), m)
                if cone_ok(q1, ds[0], ds[1]) and cone_ok(q2, ds[1], ds[2]):
                    break
            apex = (F(0), F(0))
            pcs = [{'kind': 'cone', 'apex': apex, 'd1': ds[0], 'd2': ds[1], 'q': q1},
                   {'kind': 'cone', 'apex': apex, 'd1': ds[1], 'd2': ds[2], 'q': q2}]
        out.append(('%s-%03d' % (fam, i), pcs))
    return out


# ---------------------------------------------------------------------------------------------
# evaluation points (dyadic, so the doubles the engines read are exact)
# ---------------------------------------------------------------------------------------------
def dual_points():
    S = []
    for a in range(-6, 7):
        for b in range(-6, 7):
            S.append((Fr(a, 2) + Fr(1, 32), Fr(b, 2) + Fr(3, 64)))
    for a in range(-3, 4):
        for b in range(-3, 4):
            S.append((Fr(a), Fr(b)))
    return S


def primal_points(pieces):
    xs, ys = [], []
    for p in pieces:
        if p['kind'] in ('poly', 'segment', 'point'):
            xs += [v[0] for v in p['P']]
            ys += [v[1] for v in p['P']]
        elif p['kind'] == 'cone':
            a = p['apex']
            for d in (p['d1'], p['d2']):
                xs += [a[0], a[0] + 2 * d[0]]
                ys += [a[1], a[1] + 2 * d[1]]
        else:
            xs += [-2, 2]
            ys += [-2, 2]
    lo = (min(xs) - Fr(1, 2), min(ys) - Fr(1, 2))
    hi = (max(xs) + Fr(1, 2), max(ys) + Fr(1, 2))
    X = []
    n = 12
    for i in range(n + 1):
        for j in range(n + 1):
            X.append((lo[0] + (hi[0] - lo[0]) * Fr(i, n) + Fr(1, 128),
                      lo[1] + (hi[1] - lo[1]) * Fr(j, n) + Fr(3, 256)))
    return X


def ser(v):
    if isinstance(v, Fr):
        return str(v)
    if isinstance(v, (list, tuple)):
        return [ser(x) for x in v]
    if isinstance(v, dict):
        return {k: ser(x) for k, x in v.items()}
    return v


def num(v):
    """Integer data as int, dyadic as float (exact in double)."""
    if isinstance(v, Fr):
        return int(v) if v.denominator == 1 else float(v)
    if isinstance(v, (list, tuple)):
        return [num(x) for x in v]
    return v


def build(seed=20261119, nrand=210):
    rng = random.Random(seed)
    cases = exhaustive() + random_family(rng, nrand)
    S = dual_points()
    out = []
    for cid, (name, pcs) in enumerate(cases, start=1):
        enc = encode(pcs)
        X = primal_points(pcs)
        out.append({'id': cid, 'name': name, 'pieces': ser(pcs),
                    'V': num(enc['V']), 'E': enc['E'], 'f': num(enc['f']), 'F': enc['F'],
                    'S': num(S), 'X': num(X)})
    return out


def load_pieces(js):
    def fr(v):
        return Fr(v)
    pcs = []
    for p in js:
        d = {'kind': p['kind'], 'q': tuple(fr(v) for v in p['q'])}
        if 'P' in p:
            d['P'] = [tuple(fr(c) for c in v) for v in p['P']]
        if 'apex' in p:
            d['apex'] = tuple(fr(c) for c in p['apex'])
            d['d1'] = tuple(fr(c) for c in p['d1'])
            d['d2'] = tuple(fr(c) for c in p['d2'])
        pcs.append(d)
    return {'pieces': pcs}


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--seed', type=int, default=20261119)
    ap.add_argument('--nrand', type=int, default=210)
    a = ap.parse_args()
    fx = build(a.seed, a.nrand)
    json.dump(fx, open(a.out, 'w'))
    print('%d fixtures -> %s' % (len(fx), a.out))

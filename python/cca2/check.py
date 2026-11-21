"""check -- verify a conjugate or a biconjugate against the DEFINITION, independently of MATLAB
and of the engine (the role plqCheck.m plays for the MATLAB toolbox).

    supAt(f, s)       sup_x <s,x> - f(x), computed pointwise in doubles by enumerating the KKT
                      candidates of each piece (vertices, edge stationary points, interior
                      stationary point, recession directions). No cells, no conics: nothing in
                      common with conjQ's construction.
    evalQuaPol(f, X)  f at points, +inf off the domain.
    checkConj(f, g)   g(s) == supAt(f, s) at sample dual points.
    checkBiconj(f, h) h <= f on dom f;  h convex on sampled segments;  and at each sample x,
                      h(x) == <s,x> - f*(s) for s the gradient of h at x  (f** is the sup of the
                      affine minorants of f, and the sup is attained at a subgradient).

Each check returns a list of defects (empty when the definition holds on the samples). Doubles
and tolerances throughout: this module only JUDGES results, it never produces one.
"""
import math
from fractions import Fraction as Fr

import numpy as np

INF = math.inf


# ---- the input function -------------------------------------------------------------------------
def _face(f, k):
    n = [float(x) / f.fD[k - 1] for x in f.fN[k - 1]]
    H = [[n[4], n[5]], [n[5], n[6]]]
    L = [n[7], n[8]]
    return H, L, n[9]


def _q(H, L, c, x):
    return (0.5 * (x[0] * (H[0][0] * x[0] + H[0][1] * x[1]) + x[1] * (H[1][0] * x[0] + H[1][1] * x[1]))
            + L[0] * x[0] + L[1] * x[1] + c)


def _V(f):
    return [[float(v) for v in f.Vx(i)] for i in range(1, len(f.VN) + 1)]


def _pieces(f):
    """[(H, L, c, edges, k)]: each piece's quadratic and its edges (a, b, isRay) in doubles;
    k = 0 for an edge of a chain (the edge itself is the piece)."""
    V = _V(f)
    out = []
    if f.nv == 0:
        H, L, c = _face(f, 1)
        return [(H, L, c, None, 1)]
    if f.dom['dim'] < 2:
        if not f.E:
            H, L, c = _face(f, 1)
            return [(H, L, c, [(V[0], V[0], False)], 0)]
        for j, e in enumerate(f.E):
            H, L, c = _face(f, j + 1 if len(f.fN) > 1 else 1)
            out.append((H, L, c, [(V[e[0] - 1], V[e[1] - 1], e[2] == 0)], 0))
        return out
    for k in range(1, len(f.fN) + 1):
        H, L, c = _face(f, k)
        edges = [(V[e[0] - 1], V[e[1] - 1], e[2] == 0)
                 for e, fr in zip(f.E, f.F) if fr[0] == k or fr[1] == k]
        out.append((H, L, c, edges, k))
    return out


def _inFace(f, k, x, tol=1e-9):
    """Is x in the CLOSED face k? Bounded faces by ray casting (convex or not); unbounded faces
    as the intersection of their edges' half-planes (the side F records)."""
    V = _V(f)
    own = [(e, fr) for e, fr in zip(f.E, f.F) if fr[0] == k or fr[1] == k]
    sc = max(1.0, abs(x[0]), abs(x[1]))
    # on an edge?
    for e, _ in own:
        a = V[e[0] - 1]
        b = V[e[1] - 1]
        d = [b[0] - a[0], b[1] - a[1]]
        L2 = d[0] * d[0] + d[1] * d[1]
        if L2 == 0:
            continue
        t = ((x[0] - a[0]) * d[0] + (x[1] - a[1]) * d[1]) / L2
        cr = abs(d[0] * (x[1] - a[1]) - d[1] * (x[0] - a[0])) / math.sqrt(L2)
        if cr <= tol * sc and t >= -tol and (e[2] == 0 or t <= 1 + tol):
            return True
    segs = [(V[e[0] - 1], V[e[1] - 1]) for e, _ in own if e[2] != 0]
    rays = [(V[e[0] - 1], [V[e[1] - 1][0] - V[e[0] - 1][0], V[e[1] - 1][1] - V[e[0] - 1][1]],
             fr[0] == k) for e, fr in own if e[2] == 0]
    if rays:
        # An unbounded face, convex or not: close it far away. Each ray with the face on its
        # LEFT is joined, counter-clockwise through the face, to the next ray with the face on
        # its right, by an arc of far points; the result is a (huge) polygon.
        R = 1e4 * max(1.0, max(abs(c) for v in V for c in v), abs(x[0]), abs(x[1]))
        lefts = [(a, d) for a, d, left in rays if left]
        rights = [(a, d) for a, d, left in rays if not left]
        for a, d in rays_iter(rays):
            n = math.hypot(*d)
            segs.append((a, [a[0] + R * d[0] / n, a[1] + R * d[1] / n]))
        for aL, dL in lefts:
            th0 = math.atan2(dL[1], dL[0])
            best = None
            for aR, dR in rights:
                th = (math.atan2(dR[1], dR[0]) - th0) % (2 * math.pi)
                if best is None or th < best[0]:
                    best = (th, aR, dR)
            if best is None:
                continue
            th, aR, dR = best
            nL = math.hypot(*dL)
            nR = math.hypot(*dR)
            pts = [[aL[0] + R * dL[0] / nL, aL[1] + R * dL[1] / nL]]
            m = max(1, int(math.ceil(th / (math.pi / 4))))
            for i in range(1, m):
                t = th0 + th * i / m
                pts.append([R * math.cos(t), R * math.sin(t)])
            pts.append([aR[0] + R * dR[0] / nR, aR[1] + R * dR[1] / nR])
            segs += list(zip(pts[:-1], pts[1:]))
    inside = False
    for a, b in segs:
        if (a[1] > x[1]) != (b[1] > x[1]):
            xc = a[0] + (x[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if xc > x[0]:
                inside = not inside
    return inside


def rays_iter(rays):
    return [(a, d) for a, d, _ in rays]


def evalQuaPol(f, X):
    """f at each point of X (doubles); +inf off the domain; the smaller value on a shared edge."""
    out = []
    V = _V(f)
    for x in X:
        x = [float(x[0]), float(x[1])]
        if f.nv == 0:
            H, L, c = _face(f, 1)
            out.append(_q(H, L, c, x))
            continue
        best = INF
        if f.dom['dim'] < 2:
            sc = max(1.0, abs(x[0]), abs(x[1]))
            if not f.E:
                if math.hypot(x[0] - V[0][0], x[1] - V[0][1]) <= 1e-9 * sc:
                    best = _q(*_face(f, 1), x)
            for j, e in enumerate(f.E):
                a = V[e[0] - 1]
                b = V[e[1] - 1]
                d = [b[0] - a[0], b[1] - a[1]]
                L2 = d[0] * d[0] + d[1] * d[1]
                t = ((x[0] - a[0]) * d[0] + (x[1] - a[1]) * d[1]) / L2
                cr = abs(d[0] * (x[1] - a[1]) - d[1] * (x[0] - a[0])) / math.sqrt(L2)
                if cr <= 1e-9 * sc and t >= -1e-9 and (e[2] == 0 or t <= 1 + 1e-9):
                    best = min(best, _q(*_face(f, j + 1 if len(f.fN) > 1 else 1), x))
            out.append(best)
            continue
        for k in range(1, len(f.fN) + 1):
            if _inFace(f, k, x):
                best = min(best, _q(*_face(f, k), x))
        out.append(best)
    return out


# ---- the conjugate, pointwise, from the definition -----------------------------------------------
def _supOnEdge(H, L, c, a, b, isRay, s):
    """sup over the segment [a, b] (or the ray from a through b) of <s,x> - q(x)."""
    d = [b[0] - a[0], b[1] - a[1]]

    def val(t):
        x = [a[0] + t * d[0], a[1] + t * d[1]]
        return s[0] * x[0] + s[1] * x[1] - _q(H, L, c, x)
    if d[0] == 0 and d[1] == 0:
        return val(0.0)
    al = 0.5 * (d[0] * (H[0][0] * d[0] + H[0][1] * d[1]) + d[1] * (H[1][0] * d[0] + H[1][1] * d[1]))
    g = [H[0][0] * a[0] + H[0][1] * a[1] + L[0], H[1][0] * a[0] + H[1][1] * a[1] + L[1]]
    be = (s[0] - g[0]) * d[0] + (s[1] - g[1]) * d[1]      # slope at t = 0
    sc = max(1.0, abs(be), abs(al))
    best = val(0.0)
    if isRay:
        if al < -1e-12 * sc or (abs(al) <= 1e-12 * sc and be > 1e-9 * sc):
            return INF
    else:
        best = max(best, val(1.0))
    if al > 1e-12 * sc:
        t = be / (2 * al)
        if t > 0 and (isRay or t < 1):
            best = max(best, val(t))
    return best


def supAt(f, s):
    """sup_x <s,x> - f(x) at the dual point s, in doubles (+inf when it diverges)."""
    s = [float(s[0]), float(s[1])]
    best = -INF
    for H, L, c, edges, k in _pieces(f):
        Hm = np.array(H)
        ev = np.linalg.eigvalsh(Hm)
        if edges is None:                                   # the whole plane
            if ev[0] <= 1e-12 * max(1.0, abs(ev[1])):
                # singular or indefinite: finite only in degenerate cases the engine names
                if ev[0] < -1e-12 * max(1.0, abs(ev[1])):
                    return INF
                r = np.array([s[0] - L[0], s[1] - L[1]])
                x, res, rank, _ = np.linalg.lstsq(Hm, r, rcond=None)
                if np.linalg.norm(Hm @ x - r) > 1e-9 * max(1.0, np.linalg.norm(r)):
                    return INF
                best = max(best, s[0] * x[0] + s[1] * x[1] - _q(H, L, c, x))
            else:
                x = np.linalg.solve(Hm, np.array([s[0] - L[0], s[1] - L[1]]))
                best = max(best, s[0] * x[0] + s[1] * x[1] - _q(H, L, c, x))
            continue
        for a, b, isRay in edges:
            v = _supOnEdge(H, L, c, a, b, isRay, s)
            if v == INF:
                return INF
            best = max(best, v)
        if k == 0:
            continue
        rays = [[b[0] - a[0], b[1] - a[1]] for a, b, isRay in edges if isRay]
        if len(rays) == 2:
            # directions strictly inside the recession cone
            base = edges[0][0]
            for w in np.linspace(0.02, 0.98, 49):
                d = [(1 - w) * rays[0][0] / math.hypot(*rays[0]) + w * rays[1][0] / math.hypot(*rays[1]),
                     (1 - w) * rays[0][1] / math.hypot(*rays[0]) + w * rays[1][1] / math.hypot(*rays[1])]
                if d[0] == 0 and d[1] == 0:
                    continue
                far = [base[0] + d[0], base[1] + d[1]]
                # any point of the face works as a base; use the centroid of the vertices
                vs = [e[0] for e in edges]
                cx = [sum(p[0] for p in vs) / len(vs) + 1e-3 * d[0],
                      sum(p[1] for p in vs) / len(vs) + 1e-3 * d[1]]
                if _inFace(f, k, cx):
                    far = [cx[0] + d[0], cx[1] + d[1]]
                    if _supOnEdge(H, L, c, cx, far, True, s) == INF:
                        return INF
            # A NULL DIRECTION STRICTLY INSIDE THE CONE (review 2026-11-19, finding 5): the
            # samples above can miss it, so decide it exactly from the face's rational data.
            # A singular H has a null LINE; when it lies in the closed recession cone the
            # objective is LINEAR along it (any base point works: the piece is convex and the
            # direction recedes from every point), so the sup is +inf exactly when the slope
            # is positive. At most one of the antipodal pair lies in a proper cone.
            n0 = f.fN[k - 1]
            dd = f.fD[k - 1]
            H11, H12, H22 = Fr(n0[4], dd), Fr(n0[5], dd), Fr(n0[6], dd)
            if H11 * H22 - H12 * H12 == 0 and (H11 or H12 or H22):
                nv = (-H12, H11) if H11 else (1, 0)              # spans null(H) exactly
                def vx(i):                                       # exact vertex i
                    return [Fr(f.VN[i - 1][0], f.VD[i - 1]), Fr(f.VN[i - 1][1], f.VD[i - 1])]
                ea = [e for e, fr in zip(f.E, f.F) if (fr[0] == k or fr[1] == k) and e[2] == 0][:2]
                r1 = [vx(ea[0][1])[0] - vx(ea[0][0])[0], vx(ea[0][1])[1] - vx(ea[0][0])[1]]
                r2 = [vx(ea[1][1])[0] - vx(ea[1][0])[0], vx(ea[1][1])[1] - vx(ea[1][0])[1]]

                def in2(d):
                    cr = r1[0] * r2[1] - r1[1] * r2[0]
                    if cr == 0:
                        return (r1[0] * d[1] - r1[1] * d[0] == 0
                                and (r1[0] * r2[0] + r1[1] * r2[1] < 0
                                     or r1[0] * d[0] + r1[1] * d[1] >= 0))
                    return (cr * (r1[0] * d[1] - r1[1] * d[0]) >= 0
                            and cr * (d[0] * r2[1] - d[1] * r2[0]) >= 0)
                Lx, Ly = Fr(n0[7], dd), Fr(n0[8], dd)           # exact L (L above is doubles)
                for d in (nv, [-nv[0], -nv[1]]):
                    if in2(d):
                        slope = Fr(s[0]) * d[0] + Fr(s[1]) * d[1] - (Lx * d[0] + Ly * d[1])
                        if slope > 0:
                            return INF
                        break
        if ev[0] > 1e-12 * max(1.0, abs(ev[1])):            # an interior stationary point
            x = np.linalg.solve(Hm, np.array([s[0] - L[0], s[1] - L[1]]))
            if _inFace(f, k, [float(x[0]), float(x[1])]):
                best = max(best, s[0] * x[0] + s[1] * x[1] - _q(H, L, c, x))
        elif ev[0] > -1e-12 * max(1.0, abs(ev[1])) and ev[1] > 1e-12:
            # H is PSD and SINGULAR: the stationary points, if any, form a LINE x0 + t n on
            # which the objective is constant; any point of that line inside the face attains it.
            r = np.array([s[0] - L[0], s[1] - L[1]])
            x0, *_ = np.linalg.lstsq(Hm, r, rcond=None)
            if np.linalg.norm(Hm @ x0 - r) <= 1e-9 * max(1.0, np.linalg.norm(r)):
                w, Vv = np.linalg.eigh(Hm)
                nvec = Vv[:, 0]
                ts = [0.0]
                for a, b, isRay in edges:
                    d = np.array([b[0] - a[0], b[1] - a[1]])
                    M = np.array([[nvec[0], -d[0]], [nvec[1], -d[1]]])
                    if abs(np.linalg.det(M)) > 1e-12:
                        ts.append(float(np.linalg.solve(M, np.array([a[0] - x0[0], a[1] - x0[1]]))[0]))
                ts = sorted(ts)
                cand = [ts[0] - 1.0] + [(u + v) / 2 for u, v in zip(ts[:-1], ts[1:])] + [ts[-1] + 1.0] + ts
                for t in cand:
                    x = x0 + t * nvec
                    if _inFace(f, k, [float(x[0]), float(x[1])]):
                        best = max(best, s[0] * x[0] + s[1] * x[1] - _q(H, L, c, x))
                        break
    return float(best)


# ---- the checks ------------------------------------------------------------------------------------
def _val(h, X):
    r = h.eval(X)
    return r[0] if isinstance(r, tuple) else r


def defaultDualPoints():
    g = [-3.3, -1.7, -0.6, 0.15, 0.8, 1.9, 3.1]
    return [[a, b + 0.037 * a] for a in g for b in g]


def checkConj(f, g, S=None, tol=1e-7):
    """Defects of g as the conjugate of f at the dual points S (default: a 7x7 grid)."""
    S = S or defaultDualPoints()
    gv = _val(g, S)
    bad = []
    for s, a in zip(S, gv):
        b = supAt(f, s)
        if math.isinf(a) or math.isinf(b):
            if a != b:
                # a point within rounding of the boundary of dom f* may fall either way
                near = [supAt(f, [s[0] + dx, s[1] + dy]) for dx in (-1e-6, 1e-6)
                        for dy in (-1e-6, 1e-6)]
                gn = _val(g, [[s[0] + dx, s[1] + dy] for dx in (-1e-6, 1e-6)
                              for dy in (-1e-6, 1e-6)])
                if all(math.isinf(x) == math.isinf(y) for x, y in zip(near, gn)):
                    bad.append('f*(%r): engine %r, definition %r' % (s, a, b))
            continue
        if abs(a - b) > tol * max(1.0, abs(b)):
            bad.append('f*(%r): engine %.12g, definition %.12g' % (s, a, b))
    return bad


def domainSamples(f, n=7):
    """Points spread over (a bounding box of) dom f."""
    V = _V(f)
    if not V:
        g = np.linspace(-2.6, 2.6, n)
        return [[float(a), float(b + 0.031)] for a in g for b in g]
    xs = [v[0] for v in V]
    ys = [v[1] for v in V]
    pad = 0.0 if all(e[2] != 0 for e in f.E) else 1.5 * max(1.0, max(xs) - min(xs), max(ys) - min(ys))
    gx = np.linspace(min(xs) - pad, max(xs) + pad, n + 2)[1:-1]
    gy = np.linspace(min(ys) - pad, max(ys) + pad, n + 2)[1:-1]
    return [[float(a), float(b)] for a in gx for b in gy]


def checkBiconj(f, h, X=None, tol=1e-6):
    """Defects of h as the biconjugate (closed convex envelope) of f.

    For a BOUNDED domain the points outside the vertex box are always outside conv(dom f),
    where the biconjugate is +infinity -- finding 9a of the 2026-11-19 review was invisible to
    samples inside the box, so the outside is sampled too (finding 5)."""
    X = X or domainSamples(f)
    outside = []
    if f.nv > 0 and all(e[2] != 0 for e in f.E):
        V = _V(f)
        xs = [v[0] for v in V]
        ys = [v[1] for v in V]
        lo, hi = (min(xs), min(ys)), (max(xs), max(ys))
        px = 0.5 * max(1.0, hi[0] - lo[0])
        py = 0.5 * max(1.0, hi[1] - lo[1])
        outside = [[lo[0] - px, lo[1] - py], [hi[0] + px, hi[1] + py], [lo[0] - px, hi[1] + py],
                   [hi[0] + px, lo[1] - py], [lo[0] - px, (lo[1] + hi[1]) / 2],
                   [hi[0] + px, (lo[1] + hi[1]) / 2], [(lo[0] + hi[0]) / 2, lo[1] - py],
                   [(lo[0] + hi[0]) / 2, hi[1] + py]]
        X = list(X) + outside
    fv = evalQuaPol(f, X)
    hv = _val(h, X)
    bad = []
    # (0) outside conv(dom f) the biconjugate is +infinity
    for x, a in zip(outside, hv[len(X) - len(outside):] if outside else []):
        if not (a == INF):
            bad.append('f**(%r) = %r outside conv(dom f), where it is +infinity' % (x, a))
    # (1) h <= f wherever f is finite
    for x, a, b in zip(X, hv, fv):
        if math.isfinite(b) and not (a <= b + tol * max(1.0, abs(b))):
            bad.append('f**(%r) = %.12g exceeds f = %.12g' % (x, a, b))
    # (2) h is convex on segments between sample points where it is finite
    fin = [(x, a) for x, a in zip(X, hv) if math.isfinite(a)]
    step = max(1, len(fin) // 12)
    for i in range(0, len(fin), step):
        for j in range(i + 1, len(fin), step):
            (xa, va), (xb, vb) = fin[i], fin[j]
            xm = [(xa[0] + xb[0]) / 2, (xa[1] + xb[1]) / 2]
            vm = _val(h, [xm])[0]
            if not (vm <= (va + vb) / 2 + tol * max(1.0, abs(va), abs(vb))):
                bad.append('f** not convex between %r and %r: midpoint %.12g > %.12g'
                           % (xa, xb, vm, (va + vb) / 2))
    # (3) the definition: h(x) = <s,x> - f*(s) at s = grad h(x)
    for x, a in fin:
        d = 1e-5 * max(1.0, abs(x[0]), abs(x[1]))
        nb = _val(h, [[x[0] + d, x[1]], [x[0] - d, x[1]], [x[0], x[1] + d], [x[0], x[1] - d]])
        if not all(math.isfinite(v) for v in nb):
            continue                                        # too close to the boundary of dom f**
        s = [(nb[0] - nb[1]) / (2 * d), (nb[2] - nb[3]) / (2 * d)]
        # skip a kink of h: the gradient is not defined there
        if abs(nb[0] + nb[1] - 2 * a) > 1e-7 * max(1.0, abs(a)) or \
           abs(nb[2] + nb[3] - 2 * a) > 1e-7 * max(1.0, abs(a)):
            continue
        fs = supAt(f, s)
        if not math.isfinite(fs):
            # the NUMERICAL gradient carries ~1e-11 noise (central difference of values that
            # themselves carry double-precision noise); when the true gradient sits ON the
            # boundary of dom f* (finite there) the noise can push s just outside, and the exact
            # supAt then rightly returns +inf at s. Round s off the rounding band and retest: a
            # point within rounding of the boundary is not a defect (checkConj's same policy).
            # 9 significant digits swallows the ~1e-11 noise while a genuine O(1) boundary miss
            # survives -- except for a true component 0, which read as 7.4e-12 keeps all its
            # digits (cone records 7 and 32): a component within 1e-9 of the gradient's scale is 0.
            sc = max(1.0, abs(s[0]), abs(s[1]))
            sr = [0.0 if abs(v) <= 1e-9 * sc else float('%.9g' % v) for v in s]
            if sr != s and math.isfinite(supAt(f, sr)):
                s, fs = sr, supAt(f, sr)
            else:
                bad.append('f**(%r): the gradient %r of the result is outside dom f*' % (x, s))
                continue
        lower = s[0] * x[0] + s[1] * x[1] - fs
        if abs(a - lower) > 1e-4 * max(1.0, abs(a)):
            bad.append('f**(%r) = %.10g but <s,x> - f*(s) = %.10g at its own gradient'
                       % (x, a, lower))
    return bad

"""Independent oracles for the 2026-11-19 code review (reviews/codereview-2026-11-19.md).

Written by the reviewer from the DEFINITIONS, sharing no code with cca2 (not even cca2.check):

  conj_exact(fx, s)   f*(s) = sup_x <s,x> - f(x), EXACTLY (fractions.Fraction) at a rational s.
                      A quadratic bounded above on a polyhedron attains its sup (Frank-Wolfe), so
                      the sup over each closed piece is the max over its KKT candidates: vertices,
                      the critical point of each edge, and the interior stationary set (a point,
                      or a line when H is singular). Unboundedness is decided exactly on the
                      piece's recession cone. Returns None for +inf.
  f_at(fx, x)         f(x) exactly (min over the closed pieces containing x; None off the domain).
  biconj_bracket(...) a certified bracket  lo <= f**(x) <= hi:
                        lo = <s,x> - f*(s) for any s (Fenchel-Young; s is the caller's choice);
                        hi = min sum(l_i f(x_i)) over convex combinations of sampled domain points
                             (an LP), an upper bound on the convex envelope.

A fixture ('fx') is the reviewer's own description of the input, from which the engine encoding
(V, E, f, F) is derived by plq_gen.encode:
  fx['pieces'] = [ {'q': (H11, H12, H22, L1, L2, c) as Fractions  -- q(x) = x'Hx/2 + L'x + c,
                    'kind': 'poly' | 'cone' | 'plane' | 'segment' | 'point',
                    'P': [vertices CCW] (poly) | apex + two directions (cone) ...} ]
"""
from fractions import Fraction as Fr


# ---------------------------------------------------------------------------------------------
# quadratic helpers
# ---------------------------------------------------------------------------------------------
def qval(q, x):
    H11, H12, H22, L1, L2, c = q
    return (H11 * x[0] * x[0] + 2 * H12 * x[0] * x[1] + H22 * x[1] * x[1]) / 2 \
        + L1 * x[0] + L2 * x[1] + c


def qgrad(q, x):
    H11, H12, H22, L1, L2, _ = q
    return (H11 * x[0] + H12 * x[1] + L1, H12 * x[0] + H22 * x[1] + L2)


def dHd(q, d):
    H11, H12, H22 = q[0], q[1], q[2]
    return H11 * d[0] * d[0] + 2 * H12 * d[0] * d[1] + H22 * d[1] * d[1]


def obj(q, s, x):
    return s[0] * x[0] + s[1] * x[1] - qval(q, x)


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


# ---------------------------------------------------------------------------------------------
# membership (exact)
# ---------------------------------------------------------------------------------------------
def on_segment(a, b, x):
    if cross(a, b, x) != 0:
        return False
    return min(a[0], b[0]) <= x[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= x[1] <= max(a[1], b[1])


def in_polygon(P, x):
    """Closed simple polygon P (any orientation, convex or not)."""
    n = len(P)
    for i in range(n):
        if on_segment(P[i], P[(i + 1) % n], x):
            return True
    inside = False
    for i in range(n):
        a, b = P[i], P[(i + 1) % n]
        if (a[1] > x[1]) != (b[1] > x[1]):
            xc = a[0] + (x[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
            if xc > x[0]:
                inside = not inside
    return inside


def in_cone(apex, d1, d2, x):
    """Closed cone apex + cone(d1, d2), d1 -> d2 counter-clockwise, angle < pi."""
    v = (x[0] - apex[0], x[1] - apex[1])
    return (d1[0] * v[1] - d1[1] * v[0]) >= 0 and (v[0] * d2[1] - v[1] * d2[0]) >= 0


def in_piece(p, x):
    k = p['kind']
    if k == 'plane':
        return True
    if k == 'poly':
        return in_polygon(p['P'], x)
    if k == 'cone':
        return in_cone(p['apex'], p['d1'], p['d2'], x)
    if k == 'segment':
        return on_segment(p['P'][0], p['P'][1], x)
    if k == 'point':
        return tuple(x) == tuple(p['P'][0])
    raise ValueError(k)


def f_at(fx, x):
    best = None
    for p in fx['pieces']:
        if in_piece(p, x):
            v = qval(p['q'], x)
            best = v if best is None or v < best else best
    return best


# ---------------------------------------------------------------------------------------------
# the conjugate, exactly
# ---------------------------------------------------------------------------------------------
INF = None


def _sup_segment(q, s, a, b):
    """max over the closed segment [a, b] of <s,x> - q(x): endpoints and the critical point."""
    best = max(obj(q, s, a), obj(q, s, b))
    d = (b[0] - a[0], b[1] - a[1])
    al = dHd(q, d) / 2                          # phi(t) = phi(0) + be t - al t^2
    g = qgrad(q, a)
    be = (s[0] - g[0]) * d[0] + (s[1] - g[1]) * d[1]
    if al > 0:
        t = be / (2 * al)
        if 0 < t < 1:
            best = max(best, obj(q, s, (a[0] + t * d[0], a[1] + t * d[1])))
    return best


def _sup_ray(q, s, a, d):
    """max over the ray a + t d, t >= 0; None for +inf."""
    al = dHd(q, d) / 2
    g = qgrad(q, a)
    be = (s[0] - g[0]) * d[0] + (s[1] - g[1]) * d[1]
    if al < 0 or (al == 0 and be > 0):
        return INF
    best = obj(q, s, a)
    if al > 0:
        t = be / (2 * al)
        if t > 0:
            best = max(best, obj(q, s, (a[0] + t * d[0], a[1] + t * d[1])))
    return best


def _stationary(q, s):
    """The stationary set of x -> <s,x> - q(x): ('point', x) | ('line', x0, dir) | ('plane',) |
    None (no stationary point)."""
    H11, H12, H22, L1, L2, _ = q
    r = (s[0] - L1, s[1] - L2)
    det = H11 * H22 - H12 * H12
    if det != 0:
        return ('point', ((H22 * r[0] - H12 * r[1]) / det, (H11 * r[1] - H12 * r[0]) / det))
    if H11 == 0 and H12 == 0 and H22 == 0:
        return ('plane',) if r == (0, 0) else None
    # rank 1: H = rows (H11, H12) and (H12, H22) are parallel; solve the nonzero row
    if H11 != 0 or H12 != 0:
        a, b, rr = H11, H12, r[0]
        # consistency of the second row: (H12, H22) = lam (H11, H12)
        lam = (H12 / H11) if H11 != 0 else (H22 / H12)
        if r[1] != lam * r[0]:
            return None
    else:
        a, b, rr = H12, H22, r[1]
        if r[0] != 0:
            return None
    # a x + b y = rr
    x0 = (rr / a, Fr(0)) if a != 0 else (Fr(0), rr / b)
    return ('line', x0, (-b, a))


def _line_hits_piece(p, x0, d):
    """Does the line x0 + t d meet the closed piece? Candidates: intersections with the
    piece's boundary lines, plus x0 itself."""
    cands = [x0]
    for (a, b) in _boundary_lines(p):
        # solve x0 + t d = a + u (b - a)
        e = (b[0] - a[0], b[1] - a[1])
        den = d[0] * (-e[1]) - d[1] * (-e[0])
        if den == 0:
            continue
        rx, ry = a[0] - x0[0], a[1] - x0[1]
        t = (rx * (-e[1]) - ry * (-e[0])) / den
        cands.append((x0[0] + t * d[0], x0[1] + t * d[1]))
    for c in cands:
        if in_piece(p, c):
            return c
    return None


def _boundary_lines(p):
    k = p['kind']
    if k == 'poly':
        P = p['P']
        return [(P[i], P[(i + 1) % len(P)]) for i in range(len(P))]
    if k == 'cone':
        a = p['apex']
        return [(a, (a[0] + p['d1'][0], a[1] + p['d1'][1])), (a, (a[0] + p['d2'][0], a[1] + p['d2'][1]))]
    if k == 'segment':
        return [tuple(p['P'])]
    return []


def _cone_unbounded(q, s, p):
    """On a cone piece: is sup = +inf? Exact. Directions d = (1-w) d1 + w d2, w in [0, 1]."""
    d1, d2 = p['d1'], p['d2']
    # alpha(w) = dHd(d(w)) is a quadratic in w
    A0 = dHd(q, d1)
    A1 = dHd(q, d2)
    Am = dHd(q, ((d1[0] + d2[0]) / 2, (d1[1] + d2[1]) / 2))
    # alpha(w) = c0 + c1 w + c2 w^2 through (0, A0), (1/2, Am), (1, A1)
    c0 = A0
    c2 = 2 * A0 + 2 * A1 - 4 * Am
    c1 = A1 - A0 - c2
    ws = [Fr(0), Fr(1)]
    if c2 != 0:
        wc = -c1 / (2 * c2)
        if 0 < wc < 1:
            ws.append(wc)
    mins = min(c0 + c1 * w + c2 * w * w for w in ws)
    if mins < 0:
        return True
    # alpha >= 0 on the cone; where alpha == 0 the objective is linear along d (H d = 0 needed
    # for the slope to be independent of the base point; the generator only builds cones with
    # PSD H when alpha touches 0, so H d = 0 there)
    for w in ws:
        if c0 + c1 * w + c2 * w * w == 0:
            d = ((1 - w) * d1[0] + w * d2[0], (1 - w) * d1[1] + w * d2[1])
            g = qgrad(q, p['apex'])
            be = (s[0] - g[0]) * d[0] + (s[1] - g[1]) * d[1]
            if be > 0:
                return True
    if c0 == 0 and c1 == 0 and c2 == 0:
        for d in (d1, d2):
            g = qgrad(q, p['apex'])
            if (s[0] - g[0]) * d[0] + (s[1] - g[1]) * d[1] > 0:
                return True
    return False


def sup_piece(p, s):
    q = p['q']
    k = p['kind']
    if k == 'point':
        return obj(q, s, p['P'][0])
    if k == 'segment':
        return _sup_segment(q, s, p['P'][0], p['P'][1])
    if k == 'plane':
        st = _stationary(q, s)
        H11, H12, H22 = q[0], q[1], q[2]
        psd = H11 >= 0 and H22 >= 0 and H11 * H22 - H12 * H12 >= 0
        if not psd or st is None:
            return INF
        x = st[1] if st[0] in ('point', 'line') else (Fr(0), Fr(0))
        return obj(q, s, x)
    best = None
    if k == 'poly':
        P = p['P']
        for i in range(len(P)):
            v = _sup_segment(q, s, P[i], P[(i + 1) % len(P)])
            best = v if best is None or v > best else best
    elif k == 'cone':
        if _cone_unbounded(q, s, p):
            return INF
        for d in (p['d1'], p['d2']):
            v = _sup_ray(q, s, p['apex'], d)
            if v is INF:
                return INF
            best = v if best is None or v > best else best
    st = _stationary(q, s)
    if st is not None:
        if st[0] == 'point':
            if in_piece(p, st[1]):
                best = max(best, obj(q, s, st[1]))
        elif st[0] == 'line':
            c = _line_hits_piece(p, st[1], st[2])
            if c is not None:
                best = max(best, obj(q, s, c))
        # 'plane': objective constant (q affine with slope s): any vertex already gives it
    return best


def conj_exact(fx, s):
    s = (Fr(s[0]), Fr(s[1]))
    best = None
    for p in fx['pieces']:
        v = sup_piece(p, s)
        if v is INF:
            return INF
        best = v if best is None or v > best else best
    return best


# ---------------------------------------------------------------------------------------------
# the biconjugate: a certified bracket
# ---------------------------------------------------------------------------------------------
def domain_samples(fx, n=24, box=None):
    """Points of dom f (exact Fractions): vertices, edge points, a grid inside each piece."""
    pts = set()
    for p in fx['pieces']:
        k = p['kind']
        if k in ('point',):
            pts.add(tuple(p['P'][0]))
            continue
        if k == 'segment':
            a, b = p['P']
            for i in range(n + 1):
                t = Fr(i, n)
                pts.add((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
            continue
        if k == 'poly':
            P = p['P']
            xs = [v[0] for v in P]
            ys = [v[1] for v in P]
            lo, hi = (min(xs), min(ys)), (max(xs), max(ys))
            for i in range(len(P)):
                a, b = P[i], P[(i + 1) % len(P)]
                for j in range(n + 1):
                    t = Fr(j, n)
                    pts.add((a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])))
        else:
            lo, hi = box
        for i in range(n + 1):
            for j in range(n + 1):
                x = (lo[0] + (hi[0] - lo[0]) * Fr(i, n), lo[1] + (hi[1] - lo[1]) * Fr(j, n))
                if in_piece(p, x):
                    pts.add(x)
    return sorted(pts)


def hull_upper(fx, samples, x):
    """min sum l_i f(x_i) s.t. sum l_i x_i = x, sum l_i = 1, l >= 0 (an upper bound on f**(x);
    +inf if x is not in the hull of the samples). Doubles (scipy LP)."""
    import numpy as np
    from scipy.optimize import linprog
    X = np.array([[float(a), float(b)] for (a, b) in samples])
    fv = np.array([float(f_at(fx, p)) for p in samples])
    A = np.vstack([X.T, np.ones(len(samples))])
    b = np.array([float(x[0]), float(x[1]), 1.0])
    r = linprog(fv, A_eq=A, b_eq=b, bounds=(0, None), method='highs')
    if r.status != 0:
        return float('inf')
    return float(r.fun)

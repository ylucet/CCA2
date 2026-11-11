"""alg_face_below -- is the algebraic face `lo` below one of the faces `hi` everywhere in the
part of a box inside a region? (port of algFaceBelow.m)

Per box: sample both faces, fit a rational quadratic SEPARATOR between them, and ask algRowFails
to certify "lo above it" and "hi below it" both fail; failing that, ask it to certify the direct
comparison; failing that, subdivide. True is a proof; False means "not proved".
"""
import math

import numpy as np

from . import rat_q
from .alg_branch import algBranch
from .alg_row_fails import algRowFails
from .numeric import linspace, roundHalfAway


def algFaceBelow(lo, hi, box, inBox, depth=6, factors=None, corners=None, own=None):
    factors = factors or []
    corners = corners or []
    own = [{'FP': o['FP'], 'FB': o['FB'], 'num': o['num'], 'den': o['den'], 'side': o['side'],
            'W': None} for o in (own or [])]
    if isinstance(hi, dict):
        hi = [hi]
    ctx = {'inBox': inBox, 'factors': factors, 'corners': corners, 'own': own}
    return _rec(lo, hi, list(box), ctx, depth)


def _rec(lo, hi, b, ctx, depth):
    if not ctx['inBox'](b):
        return True                                  # misses the region: nothing to prove here
    us = linspace(b[0], b[1], 5)
    vs = linspace(b[2], b[3], 5)
    X = [[u, v] for u in us for v in vs]             # meshgrid, column-major: u outer, v inner
    zl = algBranch(lo['FP'], X, lo['FB'])
    n = len(hi)
    ZH = []
    cnt = []
    for j in range(n):
        zh = algBranch(hi[j]['FP'], X, hi[j]['FB'])
        ZH.append(zh)
        cnt.append(sum(1 for a, c in zip(zl, zh)
                       if math.isfinite(a) and math.isfinite(c) and a < c))
    order = sorted(range(n), key=lambda j: -cnt[j])
    inner = 6 if depth == 0 else 3
    for j in order:
        zh = ZH[j]
        ok = [math.isfinite(a) and math.isfinite(c) and a < c for a, c in zip(zl, zh)]
        if sum(ok) < 6:
            break
        num = _separator([x for x, o in zip(X, ok) if o], [z for z, o in zip(zl, ok) if o],
                         [z for z, o in zip(zh, ok) if o])
        if num is None:
            continue
        rl = {'FP': lo['FP'], 'FB': lo['FB'], 'num': num[0:10], 'den': num[10], 'side': -1,
              'W': None}
        rh = {'FP': hi[j]['FP'], 'FB': hi[j]['FB'], 'num': num[0:10], 'den': num[10], 'side': 1,
              'W': None}
        if algRowFails([rl] + ctx['own'], b, ctx['inBox'], inner, ctx['factors'],
                       ctx['corners']) and \
           algRowFails([rh] + ctx['own'], b, ctx['inBox'], inner, ctx['factors'],
                       ctx['corners']):
            return True
    inPiece = [ctx['inBox']([x[0], x[0], x[1], x[1]]) for x in X]
    for o in ctx['own']:
        if not any(inPiece):
            break
        zo = algBranch(o['FP'], X, o['FB'])
        qo = rat_q.evalFace(o['num'], o['den'], X)
        inPiece = [p and math.isfinite(z) and o['side'] * (q - z) >= -1e-9
                   for p, z, q in zip(inPiece, zo, qo)]
    for j in order:
        zh = ZH[j]
        both = [math.isfinite(a) and math.isfinite(c) for a, c in zip(zl, zh)]
        bad = any(bo and ip and a > c + 1e-9 * max(1.0, abs(c))
                  for bo, ip, a, c in zip(both, inPiece, zl, zh))
        if sum(both) < 6 or bad:
            continue
        fr = {'FP': lo['FP'], 'FB': lo['FB'], 'num': None, 'den': None, 'side': None,
              'W': {'FP': hi[j]['FP'], 'FB': hi[j]['FB']}}
        if algRowFails([fr] + ctx['own'], b, ctx['inBox'], inner, ctx['factors'], ctx['corners']):
            return True
    if depth == 0:
        return False
    xm = (b[0] + b[1]) / 2
    ym = (b[2] + b[3]) / 2
    return (_rec(lo, hi, [b[0], xm, b[2], ym], ctx, depth - 1)
            and _rec(lo, hi, [xm, b[1], b[2], ym], ctx, depth - 1)
            and _rec(lo, hi, [b[0], xm, ym, b[3]], ctx, depth - 1)
            and _rec(lo, hi, [xm, b[1], ym, b[3]], ctx, depth - 1))


def _separator(X, zl, zh):
    """A quadratic with DYADIC coefficients strictly between zl and zh at every sample: the
    11-vector [face numerator (10), denominator], or None."""
    V = np.array([[x[0] ** 2, x[0] * x[1], x[1] ** 2, x[0], x[1], 1.0] for x in X])
    rhs = np.array([(a + c) / 2 for a, c in zip(zl, zh)])
    c, *_ = np.linalg.lstsq(V, rhs, rcond=None)
    for k in range(41):
        r = [roundHalfAway(float(v) * 2 ** k) for v in c]
        if any(abs(v) >= 2 ** 50 for v in r):
            return None
        q = V @ np.array(r, dtype=float) / 2 ** k
        if all(qi > a and qi < b for qi, a, b in zip(q, zl, zh)):
            return [0, 0, 0, 0, 2 * r[0], r[1], 2 * r[2], r[3], r[4], r[5], 2 ** k]
    return None

"""bridge_arc_cut -- the rows that cut a bridge face to the ARC of the tie lying in the cell
pair, for a STRAIGHT tie (port of bridgeArcCut.m).

A straight tie's arc ends at rational points, and the cut boundary at each end is the segment
joining the two gradients there: one rational line per end, sided by the witness. A genuinely
CONIC tie is refused by name (bridgeArcCut:conicTie) -- its cut boundaries are irrational
segments, and the caller then keeps the cut as a predicate. Returns (cells, info).
"""
from . import rat_q
from .errors import error


def bridgeArcCut(tie, P, N1, m1, N2, m2, d, wit, C1, C2, kp):
    tie = list(tie)
    if any(x != 0 for x in tie[0:3]):
        error('bridgeArcCut:conicTie',
              'the tie is a genuine CONIC, so each cut boundary is a SEGMENT of an irrational '
              'line and only the conjugate PAIR is rational. Two constructions for that have '
              'been built and MEASURED WRONG; read bridgeArcCut.m\'s header and DECISIONS.md '
              '2026-09-17 before attempting a third.')
    D, E, F = tie[3], tie[4], tie[5]
    if D == 0 and E == 0:
        error('bridgeArcCut:degenerateTie',
              'the tie has no linear part either, so it names no curve at all.')
    w = [E, -D]
    s0n = [-F * D, -F * E]
    nd = D * D + E * E
    lo = None
    hi = None                                       # rational bounds as [num, den], den > 0
    for i, r in enumerate(P):
        a = r[0:2]
        c = r[2]
        A = a[0] * w[0] + a[1] * w[1]
        B = a[0] * s0n[0] + a[1] * s0n[1] + c * nd
        if A == 0:
            if B < 0:
                error('bridgeArcCut:emptyArc',
                      'facet %d excludes the ENTIRE tie line, so this pair contains no tie point '
                      'and carries no bridge face. The pair classification and this routine '
                      'disagree, which is a caller defect.' % (i + 1))
            continue
        b = [-B, A] if A > 0 else [B, -A]           # keep the denominator positive
        if A > 0:
            if lo is None or b[0] * lo[1] > lo[0] * b[1]:
                lo = b
        else:
            if hi is None or b[0] * hi[1] < hi[0] * b[1]:
                hi = b
    if lo is not None and hi is not None and hi[0] * lo[1] < lo[0] * hi[1]:
        error('bridgeArcCut:emptyArc',
              'the facets leave no tie point at all: the arc\'s bounds cross.')
    rows = []
    for tau in (lo, hi):
        if tau is None:
            continue                                # the arc runs to infinity on that side
        row = _endpointRow(tau, w, s0n, nd, N1, m1, N2, m2, d, wit)
        if row is None:
            continue                                # the segment degenerates to a point
        rows.append(row)
    return [rows], {'kind': 'line', 'nCut': len(rows), 'nCells': 1}


def _endpointRow(tau, w, s0n, nd, N1, m1, N2, m2, d, wit):
    tn, td = tau
    seN = [s0n[0] * td + tn * w[0], s0n[1] * td + tn * w[1]]       # s_e = seN / seD
    seD = nd * td
    P1n = [a + seD * b for a, b in zip(rat_q.mv(N1, seN), m1)]     # grad g1 = P1n/(d*seD)
    P2n = [a + seD * b for a, b in zip(rat_q.mv(N2, seN), m2)]
    dr = [P2n[0] - P1n[0], P2n[1] - P1n[1]]
    if dr[0] == 0 and dr[1] == 0:
        return None                                 # the two gradients agree: no segment
    n = [-dr[1], dr[0]]
    c = rat_q.conic([0, 0, 0, d * seD * n[0], d * seD * n[1], -(n[0] * P1n[0] + n[1] * P1n[1])])
    xn, xd = wit[0:2], wit[2]
    v = (c[0] * xn[0] ** 2 + c[1] * xn[0] * xn[1] + c[2] * xn[1] ** 2
         + c[3] * xn[0] * xd + c[4] * xn[1] * xd + c[5] * xd ** 2)
    s = rat_q.sign(v)
    if s == 0:
        error('bridgeArcCut:witnessOnCut',
              'the witness lies exactly ON a cut boundary, so it names no side. A witness in '
              'the INTERIOR of the arc\'s sweep decides it; this one was chosen before the cut '
              'was known.')
    return c + [s, 0]                               # branch 0: a straight cut line IS the conic

"""envelope_cell_pairs -- classify every pair of cells of two exact conjugates (port of
envelopeCellPairs.m).

For cells i of gA and j of gB that meet in 2-D (inside the dual domain D), the tie conic
{gA_i = gB_j} either keeps one sign on the overlap ('poly': one of the two wins there) or changes
sign ('bridge': the envelope needs a bridge face). Everything is integer arithmetic.
Returns (pairs, stats); each pair is {'kind', 'ij', 'P', 'q', 'win'}.
"""
from . import rat_q
from .errors import error
from .qua_con import QuaCon


def envelopeCellPairs(gA, gB, D=None):
    D = [list(r) for r in (D or [])]
    for g in (gA, gB):
        if not isinstance(g, QuaCon):
            error('envelopeCellPairs:type',
                  'both arguments must be QuaCon (the exact conjugate of one piece); got %s. '
                  'The classification below is integer arithmetic on fN/fD and EcQ, which only '
                  'an exact conjugate carries.' % type(g).__name__)
    pairs = []
    stats = {'nPairs': 0, 'nMeet': 0, 'nSame': 0, 'nPoly': 0, 'nBridge': 0, 'nSkipCurved': 0}
    for i in range(1, gA.nf + 1):
        RA, thinA = _cellHalfPlanes(gA, i)
        if thinA:
            continue
        for j in range(1, gB.nf + 1):
            RB, thinB = _cellHalfPlanes(gB, j)
            if thinB:
                continue
            stats['nPairs'] += 1
            P = RA + RB + D
            if not rat_q.feasible2(P, True):
                continue                                  # FILTER 1: meet in 2-D
            stats['nMeet'] += 1
            q = _tieConic(gA, i, gB, j)
            if all(x == 0 for x in q):
                stats['nSame'] += 1
                stats['nPoly'] += 1
                pairs.append(_mk('poly', [i, j], P, q, 'A'))
                continue
            sgn, strict = rat_q.quadSignOnPoly(q, P)      # FILTER 2: does the tie change sign?
            if sgn != 0 and strict:
                stats['nPoly'] += 1
                pairs.append(_mk('poly', [i, j], P, q, 'A' if sgn > 0 else 'B'))
            else:
                stats['nBridge'] += 1
                pairs.append(_mk('bridge', [i, j], P, q, ''))
    return pairs, stats


def _mk(kind, ij, P, q, win):
    return {'kind': kind, 'ij': ij, 'P': [list(r) for r in P], 'q': list(q), 'win': win}


def _cellHalfPlanes(g, k):
    rows = []
    thin = False
    for e, side, _ in g.FC[k - 1]:
        cc = g.EcQ[e - 1]
        if any(x != 0 for x in cc[0:3]):
            error('envelopeCellPairs:curvedCell',
                  'cell %d of a %s has a genuinely CONIC edge (edge %d, quadratic part '
                  '[%d %d %d]). The conjugate of a CONVEX piece is polyhedral -- every cell '
                  'boundary is the affine condition that a KKT multiplier changes sign -- so a '
                  'curved cell means this is not one.'
                  % (k, type(g).__name__, e, cc[0], cc[1], cc[2]))
        if side == 0:
            thin = True                                   # a THIN face: on the curve
            continue
        rows.append([side * x for x in cc[3:6]])
    return rows, thin


def _tieConic(gA, i, gB, j):
    D = [a * gB.fD[j - 1] - b * gA.fD[i - 1] for a, b in zip(gA.fN[i - 1], gB.fN[j - 1])]
    return [D[4], 2 * D[5], D[6], 2 * D[7], 2 * D[8], 2 * D[9]]

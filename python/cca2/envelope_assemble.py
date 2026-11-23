"""envelope_assemble -- the convex envelope contributed by ONE pair of exact conjugates, as a
list of candidates (port of envelopeAssemble.m).

For every pair of cells (i of gA, j of gB) that meet, the pair's value is
sup over the overlap of <s,x> - max(gA_i, gB_j)(s). Where one function wins on the whole overlap
that is envelopeCellRegions' quadratic candidates; where the tie crosses the overlap it is the
two sides cut at the tie, the BRIDGE face along the tie (bridgeFaceAlg) and the CORNER faces
where the tie leaves the cell. Returns (env, info):
  env['quad']   : [{'num', 'den', 'con', 'pair'}]   rational candidates with their regions
  env['bridge'] : [face dicts]                      bridge / corner faces, each with 'pairIndex'
"""
from . import rat_q
from .bridge_face_alg import bridgeFaceAlg
from .corner_face_alg import cornerFaceAlg
from .envelope_cell_pairs import envelopeCellPairs
from .envelope_cell_regions import envelopeCellRegions
from .errors import error, CCA2Error
from .qua_con import QuaCon
from .rat_q import sign


def envelopeAssemble(gA, gB, D=None):
    for h in (gA, gB):
        if not isinstance(h, QuaCon):
            error('envelopeAssemble:type',
                  'both arguments must be QuaCon; got %s.' % type(h).__name__)
    D = [list(r) for r in (D or [])]                # dom f* to restrict to
    pairs, stats = envelopeCellPairs(gA, gB, D)

    info = {'nPairs': len(pairs), 'nPoly': 0, 'nTangential': 0, 'nStraight': 0, 'nCurved': 0,
            'stats': stats}
    env = {'quad': [], 'bridge': []}

    def collect(cells, k):
        for c in cells:
            env['quad'].append({'num': c['num'], 'den': c['den'], 'con': c['con'], 'pair': k})

    def corners(nA, dA, nB, dB, pr, k):
        for cf in cornerFaceAlg(nA, dA, pr['P'], pr['q'], +1) + \
                cornerFaceAlg(nB, dB, pr['P'], pr['q'], -1):
            cf['pairIndex'] = k
            env['bridge'].append(cf)

    for k0, pr in enumerate(pairs):
        k = k0 + 1
        nA, dA = gA.fN[pr['ij'][0] - 1], gA.fD[pr['ij'][0] - 1]
        nB, dB = gB.fN[pr['ij'][1] - 1], gB.fD[pr['ij'][1] - 1]

        if pr['kind'] == 'poly':
            info['nPoly'] += 1
            collect(_cellsFor(nA, dA, pr['P']) if pr['win'] == 'A'
                    else _cellsFor(nB, dB, pr['P']), k)
            continue

        # a 'bridge' pair whose tie only TOUCHES the overlap is still one winner
        sgn, strict = rat_q.quadSignOnPoly(pr['q'], pr['P'])
        if sgn != 0 and not strict:
            info['nTangential'] += 1
            collect(_cellsFor(nA, dA, pr['P']) if sgn > 0 else _cellsFor(nB, dB, pr['P']), k)
            continue

        # a STRAIGHT tie: an ordinary facet of each side polygon. The pair's value is the MAX of
        # two separate sups, one per side polygon, and each is +inf where none of its own
        # candidates covers x (an unbounded side). So the B side is its own group: under one
        # group index a finite A-side candidate would cover a point where the B side's sup is
        # +inf, and co f came out finite outside conv(dom f) (review 2026-11-19 #9a).
        if all(x == 0 for x in pr['q'][0:3]):
            info['nStraight'] += 1
            collect(_cellsFor(nA, dA, pr['P'] + [pr['q'][3:6]]), k)
            collect(_cellsFor(nB, dB, pr['P'] + [[-x for x in pr['q'][3:6]]]), len(pairs) + k)
            continue

        try:
            face, _ = bridgeFaceAlg(nA, dA, nB, dB, pr['P'])
        except CCA2Error as ME:
            if ME.identifier == 'bridgeFaceAlg:tangentialArc':
                info['nTangential'] += 1
                s2 = _sideWinner(nA, dA, nB, dB, pr['P'])
                collect(_cellsFor(nA, dA, pr['P']) if s2 > 0 else _cellsFor(nB, dB, pr['P']), k)
                continue
            if ME.identifier in ('bridgeFaceAlg:straightTie', 'bridgeFaceRat:noSeparateFace'):
                info['nStraight'] += 1
                collect(_sideCells(nA, dA, pr['P'], pr['q'], +1), k)
                collect(_sideCells(nB, dB, pr['P'], pr['q'], -1), k)
                _straightTieGuard(nA, nB, pr['q'])
                corners(nA, dA, nB, dB, pr, k)
                continue
            raise

        info['nCurved'] += 1
        face['pairIndex'] = k      # a bridge face belongs to a PAIR, and eval groups by pair
        env['bridge'].append(face)
        collect(_sideCells(nA, dA, pr['P'], pr['q'], +1), k)
        collect(_sideCells(nB, dB, pr['P'], pr['q'], -1), k)
        corners(nA, dA, nB, dB, pr, k)

    info['nQuadFaces'] = len(env['quad'])
    info['nCornerFaces'] = sum(1 for f in env['bridge'] if f['kind'] == 'corner')
    info['nBridgeFaces'] = len(env['bridge']) - info['nCornerFaces']
    return env, info


def _straightTieGuard(nA, nB, q):
    """A tie on IRRATIONAL parallel lines along which a face curves is refused by name."""
    import math
    A, B, C, D, E, F = q
    if A != 0:
        disc = 4 * D * D - 16 * A * F
        w = [-B, 2 * A]
    else:
        disc = 4 * E * E - 16 * C * F
        w = [1, 0]
    if disc <= 0 or math.isqrt(disc) ** 2 == disc:
        return                                      # double line, or rational lines
    for n in (nA, nB):
        N = [[n[4], n[5]], [n[5], n[6]]]
        if rat_q.quad(w, N, w) != 0:
            error('envelopeAssemble:irrationalStraightFacet',
                  'a cell pair ties on irrational parallel lines and a face CURVES along them, '
                  'so its sup may sit inside the tie segment -- a facet candidate on an '
                  'irrational line, which is not built.')


def _cellsFor(n, d, P):
    N = [[n[4], n[5]], [n[5], n[6]]]
    m = [n[7], n[8]]
    faces = envelopeCellRegions(N, m, n[9], d, P)
    return [{'num': f['fN'], 'den': f['fD'], 'con': _rowsToConics(f['R'])} for f in faces]


def _sideCells(n, d, P, tie, side):
    """One side's candidates, each cut to the part of its region whose MAXIMISER lies on that
    side of the tie (the tie pulled back through the maximiser map)."""
    N = [[n[4], n[5]], [n[5], n[6]]]
    m = [n[7], n[8]]
    faces = envelopeCellRegions(N, m, n[9], d, P)
    cells = []
    for f in faces:
        S = f['S']
        c = _pullback(tie, S['Mn'], S['tn'], S['d'])
        con = _rowsToConics(f['R'])
        if not all(x == 0 for x in c):
            con = con + [rat_q.conic(c) + [side * _firstNonZero(c), 0]]
        cells.append({'num': f['fN'], 'den': f['fD'], 'con': con})
    return cells


def _sideWinner(nA, dA, nB, dB, P):
    c = rat_q.diffConic(nA, dA, nB, dB)
    s, _ = rat_q.quadSignOnPoly(c, P)
    return s if s != 0 else 1      # identical faces: either is the winner, and they agree


def _rowsToConics(R):
    return [rat_q.conic([0, 0, 0] + list(r)) + [_firstNonZero(r), 0] for r in R]


def _firstNonZero(c):
    v = next((x for x in c if x != 0), None)
    return 1 if v is None else sign(v)


def _pullback(q, Mn, tn, dd):
    """The conic q(s) composed with s = (Mn*x + tn)/dd, cleared of dd^2."""
    L1 = [Mn[0][0], Mn[0][1], tn[0]]
    L2 = [Mn[1][0], Mn[1][1], tn[1]]

    def mulAff(f, g):
        return [f[0] * g[0], f[0] * g[1] + f[1] * g[0], f[1] * g[1],
                f[0] * g[2] + f[2] * g[0], f[1] * g[2] + f[2] * g[1], f[2] * g[2]]
    a = mulAff(L1, L1)
    b = mulAff(L1, L2)
    c = mulAff(L2, L2)
    out = [q[0] * x + q[1] * y + q[2] * z for x, y, z in zip(a, b, c)]
    lin = [dd * (q[3] * x + q[4] * y) for x, y in zip(L1, L2)]
    out[3] += lin[0]
    out[4] += lin[1]
    out[5] += lin[2] + q[5] * dd * dd
    return out

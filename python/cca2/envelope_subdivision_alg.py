"""envelope_subdivision_alg -- subdivide an assembled envelope when a BRIDGE (or irrational
corner) face competes with the quadratic faces (port of envelopeSubdivisionAlg.m).

The conic part is as in envelope_subdivision: tabulate the region rows as curves, build their
arrangement (cellInterior decides the live sign patterns; a curve that is irrelevant in a cell
is SKIPPED there rather than splitting it), find an interior point per cell. What is new: in a
cell an implicit face may be present, and the boundary between it and a quadratic face q is
NOT a conic. It is recorded as an ALGEBRAIC ROW "side * (q - value) >= 0" on the face, which
AlgCrv evaluates; rows and pieces that provably hold nowhere are dropped (algRowFails /
algFaceBelow certificates).

Returns (subdivision, info): subdivision is a list of faces
    {'kind': 'quad' | 'bridge', 'num', 'den', 'FP', 'FB', 'con', 'alg'}
that assembleAlgCrv turns into an AlgCrv.

Several atom pairs (a list of envelopes) are subdivided one by one and OVERLAID: co f is the MIN
over pairs, so in each overlay cell every choice of one face per pair is compared, the smallest
winning; ties between two implicit faces become rows of table 2 ("value_1 >= value_2").
"""
import math

import numpy as np

from . import rat_q
from .alg_branch import algBranch
from .alg_face_below import algFaceBelow
from .alg_row_fails import algRowFails
from .bridge_member import bridgeMember
from .cell_interior import cellInterior
from .cell_sample_point import cellSamplePoint, cellWitness
from .conic_meet import _solve2
from .envelope_rational_corners import envelopeRationalCorners
from .errors import error, CCA2Error
from .numeric import linspace, roundHalfAway
from .poly_sign_on import polySignOn
from .rat_q import sign

PAIR_STRIDE = 10 ** 6


def _pairsOf(p):
    return p if isinstance(p, list) else [p]


def envelopeSubdivisionAlg(env, dom=None, cuts=None, opts=None):
    dom = [list(r) for r in (dom or [])]
    cuts = [list(r) for r in (cuts or [])]
    opts = dict(opts or {})
    opts.setdefault('overlayBudget', 3000)
    opts.setdefault('maxCells', 3000)
    envs = env if isinstance(env, list) else [env]
    for e in envs:
        if not isinstance(e, dict) or 'quad' not in e or 'bridge' not in e:
            error('envelopeSubdivisionAlg:shape', 'env must be what envelopeAssemble returns.')
    if len(envs) > 1:
        subs = []
        infos = []
        for e in envs:
            sb, inf = envelopeSubdivisionAlg(e, dom, cuts, opts)
            subs.append(sb)
            infos.append(inf)
        return _minOverlay(subs, infos, opts['overlayBudget'])

    env = _combinePairs(envs)
    corners, bridges = _splitBridge(env['bridge'])
    C, req, pairOf = _conicTable(env['quad'], corners)
    nQuad = len(env['quad'])
    allPairs = sorted(set([p for po in pairOf for p in po] + [f['pairIndex'] for f in bridges]))

    C, cEnds = _addCornerCrossings(C, corners)
    C, bEnds = _addBridgeCurves(C, bridges)    # else a bridge's membership can change in a cell
    C, lEnds = _splitLineCurves(C)
    for cu in cuts:
        _addKey(C, rat_q.conic([0, 0, 0] + cu) + [0])
    owners = _curveOwners(C, req, nQuad, bridges, cEnds, bEnds, lEnds)

    cells = [[]]
    if dom:
        base = []
        for drow in dom:
            raw = [0, 0, 0] + drow
            key = rat_q.conic(raw) + [0]
            hit = _addKey(C, key)
            base.append([hit, sign(next(x for x in raw if x != 0))])
        cells = [base]
        while len(owners) < len(C):
            owners.append([[]])
    skips = [[]]
    nSkip = 0
    maxCells = opts['maxCells']
    for j in range(1, len(C) + 1):
        nxt = []
        nxtSk = []
        for cell, sk in zip(cells, skips):
            if any(r[0] == j for r in cell):
                nxt.append(cell)
                nxtSk.append(sk)
                continue
            liveP = _liveCell(C, cell + [[j, 1]])
            liveM = _liveCell(C, cell + [[j, -1]])
            if liveP and liveM and not _relevantHere(C, owners[j - 1], cell):
                nxt.append(cell)
                nxtSk.append(sk + [j])
                nSkip += 1
            elif liveP and liveM:
                nxt.append(cell + [[j, 1]])
                nxt.append(cell + [[j, -1]])
                nxtSk.append(sk)
                nxtSk.append(sk)
            elif liveP or liveM:
                nxt.append(cell)                    # the curve misses this cell
                nxtSk.append(sk)
        cells = nxt
        skips = nxtSk
        if len(cells) > maxCells:
            error('envelopeSubdivisionAlg:arrangementSize',
                  'the arrangement passed %d cells at curve %d of %d; envelopeEval answers for '
                  'the face list.' % (maxCells, j, len(C)))

    subdivision = []
    nAlg = nUncovered = nOverlap = nNoPoint = nBelow = 0
    cellCon = []
    faceCell = []
    for rows, sk in zip(cells, skips):
        pt, okPt = _interiorPoint(C, rows)
        if not okPt:
            nUncovered += 1
            nNoPoint += 1
            continue
        cellCon.append(_rowsToCon(C, rows))
        app = _applicable(req, rows, C, pt)
        here = [p for i in app for p in pairOf[i - 1]]
        bHere = []
        Q, Pw = _votePoints(C, rows, pt)
        rowsOut = _fullRows(C, rows, Q, sk)         # emitted; the exact tests keep the sparse rows
        for b, f in enumerate(bridges):
            if sum(bridgeMember(f['pair'], Q)[0]) * 2 > len(Q):
                here.append(f['pairIndex'])
                bHere.append(b)
        if not all(p in here for p in allPairs):
            nUncovered += 1
            continue
        quadApp = [i for i in app if i <= nQuad]
        irr = [{'FP': bridges[b]['FP'], 'FB': bridges[b]['FB']} for b in bHere]
        for j in [i for i in app if i > nQuad]:
            fc = _cornerValueFace(corners[j - nQuad - 1], pt)
            if fc is not None:
                irr.append(fc)
        if not quadApp and not irr:
            nUncovered += 1
            continue
        if len(irr) > 1:
            irr, nDrop = _dropBelow(irr, C, rows, Pw, env['quad'], quadApp)
            nBelow += nDrop
        if len(irr) > 1 and not _dominantIrrational(irr, C, rows, pt, env['quad'], quadApp):
            nOverlap += 1
        sub = _winnerFaces(env['quad'], quadApp, irr, rows, C, pt, Pw, rowsOut)
        for s in sub:
            nAlg += len(s['alg'])
            subdivision.append(s)
            faceCell.append(len(cellCon))
    info = {'cells': len(cells), 'faces': len(subdivision), 'algRows': nAlg,
            'refused': nUncovered, 'noPoint': nNoPoint, 'overlaps': nOverlap, 'below': nBelow,
            'skipped': nSkip, 'cellCon': cellCon, 'faceCell': faceCell}
    return subdivision, info


# ------------------------------------------------------------------------------------------------
def _winnerFaces(quad, quadApp, irr, rows, C, pt, Q, rowsOut):
    out = []
    base = _rowsToCon(C, rowsOut)
    if not irr:
        for j in quadApp:
            extra = _quadOver(quad, quadApp, j)
            if not _pieceLive(C, rows, extra, Q):
                continue
            out.append({'kind': 'quad', 'num': quad[j - 1]['num'], 'den': quad[j - 1]['den'],
                        'FP': None, 'FB': 0, 'con': base + extra, 'alg': []})
        return out
    for f in irr:
        algRows = [{'FP': f['FP'], 'FB': f['FB'], 'num': quad[j - 1]['num'],
                    'den': quad[j - 1]['den'], 'side': -1} for j in quadApp]
        out.append({'kind': 'bridge', 'num': [], 'den': 0, 'FP': f['FP'], 'FB': f['FB'],
                    'con': base, 'alg': algRows})
    for j in quadApp:
        extra = _quadOver(quad, quadApp, j)
        if not _pieceLive(C, rows, extra, Q):
            continue
        algRows = [{'FP': f['FP'], 'FB': f['FB'], 'num': quad[j - 1]['num'],
                    'den': quad[j - 1]['den'], 'side': 1} for f in irr]
        if _algEmpty(C, rowsOut, extra, algRows, Q):
            continue                                # the face's own rows hold nowhere
        out.append({'kind': 'quad', 'num': quad[j - 1]['num'], 'den': quad[j - 1]['den'],
                    'FP': None, 'FB': 0, 'con': base + extra, 'alg': algRows})
    return out


def _evalCurve(C, j, P):
    c = C[j - 1][0:6]
    if C[j - 1][6] != 0:
        return rat_q.pairLineEval(c, C[j - 1][6], P)
    return rat_q.evalConic(c, P)


def _fullRows(C, rows, P, skipped):
    """The cell's rows plus, for each curve neither split on nor skipped, the side the sample
    points vote for."""
    rows = [list(r) for r in rows]
    have = set(r[0] for r in rows) | set(skipped)
    for j in range(1, len(C) + 1):
        if j in have:
            continue
        v = [sign(x) for x in _evalCurve(C, j, P)]
        v = [x for x in v if x != 0]
        if not v:
            continue                                # on the curve throughout: unrecorded
        rows.append([j, sign(sum(v))])
    return rows


def _votePoints(C, rows, pt):
    Q = [pt]
    _, ok, P = cellSamplePoint(C, rows, _straightRows(C, rows), 1e3)
    if ok and len(P) >= 4:
        idx = [roundHalfAway(t) for t in linspace(1, len(P), 4)]
        Q = [pt] + [P[i - 1] for i in idx]
    return Q, [pt] + P


def _liveCell(C, rows):
    tf, info = cellInterior(C, rows)
    if not tf:
        return False
    hasBranch = any(C[r[0] - 1][6] != 0 for r in rows)
    hasCurved = any(C[r[0] - 1][6] == 0 and any(x != 0 for x in C[r[0] - 1][0:3]) for r in rows)
    if info['path'] == 'unsureTie' or (hasBranch and hasCurved):
        _, tf, _ = cellSamplePoint(C, rows, _straightRows(C, rows), 1e3)
        if not tf:                                  # the grid can miss a small cell (#9c)
            _, tf = cellWitness(C, rows)
    return tf


def _algEmpty(C, rows, conRows, algRows, Q):
    """Do the piece's algebraic rows hold NOWHERE in it? Sample first (a witness means no),
    then ask for a certificate."""
    if not algRows:
        return False
    Qk = []
    for q in Q:
        if all(r[6] * rat_q.evalConic(r[0:6], [q])[0] > 0 for r in conRows):
            Qk.append(q)
    holds = [True] * len(Qk)
    nFail = []
    for a in algRows:
        z = algBranch(a['FP'], Qk, a['FB'])
        q = rat_q.evalFace(a['num'], a['den'], Qk)
        ok = [math.isfinite(zi) and a['side'] * (qi - zi) > 1e-9 * max(1.0, abs(zi))
              for zi, qi in zip(z, q)]
        holds = [h and o for h, o in zip(holds, ok)]
        nFail.append(sum(1 for o in ok if not o))
    if any(holds):
        return False                                # a witness: not empty
    for t, a in enumerate(algRows):
        if not Qk or nFail[t] < len(Qk):
            continue                                # a lone row must fail at every point
        if _rowFails(C, rows, conRows, [a]):
            return True
    if len(algRows) > 1 or not Qk:
        o = sorted(range(len(algRows)), key=lambda t: -nFail[t])
        return _rowFails(C, rows, conRows, [algRows[t] for t in o])
    return False


_rowFailsCache = {}


def _polyKey(P):
    return tuple(sorted(P.items())) if P else ()


def _rowFails(C, rows, conRows, a):
    key = (tuple(sorted(tuple(C[r[0] - 1]) + (r[1],) for r in rows)),
           tuple(sorted(tuple(r) for r in conRows)),
           tuple((_polyKey(x['FP']), tuple(x['num']), x['den'], x['FB'], x['side']) for x in a))
    if key in _rowFailsCache:
        return _rowFailsCache[key]
    ok = False
    box, inBox, factors, corners = _regionBox(C, rows, conRows)
    if box is not None:
        try:
            ok = bool(algRowFails(a, box, inBox, 6, factors, corners))
        except CCA2Error:
            ok = False
    _rowFailsCache[key] = ok
    return ok


def _regionBox(C, rows, K):
    """The bounding box of the cell's straight rows (and branch lines), a predicate inBox(b)
    'the box b may meet the cell', the cell's rows as bivariate polynomials, and its corners."""
    H = _straightRows(C, rows)
    A = [list(h) for h in H] + [[1, 0, 1e3], [-1, 0, 1e3], [0, 1, 1e3], [0, -1, 1e3]]
    A += _branchLines(C, rows, K)
    V = []
    for i in range(len(A)):
        for j in range(i + 1, len(A)):
            Mx = [[A[i][0], A[i][1]], [A[j][0], A[j][1]]]
            if abs(Mx[0][0] * Mx[1][1] - Mx[0][1] * Mx[1][0]) < 1e-12:
                continue
            x = _solve2(Mx, [-A[i][2], -A[j][2]])
            lim = -1e-9 * max(1.0, math.sqrt(x[0] * x[0] + x[1] * x[1]))
            if all(a[0] * x[0] + a[1] * x[1] + a[2] >= lim for a in A):
                V.append(x)
    if len(V) < 3:
        return None, None, [], []
    box = [min(v[0] for v in V), max(v[0] for v in V), min(v[1] for v in V), max(v[1] for v in V)]
    rowPolys = []
    for r in rows:
        if C[r[0] - 1][6] != 0:
            continue
        rowPolys.append(_conicPoly(C[r[0] - 1][0:6], r[1]))
    for r in K:
        if len(r) >= 8 and r[7] != 0:
            continue
        rowPolys.append(_conicPoly(r[0:6], r[6]))
    lines = []
    for r in rows:
        if C[r[0] - 1][6] != 0:
            lines.append(list(C[r[0] - 1][0:6]) + [r[1], C[r[0] - 1][6]])
    lines += [list(r[0:8]) for r in K if len(r) >= 8 and r[7] != 0]

    def inBox(b):
        for P in rowPolys:
            if polySignOn(P, b, lambda _b: True, 0) == -1:
                return False
        return not _lineMisses(lines, b)
    return box, inBox, rowPolys, V


def _branchLines(C, rows, K):
    L = []
    E = [[0, 0], [1, 0], [0, 1]]

    def line(c, br, sd):
        v = [sd * x for x in rat_q.pairLineEval(c, br, E)]
        return [v[1] - v[0], v[2] - v[0], v[0]]
    for r in rows:
        br = C[r[0] - 1][6]
        if br == 0:
            continue
        L.append(line(C[r[0] - 1][0:6], br, r[1]))
    for r in K:
        if len(r) >= 8 and r[7] != 0:
            L.append(line(r[0:6], r[7], r[6]))
    out = []
    for l in L:
        if not all(math.isfinite(x) for x in l):
            continue
        nr = math.hypot(l[0], l[1])
        if nr > 0:
            out.append([x / nr for x in l])
    return out


def _lineMisses(lines, b):
    X = [[b[0], b[2]], [b[1], b[2]], [b[0], b[3]], [b[1], b[3]]]
    for r in lines:
        v = [r[6] * x for x in rat_q.pairLineEval(r[0:6], r[7], X)]
        if max(v) < -1e-7 * max([1.0] + [abs(x) for x in v]):
            return True
    return False


def _conicPoly(c, s=1):
    """The conic as a 3x3 coefficient matrix P[i][j] of x^i y^j, times s."""
    return [[s * c[5], s * c[4], s * c[2]],
            [s * c[3], s * c[1], 0],
            [s * c[0], 0, 0]]


def _pieceLive(C, rows, extra, Q):
    k0 = [all(x == 0 for x in r[0:5]) for r in extra]
    if any(k0):
        if any(r[6] * r[5] < 0 for r, z in zip(extra, k0) if z):
            return False
        extra = [r for r, z in zip(extra, k0) if not z]
    if not extra:
        return True
    for q in Q:
        if all(r[6] * rat_q.evalConic(r[0:6], [q])[0] > 1e-9 * max(abs(x) for x in r[0:6])
               for r in extra):
            return True
    fails = [sum(1 for q in Q if r[6] * rat_q.evalConic(r[0:6], [q])[0] <= 0) for r in extra]
    o = sorted(range(len(extra)), key=lambda t: -fails[t])
    for t in o[:3]:
        r = extra[t]
        if not _liveCell(C + [list(r[0:6]) + [r[7]]], rows + [[len(C) + 1, r[6]]]):
            return False
    return _liveCell(C + [list(r[0:6]) + [r[7]] for r in extra],
                     rows + [[len(C) + 1 + t, r[6]] for t, r in enumerate(extra)])


def _quadOver(quad, quadApp, j):
    """The rows 'face j is at least face i' for every other applicable quadratic face i."""
    extra = []
    for i in quadApp:
        if i == j:
            continue
        c, sgn = _diffConic(quad[j - 1], quad[i - 1])
        if c is None:
            continue
        extra.append(c + [sgn, 0])
    return extra


def _diffConic(fw, fj):
    num = [a * fj['den'] - b * fw['den'] for a, b in zip(fw['num'], fj['num'])]
    raw = [num[4], 2 * num[5], num[6], 2 * num[7], 2 * num[8], 2 * num[9]]
    if all(x == 0 for x in raw) or any(x != 0 for x in num[0:4]):
        return None, 0
    nz = next(x for x in raw if x != 0)
    return rat_q.conic(raw), sign(nz)


def _addKey(C, key):
    """Append key to C unless present; returns its 1-based index."""
    key = list(key)
    try:
        return C.index(key) + 1
    except ValueError:
        C.append(key)
        return len(C)


def _addCornerCrossings(C, corners):
    """For an IRRATIONAL corner, the rational line across which its two conjugate values swap."""
    ends = []
    for j, c in enumerate(corners):
        if 'sigma' not in c or c['sigma']['root'] == 0:
            continue                                # rational corner: one root
        n = c['src']['n']
        d = c['src']['d']
        N = [[n[4], n[5]], [n[5], n[6]]]
        m = [n[7], n[8]]
        al, be = c['sigma']['quad'][0], c['sigma']['quad'][1]
        s0n = c['sigma']['s0n']
        w = c['sigma']['w']
        na = c['sigma']['na']
        T1 = (2 * al * rat_q.quad(s0n, N, w) - be * rat_q.quad(w, N, w)
              + 2 * na * al * rat_q.dot(m, w))
        E1 = [2 * al * na * d * w[0], 2 * al * na * d * w[1], -T1]
        if E1[0] == 0 and E1[1] == 0:
            continue
        key = rat_q.conic([0, 0, 0] + E1) + [0]
        _addKey(C, key)
        ends.append((key, j))
    return C, ends


def _cornerValueFace(c, pt):
    n = c['src']['n']
    d = c['src']['d']
    s = c['s']
    z = pt[0] * s[0] + pt[1] * s[1] - ((s[0] * (n[4] * s[0] + n[5] * s[1])
                                        + s[1] * (n[5] * s[0] + n[6] * s[1])) / 2
                                       + n[7] * s[0] + n[8] * s[1] + n[9]) / d
    _, info = algBranch(c['FP'], [pt], -1, wantInfo=True)
    rr = info['allReal'][0]
    if not rr:
        return None
    gaps = [abs(r - z) for r in rr]
    gap = min(gaps)
    j = gaps.index(gap)
    if gap > 1e-6 * max(1.0, abs(z)):
        return None
    return {'FP': c['FP'], 'FB': (j + 1) - len(rr) - 1}


def _dropBelow(irr, C, rows, P, quad, quadApp):
    """Drop the implicit faces that are provably below another one on the whole cell."""
    nDrop = 0
    keep = [True] * len(irr)
    for i in range(1, len(irr)):
        for j in range(i):
            if keep[j] and irr[i]['FB'] == irr[j]['FB'] and irr[i]['FP'] == irr[j]['FP']:
                keep[i] = False
                nDrop += 1
                break
    irr = [f for f, k in zip(irr, keep) if k]
    n = len(irr)
    if n < 2:
        return irr, nDrop
    V = []
    for f in irr:
        v = algBranch(f['FP'], P, f['FB'])
        V.append([x if math.isfinite(x) else -math.inf for x in v])
    top = []
    for t in range(len(P)):
        col = [V[i][t] for i in range(n)]
        top.append(col.index(max(col)))
    wins = [top.count(i) for i in range(n)]
    order = sorted(range(n), key=lambda i: -wins[i])
    above = order[:min(2, n)]
    box, inBox, factors, corners = _regionBox(C, rows, [])
    if box is None:
        return irr, nDrop
    own = []
    for f in irr:
        own.append([{'FP': f['FP'], 'FB': f['FB'], 'num': quad[j - 1]['num'],
                     'den': quad[j - 1]['den'], 'side': -1} for j in quadApp])
    keep = [True] * n
    for i in range(n):
        for j in above:
            if j == i or not keep[j]:
                continue
            if all(a <= b + 1e-9 * max(1.0, abs(b)) for a, b in zip(V[i], V[j])) and \
                    algFaceBelow(irr[i], irr[j], box, inBox, 5, factors, corners, own[i]):
                keep[i] = False
                nDrop += 1
                break
    for i in range(n):
        if not (keep[i] and wins[i] == 0):
            continue
        others = [j for j in range(n) if keep[j] and j != i]
        if not others:
            break
        hi = [{'FP': irr[j]['FP'], 'FB': irr[j]['FB']} for j in others]
        if algFaceBelow({'FP': irr[i]['FP'], 'FB': irr[i]['FB']}, hi, box, inBox, 5, factors,
                        corners, own[i]):
            keep[i] = False
            nDrop += 1
    return [f for f, k in zip(irr, keep) if k], nDrop


def _dominantIrrational(irr, C, rows, pt, quad, quadApp):
    _, okS, P = cellSamplePoint(C, rows, _straightRows(C, rows), 1e3)
    if not okS:
        P = [pt]
    q = [-math.inf] * len(P)
    for j in quadApp:
        v = rat_q.evalFace(quad[j - 1]['num'], quad[j - 1]['den'], P)
        q = [max(a, b) for a, b in zip(q, v)]
    above = [0] * len(P)
    for f in irr:
        v = algBranch(f['FP'], P, f['FB'])
        for t in range(len(P)):
            if math.isfinite(v[t]) and v[t] > q[t] + 1e-12 * max(1.0, abs(q[t])):
                above[t] += 1
    return not any(a > 1 for a in above)


def _straightRows(C, rows):
    H = []
    for r in rows:
        c = C[r[0] - 1][0:6]
        if C[r[0] - 1][6] == 0 and all(x == 0 for x in c[0:3]):
            h = [r[1] * x for x in c[3:6]]
            nr = math.hypot(h[0], h[1])
            if nr > 0:
                h = [x / nr for x in h]
            H.append(h)
    return H


def _combinePairs(envs):
    env = {'quad': [], 'bridge': []}
    for g, e0 in enumerate(envs):
        e = envelopeRationalCorners(e0)
        off = g * PAIR_STRIDE
        for q in e['quad']:
            pr = q['pair']
            pr = [p + off for p in pr] if isinstance(pr, list) else pr + off
            env['quad'].append({'num': q['num'], 'den': q['den'], 'con': q['con'], 'pair': pr})
        for b in e['bridge']:
            b = dict(b)
            b['pairIndex'] = b['pairIndex'] + off
            env['bridge'].append(b)
    return env


def _splitBridge(bl):
    corners = []
    bridges = []
    for f in bl:
        if f['kind'] == 'corner' and f.get('con'):
            corners.append(f)
        elif 'FP' in f and 'pair' in f and 'cells' in f:      # cells: its outer bound
            bridges.append(f)
        else:
            error('envelopeSubdivisionAlg:faceKind',
                  'a face of kind %s has neither rows nor an implicit value.' % f['kind'])
    return corners, bridges


def _conicTable(faces, corners):
    C = []
    n = len(faces) + len(corners)
    req = [None] * n
    pairOf = [None] * n
    for i in range(n):
        if i < len(faces):
            con = faces[i]['con']
            pairOf[i] = _pairsOf(faces[i]['pair'])
        else:
            con = corners[i - len(faces)]['con']
            pairOf[i] = [corners[i - len(faces)]['pairIndex']]
        r = []
        dead = False
        for row in con:
            c = row[0:6]
            if all(x == 0 for x in c[0:5]):
                if not _constHolds(c[5], row[6]):
                    dead = True
                    break
                continue
            br = row[7] if len(row) >= 8 else 0
            hit = _addKey(C, list(c) + [br])
            r.append([hit, row[6]])
        req[i] = None if dead else r
    return C, req, pairOf


def _constHolds(c6, side):
    return (side > 0 and c6 >= 0) or (side < 0 and c6 <= 0) or (side == 0 and c6 == 0)


def _applicable(req, rows, C, pt):
    side = {}
    for r in rows:
        side.setdefault(r[0], r[1])
    app = []
    for i, rq in enumerate(req):
        if rq is None:
            continue
        ok = True
        for t in rq:
            s = side.get(t[0])
            if s is None:
                s = sign(_evalCurve(C, t[0], [pt])[0])       # a skipped curve: read it at pt
            if s != t[1]:
                ok = False
                break
        if ok:
            app.append(i + 1)
    return app


def _addBridgeCurves(C, bridges):
    ends = []
    for b, f in enumerate(bridges):
        for con in f['cells']:
            for row in con:
                c = row[0:6]
                if all(x == 0 for x in c[0:5]):
                    continue
                br = row[7] if len(row) >= 8 else 0
                _addKey(C, list(c) + [br])
        pr = f['pair']
        if 'N1' not in pr:
            continue                                # a rational face: exact rows
        for k in (1, 2):
            Nk, mk = (pr['N1'], pr['m1']) if k == 1 else (pr['N2'], pr['m2'])
            if all(x == 0 for r in Nk for x in r) or rat_q.detExact(Nk) != 0:
                continue
            col = next(j for j in range(2) if Nk[0][j] != 0 or Nk[1][j] != 0)
            u = [Nk[0][col], Nk[1][col]]            # the range of a rank-1 Nk
            w = [-u[1], u[0]]
            raw = [0, 0, 0, pr['d'] * w[0], pr['d'] * w[1], -(w[0] * mk[0] + w[1] * mk[1])]
            key = rat_q.conic(raw) + [0]
            _addKey(C, key)
            ends.append((key, b))
    return C, ends


def _splitLineCurves(C):
    """A conic that is a DOUBLE line (c = k * L^2) also enters as the line L itself."""
    ends = []
    for j in range(len(C)):
        c = C[j][0:6]
        if C[j][6] != 0 or all(x == 0 for x in c[0:3]):
            continue
        L = _squareRoot(c)
        if L is not None:
            key = rat_q.conic([0, 0, 0] + L) + [0]
            _addKey(C, key)
            ends.append((key, j + 1))
    return C, ends


def _curveOwners(C, req, nQuad, bridges, cEnds, bEnds, lEnds):
    """owners[j]: for each face whose region uses curve j, that face's OTHER rows. A curve
    matters in a cell only where one of its owners can still hold."""
    n = len(C)
    owners = [[] for _ in range(n)]
    for r in req:
        if not r:
            continue
        for t in r:
            owners[t[0] - 1].append([x for x in r if x[0] != t[0]])
    for f in bridges:
        for r in _bridgeBlocks(C, f):
            for u in r:
                owners[u[0] - 1].append([x for x in r if x[0] != u[0]])
    for key, cj in cEnds:
        j = C.index(key)
        r = req[nQuad + cj]
        if r is None:
            continue
        owners[j].append([x for x in r if x[0] != j + 1])
    for key, b in bEnds:
        j = C.index(key)
        for blk in _bridgeBlocks(C, bridges[b]):
            owners[j].append([x for x in blk if x[0] != j + 1])
    for key, src in lEnds:                          # after the faces: the double line's are there
        j = C.index(key)
        s = owners[src - 1] or [[]]
        owners[j] = owners[j] + s
    for j in range(n):
        if not owners[j]:
            owners[j] = [[]]
            continue
        seen = set()
        uniq = []
        for r in owners[j]:
            k = tuple(sorted(tuple(x) for x in r))
            if k in seen:
                continue
            seen.add(k)
            uniq.append(r)
        owners[j] = uniq
    return owners


def _bridgeBlocks(C, bridge):
    blocks = []
    for con in bridge['cells']:
        r = []
        for row in con:
            c = row[0:6]
            if all(x == 0 for x in c[0:5]):
                continue
            br = row[7] if len(row) >= 8 else 0
            r.append([C.index(list(c) + [br]) + 1, row[6]])
        blocks.append(r)
    return blocks


def _relevantHere(C, owners, rows):
    side = {}
    for r in rows:
        side.setdefault(r[0], r[1])
    for r in owners:
        if not r:
            return True
        clash = False
        rest = []
        for t in r:
            s = side.get(t[0])
            if s is None:
                rest.append(t)
                continue
            if s != t[1]:
                clash = True
                break
        if clash:
            continue
        if not rest or _liveCell(C, rows + rest):
            return True
    return False


def _squareRoot(c):
    """L = [p, q, r] with c = k * (p x + q y + r)^2 exactly, or None."""
    a, b, cc = c[0], c[1], c[2]
    if b * b != 4 * a * cc:
        return None                                 # quadratic part not a square
    g = math.gcd(math.gcd(abs(a), abs(cc)), abs(b))
    if g == 0:
        return None
    A, B, Cq = a // g, b // g, cc // g
    sg = 1
    if A < 0 or Cq < 0:
        A, B, Cq, sg = -A, -B, -Cq, -1
    p = math.isqrt(A)
    q = math.isqrt(Cq)
    if p * p != A or q * q != Cq:
        return None                                 # not a square over the integers
    if B < 0:
        q = -q                                      # (p x + q y)^2 has cross term 2pq
    if 2 * p * q != B:
        return None
    k = sg * g
    from fractions import Fraction
    r = Fraction(c[3], 2 * k * p) if p != 0 else Fraction(c[4], 2 * k * q)
    if r.denominator != 1:
        r2 = 2 * r
        if r2.denominator != 1:
            return None
        Lc = [2 * p, 2 * q, int(r2)]
        kk = Fraction(k, 4)
    else:
        Lc = [p, q, int(r)]
        kk = Fraction(k)
    want = [kk * Lc[0] ** 2, kk * 2 * Lc[0] * Lc[1], kk * Lc[1] ** 2, kk * 2 * Lc[0] * Lc[2],
            kk * 2 * Lc[1] * Lc[2], kk * Lc[2] ** 2]
    if all(w == x for w, x in zip(want, c)):
        return Lc
    return None


def _rowsToCon(C, rows):
    return [list(C[r[0] - 1][0:6]) + [r[1], C[r[0] - 1][6]] for r in rows]


def _holdsAll(C, rows, p):
    for r in rows:
        if r[1] * _evalCurve(C, r[0], [p])[0] <= 0:
            return False
    return True


def _interiorPoint(C, rows):
    """A point strictly inside the cell, on a coarse dyadic grid when possible. (pt, ok)."""
    H = []
    for r in rows:
        c = C[r[0] - 1][0:6]
        if C[r[0] - 1][6] == 0 and all(x == 0 for x in c[0:3]):
            H.append([r[1] * x for x in c[3:6]])
    R = 1e3
    if not H:
        cand = [0.0, 0.0]
    else:
        Hn = []
        for h in H:
            nrm = math.hypot(h[0], h[1])
            Hn.append([x / nrm for x in h] if nrm > 0 else [float(x) for x in h])
        H = Hn
        M = np.array([[h[0], h[1], -1.0] for h in H]
                     + [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], dtype=float)
        b = np.array([-h[2] for h in H] + [-R, -R, -R, -R], dtype=float)
        best = -math.inf
        cand = [0.0, 0.0]
        k = M.shape[0]
        for i in range(k - 2):
            for j in range(i + 1, k - 1):
                for l in range(j + 1, k):
                    S = M[[i, j, l], :]
                    if abs(np.linalg.det(S)) < 1e-12:
                        continue
                    z = np.linalg.solve(S, b[[i, j, l]])
                    if np.all(M @ z >= b - 1e-9) and z[2] > best:
                        best = float(z[2])
                        cand = [float(z[0]), float(z[1])]
        if best <= 0:
            return [0.0, 0.0], False
    for pas in (1, 2):
        for g in (8, 16, 64, 256, 1024):
            p = [roundHalfAway(cand[0] * g) / g, roundHalfAway(cand[1] * g) / g]
            if _holdsAll(C, rows, p):
                return p, True
        if _holdsAll(C, rows, cand):
            return cand, True
        if pas == 1:                                # a curved row excludes the centre
            cand, found, _ = cellSamplePoint(C, rows, H, R)
            if not found:
                cand, found = cellWitness(C, rows)  # the grid can miss a small cell (#9c)
            if not found:
                return [0.0, 0.0], False
    return [0.0, 0.0], False


# ---- the overlay of several pairs: a MIN of maxes ------------------------------------------------
def _keyRows(keys, con):
    rows = []
    for r in con:
        rows.append([_addKey(keys, list(r[0:6]) + [r[7]]), r[6]])
    return rows


def _mergeRows(a, b):
    side = {}
    for r in a:
        side.setdefault(r[0], r[1])
    out = [list(r) for r in a]
    for t in b:
        s = side.get(t[0])
        if s is None:
            out.append(list(t))
        elif s != t[1]:
            return None
    return out


def _uniqueStable(rows):
    seen = set()
    out = []
    for r in rows:
        k = tuple(r)
        if k not in seen:
            seen.add(k)
            out.append(list(r))
    return out


class _Samples:
    """Sample points of one overlay cell, with each face's branch values computed once."""

    def __init__(self, X):
        self.X = X
        self.Z = {}

    def branch(self, FP, FB, idx):
        key = (_polyKey(FP), FB)
        e = self.Z.setdefault(key, {})
        need = [i for i in idx if i not in e]
        if need:
            for i, v in zip(need, algBranch(FP, [self.X[i] for i in need], FB)):
                e[i] = v
        return [e[i] for i in idx]


def _conVal(r, X):
    if r[7] != 0:
        return rat_q.pairLineEval(r[0:6], r[7], X)
    return rat_q.evalConic(r[0:6], X)


def _witnessed(con, alg, ii, PP):
    """Is there a sample point strictly inside every row of the piece?"""
    X = PP.X
    if not X:
        return False
    idx = list(range(len(X)))
    for r in con:
        v = _conVal(r, [X[i] for i in idx])
        lim = 1e-9 * max(1, max(abs(x) for x in r[0:6]))
        idx = [i for i, vi in zip(idx, v) if r[6] * vi > lim]
        if not idx:
            return False
    K = 96
    if len(idx) > K:
        pick = [roundHalfAway(t) for t in linspace(1, len(idx), K)]
        idx = sorted(set(idx[p - 1] for p in pick))
    for a in alg:
        z = PP.branch(a['FP'], a['FB'], idx)
        q = rat_q.evalFace(a['num'], a['den'], [X[i] for i in idx])
        idx = [i for i, zi, qi in zip(idx, z, q)
               if math.isfinite(zi) and a['side'] * (qi - zi) > 1e-9]
        if not idx:
            return False
    for t in ii:
        z1 = PP.branch(t['FP1'], t['FB1'], idx)
        z2 = PP.branch(t['FP2'], t['FB2'], idx)
        idx = [i for i, a1, a2 in zip(idx, z1, z2)
               if math.isfinite(a1) and math.isfinite(a2) and a1 - a2 > 1e-9]
        if not idx:
            return False
    return bool(idx)


def _overlayLive(keys, rows, con, alg, ii, PP):
    if _witnessed(con, alg, ii, PP):
        return True
    k0 = [all(x == 0 for x in r[0:5]) for r in con]
    if any(r[6] * r[5] < 0 for r, z in zip(con, k0) if z):
        return False
    K = [list(k) for k in keys]
    r = _keyRows(K, [c for c, z in zip(con, k0) if not z])
    allRows = [list(t) for t in sorted(set(tuple(t) for t in rows + r))]
    if not _liveCell(K, allRows):
        return False
    if not alg:
        return True
    plain = [c[0:7] + [0] for c, z in zip(con, k0) if not z and c[7] == 0]
    return not _algEmpty(keys, rows, plain, alg, PP.X)


def _tieEmpty(keys, rows, con, ii, PP):
    """Does some tie row 'value_1 >= value_2' hold NOWHERE in the piece (so the piece is empty)?"""
    if not ii:
        return False
    K = [r for r in con if not all(x == 0 for x in r[0:5])]
    X = PP.X
    idx = list(range(len(X)))
    for r in K:
        v = _conVal(r, [X[i] for i in idx])
        lim = -1e-9 * max(1, max(abs(x) for x in r[0:6]))
        idx = [i for i, vi in zip(idx, v) if r[6] * vi >= lim]
    box = None
    for t in ii:
        if idx:
            z1 = PP.branch(t['FP1'], t['FB1'], idx)
            z2 = PP.branch(t['FP2'], t['FB2'], idx)
            if any(math.isfinite(a) and math.isfinite(b) and a > b for a, b in zip(z1, z2)):
                continue
        if box is None:
            box, inBox, factors, corners = _regionBox(keys, rows, K)
            if box is None:
                return False
        lo = {'FP': t['FP1'], 'FB': t['FB1']}
        hi = {'FP': t['FP2'], 'FB': t['FB2']}
        try:
            tf = algFaceBelow(lo, hi, box, inBox, 5, factors, corners)
        except CCA2Error:
            tf = False
        if tf:
            return True
    return False


def _leRows(f, o, con, alg, ii):
    """Append the rows saying 'face f is at most face o'."""
    if f['kind'] == 'quad' and o['kind'] == 'quad':
        c, sgn = _diffConic(o, f)                  # o >= f
        if c is not None:
            con = con + [c + [sgn, 0]]
    elif f['kind'] == 'quad':                      # quad f <= irrational o
        alg = alg + [{'FP': o['FP'], 'FB': o['FB'], 'num': f['num'], 'den': f['den'], 'side': -1}]
    elif o['kind'] == 'quad':                      # irrational f <= quad o
        alg = alg + [{'FP': f['FP'], 'FB': f['FB'], 'num': o['num'], 'den': o['den'], 'side': 1}]
    elif not (f['FP'] == o['FP'] and f['FB'] == o['FB']):
        ii = ii + [{'FP1': o['FP'], 'FB1': o['FB'], 'FP2': f['FP'], 'FB2': f['FB']}]
    return con, alg, ii


def _minOverlay(subs, infos, budget):
    nP = len(subs)
    keys = []
    cellRows = []
    covered = []
    facesOf = []
    for g in range(nP):
        cc = infos[g]['cellCon']
        cellRows.append([_keyRows(keys, c) for c in cc])
        fo = [[i for i, fc in enumerate(infos[g]['faceCell']) if fc == k + 1]
              for k in range(len(cc))]
        facesOf.append(fo)
        covered.append([bool(x) for x in fo])

    # ---- overlay the cells, one pair at a time ------------------------------------------------
    tup = [[k] for k in range(len(cellRows[0]))]
    trow = [list(r) for r in cellRows[0]]
    for g in range(1, nP):
        nt = []
        nr = []
        for a in range(len(trow)):
            for k in range(len(cellRows[g])):
                r = _mergeRows(trow[a], cellRows[g][k])
                if r is None or not _liveCell(keys, r):
                    continue
                nt.append(tup[a] + [k])
                nr.append(r)
        tup = nt
        trow = nr

    bound = 0
    for t in tup:
        G = [g for g in range(nP) if covered[g][t[g]]]
        pr = 1
        for g in G:
            pr *= len(facesOf[g][t[g]])
        bound += pr * len(G)
    if bound > budget:
        error('envelopeSubdivisionAlg:overlaySize',
              'the overlay of %d pairs would visit up to %d pieces over %d tuples, past the '
              'budget of %d; envelopeEval answers for the minimum.' % (nP, bound, len(tup), budget))

    out = []
    nUnsure = nII = nHole = nCombo = nTie = 0
    for t, rows in zip(tup, trow):
        G = [g for g in range(nP) if covered[g][t[g]]]
        if not G:
            nHole += 1
            continue
        base = []                                  # every pair's cell, holes included
        for g in range(nP):
            base += infos[g]['cellCon'][t[g]]
        _, _, PX = cellSamplePoint(keys, rows, _straightRows(keys, rows), 1e3, True)
        PP = _Samples(PX)
        combos = [{'pick': [], 'con': base, 'alg': []}]
        for g in G:
            nc = []
            for cb in combos:
                for k in facesOf[g][t[g]]:
                    fk = subs[g][k]
                    con = _uniqueStable(cb['con'] + [list(r[0:8]) for r in fk['con']])
                    alg = cb['alg'] + fk['alg']
                    if not _overlayLive(keys, rows, con, alg, [], PP):
                        continue
                    nc.append({'pick': cb['pick'] + [fk], 'con': con, 'alg': alg})
            combos = nc
        nCombo += len(combos)
        for cb in combos:
            pick = cb['pick']
            for w in range(len(G)):
                f = pick[w]
                con = cb['con']
                alg = cb['alg']
                ii = []
                for i in range(len(G)):
                    if i == w:
                        continue
                    con, alg, ii = _leRows(f, pick[i], con, alg, ii)
                con = _uniqueStable(con)
                if not _overlayLive(keys, rows, con, alg, ii, PP):
                    continue
                if not _witnessed(con, alg, ii, PP):
                    if _tieEmpty(keys, rows, con, ii, PP):
                        nTie += 1
                        continue
                    nUnsure += 1
                nII += len(ii)
                out.append({'kind': f['kind'], 'num': f['num'], 'den': f['den'], 'FP': f['FP'],
                            'FB': f['FB'], 'con': con, 'alg': alg, 'ii': ii})
    info = {'cells': len(tup), 'faces': len(out), 'algRows': sum(len(f['alg']) for f in out),
            'iiRows': nII, 'bound': bound, 'refused': nHole, 'unsure': nUnsure, 'tieEmpty': nTie,
            'combos': nCombo, 'pairCells': [x['cells'] for x in infos],
            'pairFaces': [x['faces'] for x in infos],
            'skipped': sum(x['skipped'] for x in infos)}
    return out, info

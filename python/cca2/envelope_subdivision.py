"""envelope_subdivision -- subdivide an assembled envelope (ONE atom pair) on a CONIC
subdivision, returning a QuaCon (port of envelopeSubdivision.m).

The candidates' region rows are tabulated as curves; the arrangement of those curves is built
cell by cell (cellInterior decides which sign patterns are live); in each cell the applicable
candidates are compared pairwise and the cell is refined along each difference conic, the MAX
winning. Refused by name -- so the caller can fall back -- when a bridge face competes
(envelopeSubdivision:bridgeCompetes), when a face's region is only a predicate
(:bridgeFaces), or when the arrangement is too large (:arrangementSize). Returns (g, info).
"""
import math

import numpy as np

from . import rat_q
from .assemble_qua_con_cells import assembleQuaConCells
from .bridge_member import bridgeMember
from .cell_interior import cellInterior
from .cell_sample_point import cellSamplePoint
from .envelope_eval import envelopeEval
from .envelope_rational_corners import envelopeRationalCorners
from .errors import error
from .qua_con import QuaCon
from .rat_q import sign


def _pairsOf(p):
    return p if isinstance(p, list) else [p]


def envelopeSubdivision(env, opts=None):
    opts = dict(opts or {})
    maxCurves = opts.get('maxCurves', 40)
    if not isinstance(env, dict) or 'quad' not in env or 'bridge' not in env:
        error('envelopeSubdivision:shape', 'env must be what envelopeAssemble returns.')

    env = envelopeRationalCorners(env)            # a rational corner is a value face
    corners = []
    bridges = []
    for f in env['bridge']:
        if f['kind'] == 'corner' and f.get('con'):
            corners.append(f)
        elif f.get('cells') and 'pair' in f:
            bridges.append(f)
        else:
            error('envelopeSubdivision:bridgeFaces',
                  'this envelope has a face whose region is a PREDICATE with no rows and no '
                  'outer bound (kind %s), so nothing places it in an arrangement. Keep using '
                  'envelopeEval, which answers for it.' % f['kind'])

    for f in bridges:
        if not f.get('wit'):
            continue
        x = [f['wit'][0] / f['wit'][2], f['wit'][1] / f['wit'][2]]
        if bridgeMember(f['pair'], [x])[0][0] and math.isfinite(envelopeEval(env, [x])[0][0]):
            error('envelopeSubdivision:bridgeCompetes',
                  'a bridge face owns part of the envelope\'s domain, and its value competes '
                  'with the quadratic faces there. The boundary between them is not a conic, so '
                  'this envelope cannot be subdivided on a conic subdivision -- '
                  'envelopeSubdivisionAlg can. DECISIONS.md 2026-09-19.')

    if not env['quad']:
        return QuaCon([], [], [], [], [], [], []), {'curves': 0, 'cells': 0, 'refined': 0,
                                                    'faces': 0}

    C, req, pairOf = _curveTable(env['quad'], corners)
    nQuad = len(env['quad'])
    C = _addBridgeCurves(C, bridges)              # curves only: a bridge's membership is a PREDICATE
    allPairs = sorted(set([p for po in pairOf for p in po] + [f['pairIndex'] for f in bridges]))

    if len(C) > maxCurves:
        error('envelopeSubdivision:arrangementSize',
              'the arrangement has %d curves (budget %d), and a global arrangement of that size '
              'does not finish. envelopeEval answers for the face list.' % (len(C), maxCurves))

    # ---- the arrangement, curve by curve -----------------------------------------------------
    cells = [[]]                                  # each cell: rows [curveIndex, side]
    for j in range(1, len(C) + 1):
        nxt = []
        for cell in cells:
            for sd in (1, -1):
                cand = cell + [[j, sd]]
                if cellInterior(C, cand)[0]:
                    nxt.append(cand)
        cells = nxt
    nArr = len(cells)

    out = []
    nRef = 0
    nUncovered = 0
    for cell in cells:
        app = _applicable(req, cell)
        pt, okPt = _interiorPoint(C, cell)
        here = [p for i in app for p in pairOf[i - 1]]
        bridgeHere = []
        if okPt:
            for b, f in enumerate(bridges):
                if bridgeMember(f['pair'], [pt])[0][0]:
                    here.append(f['pairIndex'])
                    bridgeHere.append(b)
        if not app and not bridgeHere:
            nUncovered += 1
            continue
        if not all(p in here for p in allPairs):
            nUncovered += 1
            continue
        if bridgeHere:
            error('envelopeSubdivision:bridgeCompetes',
                  'a bridge face owns part of a cell, and its value competes with the quadratic '
                  'faces there. The boundary between them is not a conic, so this envelope '
                  'cannot be subdivided on a conic subdivision -- envelopeEval answers for it. '
                  'DECISIONS.md 2026-09-19.')
        quadApp = [i for i in app if i <= nQuad]
        if not quadApp:
            error('envelopeSubdivision:algebraicFace',
                  'a cell is covered only by CORNER faces, whose value is affine with '
                  'quadratic-irrational coefficients and cannot be stored on a conic '
                  'subdivision.')
        sub, C = _winnerCells(C, cell, env['quad'], quadApp)
        nRef += len(sub)
        out += sub

    g = assembleQuaConCells(out)
    info = {'curves': len(C), 'cells': nArr, 'refined': nRef, 'faces': g.nf,
            'uncovered': nUncovered, 'corners': len(corners), 'bridges': len(bridges)}
    return g, info


def _branchOf(con, t):
    return con[t][7] if len(con[t]) >= 8 else 0


def _curveTable(faces, corners):
    """C: the distinct curves [conic(6), branch]; req[i]: the rows [curve, side] face i needs
    (None when a constant row of its region fails); pairOf[i]: the pairs face i covers."""
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
        for t in range(len(con)):
            c = con[t][0:6]
            if all(x == 0 for x in c[0:5]):
                if not _constHolds(c[5], con[t][6]):
                    dead = True
                    break
                continue
            key = list(c) + [_branchOf(con, t)]
            try:
                hit = C.index(key)
            except ValueError:
                C.append(key)
                hit = len(C) - 1
            r.append([hit + 1, con[t][6]])
        req[i] = None if dead else r
    return C, req, pairOf


def _constHolds(c6, side):
    return (side > 0 and c6 >= 0) or (side < 0 and c6 <= 0) or (side == 0 and c6 == 0)


def _addBridgeCurves(C, bridges):
    for f in bridges:
        for con in f['cells']:
            for r in range(len(con)):
                c = con[r][0:6]
                if all(x == 0 for x in c[0:5]):
                    continue
                key = list(c) + [_branchOf(con, r)]
                if key not in C:
                    C.append(key)
    return C


def _interiorPoint(C, rows):
    """A point strictly inside the cell (doubles): the Chebyshev centre of its straight rows,
    or a sampled point when a curved row excludes that centre. Returns (pt, ok)."""
    H = []
    for r in rows:
        c = C[r[0] - 1][0:6]
        br = C[r[0] - 1][6]
        if br == 0 and all(x == 0 for x in c[0:3]):
            H.append([r[1] * x for x in c[3:6]])
    R = 1e3
    if not H:
        cand = [0.0, 0.0]
    else:
        Hn = []
        for h in H:
            n = math.hypot(h[0], h[1])
            Hn.append([x / n for x in h] if n > 0 else [float(x) for x in h])
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
    for r in rows:
        c = C[r[0] - 1][0:6]
        br = C[r[0] - 1][6]
        if br != 0:
            v = rat_q.pairLineEval(c, br, [cand])[0]
        else:
            v = (c[0] * cand[0] ** 2 + c[1] * cand[0] * cand[1] + c[2] * cand[1] ** 2
                 + c[3] * cand[0] + c[4] * cand[1] + c[5])
        if r[1] * v <= 0:
            pt, ok, _ = cellSamplePoint(C, rows, H, R)   # a curved row excludes the centre
            return pt, ok
    return cand, True


def _applicable(req, rows):
    """The (1-based) faces whose every required row has the cell's sign."""
    side = {}
    for r in rows:
        side.setdefault(r[0], r[1])
    app = []
    for i, rq in enumerate(req):
        if rq is None:
            continue
        if all(side.get(t[0]) == t[1] for t in rq):
            app.append(i + 1)
    return app


def _winnerCells(C, rows, faces, app):
    """Refine one cell along the difference conics of its applicable faces; the MAX wins."""
    open_ = [{'rows': rows, 'win': app[0]}]
    for j in app[1:]:
        nxt = []
        for cur in open_:
            d, sgn = _diffConicOf(faces[j - 1], faces[cur['win'] - 1])   # d > 0 means j wins
            if d is None:
                nxt.append(cur)                     # identical functions
                continue
            key = list(d) + [0]
            try:
                idx = C.index(key) + 1
            except ValueError:
                C.append(key)
                idx = len(C)
            for sd in (1, -1):
                cand = {'rows': cur['rows'] + [[idx, sd * sgn]],
                        'win': j if sd > 0 else cur['win']}
                if cellInterior(C, cand['rows'])[0]:
                    nxt.append(cand)
        open_ = nxt
    out = []
    for o in open_:
        f = faces[o['win'] - 1]
        out.append({'num': f['num'], 'den': f['den'],
                    'con': [list(C[r[0] - 1][0:6]) + [r[1], C[r[0] - 1][6]] for r in o['rows']]})
    return out, C


def _diffConicOf(fj, fi):
    num = [a * fi['den'] - b * fj['den'] for a, b in zip(fj['num'], fi['num'])]
    raw = [num[4], 2 * num[5], num[6], 2 * num[7], 2 * num[8], 2 * num[9]]
    if all(x == 0 for x in raw) or any(x != 0 for x in num[0:4]):
        return None, 0
    nz = next(x for x in raw if x != 0)
    return rat_q.conic(raw), sign(nz)

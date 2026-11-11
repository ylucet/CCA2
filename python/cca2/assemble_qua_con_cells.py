"""assemble_qua_con_cells -- turn a list of cells (face function + H-form rows) into a QuaCon
(port of assembleQuaConCells.m).

A cell is a dict {'num': 10 ints, 'den': int, 'con': rows [conic(6), side] or
[conic(6), side, branch]}. Dead cells (no interior) are dropped, implied rows removed, adjacent
cells carrying one function merged, the conics tabulated, and each corner NAMED.
"""
from . import rat_q
from .cell_has_interior import cellHasInterior, satisfiesAllDouble
from .compact_cells import compactCells
from .conic_meet import conicMeet
from .errors import error, CCA2Error
from .qua_con import QuaCon


def _emptyQuaCon():
    return QuaCon([], [], [], [], [], [], [])


def assembleQuaConCells(cells):
    if not cells:
        return _emptyQuaCon()
    cells = [dict(c) for c in cells]
    for k, c in enumerate(cells):
        r = [list(x) for x in (c['con'] or [])]
        r = [x + [0] if len(x) == 7 else x for x in r]           # the pre-branch spelling
        if any(len(x) != 8 for x in r):
            error('assembleQuaConCells:conShape',
                  'cell %d has a `con` that is not k x 7 [conic sign] or k x 8 '
                  '[conic sign branch].' % (k + 1))
        c['con'] = r

    live = [True] * len(cells)
    for k, c in enumerate(cells):
        rows = c['con']
        hasBranch = any(r[7] != 0 for r in rows)
        if any(r[6] == 0 for r in rows):
            # a THIN cell: only a curve asked for on BOTH sides kills it
            if _bothSides(rows):
                live[k] = False
            continue
        keep = [True] * len(rows)
        for i, r in enumerate(rows):
            if r[7] != 0:
                continue                                         # a BRANCH row
            cs = rat_q.conicSign(r[0:6])
            if cs == 0:
                continue
            if cs == r[6]:
                keep[i] = False                                  # vacuous
            else:
                live[k] = False                                  # empty interior
                break
        if not live[k]:
            continue
        rows = [r for r, kp in zip(rows, keep) if kp]
        c['con'] = rows
        if not rows:
            continue                                             # no conditions: the whole plane
        if _bothSides(rows):
            live[k] = False
            continue
        isLine = [all(x == 0 for x in r[0:3]) and r[7] == 0 for r in rows]
        live[k] = rat_q.feasible2([[r[6] * x for x in r[3:6]]
                                   for r, l in zip(rows, isLine) if l], True)
        if live[k] and hasBranch:
            continue
        if live[k] and not all(isLine):
            live[k] = cellHasInterior([r[0:6] for r in rows],
                                      [[i + 1, r[6]] for i, r in enumerate(rows)])

    if not any(live):
        # dom is possibly a SINGLE POINT: keep the largest candidate there
        pt, den = _singlePointDomain(cells)
        if pt is not None:
            best = None
            bestD = None
            for c in cells:
                if not _cellHoldsAtPoint(c['con'], pt, den):
                    continue
                x = [(pt[0] / den, pt[1] / den)]
                if best is None or rat_q.evalFace(c['num'], c['den'], x)[0] > \
                        rat_q.evalFace(best, bestD, x)[0]:
                    best = c['num']
                    bestD = c['den']
            if best is not None:
                r1 = [den, 0, -pt[0]]
                r2 = [0, den, -pt[1]]
                con = [rat_q.conic([0, 0, 0] + r1) + [0, 0], rat_q.conic([0, 0, 0] + r2) + [0, 0]]
                cells = [{'num': best, 'den': bestD, 'con': con}]
                live = [True]

    cells = [c for c, l in zip(cells, live) if l]
    if not cells:
        return _emptyQuaCon()

    # ---- tabulate the conics ---------------------------------------------------------------
    EcQ = []
    index = {}
    FC = []
    for c in cells:
        rows = c['con']
        fc = []
        for r in rows:
            key = tuple(r[0:6])
            hit = index.get(key)
            if hit is None:
                EcQ.append(list(key))
                hit = len(EcQ)
                index[key] = hit
            fc.append([hit, r[6], r[7]])
        FC.append(fc)
    fN = []
    fD = []
    for c in cells:
        n, d = rat_q.canon(c['num'], c['den'])
        fN.append(n)
        fD.append(d)

    fN, fD, FC, _ = compactCells(fN, fD, FC, EcQ)

    used = [False] * len(EcQ)
    for fc in FC:
        for r in fc:
            used[r[0] - 1] = True
    renum = []
    acc = 0
    for u in used:
        acc += 1 if u else 0
        renum.append(acc)
    EcQ = [c for c, u in zip(EcQ, used) if u]
    for fc in FC:
        for r in fc:
            r[0] = renum[r[0] - 1]

    # ---- name the corners ------------------------------------------------------------------
    seen = set()
    Vname = []
    for rows in FC:
        for p1 in range(len(rows)):
            for p2 in range(p1 + 1, len(rows)):
                i = rows[p1][0]
                j = rows[p2][0]
                if i == j:
                    continue
                lo = min(i, j)
                hi = max(i, j)
                if (lo, hi) in seen:
                    continue
                ci = EcQ[i - 1]
                cj = EcQ[j - 1]
                if all(x == 0 for x in ci[0:3]) and all(x == 0 for x in cj[0:3]):
                    A = [[ci[3], ci[4]], [cj[3], cj[4]]]
                    if rat_q.detExact(A) == 0:
                        continue                                 # parallel: they meet nowhere
                    xn, xd = rat_q.solve2(A, [-ci[5], -cj[5]])
                    if not _cellHoldsAt(EcQ, rows, xn, xd):
                        continue
                    seen.add((lo, hi))
                    Vname.append([lo, hi, 1])
                else:
                    try:
                        P, info = conicMeet(ci, cj)
                    except CCA2Error as ME:         # only reachable under rat_q.EMULATE_2P53
                        if ME.identifier != 'ratQ:overflow':
                            raise
                        continue                    # MATLAB leaves such a corner unnamed
                    if info['degenerate']:
                        continue
                    for rt in range(len(P)):
                        if not satisfiesAllDouble(EcQ, rows, P[rt]):
                            continue
                        seen.add((lo, hi))
                        Vname.append([lo, hi, rt + 1])
                        break

    E, F = _buildIncidence(EcQ, FC, Vname)
    return QuaCon(Vname, E, EcQ, fN, fD, F, FC)


def _bothSides(rows):
    """Is some curve (conic + branch) asked for on both its strict sides?"""
    sides = {}
    for r in rows:
        key = tuple(r[0:6]) + (r[7],)
        s = sides.setdefault(key, [False, False])
        if r[6] > 0:
            s[0] = True
        elif r[6] < 0:
            s[1] = True
    return any(s[0] and s[1] for s in sides.values())


def _conicAt(c, xn, xd):
    return (c[0] * xn[0] ** 2 + c[1] * xn[0] * xn[1] + c[2] * xn[1] ** 2
            + c[3] * xn[0] * xd + c[4] * xn[1] * xd + c[5] * xd ** 2)


def _cellHoldsAt(EcQ, rows, xn, xd):
    for r in rows:
        val = _conicAt(EcQ[r[0] - 1], xn, xd)
        if r[1] == 0:
            if val != 0:
                return False                    # an EQUALITY row: the corner must be ON it
        elif r[1] * val < 0:
            return False
    return True


def _singlePointDomain(cells):
    seen = None
    for c in cells:
        rows = c['con']
        lin = [r for r in rows if all(x == 0 for x in r[0:3])]
        for a in range(len(lin)):
            for b in range(a + 1, len(lin)):
                A = [[lin[a][3], lin[a][4]], [lin[b][3], lin[b][4]]]
                if rat_q.detExact(A) == 0:
                    continue
                x, xd = rat_q.solve2(A, [-lin[a][5], -lin[b][5]])
                if not _cellHoldsAtPoint(rows, x, xd):
                    continue
                if seen is None:
                    seen = (x[0], x[1], xd)
                elif seen[0] * xd != x[0] * seen[2] or seen[1] * xd != x[1] * seen[2]:
                    return None, None           # two distinct points: not a single-point domain
    if seen is None:
        return None, None
    return [seen[0], seen[1]], seen[2]


def _cellHoldsAtPoint(rows, xn, xd):
    for r in rows:
        v = _conicAt(r[0:6], xn, xd)
        if r[6] == 0:
            if v != 0:
                return False
        elif r[6] * v < 0:
            return False
    return True


def _buildIncidence(EcQ, FC, Vname):
    ne = len(EcQ)
    E = [[0, 0, 0] for _ in range(ne)]
    F = [[0, 0] for _ in range(ne)]
    for j in range(1, ne + 1):
        vs = [i + 1 for i, v in enumerate(Vname) if v[0] == j or v[1] == j]
        if len(vs) >= 2:
            E[j - 1] = [vs[0], vs[1], 1]        # two named corners: a segment between them
        elif len(vs) == 1:
            E[j - 1] = [vs[0], 0, 0]            # one named corner: unbounded from it
        for k, fc in enumerate(FC):
            r = next((x for x in fc if x[0] == j), None)
            if r is None:
                continue
            if r[1] > 0:
                if F[j - 1][0] == 0:
                    F[j - 1][0] = k + 1
            elif r[1] < 0:
                if F[j - 1][1] == 0:
                    F[j - 1][1] = k + 1
    return E, F

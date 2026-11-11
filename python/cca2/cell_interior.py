"""cell_interior -- interior test for a cell whose rows may name ONE LINE of a degenerate conic
(port of cellInterior.m).

C is a table of [conic(6), branch] rows; rows is a list of [index (1-based) into C, side].
Lines are decided exactly; curved rows by cellHasInterior's witness search; a cell of half-planes
(rational or branch lines) by a double-precision margin LP whose TIES are settled exactly in
Q(t) / Q(t1,t2) by a Farkas certificate. Returns (tf, info).
"""
import math

import numpy as np

from . import rat_q
from .cell_has_interior import cellHasInterior
from .errors import error


def cellInterior(C, rows):
    if any(len(c) != 7 for c in C):
        error('cellInterior:shape', 'C must be k x 7 -- [conic(1:6), branch].')
    info = {'path': 'rational', 'margin': math.inf}
    if not rows:
        return True, info
    lin = []
    curved = []
    sides = []
    H = []
    src = []
    for t, row in enumerate(rows):
        c = C[row[0] - 1][0:6]
        br = C[row[0] - 1][6]
        sd = row[1]
        if all(x == 0 for x in c[0:5]):
            error('cellInterior:constant',
                  'row %d is a CONSTANT, not a curve: %s. Such a row is vacuous when it holds '
                  'and empties the cell when it does not, and must be resolved before the '
                  'arrangement is built.' % (t + 1, c))
        if br != 0:
            H.append([sd * x for x in _branchLine(c, br)])
            src.append(t)
        elif all(x == 0 for x in c[0:3]):
            lin.append([sd * x for x in c[3:6]])
            H.append([float(sd * x) for x in c[3:6]])
            src.append(t)
        else:
            curved.append(list(c))
            sides.append(sd)
    if lin and not rat_q.feasible2(lin, True):
        return False, info
    if curved:
        tab = curved + [[0, 0, 0] + l for l in lin]
        rr = [[i + 1, s] for i, s in enumerate(sides)] + \
             [[len(curved) + i + 1, 1] for i in range(len(lin))]
        if not cellHasInterior(tab, rr):
            return False, info
    if len(H) < 2:
        return True, info
    m = _marginLP(H)
    info['margin'] = m
    info['path'] = 'margin'
    tol = 1e-9 * max(1.0, max(abs(x) for r in H for x in r))
    if m > tol:
        return True, info
    if m < -tol:
        return False, info
    info['path'] = 'exactTie'
    empty, skipped = _certificate(C, rows, src)
    if empty:
        return False, info
    if skipped:
        info['path'] = 'unsureTie'                 # kept, because keeping only costs a cell
    return True, info


def _certificate(C, rows, src):
    empty = False
    skipped = False
    m = len(src)
    if m < 2:
        return empty, skipped
    for i in range(m - 1):
        for j in range(i + 1, m):
            e, sk = _tryOne(C, rows, [src[i], src[j]])
            empty = empty or e
            skipped = skipped or sk
            if empty:
                return empty, skipped
            for k in range(j + 1, m):
                e, sk = _tryOne(C, rows, [src[i], src[j], src[k]])
                empty = empty or e
                skipped = skipped or sk
                if empty:
                    return empty, skipped
    return empty, skipped


def _uniqueRows(M):
    return sorted(set(tuple(r) for r in M))


def _tryOne(C, rows, idx):
    """idx: 0-based positions into rows."""
    conics = [C[rows[t][0] - 1][0:6] for t in idx if C[rows[t][0] - 1][6] != 0]
    quad = []
    iv = (-1, 1)
    ref = 0
    if conics:
        nu = len(_uniqueRows(conics))
        if nu == 2:
            return _tryTwoFields(C, rows, idx)
        if nu > 2:
            return False, True
        pd = rat_q.pairLineData(conics[0])
        quad = pd['quad']
        for t in idx:
            if C[rows[t][0] - 1][6] != 0:
                ref = C[rows[t][0] - 1][6]
                break
        if len(quad) >= 3:
            iv = pd['iv'][ref - 1]
    A1 = []
    A2 = []
    Cc = []
    for t in idx:
        a, c, ok = _exactHalfPlane(C, rows, t, ref)
        if not ok:
            return False, True
        A1.append(a[0])
        A2.append(a[1])
        Cc.append(c)
    if not quad:
        quad = [1, 0]                               # no branch row: t never occurs
    if len(quad) == 2:
        quad = [1, 0]                               # a degenerate pair: every slope is rational
    if len(A1) == 2:
        cr = _cross2(A1[0], A2[0], A1[1], A2[1], quad)
        if rat_q.qsSign(cr, quad, iv) != 0:
            return False, False
        dt = rat_q.qsAdd(rat_q.qsMul(A1[0], A1[1], quad), rat_q.qsMul(A2[0], A2[1], quad))
        if rat_q.qsSign(dt, quad, iv) >= 0:
            return False, False
        n2 = rat_q.qsAdd(rat_q.qsMul(A1[1], A1[1], quad), rat_q.qsMul(A2[1], A2[1], quad))
        lam = [n2, _qneg(dt)]
    else:
        lam = [_cross2(A1[1], A2[1], A1[2], A2[2], quad),
               _cross2(A1[2], A2[2], A1[0], A2[0], quad),
               _cross2(A1[0], A2[0], A1[1], A2[1], quad)]
        sg = [rat_q.qsSign(l, quad, iv) for l in lam]
        if all(s <= 0 for s in sg):
            lam = [_qneg(l) for l in lam]
            sg = [-s for s in sg]
        if not all(s >= 0 for s in sg):
            return False, False                     # no nonnegative combination annihilates a
        if all(s == 0 for s in sg):
            return False, False                     # the three normals are all parallel
    tot = [0, 0, 1]
    for k in range(len(Cc)):
        tot = rat_q.qsAdd(tot, rat_q.qsMul(lam[k], Cc[k], quad))
    return rat_q.qsSign(tot, quad, iv) <= 0, False


def _exactHalfPlane(C, rows, t, ref):
    a = [[0, 0, 1], [0, 0, 1]]
    c = [0, 0, 1]
    con = C[rows[t][0] - 1][0:6]
    br = C[rows[t][0] - 1][6]
    sd = rows[t][1]
    if br == 0:
        if not all(x == 0 for x in con[0:3]):
            return a, c, False                      # a conic is not a half-plane
        return [[sd * con[3], 0, 1], [sd * con[4], 0, 1]], [sd * con[5], 0, 1], True
    pd = rat_q.pairLineData(con)
    if pd['kind'] == 'crossing' and pd['vertical'] and br == 2:
        p = pd['apex']                              # the VERTICAL line x1 = p1, rational
        return [[sd * p[2], 0, 1], [0, 0, 1]], [-sd * p[0], 0, 1], True
    T = _rootElement(pd['quad'], br, ref)           # the slope or level, in Q(t_ref)
    if T is None:
        return a, c, False
    if pd['kind'] == 'parallel':
        a = [[sd * pd['form'][0], 0, 1], [sd * pd['form'][1], 0, 1]]
        c = rat_q.qsMul([-sd, 0, 1], T, pd['quad'])
    else:
        p = pd['apex']
        a = [rat_q.qsMul([-sd * p[2], 0, 1], T, pd['quad']), [sd * p[2], 0, 1]]
        c = rat_q.qsAdd([-sd * p[1], 0, 1], rat_q.qsMul([sd * p[0], 0, 1], T, pd['quad']))
    return a, c, True


def _tryTwoFields(C, rows, idx):
    cs = [C[rows[t][0] - 1][0:6] for t in idx if C[rows[t][0] - 1][6] != 0]
    u = [list(r) for r in _uniqueRows(cs)]
    pd1 = rat_q.pairLineData(u[0])
    pd2 = rat_q.pairLineData(u[1])
    Q1 = pd1['quad']
    Q2 = pd2['quad']
    if len(Q1) != 3 or len(Q2) != 3:
        return False, True
    ref1 = 0
    ref2 = 0
    for t in idx:
        b = C[rows[t][0] - 1][6]
        if b == 0:
            continue
        if C[rows[t][0] - 1][0:6] == u[0] and ref1 == 0:
            ref1 = b
        if C[rows[t][0] - 1][0:6] == u[1] and ref2 == 0:
            ref2 = b
    iv1 = pd1['iv'][ref1 - 1]
    iv2 = pd2['iv'][ref2 - 1]
    A1 = []
    A2 = []
    Cc = []
    for t in idx:
        b = C[rows[t][0] - 1][6]
        if b == 0:
            a, c, ok = _exactHalfPlane(C, rows, t, 0)
            slot = 0
        elif C[rows[t][0] - 1][0:6] == u[0]:
            a, c, ok = _exactHalfPlane(C, rows, t, ref1)
            slot = 1
        else:
            a, c, ok = _exactHalfPlane(C, rows, t, ref2)
            slot = 2
        if not ok:
            return False, True
        A1.append(_embed(a[0], slot))
        A2.append(_embed(a[1], slot))
        Cc.append(_embed(c, slot))
    n = len(A1)
    if n == 2:
        cr = _cross2two(A1[0], A2[0], A1[1], A2[1], Q1, Q2)
        sg, ok = rat_q.qs2Sign(cr, Q1, iv1, Q2, iv2)
        if not ok:
            return False, True
        if sg != 0:
            return False, False
        dt = rat_q.qs2Add(rat_q.qs2Mul(A1[0], A1[1], Q1, Q2), rat_q.qs2Mul(A2[0], A2[1], Q1, Q2))
        sg, ok = rat_q.qs2Sign(dt, Q1, iv1, Q2, iv2)
        if not ok:
            return False, True
        if sg >= 0:
            return False, False
        n2 = rat_q.qs2Add(rat_q.qs2Mul(A1[1], A1[1], Q1, Q2), rat_q.qs2Mul(A2[1], A2[1], Q1, Q2))
        lam = [n2, _qneg2(dt)]
    else:
        lam = [_cross2two(A1[1], A2[1], A1[2], A2[2], Q1, Q2),
               _cross2two(A1[2], A2[2], A1[0], A2[0], Q1, Q2),
               _cross2two(A1[0], A2[0], A1[1], A2[1], Q1, Q2)]
        sg = [0, 0, 0]
        for k in range(3):
            sg[k], ok = rat_q.qs2Sign(lam[k], Q1, iv1, Q2, iv2)
            if not ok:
                return False, True
        if all(s <= 0 for s in sg):
            lam = [_qneg2(l) for l in lam]
            sg = [-s for s in sg]
        if not all(s >= 0 for s in sg):
            return False, False
        if all(s == 0 for s in sg):
            return False, False
    tot = [0, 0, 0, 0, 1]
    for k in range(len(Cc)):
        tot = rat_q.qs2Add(tot, rat_q.qs2Mul(lam[k], Cc[k], Q1, Q2))
    sg, ok = rat_q.qs2Sign(tot, Q1, iv1, Q2, iv2)
    if not ok:
        return False, True
    return sg <= 0, False


def _embed(u, slot):
    if slot == 1:
        w = [u[0], u[1], 0, 0, u[2]]
    elif slot == 2:
        w = [u[0], 0, u[1], 0, u[2]]
    else:
        if u[1] != 0:
            error('cellInterior:embed', 'a rational row cannot carry a t term.')
        w = [u[0], 0, 0, 0, u[2]]
    return rat_q.qs2Red(w)


def _cross2two(u1, u2, v1, v2, Q1, Q2):
    return rat_q.qs2Add(rat_q.qs2Mul(u1, v2, Q1, Q2), _qneg2(rat_q.qs2Mul(u2, v1, Q1, Q2)))


def _qneg2(u):
    return [-u[0], -u[1], -u[2], -u[3], u[4]]


def _rootElement(quad, br, ref):
    if len(quad) == 2:
        return [-quad[1], 0, quad[0]]               # degenerate: the slope is rational
    if len(quad) != 3:
        return None
    if br == ref:
        T = [0, 1, 1]
    else:
        T = [-quad[1], -quad[0], quad[0]]
    return rat_q.qsRed(T)


def _cross2(u1, u2, v1, v2, Q):
    return rat_q.qsAdd(rat_q.qsMul(u1, v2, Q), _qneg(rat_q.qsMul(u2, v1, Q)))


def _qneg(u):
    return [-u[0], -u[1], u[2]]


def _branchLine(c, b):
    v = rat_q.pairLineEval(c, b, [[0, 0], [1, 0], [0, 1]])
    L = [v[1] - v[0], v[2] - v[0], v[0]]
    n = math.hypot(L[0], L[1])
    if n > 0:
        L = [x / n for x in L]
    return L


def _marginLP(H):
    """Largest t with H*[x;1] >= t inside a big box, by vertex enumeration (doubles)."""
    R = 1e4
    Hn = []
    for r in H:
        n = math.hypot(r[0], r[1])
        Hn.append([x / n for x in r] if n > 0 else [float(x) for x in r])
    M = np.array([[r[0], r[1], -1.0] for r in Hn]
                 + [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]], dtype=float)
    b = np.array([-r[2] for r in Hn] + [-R, -R, -R, -R], dtype=float)
    m = -math.inf
    k = M.shape[0]
    for i in range(k - 2):
        for j in range(i + 1, k - 1):
            for l in range(j + 1, k):
                S = M[[i, j, l], :]
                if abs(np.linalg.det(S)) < 1e-12:
                    continue
                z = np.linalg.solve(S, b[[i, j, l]])
                if np.all(M @ z >= b - 1e-9) and z[2] > m:
                    m = float(z[2])
    if math.isinf(m):
        m = -1.0                                    # no vertex satisfies all: empty
    return m

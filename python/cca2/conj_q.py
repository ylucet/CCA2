"""conj_q -- the Legendre-Fenchel conjugate computed EXACTLY over the rationals (port of conjQ.m).

    f*(s) = sup_x <s,x> - f(x)

with every coefficient an exact rational and every decision an integer predicate. Input: an exact
QuaPol. Output: (g, dom) with g a QuaCon (exact faces, primitive integer edge conics, NAMED
vertices) and dom a list of integer half-planes [a b c] (a s1 + b s2 + c >= 0) whose intersection
is dom f* -- empty when f has no unbounded piece, and the impossible row [0 0 -1] when dom f* is
empty.

EXACT OR A NAMED REFUSAL, NEVER A FALLBACK. conjQ.m's header has the measured reasons.

The case analysis per piece (conjPieceQ) is two questions asked in order: can the sup be finite
at all (a condition on the recession cone), and is H positive definite (then the maximiser can be
interior: the KKT active set, caseB) or not (then it lies on the boundary: caseD). A multi-face
input is the fold of its pieces through maxQ, because conjugation turns a union into a max.
"""
from fractions import Fraction
from math import gcd

from . import rat_q
from .assemble_qua_con_cells import assembleQuaConCells
from .errors import error
from .max_q import maxQ
from .qua_con import QuaCon
from .qua_pol import QuaPol
from .rat_q import sign, dot, mv, quad


def assertNoClosedEdges(E, who):
    bad = [j + 1 for j, e in enumerate(E) if e[2] == 2]
    if bad:
        error('RatCon:closedEdge',
              '%s cannot yet handle a CLOSED edge (E(:,3) = 2): its edge handling assumes every '
              'edge names two endpoints. Edge(s): %s.' % (who, bad))


def conjQ(obj):
    """Returns (g, dom)."""
    if not isinstance(obj, QuaPol):
        error('PLQ:conjQ:input',
              'conjQ takes a QuaPol (quadratic on polyhedral); got a %s.' % type(obj).__name__)
    assertNoClosedEdges(obj.E, 'conjQ')
    obj.assertOperable()
    obj.assertExact()

    nf = len(obj.fN)

    # ---- a CHAIN carrying one function per edge -------------------------------------------------
    if obj.dom['dim'] < 2 and nf > 1:
        if nf != obj.ne:
            error('PLQ:conjQ:ambiguousChain',
                  'a domain of dimension < 2 carrying %d function rows on %d edges: the '
                  'convention is one row per edge, so which edge takes which is undefined.'
                  % (nf, obj.ne))
        g = None
        dom = []
        for j in range(1, obj.ne + 1):
            e = obj.E[j - 1]
            ej = QuaPol([obj.Vx(e[0]), obj.Vx(e[1])], [[1, 2, e[2]]], [obj.fx(j)], [[0, 0]])
            gj, dj = conjPieceQ(ej, 1)
            g = gj if g is None else maxQ(g, gj)
            dom = dom + dj
        return g, dom

    g, dom = conjPieceQ(obj, 1)
    for k in range(2, nf + 1):
        gk, dk = conjPieceQ(obj, k)
        g = maxQ(g, gk)
        dom = dom + dk
    return g, dom


# ==================================================================================================
def conjPieceQ(obj, k):
    """The exact conjugate of ONE piece -- q_k restricted to its own face P_k. Returns (g, dom)."""
    fN, fD = obj.faceQ(k)
    Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
    D = rat_q.detExact(Hn)
    dom = []

    if obj.nv == 0:
        return caseAFullDomain(fN, fD), dom       # the whole plane: no domain constraints at all

    Vi, vd = rat_q.combineDen(obj.VN, obj.VD)

    # ---- a domain of dimension < 2: a needle, a segment, a ray, or a chain of edges ------------
    if obj.dom['dim'] < 2:
        if len(obj.fN) > 1:
            error('PLQ:conjQ:ambiguousChain',
                  'a domain of dimension < 2 carrying %d different functions has no F to say '
                  'which edge takes which, so the mapping is ambiguous. One function only.'
                  % len(obj.fN))
        sh = degenerateShape(obj, Vi)
        okFinite, dom, _ = recessionConditions(sh, Hn, [fN[7], fN[8]], fD, Vi, vd)
        if not okFinite:
            return emptyConjugate(), [[0, 0, -1]]
        return caseDBoundaryMax(fN, fD, Vi, vd, sh, dom), dom

    sh = pieceShape(obj, k, Vi, vd)
    if len(sh['vs']) == 0:
        error('PLQ:conjQ:noVertex',
              'piece %d has no vertex, so it is a half-plane or a slab and its conjugate is '
              'supported on a line or a point -- a domain of dimension < 2, which this type '
              'cannot carry.' % k)

    # ---- a NON-CONVEX piece is SPLIT, not refused ----------------------------------------------
    if not sh['convex']:
        if not sh['bounded']:
            return conjUnboundedNonConvex(obj, k)
        tris = earClip(Vi, boundaryCycle(sh, Vi))
        g = None
        for tri in tris:
            sht = triangleShape(Vi, vd, tri)
            if Hn[0][0] > 0 and D > 0:
                gt = caseBConvexOnPiece(fN, fD, Vi, vd, sht, [])
            else:
                gt = caseDBoundaryMax(fN, fD, Vi, vd, sht, [])
            g = gt if g is None else maxQ(g, gt)
        return g, dom

    okFinite, dom, why = recessionConditions(sh, Hn, [fN[7], fN[8]], fD, Vi, vd)
    if not okFinite:
        # dom f* is EMPTY: a QuaCon with zero faces evaluates to +infinity everywhere.
        assert why
        return emptyConjugate(), [[0, 0, -1]]     # 0 >= 1: no s at all

    gT, isThin = linealityConj(obj, k, sh, Hn, fN, fD, Vi, vd)
    if isThin:
        return gT, dom                            # dom f* is a LINE (or part of one)

    if Hn[0][0] > 0 and D > 0:
        g = caseBConvexOnPiece(fN, fD, Vi, vd, sh, dom)     # strictly convex: KKT active set
    else:
        g = caseDBoundaryMax(fN, fD, Vi, vd, sh, dom)       # everything else: max over the boundary
    return g, dom


def emptyConjugate():
    """The function that is +infinity everywhere: a QuaCon with NO faces."""
    return QuaCon([], [], [], [], [], [], [])


def degenerateShape(obj, Vi):
    """The shape record for a domain of dimension < 2: every edge belongs to the one piece and
    there are no bounding half-planes (the segment IS the candidate set)."""
    sh = {'vs': [], 'ed': [], 'rays': [], 'hp': [], 'bounded': True, 'convex': True}
    for e in obj.E:
        a, b = e[0], e[1]
        d = [Vi[b - 1][0] - Vi[a - 1][0], Vi[b - 1][1] - Vi[a - 1][1]]
        if e[2] != 0:
            sh['vs'] += [a, b]
            sh['ed'].append({'a': a, 'b': b, 'd': d, 'isRay': False})
        else:
            sh['vs'] += [a]
            sh['ed'].append({'a': a, 'b': 0, 'd': d, 'isRay': True})
            sh['rays'].append(d)
            sh['bounded'] = False
    if not sh['ed']:
        sh['vs'] = list(range(1, len(obj.VN) + 1))   # a NEEDLE: the single vertex, no edge at all
    sh['vs'] = sorted(set(sh['vs']))
    return sh


def caseAFullDomain(fN, fD):
    """The conjugate of f(x) = 1/2 x'Hx + L'x + kappa over ALL of R^2, exactly."""
    Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
    Ln = [fN[7], fN[8]]
    D = rat_q.detExact(Hn)

    # Sylvester's criterion on exact integers: H > 0 iff Hn(1,1) > 0 and det(Hn) > 0.
    if not (Hn[0][0] > 0 and D > 0):
        # A direction of NEGATIVE curvature: f* is +infinity everywhere.
        if not rat_q.isPSD2(Hn):
            return emptyConjugate()
        # Otherwise H is PSD and singular: dom f* is THIN (a line, or a point when H = 0),
        # carried by the H-form's equality side.
        if all(x == 0 for r in Hn for x in r):
            # q is AFFINE on the plane: f*(s) = -kappa at the single point s = L.
            con = [rat_q.conic([0, 0, 0, fD, 0, -Ln[0]]) + [sgnOf([fD, 0, -Ln[0]], 0)],
                   rat_q.conic([0, 0, 0, 0, fD, -Ln[1]]) + [sgnOf([0, fD, -Ln[1]], 0)]]
            return assembleQuaConCells([{'num': [0, 0, 0, 0, 0, 0, 0, 0, 0, -fN[9]],
                                         'den': fD, 'con': con}])
        # q is positive semidefinite and SINGULAR: finite exactly on the line u.n = 0, u = fD*s-Ln.
        if Hn[0][0] != 0:
            nl = [-Hn[0][1], Hn[0][0]]              # Hn*nl = 0 because det(Hn) = 0
        else:
            nl = [1, 0]                             # Hn(1,1) = 0 forces Hn(1,2) = 0 here
        m = [Hn[0][0], Hn[0][1]]
        if all(x == 0 for x in m):
            m = [Hn[0][1], Hn[1][1]]
        mHm = quad(m, Hn, m)
        P = [fD * m[0], fD * m[1], -dot(Ln, m)]
        den = 2 * mHm * fD
        num = [0, 0, 0, 0, 2 * P[0] ** 2, 2 * P[0] * P[1], 2 * P[1] ** 2,
               2 * P[0] * P[2], 2 * P[1] * P[2], P[2] ** 2 - 2 * mHm * fN[9]]
        row = [fD * nl[0], fD * nl[1], -dot(Ln, nl)]
        return assembleQuaConCells([{'num': num, 'den': den,
                                     'con': [rat_q.conic([0, 0, 0] + row) + [sgnOf(row, 0)]]}])

    A = [[Hn[1][1], -Hn[0][1]], [-Hn[0][1], Hn[0][0]]]      # adjugate of a 2x2, exactly
    ALn = mv(A, Ln)
    den = 2 * fD * D
    c5 = 2 * fD ** 2 * A[0][0]
    c6 = 2 * fD ** 2 * A[0][1]
    c7 = 2 * fD ** 2 * A[1][1]
    c8 = -2 * fD * ALn[0]
    c9 = -2 * fD * ALn[1]
    c10 = dot(Ln, ALn) - 2 * D * fN[9]
    gN, gD = rat_q.canon([0, 0, 0, 0, c5, c6, c7, c8, c9, c10], den)
    return QuaCon([], [], [], [gN], [gD], [], [[]])


def caseBConvexOnPiece(fN, fD, Vi, vd, sh, dom):
    """The conjugate of a STRICTLY CONVEX quadratic over one polyhedral piece, bounded or not:
    the KKT active-set subdivision -- one AFFINE cell per vertex (its normal cone shifted by
    grad q(v)), one QUADRATIC cell per edge, one interior cell."""
    Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
    Ln = [fN[7], fN[8]]
    kn = fN[9]
    D = rat_q.detExact(Hn)
    if not (Hn[0][0] > 0 and D > 0):
        error('PLQ:conjQ:notStrictlyConvex',
              'caseB needs a strictly convex Q; leading minor %d, determinant %d.'
              % (Hn[0][0], D))

    denV = 2 * fD * vd ** 2
    cells = []

    # ---- one cell per VERTEX: its normal cone, shifted by grad q(v) ---------------------------
    for vi in sh['vs']:
        v = Vi[vi - 1]
        gv = [a + vd * b for a, b in zip(mv(Hn, v), Ln)]            # over fD*vd
        qv = quad(v, Hn, v) + 2 * vd * dot(Ln, v) + 2 * vd ** 2 * kn
        num = [0, 0, 0, 0, 0, 0, 0, 2 * fD * vd * v[0], 2 * fD * vd * v[1], -qv]
        rows = [list(r) for r in dom]
        for ed in dirsAtVertex(sh, vi):
            # <s - grad q(v), e> <= 0, cleared by fD*vd, written as a >= 0 half-plane
            rows.append([-fD * vd * ed[0], -fD * vd * ed[1], dot(gv, ed)])
        con, ok = toCon(rows)
        if ok:
            cells.append({'num': num, 'den': denV, 'con': con})

    # ---- one cell per EDGE of positive curvature ----------------------------------------------
    for j, e in enumerate(sh['ed']):
        a = Vi[e['a'] - 1]
        d = e['d']
        alpha = quad(d, Hn, d)
        if alpha <= 0:
            error('PLQ:conjQ:degenerateEdge',
                  'edge %d has d\'Qd = %d <= 0 under a strictly convex Q.' % (j + 1, alpha))
        ga = [x + vd * y for x, y in zip(mv(Hn, a), Ln)]
        qa = quad(a, Hn, a) + 2 * vd * dot(Ln, a) + 2 * vd ** 2 * kn
        Td = [fD * vd * d[0], fD * vd * d[1], -dot(ga, d)]

        den = 2 * alpha * fD * vd ** 2
        lin = [2 * alpha * fD * vd * a[0], 2 * alpha * fD * vd * a[1]]
        num = [0, 0, 0, 0, 2 * Td[0] ** 2, 2 * Td[0] * Td[1], 2 * Td[1] ** 2,
               lin[0] + 2 * Td[0] * Td[2], lin[1] + 2 * Td[1] * Td[2], Td[2] ** 2 - alpha * qa]

        rows = [list(r) for r in dom] + [Td]                         # 0 <= t*
        if not e['isRay']:
            rows.append([-Td[0], -Td[1], alpha - Td[2]])
        # the multiplier condition, cleared by alpha > 0
        n = outwardNormalOf(sh, j)
        Hdn = quad(d, Hn, n)
        rows.append([alpha * fD * vd * n[0] - Td[0] * Hdn,
                     alpha * fD * vd * n[1] - Td[1] * Hdn,
                     -alpha * dot(ga, n) - Td[2] * Hdn])
        con, ok = toCon(rows)
        if ok:
            cells.append({'num': num, 'den': den, 'con': con})

    # ---- the INTERIOR cell ----------------------------------------------------------------------
    A = [[Hn[1][1], -Hn[0][1]], [-Hn[0][1], Hn[0][0]]]
    ALn = mv(A, Ln)
    numI = [0, 0, 0, 0, 2 * fD ** 2 * A[0][0], 2 * fD ** 2 * A[0][1], 2 * fD ** 2 * A[1][1],
            -2 * fD * ALn[0], -2 * fD * ALn[1], dot(Ln, ALn) - 2 * D * kn]
    rows = [list(r) for r in dom]
    for hp in sh['hp']:
        # hp is [p q c] meaning p*x + q*y + c >= 0 at actual points; substitute x = H^-1(s - L)
        p = hp[0:2]
        Ap = [A[0][0] * p[0] + A[1][0] * p[1], A[0][1] * p[0] + A[1][1] * p[1]]   # A' * p
        rows.append([fD * Ap[0], fD * Ap[1], -dot(Ap, Ln) + D * hp[2]])
    con, ok = toCon(rows)
    if ok:
        cells.append({'num': numI, 'den': 2 * fD * D, 'con': con})

    return assembleQuaConCells(cells)


def outwardNormalOf(sh, j):
    """The OUTWARD normal of edge j (0-based position in sh['ed']), read from the inward normal
    pieceShape already decided -- never re-derived (a half-plane piece cannot orient it)."""
    e = sh['ed'][j]
    if e.get('nIn'):
        return [-e['nIn'][0], -e['nIn'][1]]
    d = e['d']
    return [d[1], -d[0]]


def toCon(rows):
    """[a b c] half-planes -> H-form rows [conic6, side]; a vacuous row is dropped, an
    unsatisfiable one makes the cell EMPTY (ok False)."""
    con = []
    for r in rows:
        if r[0] == 0 and r[1] == 0:
            if r[2] < 0:
                return con, False
            continue                                   # vacuous
        con.append(rat_q.conic([0, 0, 0] + list(r)) + [sgnOf(r, +1)])
    return con, True


def sgnOf(row, want):
    """The H-form sign for `row >= 0` (want = +1) or `row <= 0` (want = -1) AFTER rat_q.conic
    has normalised the row's overall sign; an equality (want = 0) is symmetric."""
    if want == 0:
        return 0
    r = [0, 0, 0] + list(row) if len(row) == 3 else list(row)
    nz = next((x for x in r if x != 0), None)
    if nz is None:
        error('ratQ:zeroConic', 'a constraint with all-zero coefficients names no half-plane.')
    return -want if nz < 0 else want


def pieceShape(obj, k, Vi, vd):
    """Everything the conjugate needs to know about ONE piece's geometry, bounded or not:
      vs      the piece's true vertex indices (1-based; a ray's direction marker is NOT one)
      ed      edges: a (base vertex), b (other vertex, 0 for a ray), d (integer direction),
              isRay, j (edge index in obj.E), nIn (inward normal)
      rays    the integer extreme recession directions
      hp      half-planes [a b c], a*x + b*y + c >= 0 on the piece, in ACTUAL coordinates
      bounded, convex
    The SIDE of each edge is decided at numerators (where the vertices live); the row is emitted
    for actual points. Mixing the two rescales the offset into a PARALLEL line."""
    own = [j + 1 for j, r in enumerate(obj.F) if r[0] == k or r[1] == k]
    if not own:
        error('PLQ:conjQ:noFace', 'face %d has no edges.' % k)

    sh = {'vs': [], 'ed': [], 'rays': [], 'hp': [], 'bounded': True, 'convex': True}
    for jg in own:
        e = obj.E[jg - 1]
        a, b = e[0], e[1]
        d = [Vi[b - 1][0] - Vi[a - 1][0], Vi[b - 1][1] - Vi[a - 1][1]]
        if e[2] != 0:
            sh['vs'] += [a, b]
            sh['ed'].append({'a': a, 'b': b, 'd': d, 'isRay': False, 'j': jg, 'nIn': None})
        else:
            sh['vs'] += [a]
            sh['ed'].append({'a': a, 'b': 0, 'd': d, 'isRay': True, 'j': jg, 'nIn': None})
            sh['rays'].append(d)
    sh['vs'] = sorted(set(sh['vs']))
    sh['bounded'] = not sh['rays']

    # ---- the half-planes -----------------------------------------------------------------------
    for ed in sh['ed']:
        p = Vi[ed['a'] - 1]
        n = [ed['d'][1], -ed['d'][0]]                 # a normal to this edge
        c0 = -dot(n, p)
        sgn = 0
        for vi in sh['vs']:
            v = Vi[vi - 1]
            t = n[0] * v[0] + n[1] * v[1] + c0
            if t != 0:
                sgn = sign(t)
                break
        if sgn == 0:
            for r in sh['rays']:
                t = dot(n, r)
                if t != 0:
                    sgn = sign(t)
                    break
        if sgn == 0:
            # EVERY vertex and recession direction lies ON this edge's line: a piece bounded BY
            # A LINE, i.e. a half-plane. Only F knows which side; n = (d2, -d1) is the RIGHT normal.
            jg = ed['j']
            if obj.F[jg - 1][1] == k:
                sgn = 1
            elif obj.F[jg - 1][0] == k:
                sgn = -1
            else:
                error('PLQ:conjQ:badSubdivision',
                      'edge %d lists faces [%d %d], neither of which is face %d.'
                      % (jg, obj.F[jg - 1][0], obj.F[jg - 1][1], k))

        # IS THE PIECE CONVEX? Recorded rather than raised: a non-convex face is SPLIT.
        for vi in sh['vs']:
            v = Vi[vi - 1]
            if sgn * (n[0] * v[0] + n[1] * v[1] + c0) < 0:
                sh['convex'] = False
        for r in sh['rays']:
            if sgn * dot(n, r) < 0:
                sh['convex'] = False

        sh['hp'].append([sgn * vd * n[0], sgn * vd * n[1], sgn * c0])
        ed['nIn'] = [sgn * n[0], sgn * n[1]]
    return sh


def dirsAtVertex(sh, i):
    """The edge directions leaving vertex index i. A ray contributes at its BASE only."""
    dirs = []
    for e in sh['ed']:
        if e['a'] == i:
            dirs.append(list(e['d']))
        elif not e['isRay'] and e['b'] == i:
            dirs.append([-e['d'][0], -e['d'][1]])
    return dirs


def recessionConditions(sh, Hn, Ln, fD, Vi, vd):
    """Can the sup be finite on this piece, and under which LINEAR conditions on s?
    Returns (ok, dom, why). Along a recession direction d the objective carries
    t*slope(x) - (t^2/2) d'Hd:  d'Hd > 0 no condition;  d'Hd < 0 dom f* EMPTY;  d'Hd = 0 finite
    exactly when the slope is <= 0 for EVERY x of the piece -- an LP over the polyhedron."""
    ok = True
    dom = []
    why = ''
    R = sh['rays']
    if not R:
        return ok, dom, why

    # ---- IS THE RECESSION CONE A HALF-PLANE? --------------------------------------------------
    if len(sh['hp']) >= 1:
        Np = sorted(set(tuple(r) for r in sign2([hp[0:2] for hp in sh['hp']])))
        if len(Np) == 1:                           # one distinct inward normal: a half-plane
            p = list(Np[0])
            if not rat_q.isPSD2(Hn):
                return False, dom, ('the piece recedes over a HALF-PLANE of directions, on which '
                                    'some direction has negative curvature (H is not positive '
                                    'semidefinite)')
            if rat_q.detExact(Hn) == 0:
                if any(x != 0 for r in Hn for x in r):
                    d0 = [-Hn[0][1], Hn[0][0]] if Hn[0][0] != 0 else [1, 0]
                    cand = [d0, [-d0[0], -d0[1]]]
                else:
                    cand = [[1, 0], [-1, 0], [0, 1], [0, -1]]   # H = 0: every direction is null
                for d in cand:
                    if dot(p, d) < 0:
                        continue
                    row = [-fD * vd * d[0], -fD * vd * d[1], vd * dot(Ln, d)]
                    if all(x == 0 for x in row):
                        continue
                    dom.append(row)
            return ok, dom, why

    for r in R:
        if quad(r, Hn, r) < 0:
            return False, dom, ('the piece recedes along (%d,%d), where the objective rises '
                                'without bound for every s' % (r[0], r[1]))
    # THE CONE, NOT ONLY ITS GENERATORS: a*l1^2 + 2b*l1*l2 + c*l2^2 >= 0 for every l >= 0.
    if len(R) >= 2:
        a = quad(R[0], Hn, R[0])
        b = quad(R[0], Hn, R[1])
        c = quad(R[1], Hn, R[1])
        if b < 0 and b * b - a * c > 0:
            return False, dom, ('some direction strictly inside the recession cone has negative '
                                'curvature')

    for d in R:
        if quad(d, Hn, d) != 0:
            continue
        Hd = mv(Hn, d)
        # inf over the piece of <Hd, x>: -infinity if it decreases along some recession direction
        for q in R:
            if dot(q, Hd) < 0:
                return False, dom, ('the objective\'s slope along the null recession direction '
                                    '(%d,%d) is unbounded above on the piece' % (d[0], d[1]))
        mu = None
        for vi in sh['vs']:
            t = dot(Vi[vi - 1], Hd)
            if mu is None or t < mu:
                mu = t                             # over vd, like Vi itself
        # <s,d> - (Ln.d)/fD <= mu/vd, cleared by fD*vd > 0 and written as a >= 0 half-plane
        row = [-fD * vd * d[0], -fD * vd * d[1], fD * mu + vd * dot(Ln, d)]
        if all(x == 0 for x in row):
            continue
        dom.append(row)
    return ok, dom, why


def caseDBoundaryMax(fN, fD, Vi, vd, sh, dom):
    """The conjugate of a quadratic that is NOT positive definite over one polyhedral piece.

    The sup is attained on the BOUNDARY, so q* is the largest of finitely many candidates: every
    vertex value <s,v> - q(v) (affine), and for every edge of positive curvature its clamped 1-D
    maximum (quadratic, available where the clamp is inactive -- an AFFINE condition). Classify
    first (a LINE arrangement, decided by feasible2), split last (one conic per pair of surviving
    candidates)."""
    Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
    Ln = [fN[7], fN[8]]
    kn = fN[9]
    m = len(sh['vs'])

    denV = 2 * fD * vd ** 2
    qv = []
    numV = []
    for vi in sh['vs']:
        v = Vi[vi - 1]
        q = quad(v, Hn, v) + 2 * vd * dot(Ln, v) + 2 * vd ** 2 * kn
        qv.append(q)
        numV.append([0, 0, 0, 0, 0, 0, 0, 2 * fD * vd * v[0], 2 * fD * vd * v[1], -q])

    Ed = []
    for e in sh['ed']:
        d = e['d']
        alpha = quad(d, Hn, d)
        if alpha <= 0:
            continue
        a = Vi[e['a'] - 1]
        ga = [x + vd * y for x, y in zip(mv(Hn, a), Ln)]
        qa = quad(a, Hn, a) + 2 * vd * dot(Ln, a) + 2 * vd ** 2 * kn
        Td = [fD * vd * d[0], fD * vd * d[1], -dot(ga, d)]
        lin = [2 * alpha * fD * vd * a[0], 2 * alpha * fD * vd * a[1]]
        Ed.append({'Td': Td, 'alpha': alpha,
                   'num': [0, 0, 0, 0, 2 * Td[0] ** 2, 2 * Td[0] * Td[1], 2 * Td[1] ** 2,
                           lin[0] + 2 * Td[0] * Td[2], lin[1] + 2 * Td[1] * Td[2],
                           Td[2] ** 2 - alpha * qa],
                   'den': 2 * alpha * fD * vd ** 2, 'isRay': e['isRay']})

    # ---- the LINE arrangement: which vertex leads, and which edges are available ---------------
    cells0 = []
    for i in range(m):
        rows = [list(r) for r in dom]
        vI = Vi[sh['vs'][i] - 1]
        for j in range(m):
            if j == i:
                continue
            vJ = Vi[sh['vs'][j] - 1]
            r = [2 * fD * vd * (vI[0] - vJ[0]), 2 * fD * vd * (vI[1] - vJ[1]), -(qv[i] - qv[j])]
            if all(x == 0 for x in r):
                error('PLQ:conjQ:repeatedVertex',
                      'piece vertices %d and %d coincide.' % (sh['vs'][i], sh['vs'][j]))
            rows.append(r)
        if rat_q.feasible2(rows, True):
            cells0.append({'lin': rows, 'vi': i, 'ei': []})

    for e, ed in enumerate(Ed):
        Td = ed['Td']
        alpha = ed['alpha']
        neg = [-Td[0], -Td[1], -Td[2]]
        if ed['isRay']:
            # clamped at the base only: available on t* >= 0, unavailable on t* <= 0
            regimes = [[neg], [Td]]
            avail = [False, True]
        else:
            regimes = [[neg],
                       [Td, [-Td[0], -Td[1], alpha - Td[2]]],
                       [[Td[0], Td[1], Td[2] - alpha]]]
            avail = [False, True, False]
        nxt = []
        for c in cells0:
            for r in range(len(regimes)):
                lin = c['lin'] + regimes[r]
                if not rat_q.feasible2(lin, True):
                    continue
                ei = c['ei'] + [e] if avail[r] else c['ei']
                nxt.append({'lin': lin, 'vi': c['vi'], 'ei': ei})
        cells0 = nxt

    # ---- inside each linear cell, the winner among its few candidates --------------------------
    out = []
    for c in cells0:
        base, ok = toCon(c['lin'])
        if not ok:
            continue
        cand = [{'num': numV[c['vi']], 'den': denV}]
        for e in c['ei']:
            cand.append({'num': Ed[e]['num'], 'den': Ed[e]['den']})
        for i in range(len(cand)):
            con = [list(r) for r in base]
            for j in range(len(cand)):
                if j == i:
                    continue
                dn, _ = rat_q.sub(cand[i]['num'], cand[i]['den'], cand[j]['num'], cand[j]['den'])
                if all(x == 0 for x in dn):
                    continue
                Rr = [dn[4], 2 * dn[5], dn[6], 2 * dn[7], 2 * dn[8], 2 * dn[9]]
                con.append(rat_q.conic(Rr) + [sgnOf(Rr, +1)])
            out.append({'num': cand[i]['num'], 'den': cand[i]['den'], 'con': con})

    return assembleQuaConCells(out)


def boundaryCycle(sh, Vi):
    """The face's boundary as an ORDERED cycle of vertex indices, counter-clockwise, walked from
    the edge list (so no orientation convention is consulted)."""
    m = len(sh['ed'])
    adj = {}
    for e in sh['ed']:
        a, b = e['a'], e['b']
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, []).append(a)
    cyc = [sh['ed'][0]['a'], sh['ed'][0]['b']]
    for t in range(2, m):
        nxt = [v for v in adj[cyc[t - 1]] if v != cyc[t - 2]]
        if not nxt:
            error('PLQ:conjQ:notSimple', 'the face boundary is not a single closed cycle.')
        cyc.append(nxt[0])
    if len(set(cyc)) != m:
        error('PLQ:conjQ:notSimple',
              'the face boundary revisits a vertex, so it is not a simple polygon.')
    if twiceArea(Vi, cyc) < 0:
        cyc = cyc[::-1]
    return cyc


def twiceArea(Vi, cyc):
    """Twice the SIGNED area of the polygon, an exact integer (only its sign is used)."""
    A2 = 0
    m = len(cyc)
    for i in range(m):
        p = Vi[cyc[i] - 1]
        q = Vi[cyc[(i + 1) % m] - 1]
        A2 += p[0] * q[1] - p[1] * q[0]
    return A2


def earClip(Vi, cyc):
    """Triangulate a simple counter-clockwise polygon EXACTLY by ear clipping: every decision is
    the sign of an integer cross product. Returns a list of [a, b, c] vertex-index triples."""
    cyc = dropCollinear(Vi, list(cyc))
    m = len(cyc)
    if m < 3:
        error('PLQ:conjQ:degenerate', 'the face collapses to fewer than three corners.')
    tris = []
    guard = 0
    while len(cyc) > 3:
        guard += 1
        if guard > m + 2:
            error('PLQ:conjQ:noEar',
                  'no ear was found on a polygon of %d corners, which cannot happen for a simple '
                  'polygon -- so the boundary is self-intersecting.' % len(cyc))
        n = len(cyc)
        clipped = False
        for i in range(n):
            a = cyc[(i - 1) % n]
            b = cyc[i]
            c = cyc[(i + 1) % n]
            if not isEar(Vi, cyc, a, b, c):
                continue
            tris.append([a, b, c])
            del cyc[i]
            clipped = True
            break
        if not clipped:
            error('PLQ:conjQ:noEar',
                  'no ear was found, so the boundary is self-intersecting or not '
                  'counter-clockwise.')
    tris.append(list(cyc))
    return tris


def isEar(Vi, cyc, a, b, c):
    if cross2(Vi, a, b, c) <= 0:
        return False                              # reflex, or collinear: not an ear
    for v in cyc:
        if v == a or v == b or v == c:
            continue
        if inTriangle(Vi, a, b, c, v):
            return False                          # another corner is inside it
    return True


def cross2(Vi, a, b, c):
    """z of (b-a) x (c-b): positive means a LEFT turn at b."""
    p, q, r = Vi[a - 1], Vi[b - 1], Vi[c - 1]
    return (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0])


def inTriangle(Vi, a, b, c, v):
    """Is vertex v inside the counter-clockwise triangle (a,b,c), boundary INCLUDED?"""
    return (cross2v(Vi, a, b, v) >= 0 and cross2v(Vi, b, c, v) >= 0
            and cross2v(Vi, c, a, v) >= 0)


def cross2v(Vi, a, b, v):
    """z of (b-a) x (v-a): positive means v is to the LEFT of the directed edge a->b."""
    p, q, r = Vi[a - 1], Vi[b - 1], Vi[v - 1]
    return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])


def dropCollinear(Vi, cyc):
    changed = True
    while changed and len(cyc) > 3:
        changed = False
        n = len(cyc)
        for i in range(n):
            a = cyc[(i - 1) % n]
            b = cyc[i]
            c = cyc[(i + 1) % n]
            if cross2(Vi, a, b, c) == 0:
                del cyc[i]
                changed = True
                break
    return cyc


def triangleShape(Vi, vd, tri):
    """The shape record for one triangle of the split, in the form pieceShape produces."""
    sh = {'vs': sorted(tri), 'ed': [], 'rays': [], 'hp': [], 'bounded': True, 'convex': True}
    for i in range(3):
        a = tri[i]
        b = tri[(i + 1) % 3]
        o = tri[(i + 2) % 3]
        d = [Vi[b - 1][0] - Vi[a - 1][0], Vi[b - 1][1] - Vi[a - 1][1]]
        sh['ed'].append({'a': a, 'b': b, 'd': d, 'isRay': False})
        n = [d[1], -d[0]]
        c0 = -dot(n, Vi[a - 1])
        t = n[0] * Vi[o - 1][0] + n[1] * Vi[o - 1][1] + c0
        if t == 0:
            error('PLQ:conjQ:degenerate', 'a split triangle has zero area.')
        s = sign(t)
        sh['hp'].append([s * vd * n[0], s * vd * n[1], s * c0])
    return sh


def sign2(N):
    """Each row reduced to a primitive integer direction, SIGN PRESERVED (a normal has a side)."""
    P = []
    for r in N:
        g = gcd(abs(r[0]), abs(r[1]))
        if g == 0:
            P.append(list(r))
        else:
            P.append([r[0] // g, r[1] // g])
    return P


# ------------------------------------------------------------------------------------------------
def conjUnboundedNonConvex(obj, k):
    """The conjugate of q_k on an UNBOUNDED NON-CONVEX face, by a cover with convex parts: the
    face is a boundary graph swept along one direction r, i.e. a union of half-strips and two
    wedges, and f* is the max of the parts' conjugates."""
    Vi, _ = rat_q.combineDen(obj.VN, obj.VD)
    own = [j + 1 for j, r in enumerate(obj.F) if r[0] == k or r[1] == k]
    segs = []
    rays = []                                      # [a b] vertex indices
    nIn = []
    for e in own:
        a, b = obj.E[e - 1][0], obj.E[e - 1][1]
        d = [Vi[b - 1][0] - Vi[a - 1][0], Vi[b - 1][1] - Vi[a - 1][1]]
        if obj.F[e - 1][0] == k:
            n = [-d[1], d[0]]
        else:
            n = [d[1], -d[0]]                      # the face is on n's side
        nIn.append(n)
        if obj.E[e - 1][2] == 0:
            rays.append([a, b])
        else:
            segs.append([a, b])
    if len(rays) != 2:
        error('PLQ:conjQ:nonConvexPiece',
              'face %d is non-convex and unbounded with %d boundary rays; the cover needs '
              'exactly two.' % (k, len(rays)))
    dA = [Vi[rays[0][1] - 1][i] - Vi[rays[0][0] - 1][i] for i in range(2)]
    dB = [Vi[rays[1][1] - 1][i] - Vi[rays[1][0] - 1][i] for i in range(2)]
    r = None
    for c in ((1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2)):
        cand = [c[0] * dA[0] + c[1] * dB[0], c[0] * dA[1] + c[1] * dB[1]]
        if all(dot(n, cand) > 0 for n in nIn):
            g = max(1, gcd(cand[0], cand[1]))
            r = [cand[0] // g, cand[1] // g]
            break
    if r is None:
        error('PLQ:conjQ:nonConvexPiece',
              'face %d is non-convex and unbounded, and no direction between its rays points '
              'into it from every edge, so it is not swept by one direction and the cover does '
              'not apply.' % k)
    q = obj.fx(k)
    # Every new vertex is an old one plus the INTEGER direction r, so the parts stay exact.
    parts = []
    for s in segs:
        p = obj.Vx(s[0])
        w = obj.Vx(s[1])
        parts.append(onePart([p, w, [p[0] + r[0], p[1] + r[1]], [w[0] + r[0], w[1] + r[1]]],
                             [[1, 2, 1], [1, 3, 0], [2, 4, 0]],
                             [(p[0] + w[0]) / 2 + r[0], (p[1] + w[1]) / 2 + r[1]], q))
    for i in range(2):
        A = obj.Vx(rays[i][0])
        B = obj.Vx(rays[i][1])
        parts.append(onePart([A, B, [A[0] + r[0], A[1] + r[1]]], [[1, 2, 0], [1, 3, 0]],
                             [(A[0] + B[0]) / 2 + r[0], (A[1] + B[1]) / 2 + r[1]], q))
    g = None
    dom = []
    for part in parts:
        gi, di = conjPieceQ(part, 1)
        g = gi if g is None else maxQ(g, gi)
        dom = dom + di
    return g, dom


def onePart(V, E, ip, q):
    """A one-face QuaPol on V/E with F set so the face is on the side of each edge that holds the
    interior point ip -- column 1 is the face on an edge's LEFT."""
    F = []
    for e in E:
        a = V[e[0] - 1]
        d = [V[e[1] - 1][0] - a[0], V[e[1] - 1][1] - a[1]]
        if d[0] * (ip[1] - a[1]) - d[1] * (ip[0] - a[0]) > 0:
            F.append([1, 0])
        else:
            F.append([0, 1])
    return QuaPol(V, E, [q], F)


# ------------------------------------------------------------------------------------------------
def linealityConj(obj, k, sh, Hn, fN, fD, Vi, vd):
    """The conjugate of a convex piece that contains a whole LINE along which q has no curvature
    (a half-plane or a slab with H d = 0): dom f* is THIN -- the line <s,d> = <L,d> -- and on it
    f* is the conjugate of f restricted to the cross-section. Returns (g, isThin).

    PORT NOTE: MATLAB hands the cross-section the NUMERATOR row fN (dropping fD), which is the
    function only when fD = 1; the port passes fN/fD. See PORT_NOTES.md.
    """
    if len(sh['hp']) < 1 or (sh['rays'] and len(sh['rays']) != 2):
        return None, False
    N = [hp[0:2] for hp in sh['hp']]
    g0 = max(1, gcd(abs(N[0][0]), abs(N[0][1])))
    p = [N[0][0] // g0, N[0][1] // g0]
    for i in range(1, len(N)):
        if N[i][0] * p[1] - N[i][1] * p[0] != 0:
            return None, False                       # two directions: no line in the piece
    d = [-p[1], p[0]]
    if any(x != 0 for x in mv(Hn, d)):
        return None, False                           # curvature along d
    # the cross-section through a vertex v, along p: x = v + t p, each row a*x + c >= 0 a bound on t
    v = [Fraction(x, vd) for x in Vi[sh['vs'][0] - 1]]
    lo = None
    hi = None
    for hp in sh['hp']:
        a = hp[0:2]
        c = hp[2]
        coef = dot(a, p)
        rhs = -(dot(a, v) + c)                       # coef * t >= rhs
        if coef > 0:
            b = Fraction(rhs) / coef
            lo = b if lo is None else max(lo, b)
        else:
            b = Fraction(rhs) / coef
            hi = b if hi is None else min(hi, b)
    row = [Fraction(x, fD) for x in fN]

    def at(t):
        return [v[0] + t * p[0], v[1] + t * p[1]]
    if lo is not None and hi is not None:
        if lo == hi:
            Pp = QuaPol([at(lo)], [], [row], [])
        else:
            Pp = QuaPol([at(lo), at(hi)], [[1, 2, 1]], [row], [[0, 0]])
    elif lo is not None:
        Pp = QuaPol([at(lo), at(lo + 1)], [[1, 2, 0]], [row], [[0, 0]])
    elif hi is not None:
        Pp = QuaPol([at(hi), at(hi - 1)], [[1, 2, 0]], [row], [[0, 0]])
    else:
        return None, False                           # no bound across d: the plane, caseA's
    gP, _ = conjPieceQ(Pp, 1)
    Ln = [fN[7], fN[8]]
    eq = [fD * d[0], fD * d[1], -dot(Ln, d)]         # fD <s,d> - <Ln,d> = 0
    eqRow = rat_q.conic([0, 0, 0] + eq) + [sgnOf(eq, 0), 0]
    cells = []
    for f in range(gP.nf):
        con = [list(gP.EcQ[r[0] - 1]) + [r[1], r[2]] for r in gP.FC[f]]
        cells.append({'num': gP.fN[f], 'den': gP.fD[f], 'con': con + [eqRow]})
    return assembleQuaConCells(cells), True

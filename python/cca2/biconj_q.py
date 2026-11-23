"""biconj_q -- the biconjugate f** = cl co f (the closed convex envelope), computed EXACTLY over
the rationals (port of biconjQ.m).

DIRECT, not conj(conj(f)) -- ALGORITHM.md has the principle and the 436 s measurement behind it.
Returns an AlgCon: a QuaCon (rational faces), an AlgCrv (a bridge or tie face), or an EnvelopeMin
(eval only) where the subdivision is declined; kind() says which.

What is decided here, in order:
  * a domain of dimension < 2 -> thinEnvelope / the chain envelopes;
  * f CONVEX on its own domain -> f itself, as a QuaCon;
  * two convex pieces -> the assembled envelope of their conjugates (may decline);
  * a piece that curves UPWARD along one of its own edges -> the atoms route (biconjAtoms and the
    subdivision behind it);
  * otherwise every piece is EDGE-CONCAVE and co f is the affine interpolant over the LOWER HULL
    of the lifted vertices (and lifted recession rays): a hull of rational points in R^3.
"""
import os
from math import gcd

from . import rat_q
from .assemble_qua_con_cells import assembleQuaConCells
from .conj_q import conjQ, assertNoClosedEdges
from .errors import error, CCA2Error
from .is_convex_q import isConvexQ
from .qua_pol import QuaPol
from .rat_q import sign, dot, mv, quad, lcm

# MATLAB's `global BICONJ_ASSERT`: >= 1 refutes a caller's fIsConvex flag instead of trusting it.
BICONJ_ASSERT = None


def _assertLevel():
    if BICONJ_ASSERT is not None:
        return BICONJ_ASSERT
    try:
        return int(os.environ.get('BICONJ_ASSERT', '0') or 0)
    except ValueError:
        return 0


def biconjQ(obj, opts=None):
    if not isinstance(obj, QuaPol):
        error('PLQ:biconjQ:input',
              'biconjQ takes a QuaPol (quadratic on polyhedral); got a %s.' % type(obj).__name__)
    assertNoClosedEdges(obj.E, 'biconjQ')
    obj.assertOperable()
    obj.assertExact()
    obj.assertContinuous()

    nf = len(obj.fN)

    # ---- a domain of dimension < 2 ---------------------------------------------------------------
    if obj.dom['dim'] < 2:
        if isChainSubdivision(obj):
            assertChainEdgeConcave(obj)
            if isChain(obj):
                return concaveEnvelope(obj)         # non-collinear: a 2-D hull
            return collinearChainEnvelope(obj)      # collinear: a 1-D hull along its own span
        return thinEnvelope(obj)

    def Hof(k):
        r = obj.fN[k]
        return [[r[4], r[5]], [r[5], r[6]]]

    allPSD = all(rat_q.isPSD2(Hof(k)) for k in range(nf))

    flagged = obj.fIsConvex is True
    if flagged:
        assertFlagNotRefuted(obj)
    if allPSD and (nf == 1 or flagged or isConvexQ(obj)[0]):
        if obj.nv == 0:
            return assembleQuaConCells([{'num': obj.fN[0], 'den': obj.fD[0], 'con': []}])
        return quaPolAsQuaCon(obj)
    if obj.fIsConvex is True and not allPSD:
        error('PLQ:biconjQ:notConvex',
              'fIsConvex is asserted true but some piece has a Hessian that is not positive '
              'semidefinite, which is a free necessary condition for convexity.')

    # ---- two convex pieces: the assembled envelope -----------------------------------------------
    assemblyWhy = ''
    if allPSD and nf == 2:
        try:
            return assembledEnvelope(obj)
        except CCA2Error as ME:
            if not isDeclined(ME):
                raise
            assemblyWhy = ' The assembly route declined first: %s -- %s' % (
                ME.identifier, ME.message[:200])

    Vi, vd = rat_q.combineDen(obj.VN, obj.VD)
    for k in range(1, nf + 1):
        Hk = Hof(k - 1)
        if not any(r[0] == k or r[1] == k for r in obj.F):
            # A face owning no edge is the whole plane, and it reaches this loop only because it
            # is not convex (a convex single piece took the short-circuit above): f then runs
            # DOWNWARD at infinity, has no affine minorant, and co f is -infinity everywhere --
            # the same answer as the cone case below, so the same identifier (review #17).
            error('PLQ:biconjQ:minusInfinity',
                  'piece %d owns no edge, so it is the whole plane, and it is not convex (a '
                  'convex single piece would have taken the short-circuit above): f runs '
                  'DOWNWARD at infinity, has no affine minorant, and co f is -infinity '
                  'everywhere. That is the right answer and there is nowhere to put it.' % k)
        vsk, rays, bnd = pieceGeometryB(obj, k)
        if not bnd:
            Lk = [obj.fN[k - 1][7], obj.fN[k - 1][8]]
            diverges = False
            for d in rays:
                a = quad(d, Hk, d)
                if a < 0:
                    diverges = True
                    break
                if a == 0:
                    for vi in vsk:
                        v = Vi[vi - 1]
                        rate = dot([x + vd * y for x, y in zip(mv(Hk, v), Lk)], d)
                        if rate < 0:
                            diverges = True
                            break
                    if diverges:
                        break
            if not diverges and len(rays) >= 2:
                a = quad(rays[0], Hk, rays[0])
                b = quad(rays[0], Hk, rays[1])
                c = quad(rays[1], Hk, rays[1])
                if b < 0 and b * b - a * c > 0:
                    diverges = True
            if diverges:
                error('PLQ:biconjQ:minusInfinity',
                      'piece %d is unbounded and q runs DOWNWARD along one of its recession '
                      'directions, so co f is -infinity there. That is the right answer and '
                      'there is nowhere to put it.' % k)

        own = [j + 1 for j, r in enumerate(obj.F) if r[0] == k or r[1] == k]
        for j in own:
            a, b = obj.E[j - 1][0], obj.E[j - 1][1]
            va = obj.Vx(a)
            vb = obj.Vx(b)
            dn, _ = rat_q.fromDouble([vb[0] - va[0], vb[1] - va[1]])
            if quad(dn, Hk, dn) > 0:
                if allPSD:
                    error('PLQ:biconjQ:convexPieceCoupling',
                          'every piece is convex but f is not, so co f must BRIDGE the gradient '
                          'jump across a shared edge (piece %d curves upward along one of its '
                          'own edges, so its vertex values do not determine the bridge). The '
                          'bridge face carries the square root of a quadratic, which no CCA2 '
                          'type can hold -- see DECISIONS.md 2026-09-05.%s' % (k, assemblyWhy))
                return _atomsRoute(obj, opts)

        if bnd:
            if len(vsk) < 3:
                error('PLQ:biconjQ:degenerate',
                      'piece %d has %d vertices; a 2-D bounded face needs at least three.'
                      % (k, len(vsk)))
        elif len(vsk) + len(rays) < 3:
            error('PLQ:biconjQ:degenerate',
                  'piece %d has %d vertices and %d rays; a 2-D face needs three generators.'
                  % (k, len(vsk), len(rays)))

    return concaveEnvelope(obj)


def _atomsRoute(obj, opts):
    """A piece curves upward along one of its own edges: the envelope couples atoms."""
    from .biconj_atoms import biconjAtoms
    from .envelope_prune import envelopePrune
    h = biconjAtoms(obj)
    if opts is None:
        opts = {}
    if len(h.pairs) == 1 and not h.tris:
        try:
            h = subdividePair(envelopePrune(h.pairs[0]), hullRows(obj))
        except CCA2Error as ME:
            if not isDeclined(ME):
                raise
    elif not h.tris:
        try:
            from .envelope_subdivision_alg import envelopeSubdivisionAlg
            from .assemble_alg_crv import assembleAlgCrv
            envs = [envelopePrune(p) for p in h.pairs]
            h = assembleAlgCrv(envelopeSubdivisionAlg(envs, hullRows(obj), [], opts)[0])
        except CCA2Error as ME:
            if not isDeclined(ME):
                raise
    else:
        try:
            from .envelope_subdivision_alg import envelopeSubdivisionAlg
            from .assemble_alg_crv import assembleAlgCrv
            from .triangle_face import triangleFace
            from .cut_triangles import cutTriangles
            envs = [envelopePrune(p) for p in h.pairs]
            tfs = [triangleFace(t) for t in h.tris]
            if len(envs) == 1:
                envs = envs[0]
            sub = envelopeSubdivisionAlg(envs, hullRows(obj), [], opts)[0]
            h = assembleAlgCrv(cutTriangles(sub, tfs)[0])
        except CCA2Error as ME:
            if not isDeclined(ME):
                raise
    return h


def assembledEnvelope(obj):
    from .envelope_assemble import envelopeAssemble
    from .envelope_prune import envelopePrune
    pieces = splitPieces(obj)
    gA, _ = conjQ(pieces[0])
    gB, _ = conjQ(pieces[1])
    return subdividePair(envelopePrune(envelopeAssemble(gA, gB)[0]), hullRows(obj))


def subdividePair(env, dom):
    from .envelope_subdivision import envelopeSubdivision
    from .envelope_subdivision_alg import envelopeSubdivisionAlg
    from .assemble_alg_crv import assembleAlgCrv
    if hasIrrationalFace(env):
        try:
            m, _ = envelopeSubdivisionAlg(env, dom)
            if any(x['kind'] == 'bridge' for x in m):
                return assembleAlgCrv(m)
        except CCA2Error as ME:
            if ME.identifier not in ('envelopeSubdivisionAlg:faceKind',
                                     'envelopeSubdivisionAlg:arrangementSize'):
                raise
            return unsubdivided(env)
    try:
        return envelopeSubdivision(env)[0]
    except CCA2Error as ME:
        if ME.identifier == 'envelopeSubdivision:bridgeCompetes':
            try:
                return assembleAlgCrv(envelopeSubdivisionAlg(env)[0])
            except CCA2Error as ME2:
                if ME2.identifier != 'envelopeSubdivisionAlg:faceKind':
                    raise
                return unsubdivided(env)
        if ME.identifier in ('envelopeSubdivision:bridgeFaces',
                             'envelopeSubdivision:arrangementSize'):
            return unsubdivided(env)
        raise


def convhull2(pts):
    """MATLAB ``convhull(x, y)``: indices (1-based) of the hull, counter-clockwise, CLOSED (first
    index repeated at the end), starting at the lowest-x point (lowest y among ties). Exact:
    the points are Fractions / ints. Collinear boundary points are left out."""
    n = len(pts)
    order = sorted(range(n), key=lambda i: (pts[i][0], pts[i][1]))

    def cross(o, a, b):
        return ((pts[a][0] - pts[o][0]) * (pts[b][1] - pts[o][1])
                - (pts[a][1] - pts[o][1]) * (pts[b][0] - pts[o][0]))
    lower = []
    for i in order:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], i) <= 0:
            lower.pop()
        lower.append(i)
    upper = []
    for i in reversed(order):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], i) <= 0:
            upper.pop()
        upper.append(i)
    hull = lower[:-1] + upper[:-1]
    return [i + 1 for i in hull] + [hull[0] + 1]


def hullRows(obj):
    """The half-planes of conv(dom f) for a BOUNDED domain, as primitive integer rows; empty when
    the domain is unbounded."""
    R = []
    if any(e[2] == 0 for e in obj.E):
        return R
    V = [obj.Vx(i) for i in range(1, len(obj.VN) + 1)]
    k = convhull2(V)
    for i in range(len(k) - 1):
        a, b = k[i], k[i + 1]
        pN, pD = obj.VN[a - 1], obj.VD[a - 1]
        qN, qD = obj.VN[b - 1], obj.VD[b - 1]
        n = [-(qN[1] * pD - pN[1] * qD), qN[0] * pD - pN[0] * qD]
        r = [n[0] * pD, n[1] * pD, -dot(n, pN)]
        g = gcd(gcd(r[0], r[1]), r[2])
        R.append([x // g for x in r])
    return R


def hasIrrationalFace(env):
    from .envelope_rational_corners import envelopeRationalCorners
    env = envelopeRationalCorners(env)
    return bool(env['bridge'])


def isDeclined(ME):
    return ME.identifier.startswith(('envelope', 'bridgeFace', 'bridgeArcCut', 'cornerFace',
                                     'triangleFace', 'cutTriangles'))


def unsubdivided(env):
    from .envelope_min import EnvelopeMin
    return EnvelopeMin([env], [], [])


def splitPieces(obj):
    """The two pieces of a two-face QuaPol as one-face QuaPols."""
    P = []
    for k in (1, 2):
        own = [j for j, r in enumerate(obj.F) if r[0] == k or r[1] == k]
        vs = sorted(set(v for j in own for v in obj.E[j][0:2] if v > 0))
        idx = {v: i + 1 for i, v in enumerate(vs)}
        V = [obj.Vx(v) for v in vs]
        E = [[idx[obj.E[j][0]], idx[obj.E[j][1]], obj.E[j][2]] for j in own]
        F = [[1 if x == k else 0 for x in obj.F[j]] for j in own]
        P.append(QuaPol(V, E, [obj.fx(k)], F))
    return P


def assertFlagNotRefuted(obj):
    if _assertLevel() < 1:
        return
    tf, why = isConvexQ(obj)
    if not tf:
        error('PLQ:biconjQ:notConvex',
              'fIsConvex is asserted true but the exact test refutes it: %s. (BICONJ_ASSERT is '
              'on; unset it to trust the flag again.)' % why)


def concaveEnvelope(obj):
    """co f for EDGE-CONCAVE pieces: the lower convex hull of the lifted vertices {(v, q(v))} and
    lifted recession rays, every predicate an integer sign test."""
    Vi, vd = rat_q.combineDen(obj.VN, obj.VD)
    nf = len(obj.fN)
    chain = isChain(obj)

    Lden = 1
    for k in range(nf):
        Lden = lcm(Lden, obj.fD[k])
    Dl = 2 * Lden * vd ** 2

    lifted = {}
    R3 = []
    for k in range(1, nf + 1):
        fN, fD = obj.faceQ(k)
        Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
        Ln = [fN[7], fN[8]]
        kn = fN[9]
        scale = Lden // fD
        if chain:
            e = obj.E[k - 1]
            if e[2] != 0:
                vsk = [e[0], e[1]]
                raysk = []
            else:
                vsk = [e[0]]
                raysk = [[Vi[e[1] - 1][0] - Vi[e[0] - 1][0], Vi[e[1] - 1][1] - Vi[e[0] - 1][1]]]
        else:
            vsk, raysk, _ = pieceGeometryB(obj, k)
        for d in raysk:
            rho = None
            for vi in vsk:
                v = Vi[vi - 1]
                ri = dot([x + vd * y for x, y in zip(mv(Hn, v), Ln)], d)
                if rho is None or ri < rho:
                    rho = ri
            R3.append([2 * Lden * vd * d[0], 2 * Lden * vd * d[1], 2 * scale * rho])
        for vi in vsk:
            v = Vi[vi - 1]
            qv = scale * (quad(v, Hn, v) + 2 * vd * dot(Ln, v) + 2 * vd ** 2 * kn)
            if vi not in lifted or qv < lifted[vi]:
                lifted[vi] = qv

    vs = sorted(lifted)
    m = len(vs)
    R3 = [list(r) for r in sorted(set(tuple(r) for r in R3))]
    nr = len(R3)
    if m + nr < 3 or m < 1:
        error('PLQ:biconjQ:degenerate',
              'the domain has %d vertices and %d recession rays; a 2-D one needs three '
              'generators.' % (m, nr))
    Q3 = [[2 * Lden * vd * Vi[v - 1][0], 2 * Lden * vd * Vi[v - 1][1], lifted[v]] for v in vs]
    G = Q3 + R3
    isP = [True] * m + [False] * nr
    ng = m + nr

    facets = set()                                # rows [a b c d] with a X + b Y + c Z = d, c > 0
    for i in range(ng):
        for j in range(i + 1, ng):
            for k in range(j + 1, ng):
                idx = [i, j, k]
                pts = [t for t in idx if isP[t]]
                if not pts:
                    continue                      # three rays fix no offset: not a facet
                base = pts[0]
                dirs = []
                for t in idx:
                    if t == base:
                        continue
                    if isP[t]:
                        dirs.append([G[t][c] - G[base][c] for c in range(3)])
                    else:
                        dirs.append(G[t])
                u, w = dirs
                n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2],
                     u[0] * w[1] - u[1] * w[0]]
                if n[2] == 0:
                    continue                      # vertical plane: not the graph of a function
                if n[2] < 0:
                    n = [-x for x in n]           # orient upward, so "above" is n.p >= d
                d = dot(n, G[base])
                below = any(dot(n, q) - d < 0 for q in Q3)
                if not below:
                    below = any(dot(n, r) < 0 for r in R3)
                if below:
                    continue                      # undercut by a vertex or at infinity
                g = gcd(gcd(abs(n[0]), abs(n[1])), gcd(abs(n[2]), abs(d)))
                if g == 0:
                    g = 1
                facets.add((n[0] // g, n[1] // g, n[2] // g, d // g))
    facets = sorted(facets)
    if not facets:
        error('PLQ:biconjQ:noFacet',
              'no lower-hull facet was found, which cannot happen for three non-collinear '
              'vertices.')
    nfac = len(facets)

    dom = hullHalfPlanes(Vi, vd, vs, [r[0:2] for r in R3])

    def face(t):
        a, b, c, d = facets[t]
        return [0, 0, 0, 0, 0, 0, 0, -a * Dl, -b * Dl, d], Dl * c

    cells = []
    for t in range(nfac):
        num, den = face(t)
        rows = [list(r) for r in dom]
        for s in range(nfac):
            if s == t:
                continue
            n2, d2 = face(s)
            dn, _ = rat_q.sub(num, den, n2, d2)
            if all(x == 0 for x in dn):
                continue
            rows.append([2 * dn[7], 2 * dn[8], 2 * dn[9]])
        if not rat_q.feasible2(rows, True):
            continue                              # this facet never wins inside P
        con = [rat_q.conic([0, 0, 0] + r) + [sgnOfB(r, +1)] for r in rows]
        cells.append({'num': num, 'den': den, 'con': con})
    return assembleQuaConCells(cells)


def hullHalfPlanes(Vi, vd, vs, dirs=None):
    """The half-planes of conv(vertices) + cone(dirs), as primitive integer rows on ACTUAL points."""
    dirs = dirs or []
    rows = []
    m = len(vs)
    r = len(dirs)
    for a in range(m):
        for b in range(m + r):
            p = Vi[vs[a] - 1]
            if b < m:
                if a == b:
                    continue
                q = Vi[vs[b] - 1]
                n = [q[1] - p[1], -(q[0] - p[0])]
            else:
                e = dirs[b - m]
                n = [e[1], -e[0]]
            if n[0] == 0 and n[1] == 0:
                continue                          # coincident vertices, or a null direction
            c0 = -dot(n, p)
            sgn = 0
            ok = True
            for i in range(m):
                v = Vi[vs[i] - 1]
                t = n[0] * v[0] + n[1] * v[1] + c0
                if t == 0:
                    continue
                if sgn == 0:
                    sgn = sign(t)
                elif sign(t) != sgn:
                    ok = False
                    break
            if ok:
                for i in range(r):
                    t = dot(n, dirs[i])
                    if t == 0:
                        continue
                    if sgn == 0:
                        sgn = sign(t)
                    elif sign(t) != sgn:
                        ok = False
                        break
            if not ok or sgn == 0:
                continue                          # generators on both sides: not a facet
            row = [sgn * vd * n[0], sgn * vd * n[1], sgn * c0]
            g = gcdRow(row)
            row = [x // g for x in row]
            if row not in rows:
                rows.append(row)
    return rows


def gcdRow(row):
    g = 0
    for x in row:
        if x != 0:
            g = gcd(g, abs(x))
    return g if g != 0 else 1


def polygonHalfPlanes(obj, Vi, vd, vs, k, rays=None):
    rows = []
    for j in [j for j, r in enumerate(obj.F) if r[0] == k or r[1] == k]:
        p = Vi[obj.E[j][0] - 1]
        q = Vi[obj.E[j][1] - 1]
        n = [q[1] - p[1], -(q[0] - p[0])]
        c0 = -dot(n, p)
        sgn = 0
        for vi in vs:
            v = Vi[vi - 1]
            t = n[0] * v[0] + n[1] * v[1] + c0
            if t != 0:
                sgn = sign(t)
                break
        if sgn == 0 and rays is not None:
            for ray in rays:
                t = dot(n, ray)
                if t != 0:
                    sgn = sign(t)
                    break
        if sgn == 0:
            if obj.F[j][1] == k:
                sgn = 1
            elif obj.F[j][0] == k:
                sgn = -1
            else:
                continue                          # genuinely degenerate: skip the edge
        rows.append([sgn * vd * n[0], sgn * vd * n[1], sgn * c0])
    return rows


def sgnOfB(row, want):
    r = [0, 0, 0, row[0], row[1], row[2]]
    nz = next((x for x in r if x != 0), None)
    if nz is None:
        error('ratQ:zeroConic', 'a constraint with all-zero coefficients names no half-plane.')
    return -want if nz < 0 else want


def quaPolAsQuaCon(obj):
    """f is convex, so co f = f: the same faces, as a QuaCon. A NON-CONVEX piece is refused:
    its region is not an intersection of half-planes, so it cannot be one face, and this route
    only returns f unchanged -- it does not compute the envelope of a non-convex piece."""
    Vi, vd = rat_q.combineDen(obj.VN, obj.VD)
    cells = []
    for k in range(1, len(obj.fN) + 1):
        vs, rays, _ = pieceGeometryB(obj, k)
        rows = polygonHalfPlanes(obj, Vi, vd, vs, k, rays)
        for r in rows:
            for vi in vs:
                v = Vi[vi - 1]
                if r[0] * v[0] + r[1] * v[1] + vd * r[2] < 0:
                    error('PLQ:biconjQ:nonConvexPiece',
                          'piece %d is NOT CONVEX, so its region is not an intersection of '
                          'half-planes and cannot be one face of a QuaCon. This route only '
                          'returns f unchanged (one face per piece); the convex envelope of a '
                          'non-convex piece is a different function (on the convex hull of the '
                          'piece), and this route does not compute it, so it declines.' % k)
            for q in rays:
                if r[0] * q[0] + r[1] * q[1] < 0:
                    error('PLQ:biconjQ:nonConvexPiece',
                          'piece %d recedes along (%d,%d), leaving one of its own half-planes.'
                          % (k, q[0], q[1]))
        con = [rat_q.conic([0, 0, 0] + r) + [sgnOfB(r, +1)] for r in rows]
        cells.append({'num': obj.fN[k - 1], 'den': obj.fD[k - 1], 'con': con})
    return assembleQuaConCells(cells)


def _allFzero(obj):
    return bool(obj.F) and not any(x != 0 for r in obj.F for x in r)


def isChainSubdivision(obj):
    return _allFzero(obj) and len(obj.E) >= 2 and obj.nv >= 3


def _collinearExact(obj):
    Vi, _ = rat_q.combineDen(obj.VN, obj.VD)
    A = Vi[0]
    d = None
    for p in Vi[1:]:
        c = [p[0] - A[0], p[1] - A[1]]
        if c[0] != 0 or c[1] != 0:
            d = c
            break
    if d is None:
        return True
    return all(d[0] * (p[1] - A[1]) - d[1] * (p[0] - A[0]) == 0 for p in Vi)


def isChain(obj):
    """A chain of edges that is NOT collinear (decided exactly; MATLAB screens with sqrt(eps))."""
    if not _allFzero(obj):
        return False
    if len(obj.E) < 2:
        return False
    if obj.nv < 3:
        return False
    return not _collinearExact(obj)


def assertChainEdgeConcave(obj):
    Vi, _ = rat_q.combineDen(obj.VN, obj.VD)
    for j in range(1, len(obj.E) + 1):
        fN, _ = obj.faceQ(j)
        Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
        e = obj.E[j - 1]
        d = [Vi[e[1] - 1][0] - Vi[e[0] - 1][0], Vi[e[1] - 1][1] - Vi[e[0] - 1][1]]
        if quad(d, Hn, d) > 0:
            error('PLQ:biconjQ:notImplemented',
                  'chain edge %d curves UPWARD along itself, so its lifted arc dips below its '
                  'own chord and the envelope is not determined by the vertex values.' % j)


def pieceGeometryB(obj, k):
    """(vs, rays, bounded): the vertices (1-based, sorted), integer ray directions of face k."""
    own = [j for j, r in enumerate(obj.F) if r[0] == k or r[1] == k]
    if not own:
        error('PLQ:biconjQ:noFace', 'face %d has no edges.' % k)
    Vi, _ = rat_q.combineDen(obj.VN, obj.VD)
    vs = []
    rays = []
    for j in own:
        e = obj.E[j]
        if e[2] != 0:
            vs += [e[0], e[1]]
        else:
            vs += [e[0]]
            rays.append([Vi[e[1] - 1][0] - Vi[e[0] - 1][0], Vi[e[1] - 1][1] - Vi[e[0] - 1][1]])
    vs = sorted(set(vs))
    return vs, rays, not rays


def collinearChainEnvelope(obj):
    """A chain whose vertices all lie on ONE line: the 1-D lower hull along that line, as thin
    cells (an equality row for the line, two inequalities for the span)."""
    Vi, vd = rat_q.combineDen(obj.VN, obj.VD)
    ne = len(obj.E)
    A = Vi[0]
    d = None
    for p in Vi[1:]:
        cand = [p[0] - A[0], p[1] - A[1]]
        if cand[0] != 0 or cand[1] != 0:
            d = cand
            break
    if d is None:
        error('PLQ:biconjQ:degenerate', 'every vertex of the chain is the same point.')
    g = gcd(abs(d[0]), abs(d[1]))
    if g > 1:
        d = [d[0] // g, d[1] // g]

    Lden = 1
    for j in range(ne):
        Lden = lcm(Lden, obj.fD[j])
    Dl = 2 * Lden * vd ** 2
    lifted = {}
    for j in range(1, ne + 1):
        fN, fD = obj.faceQ(j)
        Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
        Ln = [fN[7], fN[8]]
        kn = fN[9]
        scale = Lden // fD
        for e in range(2):
            iv = obj.E[j - 1][e]
            v = Vi[iv - 1]
            qv = scale * (quad(v, Hn, v) + 2 * vd * dot(Ln, v) + 2 * vd ** 2 * kn)
            if iv not in lifted or qv < lifted[iv]:
                lifted[iv] = qv

    vs = sorted(lifted)
    TZ = [(dot([Vi[v - 1][0] - A[0], Vi[v - 1][1] - A[1]], d), lifted[v]) for v in vs]
    # MATLAB: [T, ord] = sort(T) -- stable, so ties keep the vertex order
    TZ = [TZ[i] for i in sorted(range(len(TZ)), key=lambda i: TZ[i][0])]
    T = [t for t, _ in TZ]
    Z = [z for _, z in TZ]
    st = []
    for i in range(len(T)):
        while len(st) >= 2:
            o = st[-2]
            a1 = st[-1]
            cr = (T[a1] - T[o]) * (Z[i] - Z[o]) - (Z[a1] - Z[o]) * (T[i] - T[o])
            if cr <= 0:
                st.pop()
            else:
                break
        st.append(i)
    if len(st) < 2:
        error('PLQ:biconjQ:degenerate',
              'the collinear chain produced fewer than two hull points.')

    n0 = [d[1], -d[0]]                              # normal to the line
    c0 = -dot(n0, A)
    onLine = [vd * n0[0], vd * n0[1], c0]
    Ad = dot(A, d)
    cells = []
    for k in range(len(st) - 1):
        p = st[k]
        q = st[k + 1]
        DT = T[q] - T[p]
        DZ = Z[q] - Z[p]
        num = [0, 0, 0, 0, 0, 0, 0, DZ * vd * d[0], DZ * vd * d[1], Z[p] * DT - DZ * (Ad + T[p])]
        den = Dl * DT
        lo = [vd * d[0], vd * d[1], -(Ad + T[p])]
        hi = [-vd * d[0], -vd * d[1], Ad + T[q]]
        con = [rat_q.conic([0, 0, 0] + onLine) + [0],
               rat_q.conic([0, 0, 0] + lo) + [sgnOfB(lo, +1)],
               rat_q.conic([0, 0, 0] + hi) + [sgnOfB(hi, +1)]]
        cells.append({'num': num, 'den': den, 'con': con})
    return assembleQuaConCells(cells)


def thinEnvelope(obj):
    """The envelope on a NEEDLE (a point) or a single SEGMENT: a thin QuaCon."""
    fN, fD = obj.faceQ(1)
    Hn = [[fN[4], fN[5]], [fN[5], fN[6]]]
    Ln = [fN[7], fN[8]]
    kn = fN[9]
    Vi, vd = rat_q.combineDen(obj.VN, obj.VD)

    if not obj.E:
        p = Vi[0]
        qp = quad(p, Hn, p) + 2 * vd * dot(Ln, p) + 2 * vd ** 2 * kn     # over 2*fD*vd^2
        r1 = [vd, 0, -p[0]]
        r2 = [0, vd, -p[1]]
        con = [rat_q.conic([0, 0, 0] + r1) + [0], rat_q.conic([0, 0, 0] + r2) + [0]]
        return assembleQuaConCells([{'num': [0, 0, 0, 0, 0, 0, 0, 0, 0, qp],
                                     'den': 2 * fD * vd ** 2, 'con': con}])

    if len(obj.E) != 1:
        error('PLQ:biconjQ:chainNotImplemented',
              'a chain of %d edges reached thinEnvelope, which handles only a single segment. '
              'Chains are routed before this point, so this is a defect in that routing.'
              % len(obj.E))
    if obj.E[0][2] == 0:
        error('PLQ:biconjQ:unbounded',
              'the domain is a RAY. Where q curves down along it the envelope runs to '
              '-infinity, which is a correct answer with nowhere to be stored.')

    a = Vi[obj.E[0][0] - 1]
    b = Vi[obj.E[0][1] - 1]
    d = [b[0] - a[0], b[1] - a[1]]
    alpha = quad(d, Hn, d)
    n = [d[1], -d[0]]
    c0 = -dot(n, a)
    rows = [vd * n[0], vd * n[1], c0]                                   # ON the line
    lo = [vd * d[0], vd * d[1], -dot(d, a)]
    hi = [-vd * d[0], -vd * d[1], dot(d, b)]
    con = [rat_q.conic([0, 0, 0] + rows) + [0],
           rat_q.conic([0, 0, 0] + lo) + [sgnOfB(lo, +1)],
           rat_q.conic([0, 0, 0] + hi) + [sgnOfB(hi, +1)]]
    if alpha >= 0:
        num, den = fN, fD
    else:
        qa = quad(a, Hn, a) + 2 * vd * dot(Ln, a) + 2 * vd ** 2 * kn
        qb = quad(b, Hn, b) + 2 * vd * dot(Ln, b) + 2 * vd ** 2 * kn
        dd = dot(d, d)
        num = [0, 0, 0, 0, 0, 0, 0, (qb - qa) * vd * d[0], (qb - qa) * vd * d[1],
               qa * dd - (qb - qa) * dot(a, d)]
        den = 2 * fD * vd ** 2 * dd
    return assembleQuaConCells([{'num': num, 'den': den, 'con': con}])

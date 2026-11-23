"""biconj_atoms -- the envelope of a QuaPol as a MIN over atom pairs (port of biconjAtoms.m).

f is cut into ATOMS -- convex pieces, the edges and rays of the non-convex pieces along which q
curves upward, and the vertices nothing else holds -- each a small exact QuaPol with an exact
conjugate. co f is then the min over PAIRS of atoms of the pair's assembled envelope, plus the
contact TRIANGLES where three dual forms tie at one isolated dual point.

Which pairs are needed is decided here in DOUBLES (the forms of f*, their ties, the arcs on which
two forms are both active): that screen only prunes the list of pairs -- every pair kept is then
assembled exactly -- so a tolerance here costs at most an unnecessary pair, never a wrong value.
Returns an EnvelopeMin(pairs, tris, atoms); with mode 'triangles', just the contact triangles.
"""
import math

import numpy as np

from . import rat_q
from .conj_q import conjQ
from .envelope_min import EnvelopeMin
from .errors import error
from .qua_pol import QuaPol

INF = math.inf


def biconjAtoms(f, mode=None):
    if not isinstance(f, QuaPol):
        error('biconjAtoms:input', 'expects a QuaPol; got %s.' % type(f).__name__)
    if mode is not None and mode != 'triangles':
        error('biconjAtoms:mode', 'the only mode is \'triangles\'; got %r.' % (mode,))
    from .envelope_assemble import envelopeAssemble

    D = []
    if any(e[2] == 0 for e in f.E):
        _, D = conjQ(f)
        D = [list(r) for r in sorted(set(tuple(r) for r in D))]
        if any(r[0] == 0 and r[1] == 0 and r[2] < 0 for r in D):
            error('biconjAtoms:emptyDomain',
                  'dom f* is empty: f has no affine minorant, so co f is -infinity everywhere.')

    atoms, forms = collectAtoms(f)
    if len(atoms) == 1 and D:
        g1, _ = conjQ(atoms[0]['q'])
        tris = [_dropForms(t) for t in tripleFaces(forms, D)]
        return EnvelopeMin([envelopeAssemble(g1, g1, D)[0]], tris, atoms)
    if len(atoms) < 2:
        error('biconjAtoms:oneAtom',
              'the input reduces to %d atom(s); a single atom is convex, so biconjQ\'s convex '
              'path is the one to take.' % len(atoms))

    tris = tripleFaces(forms, D)
    ab, tris = activeAtomPairs(atoms, forms, tris, D)
    if mode is not None:
        return tris

    g = {}
    for a in sorted(set(x for r in ab for x in r)):
        g[a], _ = conjQ(atoms[a - 1]['q'])
    pairs = [envelopeAssemble(g[r[0]], g[r[1]], D)[0] for r in ab]
    return EnvelopeMin(pairs, tris, atoms)


def _dropForms(t):
    t = dict(t)
    t.pop('forms', None)
    return t


# ------------------------------------------------------------------------------------------------
def activeAtomPairs(atoms, forms, tris, D):
    nA = len(atoms)
    holds = formHolders(atoms, forms)
    req = []                                   # each: {'pairs': [[A, B], ...], 'any': atoms}
    nF = len(forms)
    for i in range(nF):
        for j in range(i + 1, nF):
            if sweepsNoArea(forms[i], forms[j]):
                continue
            if not tieArc(i, j, forms, D):
                continue
            req.append(requirement([holds[i], holds[j]], nA))
    for i in range(nF):
        if forms[i]['kind'] == 'piece' or (D and forms[i]['kind'] == 'edge'):
            req.append({'pairs': [], 'any': list(holds[i])})
    keepT = [True] * len(tris)
    for t, tr in enumerate(tris):
        if not triangleActive(tr, forms, D):
            keepT[t] = False
            continue
        r = requirement([holds[c - 1] for c in tr['forms']], nA)
        if r['pairs'] or r['any']:
            keepT[t] = False                   # two atoms hold all three: the pair attains it
            req.append(r)
    tris = [_dropForms(t) for t, k in zip(tris, keepT) if k]
    if not req:
        ab = [[a, b] for a in range(1, nA + 1) for b in range(a + 1, nA + 1)]
        return ab, tris                        # no arc tie found at all: keep everything (safe)
    return coverRequirements(req, nA), tris


def requirement(H, nA):
    all1 = set(H[0])
    for h in H[1:]:
        all1 &= set(h)
    P = []
    for A in range(1, nA + 1):
        for B in range(A + 1, nA + 1):
            if all((A in h) or (B in h) for h in H):
                P.append([A, B])
    return {'pairs': P, 'any': sorted(all1)}


def coverRequirements(req, nA):
    cand = [[a, b] for a in range(1, nA + 1) for b in range(a + 1, nA + 1)]
    sat = [[(c in r['pairs']) or (c[0] in r['any']) or (c[1] in r['any']) for c in cand]
           for r in req]
    if any(not any(row) for row in sat):
        error('biconjAtoms:uncovered',
              'a tie no atom pair holds -- the form holders are inconsistent.')
    ab = []
    open_ = [True] * len(req)
    while any(open_):
        score = [sum(1 for r in range(len(req)) if open_[r] and sat[r][c])
                 for c in range(len(cand))]
        c = score.index(max(score))
        ab.append(cand[c])
        open_ = [o and not sat[r][c] for r, o in enumerate(open_)]
    return sorted(ab)


def formHolders(atoms, forms):
    return [[a + 1 for a, at in enumerate(atoms) if atomHolds(at, fm)] for fm in forms]


def _holdsValues(at, fm):
    """An atom that holds a form's support geometrically holds the FORM only if it carries the
    same values there (lsc input: across a jump, two faces on one edge or vertex are different
    forms). Atoms and forms without exact data (hand-built in tests) are judged by geometry."""
    if at.get('fq') is None or fm.get('fq') is None or fm.get('sup') is None:
        return True
    return _agree(at['fq'], fm['fq'], fm['sup'])


def _norm(v):
    return math.hypot(float(v[0]), float(v[1]))


def atomHolds(at, fm):
    return _atomHoldsGeom(at, fm) and _holdsValues(at, fm)


def _atomHoldsGeom(at, fm):
    tol = 1e-12
    if fm['kind'] == 'vertex':
        pts = [fm['u']]
        dirs = []
    elif fm['kind'] == 'edge':
        if math.isinf(fm['tmax']):
            pts = [fm['u']]
            dirs = [fm['d']]
        else:
            pts = [fm['u'], [fm['u'][0] + fm['d'][0], fm['u'][1] + fm['d'][1]]]
            dirs = []
    else:
        return at['kind'] == 'piece' and _sameRegion(at['region'], fm['poly'])
    k = at['kind']
    if k == 'vertex':
        return not dirs and all(_norm([p[0] - at['u'][0], p[1] - at['u'][1]]) < tol for p in pts)
    if k == 'edge':
        return not dirs and all(onSegment(p, at['u'], at['w']) for p in pts)
    if k == 'ray':
        dr = at['dr']                          # a ray runs THROUGH its second vertex
        return (all(onRay(p, at['u'], dr) for p in pts)
                and all(parallelSame(d, dr) for d in dirs))
    P = at['region']
    tf = all(inPiece(p, P) for p in pts)
    if tf and dirs:
        if not isinstance(P, dict):
            return False
        R = P['rows']
        m = max(1.0, max(abs(x) for r in R for x in r))
        for d in dirs:
            nd = max(1.0, _norm(d))
            if not all(r[0] * d[0] + r[1] * d[1] >= -1e-9 * m * nd for r in R):
                return False
    return tf


def _sameRegion(P, Q):
    if isinstance(P, dict) != isinstance(Q, dict):
        return False
    if isinstance(P, dict):
        return P['rows'] == Q['rows']
    return P == Q


def onSegment(p, u, w):
    d = [w[0] - u[0], w[1] - u[1]]
    L2 = d[0] * d[0] + d[1] * d[1]
    t = ((p[0] - u[0]) * d[0] + (p[1] - u[1]) * d[1]) / L2
    return (t >= -1e-12 and t <= 1 + 1e-12
            and _norm([u[0] + t * d[0] - p[0], u[1] + t * d[1] - p[1]]) <= 1e-12 * max(1.0, _norm(p)))


def onRay(p, u, d):
    t = ((p[0] - u[0]) * d[0] + (p[1] - u[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1])
    return t >= -1e-12 and \
        _norm([u[0] + t * d[0] - p[0], u[1] + t * d[1] - p[1]]) <= 1e-12 * max(1.0, _norm(p))


def parallelSame(a, b):
    return (abs(a[0] * b[1] - a[1] * b[0]) <= 1e-12 * _norm(a) * _norm(b)
            and a[0] * b[0] + a[1] * b[1] > 0)


def sweepsNoArea(f1, f2):
    k = (f1['kind'], f2['kind'])
    if 'piece' in k:
        return False
    if k == ('vertex', 'vertex'):
        return True
    if f1['kind'] == 'vertex':
        f1, f2 = f2, f1                        # f1 is an edge now

    def onLine(p):
        return abs(f1['d'][0] * (p[1] - f1['u'][1]) - f1['d'][1] * (p[0] - f1['u'][0])) <= \
            1e-12 * _norm(f1['d']) * max(1.0, _norm([p[0] - f1['u'][0], p[1] - f1['u'][1]]))
    if f2['kind'] == 'vertex':
        return onLine(f2['u'])
    return onLine(f2['u']) and onLine([f2['u'][0] + f2['d'][0], f2['u'][1] + f2['d'][1]])


def tieArc(i, j, forms, D):
    """Is there an arc of the tie {g_i = g_j} on which both forms are active (attain f*)?"""
    Q = conicOf(forms[i], forms[j])
    for par in conicParams(Q):
        bp = []
        for c in range(len(forms)):
            for r in domRows(forms[c]):
                bp += realRoots(compose(par, [0, 0, 0] + list(r)))
            if c != i and c != j:
                bp += realRoots(compose(par, conicOf(forms[c], forms[i])))
        for r in D:
            bp += realRoots(compose(par, [0, 0, 0] + list(r)))
        if par['kind'] == 'rat':
            bp += realRoots(par['Dn'])
        bp = sorted(set(bp))
        if not bp:
            samp = [0.0]
        else:
            samp = ([bp[0] - 1 - abs(bp[0])]
                    + [(a + b) / 2 for a, b in zip(bp[:-1], bp[1:])]
                    + [bp[-1] + 1 + abs(bp[-1])])
        for t in samp:
            s, ok = pointAt(par, t)
            if ok and bothActive(s, i, j, forms, D):
                return True
    return False


def bothActive(s, i, j, forms, D):
    if not (inDom(s, domRows(forms[i]), -1) and inDom(s, domRows(forms[j]), -1)
            and inDom(s, D, -1)):
        return False
    vi = gval(forms[i], s)
    for c in range(len(forms)):
        if c == i or c == j or not inDom(s, domRows(forms[c]), +1):
            continue
        if gval(forms[c], s) > vi + 1e-9 * max(1.0, abs(vi)):
            return False
    return True


def triangleActive(tr, forms, D):
    s = tr['s']
    tf = inDom(s, D, -1)
    v = tr['c']
    for c in range(len(forms)):
        if (c + 1) in tr['forms']:
            continue
        if not tf:
            return tf
        if inDom(s, domRows(forms[c]), +1) and \
                gval(forms[c], s) > v + 1e-9 * max(1.0, abs(v)):
            tf = False
    return tf


def inDom(s, R, strictness):
    if not R:
        return True
    tol = 1e-9 * max(1.0, max(abs(x) for r in R for x in r)) * max(1.0, _norm(s))
    return all(r[0] * s[0] + r[1] * s[1] + r[2] >= -strictness * tol for r in R)


def domRows(fm):
    """The half-planes of s on which a form is ATTAINED (its cell of dom f*), in doubles."""
    if 'domRows' in fm:
        return fm['domRows']
    if fm['kind'] == 'vertex':
        R = []
    elif fm['kind'] == 'edge':
        al2 = 2 * fm['al']
        r = [fm['d'][0] / al2, fm['d'][1] / al2, -fm['be'] / al2]        # t(s) >= 0
        R = [r]
        if math.isfinite(fm['tmax']):
            R.append([-r[0], -r[1], fm['tmax'] - r[2]])
    else:
        Hi = np.linalg.inv(np.array(fm['H'], dtype=float))
        P = fm['poly']                          # m(s) = Hi (s - L) in the piece
        if isinstance(P, dict):
            Rm = [list(r) for r in P['rows']]
        else:
            n = len(P)
            Rm = []
            for e in range(n):
                a = P[e]
                b = P[(e + 1) % n]
                nr = [-(b[1] - a[1]), b[0] - a[0]]
                Rm.append([nr[0], nr[1], -(nr[0] * a[0] + nr[1] * a[1])])
        Rm = np.array(Rm, dtype=float).reshape(-1, 3)
        L = np.array(fm['L'], dtype=float)
        A = Rm[:, 0:2] @ Hi
        R = [[float(A[k, 0]), float(A[k, 1]), float(Rm[k, 2] - A[k, :] @ L)]
             for k in range(Rm.shape[0])]
    fm['domRows'] = R
    return R


def _roots(c):
    """MATLAB roots(c) as a list of complex numbers."""
    c = [float(x) for x in c]
    while c and c[0] == 0:
        c = c[1:]
    if len(c) < 2:
        return []
    return [complex(z) for z in np.roots(c)]


def conicParams(Q):
    """Rational / linear parametrisations of the real points of the conic Q (doubles)."""
    P = []
    sc = max(abs(x) for x in Q)
    if sc == 0:
        return P
    Q = [x / sc for x in Q]
    A = np.array([[Q[0], Q[1] / 2], [Q[1] / 2, Q[2]]])
    g = np.array([Q[3], Q[4]])
    f0 = Q[5]
    if np.linalg.norm(A, 2) <= 1e-12:
        if np.linalg.norm(g) <= 1e-12:
            return P
        gg = float(g @ g)
        P.append(lineP([-f0 * g[0] / gg, -f0 * g[1] / gg], [-g[1], g[0]]))
        return P
    M = np.array([[A[0, 0], A[0, 1], g[0] / 2], [A[1, 0], A[1, 1], g[1] / 2],
                  [g[0] / 2, g[1] / 2, f0]])
    if abs(np.linalg.det(M)) > 1e-10 * max(1.0, np.linalg.norm(M, 2)) ** 3:
        p0 = realPoint(A, g, f0)
        if p0 is None:
            return P
        G = 2 * A @ p0 + g
        Dn = [A[1, 1], 2 * A[0, 1], A[0, 0]]
        N = [G[1], G[0]]
        P.append({'kind': 'rat',
                  'S1': padd([p0[0] * x for x in Dn], [-x for x in N]),
                  'S2': padd([p0[1] * x for x in Dn], [-N[0], -N[1], 0.0]),
                  'Dn': [float(x) for x in Dn], 'p': None, 'w': None})
        return P
    dA = float(np.linalg.det(A))
    if dA > 1e-12:
        return P                                # a single point
    if dA < -1e-12:                             # two crossing lines
        c0 = -np.linalg.solve(A, g) / 2
        l, V = np.linalg.eigh(A)
        o = np.argsort(-l, kind='stable')
        V = V[:, o]
        l = l[o]
        for sg in (1, -1):
            w = V[:, 0] / math.sqrt(l[0]) + sg * V[:, 1] / math.sqrt(-l[1])
            P.append(lineP([float(c0[0]), float(c0[1])], [float(w[0]), float(w[1])]))
        return P
    l, V = np.linalg.eigh(A)
    k = int(np.argmax(np.abs(l)))
    lk = float(l[k])
    v = V[:, k]
    vp = [-float(v[1]), float(v[0])]
    for r in _roots([lk, float(g @ v), f0]):    # parallel lines v.s = r
        if abs(r.imag) > 1e-9 * max(1.0, abs(r)):
            continue
        P.append(lineP([r.real * float(v[0]), r.real * float(v[1])], vp))
    return P


def lineP(p0, w):
    return {'kind': 'line', 'S1': None, 'S2': None, 'Dn': None,
            'p': [float(p0[0]), float(p0[1])], 'w': [float(w[0]), float(w[1])]}


def _null(A):
    """MATLAB null(A): an orthonormal basis of the null space, as columns."""
    U, S, Vt = np.linalg.svd(A)
    tol = max(A.shape) * (S[0] if len(S) else 0.0) * np.finfo(float).eps
    r = int(np.sum(S > tol))
    return Vt[r:, :].T


def realPoint(A, g, f0):
    l, V = np.linalg.eigh(A)
    base = [np.zeros(2)]
    if abs(np.linalg.det(A)) > 1e-12:
        base.append(-np.linalg.solve(A, g) / 2)
    dirs = np.hstack([V, np.eye(2), _null(A)])
    for p in base:
        for k in range(dirs.shape[1]):
            w = dirs[:, k]
            c = [float(w @ A @ w), float((2 * A @ p + g) @ w), float(p @ A @ p + g @ p + f0)]
            t = [z.real for z in _roots(c) if abs(z.imag) <= 1e-9 * max(1.0, abs(z))]
            if t:
                # the root of SMALLEST magnitude: along a parabola's axis w'Aw is round-off
                # (~1e-19), and its other root is a spurious ~1e34 point that wrecks the
                # parametrisation -- tieArc then missed real tie arcs (review 2026-11-19 #9c)
                return p + min(t, key=abs) * w
    return None


def _conv(a, b):
    out = [0.0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def compose(par, Q):
    """The polynomial in the parameter t that Q becomes along the parametrised curve."""
    if par['kind'] == 'line':
        A = np.array([[Q[0], Q[1] / 2], [Q[1] / 2, Q[2]]], dtype=float)
        g = np.array([Q[3], Q[4]], dtype=float)
        p = np.array(par['p'])
        w = np.array(par['w'])
        return [float(w @ A @ w), float((2 * A @ p + g) @ w), float(p @ A @ p + g @ p + Q[5])]
    S1, S2, Dn = par['S1'], par['S2'], par['Dn']
    return padd(padd(padd([Q[0] * x for x in _conv(S1, S1)], [Q[1] * x for x in _conv(S1, S2)]),
                     [Q[2] * x for x in _conv(S2, S2)]),
                padd(_conv(padd([Q[3] * x for x in S1], [Q[4] * x for x in S2]), Dn),
                     [Q[5] * x for x in _conv(Dn, Dn)]))


def realRoots(c):
    c = trimLead(c)
    if len(c) < 2:
        return []
    # non-finite coefficients make the eigenvalue run produce NaN roots; they are not
    # breakpoints (review 2026-11-19 #6), and sampling at them silently degrades tieArc
    return [z.real for z in _roots(c) if math.isfinite(z.real) and abs(z.imag) <= 1e-6 * max(1.0, abs(z))]


def _polyval(p, t):
    v = 0.0
    for c in p:
        v = v * t + c
    return v


def pointAt(par, t):
    if par['kind'] == 'line':
        return [par['p'][0] + t * par['w'][0], par['p'][1] + t * par['w'][1]], True
    dn = _polyval(par['Dn'], t)
    ok = abs(dn) > 1e-12 * max(1.0, max(abs(x) for x in par['Dn'])) * max(1.0, t * t)
    if dn == 0:
        return [math.nan, math.nan], False
    return [_polyval(par['S1'], t) / dn, _polyval(par['S2'], t) / dn], ok


# ------------------------------------------------------------------------------------------------
def collectAtoms(f):
    """The atoms of f and the dual FORMS (vertex / edge / piece candidates of f*)."""
    atoms = []
    forms = []
    cvx = []                                   # the convex-piece atoms' regions (see inPiece)
    edgeAt = []                                # [u w] of each edge atom, for dedup
    rayAt = []                                 # [u d] of each RAY atom, d unit
    vtx = []                                   # every vertex that gets a form: [point, value,
    #                                            face, vertex index]
    Vf = [[float(x) for x in f.Vx(i)] for i in range(1, len(f.VN) + 1)]

    for k in range(1, f.nf + 1):
        n, d = f.faceQ(k)
        H = [[n[4] / d, n[5] / d], [n[5] / d, n[6] / d]]
        L = [n[7] / d, n[8] / d]
        k0 = n[9] / d
        isCvx = n[4] >= 0 and n[6] >= 0 and n[4] * n[6] - n[5] ** 2 >= 0   # Sylvester on integers

        def q(x, H=H, L=L, k0=k0):
            return (0.5 * ((x[0] * H[0][0] + x[1] * H[1][0]) * x[0]
                           + (x[0] * H[0][1] + x[1] * H[1][1]) * x[1])
                    + x[0] * L[0] + x[1] * L[1] + k0)

        own = [j + 1 for j, r in enumerate(f.F) if r[0] == k or r[1] == k]
        if any(f.E[j - 1][2] == 0 for j in own):
            ed = []                             # rows [u, w-u, isRay] + vertex indices
            for j in own:
                iu, iw = f.E[j - 1][0], f.E[j - 1][1]
                u = Vf[iu - 1]
                w = Vf[iw - 1]
                ed.append({'u': u, 'dv': [w[0] - u[0], w[1] - u[1]],
                           'ray': f.E[j - 1][2] == 0, 'iu': iu, 'iw': iw})
            poly = {'rows': pieceRows(ed)}
            eidx = own
        else:
            order, eidx = piecePolygon(f, k)
            poly = [Vf[i - 1] for i in order]
            m = len(order)
            ed = []
            for i in range(m):
                u = poly[i]
                w = poly[(i + 1) % m]
                ed.append({'u': u, 'dv': [w[0] - u[0], w[1] - u[1]], 'ray': False,
                           'iu': order[i], 'iw': order[(i + 1) % m]})
        fq = (n, d)                             # this face's exact quadratic (lsc checks)
        for e in ed:
            addVertex(vtx, e['u'], q(e['u']), k, e['iu'], fq, f)
            if not e['ray']:
                w = [e['u'][0] + e['dv'][0], e['u'][1] + e['dv'][1]]
                addVertex(vtx, w, q(w), k, e['iw'], fq, f)

        if isCvx:
            cvx.append((poly, fq))
            atoms.append({'kind': 'piece', 'q': singlePiece(f, k, eidx), 'region': poly,
                          'fq': fq})
            if n[4] * n[6] - n[5] ** 2 > 0:
                Hn = np.array(H)
                Hi = np.linalg.inv(Hn)
                Ln = np.array(L)
                HiL = Hi @ Ln
                Gn, Gd = formQ('piece', f, k)
                forms.append(mkForm('piece', Hi.tolist(), [-float(HiL[0]), -float(HiL[1])],
                                    float(0.5 * Ln @ Hi @ Ln - k0), None, None, None, None,
                                    H, L, poly, Gn, Gd, fq, None))

        for e in ed:
            u = e['u']
            dv = e['dv']
            isRay = e['ray']
            al = 0.5 * ((dv[0] * H[0][0] + dv[1] * H[1][0]) * dv[0]
                        + (dv[0] * H[0][1] + dv[1] * H[1][1]) * dv[1])
            if al <= 0:
                continue                        # a ray of linear growth only shapes dom f* (D)
            be = ((u[0] * H[0][0] + u[1] * H[1][0] + L[0]) * dv[0]
                  + (u[0] * H[0][1] + u[1] * H[1][1] + L[1]) * dv[1])
            ga = q(u)
            Gn, Gd = formQ('edge', f, k, e['iu'], e['iw'])
            sup = _edgeSupport(f, e['iu'], e['iw'], isRay)
            a2 = 2 * al
            fm = mkForm('edge',
                        [[dv[0] * dv[0] / a2, dv[0] * dv[1] / a2],
                         [dv[1] * dv[0] / a2, dv[1] * dv[1] / a2]],
                        [u[0] - be * dv[0] / a2, u[1] - be * dv[1] / a2],
                        -ga + be ** 2 / (4 * al), u, dv, al, be, None, None, None, Gn, Gd,
                        fq, sup)
            if isRay:
                fm['tmax'] = INF
            forms.append(fm)
            if isCvx:
                continue                        # the convex piece's atom already holds it
            if isRay:
                nd = _norm(dv)
                key = [u[0], u[1], dv[0] / nd, dv[1] / nd]     # a ray is its base and direction
                # the same ray of another face is the same atom only where the two faces AGREE
                # on it (lsc input: across a jump both values are kept, and the min is theirs)
                if any(all(abs(a - b) < 1e-12 for a, b in zip(r, key)) and _agree(rq, fq, sup)
                       for r, rq in rayAt):
                    continue
                rayAt.append((key, fq))
                atoms.append({'kind': 'ray', 'u': u, 'w': None, 'dr': dv,
                              'q': QuaPol([f.Vx(e['iu']), f.Vx(e['iw'])], [[1, 2, 0]],
                                          [f.fx(k)], [[0, 0]]), 'fq': fq, 'sup': sup})
                continue
            w = [u[0] + dv[0], u[1] + dv[1]]
            pr = sorted([u, w])
            key = [pr[0][0], pr[0][1], pr[1][0], pr[1][1]]
            if any(all(abs(a - b) < 1e-12 for a, b in zip(r, key)) and _agree(rq, fq, sup)
                   for r, rq in edgeAt):
                continue
            edgeAt.append((key, fq))
            atoms.append({'kind': 'edge', 'u': u, 'w': w,
                          'q': QuaPol([f.Vx(e['iu']), f.Vx(e['iw'])], [[1, 2, 1]],
                                      [f.fx(k)], [[0, 0]]), 'fq': fq, 'sup': sup})

    for p, val, face, iv, fqv, vex in vtx:
        # the vertex's value is the MIN over the faces at it (lsc input; = the common value when
        # f is continuous there), so an atom covers it only where its own face takes that value
        Gn, Gd = formQ('vertex', f, face, iv)
        sup = [f.Vx(iv)]
        forms.append(mkForm('vertex', [[0.0, 0.0], [0.0, 0.0]], [p[0], p[1]], -val, p, None,
                            None, None, None, None, None, Gn, Gd, fqv, sup))
        covered = any(inPiece(p, P) and _agree(pq, fqv, sup) for P, pq in cvx)
        for r, rq in edgeAt:
            covered = covered or ((_norm([p[0] - r[0], p[1] - r[1]]) < 1e-12
                                   or _norm([p[0] - r[2], p[1] - r[3]]) < 1e-12)
                                  and _agree(rq, fqv, sup))
        for r, rq in rayAt:
            covered = covered or (_norm([p[0] - r[0], p[1] - r[1]]) < 1e-12   # a ray holds its base
                                  and _agree(rq, fqv, sup))
        if not covered:
            atoms.append({'kind': 'vertex', 'u': p,
                          'q': QuaPol([f.Vx(iv)], [], [[0, 0, 0, 0, 0, 0, 0, 0, 0, vex]], []),
                          'fq': fqv})

    keep = []
    for at in atoms:
        if at['kind'] == 'edge' and any(inPiece(at['u'], P) and inPiece(at['w'], P)
                                        and _agree(pq, at['fq'], at['sup']) for P, pq in cvx):
            continue
        keep.append(at)
    return keep, forms


def mkForm(kind, A, b, c, u, d, al, be, H, L, poly, Gn, Gd, fq=None, sup=None):
    # fq: the exact quadratic (num, den) of the face the form comes from; sup: the exact points
    # that decide whether another face takes the same values on the form's support (_agree)
    return {'kind': kind, 'A': A, 'b': b, 'c': c, 'u': u, 'd': d, 'al': al, 'be': be,
            'H': H, 'L': L, 'poly': poly, 'tmax': 1, 'Gn': Gn, 'Gd': Gd, 'fq': fq, 'sup': sup}


def addVertex(vtx, p, v, k, iv, fq, f):
    """One entry per vertex POSITION, carrying the SMALLEST value of the faces at it (lsc input:
    f at a shared vertex is the min over the closed pieces holding it; ties keep the first face,
    so a continuous f is unchanged)."""
    vex = rat_q.evalFaceExact(fq[0], fq[1], f.Vx(iv))
    for t, r in enumerate(vtx):
        if abs(r[0][0] - p[0]) < 1e-12 and abs(r[0][1] - p[1]) < 1e-12:
            if vex < r[5]:
                vtx[t] = (list(p), v, k, iv, fq, vex)
            return
    vtx.append((list(p), v, k, iv, fq, vex))


def _edgeSupport(f, iu, iw, isRay):
    """Three exact points of an edge (both ends and the midpoint) or of a ray (base, the point
    through which it runs, and one more): a quadratic difference vanishing at them vanishes on
    the whole edge or ray (it is a univariate quadratic along the line)."""
    u = f.Vx(iu)
    w = f.Vx(iw)
    if isRay:
        return [u, w, (2 * w[0] - u[0], 2 * w[1] - u[1])]
    return [u, w, ((u[0] + w[0]) / 2, (u[1] + w[1]) / 2)]


def _agree(fa, fb, pts):
    """Do the exact face quadratics fa and fb take the same value at every point of pts?

    This is the exact continuity test behind every place where one atom stands in for another
    face's values (review 2026-11-19 #10, #108): True for a continuous f, so those inputs take
    the old path unchanged."""
    if fa is fb or fa == fb:
        return True
    return all(rat_q.evalFaceExact(fa[0], fa[1], P) == rat_q.evalFaceExact(fb[0], fb[1], P)
               for P in pts)


def formQ(kind, f, k, iu=None, iw=None):
    """The form's conic [a b c d e f] over a common denominator, EXACTLY: (Gn, Gd)."""
    from fractions import Fraction as Fr
    n, d = f.faceQ(k)
    H = [[Fr(n[4], d), Fr(n[5], d)], [Fr(n[5], d), Fr(n[6], d)]]
    L = [Fr(n[7], d), Fr(n[8], d)]
    k0 = Fr(n[9], d)

    def quadf(a, b):
        return sum(a[i] * H[i][j] * b[j] for i in range(2) for j in range(2))

    def qv(u):
        return Fr(1, 2) * quadf(u, u) + L[0] * u[0] + L[1] * u[1] + k0
    if kind == 'vertex':
        u = f.Vx(iu)
        c = [Fr(0), Fr(0), Fr(0), u[0], u[1], -qv(u)]
    elif kind == 'edge':
        u = f.Vx(iu)
        w = f.Vx(iw)
        dv = [w[0] - u[0], w[1] - u[1]]
        al = Fr(1, 2) * quadf(dv, dv)
        Hu = [H[0][0] * u[0] + H[0][1] * u[1] + L[0], H[1][0] * u[0] + H[1][1] * u[1] + L[1]]
        be = Hu[0] * dv[0] + Hu[1] * dv[1]
        i4 = 1 / (4 * al)
        i2 = 1 / (2 * al)
        c = [dv[0] * dv[0] * i4, 2 * dv[0] * dv[1] * i4, dv[1] * dv[1] * i4,
             u[0] - be * dv[0] * i2, u[1] - be * dv[1] * i2, be * be * i4 - qv(u)]
    else:
        dt = H[0][0] * H[1][1] - H[0][1] * H[0][1]
        Hi = [[H[1][1] / dt, -H[0][1] / dt], [-H[0][1] / dt, H[0][0] / dt]]
        HiL = [Hi[0][0] * L[0] + Hi[0][1] * L[1], Hi[1][0] * L[0] + Hi[1][1] * L[1]]
        c = [Fr(1, 2) * Hi[0][0], Hi[0][1], Fr(1, 2) * Hi[1][1], -HiL[0], -HiL[1],
             Fr(1, 2) * (L[0] * HiL[0] + L[1] * HiL[1]) - k0]
    D = 1
    for x in c:
        D = rat_q.lcm(D, x.denominator)
    N = [int(x * D) for x in c]
    return rat_q.canon(N, D)


def piecePolygon(f, k):
    """(order, eidx): the vertex indices of face k walked with the face on the LEFT, and the
    face's edge indices."""
    eidx = [j + 1 for j, r in enumerate(f.F) if r[0] == k or r[1] == k]
    E = []
    for j in eidx:
        a, b = f.E[j - 1][0], f.E[j - 1][1]
        if f.F[j - 1][0] != k:
            a, b = b, a                         # orient every edge with face k on its LEFT
        E.append([a, b])
    order = [E[0][0]]
    nxt = E[0][1]
    while nxt != order[0]:
        order.append(nxt)
        cand = [e[1] for e in E if e[0] == nxt]
        if len(cand) != 1:
            error('biconjAtoms:subdivision', 'face %d is not a simple polygon.' % k)
        nxt = cand[0]
    return order, eidx


def singlePiece(f, k, eidx):
    vs = sorted(set(v for j in eidx for v in f.E[j - 1][0:2]))
    idx = {v: i + 1 for i, v in enumerate(vs)}
    E = [[idx[f.E[j - 1][0]], idx[f.E[j - 1][1]], f.E[j - 1][2]] for j in eidx]
    F = [[1 if x == k else 0 for x in f.F[j - 1]] for j in eidx]
    return QuaPol([f.Vx(v) for v in vs], E, [f.fx(k)], F)


def inPiece(p, P):
    if isinstance(P, dict):
        R = P['rows']
        m = max(1.0, max(abs(x) for r in R for x in r))
        return all(r[0] * p[0] + r[1] * p[1] + r[2] >= -1e-9 * m * max(1.0, _norm(p)) for r in R)
    return inPoly(p, P)


def pieceRows(ed):
    rays = [e['dv'] for e in ed if e['ray']]
    ip = [sum(e['u'][0] for e in ed) / len(ed), sum(e['u'][1] for e in ed) / len(ed)]
    for r in rays:
        nr = _norm(r)
        ip[0] += r[0] / nr
        ip[1] += r[1] / nr
    R = []
    for e in ed:
        nrm = [-e['dv'][1], e['dv'][0]]
        row = [nrm[0], nrm[1], -(nrm[0] * e['u'][0] + nrm[1] * e['u'][1])]
        if row[0] * ip[0] + row[1] * ip[1] + row[2] < 0:
            row = [-x for x in row]
        R.append(row)
    return R


def inPoly(p, P):
    n = len(P)
    m = max(1.0, max(abs(x) for r in P for x in r)) ** 2
    for i in range(n):
        a = P[i]
        b = P[(i + 1) % n]
        if (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]) < -1e-9 * m:
            return False
    return True


# ------------------------------------------------------------------------------------------------
def tripleFaces(forms, D):
    """The contact TRIANGLES: three forms tying at one isolated dual point s0, whose three contact
    points span a triangle of positive area. Each carries the affine face <s0, x> - g(s0)."""
    tris = []
    n = len(forms)
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                F3 = [forms[a], forms[b], forms[c]]
                for s0 in commonTies(F3):
                    M = []
                    ok = True
                    for t in range(3):
                        m, ok1 = contact(F3[t], s0)
                        M.append(m)
                        ok = ok and ok1
                    if not ok:
                        continue
                    ar = ((M[1][0] - M[0][0]) * (M[2][1] - M[0][1])
                          - (M[1][1] - M[0][1]) * (M[2][0] - M[0][0]))
                    if abs(ar) <= 1e-9 * max(1.0, max(abs(x) for r in M for x in r)) ** 2:
                        continue
                    tris.append({'T': M, 's': list(s0), 'c': gval(F3[0], s0),
                                 'G': normalCone(D, s0), 'forms': [a + 1, b + 1, c + 1],
                                 'Gn': [fm['Gn'] for fm in F3 if fm['Gn']],
                                 'Gd': [fm['Gd'] for fm in F3 if fm['Gn']]})
    return tris


def normalCone(D, s0):
    """Generators (columns, as a list of [g1, g2]) of the normal cone of dom f* at s0."""
    if not D:
        return []
    r = [d[0] * s0[0] + d[1] * s0[1] + d[2] for d in D]
    tol = 1e-9 * max(1.0, max(abs(x) for d in D for x in d))
    if any(x < -tol for x in r):
        return []                               # s0 is outside D: the triangle alone
    return [[-d[0], -d[1]] for d, x in zip(D, r) if abs(x) <= tol * max(1.0, _norm(s0))]


def gval(fm, s):
    A = fm['A']
    return (0.5 * (s[0] * (A[0][0] * s[0] + A[0][1] * s[1]) + s[1] * (A[1][0] * s[0] + A[1][1] * s[1]))
            + fm['b'][0] * s[0] + fm['b'][1] * s[1] + fm['c'])


def contact(fm, s):
    if fm['kind'] == 'vertex':
        return list(fm['u']), True
    if fm['kind'] == 'edge':
        t = (s[0] * fm['d'][0] + s[1] * fm['d'][1] - fm['be']) / (2 * fm['al'])
        m = [fm['u'][0] + t * fm['d'][0], fm['u'][1] + t * fm['d'][1]]
        return m, (t >= -1e-9 and t <= fm['tmax'] + 1e-9)        # tmax = inf on a ray
    m = np.linalg.solve(np.array(fm['H'], dtype=float),
                        np.array([s[0] - fm['L'][0], s[1] - fm['L'][1]], dtype=float))
    m = [float(m[0]), float(m[1])]
    return m, inPiece(m, fm['poly'])


def commonTies(F3):
    """The isolated points where the three forms tie pairwise (doubles, Newton-polished)."""
    c1 = conicOf(F3[0], F3[1])
    c2 = conicOf(F3[0], F3[2])
    S = []
    if all(abs(x) < 1e-14 for x in c1) or all(abs(x) < 1e-14 for x in c2):
        return S
    for swap in (False, True):
        if swap:
            p1 = [c1[i] for i in (2, 1, 0, 4, 3, 5)]
            p2 = [c2[i] for i in (2, 1, 0, 4, 3, 5)]
        else:
            p1, p2 = c1, c2
        cand = conicMeetDouble(p1, p2)
        if swap:
            cand = [[s[1], s[0]] for s in cand]
        S += cand
        if cand:
            break
    out = []
    kept = []
    for s in S:
        s = list(s)
        for _ in range(30):
            r = [conicAt(c1, s), conicAt(c2, s)]
            J = [conicGrad(c1, s), conicGrad(c2, s)]
            if _rcond2(J) < 1e-14:
                break
            ds = _solve2(J, [-r[0], -r[1]])
            s = [s[0] + ds[0], s[1] + ds[1]]
            if _norm(ds) <= 1e-15 * max(1.0, _norm(s)):
                break
        sc = max(1.0, _norm(s)) ** 2 * max(1.0, max(abs(x) for x in c1 + c2))
        keep = (math.isfinite(s[0]) and math.isfinite(s[1])
                and abs(conicAt(c1, s)) < 1e-8 * sc and abs(conicAt(c2, s)) < 1e-8 * sc)
        for sj, kj in zip(out, kept):
            if kj and _norm([sj[0] - s[0], sj[1] - s[1]]) < 1e-9 * max(1.0, _norm(s)):
                keep = False
        out.append(s)
        kept.append(keep)
    return [s for s, k in zip(out, kept) if k]


def _rcond2(J):
    from .conic_meet import _rcond2 as rc
    return rc(J)


def _solve2(J, r):
    from .conic_meet import _solve2 as sv
    return sv(J, r)


def conicOf(fa, fb):
    A = [[fa['A'][i][j] - fb['A'][i][j] for j in range(2)] for i in range(2)]
    b = [fa['b'][0] - fb['b'][0], fa['b'][1] - fb['b'][1]]
    return [A[0][0] / 2, A[0][1], A[1][1] / 2, b[0], b[1], fa['c'] - fb['c']]


def conicAt(c, s):
    return (c[0] * s[0] ** 2 + c[1] * s[0] * s[1] + c[2] * s[1] ** 2
            + c[3] * s[0] + c[4] * s[1] + c[5])


def conicGrad(c, s):
    return [2 * c[0] * s[0] + c[1] * s[1] + c[3], c[1] * s[0] + 2 * c[2] * s[1] + c[4]]


def conicMeetDouble(c1, c2):
    """Real common points of two double conics, by eliminating s2 (resultant in s1)."""
    A1 = c1[2]
    B1 = [c1[1], c1[4]]
    C1 = [c1[0], c1[3], c1[5]]
    A2 = c2[2]
    B2 = [c2[1], c2[4]]
    C2 = [c2[0], c2[3], c2[5]]
    if A1 == 0 and A2 == 0:
        R = padd(_conv(B1, C2), [-x for x in _conv(B2, C1)])
    else:
        X = padd([A1 * x for x in C2], [-A2 * x for x in C1])
        Y = padd([A1 * x for x in B2], [-A2 * x for x in B1])
        Z = padd(_conv(B1, C2), [-x for x in _conv(B2, C1)])
        R = padd(_conv(X, X), [-x for x in _conv(Y, Z)])
    R = trimLead(R)
    S = []
    if not R or all(x == 0 for x in R):
        return S
    r = [z.real for z in _roots(R) if abs(z.imag) <= 1e-7 * max(1.0, abs(z))]
    for s1 in r:
        a = A1
        bb = _polyval(B1, s1)
        cc = _polyval(C1, s1)
        if abs(a) > 0:
            s2c = _roots([a, bb, cc])
        elif bb != 0:
            s2c = [complex(-cc / bb)]
        else:
            s2c = []
        a2 = A2
        b2 = _polyval(B2, s1)
        cc2 = _polyval(C2, s1)
        if not s2c or (abs(a) == 0 and bb == 0):
            if abs(a2) > 0:
                s2c = _roots([a2, b2, cc2])
            elif b2 != 0:
                s2c = [complex(-cc2 / b2)]
        for z in s2c:
            if abs(z.imag) <= 1e-7 * max(1.0, abs(z)):
                S.append([s1, z.real])
    return S


def padd(a, b):
    n = max(len(a), len(b))
    a = [0.0] * (n - len(a)) + [float(x) for x in a]
    b = [0.0] * (n - len(b)) + [float(x) for x in b]
    return [x + y for x, y in zip(a, b)]


def trimLead(p):
    if not p:
        return []
    m = max(1.0, max(abs(x) for x in p))
    i = next((k for k, x in enumerate(p) if abs(x) > 1e-12 * m), None)
    if i is None:
        return []
    return list(p[i:])

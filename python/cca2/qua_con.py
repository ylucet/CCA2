"""QuaCon -- a quadratic function on a CONIC subdivision, stored EXACTLY (port of QuaCon.m).

What is stored, and why each is the shape it is (QuaCon.m's header has the full argument):

  fN, fD   the face functions: row k is a numerator in the 10-wide weighted cubic basis over the
           positive integer fD[k], in rat_q.canon form.
  EcQ      the edge conics, rat_q.conic form (primitive, first nonzero entry positive), so two
           edges lie on one curve exactly when the two integer vectors are equal.
  Vname    the vertices as NAMES [edgeA, edgeB, rootIdx] (all 1-based): the rootIdx-th real
           intersection, in conicMeet's canonical order, of those two edge conics. No coordinate
           is ever stored -- a vertex of f* is generically irrational.
  FC       the faces in H-form: FC[k] is a list of [edgeIndex, side, branch]; the face is
           { x : side * Q_edge(x) >= 0 for every row }, side 0 meaning ON the curve (a thin face),
           branch 1 or 2 naming ONE line of a degenerate conic (0 = the whole conic).
  E, F     the combinatorial half: E rows [v1, v2, isSegment], F rows [face on +side, face on
           -side], 1-based with 0 = none.

This is the type `conj` returns. Only `eval` and `vertexCoords` produce doubles, on purpose.
"""
import math

from . import rat_q
from .conic_meet import conicMeet
from .errors import error, CCA2Error


class AlgCon:
    """Root of the result lattice: every `biconj` result is an AlgCon (QuaCon, AlgCrv or
    EnvelopeMin). kind() names the actual class, never a stored flag."""

    def kind(self):
        return type(self).__name__

    def isSubdivided(self):
        return True


class QuaCon(AlgCon):
    def __init__(self, Vname, E, EcQ, fN, fD, F, FC):
        Vname = [list(r) for r in Vname]
        E = [list(r) for r in E]
        EcQ = [list(r) for r in EcQ]
        fN = [list(r) for r in fN]
        fD = list(fD)
        F = [list(r) for r in F]
        FC = [[list(r) for r in R] for R in FC]
        if any(len(r) != 3 for r in Vname):
            error('QuaCon:size', 'Vname must be nv x 3 [edgeA edgeB rootIdx].')
        nv = len(Vname)
        ne = len(E)
        nf = len(fN)
        if any(len(r) != 3 for r in E) or len(EcQ) != ne or len(F) != ne or \
                any(len(r) != 2 for r in F):
            error('QuaCon:size',
                  'E must be ne x 3 and EcQ / F must have one row per edge; got E %d rows, '
                  'EcQ %d rows, F %d rows.' % (len(E), len(EcQ), len(F)))
        if len(fD) != nf or len(FC) != nf:
            error('QuaCon:size',
                  'fD and FC need one entry per face: %d faces, %d denominators, %d constraint '
                  'lists.' % (nf, len(fD), len(FC)))
        if nv > 0 and ne > 0 and max(max(r[0], r[1]) for r in E) > nv:
            error('QuaCon:badEdge', 'E names vertex %d but there are only %d vertices.'
                  % (max(max(r[0], r[1]) for r in E), nv))

        # ---- the exact face functions, canonically ------------------------------------------
        for k in range(nf):
            fN[k], fD[k] = rat_q.canon(fN[k], fD[k])
            if any(x != 0 for x in fN[k][0:4]):
                error('QuaCon:cubicFace',
                      'face %d has a cubic numerator. A QuaCon holds the QUADRATIC output of '
                      'conj; a cubic has no conic level set, so its edges could not be stored '
                      'here at all (see rat_q.diffConic).' % (k + 1))

        # ---- the exact edge conics, canonically ----------------------------------------------
        for j in range(ne):
            EcQ[j] = rat_q.conic(EcQ[j])           # raises ratQ:zeroConic on an all-zero row

        # ---- the vertex NAMES, and that each one can actually be realised --------------------
        for i in range(nv):
            a, b, r = Vname[i]
            if a < 1 or b < 1 or a > ne or b > ne:
                error('QuaCon:badVertexName',
                      'vertex %d names edges %d and %d, but there are %d edges.' % (i + 1, a, b, ne))
            if a == b:
                error('QuaCon:badVertexName',
                      'vertex %d names edge %d twice. A vertex is where TWO curves meet; one '
                      'curve does not name a point.' % (i + 1, a))
            P, info = conicMeet(EcQ[a - 1], EcQ[b - 1])
            if info['degenerate']:
                error('QuaCon:degenerateVertex',
                      'vertex %d names edges %d and %d, whose conics share a component, so they '
                      'meet in infinitely many points and name no vertex.' % (i + 1, a, b))
            if r < 1 or r > len(P):
                error('QuaCon:badRootIndex',
                      'vertex %d asks for real intersection %d of edges %d and %d, which meet '
                      'in %d real points.' % (i + 1, r, a, b, len(P)))

        # ---- the H-form face constraints -----------------------------------------------------
        for k in range(nf):
            Ck = FC[k]
            Ck = [r + [0] if len(r) == 2 else r for r in Ck]     # the pre-branch spelling
            if any(len(r) != 3 for r in Ck):
                error('QuaCon:badConstraint',
                      'FC{%d} must be m x 3 of [edgeIndex, sign, branch] (or m x 2, read as '
                      'branch 0).' % (k + 1))
            if Ck:
                if any(r[0] < 1 or r[0] > ne for r in Ck):
                    error('QuaCon:badConstraint',
                          'FC{%d} names an edge outside 1..%d.' % (k + 1, ne))
                if not all(r[1] in (-1, 0, 1) for r in Ck):
                    error('QuaCon:badConstraint',
                          'FC{%d} carries a side that is not -1, 0 or +1. The side says WHICH '
                          'side of the conic the face lies on and cannot be omitted; 0 means ON '
                          'the curve.' % (k + 1))
                if not all(r[2] in (0, 1, 2) for r in Ck):
                    error('QuaCon:badConstraint',
                          'FC{%d} carries a branch that is not 0, 1 or 2. A branch NAMES one '
                          'line of a degenerate conic; 0 means the whole conic.' % (k + 1))
                for ri, r in enumerate(Ck):
                    if r[2] == 0:
                        continue
                    try:
                        rat_q.pairLineData(EcQ[r[0] - 1])
                    except CCA2Error as ME:
                        error('QuaCon:badConstraint',
                              'FC{%d} row %d names branch %d of edge %d, but that conic is not '
                              'two distinct crossing lines (%s).'
                              % (k + 1, ri + 1, r[2], r[0], ME.identifier))
            FC[k] = Ck

        self.nv = nv
        self.ne = ne
        self.nf = nf
        self.E = E
        self.F = F
        self.fN = fN
        self.fD = fD
        self.EcQ = EcQ
        self.Vname = Vname
        self.FC = FC
        self.fIsConvex = None

    # ------------------------------------------------------------------------------------------
    def vertexCoords(self):
        """Realise the named vertices as double coordinates, ON DEMAND (never stored)."""
        V = []
        for a, b, r in self.Vname:
            P, _ = conicMeet(self.EcQ[a - 1], self.EcQ[b - 1])
            V.append(P[r - 1])
        return V

    def sameEdgeCurve(self, j1, j2):
        return self.EcQ[j1 - 1] == self.EcQ[j2 - 1]

    def edgeKind(self, j):
        return rat_q.conicKind(self.EcQ[j - 1])

    def eval(self, X):
        """Value at each point of X (list of [x, y]) and the 1-based face that supplied it
        (0 where +inf). Doubles on purpose: point location screens with a tolerance, and a
        point takes the face it is MOST inside."""
        n = len(X)
        val = [math.inf] * n
        idx = [0] * n
        if n == 0:
            return val, idx
        best = [math.inf] * n
        xmax = [max(1.0, abs(float(p[0])), abs(float(p[1]))) for p in X]
        tol = math.sqrt(2.0 ** -52)
        Xf = [[float(p[0]), float(p[1])] for p in X]
        for k in range(self.nf):
            Ck = self.FC[k]
            viol = [-math.inf] * n
            for e, side, br in Ck:
                c = self.EcQ[e - 1]
                if br == 0:
                    q = rat_q.evalConic(c, Xf)
                    cm = max(abs(x) for x in c)
                    scale = [max(1.0, cm * xm * xm) for xm in xmax]
                else:
                    q = rat_q.pairLineEval(c, br, Xf)
                    scale = xmax
                if side == 0:
                    viol = [max(v, abs(qi) / s) for v, qi, s in zip(viol, q, scale)]
                else:
                    viol = [max(v, -side * qi / s) for v, qi, s in zip(viol, q, scale)]
            take = [i for i in range(n) if viol[i] <= tol and viol[i] < best[i]]
            if take:
                vals = rat_q.evalFace(self.fN[k], self.fD[k], [Xf[i] for i in take])
                for i, v in zip(take, vals):
                    val[i] = v
                    idx[i] = k + 1
                    best[i] = viol[i]
        return val, idx

    def __repr__(self):
        return 'QuaCon(%d vertices (named), %d edges, %d faces)' % (self.nv, self.ne, self.nf)

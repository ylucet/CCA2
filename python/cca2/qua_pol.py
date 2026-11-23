"""QuaPol -- a piecewise QUAdratic function on a POLyhedral subdivision: the INPUT type
(port of the exact subset of QuaPol.m that conj / biconj reach).

    QuaPol(f)              one quadratic on the whole plane
    QuaPol(V, E, f, F)     V  nv x 2 vertices
                           E  ne x 3 rows [v1 v2 flag], flag 1 = segment, 0 = ray from V[v1]
                              THROUGH V[v2]; vertex indices are 1-based
                           f  one row per face (6 or 10 wide, weighted basis):
                              (c5 x^2 + 2 c6 xy + c7 y^2)/2 + c8 x + c9 y + c10
                           F  ne x 2 rows [left face, right face] of each edge, 1-based, 0 = outside

Numbers may be ints, Fractions (taken exactly) or floats (reconstructed as the rational the
caller wrote, as MATLAB does; a float that is not one makes the object INEXACT, and the exact
operators then refuse it by name -- QuaPol:notExact).
"""
from fractions import Fraction
import math

from . import rat_q
from .errors import error, CCA2Error


class QuaPol:
    def __init__(self, *args):
        if len(args) == 1:
            f = args[0]
            V, E, F = [], [], []
        elif len(args) == 4:
            V, E, f, F = args
        else:
            error('QuaPol:nargin',
                  'Invalid input arguments in QuaPol constructor; need 4 arguments but got %d'
                  % len(args))
        if f and not isinstance(f[0], (list, tuple)):
            f = [f]
        self.f = [[0] * (10 - len(r)) + list(r) for r in f]
        self.V = [list(r) for r in V]
        self.E = [[int(x) for x in r] for r in E]
        self.F = [[int(x) for x in r] for r in F]
        self.fIsConvex = None
        if len(args) == 1:
            self.nv, self.ne, self.nf = 0, 0, 1
        else:
            nv, ne, nf = len(self.V), len(self.E), len(self.f)
            ok = [all(len(r) == 2 for r in self.V),
                  all(len(r) == 3 for r in self.E) and ne == len(self.F),
                  all(len(r) in (6, 10) for r in f) and nf > 0,
                  all(len(r) == 2 for r in self.F) and ne == len(self.F)]
            if not all(ok):
                error('QuaPol:size',
                      'Invalid matrix size in QuaPol constructor. Validation vector: '
                      '[V, E, f, F]=%s' % ([int(x) for x in ok],))
            self.nv, self.ne, self.nf = nv, ne, nf
            if (nv, ne) in ((1, 0), (2, 1)):
                self.nf = 0                          # needle, or one segment / ray
            if self.F and max(max(r) for r in self.F) == 0:
                self.nf = 0                          # a domain of dimension < 2: chain of edges
            if self.nf > 0:
                if ne > 1 and max(max(r) for r in self.F) != self.nf:
                    error('QuaPol:badF',
                          'max index in F must be equal to the number of faces')
                if nv > 1 and max(max(r[0], r[1]) for r in self.E) != nv:
                    error('QuaPol:badE', 'max index in E must be equal to size(V,1)')
        self._dom = None
        self._readExactly()

    # ---- the exact reading of the caller's numbers -------------------------------------------
    def _readExactly(self):
        try:
            N, D = [], []
            for row in self.f:
                n, d = rat_q.fromDouble(row)
                N.append(n)
                D.append(d)
            NV, DV = [], []
            for row in self.V:
                n, d = rat_q.fromDouble(row)
                NV.append(n)
                DV.append(d)
            self.fN, self.fD, self.VN, self.VD = N, D, NV, DV
        except CCA2Error as e:
            if e.identifier not in ('ratQ:notExact', 'ratQ:overflow', 'ratQ:notFinite'):
                raise
            self.fN = self.fD = self.VN = self.VD = None

    def isExact(self):
        return self.fN is not None and (self.nv == 0 or self.VN is not None)

    def assertExact(self):
        if not self.isExact():
            error('PLQ:QuaPol:notExact',
                  'this QuaPol does not carry exact coefficients, so it cannot be used on the '
                  'exact path. It was built from COMPUTED doubles rather than from data a caller '
                  'wrote down. Computing exactly on inexact input yields exactly the wrong '
                  'number, which carries no warning.')

    def faceQ(self, k):
        """Face k (1-based), exactly: (numerator row, denominator)."""
        self.assertExact()
        return list(self.fN[k - 1]), self.fD[k - 1]

    def vertexQ(self, i):
        self.assertExact()
        return list(self.VN[i - 1]), self.VD[i - 1]

    def Vx(self, i):
        """Vertex i (1-based) as exact Fractions."""
        return [Fraction(x, self.VD[i - 1]) for x in self.VN[i - 1]]

    def fx(self, k):
        """Function row k (1-based) as exact Fractions."""
        return [Fraction(x, self.fD[k - 1]) for x in self.fN[k - 1]]

    # ---- the domain ----------------------------------------------------------------------------
    @property
    def dom(self):
        if self._dom is None:
            self._dom = self._createDom()
        return self._dom

    def _createDom(self):
        dom = {'dim': 0, 'isConvex': False}
        if self.nv == 0:
            if self.nf != 1:
                error('QuaPol:dom', 'a function with no vertex must have exactly 1 face')
            return {'dim': 2, 'isConvex': True}
        if self.nv == 1:
            dom['isConvex'] = True
            return dom
        if self.nf == 0:
            dom['dim'] = 1
            if self.nv == 2 and self.ne == 1:
                dom['isConvex'] = True
            elif isCollinear(self.V, math.sqrt(2.0 ** -52)):
                dom['isConvex'] = True
            return dom
        dom['dim'] = 2
        if not any(x == 0 for r in self.F for x in r):
            dom['isConvex'] = True
            return dom
        from .domain_is_convex_d import domainIsConvexD
        dom['isConvex'] = domainIsConvexD(self.V, self.E, self.F)
        return dom

    def isDomBounded(self):
        if self.ne == 0:
            return self.nv == 1
        return all(r[2] != 0 for r in self.E)

    # ---- guards ----------------------------------------------------------------------------------
    def degree(self):
        if not self.f:
            return 0
        tol = math.sqrt(2.0 ** -52)

        def nrm(cols):
            return math.sqrt(sum(float(r[c]) ** 2 for r in self.f for c in cols))
        # MATLAB takes the 2-norm of the block; any norm separates zero from nonzero the same way
        if nrm(range(0, 4)) > tol:
            return 3
        if nrm(range(4, 7)) > tol:
            return 2
        if nrm(range(7, 9)) > tol:
            return 1
        return 0

    def assertOperable(self):
        if self.degree() > 2:
            error('PLQ:op:unsupportedType',
                  'Convex-analysis operators require a quadratic numerator (degree<=2); got '
                  'degree %d. Cubic is allowed for isConvex only.' % self.degree())

    def isContinuous(self):
        """(tf, why): is f CONTINUOUS across every shared edge? Exact, three points per edge.

        d = f_L - f_R restricted to the edge LINE is a univariate quadratic, and it vanishes on
        the whole edge iff it vanishes at three points of it: both endpoints and the midpoint
        (the same condition (b) three-point test as MATLAB isConvexQ.m, whose header carries the
        reference). A degenerate edge (both endpoints at one vertex) leaves one point, and the
        midpoint formula below evaluates to that point, so the single test that remains decides
        d on it.

        A discontinuous f is NOT refused (#108): conj and biconj read it as LOWER SEMICONTINUOUS
        -- at a shared point f is the min over the closed pieces holding it. This is the exact
        test the continuity-dependent shortcuts use to decide their path (isConvexQ's convex
        shortcut runs the same three-point test; biconjAtoms applies it per edge / vertex).
        """
        self.assertExact()
        if self.nf < 2:
            return True, ''
        fN, fD = self.fN, self.fD
        for j, row in enumerate(self.F):
            iL, iR = row
            if iL == 0 or iR == 0 or iL == iR:
                continue
            dn = [a * fD[iR - 1] - b * fD[iL - 1] for a, b in zip(fN[iL - 1], fN[iR - 1])]
            e = self.E[j]
            P1 = self.VN[e[0] - 1]; D1 = self.VD[e[0] - 1]
            P2 = self.VN[e[1] - 1]; D2 = self.VD[e[1] - 1]

            def val(P, D):
                # 2*D^2 * ((n5 x^2)/2 + n6 xy + (n7 y^2)/2 + n8 x + n9 y + n10) at (x,y) = P/D:
                # all integers, zero exactly when the jump value is (isConvexQ.m evalNum).
                return (dn[4] * P[0] * P[0] + 2 * dn[5] * P[0] * P[1] + dn[6] * P[1] * P[1]
                        + 2 * dn[7] * P[0] * D + 2 * dn[8] * P[1] * D + 2 * dn[9] * D * D)
            if (val(P1, D1) != 0 or val(P2, D2) != 0
                    or val((P1[0] * D2 + P2[0] * D1, P1[1] * D2 + P2[1] * D1), 2 * D1 * D2) != 0):
                return False, ('f is discontinuous across edge %d: faces %d and %d do not agree '
                               'on the edge' % (j + 1, iL, iR))
        return True, ''

    # ---- the operators ---------------------------------------------------------------------------
    def conj(self, engine='exact'):
        """The Legendre-Fenchel conjugate f*(s) = sup_x <s,x> - f(x), exactly, as a QuaCon."""
        self.assertOperable()
        e = engine.lower()
        if e == 'exact':
            from .conj_q import conjQ
            return conjQ(self)[0]
        if e in ('cplq', 'symbolic'):
            error('PLQ:conj:engine',
                  'the legacy symbolic engine \'%s\' is MATLAB-only; the Python port carries '
                  'the exact engine.' % engine)
        if e in ('pqp', 'graph'):
            error('PLQ:conj:engine',
                  'Conjugate engine \'%s\' is not implemented yet; use \'exact\'.' % engine)
        error('PLQ:conj:engine', 'Unknown conjugate engine \'%s\'.' % engine)

    def biconj(self, engine='exact', opts=None):
        """The biconjugate f** = cl co f (the closed convex envelope), exactly, as an AlgCon."""
        self.assertOperable()
        e = engine.lower()
        if e == 'exact':
            from .biconj_q import biconjQ
            return biconjQ(self, opts)
        if e in ('cplq', 'symbolic'):
            error('PLQ:biconj:engine',
                  'the legacy symbolic engine \'%s\' is MATLAB-only; the Python port carries '
                  'the exact engine.' % engine)
        error('PLQ:biconj:engine', 'Unknown biconjugate engine \'%s\'.' % engine)

    def __repr__(self):
        return 'QuaPol(nv=%d, ne=%d, nf=%d)' % (self.nv, self.ne, self.nf)


def isCollinear(V, tol):
    """Are all the points of V on one line (to tol)? Doubles: a property of the legacy `dom`."""
    if len(V) < 3:
        return True
    p0 = [float(x) for x in V[0]]
    d = None
    for p in V[1:]:
        u = [float(p[0]) - p0[0], float(p[1]) - p0[1]]
        if math.hypot(u[0], u[1]) > tol:
            d = u
            break
    if d is None:
        return True
    n = math.hypot(d[0], d[1])
    for p in V[1:]:
        u = [float(p[0]) - p0[0], float(p[1]) - p0[1]]
        if abs(d[0] * u[1] - d[1] * u[0]) / n > tol * max(1.0, math.hypot(u[0], u[1])):
            return False
    return True

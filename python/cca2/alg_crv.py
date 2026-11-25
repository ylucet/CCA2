"""AlgCrv -- an AlgCon whose faces may be IMPLICIT (a named root of a polynomial) and whose rows
may compare such values (port of AlgCrv.m).

On top of QuaCon's tables (EcQ conics, FC rows, fN/fD rational faces) it carries:
  FP, FB   face k's implicit polynomial (poly_n dict in x, y, z) and branch index, or None/0
  EcA      table 1: rows "rational face  vs  algebraic value":  side * (q - value) >= 0
  EcI      table 2: rows "value_1 >= value_2" between two algebraic values
  Tri, FT  table 3: contact triangles; FT[k] > 0 makes face k the affine face of triangle FT[k]
An FC row is [index, side, branch, table]; table 0 is EcQ (branch names a line of a pair).
eval takes, at each point, the LARGEST value among the faces whose rows hold there.
"""
import math

from . import rat_q
from .alg_branch import algBranch, fibreVanishes
from .errors import error
from .qua_con import AlgCon

_R2 = math.sqrt(2.0)
_DIRS = [[1, 0], [-1, 0], [0, 1], [0, -1],
         [1 / _R2, 1 / _R2], [1 / _R2, -1 / _R2], [-1 / _R2, 1 / _R2], [-1 / _R2, -1 / _R2]]


class AlgCrv(AlgCon):
    def __init__(self, Vname, E, EcQ, fN, fD, F, FC, FP=None, FB=None, EcA=None, EcI=None,
                 Tri=None, FT=None):
        self.Vname = [list(r) for r in (Vname or [])]
        self.E = [list(r) for r in (E or [])]
        self.EcQ = [list(r) for r in (EcQ or [])]
        self.fN = [list(r) for r in (fN or [])]
        self.fD = list(fD or [])
        self.F = [list(r) for r in (F or [])]
        self.FC = [[list(r) for r in R] for R in FC]
        self.nv = len(self.Vname)
        self.nf = len(self.FC)
        self.ne = len(self.EcQ)
        self.FP = list(FP) if FP is not None else [None] * self.nf
        self.FB = list(FB) if FB is not None else [0] * self.nf
        self.EcA = list(EcA or [])
        self.EcI = list(EcI or [])
        self.Tri = list(Tri or [])
        self.FT = list(FT) if FT is not None else [0] * self.nf
        self.fIsConvex = None

    def eval(self, X):
        """(val, idx): the value at each point and the 1-based face that supplied it.

        A CONTACT TRIANGLE decides where it holds: its plane <s,x> - f*(s) is an affine minorant
        of f** that equals f** on the whole closed triangle, so no other face can be right there
        with a different value. Cutting the triangle out of a face leaves closed remnants that
        can meet the triangle on a ZERO-AREA set (a segment of its edge); such a remnant carries
        its pair's value, which for a discontinuous f (#108 residual, D-132, D-111) is larger
        than f** on that edge and would win the max."""
        n = len(X)
        Xf = [[float(p[0]), float(p[1])] for p in X]
        val = [math.inf] * n
        idx = [0] * n
        tri = [False] * n
        for k in range(self.nf):
            M = self._faceClosurePoints(k, Xf)
            pts = [i for i in range(n) if M[i]]
            if not pts:
                continue
            v = self._faceValues(k, [Xf[i] for i in pts])
            isTri = self.FT[k] > 0
            for i, vi in zip(pts, v):
                if math.isnan(vi):
                    continue
                if tri[i] and not isTri:
                    continue                    # a triangle already decides this point
                if (isTri and not tri[i]) or math.isinf(val[i]) or vi > val[i]:
                    val[i] = vi
                    idx[i] = k + 1
                    tri[i] = tri[i] or isTri
        return val, idx

    def _faceClosurePoints(self, k, X):
        """_facePoints, with a row that is UNDEFINED at a point decided by the face's closure.

        An algebraic row compares a value that is a root of a polynomial in z; where the whole
        fibre vanishes (a BASE POINT: every z is a root, e.g. a 0/0 of a bridge between two atoms
        that share a vertex) the row says nothing, and reading it as "fails" left a hole: D-153's
        domain corner (0,0), where the face around it is worth -2 = f**, was +inf (#108 residual).
        There the face holds iff it holds at points arbitrarily close: its rows are read at x
        nudged by 1e-7 in eight directions (envelopeEval's nudge), and the value is still the
        face's own at x."""
        M, U = self._facePoints(k, X, undecided=True)
        for i in range(len(X)):
            if M[i] or not U[i]:
                continue
            h = 1e-7 * max(1.0, abs(X[i][0]), abs(X[i][1]))
            Y = [[X[i][0] + h * d[0], X[i][1] + h * d[1]] for d in _DIRS]
            if any(self._facePoints(k, Y)):
                M[i] = True
        return M

    def _facePoints(self, k, X, undecided=False):
        M = [True] * len(X)
        U = [False] * len(X)                    # held but for rows undefined at a base point
        for row in self.FC[k]:
            live = [i for i in range(len(X)) if M[i] or U[i]]
            if not live:
                return (M, U) if undecided else M
            P = [X[i] for i in live]
            tbl = row[3] if len(row) >= 4 else 0
            if tbl == 3:
                T = self.Tri[row[0] - 1]['T']
                i0 = row[1]
                j0 = i0 % 3 + 1
                o0 = 6 - i0 - j0
                Ti, Tj, To = T[i0 - 1], T[j0 - 1], T[o0 - 1]
                e = [Tj[0] - Ti[0], Tj[1] - Ti[1]]
                so = rat_q.sign(e[0] * (To[1] - Ti[1]) - e[1] * (To[0] - Ti[0]))
                ne = max(1.0, math.hypot(e[0], e[1]))
                ok = []
                for p in P:
                    cr = (e[0] * (p[1] - Ti[1]) - e[1] * (p[0] - Ti[0])) * so
                    ok.append(not (row[2] * cr < -1e-9 * ne * max(1.0, math.hypot(p[0], p[1]))))
            elif tbl == 2:
                r = self.EcI[row[0] - 1]
                v1 = algBranch(r['FP1'], P, r['FB1'])
                v2 = algBranch(r['FP2'], P, r['FB2'])
                ok = [math.isfinite(a) and math.isfinite(b) and not (row[1] * (a - b) < -1e-9)
                      for a, b in zip(v1, v2)]
                und = [(not math.isfinite(a) and fibreVanishes(r['FP1'], p))
                       or (not math.isfinite(b) and fibreVanishes(r['FP2'], p))
                       for a, b, p in zip(v1, v2, P)]
            elif tbl == 1:
                a = self.EcA[row[0] - 1]
                v = algBranch(a['FP'], P, a['FB'])
                q = rat_q.evalFace(a['num'], a['den'], P)
                ok = [math.isfinite(vi) and not (row[1] * (qi - vi) < -1e-9)
                      for vi, qi in zip(v, q)]
                und = [not math.isfinite(vi) and fibreVanishes(a['FP'], p) for vi, p in zip(v, P)]
            else:
                c = self.EcQ[row[0] - 1]
                br = row[2] if len(row) >= 3 else 0
                if br != 0:
                    v = rat_q.pairLineEval(c, br, P)
                else:
                    v = rat_q.evalConic(c, P)
                if row[1] == 0:
                    ok = [not (abs(vi) > 1e-9) for vi in v]          # an EQUALITY row: thin face
                else:
                    ok = [not (row[1] * vi < -1e-9) for vi in v]
            if tbl not in (1, 2):
                und = [False] * len(live)
            for i, o, u in zip(live, ok, und):
                if M[i]:
                    M[i] = o
                    U[i] = (not o) and u
                else:
                    U[i] = o or u               # an undecided point stays undecided, or fails
        return (M, U) if undecided else M

    def _faceValues(self, k, X):
        if self.FT[k] > 0:
            t = self.Tri[self.FT[k] - 1]
            return [x[0] * t['s'][0] + x[1] * t['s'][1] - t['c'] for x in X]
        if self.FP[k]:
            return algBranch(self.FP[k], X, self.FB[k])
        return rat_q.evalFace(self.fN[k], self.fD[k], X)

    def __repr__(self):
        return 'AlgCrv(%d faces, %d conics, %d algebraic rows)' % (self.nf, self.ne, len(self.EcA))

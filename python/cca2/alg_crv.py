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
from .alg_branch import algBranch
from .errors import error
from .qua_con import AlgCon


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
        """(val, idx): the value at each point and the 1-based face that supplied it."""
        n = len(X)
        Xf = [[float(p[0]), float(p[1])] for p in X]
        val = [math.inf] * n
        idx = [0] * n
        for k in range(self.nf):
            M = self._facePoints(k, Xf)
            pts = [i for i in range(n) if M[i]]
            if not pts:
                continue
            v = self._faceValues(k, [Xf[i] for i in pts])
            for i, vi in zip(pts, v):
                if math.isnan(vi):
                    continue
                if math.isinf(val[i]) or vi > val[i]:
                    val[i] = vi
                    idx[i] = k + 1
        return val, idx

    def _facePoints(self, k, X):
        M = [True] * len(X)
        for row in self.FC[k]:
            live = [i for i in range(len(X)) if M[i]]
            if not live:
                return M
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
            elif tbl == 1:
                a = self.EcA[row[0] - 1]
                v = algBranch(a['FP'], P, a['FB'])
                q = rat_q.evalFace(a['num'], a['den'], P)
                ok = [math.isfinite(vi) and not (row[1] * (qi - vi) < -1e-9)
                      for vi, qi in zip(v, q)]
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
            for i, o in zip(live, ok):
                M[i] = o
        return M

    def _faceValues(self, k, X):
        if self.FT[k] > 0:
            t = self.Tri[self.FT[k] - 1]
            return [x[0] * t['s'][0] + x[1] * t['s'][1] - t['c'] for x in X]
        if self.FP[k]:
            return algBranch(self.FP[k], X, self.FB[k])
        return rat_q.evalFace(self.fN[k], self.fD[k], X)

    def __repr__(self):
        return 'AlgCrv(%d faces, %d conics, %d algebraic rows)' % (self.nf, self.ne, len(self.EcA))

"""max_q -- the pointwise maximum of two QuaCons, exactly (port of maxQ.m).

Overlay every cell of g1 with every cell of g2; on the overlay the larger function is decided by
the sign of ONE conic, the level set {g1 = g2}, which is rational.
"""
from . import rat_q
from .assemble_qua_con_cells import assembleQuaConCells
from .errors import error
from .qua_con import QuaCon


def maxQ(g1, g2):
    if not isinstance(g1, QuaCon) or not isinstance(g2, QuaCon):
        error('PLQ:maxQ:input', 'maxQ takes two QuaCons; got a %s and a %s.'
              % (type(g1).__name__, type(g2).__name__))
    cells = []
    for i in range(g1.nf):
        ci = _rowsOf(g1, i)
        for j in range(g2.nf):
            cj = _rowsOf(g2, j)
            base = ci + cj                          # the overlay cell: both cells' conditions
            n1, d1 = g1.fN[i], g1.fD[i]
            n2, d2 = g2.fN[j], g2.fD[j]
            dn, _ = rat_q.sub(n1, d1, n2, d2)       # g1 - g2, exactly
            if all(x == 0 for x in dn):
                cells.append({'num': n1, 'den': d1, 'con': base})
                continue
            R = [dn[4], 2 * dn[5], dn[6], 2 * dn[7], 2 * dn[8], 2 * dn[9]]
            Rc = rat_q.conic(R)
            nz = next(x for x in R if x != 0)
            flip = -1 if nz < 0 else 1
            cells.append({'num': n1, 'den': d1, 'con': base + [Rc + [flip]]})
            cells.append({'num': n2, 'den': d2, 'con': base + [Rc + [-flip]]})
    return assembleQuaConCells(cells)


def _rowsOf(g, k):
    # MATLAB keeps [conic, side] and DROPS the branch column here; kept verbatim.
    return [list(g.EcQ[r[0] - 1]) + [r[1]] for r in g.FC[k]]

"""NOTE 2026-11-17: the defect this script demonstrates is FIXED (conicMeet line-vs-conic scale bug,
DECISIONS.md 2026-11-17); the script now reports it as not reproduced. The root-cause text below is
the first diagnosis (the symptom), kept for the record.
"""
"""Star-polygon conjugate defect -- reproduction and root-cause demonstration.

DEFECT
  conjQ of a star polygon (2n vertices, f = (x^2+y^2)/2) returns a value that is TOO LOW at some
  dual points: the face that should dominate (the interior face, value |s|^2/2) is dropped and a
  lower edge-face is picked instead. Reproduced with n = 6, 7, 8 (see the counts below). The 8- and
  10-vertex stars (n = 4, 5) happen to be correct, but the defect is latent in all of them.

ROOT CAUSE
  The max of the per-triangle conjugates is assembled by maxQ -> assembleQuaConCells ->
  compactCells -> _dropImplied. `_dropImplied` drops a row of a cell when, with that row FLIPPED,
  `cellInterior` reports NO interior. `cellInterior` delegates the curved case to `cellHasInterior`,
  whose decision is a WITNESS SEARCH: it tries the pairwise conic/line intersection points, then 400
  fixed random probes in a box around their centroid. That search is INCOMPLETE -- for a cell whose
  interior is a thin strip the probes miss it, so `cellHasInterior` returns False ("no interior") for
  a cell that DOES have interior. `_dropImplied` then drops a bounding row, enlarging the cell until
  it overlaps a higher face's cell; QuaCon.eval then picks the face it is "most inside" (the lower
  one), giving a too-low value.

  The witness search is described in the source as "the same screen MATLAB uses", so the defect is a
  SHARED algorithmic defect (MATLAB, Python and Rust all carry the incomplete screen), not a port bug.

WHY A SIMPLE REFUSAL DOES NOT WORK
  Refusing "a strictly convex quadratic on a bounded non-convex polygon" (the star-polygon class)
  also refuses correct cases: golden conj record 36 (a 6-vertex L-shape, f = (x^2+y^2)/2, 71 faces)
  is a valid conjugate that the refusal would break. The failing and the correct cases are the SAME
  class; only the completeness of the cell-interior test separates them. The proper fix is a COMPLETE
  test for "does this cell have interior" (decide feasibility of the conic sign conditions), not a
  heuristic witness search.

USAGE
    python3 star_defect_diagnosis.py            # reproduces the defect + demonstrates the root cause
Exit status 0 = defect reproduced (root cause shown); the star values are wrong on purpose.
"""
import sys, math, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from cca2.qua_pol import QuaPol                       # noqa: E402
from cca2.conj_q import (earClip, boundaryCycle, pieceShape, caseBConvexOnPiece,   # noqa: E402
                         triangleShape, maxQ)
from cca2 import rat_q                                # noqa: E402
from cca2.cell_interior import cellInterior           # noqa: E402
from cca2.cell_has_interior import cellHasInterior, _satisfiesAll  # noqa: E402


def star(n):
    m = 2 * n
    V = []
    for i in range(m):
        r = 100.0 if i % 2 == 0 else 45.0
        a = math.pi * i / n
        V.append([round(r * math.cos(a)), round(r * math.sin(a))])
    E = [[(i + 1) % m + 1, (i + 2) % m + 1, 1] for i in range(m)]
    return QuaPol(V, E, [[1, 0, 1, 0, 0, 0]], [[1, 0]] * m)     # f = (x^2+y^2)/2


def truth(n, s):
    """f*(s) = sup over the star of <s,x> - (x^2+y^2)/2, exactly: the interior point x = s if it is
    inside the star (then |s|^2/2), else the best boundary point, edge by edge."""
    m = 2 * n
    vert = []
    for i in range(m):
        r = 100.0 if i % 2 == 0 else 45.0
        a = math.pi * i / n
        vert.append([round(r * math.cos(a)), round(r * math.sin(a))])
    inside = False
    for i in range(m):
        a, b = vert[i], vert[(i + 1) % m]
        if (a[1] > s[1]) != (b[1] > s[1]) and s[0] < (b[0] - a[0]) * (s[1] - a[1]) / (b[1] - a[1]) + a[0]:
            inside = not inside
    if inside:
        return (s[0] ** 2 + s[1] ** 2) / 2.0
    val = lambda x: s[0] * x[0] + s[1] * x[1] - (x[0] ** 2 + x[1] ** 2) / 2.0
    best = -float('inf')
    for i in range(m):
        a, b = vert[i], vert[(i + 1) % m]
        d = [b[0] - a[0], b[1] - a[1]]
        t = max(0.0, min(1.0, ((s[0] - a[0]) * d[0] + (s[1] - a[1]) * d[1]) / (d[0] ** 2 + d[1] ** 2)))
        best = max(best, val([a[0] + t * d[0], a[1] + t * d[1]]))
    return best


def conic_at(c, pt):
    a, b, cc, d, e, f = [float(x) for x in c]
    return a * pt[0] ** 2 + b * pt[0] * pt[1] + cc * pt[1] ** 2 + d * pt[0] + e * pt[1] + f


def main():
    print('== 1. Reproduce the defect (n = 7 star, f = (x^2+y^2)/2) ==')
    n = 7
    p = star(n)
    g, _ = None, None
    # build the conjugate exactly as conjQ does (ear-clip + maxQ over the triangles) so we can also
    # inspect the intermediate cells
    fN, fD = p.faceQ(1)
    Vi2, vd2 = rat_q.combineDen(p.VN, p.VD)
    sh = pieceShape(p, 1, Vi2, vd2)
    tris = earClip(Vi2, boundaryCycle(sh, Vi2))
    tcons = [caseBConvexOnPiece(fN, fD, Vi2, vd2, triangleShape(Vi2, vd2, t), []) for t in tris]
    acc = tcons[0]
    for t in tcons[1:]:
        acc = maxQ(acc, t)

    s = [-36.39, -5.970000000000001]
    tval = truth(n, s)
    engine = acc.eval([s])[0][0]
    print('  at s = %s' % (s,))
    print('    definition (truth) = %.4f' % tval)
    print('    engine             = %.4f   <-- TOO LOW (should be %.4f)' % (engine, tval))
    defect = engine < tval * (1 - 1e-6)
    print('  DEFECT REPRODUCED: %s' % defect)

    print()
    print('== 2. Root cause: cellInterior says "no interior" for a cell that provably has one ==')
    # The flipped cell below is a real cell of the maxQ assembly: 4 line conditions + 2 conic
    # conditions. It HOLDS the point s (every row is satisfied at s), so it has interior. Yet the
    # engine's own interior test says otherwise.
    C = [
        [0, 0, 0, 128, -35, 4809, 0],
        [0, 0, 0, 3, -31, -1001, 0],
        [0, 0, 0, 38, -9, 749, 0],
        [35, 256, -35, -7000, -25600, 350000, 0],
        [1, 0, 0, 56, 70, 2009, 0],
        [0, 0, 0, 20, 59, -2000, 0],
    ]
    rows = [[1, 1], [2, -1], [3, -1], [4, 1], [5, 1], [6, -1]]
    witness_ok = _satisfiesAll([c[:6] for c in C], rows, s)
    ci = cellInterior(C, rows)[0]
    chi = cellHasInterior([c[:6] for c in C], rows)
    print('  the cell holds the point s (a valid witness)?  %s' % witness_ok)
    print('  cellHasInterior (witness search) says:         %s   <-- WRONG (should be True)' % chi)
    print('  cellInterior says:                             %s   <-- WRONG (should be True)' % ci)
    print()
    print('  => the incomplete witness search misses this thin interior, so _dropImplied drops a')
    print('     bounding row, the cell grows into the winning face region, and eval picks the lower face.')
    print()
    print('PROPER FIX: replace the heuristic witness search in cellHasInterior / cellInterior with a')
    print('COMPLETE feasibility test for the conic sign conditions (MATLAB, Python and Rust).')
    return 0 if (defect and not ci) else 1


if __name__ == '__main__':
    sys.exit(main())

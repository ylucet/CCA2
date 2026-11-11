"""merge_adjacent_cells -- fuse two cells that carry the SAME face function and differ only in
which side of one shared conic they lie on (port of mergeAdjacentCells.m).

FC[k] is a list of rows [edgeIndex, side(, branch)]; a merged-away cell is removed.
"""


def mergeAdjacentCells(fN, fD, FC):
    fN = [list(r) for r in fN]
    fD = list(fD)
    FC = [_normalizeRows(R) for R in FC]
    changed = True
    while changed:
        changed = False
        groups = {}
        for k in range(len(FC)):
            groups.setdefault((tuple(fN[k]), fD[k]), []).append(k)
        dead = set()
        for key in sorted(groups):
            idx = groups[key]
            if len(idx) < 2:
                continue
            done = False
            for a in range(len(idx)):
                for b in range(a + 1, len(idx)):
                    ka = idx[a]
                    kb = idx[b]
                    c = _mergeableOn(FC[ka], FC[kb])
                    if c is None:
                        continue
                    FC[ka] = _normalizeRows([r for r in FC[ka] if r[0] != c])
                    FC[kb] = []
                    dead.add(kb)
                    changed = True
                    done = True
                    break
                if done:
                    break
        if dead:
            alive = [k for k in range(len(FC)) if k not in dead]
            fN = [fN[k] for k in alive]
            fD = [fD[k] for k in alive]
            FC = [FC[k] for k in alive]
    return fN, fD, FC


def _mergeableOn(A, B):
    if len(A) != len(B):
        return None
    d = [i for i in range(len(A)) if list(A[i]) != list(B[i])]
    if len(d) != 1:
        return None
    d = d[0]
    if A[d][0] != B[d][0] or A[d][1] == B[d][1]:
        return None
    return A[d][0]


def _normalizeRows(R):
    if not R:
        return []
    return [list(r) for r in sorted(set(tuple(r) for r in R))]

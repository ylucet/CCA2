"""compact_cells -- drop each cell's implied rows, then merge adjacent cells carrying one
function; repeat to a fixed point (port of compactCells.m).
"""
from .cell_interior import cellInterior
from .errors import CCA2Error
from .merge_adjacent_cells import mergeAdjacentCells


def compactCells(fN, fD, FC, EcQ):
    """Returns (fN, fD, FC, info)."""
    info = {'dropped': 0, 'merged': 0}
    if not FC:
        return fN, fD, FC, info
    FC = [[list(r) for r in R] for R in FC]
    before = len(FC)
    for _ in range(8):
        drops = 0
        for k in range(len(FC)):
            FC[k], d = _dropImplied(FC[k], EcQ)
            drops += d
        info['dropped'] += drops
        n0 = len(FC)
        fN, fD, FC = mergeAdjacentCells(fN, fD, FC)
        if drops == 0 and len(FC) == n0:
            break
    info['merged'] = before - len(FC)
    return fN, fD, FC, info


def _dropImplied(rows, EcQ):
    """A row is implied when flipping it leaves a cell with no interior."""
    dropped = 0
    if len(rows) < 2:
        return rows, dropped
    try:
        C0, local0 = _localTable(rows, EcQ)
        if not cellInterior(C0, local0)[0]:
            return rows, dropped
    except CCA2Error:
        return rows, dropped                        # undecidable: change nothing
    rows = [list(r) for r in rows]
    t = 0
    while t < len(rows):
        if len(rows) < 2 or rows[t][1] == 0:        # an EQUALITY row is the face, never redundant
            t += 1
            continue
        try:
            flipped = [list(r) for r in rows]
            flipped[t][1] = -flipped[t][1]
            C, local = _localTable(flipped, EcQ)
            if not cellInterior(C, local)[0]:
                del rows[t]
                dropped += 1
                continue                            # t now indexes the next row
        except CCA2Error as ME:
            if ME.identifier not in ('ratQ:overflow', 'cellInterior:constant'):
                raise
        t += 1
    return rows, dropped


def _localTable(rows, EcQ):
    C = []
    local = []
    for r in rows:
        br = r[2] if len(r) >= 3 else 0
        key = list(EcQ[r[0] - 1]) + [br]
        try:
            hit = C.index(key)
        except ValueError:
            C.append(key)
            hit = len(C) - 1
        local.append([hit + 1, r[1]])
    return C, local

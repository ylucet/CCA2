"""The exact CAD fallback of cellHasInterior (review 2026-11-19 finding #2).

_cadInterior is a complete, exact decision for open cells (see its docstring); these tests pin
(a) the two reported mechanisms, (b) agreement with a dense exact rational grid on random cells,
including conic pairs that share a line factor and vertical lines, and (c) that the fallback is
what decides the far cells (mutation: with the fallback off they are declared empty).
"""
import random
from fractions import Fraction as Fr

from cca2 import cell_has_interior as CHI
from cca2 import rat_q


def _strip(vertical):
    c = [[0, 0, 0, 1, 0, -10], [0, 0, 0, 2, 0, -21], [1, 0, 1, -200, 0, 9999]]
    if not vertical:
        c = [[0, 0, 0, r[4], r[3], r[5]] for r in c[:2]] + [[1, 0, 1, 0, -200, 9999]]
    return [rat_q.conic(r) for r in c], [[1, 1], [2, -1], [3, 1]]


def test_far_strips_vertical_and_horizontal():
    for vertical in (True, False):
        E, rows = _strip(vertical)
        assert CHI.cellHasInterior(E, rows)


def test_mutation_fallback_off_is_detected(monkeypatch):
    """With the CAD fallback switched off these cells are declared empty: the tests above depend
    on it (mutation M4 of reviews/codereview-2026-11-19.md)."""
    monkeypatch.setattr(CHI, '_cadInterior', lambda E, r: False)
    E, rows = _strip(True)
    assert not CHI.cellHasInterior(E, rows)


def _prod(l1, l2):
    a, b = l1, l2     # (a0 x + a1 y + a2)(b0 x + b1 y + b2)
    return [a[0] * b[0], a[0] * b[1] + a[1] * b[0], a[1] * b[1], a[0] * b[2] + a[2] * b[0],
            a[1] * b[2] + a[2] * b[1], a[2] * b[2]]


def _rand_conic(rng):
    k = rng.random()
    L = lambda: [rng.randint(-2, 2), rng.randint(-2, 2), rng.randint(-4, 4)]   # noqa: E731
    if k < 0.2:
        return [0, 0, 0] + L()
    if k < 0.4:
        return [0, 0, 0, rng.choice([-1, 1]), 0, rng.randint(-4, 4)]           # vertical line
    if k < 0.75:
        return _prod(L(), L())
    return [rng.randint(-2, 2) for _ in range(6)]


def _grid_has(cs, sg):
    pts = [Fr(i, 4) + Fr(1, 97) for i in range(-48, 49)]
    for x in pts:
        for y in pts:
            if all(s * (c[0] * x * x + c[1] * x * y + c[2] * y * y + c[3] * x + c[4] * y + c[5]) > 0
                   for c, s in zip(cs, sg)):
                return True
    return False


def test_random_cells_against_exact_grid():
    rng = random.Random(11)
    n_true = n_shared = 0
    for it in range(400):
        base = [_rand_conic(rng) for _ in range(rng.randint(2, 4))]
        if rng.random() < 0.5 and len(base) >= 2:        # force a shared line factor
            L = [rng.randint(-2, 2), rng.randint(1, 2), rng.randint(-3, 3)]
            base[0] = _prod(L, [rng.randint(-2, 2), rng.randint(-2, 2), rng.randint(-3, 3)])
            base[1] = _prod(L, [rng.randint(-2, 2), rng.randint(-2, 2), rng.randint(-3, 3)])
            n_shared += 1
        base = [c for c in base if any(t != 0 for t in c[:5])]
        if len(base) < 2:
            continue
        E = [rat_q.conic(c) for c in base]
        sg = [rng.choice([-1, 1]) for _ in E]
        rows = [[i + 1, s] for i, s in enumerate(sg)]
        pt = CHI._cadPoint(E, rows)
        got = pt is not None
        if _grid_has(E, sg):
            assert got, (E, sg)                           # completeness: never a miss
            n_true += 1
        if got:                                           # soundness: the point is a strict witness
            assert all(s * (c[0] * pt[0] ** 2 + c[1] * pt[0] * pt[1] + c[2] * pt[1] ** 2 + c[3] * pt[0]
                            + c[4] * pt[1] + c[5]) > 0 for c, s in zip(E, sg)), (E, sg, pt)
    assert n_true > 20 and n_shared > 100

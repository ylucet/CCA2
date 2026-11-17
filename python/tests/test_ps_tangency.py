"""The perfect-square tangency certificate (psTangency, 2026-11-16, s-symmetric).

A quad-vs-bridge row whose bridge branch exceeds the quad by a perfect square holds only on the
square's zero set -- a curve -- so the face owns no area. The sign of the square root is
region-dependent, so the certificate is symmetric in the two branches P+ / P-:

    upper, side > 0:  P+ OR  P- a NEGATIVE square (<= 0);
    upper, side < 0:  P+ AND P- POSITIVE squares (>= 0);
    lower, side > 0:  P+ AND P- NEGATIVE squares;
    lower, side < 0:  P+ OR  P- a POSITIVE square.

The hand-checked rows below pin the two corners of the rule: the diagonalKink(c = -2) tangency
(which must certify) and the one-sign draft's false positive (which must be refused). The many
fixture rows the certificate newly takes are pinned by the MATLAB trace replay
(test_trace_subdivision.py, trace_biconj.jsonl).
"""
from cca2.alg_row_fails import _psTangency


def _fp(A2, A1, A0):
    """A2, A1, A0 as dicts {(x, y): c} -> the face polynomial {(x, y, z): c}."""
    FP = {}
    for (x, y), c in A2.items():
        FP[(x, y, 2)] = FP.get((x, y, 2), 0) + c
    for (x, y), c in A1.items():
        FP[(x, y, 1)] = FP.get((x, y, 1), 0) + c
    for (x, y), c in A0.items():
        FP[(x, y, 0)] = FP.get((x, y, 0), 0) + c
    return FP


# diagonalKink(c = -2), top row: D = 16 (x - 1)^2 (y - 1)^2, s = 4 (x - 1) (y - 1),
# P+ = -4 (x + y - 1)^2 (a NEGATIVE square), side > 0 -- certifies (the face owns no area).
FIX = _fp(
    A2={(0, 0): 1},
    A1={(2, 0): -2, (1, 0): -4, (0, 2): -2, (0, 0): 4},
    A0={(4, 0): 1, (3, 0): 4, (2, 2): -2, (2, 1): 8, (2, 0): -4,
        (1, 2): 12, (1, 1): -16, (0, 4): 1, (0, 2): -8, (0, 1): 8})
FIX_ROW = {'FP': FIX, 'num': [0, 0, 0, 0, 0, 0, 0, 2, 0, -1], 'den': 1, 'FB': -1, 'side': 1}


def test_upper_side_positive_negative_square_certifies():
    assert _psTangency(FIX_ROW) is True


# (z - (1 - x - y))(z + 1) against q = 0: D = (x + y - 2)^2, P+ = 4 (a positive square),
# P- = 4 (x + y - 1) (NOT a square). The one-sign draft read "P+ a positive square" as enough;
# the s-symmetric rule needs BOTH branches positive squares for upper / side < 0 -- refused.
FAIL_FP = _fp(A2={(0, 0): 1}, A1={(1, 0): 1, (0, 1): 1}, A0={(1, 0): 1, (0, 1): 1, (0, 0): -1})


def _fail_row(side):
    return {'FP': FAIL_FP, 'num': [0] * 10, 'den': 1, 'FB': -1, 'side': side}


def test_upper_side_negative_is_refused_when_only_one_branch_is_a_square():
    assert _psTangency(_fail_row(-1)) is False


def test_upper_side_positive_is_refused_when_not_a_negative_square():
    assert _psTangency(_fail_row(1)) is False

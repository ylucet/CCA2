"""AlgCrv: a row undefined at a BASE POINT of its polynomial is decided by the face's closure
(#108 residual: D-153's domain corner (0,0) was +inf). Mirrors MATLAB AlgCrvTest's
aROWUndefinedAtABASEPointIsDecidedByTheFacesCLOSURE and Rust tests/alg_crv_base_point.rs."""
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from cca2.alg_crv import AlgCrv          # noqa: E402


def test_a_row_undefined_at_a_base_point_is_decided_by_the_closure():
    # the face 0 <= y <= x: x >= 0, y >= 0 and 1 - z >= 0 with x z - y = 0 (z = y/x); at the
    # corner every z is a root of the fibre, so the row says nothing there
    eca = [{'FP': {(1, 0, 1): 1, (0, 1, 0): -1}, 'FB': 1, 'num': [0] * 9 + [1], 'den': 1}]
    g = AlgCrv([], [], [[0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0]], [[0] * 9 + [5]], [1], [],
               [[[1, 1, 0, 0], [2, 1, 0, 0], [1, 1, 0, 1]]], [None], [0], eca)
    v, _ = g.eval([[0.0, 0.0], [0.5, 0.25], [0.25, 0.5], [-1e-3, 0.0]])
    assert v[0] == 5 and v[1] == 5, v                     # the corner is in the closure
    assert math.isinf(v[2]) and math.isinf(v[3]), v       # outside the face: still a hole

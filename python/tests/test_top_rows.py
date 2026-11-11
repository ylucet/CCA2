"""The rows against the TOP root of a face of degree >= 3 in z (topRowPolys / topFaceRow, ported
2026-11-11), and the double-root rule of algBranch. Hand-checkable cases: G = z^3 - 8 has the top
root 2 everywhere."""
from cca2.alg_branch import algBranch
from cca2.alg_row_fails import algRowFails

G = {(0, 0, 3): 1, (0, 0, 0): -8}                       # z^3 - 8: top root 2
BOX = [0.0, 1.0, 0.0, 1.0]


def _yes(_b):
    return True


def _row(q, side):
    return {'FP': G, 'FB': -1, 'num': [0] * 9 + [q], 'den': 1, 'side': side, 'W': None}


def test_a_constant_below_the_top_root_fails_q_above_it():
    assert algRowFails([_row(0, 1)], BOX, _yes, 3) is True      # "0 >= z_top" holds nowhere
    assert algRowFails([_row(0, -1)], BOX, _yes, 3) is False    # "0 <= z_top" holds: no proof


def test_a_constant_above_the_top_root_fails_q_below_it():
    assert algRowFails([_row(3, -1)], BOX, _yes, 3) is True     # "3 <= z_top" holds nowhere
    assert algRowFails([_row(3, 1)], BOX, _yes, 3) is False


def test_a_face_below_a_top_root_is_proved_below():
    lin = {(0, 0, 1): 1, (0, 0, 0): -1}                        # z = 1
    fr = {'FP': lin, 'FB': 1, 'num': None, 'den': None, 'side': None,
          'W': {'FP': G, 'FB': -1}}
    assert algRowFails([fr], BOX, _yes, 3) is True
    above = {(0, 0, 1): 1, (0, 0, 0): -3}                      # z = 3: above the top root
    fr['FP'] = above
    assert algRowFails([fr], BOX, _yes, 3) is False


# the poly2 kink's cubic bridge face (kinkSubdivisionTest); its fibre at (1, 2) is
# (z - 1)^2 (169 z + 1703): the face touches its neighbour there at z = 1
KINK = {(0, 0, 0): 576, (0, 0, 1): 5360, (0, 0, 2): 1885, (0, 0, 3): 169, (0, 1, 0): 3072,
        (0, 1, 1): -1456, (0, 1, 2): -364, (0, 2, 0): -1984, (0, 2, 1): -658, (0, 2, 2): -104,
        (0, 3, 0): -80, (0, 3, 1): -228, (0, 4, 0): -95, (0, 4, 1): 16, (0, 5, 0): 16,
        (1, 0, 0): 192, (1, 0, 1): 8216, (1, 0, 2): 1482, (1, 1, 0): 3200, (1, 1, 1): -3248,
        (1, 1, 2): -312, (1, 2, 0): -1720, (1, 2, 1): -454, (1, 3, 0): 112, (1, 3, 1): 96,
        (1, 4, 0): 136, (2, 0, 0): -5920, (2, 0, 1): 918, (2, 0, 2): -234, (2, 1, 0): 944,
        (2, 1, 1): -1344, (2, 2, 0): 922, (2, 2, 1): 216, (2, 3, 0): 508, (3, 0, 0): -8504,
        (3, 0, 1): -1764, (3, 1, 0): 1488, (3, 1, 1): 216, (3, 2, 0): 990, (4, 0, 0): -2691,
        (4, 0, 1): 81, (4, 1, 0): 972, (5, 0, 0): 378}


def test_a_double_root_split_off_the_real_axis_still_counts():
    # At the box corner (1 - 1 ulp, 2 + 1 ulp) numpy splits the double root 1 into a conjugate pair with
    # imaginary part ~4e-8 and MATLAB's LAPACK into two reals; exactly, the three roots are real.
    # The top branch is 1 either way -- not -10.08, which made the kink's below-test skip the
    # direct comparison MATLAB proves (334 faces against 317).
    z, info = algBranch(KINK, [[0.9999999999999998, 2.000000000000001]], -1,
                         wantInfo=True)
    assert abs(z[0] - 1) < 1e-6 and info['nreal'][0] == 3
    # a genuine complex pair is still dropped: z^3 + z = z (z^2 + 1)
    z, info = algBranch({(0, 0, 3): 1, (0, 0, 1): 1}, [[0.0, 0.0]], -1, wantInfo=True)
    assert z[0] == 0 and info['nreal'][0] == 1

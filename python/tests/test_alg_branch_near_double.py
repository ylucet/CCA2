"""algBranch at NEAR-DOUBLE roots (#108 residual, D-153): mirrors MATLAB algBranchTest's
aDOUBLERootRealOnlyOnALINEIsNotABranch / aDOUBLERootOfTwoCROSSINGBranchesIsKept and Rust
tests/alg_branch_near_double.rs."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

from cca2.alg_branch import algBranch                  # noqa: E402

# p = (z - 5)((z - 10)^2 + y^2): the double root 10 is real ONLY on the line y = 0
ISO = {(0, 0, 3): 1, (0, 0, 2): -25, (0, 0, 1): 200, (0, 2, 1): 1, (0, 0, 0): -500, (0, 2, 0): -5}
# p = (z - 5)((z - 10)^2 - y^2): two real branches 10 +- y CROSSING on y = 0
CRS = {(0, 0, 3): 1, (0, 0, 2): -25, (0, 0, 1): 200, (0, 2, 1): -1, (0, 0, 0): -500, (0, 2, 0): 5}


def test_a_double_root_real_only_on_a_line_is_not_a_branch():
    # f** is closed: its value is the limit from around, and off the line the pair is complex on
    # both sides. Before the fix the top branch was 10 on the line and, through the double-root
    # noise rule, in a ~1e-6 strip next to it.
    for y in (0.0, 1e-9, 1e-7, 1e-6):
        z, info = algBranch(ISO, [[0.3, y]], -1, wantInfo=True)
        assert abs(z[0] - 5) < 1e-9, (y, z, info)
        assert info['nreal'][0] == 1, (y, info)


def test_a_double_root_of_two_crossing_branches_is_kept():
    z, info = algBranch(CRS, [[0.3, 0.0]], -1, wantInfo=True)
    assert abs(z[0] - 10) < 1e-6, z
    assert info['nreal'][0] == 3, info


# the poly2 kink's face at its jump vertex (1,2): 169 z^3 + 1365 z^2 - 3237 z + 1703 =
# (z - 1)^2 (169 z + 1703); f** = 1 is the double root, which rounding returned as 1.0000000256
KINK = {(0, 0, 3): 169, (0, 0, 2): 1365, (0, 0, 1): -3237, (0, 0, 0): 1703}
# (z - 1)^3 (z + 1): a triple root, which rounding returned as 1.0000013
TRIPLE = {(0, 0, 4): 1, (0, 0, 3): -2, (0, 0, 1): 2, (0, 0, 0): -1}


def test_a_genuine_multiple_root_is_returned_to_full_precision():
    z, info = algBranch(KINK, [[0.3, 0.7]], -1, wantInfo=True)
    assert abs(z[0] - 1) <= 4e-16 and info['nreal'][0] == 3, (z, info)
    z, info = algBranch(TRIPLE, [[0.3, 0.7]], -1, wantInfo=True)
    assert abs(z[0] - 1) <= 4e-16 and info['nreal'][0] == 4, (z, info)

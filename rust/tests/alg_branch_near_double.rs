//! alg_branch at NEAR-DOUBLE roots (#108 residual, D-153): mirrors
//! python/tests/test_alg_branch_near_double.py and MATLAB algBranchTest.

use cca2::alg_branch::alg_branch_info;
use cca2::poly_n::PolyN;
use num_bigint::BigInt;

fn poly(t: &[([u32; 3], i64)]) -> PolyN {
    t.iter().map(|(e, c)| (e.to_vec(), BigInt::from(*c))).collect()
}

/// p = (z - 5)((z - 10)^2 + y^2): the double root 10 is real ONLY on the line y = 0, so it is not
/// the value of any branch there (f** is closed: the limit from around). Before the fix the top
/// branch was 10 on the line and, through the double-root noise rule, in a ~1e-6 strip next to it.
#[test]
fn a_double_root_real_only_on_a_line_is_not_a_branch() {
    let p = poly(&[([0, 0, 3], 1), ([0, 0, 2], -25), ([0, 0, 1], 200), ([0, 2, 1], 1), ([0, 0, 0], -500), ([0, 2, 0], -5)]);
    for y in [0.0, 1e-9, 1e-7, 1e-6] {
        let (z, all) = alg_branch_info(&p, &[[0.3, y]], -1).unwrap();
        assert!((z[0] - 5.0).abs() < 1e-9, "y = {y}: {:?} {:?}", z, all);
        assert_eq!(all[0].len(), 1, "y = {y}: {:?}", all);
    }
}

/// The control: p = (z - 5)((z - 10)^2 - y^2), two REAL branches crossing on y = 0: kept.
#[test]
fn a_double_root_of_two_crossing_branches_is_kept() {
    let p = poly(&[([0, 0, 3], 1), ([0, 0, 2], -25), ([0, 0, 1], 200), ([0, 2, 1], -1), ([0, 0, 0], -500), ([0, 2, 0], 5)]);
    let (z, all) = alg_branch_info(&p, &[[0.3, 0.0]], -1).unwrap();
    assert!((z[0] - 10.0).abs() < 1e-6, "{:?}", z);
    assert_eq!(all[0].len(), 3, "{:?}", all);
}

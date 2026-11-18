//! The perfect-square tangency certificate (psTangency, 2026-11-16, s-symmetric).
//!
//! A quad-vs-bridge row whose bridge branch exceeds the quad by a perfect square holds only on the
//! square's zero set -- a curve -- so the face owns no area. The sign of the square root is
//! region-dependent, so the certificate is symmetric in the two branches P+ / P-:
//!
//!     upper, side > 0:  P+ OR  P- a NEGATIVE square (<= 0);
//!     upper, side < 0:  P+ AND P- POSITIVE squares (>= 0);
//!     lower, side > 0:  P+ AND P- NEGATIVE squares;
//!     lower, side < 0:  P+ OR  P- a POSITIVE square.
//!
//! The hand-checked rows below pin the two corners of the rule: the diagonalKink(c = -2) tangency
//! (which must certify) and the one-sign draft's false positive (which must be refused). The many
//! fixture rows the certificate newly takes are pinned by the parity replay (tests/parity.rs,
//! biconj.jsonl #69/#95/#110 now compare the face count too).
//! Mirror of python/tests/test_ps_tangency.py.

use cca2::alg_row_fails::{ps_tangency, RowSpec};
use cca2::poly_n::PolyN;
use cca2::rat_q::Z;
use num_traits::Zero;
use std::collections::BTreeMap;
use std::sync::Arc;

type Biv = BTreeMap<(i64, i64), i64>;

fn terms(v: &[(i64, i64, i64)]) -> Biv {
    v.iter().map(|(x, y, c)| ((*x, *y), *c)).collect()
}

/// A2, A1, A0 as dicts {(x, y): c} -> the face polynomial {(x, y, z): c}.
fn fp(a2: Biv, a1: Biv, a0: Biv) -> PolyN {
    let mut p: PolyN = BTreeMap::new();
    for (k, c) in a2 {
        *p.entry(vec![k.0 as u32, k.1 as u32, 2]).or_insert_with(Z::zero) += c;
    }
    for (k, c) in a1 {
        *p.entry(vec![k.0 as u32, k.1 as u32, 1]).or_insert_with(Z::zero) += c;
    }
    for (k, c) in a0 {
        *p.entry(vec![k.0 as u32, k.1 as u32, 0]).or_insert_with(Z::zero) += c;
    }
    p
}

fn row(fp: PolyN, num: [i64; 10], den: i64, fb: i64, side: i64) -> RowSpec {
    let num: Vec<Z> = num.iter().map(|x| Z::from(*x)).collect();
    RowSpec::Quad { fp: Arc::new(fp), fb, num, den: Z::from(den), side }
}

// diagonalKink(c = -2), top row: D = 16 (x - 1)^2 (y - 1)^2, s = 4 (x - 1) (y - 1),
// P+ = -4 (x + y - 1)^2 (a NEGATIVE square), side > 0 -- certifies (the face owns no area).
#[test]
fn upper_side_positive_negative_square_certifies() {
    let f = fp(
        terms(&[(0, 0, 1)]),
        terms(&[(2, 0, -2), (1, 0, -4), (0, 2, -2), (0, 0, 4)]),
        terms(&[(4, 0, 1), (3, 0, 4), (2, 2, -2), (2, 1, 8), (2, 0, -4), (1, 2, 12), (1, 1, -16), (0, 4, 1), (0, 2, -8), (0, 1, 8)]),
    );
    assert!(ps_tangency(&row(f, [0, 0, 0, 0, 0, 0, 0, 2, 0, -1], 1, -1, 1)));
}

// (z - (1 - x - y))(z + 1) against q = 0: D = (x + y - 2)^2, P+ = 4 (a positive square),
// P- = 4 (x + y - 1) (NOT a square). The one-sign draft read "P+ a positive square" as enough;
// the s-symmetric rule needs BOTH branches positive squares for upper / side < 0 -- refused.
fn fail_fp() -> PolyN {
    fp(terms(&[(0, 0, 1)]), terms(&[(1, 0, 1), (0, 1, 1)]), terms(&[(1, 0, 1), (0, 1, 1), (0, 0, -1)]))
}

#[test]
fn upper_side_negative_is_refused_when_only_one_branch_is_a_square() {
    assert!(!ps_tangency(&row(fail_fp(), [0; 10], 1, -1, -1)));
}

#[test]
fn upper_side_positive_is_refused_when_not_a_negative_square() {
    assert!(!ps_tangency(&row(fail_fp(), [0; 10], 1, -1, 1)));
}

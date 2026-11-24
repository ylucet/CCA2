//! AlgCrv: a row undefined at a BASE POINT of its polynomial is decided by the face's closure
//! (#108 residual: D-153's domain corner (0,0) was +inf). Mirrors
//! python/tests/test_alg_crv_base_point.py and MATLAB AlgCrvTest.

use cca2::alg_con::{AlgCrv, EcA};
use num_bigint::BigInt;
use std::sync::Arc;

fn z(v: &[i64]) -> Vec<BigInt> {
    v.iter().map(|x| BigInt::from(*x)).collect()
}

/// The face 0 <= y <= x: x >= 0, y >= 0 and 1 - z >= 0 with x z - y = 0 (z = y/x). At the corner
/// every z is a root of the fibre: the row is undefined there, not failed.
#[test]
fn a_row_undefined_at_a_base_point_is_decided_by_the_closure() {
    let fp: cca2::poly_n::PolyN = [(vec![1u32, 0, 1], BigInt::from(1)), (vec![0u32, 1, 0], BigInt::from(-1))].into_iter().collect();
    let g = AlgCrv {
        nf: 1,
        ecq: vec![z(&[0, 0, 0, 1, 0, 0]), z(&[0, 0, 0, 0, 1, 0])],
        fn_: vec![z(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 5])],
        fd: vec![BigInt::from(1)],
        fc: vec![vec![[1, 1, 0, 0], [2, 1, 0, 0], [1, 1, 0, 1]]],
        fp: vec![None],
        fb: vec![0],
        eca: vec![EcA { fp: Arc::new(fp), fb: 1, num: z(&[0, 0, 0, 0, 0, 0, 0, 0, 0, 1]), den: BigInt::from(1) }],
        eci: vec![],
        tri: vec![],
        ft: vec![0],
    };
    let (v, _) = g.eval(&[[0.0, 0.0], [0.5, 0.25], [0.25, 0.5], [-1e-3, 0.0]]).unwrap();
    assert!(v[0] == 5.0 && v[1] == 5.0, "the corner is in the closure: {:?}", v);
    assert!(v[2].is_infinite() && v[3].is_infinite(), "outside the face: still a hole: {:?}", v);
}

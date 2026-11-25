//! Two BRIDGE faces never own one cell (DESIGN.md "Two bridge faces never share a cell").
//!
//! A point x in a bridge region has its base point in the subdifferential of the pair's envelope at
//! x, a single point for almost every x; two bridges of one envelope have tie arcs meeting in
//! finitely many points, so their regions meet in a null set. Checked on the diagonal kink (c = -2,
//! two bridges) two independent ways: the defining predicate `bridge_member` on random points, and
//! the subdivision's own count `SubInfo::bridge_overlaps`.
//! Mirror of python/tests/test_two_bridges.py.
use cca2::bridge::bridge_member;
use cca2::conj_q::conj_q;
use cca2::envelope::{envelope_assemble, envelope_prune, Env, FaceKind};
use cca2::envelope_subdivision_alg::esa;
use cca2::alg_con::SubKind;
use cca2::qua_pol::QuaPol;
use cca2::rat_q::{qi, Num};

fn n(v: &[i64]) -> Vec<Num> {
    v.iter().map(|&x| Num::Q(qi(x))).collect()
}

fn tri(v: [[i64; 2]; 3], c: i64) -> QuaPol {
    let mut q = [0, 0, 0, 0, 2, 0, 2, 0, 0, 0];
    q[7] -= c;
    q[8] += c;
    QuaPol::new(v.iter().map(|p| n(p)).collect(), vec![[1, 2, 1], [2, 3, 1], [3, 1, 1]], vec![n(&q)], vec![[1, 0], [1, 0], [1, 0]]).unwrap()
}

fn kink(c: i64) -> Env {
    let (ga, _) = conj_q(&tri([[0, 0], [1, 0], [1, 1]], 0)).unwrap();
    let (gb, _) = conj_q(&tri([[0, 0], [1, 1], [0, 1]], c)).unwrap();
    envelope_prune(&envelope_assemble(&ga, &gb, &[]).unwrap().0).unwrap()
}

#[test]
fn two_bridge_regions_meet_in_no_sampled_point() {
    let env = kink(-2);
    let bridges: Vec<_> = env.bridge.iter().filter(|f| f.kind != FaceKind::Corner && f.pair.is_some()).collect();
    assert_eq!(bridges.len(), 2, "the fixture must have two bridges, or the test is vacuous");
    assert!(bridges.iter().all(|b| b.kind == FaceKind::Algebraic), "membership exactly the sweep");
    // a fixed LCG: no rand dependency, and the same points every run
    let mut s: u64 = 20261125;
    let mut u = || {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let x: Vec<[f64; 2]> = (0..20000).map(|_| [u() * 5.0 - 2.0, u() * 5.0 - 2.0]).collect();
    let m: Vec<Vec<bool>> = bridges.iter().map(|b| bridge_member(b.pair.as_ref().unwrap(), &x).unwrap()).collect();
    for mi in &m {
        assert!(mi.iter().filter(|t| **t).count() > 300, "each bridge must own a sizeable sample");
    }
    let both = (0..x.len()).filter(|&k| m[0][k] && m[1][k]).count();
    assert_eq!(both, 0, "{} sampled points lie in both bridge regions", both);
}

#[test]
fn the_subdivision_puts_no_cell_in_two_bridges() {
    let env = kink(-2);
    let (sub, info) = esa(std::slice::from_ref(&env), &[], &[], &Default::default()).unwrap();
    assert!(sub.iter().any(|f| f.kind == SubKind::Bridge));
    assert_eq!(info.bridge_overlaps, 0);
}

//! Failing tests for the findings of the 2026-11-19 code review (reviews/codereview-2026-11-19.md).
//!
//! Each test asserts the CORRECT behaviour and FAILS against the reviewed code, so each is
//! `#[ignore]`d with the finding number; `cargo test --release --test review_2026_11_19 -- --ignored`
//! runs them. Remove the `#[ignore]` when the finding is fixed.
use cca2::qua_pol::QuaPol;
use cca2::rat_q::{qi, Num, Z};

fn n(v: &[i64]) -> Vec<Num> {
    v.iter().map(|&x| Num::Q(qi(x))).collect()
}

fn pol(v: &[[i64; 2]], e: &[[usize; 3]], f: &[[i64; 10]], ff: &[[usize; 2]]) -> QuaPol {
    QuaPol::new(v.iter().map(|p| n(p)).collect(), e.to_vec(), f.iter().map(|r| n(r)).collect(), ff.to_vec()).unwrap()
}

/// Finding 1 (CRITICAL): f(x,y) = y^2/2 on the cone {x >= |y|}. Along x = (t,0) the objective is
/// s1*t, so f*(s) = +inf for s1 > 0; recession_conditions only tests the two edge rays.
#[test]
#[ignore = "review 2026-11-19 #1: conj finite where f* = +inf (interior null recession direction)"]
fn conj_cone_with_interior_null_direction() {
    let f = pol(&[[0, 0], [1, -1], [1, 1]], &[[1, 2, 0], [1, 3, 0]], &[[0, 0, 0, 0, 0, 0, 1, 0, 0, 0]], &[[1, 0], [0, 1]]);
    let g = cca2::conj_q::conj_q(&f).unwrap().0;
    let (v, _) = g.eval(&[[1.0, 0.0], [0.5, 0.3], [0.1, 0.0], [-1.0, 0.0]]).unwrap();
    assert!(v[0].is_infinite() && v[1].is_infinite() && v[2].is_infinite(), "got {:?}", v);
    assert_eq!(v[3], 0.0);
}

/// Finding 2 (MAJOR): {10 <= x <= 10.5} outside the circle |x - (100,0)| = 1 has interior
/// (e.g. (10.25, 0)); the CAD fallback skips vertical lines and calls it empty.
#[test]
#[ignore = "review 2026-11-19 #2: cell_has_interior misses a far vertical strip"]
fn cell_has_interior_far_vertical_strip() {
    let z = |v: &[i64]| v.iter().map(|&x| Z::from(x)).collect::<Vec<Z>>();
    let ecq = vec![z(&[0, 0, 0, 1, 0, -10]), z(&[0, 0, 0, 2, 0, -21]), z(&[1, 0, 1, -200, 0, 9999])];
    assert!(cca2::cells::cell_has_interior(&ecq, &[(1, 1), (2, -1), (3, 1)]).unwrap());
}

/// Finding 2, second mechanism: degenerate conic rows sharing the line x + y + 3 = 0 (conic_meet
/// reports them degenerate and returns no points, so the other lines' crossings are not CAD
/// abscissae). Reached inside a real conj call (python/tests/review/hunt_cell_interior.py);
/// (-58.4253, -109.3822) satisfies every row strictly.
#[test]
#[ignore = "review 2026-11-19 #2: cell_has_interior misses a cell bounded by line-pair conics"]
fn cell_has_interior_shared_line_factor() {
    let z = |v: &[i64]| v.iter().map(|&x| Z::from(x)).collect::<Vec<Z>>();
    let ecq = vec![
        z(&[101, 76, -25, 3378, 3000, 9225]), z(&[91, 38, -53, -174, -606, -1341]),
        z(&[64, 37, -27, 2386, 2113, 6582]), z(&[1221, 154, -1067, -42386, -49250, -138147]),
        z(&[0, 0, 0, 3, -5, 65]), z(&[0, 0, 0, 3, -5, -31]), z(&[0, 0, 0, -18, 8, -119]),
        z(&[0, 0, 0, 9, -4, 313]), z(&[0, 0, 0, 5, -2, 169]),
    ];
    let rows = [(1, 1), (2, -1), (3, -1), (4, 1), (5, 1), (6, 1), (7, 1), (8, 1), (9, 1)];
    assert!(cca2::cells::cell_has_interior(&ecq, &rows).unwrap());
}

/// Finding 6 (MAJOR, Rust only): biconj panics -- `sort_by(partial_cmp().unwrap_or(Equal))` on a
/// breakpoint list holding NaN in biconj_atoms::tie_arc ("user-provided comparison function does
/// not correctly implement a total order", Rust >= 1.81). Input: random differential case R4-017,
/// a continuous three-piece fan on a convex pentagon. Python and MATLAB return a result.
#[test]
#[ignore = "review 2026-11-19 #6: panic in biconj_atoms::tie_arc on NaN breakpoints"]
fn biconj_does_not_panic_on_r4_017() {
    let f = pol(
        &[[-3, 1], [-2, -2], [-1, 0], [0, 3], [-2, 3]],
        &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1], [4, 5, 1], [5, 1, 1]],
        &[[0, 0, 0, 0, 2, 2, 1, 1, -1, -2], [0, 0, 0, 0, 2, 0, -7, 2, -1, -1], [0, 0, 0, 0, 10, -10, 5, 24, -25, 17]],
        &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 3], [3, 0], [3, 0]],
    );
    let r = std::panic::catch_unwind(|| cca2::biconj_q::biconj_q(&f, &Default::default()).map(|h| h.kind().to_string()));
    assert!(r.is_ok(), "biconj_q panicked");
}

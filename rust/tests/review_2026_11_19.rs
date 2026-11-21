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
/// not correctly implement a total order", Rust >= 1.81). Five random-differential inputs hit it
/// (R4-017, R4-080, R1-091, R4-115, R4-171; Python answers or refuses by name on all five);
/// each must answer or refuse by name, never panic.
#[test]
fn biconj_does_not_panic_on_r4_017() {
    let cases: Vec<(&str, &[[i64; 2]], &[[usize; 3]], &[[i64; 10]], &[[usize; 2]])> = vec![
        ("R4-017",
         &[[-3, 1], [-2, -2], [-1, 0], [0, 3], [-2, 3]],
         &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1], [4, 5, 1], [5, 1, 1]],
         &[[0, 0, 0, 0, 2, 2, 1, 1, -1, -2], [0, 0, 0, 0, 2, 0, -7, 2, -1, -1], [0, 0, 0, 0, 10, -10, 5, 24, -25, 17]],
         &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 3], [3, 0], [3, 0]]),
        ("R4-080",
         &[[-3, -2], [0, -3], [0, -1], [-2, 1], [-3, -1]],
         &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1], [4, 5, 1], [5, 1, 1]],
         &[[0, 0, 0, 0, -1, -2, 1, -1, 0, 1], [0, 0, 0, 0, -3, -1, 13, 0, 12, 7], [0, 0, 0, 0, -9, 6, 9, -4, 25, 14]],
         &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 3], [3, 0], [3, 0]]),
        ("R1-091",
         &[[-3, -3], [0, 0], [0, 2], [-2, 0]],
         &[[1, 2, 1], [2, 3, 1], [3, 4, 1], [4, 1, 1]],
         &[[0, 0, 0, 0, 2, 0, -1, -1, 0, -2]],
         &[[1, 0], [1, 0], [1, 0], [1, 0]]),
        ("R4-115",
         &[[-3, -2], [-1, -3], [3, -1], [3, 1], [-2, 3]],
         &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1], [4, 5, 1], [5, 1, 1]],
         &[[0, 0, 0, 0, 1, 2, 2, -2, 1, 1], [0, 0, 0, 0, -1, 9, -10, 9, -20, -17], [0, 0, 0, 0, 11, 3, -34, 6, -32, -20]],
         &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 3], [3, 0], [3, 0]]),
        ("R4-171",
         &[[-1, 0], [1, -2], [2, -1], [3, 1], [1, 2]],
         &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1], [4, 5, 1], [5, 1, 1]],
         &[[0, 0, 0, 0, 1, 0, 1, 2, -2, 0], [0, 0, 0, 0, 3, 3, 1, 5, 4, 2], [0, 0, 0, 0, 3, 1, 17, 6, -2, 3]],
         &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 3], [3, 0], [3, 0]]),
    ];
    for (name, v, e, f, ff) in cases {
        let f = pol(v, e, f, ff);
        let r = std::panic::catch_unwind(|| cca2::biconj_q::biconj_q(&f, &Default::default()).map(|h| h.kind().to_string()));
        assert!(r.is_ok(), "biconj_q panicked on {}", name);
    }
}

/// Finding 8 (CRITICAL): f = y^2/2 on the L-shape (0,0),(2,0),(2,1),(1,1),(1,2),(0,2). At
/// s = (2.53125, -1.453125) the sup is 2*s1 = 5.0625 (at (2,0)); all three engines return
/// 6.5137939453125 = (2 s1 + s2)^2 / 2, a face used outside its availability range.
#[test]
#[ignore = "review 2026-11-19 #8: conj too high on the L-shape with y^2/2"]
fn conj_l_shape_y2() {
    let v = [[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]];
    let e: Vec<[usize; 3]> = (0..6).map(|i| [i + 1, (i + 1) % 6 + 1, 1]).collect();
    let f = pol(&v, &e, &[[0, 0, 0, 0, 0, 0, 1, 0, 0, 0]], &[[1, 0]; 6]);
    let (val, _) = cca2::conj_q::conj_q(&f).unwrap().0.eval(&[[2.53125, -1.453125]]).unwrap();
    assert!((val[0] - 5.0625).abs() < 1e-12, "got {}", val[0]);
}

/// Finding 9a (CRITICAL): f = (y^2 - x^2)/2 on the unit square; at (-0.49, 0.01), outside the
/// square, f** = +inf. Rust (and Python, MATLAB) return 0.24616241455078125.
#[test]
fn biconj_infinite_outside_the_square() {
    let f = pol(&[[0, 0], [1, 0], [1, 1], [0, 1]], &[[1, 2, 1], [2, 3, 1], [3, 4, 1], [4, 1, 1]],
                &[[0, 0, 0, 0, -1, 0, 1, 0, 0, 0]], &[[1, 0]; 4]);
    let v = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap().eval(&[[-0.4921875, 0.01171875]]).unwrap();
    assert!(v[0].is_infinite(), "got {}", v[0]);
}

/// Finding 9b (CRITICAL): f = -xy + y^2/2 on the L-shape. Fenchel-Young against the engine's OWN
/// conj: f**(x) >= <s,x> - f*(s) = 1.04606 at the s below; the engine's f**(x) is 0.40871.
#[test]
fn biconj_fenchel_young_l_shape() {
    let v = [[0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2]];
    let e: Vec<[usize; 3]> = (0..6).map(|i| [i + 1, (i + 1) % 6 + 1, 1]).collect();
    let f = pol(&v, &e, &[[0, 0, 0, 0, 0, -1, 1, 0, 0, 0]], &[[1, 0]; 6]);
    let x = [0.2578125, 1.76171875];
    let s = [-1.9484626037100141, 1.6789473678974558];
    let g = cca2::conj_q::conj_q(&f).unwrap().0.eval(&[s]).unwrap().0[0];
    let h = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap().eval(&[x]).unwrap()[0];
    assert!(h >= s[0] * x[0] + s[1] * x[1] - g - 1e-9, "f** = {} < <s,x> - f*(s) = {}", h, s[0] * x[0] + s[1] * x[1] - g);
}

/// Finding 9c/9d (CRITICAL): two continuous pieces on a quadrilateral (R3-072). x is in piece 1,
/// f(x) = -12.29, so f**(x) is finite; the certified bracket is [-14.68169, -14.67993]. At
/// (-2.49, 1.85) a convex combination of three domain points gives -18.455, so f** <= that.
/// The tie-arc screen parametrised the parabola tie through a spurious 1e34 point and dropped
/// the needed atom pair (+inf, and too high where another pair covered).
#[test]
fn biconj_finite_inside_the_domain_r3_072() {
    let f = pol(&[[-3, 2], [-1, 0], [2, 0], [-2, 3]], &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1]],
                &[[0, 0, 0, 0, 0, 2, -2, 2, 1, 1], [0, 0, 0, 0, 4, 9, 8, -4, -8, 5]],
                &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 0]]);
    let h = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap();
    let v = h.eval(&[[-2.4921875, 1.51171875], [-2.49, 1.85]]).unwrap();
    assert!(v[0] >= -14.681688881664497 - 1e-9 && v[0] <= -14.679931640625, "got {}", v[0]);
    assert!(v[1] <= -18.455, "got {}", v[1]);
}

/// Finding 9c (CRITICAL), small-cell mechanism: two continuous pieces on a two-cone fan (R6-075).
/// x is in piece 2, so f**(x) is finite; certified bracket [1.80731, 1.86682]. The subdivision
/// dropped the live unit-wide cell holding x: its interior-point grid spans a 1e3 box.
#[test]
fn biconj_finite_inside_the_domain_r6_075() {
    let f = pol(&[[0, 0], [1, 1], [0, 1], [-1, 2]], &[[1, 2, 0], [1, 3, 0], [1, 4, 0]],
                &[[0, 0, 0, 0, 2, 1, 2, 1, 1, 0], [0, 0, 0, 0, 0, -1, 2, 3, 1, 0]],
                &[[1, 0], [2, 1], [0, 2]]);
    let v = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap().eval(&[[-0.4088541666666667, 1.1783854166666667]]).unwrap();
    assert!(v[0] >= 1.8073059721130218 - 1e-9 && v[0] <= 1.8668169186760986, "got {}", v[0]);
}

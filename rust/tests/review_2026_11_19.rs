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

/// Finding 10 (MAJOR) / #108: a discontinuous input is read LOWER SEMICONTINUOUS (at a shared
/// point, the min over the closed pieces holding it). The 30 D fixtures of the review generator
/// must give the values of python/tests/review/d_family_ref.json -- the Python engine's, graded 0
/// wrong by compare_engines.py (conj exact, biconj inside the certified bracket) -- or the same
/// named refusal. Before #108 biconj was above a certified upper bound on 12 of them (D-160: 14.36
/// against -4.89), and #96 refused them all as PLQ:notContinuous.
#[test]
fn discontinuous_input_lsc_values() {
    check_d_family(|name| !D_SLOW.contains(&name));
}

/// #108 residual (fixed 2026-11-24): the AlgCrv took the max over the faces holding at a point, and
/// cutting the contact triangle out of a pair's face left a ZERO-AREA remnant on the shared diagonal
/// carrying the pair's value: 11 at the vertex (-2,2), where lsc f = min(11, -7) = -7, and 9.875 at
/// (-13/8, 11/8) (lsc f = -4.015625). Now the triangle decides where it holds (its plane is f** on
/// the closed triangle).
#[test]
fn discontinuous_input_shared_vertex_d132() {
    let f = pol(&[[-2, 2], [-1, -2], [1, -3], [2, 2]], &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1]],
                &[[0, 0, 0, 0, -1, -1, 2, -2, 1, -1], [0, 0, 0, 0, -1, 0, -1, 0, -2, 1]],
                &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 0]]);
    let v = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap().eval(&[[-2.0, 2.0], [-1.625, 1.375]]).unwrap();
    assert!((v[0] + 7.0).abs() < 1e-9, "f** at the vertex (-2,2): {}", v[0]);
    assert!(v[1] <= -4.015625 + 1e-9, "f** at (-13/8, 11/8) above lsc f = -4.015625: {}", v[1]);
}

/// Same cause on D-111: 8.5 at the jump vertex (-2,3) (lsc f = 3) and above lsc f on the diagonal
/// next to it, through zero-area remnants of the cut contact triangle.
#[test]
fn discontinuous_input_shared_vertex_d111() {
    let f = pol(&[[-2, 3], [1, 1], [3, 2], [0, 3]], &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1]],
                &[[0, 0, 0, 0, 0, -2, -1, -2, -1, 0], [0, 0, 0, 0, -1, -2, -2, -2, 0, -2]],
                &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 0]]);
    let x = [[-2.0, 3.0], [-1.375, 2.875], [-0.75, 2.75], [-0.125, 2.625]];
    let lsc: [f64; 4] = [3.0, -0.5546875, -4.21875, -7.9921875];
    let v = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap().eval(&x).unwrap();
    assert!((v[0] - 3.0).abs() < 1e-9, "f** at the vertex (-2,3): {}", v[0]);
    for i in 0..x.len() {
        assert!(v[i] <= lsc[i] + 1e-7 * lsc[i].abs().max(1.0), "f** at {:?} above lsc f = {}: {}", x[i], lsc[i], v[i]);
    }
}

/// #108 residual on D-153 (slow): f** <= lsc f at its vertices and next to them -- 10.2578 for
/// 9.5177 on the hull edge from the jump vertex (-3,-2) (a double root real only on that edge's line
/// counted as a branch), 7.5000017 for 7.5 at (-3,-2) (a double root ~sqrt(eps) off), +inf at the
/// corner (0,0) (a row with a 0/0 there read as failed).
#[test]
#[ignore = "slow: D-153 (run with -- --ignored)"]
fn discontinuous_input_shared_vertex_d153() {
    let f = pol(&[[-3, -2], [0, 0], [2, 3], [-2, 3]], &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1]],
                &[[0, 0, 0, 0, 1, 1, -1, -1, 1, -2], [0, 0, 0, 0, 2, -1, 1, -2, 2, 2]],
                &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 0]]);
    let x = [[-3.0, -2.0], [-2.875, -1.375], [0.0, 0.0], [2.0, 3.0], [0.0, 3.0], [0.5, 3.0]];
    // lsc f there, and the value of f** where it is known in closed form
    let lsc: [f64; 6] = [7.5, 10.2578125, -2.0, 2.5, 12.5, 10.25];
    let v = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap().eval(&x).unwrap();
    for i in 0..x.len() {
        assert!(v[i] <= lsc[i] + 1e-7 * lsc[i].abs().max(1.0), "f** at {:?} above lsc f = {}: {}", x[i], lsc[i], v[i]);
    }
    assert!((v[0] - 7.5).abs() <= 1e-12 && (v[2] + 2.0).abs() <= 1e-12 && (v[3] - 2.5).abs() <= 1e-12, "{:?}", v);
    // on the hull edges f** is the 1-D envelope of lsc f along them: the tangent from the jump
    // vertex to the convex piece's parabola -- 7.5 + (9 + sqrt 51) t on edge (-3,-2)-(-2,3) at t = 1/8,
    // 12.5 - 5x on y = 3 (tangent at x = 0, a triple root of the face there)
    assert!((v[1] - (7.5 + (9.0 + 51f64.sqrt()) / 8.0)).abs() <= 1e-9, "{:?}", v);
    assert!((v[4] - 12.5).abs() <= 1e-12 && (v[5] - 10.0).abs() <= 1e-9, "{:?}", v);
}

/// #108 residual, the two-convex-pieces route: the poly2 kink at c = -2 jumps across its diagonal;
/// at the jump vertex (1,2) lsc f = min(5, 1) = 1 = f**, a DOUBLE root of the face there, which
/// rounding returned as 1.0000000256.
#[test]
#[ignore = "slow: the poly2 kink (run with -- --ignored)"]
fn discontinuous_kink_jump_vertex() {
    let f = pol(&[[-2, 0], [2, 0], [1, 2], [-1, 2]], &[[1, 2, 1], [2, 3, 1], [3, 4, 1], [4, 1, 1], [1, 3, 1]],
                &[[0, 0, 0, 0, 2, 0, 2, 0, 0, 0], [0, 0, 0, 0, 2, 0, 2, 0, -2, 0]],
                &[[1, 0], [1, 0], [2, 0], [2, 0], [2, 1]]);
    let v = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap()
        .eval(&[[1.0, 2.0], [-1.0, 2.0], [-2.0, 0.0], [2.0, 0.0]]).unwrap();
    assert!((v[0] - 1.0).abs() <= 1e-14, "f** at the jump vertex (1,2): {:?}", v);
    assert!((v[1] - 1.0).abs() <= 1e-12 && (v[2] - 4.0).abs() <= 1e-12 && (v[3] - 4.0).abs() <= 1e-12, "{:?}", v);
}

/// The D fixtures that take minutes (in Python; the slowest in Rust too): `-- --ignored`.
const D_SLOW: [&str; 5] = ["D-069", "D-111", "D-153", "D-188", "D-209"];

#[test]
#[ignore = "slow: the minutes-long D fixtures of #108 (run with -- --ignored)"]
fn discontinuous_input_lsc_values_slow() {
    check_d_family(|name| D_SLOW.contains(&name));
}

fn check_d_family(keep: impl Fn(&str) -> bool) {
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let fixes: Vec<serde_json::Value> = serde_json::from_str(
        &std::fs::read_to_string(root.join("../python/tests/review/d_family_ref.json")).unwrap())
        .unwrap();
    assert_eq!(fixes.len(), 30, "the D family has 30 fixtures");
    let num = |c: &serde_json::Value| -> f64 {
        match c.as_str() {
            Some("inf") => f64::INFINITY,
            Some("-inf") => f64::NEG_INFINITY,
            Some(o) => panic!("unexpected value {}", o),
            None => c.as_f64().unwrap(),
        }
    };
    let pts = |v: &serde_json::Value| -> Vec<[f64; 2]> {
        v.as_array().unwrap().iter().map(|p| [p[0].as_f64().unwrap(), p[1].as_f64().unwrap()]).collect()
    };
    let same = |got: f64, want: f64, tol: f64| -> bool {
        if want.is_infinite() { got == want } else { (got - want).abs() <= tol * want.abs().max(1.0) }
    };
    for x in &fixes {
        let name = x["name"].as_str().unwrap();
        if !keep(name) {
            continue;
        }
        let v: Vec<Vec<Num>> = x["V"].as_array().unwrap().iter().map(|r| {
            r.as_array().unwrap().iter().map(|c| Num::Q(qi(c.as_i64().unwrap()))).collect()
        }).collect();
        let e: Vec<[usize; 3]> = x["E"].as_array().unwrap().iter().map(|r| {
            [r[0].as_u64().unwrap() as usize, r[1].as_u64().unwrap() as usize, r[2].as_u64().unwrap() as usize]
        }).collect();
        let f: Vec<Vec<Num>> = x["f"].as_array().unwrap().iter().map(|r| {
            r.as_array().unwrap().iter().map(|c| Num::Q(qi(c.as_i64().unwrap()))).collect()
        }).collect();
        let ff: Vec<[usize; 2]> = x["F"].as_array().unwrap().iter().map(|r| {
            [r[0].as_u64().unwrap() as usize, r[1].as_u64().unwrap() as usize]
        }).collect();
        let p = QuaPol::new(v, e, f, ff).unwrap();
        assert!(!p.is_continuous().unwrap().0, "{}: the D family jumps across its diagonal", name);
        let (g, _) = cca2::conj_q::conj_q(&p).unwrap();
        let (gv, _) = g.eval(&pts(&x["S"])).unwrap();
        for (k, (a, w)) in gv.iter().zip(x["conj"].as_array().unwrap()).enumerate() {
            assert!(same(*a, num(w), 1e-9), "{}: conj at S[{}] = {} vs {}", name, k, a, num(w));
        }
        let h = cca2::biconj_q::biconj_q(&p, &Default::default());
        if let Some(err) = x.get("biconj_err") {
            let id = h.err().map(|e| e.id).unwrap_or_else(|| "answered".to_string());
            assert_eq!(id, err.as_str().unwrap(), "{}: biconj", name);
            continue;
        }
        let hv = h.unwrap().eval(&pts(&x["X"])).unwrap();
        for (k, (a, w)) in hv.iter().zip(x["biconj"].as_array().unwrap()).enumerate() {
            assert!(same(*a, num(w), 1e-7), "{}: biconj at X[{}] = {} vs {}", name, k, a, num(w));
        }
    }
}

/// Finding 10 guard. R3-072 (a continuous 2-piece split of a quadrilateral: the two quadratics
/// differ by a multiple of the shared edge's own linear form) must PASS the continuity check, so
/// conjQ and biconjQ still answer it.
#[test]
fn continuous_input_not_refused() {
    let f = pol(&[[-3, 2], [-1, 0], [2, 0], [-2, 3]],
                &[[1, 2, 1], [2, 3, 1], [3, 1, 1], [3, 4, 1], [4, 1, 1]],
                &[[0, 0, 0, 0, 0, 2, -2, 2, 1, 1], [0, 0, 0, 0, 4, 9, 8, -4, -8, 5]],
                &[[1, 0], [1, 0], [1, 2], [2, 0], [2, 0]]);
    cca2::conj_q::conj_q(&f).unwrap();
    // review #9c (fixed by #108's realPoint change): f** was +inf at this point of piece 1, where
    // f = -12.29; the value is inside compare_engines.py's certified bracket, same in Python.
    let v = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap().eval(&[[-2.4921875, 1.51171875]]).unwrap();
    assert!((v[0] - (-14.681688881664499)).abs() < 1e-9, "got {}", v[0]);
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

/// Finding 17 (MINOR): ONE identifier for every f** = -infinity refusal. f = -x^2/2 - y^2/2 on
/// the WHOLE PLANE (no vertices, no edges): not convex, so it did not take the convex
/// short-circuit, and f has no affine minorant. It used to refuse as PLQ:biconjQ:noFace while the
/// cone case uses PLQ:biconjQ:minusInfinity for the same answer.
#[test]
fn biconj_whole_plane_refused_as_minus_infinity() {
    let f = QuaPol::new(vec![], vec![], vec![n(&[0, 0, 0, 0, -1, 0, -1, 0, 0, 0])], vec![]).unwrap();
    let e = cca2::biconj_q::biconj_q(&f, &Default::default()).unwrap_err();
    assert_eq!(e.id, "PLQ:biconjQ:minusInfinity", "got {:?}", e.id);
}

/// Finding 17 (MINOR), fan case: f = -x^2 - y^2 on the cone apex (0,0), rays (1,0), (0,1): dom f*
/// is empty (golden record trace_biconj 774). It used to refuse as biconjAtoms:emptyDomain, the
/// third identifier for the same answer.
#[test]
fn biconj_atoms_fan_refused_as_minus_infinity() {
    let f = pol(&[[0, 0], [1, 0], [0, 1]], &[[1, 2, 0], [1, 3, 0]],
                &[[0, 0, 0, 0, -2, 0, -2, 0, 0, 0]], &[[1, 0], [0, 1]]);
    let e = cca2::biconj_atoms::biconj_atoms(&f).unwrap_err();
    assert_eq!(e.id, "PLQ:biconjQ:minusInfinity", "got {:?}", e.id);
}

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
#[ignore = "review 2026-11-19 #9a: biconj finite outside conv(dom f)"]
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
#[ignore = "review 2026-11-19 #9b: biconj below <s,x> - f*(s)"]
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

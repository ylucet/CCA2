//! Star polygons (2n vertices, f = (x^2+y^2)/2): the conjugate must equal the DEFINITION
//! f*(s) = sup_star <s,x> - f(x), computed independently (interior point x = s, else the best point of
//! the boundary, edge by edge) -- not another engine. Regression for the line-vs-conic intersection
//! that dropped roots (DECISIONS.md 2026-11-17: wrong at n = 6..10 before the fix).
//!
//!     cargo test --release --test star                      # n = 6..10, ~15 s
//!     cargo test --release --test star -- --ignored         # n = 11..20, hours (fold only)
use cca2::qua_pol::QuaPol;
use cca2::rat_q::Num;

fn star(n: usize) -> QuaPol {
    let mut v: Vec<Vec<Num>> = vec![];
    for i in 0..2 * n {
        let r = if i % 2 == 0 { 100.0 } else { 45.0 };
        let a = std::f64::consts::PI * i as f64 / n as f64;
        v.push(vec![Num::Q(cca2::rat_q::qi((r * a.cos()).round() as i64)), Num::Q(cca2::rat_q::qi((r * a.sin()).round() as i64))]);
    }
    let m = 2 * n;
    let e: Vec<[usize; 3]> = (0..m).map(|i| [i + 1, (i + 1) % m + 1, 1]).collect();
    let f = vec![[0, 0, 0, 0, 1, 0, 1, 0, 0, 0].iter().map(|&x| Num::Q(cca2::rat_q::qi(x))).collect::<Vec<_>>()];
    let ff: Vec<[usize; 2]> = (0..m).map(|_| [1, 0]).collect();
    QuaPol::new(v, e, f, ff).unwrap()
}

/// The definition: f*(s) = sup over the polygon of <s,x> - |x|^2 (concave in x: the interior point
/// s/2 if it lies in the polygon, else the best point of the boundary, edge by edge).
fn truth(n: usize, s: [f64; 2]) -> f64 {
    let m = 2 * n;
    let vert: Vec<[f64; 2]> = (0..m)
        .map(|i| {
            let r = if i % 2 == 0 { 100.0 } else { 45.0 };
            let a = std::f64::consts::PI * i as f64 / n as f64;
            [(r * a.cos()).round(), (r * a.sin()).round()]
        })
        .collect();
    let h = s;
    let mut inside = false;
    for i in 0..m {
        let (a, b) = (vert[i], vert[(i + 1) % m]);
        if (a[1] > h[1]) != (b[1] > h[1]) && h[0] < (b[0] - a[0]) * (h[1] - a[1]) / (b[1] - a[1]) + a[0] {
            inside = !inside;
        }
    }
    if inside {
        return (s[0] * s[0] + s[1] * s[1]) / 2.0;
    }
    let val = |x: [f64; 2]| s[0] * x[0] + s[1] * x[1] - (x[0] * x[0] + x[1] * x[1]) / 2.0;
    let mut best = f64::NEG_INFINITY;
    for i in 0..m {
        let (a, b) = (vert[i], vert[(i + 1) % m]);
        let d = [b[0] - a[0], b[1] - a[1]];
        // d/dt val(a + t d) = <s - a, d> - t |d|^2
        let t = (((s[0] - a[0]) * d[0] + (s[1] - a[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1])).clamp(0.0, 1.0);
        best = best.max(val([a[0] + t * d[0], a[1] + t * d[1]]));
    }
    best
}


fn check(n: usize) {
    cca2::par::set_tree_reduce(false);
    let g = cca2::conj_q::conj_q(&star(n)).unwrap().0;
    let mut pts = vec![];
    for i in -12..=12 {
        for j in -12..=12 {
            pts.push([i as f64 * 7.3 + 0.11, j as f64 * 5.9 - 0.07]);
        }
    }
    let (v, _) = g.eval(&pts).unwrap();
    let bad: Vec<String> = pts
        .iter()
        .zip(v.iter())
        .filter(|(p, x)| (**x - truth(n, **p)).abs() > 1e-6 * 1f64.max(truth(n, **p).abs()))
        .map(|(p, x)| format!("s={:?} engine {} definition {}", p, x, truth(n, *p)))
        .collect();
    assert!(bad.is_empty(), "n={}: {} of {} points wrong, e.g. {}", n, bad.len(), pts.len(), bad[0]);
}

#[test]
fn stars_6_to_10() {
    for n in 6..=10 {
        check(n);
    }
}

#[test]
#[ignore]
fn stars_11_to_20() {
    for n in 11..=20 {
        check(n);
    }
}

//! conjstress -- the reduction at a size the benchmark fixtures do not reach: conjQ of a star polygon
//! (2n vertices, f = x^2 + y^2) is the maximum of the conjugates of 2n - 2 triangles. Compares the
//! left fold with the balanced tree (same FUNCTION, different face partition) and times both.
//!
//!     conjstress [n ...]        default n = 3 4 5 6   (CCA2_THREADS selects the pool size)
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use cca2::qua_pol::QuaPol;
use cca2::rat_q::Num;
use std::time::Instant;

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

fn timed(p: &QuaPol, tree: bool, reps: usize) -> (f64, cca2::qua_con::QuaCon) {
    cca2::par::set_tree_reduce(tree);
    let mut best = f64::INFINITY;
    let mut g = None;
    for _ in 0..reps {
        cca2::conic_meet::clear_memo();
        cca2::alg_row_fails::clear_memos();
        cca2::rat_q::clear_memos();
        let t0 = Instant::now();
        let r = cca2::conj_q::conj_q(p).unwrap().0;
        best = best.min(t0.elapsed().as_secs_f64());
        g = Some(r);
    }
    (best, g.unwrap())
}

fn main() {
    let ns: Vec<usize> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let ns = if ns.is_empty() { vec![3, 4, 5, 6] } else { ns };
    eprintln!("{}  threads {}", cca2::flint::status(), if cca2::par::enabled() { std::env::var("CCA2_THREADS").unwrap_or("all".into()) } else { "1 (serial)".into() });
    for n in ns {
        let p = star(n);
        let (tf, gf) = timed(&p, false, 3);
        let (tt, gt) = timed(&p, true, 3);
        // same function? sample both on a grid of dual points
        let mut pts = vec![];
        for i in -12..=12 {
            for j in -12..=12 {
                pts.push([i as f64 * 7.3 + 0.11, j as f64 * 5.9 - 0.07]);
            }
        }
        let (vf, _) = gf.eval(&pts).unwrap();
        let (vt, _) = gt.eval(&pts).unwrap();
        let worst = vf.iter().zip(vt.iter()).map(|(a, b)| (a - b).abs() / 1f64.max(a.abs())).fold(0.0, f64::max);
        let bad: Vec<usize> = (0..pts.len()).filter(|&k| (vf[k] - truth(n, pts[k])).abs() > 1e-6 * 1f64.max(truth(n, pts[k]).abs())).collect();
        eprintln!("  fold wrong at {} of {} points; first few:", bad.len(), pts.len());
        for &k in bad.iter().take(4) { eprintln!("    s={:?} fold {} tree {} truth {}", pts[k], vf[k], vt[k], truth(n, pts[k])); }
        let err = |v: &Vec<f64>| pts.iter().zip(v.iter()).map(|(p, x)| (truth(n, *p) - x).abs() / 1f64.max(truth(n, *p).abs())).fold(0.0, f64::max);
        println!("   error vs the definition: fold {:.2e}  tree {:.2e}", err(&vf), err(&vt));
        println!("n={} triangles={} fold {:.4} s ({} faces)  tree {:.4} s ({} faces)  max rel diff of values {:.2e}", n, 2 * n - 2, tf, gf.nf, tt, gt.nf, worst);
    }
}

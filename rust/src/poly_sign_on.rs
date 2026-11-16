//! poly_sign_on -- the sign of a bivariate polynomial on a box, by Bernstein enclosures
//! (port of poly_sign_on.py and its numba kernel _poly_sign_nb.py).

use std::cell::RefCell;
use std::collections::HashMap;
use std::sync::Arc;

pub type FMat = Vec<Vec<f64>>;

struct Tables {
    cb: Vec<Vec<f64>>,
    m1: Vec<Vec<f64>>,
    m2: Vec<Vec<f64>>,
}

fn binom(n: usize) -> Vec<Vec<f64>> {
    let mut c = vec![vec![0.0; n]; n];
    for i in 0..n {
        for k in 0..=i {
            c[i][k] = comb(i, k) as f64;
        }
    }
    c
}

pub fn comb(n: usize, k: usize) -> u128 {
    let mut r: u128 = 1;
    for i in 0..k {
        r = r * (n - i) as u128 / (i + 1) as u128;
    }
    r
}

fn conv_mat(n1: usize) -> Vec<Vec<f64>> {
    let cb = binom(n1);
    let n = n1 - 1;
    let mut m = vec![vec![0.0; n1]; n1];
    for i in 0..n1 {
        for k in 0..=i {
            m[i][k] = cb[i][k] / cb[n][k];
        }
    }
    m
}

thread_local! {
    static TABLES: RefCell<HashMap<(usize, usize), Arc<Tables>>> = RefCell::new(HashMap::new());
}

fn tables(n1: usize, n2: usize) -> Arc<Tables> {
    if let Some(t) = TABLES.with(|m| m.borrow().get(&(n1, n2)).cloned()) {
        return t;
    }
    let t = Arc::new(Tables { cb: binom(n1.max(n2)), m1: conv_mat(n1), m2: conv_mat(n2) });
    TABLES.with(|m| m.borrow_mut().insert((n1, n2), t.clone()));
    t
}

/// Work space of box_sign, kept per thread so that a call allocates nothing. The matrices are
/// row-major, entry (i, j) at i * columns + j.
#[derive(Default)]
struct Scratch {
    pw: Vec<f64>,
    sx: Vec<f64>,
    sy: Vec<f64>,
    t1: Vec<f64>,
    t2: Vec<f64>,
}

thread_local! {
    static SCRATCH: RefCell<Scratch> = RefCell::new(Scratch::default());
}

/// s[k][i] = cb[i][k] * x0^(i-k) * h^k for i >= k, zero below (the Taylor shift to the box).
fn shift(s: &mut Vec<f64>, pw: &mut Vec<f64>, n: usize, x0: f64, h: f64, cb: &[Vec<f64>]) {
    pw.clear();
    pw.extend((0..n).map(|d| x0.powf(d as f64)));
    s.clear();
    s.resize(n * n, 0.0);
    let mut hk = 1.0;
    for k in 0..n {
        for i in k..n {
            s[k * n + i] = cb[i][k] * pw[i - k] * hk;
        }
        hk *= h;
    }
}

/// +1 / -1 when every Bernstein coefficient of P on box b clears the margin, else 0.
pub fn box_sign(p: &FMat, b: &[f64; 4]) -> i32 {
    let n1 = p.len();
    let n2 = p[0].len();
    let t = tables(n1, n2);
    SCRATCH.with(|sc| {
        let sc = &mut *sc.borrow_mut();
        shift(&mut sc.sx, &mut sc.pw, n1, b[0], b[1] - b[0], &t.cb);
        shift(&mut sc.sy, &mut sc.pw, n2, b[2], b[3] - b[2], &t.cb);
        let (sx, sy) = (&sc.sx, &sc.sy);
        let (tt, a) = (&mut sc.t1, &mut sc.t2);
        tt.clear();
        tt.resize(n1 * n2, 0.0);
        a.clear();
        a.resize(n1 * n2, 0.0);
        // tt = sx * p
        for i in 0..n1 {
            for j in 0..n2 {
                let mut acc = 0.0;
                for k in 0..n1 {
                    acc += sx[i * n1 + k] * p[k][j];
                }
                tt[i * n2 + j] = acc;
            }
        }
        // a = tt * sy'
        for i in 0..n1 {
            for j in 0..n2 {
                let mut acc = 0.0;
                for k in 0..n2 {
                    acc += tt[i * n2 + k] * sy[j * n2 + k];
                }
                a[i * n2 + j] = acc;
            }
        }
        // u = m1 * a   (into tt)
        let u = tt;
        for i in 0..n1 {
            let mi = &t.m1[i];
            for j in 0..n2 {
                let mut acc = 0.0;
                for k in 0..n1 {
                    acc += mi[k] * a[k * n2 + j];
                }
                u[i * n2 + j] = acc;
            }
        }
        // bb = u * m2'  (into a)
        let bb = a;
        let mut m = 0.0f64;
        for i in 0..n1 {
            for j in 0..n2 {
                let mj = &t.m2[j];
                let mut acc = 0.0;
                for k in 0..n2 {
                    acc += u[i * n2 + k] * mj[k];
                }
                bb[i * n2 + j] = acc;
                if acc.abs() > m {
                    m = acc.abs();
                }
            }
        }
        let tol = 1e-10 * m.max(2f64.powi(-52));
        let mut pos = true;
        let mut neg = true;
        for &v in bb.iter() {
            if !(v > tol) {
                pos = false;
            }
            if !(v < -tol) {
                neg = false;
            }
        }
        if pos {
            1
        } else if neg {
            -1
        } else {
            0
        }
    })
}

/// polySignOn: Some(+1/-1/0) on the boxes that meet the region, None when none does.
pub fn poly_sign_on(p: &FMat, b: [f64; 4], in_box: &dyn Fn(&[f64; 4]) -> bool, depth: usize) -> Option<i32> {
    if !in_box(&b) {
        return None;
    }
    let s = box_sign(p, &b);
    if s != 0 {
        return Some(s);
    }
    if depth == 0 {
        return Some(0);
    }
    let xm = (b[0] + b[1]) / 2.0;
    let ym = (b[2] + b[3]) / 2.0;
    let q: Vec<i32> = [[b[0], xm, b[2], ym], [xm, b[1], b[2], ym], [b[0], xm, ym, b[3]], [xm, b[1], ym, b[3]]]
        .iter()
        .filter_map(|bb| poly_sign_on(p, *bb, in_box, depth - 1))
        .collect();
    if q.is_empty() {
        return None;
    }
    if q.iter().all(|v| *v == 1) {
        return Some(1);
    }
    if q.iter().all(|v| *v == -1) {
        return Some(-1);
    }
    Some(0)
}

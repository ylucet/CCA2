//! poly_sign_on -- the sign of a bivariate polynomial on a box, by Bernstein enclosures
//! (port of poly_sign_on.py and its numba kernel _poly_sign_nb.py).

use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

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
    static TABLES: RefCell<HashMap<(usize, usize), Rc<Tables>>> = RefCell::new(HashMap::new());
}

fn tables(n1: usize, n2: usize) -> Rc<Tables> {
    if let Some(t) = TABLES.with(|m| m.borrow().get(&(n1, n2)).cloned()) {
        return t;
    }
    let t = Rc::new(Tables { cb: binom(n1.max(n2)), m1: conv_mat(n1), m2: conv_mat(n2) });
    TABLES.with(|m| m.borrow_mut().insert((n1, n2), t.clone()));
    t
}

fn shift(n: usize, x0: f64, h: f64, cb: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut s = vec![vec![0.0; n]; n];
    let mut hk = 1.0;
    for k in 0..n {
        for i in k..n {
            s[k][i] = cb[i][k] * x0.powf((i - k) as f64) * hk;
        }
        hk *= h;
    }
    s
}

/// +1 / -1 when every Bernstein coefficient of P on box b clears the margin, else 0.
pub fn box_sign(p: &FMat, b: &[f64; 4]) -> i32 {
    let n1 = p.len();
    let n2 = p[0].len();
    let t = tables(n1, n2);
    let sx = shift(n1, b[0], b[1] - b[0], &t.cb);
    let sy = shift(n2, b[2], b[3] - b[2], &t.cb);
    let mut tt = vec![vec![0.0; n2]; n1];
    for i in 0..n1 {
        for j in 0..n2 {
            let mut acc = 0.0;
            for k in 0..n1 {
                acc += sx[i][k] * p[k][j];
            }
            tt[i][j] = acc;
        }
    }
    let mut a = vec![vec![0.0; n2]; n1];
    for i in 0..n1 {
        for j in 0..n2 {
            let mut acc = 0.0;
            for k in 0..n2 {
                acc += tt[i][k] * sy[j][k];
            }
            a[i][j] = acc;
        }
    }
    let mut u = vec![vec![0.0; n2]; n1];
    for i in 0..n1 {
        for j in 0..n2 {
            let mut acc = 0.0;
            for k in 0..n1 {
                acc += t.m1[i][k] * a[k][j];
            }
            u[i][j] = acc;
        }
    }
    let mut m = 0.0f64;
    let mut bb = vec![vec![0.0; n2]; n1];
    for i in 0..n1 {
        for j in 0..n2 {
            let mut acc = 0.0;
            for k in 0..n2 {
                acc += u[i][k] * t.m2[j][k];
            }
            bb[i][j] = acc;
            if acc.abs() > m {
                m = acc.abs();
            }
        }
    }
    let tol = 1e-10 * m.max(2f64.powi(-52));
    let mut pos = true;
    let mut neg = true;
    for r in &bb {
        for &v in r {
            if !(v > tol) {
                pos = false;
            }
            if !(v < -tol) {
                neg = false;
            }
        }
    }
    if pos {
        1
    } else if neg {
        -1
    } else {
        0
    }
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

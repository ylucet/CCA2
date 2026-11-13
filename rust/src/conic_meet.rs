//! conic_meet -- the real intersection points of two integer conics (port of conic_meet.py).
//! The elimination is exact (the resultant quartic is built in integers); the POINTS are doubles,
//! returned in canonical (sorted) order, which is what makes "root r of conics i and j" a NAME.

use crate::bail;
use crate::errors::R;
use crate::numeric::{rcond2, roots, EPS};
use crate::rat_q::*;
use num_traits::{Signed, Zero};
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

#[derive(Clone, Debug, Default)]
pub struct MeetInfo {
    pub degenerate: bool,
    pub shear: i64,
    pub resid: Vec<f64>,
    pub sep: f64,
}

pub type Meet = (Vec<[f64; 2]>, MeetInfo);

thread_local! {
    static MEMO: RefCell<HashMap<(Vec<Z>, Vec<Z>), Rc<Meet>>> = RefCell::new(HashMap::new());
}

pub fn clear_memo() {
    MEMO.with(|m| m.borrow_mut().clear());
}

/// Returns (P, info): P a list of [x, y] doubles in canonical (sorted) order. Memoised.
pub fn conic_meet(a: &[Z], b: &[Z]) -> R<Rc<Meet>> {
    let key = (a.to_vec(), b.to_vec());
    if let Some(h) = MEMO.with(|m| m.borrow().get(&key).cloned()) {
        return Ok(h);
    }
    let h = Rc::new(conic_meet_raw(&conic(a)?, &conic(b)?, 1e-9, 8)?);
    MEMO.with(|m| {
        let mut m = m.borrow_mut();
        if m.len() > 50000 {
            m.clear();
        }
        m.insert(key, h.clone());
    });
    Ok(h)
}

fn conic_meet_raw(a: &[Z], b: &[Z], tol: f64, polish: usize) -> R<Meet> {
    let mut info = MeetInfo { sep: f64::INFINITY, ..Default::default() };
    if a == b {
        info.degenerate = true;
        return Ok((vec![], info));
    }
    if all_zero(&a[0..5]) || all_zero(&b[0..5]) {
        bail!("conicMeet:constant", "one argument has no x or y in it, so it is a CONSTANT and not a curve.");
    }
    let mut p: Vec<[f64; 2]>;
    if all_zero(&a[0..3]) || all_zero(&b[0..3]) {
        let (l, c) = if all_zero(&a[0..3]) { (a, b) } else { (b, a) };
        p = line_meets_conic(l, c, tol);
    } else {
        let (pp, deg, k) = general_meet(a, b, tol, polish)?;
        info.shear = k;
        if deg {
            info.degenerate = true;
            return Ok((vec![], info));
        }
        p = pp;
    }
    if !p.is_empty() {
        let mut keyed: Vec<[f64; 2]> = p.iter().map(|q| [round12(q[0]), round12(q[1])]).collect();
        keyed.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        let af: Vec<f64> = a.iter().map(to_f64).collect();
        let bf: Vec<f64> = b.iter().map(to_f64).collect();
        let mut out = Vec::new();
        let mut resid = Vec::new();
        for q in keyed {
            let q = newton(&af, &bf, q, polish, 1e-12);
            resid.push((eval_at(&af, q).abs() / conic_scale(&af, q)).max(eval_at(&bf, q).abs() / conic_scale(&bf, q)));
            out.push(q);
        }
        if out.len() > 1 {
            let mut sep = f64::INFINITY;
            for i in 0..out.len() - 1 {
                for j in i + 1..out.len() {
                    sep = sep.min((out[i][0] - out[j][0]).hypot(out[i][1] - out[j][1]));
                }
            }
            info.sep = sep;
        }
        info.resid = resid;
        p = out;
    }
    Ok((p, info))
}

/// MATLAB round(x, 12): an ORDERING key only.
fn round12(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let y = x.abs() * 1e12;
    if y >= 2f64.powi(52) {
        return x;
    }
    let r = (y + 0.5).floor() / 1e12;
    if x >= 0.0 {
        r
    } else {
        -r
    }
}

fn line_meets_conic(l: &[Z], c: &[Z], tol: f64) -> Vec<[f64; 2]> {
    let (d, e, f) = (&l[3], &l[4], &l[5]);
    let (p0, sc) = if d.abs() >= e.abs() { ([-f, zi(0)], d.clone()) } else { ([zi(0), -f], e.clone()) };
    let u = [-e, d.clone()];
    let q2 = &c[0] * &u[0] * &u[0] + &c[1] * &u[0] * &u[1] + &c[2] * &u[1] * &u[1];
    let q1 = zi(2) * &c[0] * &p0[0] * &u[0]
        + &c[1] * (&p0[0] * &u[1] + &p0[1] * &u[0])
        + zi(2) * &c[2] * &p0[1] * &u[1]
        + &sc * (&c[3] * &u[0] + &c[4] * &u[1]);
    let q0 = &c[0] * &p0[0] * &p0[0] + &c[1] * &p0[0] * &p0[1] + &c[2] * &p0[1] * &p0[1]
        + &sc * (&c[3] * &p0[0] + &c[4] * &p0[1])
        + &sc * &sc * &c[5];
    let t = real_roots_exact(&[q2, q1, q0]);
    let p0f = [q_to_f64(&Q::new(p0[0].clone(), sc.clone())), q_to_f64(&Q::new(p0[1].clone(), sc.clone()))];
    let uf = [to_f64(&u[0]), to_f64(&u[1])];
    let mut p: Vec<[f64; 2]> = t
        .iter()
        .map(|ti| [p0f[0] + ti * uf[0], p0f[1] + ti * uf[1]])
        .filter(|q| q[0].is_finite() && q[1].is_finite())
        .collect();
    if p.len() == 2 && (p[0][0] - p[1][0]).hypot(p[0][1] - p[1][1]) <= tol * 1f64.max(p[0][0].hypot(p[0][1])) {
        p.truncate(1);
    }
    p
}

fn general_meet(a: &[Z], b: &[Z], tol: f64, polish: usize) -> R<(Vec<[f64; 2]>, bool, i64)> {
    let k = shear_for(a, b)?;
    let as_ = shear_conic(a, k);
    let bs = shear_conic(b, k);
    let q = resultant_y(&as_, &bs);
    if all_zero(&q) {
        return Ok((vec![], true, k));
    }
    let xs = real_roots_exact(&q);
    let af: Vec<f64> = as_.iter().map(to_f64).collect();
    let bf: Vec<f64> = bs.iter().map(to_f64).collect();
    let mut p = Vec::new();
    for x in xs {
        let ys = real_roots_f(&[af[2], af[1] * x + af[4], af[0] * x * x + af[3] * x + af[5]]);
        for y in ys {
            let pt = newton(&af, &bf, [x, y], polish, 1e-14);
            if eval_at(&bf, pt).abs() <= tol * conic_scale(&bf, pt) && eval_at(&af, pt).abs() <= tol * conic_scale(&af, pt) {
                p.push(pt);
            }
        }
    }
    if !p.is_empty() {
        let kf = k as f64;
        let p2: Vec<[f64; 2]> = p.iter().map(|q| [q[0] + kf * q[1], q[1]]).collect();
        p = dedupe(&p2, tol);
    }
    Ok((p, false, k))
}

fn shear_for(a: &[Z], b: &[Z]) -> R<i64> {
    for k in 0..5i64 {
        let kz = zi(k);
        if !(&a[0] * &kz * &kz + &a[1] * &kz + &a[2]).is_zero() && !(&b[0] * &kz * &kz + &b[1] * &kz + &b[2]).is_zero() {
            return Ok(k);
        }
    }
    bail!("conicMeet:noShear", "no integer shear in 0..4 gives both conics a y^2 term.")
}

fn shear_conic(c: &[Z], k: i64) -> Vec<Z> {
    let k = zi(k);
    vec![
        c[0].clone(),
        zi(2) * &c[0] * &k + &c[1],
        &c[0] * &k * &k + &c[1] * &k + &c[2],
        c[3].clone(),
        &c[3] * &k + &c[4],
        c[5].clone(),
    ]
}

/// Resultant in y of the two conics, as a quartic in x (5 integer coefficients).
fn resultant_y(a: &[Z], b: &[Z]) -> Vec<Z> {
    let a1 = vec![a[2].clone()];
    let b1 = vec![a[1].clone(), a[4].clone()];
    let c1 = vec![a[0].clone(), a[3].clone(), a[5].clone()];
    let a2 = vec![b[2].clone()];
    let b2 = vec![b[1].clone(), b[4].clone()];
    let c2 = vec![b[0].clone(), b[3].clone(), b[5].clone()];
    let z = vec![zi(0)];
    let m = vec![
        vec![a1.clone(), b1.clone(), c1.clone(), z.clone()],
        vec![z.clone(), a1, b1, c1],
        vec![a2.clone(), b2.clone(), c2.clone(), z.clone()],
        vec![z, a2, b2, c2],
    ];
    let q = trim_poly(&poly_det(&m));
    let mut out = vec![zi(0); 5usize.saturating_sub(q.len())];
    out.extend(q);
    out
}

fn poly_det(m: &[Vec<Vec<Z>>]) -> Vec<Z> {
    let n = m.len();
    if n == 1 {
        return m[0][0].clone();
    }
    let mut d = vec![zi(0)];
    for j in 0..n {
        if all_zero(&m[0][j]) {
            continue;
        }
        let sub: Vec<Vec<Vec<Z>>> = m[1..]
            .iter()
            .map(|row| {
                let mut r = row[..j].to_vec();
                r.extend_from_slice(&row[j + 1..]);
                r
            })
            .collect();
        let mut term = conv(&m[0][j], &poly_det(&sub));
        if j % 2 == 1 {
            term = neg(&term);
        }
        d = poly_add1(&d, &term);
    }
    d
}

fn conv(p: &[Z], q: &[Z]) -> Vec<Z> {
    let mut out = vec![zi(0); p.len() + q.len() - 1];
    for (i, x) in p.iter().enumerate() {
        for (j, y) in q.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

fn trim_poly(p: &[Z]) -> Vec<Z> {
    match p.iter().position(|x| !x.is_zero()) {
        None => vec![zi(0)],
        Some(i) => p[i..].to_vec(),
    }
}

/// Distinct real roots of a real polynomial (highest first), sorted. Doubles.
pub fn real_roots_f(p: &[f64]) -> Vec<f64> {
    let i = match p.iter().position(|x| *x != 0.0) {
        None => return vec![],
        Some(i) => i,
    };
    let p = &p[i..];
    if p.len() <= 1 {
        return vec![];
    }
    let m = p.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    let pf: Vec<f64> = p.iter().map(|x| x / m).collect();
    let rr = roots(&pf);
    let mut re: Vec<f64> = rr.iter().filter(|z| z.im.abs() <= 1e-7 * 1f64.max(z.re.abs())).map(|z| z.re).collect();
    re.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    if re.len() > 1 {
        let mut keep = vec![re[0]];
        for i in 1..re.len() {
            if re[i] - re[i - 1] > 1e-10 * 1f64.max(re[i].abs()) {
                keep.push(re[i]);
            }
        }
        re = keep;
    }
    re
}

fn real_roots_exact(p: &[Z]) -> Vec<f64> {
    let s = squarefree_part(p);
    let m = s.iter().map(|x| x.abs()).max().unwrap_or_else(|| zi(0));
    if m.is_zero() {
        return vec![];
    }
    // float(x)/m in Python divides two floats: convert each, then divide
    let mf = to_f64(&m);
    let pf: Vec<f64> = s.iter().map(|x| to_f64(x) / mf).collect();
    // same as _realRoots(float list) except the scaling already done (m of a scaled list is 1)
    real_roots_f(&pf)
}

fn squarefree_part(p: &[Z]) -> Vec<Z> {
    let p = trim_poly(p);
    if p.len() <= 2 {
        return p;
    }
    let g = poly_gcd_int(&p, &polyder(&p));
    if g.len() <= 1 {
        return p;
    }
    match poly_div_exact1(&p, &g) {
        Some(s) if !s.is_empty() => s,
        _ => p,
    }
}

fn poly_gcd_int(p: &[Z], q: &[Z]) -> Vec<Z> {
    let mut p = primitive_part(&trim_poly(p));
    let mut q = primitive_part(&trim_poly(q));
    while q.len() > 1 || (q.len() == 1 && !q[0].is_zero()) {
        let r = pseudo_rem(&p, &q);
        p = q;
        q = primitive_part(&trim_poly(&r));
        if q.len() == 1 && q[0].is_zero() {
            break;
        }
    }
    primitive_part(&p)
}

fn pseudo_rem(p: &[Z], q: &[Z]) -> Vec<Z> {
    let mut r = trim_poly(p);
    let q = trim_poly(q);
    let lq = q[0].clone();
    while r.len() >= q.len() && r.iter().any(|x| !x.is_zero()) {
        let shift = r.len() - q.len();
        let mut qq = q.clone();
        qq.extend(std::iter::repeat(zi(0)).take(shift));
        let r0 = r[0].clone();
        let nr: Vec<Z> = r.iter().zip(qq.iter()).map(|(a, b)| &lq * a - &r0 * b).collect();
        r = trim_poly(&nr[1..]);
        r = primitive_part(&r);
        if r.len() == 1 && r[0].is_zero() {
            break;
        }
    }
    r
}

fn primitive_part(p: &[Z]) -> Vec<Z> {
    let p = trim_poly(p);
    let g = gcd_all(&p);
    let mut p: Vec<Z> = if g > zi(1) { p.iter().map(|x| x / &g).collect() } else { p };
    if p[0].is_negative() {
        p = neg(&p);
    }
    p
}

/// J \ r for a 2x2 system, by elimination with partial pivoting (as LAPACK does).
pub fn solve2f(j: &[[f64; 2]; 2], r: [f64; 2]) -> [f64; 2] {
    let (mut a, mut b, mut c, mut d) = (j[0][0], j[0][1], j[1][0], j[1][1]);
    let (mut r0, mut r1) = (r[0], r[1]);
    if c.abs() > a.abs() {
        std::mem::swap(&mut a, &mut c);
        std::mem::swap(&mut b, &mut d);
        std::mem::swap(&mut r0, &mut r1);
    }
    if a == 0.0 {
        return [f64::NAN, f64::NAN];
    }
    let m = c / a;
    let d2 = d - m * b;
    let y2 = r1 - m * r0;
    if d2 == 0.0 {
        return [f64::NAN, f64::NAN];
    }
    let x2 = y2 / d2;
    let x1 = (r0 - b * x2) / a;
    [x1, x2]
}

fn newton(a: &[f64], b: &[f64], mut p: [f64; 2], steps: usize, rc: f64) -> [f64; 2] {
    for _ in 0..steps {
        let j = [grad(a, p), grad(b, p)];
        let r = [eval_at(a, p), eval_at(b, p)];
        if !j.iter().flatten().all(|v| v.is_finite()) || rcond2(&j) < rc {
            break;
        }
        let st = solve2f(&j, r);
        if !st.iter().all(|v| v.is_finite()) {
            break;
        }
        p = [p[0] - st[0], p[1] - st[1]];
        if st[0].hypot(st[1]) <= EPS * 1f64.max(p[0].hypot(p[1])) {
            break;
        }
    }
    p
}

#[inline]
pub fn eval_at(c: &[f64], p: [f64; 2]) -> f64 {
    c[0] * p[0].powi(2) + c[1] * p[0] * p[1] + c[2] * p[1].powi(2) + c[3] * p[0] + c[4] * p[1] + c[5]
}

#[inline]
fn conic_scale(c: &[f64], p: [f64; 2]) -> f64 {
    let cm = c.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    1f64.max(cm * 1f64.max(p[0].abs().max(p[1].abs())).powi(2))
}

fn grad(c: &[f64], p: [f64; 2]) -> [f64; 2] {
    [2.0 * c[0] * p[0] + c[1] * p[1] + c[3], c[1] * p[0] + 2.0 * c[2] * p[1] + c[4]]
}

fn dedupe(p: &[[f64; 2]], tol: f64) -> Vec<[f64; 2]> {
    let mut keep = vec![true; p.len()];
    for i in 1..p.len() {
        for j in 0..i {
            if keep[j] && (p[i][0] - p[j][0]).hypot(p[i][1] - p[j][1]) <= tol * 1f64.max(p[j][0].hypot(p[j][1])) {
                keep[i] = false;
                break;
            }
        }
    }
    p.iter().zip(keep).filter(|(_, k)| *k).map(|(q, _)| *q).collect()
}

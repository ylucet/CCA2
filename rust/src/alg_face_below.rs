//! alg_face_below -- is the algebraic face `lo` below one of the faces `hi` everywhere in the
//! part of a box inside a region? (port of alg_face_below.py). True is a proof.

use crate::alg_branch::alg_branch;
use crate::alg_con::AlgRow;
use crate::alg_row_fails::{alg_row_fails_box, Region, RowCache, RowSpec};
use crate::errors::R;
use crate::numeric::{linspace, round_half_away};
use crate::poly_n::PolyN;
use crate::rat_q::*;
use std::sync::Arc;

#[derive(Clone)]
pub struct BelowFace {
    pub fp: Arc<PolyN>,
    pub fb: i64,
}

struct Ctx<'a> {
    reg: &'a Region,
    own: Vec<RowSpec>,
    /// the per-cell cache of processed rows (the caller's), or None
    rows: Option<&'a RowCache>,
    /// the depth the recursion started at: the two top levels fork their sub-boxes
    top: usize,
}

/// `rows`: the caller's per-cell RowCache (drop_below), or None.
pub fn alg_face_below(lo: &BelowFace, hi: &[BelowFace], reg: &Region, depth: usize, own: &[AlgRow], rows: Option<&RowCache>) -> R<bool> {
    let own: Vec<RowSpec> = own.iter().map(|o| RowSpec::Quad { fp: o.fp.clone(), fb: o.fb, num: o.num.clone(), den: o.den.clone(), side: o.side }).collect();
    let ctx = Ctx { reg, own, rows, top: depth };
    rec(lo, hi, reg.bx, &ctx, depth)
}

fn rec(lo: &BelowFace, hi: &[BelowFace], b: [f64; 4], ctx: &Ctx, depth: usize) -> R<bool> {
    let in_box = |bb: &[f64; 4]| ctx.reg.in_box(bb);
    if !in_box(&b) {
        return Ok(true);
    }
    let us = linspace(b[0], b[1], 5);
    let vs = linspace(b[2], b[3], 5);
    let mut x: Vec<[f64; 2]> = vec![];
    for u in &us {
        for v in &vs {
            x.push([*u, *v]);
        }
    }
    let zl = alg_branch(&lo.fp, &x, lo.fb)?;
    let n = hi.len();
    let mut zh_all = vec![];
    let mut cnt = vec![];
    for h in hi {
        let zh = alg_branch(&h.fp, &x, h.fb)?;
        cnt.push(zl.iter().zip(zh.iter()).filter(|(a, c)| a.is_finite() && c.is_finite() && a < c).count());
        zh_all.push(zh);
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&j| std::cmp::Reverse(cnt[j]));
    let inner = if depth == 0 { 6 } else { 3 };
    let corners = &ctx.reg.corners;
    let factors = &ctx.reg.factors;
    for &j in &order {
        let zh = &zh_all[j];
        let ok: Vec<bool> = zl.iter().zip(zh.iter()).map(|(a, c)| a.is_finite() && c.is_finite() && a < c).collect();
        if ok.iter().filter(|o| **o).count() < 6 {
            break;
        }
        let xs: Vec<[f64; 2]> = x.iter().zip(ok.iter()).filter(|(_, o)| **o).map(|(p, _)| *p).collect();
        let ls: Vec<f64> = zl.iter().zip(ok.iter()).filter(|(_, o)| **o).map(|(p, _)| *p).collect();
        let hs: Vec<f64> = zh.iter().zip(ok.iter()).filter(|(_, o)| **o).map(|(p, _)| *p).collect();
        let num = match separator(&xs, &ls, &hs) {
            None => continue,
            Some(n) => n,
        };
        let rl = RowSpec::Quad { fp: lo.fp.clone(), fb: lo.fb, num: num[0..10].to_vec(), den: num[10].clone(), side: -1 };
        let rh = RowSpec::Quad { fp: hi[j].fp.clone(), fb: hi[j].fb, num: num[0..10].to_vec(), den: num[10].clone(), side: 1 };
        let mut al = vec![rl];
        al.extend(ctx.own.iter().cloned());
        let mut ah = vec![rh];
        ah.extend(ctx.own.iter().cloned());
        if alg_row_fails_box(&al, b, &in_box, inner, factors, corners, ctx.rows)? && alg_row_fails_box(&ah, b, &in_box, inner, factors, corners, ctx.rows)? {
            return Ok(true);
        }
    }
    let mut in_piece: Vec<bool> = x.iter().map(|p| in_box(&[p[0], p[0], p[1], p[1]])).collect();
    for o in &ctx.own {
        if !in_piece.iter().any(|v| *v) {
            break;
        }
        if let RowSpec::Quad { fp, fb, num, den, side } = o {
            let zo = alg_branch(fp, &x, *fb)?;
            for (i, p) in x.iter().enumerate() {
                let q = eval_face(num, den, p[0], p[1]);
                in_piece[i] = in_piece[i] && zo[i].is_finite() && (*side as f64) * (q - zo[i]) >= -1e-9;
            }
        }
    }
    for &j in &order {
        let zh = &zh_all[j];
        let both: Vec<bool> = zl.iter().zip(zh.iter()).map(|(a, c)| a.is_finite() && c.is_finite()).collect();
        let bad = (0..x.len()).any(|i| both[i] && in_piece[i] && zl[i] > zh[i] + 1e-9 * 1f64.max(zh[i].abs()));
        if both.iter().filter(|o| **o).count() < 6 || bad {
            continue;
        }
        let fr = RowSpec::Face { fp: lo.fp.clone(), fb: lo.fb, wfp: hi[j].fp.clone(), wfb: hi[j].fb };
        let mut al = vec![fr];
        al.extend(ctx.own.iter().cloned());
        if alg_row_fails_box(&al, b, &in_box, inner, factors, corners, ctx.rows)? {
            return Ok(true);
        }
    }
    if depth == 0 {
        return Ok(false);
    }
    let xm = (b[0] + b[1]) / 2.0;
    let ym = (b[2] + b[3]) / 2.0;
    let sub = [[b[0], xm, b[2], ym], [xm, b[1], b[2], ym], [b[0], xm, ym, b[3]], [xm, b[1], ym, b[3]]];
    crate::par::all_of(4, depth + 2 > ctx.top, |i| rec(lo, hi, sub[i], ctx, depth - 1))
}

/// A quadratic with DYADIC coefficients strictly between zl and zh at every sample.
fn separator(x: &[[f64; 2]], zl: &[f64], zh: &[f64]) -> Option<Vec<Z>> {
    let v: Vec<[f64; 6]> = x.iter().map(|p| [p[0] * p[0], p[0] * p[1], p[1] * p[1], p[0], p[1], 1.0]).collect();
    let rhs: Vec<f64> = zl.iter().zip(zh.iter()).map(|(a, c)| (a + c) / 2.0).collect();
    let c = lstsq(&v, &rhs);
    for k in 0..41 {
        let s = 2f64.powi(k);
        let r: Vec<i64> = c.iter().map(|v| round_half_away(v * s)).collect();
        if r.iter().any(|v| (v.abs() as f64) >= 2f64.powi(50)) {
            return None;
        }
        let rf: Vec<f64> = r.iter().map(|x| *x as f64).collect();
        let ok = v.iter().enumerate().all(|(i, row)| {
            let q = (0..6).map(|t| row[t] * rf[t]).sum::<f64>() / s;
            q > zl[i] && q < zh[i]
        });
        if ok {
            return Some(vec![zi(0), zi(0), zi(0), zi(0), zi(2 * r[0]), zi(r[1]), zi(2 * r[2]), zi(r[3]), zi(r[4]), zi(r[5]), Z::from(1i64 << k)]);
        }
    }
    None
}

/// Minimum-norm least squares by one-sided Jacobi SVD (numpy.linalg.lstsq, rcond = eps*max(m,n)).
// QR by Givens rotations: the column index addresses rows p and q at once, so it is load-bearing
#[allow(clippy::needless_range_loop)]
fn lstsq(a: &[[f64; 6]], b: &[f64]) -> Vec<f64> {
    let m = a.len();
    let n = 6;
    let mut u: Vec<Vec<f64>> = (0..n).map(|j| (0..m).map(|i| a[i][j]).collect()).collect(); // columns
    let mut v: Vec<Vec<f64>> = (0..n).map(|j| (0..n).map(|i| if i == j { 1.0 } else { 0.0 }).collect()).collect();
    for _sweep in 0..60 {
        let mut off = 0.0;
        for p in 0..n {
            for q in p + 1..n {
                let alpha: f64 = u[p].iter().map(|x| x * x).sum();
                let beta: f64 = u[q].iter().map(|x| x * x).sum();
                let gamma: f64 = u[p].iter().zip(u[q].iter()).map(|(x, y)| x * y).sum();
                if gamma.abs() <= 1e-300 || gamma.abs() <= f64::EPSILON * (alpha * beta).sqrt() {
                    continue;
                }
                off += gamma.abs();
                let zeta = (beta - alpha) / (2.0 * gamma);
                let t = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
                let t = if zeta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for i in 0..m {
                    let (x, y) = (u[p][i], u[q][i]);
                    u[p][i] = c * x - s * y;
                    u[q][i] = s * x + c * y;
                }
                for i in 0..n {
                    let (x, y) = (v[p][i], v[q][i]);
                    v[p][i] = c * x - s * y;
                    v[q][i] = s * x + c * y;
                }
            }
        }
        if off == 0.0 {
            break;
        }
    }
    let sig: Vec<f64> = (0..n).map(|j| u[j].iter().map(|x| x * x).sum::<f64>().sqrt()).collect();
    let smax = sig.iter().cloned().fold(0.0, f64::max);
    let tol = f64::EPSILON * (m.max(n) as f64) * smax;
    let mut x = vec![0.0; n];
    for j in 0..n {
        if sig[j] <= tol || sig[j] == 0.0 {
            continue;
        }
        let ub: f64 = u[j].iter().zip(b.iter()).map(|(p, q)| p * q).sum::<f64>() / (sig[j] * sig[j]);
        for i in 0..n {
            x[i] += v[j][i] * ub;
        }
    }
    x
}

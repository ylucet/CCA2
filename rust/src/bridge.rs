//! Bridge faces (ports of bridge_phi.py, bridge_witness.py, bridge_arc_cut.py, bridge_member.py,
//! bridge_face_alg.py, bridge_face_rat.py).

use crate::bail;
use crate::envelope::{AlgPair, BFace, FaceKind, PairData, RatPair};
use crate::errors::R;
use crate::numeric::{convf, linspace, polyval_f, roots, round_half_away, EPS};
use crate::poly_n::{self as pn, PolyN};
use crate::rat_q::*;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeMap;
use std::sync::Arc;

// ==========================================================================================
// bridge_phi (ascending powers of mu)
// ==========================================================================================
fn conv(a: &[Z], b: &[Z]) -> Vec<Z> {
    let mut out = vec![Z::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

fn vadd(a: &[Z], b: &[Z]) -> Vec<Z> {
    let n = a.len().max(b.len());
    (0..n).map(|i| a.get(i).cloned().unwrap_or_else(Z::zero) + b.get(i).cloned().unwrap_or_else(Z::zero)).collect()
}

fn sc(a: &[Z], k: &Z) -> Vec<Z> {
    a.iter().map(|x| x * k).collect()
}

pub struct Phi {
    pub phi: Vec<Z>,
    pub dt: Vec<Z>,
    pub u1: Vec<Z>,
    pub u2: Vec<Z>,
}

pub fn bridge_phi(pr: &AlgPair, xn: &[Z], xd: &Z) -> Phi {
    let (n1, m1, n2, m2, dl, dlv, kp, d) = (&pr.n1, &pr.m1, &pr.n2, &pr.m2, &pr.dlm, &pr.dl, &pr.kp, &pr.d);
    let e11 = vec![n1[0][0].clone(), &n1[0][0] - &n2[0][0]];
    let e12 = vec![n1[0][1].clone(), &n1[0][1] - &n2[0][1]];
    let e22 = vec![n1[1][1].clone(), &n1[1][1] - &n2[1][1]];
    let dt = vadd(&conv(&e11, &e22), &neg(&conv(&e12, &e12)));
    let a11 = e22;
    let a12 = neg(&e12);
    let a22 = e11;
    let w1 = vec![d * &xn[0] - xd * &m1[0], -xd * (&m1[0] - &m2[0])];
    let w2 = vec![d * &xn[1] - xd * &m1[1], -xd * (&m1[1] - &m2[1])];
    let u1 = vadd(&conv(&a11, &w1), &conv(&a12, &w2));
    let u2 = vadd(&conv(&a12, &w1), &conv(&a22, &w2));
    let mut phi = vadd(&vadd(&sc(&conv(&u1, &u1), &dl[0][0]), &sc(&conv(&u1, &u2), &(zi(2) * &dl[0][1]))), &sc(&conv(&u2, &u2), &dl[1][1]));
    phi = vadd(&phi, &sc(&conv(&dt, &vadd(&sc(&u1, &dlv[0]), &sc(&u2, &dlv[1]))), &(zi(2) * xd)));
    phi = vadd(&phi, &sc(&conv(&dt, &dt), &(zi(2) * kp * xd * xd)));
    Phi { phi, dt, u1, u2 }
}

fn vaddf(a: &[f64], b: &[f64]) -> Vec<f64> {
    let n = a.len().max(b.len());
    (0..n).map(|i| a.get(i).cloned().unwrap_or(0.0) + b.get(i).cloned().unwrap_or(0.0)).collect()
}
fn scf(a: &[f64], k: f64) -> Vec<f64> {
    a.iter().map(|x| x * k).collect()
}

/// The pair's integer data as doubles (bridgeMember calls this per point).
#[derive(Debug)]
#[allow(dead_code)]
pub struct AlgPairF {
    n1: [[f64; 2]; 2],
    m1: [f64; 2],
    n2: [[f64; 2]; 2],
    m2: [f64; 2],
    dl: [[f64; 2]; 2],
    dlv: [f64; 2],
    kp: f64,
    d: f64,
    // the integer-exact pieces of the family that do not depend on x
    dt: Vec<f64>,
    a11: Vec<f64>,
    a12: Vec<f64>,
    a22: Vec<f64>,
}

impl AlgPairF {
    pub fn new(pr: &AlgPair) -> AlgPairF {
        let f = |m: &Vec<Vec<Z>>| [[to_f64(&m[0][0]), to_f64(&m[0][1])], [to_f64(&m[1][0]), to_f64(&m[1][1])]];
        let ph = bridge_phi(pr, &[zi(0), zi(0)], &zi(1));
        let e11 = vec![pr.n1[0][0].clone(), &pr.n1[0][0] - &pr.n2[0][0]];
        let e12 = vec![pr.n1[0][1].clone(), &pr.n1[0][1] - &pr.n2[0][1]];
        let e22 = vec![pr.n1[1][1].clone(), &pr.n1[1][1] - &pr.n2[1][1]];
        AlgPairF {
            n1: f(&pr.n1),
            m1: [to_f64(&pr.m1[0]), to_f64(&pr.m1[1])],
            n2: f(&pr.n2),
            m2: [to_f64(&pr.m2[0]), to_f64(&pr.m2[1])],
            dl: f(&pr.dlm),
            dlv: [to_f64(&pr.dl[0]), to_f64(&pr.dl[1])],
            kp: to_f64(&pr.kp),
            d: to_f64(&pr.d),
            dt: ph.dt.iter().map(to_f64).collect(),
            a11: e22.iter().map(to_f64).collect(),
            a12: e12.iter().map(|x| -to_f64(x)).collect(),
            a22: e11.iter().map(to_f64).collect(),
        }
    }

    /// bridgePhi at the double point x (xd = 1), in doubles.
    pub fn phi(&self, x: [f64; 2]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        let w1 = vec![self.d * x[0] - self.m1[0], -(self.m1[0] - self.m2[0])];
        let w2 = vec![self.d * x[1] - self.m1[1], -(self.m1[1] - self.m2[1])];
        let u1 = vaddf(&convf(&self.a11, &w1), &convf(&self.a12, &w2));
        let u2 = vaddf(&convf(&self.a12, &w1), &convf(&self.a22, &w2));
        let dt = &self.dt;
        let mut phi = vaddf(&vaddf(&scf(&convf(&u1, &u1), self.dl[0][0]), &scf(&convf(&u1, &u2), 2.0 * self.dl[0][1])), &scf(&convf(&u2, &u2), self.dl[1][1]));
        phi = vaddf(&phi, &scf(&convf(dt, &vaddf(&scf(&u1, self.dlv[0]), &scf(&u2, self.dlv[1]))), 2.0));
        phi = vaddf(&phi, &scf(&convf(dt, dt), 2.0 * self.kp));
        (phi, dt.clone(), u1, u2)
    }
}

// ==========================================================================================
// bridge_member
// ==========================================================================================
fn cabs(z: crate::numeric::Cx) -> f64 {
    z.re.hypot(z.im)
}

fn pairf_of(a: &AlgPair) -> (Arc<AlgPairF>, Arc<AlgPairF>) {
    a.fcache
        .get_or_init(|| {
            let f = Arc::new(AlgPairF::new(a));
            // the re-based family near mu = -1 (nearArtefact)
            let nn2: Vec<Vec<Z>> = a.n2.iter().zip(a.n1.iter()).map(|(ra, rb)| ra.iter().zip(rb.iter()).map(|(x, y)| zi(2) * x - y).collect()).collect();
            let mm2: Vec<Z> = a.m2.iter().zip(a.m1.iter()).map(|(x, y)| zi(2) * x - y).collect();
            let b = AlgPair { d: a.d.clone(), n1: a.n2.clone(), m1: a.m2.clone(), n2: nn2, m2: mm2, dlm: a.dlm.clone(), dl: a.dl.clone(), kp: a.kp.clone(), p: a.p.clone(), fcache: Default::default() };
            (f, Arc::new(AlgPairF::new(&b)))
        })
        .clone()
}

/// Is the (double) point in the bridge face's region? One verdict per point.
pub fn bridge_member(pr: &Arc<PairData>, x: &[[f64; 2]]) -> R<Vec<bool>> {
    match &**pr {
        PairData::Rat(p) => {
            let bf: Vec<f64> = p.b.iter().map(to_f64).collect();
            let cf: Vec<f64> = p.c.iter().map(to_f64).collect();
            let scale = to_f64(&p.scale);
            let (rq, rl, rc, rd) = (to_f64(&p.rq), to_f64(&p.rl), to_f64(&p.rc), to_f64(&p.rd));
            let (u, w, nu) = ([to_f64(&p.u[0]), to_f64(&p.u[1])], [to_f64(&p.w[0]), to_f64(&p.w[1])], to_f64(&p.nu));
            let pf: Vec<[f64; 3]> = p.p.iter().map(|r| [to_f64(&r[0]), to_f64(&r[1]), to_f64(&r[2])]).collect();
            let mut tf = vec![];
            for xx in x {
                let bv = xx[0] * bf[0] + xx[1] * bf[1] + bf[2];
                let cv = xx[0] * cf[0] + xx[1] * cf[1] + cf[2];
                let mut ok = cv * scale < 0.0;
                if ok && !pf.is_empty() {
                    let t = -bv / (2.0 * cv);
                    let r = (rq * t * t + rl * t + rc) / rd;
                    let s = [(t * u[0] + r * w[0]) / nu, (t * u[1] + r * w[1]) / nu];
                    let tol = EPS.sqrt() * 1f64.max(s[0].abs()).max(s[1].abs());
                    ok = pf.iter().all(|q| s[0] * q[0] + s[1] * q[1] + q[2] >= -tol);
                }
                tf.push(ok);
            }
            Ok(tf)
        }
        PairData::Alg(a) => {
            let (f, g) = pairf_of(a);
            let pf: Vec<[f64; 3]> = a.p.iter().map(|r| [to_f64(&r[0]), to_f64(&r[1]), to_f64(&r[2])]).collect();
            let mut tf = vec![false; x.len()];
            for (i, xx) in x.iter().enumerate() {
                let (phi, dt, u1, u2) = f.phi(*xx);
                let p: Vec<f64> = phi.iter().rev().cloned().collect();
                let mx = p.iter().fold(0.0f64, |m, v| m.max(v.abs()));
                let nz = match p.iter().position(|v| v.abs() > EPS * mx) {
                    None => continue,
                    Some(k) => k,
                };
                let p = &p[nz..];
                if p.len() < 2 {
                    continue;
                }
                let r: Vec<f64> = roots(p).into_iter().filter(|z| z.im.abs() <= 1e-9 * 1f64.max(cabs(*z))).map(|z| z.re).filter(|v| *v >= -1.0 - 1e-12 && *v <= 1e-12).collect();
                let dscale = dt.iter().fold(0.0f64, |m, v| m.max(v.abs()));
                let dtr: Vec<f64> = dt.iter().rev().cloned().collect();
                let u1r: Vec<f64> = u1.iter().rev().cloned().collect();
                let u2r: Vec<f64> = u2.iter().rev().cloned().collect();
                for &rj in &r {
                    let mu = rj.max(-1.0).min(0.0);
                    let d = polyval_f(&dtr, mu);
                    if d.abs() <= 1e-9 * 1f64.max(dscale) {
                        continue;
                    }
                    let s = [polyval_f(&u1r, mu) / d, polyval_f(&u2r, mu) / d];
                    if !pf.is_empty() {
                        let lim = -1e-9 * 1f64.max(s[0].abs()).max(s[1].abs());
                        if pf.iter().any(|q| q[0] * s[0] + q[1] * s[1] + q[2] < lim) {
                            continue;
                        }
                    }
                    tf[i] = true;
                    break;
                }
                if tf[i] || !r.iter().any(|v| (v + 1.0).abs() <= 1e-6) || polyval_f(&dtr, -1.0) != 0.0 {
                    continue;
                }
                if near_artefact(&g, &pf, *xx) {
                    tf[i] = true;
                }
            }
            Ok(tf)
        }
    }
}

fn near_artefact(g: &AlgPairF, pf: &[[f64; 3]], x: [f64; 2]) -> bool {
    let (phi, dt, u1, u2) = g.phi(x);
    let p: Vec<f64> = phi.iter().rev().cloned().collect();
    let nz = match p.iter().position(|v| *v != 0.0) {
        None => return false,
        Some(k) => k,
    };
    let p = &p[nz..];
    if p.len() < 2 {
        return false;
    }
    let r: Vec<f64> = roots(p)
        .into_iter()
        .filter(|z| z.im.abs() <= 1e-9 * 1e-12f64.max(cabs(*z)))
        .map(|z| z.re)
        .filter(|v| *v >= -1e-12 && *v <= 1.0 + 1e-12 && v.abs() <= 1e-6)
        .collect();
    let dtr: Vec<f64> = dt.iter().rev().cloned().collect();
    let u1r: Vec<f64> = u1.iter().rev().cloned().collect();
    let u2r: Vec<f64> = u2.iter().rev().cloned().collect();
    for rj in r {
        let d = polyval_f(&dtr, rj);
        if d == 0.0 {
            continue;
        }
        let sj = [polyval_f(&u1r, rj) / d, polyval_f(&u2r, rj) / d];
        if !(sj[0].is_finite() && sj[1].is_finite()) {
            continue;
        }
        if !pf.is_empty() {
            let lim = -1e-9 * 1f64.max(sj[0].abs()).max(sj[1].abs());
            if pf.iter().any(|q| q[0] * sj[0] + q[1] * sj[1] + q[2] < lim) {
                continue;
            }
        }
        return true;
    }
    false
}

// ==========================================================================================
// bridge_witness
// ==========================================================================================
fn sign_at_rational(c: &[Z], xn: &[Z], xd: &Z) -> i32 {
    sign(&(&c[0] * &xn[0] * &xn[0] + &c[1] * &xn[0] * &xn[1] + &c[2] * &xn[1] * &xn[1] + &c[3] * &xn[0] * xd + &c[4] * &xn[1] * xd + &c[5] * xd * xd))
}

fn rev(a: &[Z]) -> Vec<Z> {
    a.iter().rev().cloned().collect()
}

fn member_exact_arc(pr: &AlgPair, xn: &[Z], xd: &Z) -> R<(bool, usize)> {
    let ph = bridge_phi(pr, xn, xd);
    let phr = poly_trim(&rev(&ph.phi));
    if phr.len() <= 1 {
        return Ok((false, 0));
    }
    let iv = isolate_roots(&phr)?;
    let mut n = 0;
    let (m1, z0) = (qi(-1), qi(0));
    for i in &iv {
        let lo = if i.0 > m1 { i.0.clone() } else { m1.clone() };
        let hi = if i.1 < z0 { i.1.clone() } else { z0.clone() };
        if lo > hi {
            continue;
        }
        if count_roots_in_q(&phr, &lo, &hi) != 1 && !polyval_q(&phr, &lo).is_zero() {
            continue;
        }
        n += 1;
        let sd = sign_at(&rev(&ph.dt), &phr, i)?;
        if sd == 0 {
            continue;
        }
        let mut ok = true;
        for r in &pr.p {
            let q = vadd(&vadd(&sc(&ph.u1, &r[0]), &sc(&ph.u2, &r[1])), &sc(&ph.dt, &(&r[2] * xd)));
            let sq = sign_at(&rev(&q), &phr, i)?;
            if sq * sd <= 0 {
                ok = false;
                break;
            }
        }
        if ok {
            return Ok((true, n));
        }
    }
    Ok((false, n))
}

fn member_exact(pr: &AlgPair, xn: &[Z], xd: &Z) -> R<(bool, usize)> {
    let ph = bridge_phi(pr, xn, xd);
    let phr = poly_trim(&rev(&ph.phi));
    if phr.len() <= 1 {
        return Ok((false, 0));
    }
    let n = count_roots_in_q(&phr, &qi(-1), &qi(0));
    Ok((n > 0 || polyval(&phr, &zi(-1)).is_zero(), n))
}

fn in_p(s: [f64; 2], p: &[[f64; 3]]) -> bool {
    p.iter().all(|r| s[0] * r[0] + s[1] * r[1] + r[2] >= 0.0)
}

fn arc_extent(s0: [f64; 2], w: [f64; 2], p: &[[f64; 3]]) -> (f64, f64) {
    let mut lo = f64::NEG_INFINITY;
    let mut hi = f64::INFINITY;
    for r in p {
        let a = r[0] * w[0] + r[1] * w[1];
        let b = r[0] * s0[0] + r[1] * s0[1] + r[2];
        if a.abs() <= 0.0 {
            if b < 0.0 {
                return (f64::NAN, f64::NAN);
            }
            continue;
        }
        if a > 0.0 {
            lo = lo.max(-b / a);
        } else {
            hi = hi.min(-b / a);
        }
    }
    if lo >= hi {
        return (f64::NAN, f64::NAN);
    }
    let span = 24.0;
    if lo.is_infinite() {
        lo = if hi.is_infinite() { -12.0 } else { hi - span };
    }
    if hi.is_infinite() {
        hi = if lo.is_infinite() { 12.0 } else { lo + span };
    }
    (lo, hi)
}

fn conic_points(dl: &[[f64; 2]; 2], dlv: [f64; 2], kp: f64, ctr: f64, l: f64, n: usize) -> Vec<[f64; 2]> {
    let mut s = vec![];
    for t in linspace(-l, l, n) {
        let s1 = ctr + t;
        let a = dl[1][1];
        let b = 2.0 * dl[0][1] * s1 + 2.0 * dlv[1];
        let c = dl[0][0] * s1 * s1 + 2.0 * dlv[0] * s1 + 2.0 * kp;
        if a.abs() < EPS * 1f64.max(b.abs()) {
            if b.abs() > 0.0 {
                s.push([s1, -c / b]);
            }
            continue;
        }
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            continue;
        }
        let r = disc.sqrt();
        s.push([s1, (-b + r) / (2.0 * a)]);
        if r > 0.0 {
            s.push([s1, (-b - r) / (2.0 * a)]);
        }
    }
    s.into_iter().filter(|p| p[0].is_finite() && p[1].is_finite() && p[0].abs() < 1e6 && p[1].abs() < 1e6).collect()
}

fn tie_sample(dlm: &[Vec<Z>], dl: &[Z], kp: &Z, p: &[Vec<Z>]) -> R<Vec<[f64; 2]>> {
    let pf: Vec<[f64; 3]> = p.iter().map(|r| [to_f64(&r[0]), to_f64(&r[1]), to_f64(&r[2])]).collect();
    if dlm.iter().flatten().all(|x| x.is_zero()) {
        let (d, e) = (&dl[0], &dl[1]);
        if d.is_zero() && e.is_zero() {
            return Ok(vec![]);
        }
        let w = [to_f64(e), -to_f64(d)];
        let nn = to_f64(&(d * d + e * e));
        // -kp * D / nn: an int product, then int / float
        let s0 = [to_f64(&(-kp * d)) / nn, to_f64(&(-kp * e)) / nn];
        let (lo, hi) = arc_extent(s0, w, &pf);
        if lo.is_nan() {
            return Ok(vec![]);
        }
        let taus = linspace(lo, hi, 43);
        let s: Vec<[f64; 2]> = taus[1..42].iter().map(|t| [s0[0] + w[0] * t, s0[1] + w[1] * t]).filter(|s| s[0].is_finite() && s[1].is_finite()).collect();
        return Ok(s);
    }
    let df = [[to_f64(&dlm[0][0]), to_f64(&dlm[0][1])], [to_f64(&dlm[1][0]), to_f64(&dlm[1][1])]];
    let dfv = [to_f64(&dl[0]), to_f64(&dl[1])];
    let kpf = to_f64(kp);
    let (l, ctr);
    if !det_exact(dlm)?.is_zero() {
        let sol = crate::conic_meet::solve2f(&df, dfv);
        let c0 = [-sol[0], -sol[1]];
        let cdc = c0[0] * (df[0][0] * c0[0] + df[0][1] * c0[1]) + c0[1] * (df[1][0] * c0[0] + df[1][1] * c0[1]);
        let v0 = cdc / 2.0 + (dfv[0] * c0[0] + dfv[1] * c0[1]) + kpf;
        let (e1, e2) = eig2sym(df[0][0], df[0][1], df[1][1]);
        let ev = e1.abs().max(e2.abs());
        l = 4.0 * ((2.0 * v0).abs() / ev.max(EPS)).sqrt() + 4.0;
        ctr = c0[0];
    } else {
        l = 20.0;
        ctr = 0.0;
    }
    let mut s: Vec<[f64; 2]> = vec![];
    for l2 in [l, 4.0 * l, 20.0 * l, 100.0 * l, 500.0 * l] {
        s = conic_points(&df, dfv, kpf, ctr, l2, 2001);
        if p.is_empty() || s.is_empty() {
            break;
        }
        let keep: Vec<[f64; 2]> = s.iter().cloned().filter(|q| in_p(*q, &pf)).collect();
        if keep.is_empty() {
            continue;
        }
        s = keep;
        if s.len() >= 5 {
            let lo = s.iter().map(|q| q[0]).fold(f64::INFINITY, f64::min);
            let hi = s.iter().map(|q| q[0]).fold(f64::NEG_INFINITY, f64::max);
            let pad = 1e-9f64.max((hi - lo) / 20.0);
            let s2 = conic_points(&df, dfv, kpf, (lo + hi) / 2.0, (hi - lo) / 2.0 + pad, 2001);
            if !s2.is_empty() {
                let k2: Vec<[f64; 2]> = s2.into_iter().filter(|q| in_p(*q, &pf)).collect();
                if !k2.is_empty() {
                    s = k2;
                }
            }
        }
        break;
    }
    if !s.is_empty() {
        let mut mx = 0.0;
        let mut my = 0.0;
        for q in &s {
            mx += q[0];
            my += q[1];
        }
        mx /= s.len() as f64;
        my /= s.len() as f64;
        s.sort_by(|a, b| {
            let da = ((a[0] - mx).powi(2) + (a[1] - my).powi(2)).sqrt();
            let db = ((b[0] - mx).powi(2) + (b[1] - my).powi(2)).sqrt();
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    Ok(s)
}

/// Eigenvalues (ascending) of a symmetric 2x2.
pub fn eig2sym(a: f64, b: f64, c: f64) -> (f64, f64) {
    let t = (a + c) / 2.0;
    let d = ((a - c) / 2.0).hypot(b);
    (t - d, t + d)
}

pub fn bridge_witness(pr: &AlgPair, c1: &[Z], c2: &[Z]) -> R<Vec<[Z; 3]>> {
    let p = &pr.p;
    let s = tie_sample(&pr.dlm, &pr.dl, &pr.kp, p)?;
    let n1: Vec<Vec<f64>> = pr.n1.iter().map(|r| r.iter().map(to_f64).collect()).collect();
    let n2: Vec<Vec<f64>> = pr.n2.iter().map(|r| r.iter().map(to_f64).collect()).collect();
    let m1: Vec<f64> = pr.m1.iter().map(to_f64).collect();
    let m2: Vec<f64> = pr.m2.iter().map(to_f64).collect();
    let d2 = 2.0 * to_f64(&pr.d);
    let mut tries = 0;
    let mut w: Vec<[Z; 3]> = vec![];
    for q in &s {
        let x0 = [
            ((n1[0][0] * q[0] + n1[0][1] * q[1] + m1[0]) + (n2[0][0] * q[0] + n2[0][1] * q[1] + m2[0])) / d2,
            ((n1[1][0] * q[0] + n1[1][1] * q[1] + m1[1]) + (n2[1][0] * q[0] + n2[1][1] * q[1] + m2[1])) / d2,
        ];
        if !(x0[0].is_finite() && x0[1].is_finite()) {
            continue;
        }
        for xd_try in [1i64, 2, 4, 8, 16, 32, 64] {
            tries += 1;
            let xn = [Z::from(round_half_away(x0[0] * xd_try as f64)), Z::from(round_half_away(x0[1] * xd_try as f64))];
            let xd = zi(xd_try);
            if (!c1.is_empty() && sign_at_rational(c1, &xn, &xd) == 0) || (!c2.is_empty() && sign_at_rational(c2, &xn, &xd) == 0) {
                continue;
            }
            let (in_r, _) = if p.is_empty() { member_exact(pr, &xn, &xd)? } else { member_exact_arc(pr, &xn, &xd)? };
            if in_r {
                let row = [xn[0].clone(), xn[1].clone(), xd];
                if !w.contains(&row) {
                    w.push(row);
                }
                break;
            }
        }
        if w.len() >= 12 {
            break;
        }
    }
    if !w.is_empty() {
        let mut order: Vec<usize> = (0..w.len()).collect();
        order.sort_by(|&i, &j| {
            let ki = (w[i][2].clone(), w[i][0].abs().max(w[i][1].abs()));
            let kj = (w[j][2].clone(), w[j][0].abs().max(w[j][1].abs()));
            ki.cmp(&kj)
        });
        return Ok(order.into_iter().map(|i| w[i].clone()).collect());
    }
    if !p.is_empty() && !s.is_empty() && s.len() <= 2 {
        bail!("bridgeFaceAlg:tangentialArc", "the tie meets this pair in at most a POINT ({} sampled tie points lie in it).", s.len());
    }
    bail!("bridgeFaceAlg:noWitness", "no interior point of the bridge region could be VERIFIED, after {} attempts over {} tie points.", tries, s.len())
}

// ==========================================================================================
// bridge_arc_cut
// ==========================================================================================
pub fn bridge_arc_cut(tie: &[Z], p: &[Vec<Z>], pr: &AlgPair, wit: &[Z; 3]) -> R<Vec<Vec<Vec<Z>>>> {
    if !all_zero(&tie[0..3]) {
        bail!("bridgeArcCut:conicTie", "the tie is a genuine CONIC, so each cut boundary is a SEGMENT of an irrational line.");
    }
    let (d, e, f) = (&tie[3], &tie[4], &tie[5]);
    if d.is_zero() && e.is_zero() {
        bail!("bridgeArcCut:degenerateTie", "the tie has no linear part either, so it names no curve at all.");
    }
    let w = vec![e.clone(), -d];
    let s0n = vec![-f * d, -f * e];
    let nd = d * d + e * e;
    let mut lo: Option<[Z; 2]> = None;
    let mut hi: Option<[Z; 2]> = None;
    for (i, r) in p.iter().enumerate() {
        let a = &r[0..2];
        let c = &r[2];
        let aa = &a[0] * &w[0] + &a[1] * &w[1];
        let bb = &a[0] * &s0n[0] + &a[1] * &s0n[1] + c * &nd;
        if aa.is_zero() {
            if bb.is_negative() {
                bail!("bridgeArcCut:emptyArc", "facet {} excludes the ENTIRE tie line.", i + 1);
            }
            continue;
        }
        let b = if aa.is_positive() { [-&bb, aa.clone()] } else { [bb.clone(), -&aa] };
        if aa.is_positive() {
            if lo.as_ref().map_or(true, |l| &b[0] * &l[1] > &l[0] * &b[1]) {
                lo = Some(b);
            }
        } else if hi.as_ref().map_or(true, |h| &b[0] * &h[1] < &h[0] * &b[1]) {
            hi = Some(b);
        }
    }
    if let (Some(l), Some(h)) = (&lo, &hi) {
        if &h[0] * &l[1] < &l[0] * &h[1] {
            bail!("bridgeArcCut:emptyArc", "the facets leave no tie point at all: the arc's bounds cross.");
        }
    }
    let mut rows = vec![];
    for tau in [lo, hi].iter().flatten() {
        if let Some(row) = endpoint_row(tau, &w, &s0n, &nd, pr, wit)? {
            rows.push(row);
        }
    }
    Ok(vec![rows])
}

fn endpoint_row(tau: &[Z; 2], w: &[Z], s0n: &[Z], nd: &Z, pr: &AlgPair, wit: &[Z; 3]) -> R<Option<Vec<Z>>> {
    let (tn, td) = (&tau[0], &tau[1]);
    let sen = vec![&s0n[0] * td + tn * &w[0], &s0n[1] * td + tn * &w[1]];
    let sed = nd * td;
    let p1n: Vec<Z> = mv(&pr.n1, &sen).iter().zip(pr.m1.iter()).map(|(a, b)| a + &sed * b).collect();
    let p2n: Vec<Z> = mv(&pr.n2, &sen).iter().zip(pr.m2.iter()).map(|(a, b)| a + &sed * b).collect();
    let dr = [&p2n[0] - &p1n[0], &p2n[1] - &p1n[1]];
    if dr[0].is_zero() && dr[1].is_zero() {
        return Ok(None);
    }
    let n = [-&dr[1], dr[0].clone()];
    let c = conic(&[zi(0), zi(0), zi(0), &pr.d * &sed * &n[0], &pr.d * &sed * &n[1], -(&n[0] * &p1n[0] + &n[1] * &p1n[1])])?;
    let s = sign_at_rational(&c, &wit[0..2], &wit[2]);
    if s == 0 {
        bail!("bridgeArcCut:witnessOnCut", "the witness lies exactly ON a cut boundary, so it names no side.");
    }
    let mut row = c;
    row.push(zi(s as i64));
    row.push(zi(0));
    Ok(Some(row))
}

// ==========================================================================================
// bridge_face_alg
// ==========================================================================================
type MuPoly = Vec<PolyN>;

fn mu_mul(a: &MuPoly, b: &MuPoly) -> MuPoly {
    let mut c: MuPoly = vec![PolyN::new(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_empty() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if y.is_empty() {
                continue;
            }
            c[i + j] = pn::add(&c[i + j], &pn::mul(x, y));
        }
    }
    c
}

fn mu_add(a: &MuPoly, b: &MuPoly) -> MuPoly {
    let n = a.len().max(b.len());
    (0..n)
        .map(|i| {
            let mut c = PolyN::new();
            if i < a.len() {
                c = pn::add(&c, &a[i]);
            }
            if i < b.len() {
                c = pn::add(&c, &b[i]);
            }
            c
        })
        .collect()
}

fn mu_scale(a: &MuPoly, k: &Z) -> MuPoly {
    a.iter().map(|x| pn::scale(x, k)).collect()
}

pub fn implicit_poly(d: &Z, n1: &[Vec<Z>], m1: &[Z], q1: &Z, n2: &[Vec<Z>], m2: &[Z], q2: &Z) -> R<PolyN> {
    let x = pn::mono(&[1, 0, 0], zi(1));
    let y = pn::mono(&[0, 1, 0], zi(1));
    let z = pn::mono(&[0, 0, 1], zi(1));
    let u = pn::mono(&[0, 0, 0], zi(1));
    let num = |v: &[Z]| -> MuPoly { v.iter().map(|c| pn::scale(&u, c)).collect() };
    let dl: Vec<Vec<Z>> = (0..2).map(|i| (0..2).map(|j| &n1[i][j] - &n2[i][j]).collect()).collect();
    let e11 = vec![n1[0][0].clone(), dl[0][0].clone()];
    let e12 = vec![n1[0][1].clone(), dl[0][1].clone()];
    let e22 = vec![n1[1][1].clone(), dl[1][1].clone()];
    let dt: Vec<Z> = conv(&e11, &e22).iter().zip(conv(&e12, &e12).iter()).map(|(a, b)| a - b).collect();
    let a11 = e22.clone();
    let a12 = neg(&e12);
    let a22 = e11.clone();
    let c1 = [m1[0].clone(), &m1[0] - &m2[0]];
    let c2 = [m1[1].clone(), &m1[1] - &m2[1]];
    let kk = [q1.clone(), q1 - q2];
    let w1: MuPoly = vec![pn::add(&pn::scale(&x, d), &pn::scale(&u, &-&c1[0])), pn::scale(&u, &-&c1[1])];
    let w2: MuPoly = vec![pn::add(&pn::scale(&y, d), &pn::scale(&u, &-&c2[0])), pn::scale(&u, &-&c2[1])];
    let term1 = mu_scale(&mu_mul(&vec![pn::add(&pn::scale(&z, d), &pn::scale(&u, &kk[0])), pn::scale(&u, &kk[1])], &num(&dt)), &zi(2));
    let term2 = mu_add(
        &mu_add(&mu_mul(&num(&a11), &mu_mul(&w1, &w1)), &mu_scale(&mu_mul(&num(&a12), &mu_mul(&w1, &w2)), &zi(2))),
        &mu_mul(&num(&a22), &mu_mul(&w2, &w2)),
    );
    let mut psi = mu_add(&term1, &mu_scale(&term2, &zi(-1)));
    let mut last = 0;
    for (i, p) in psi.iter().enumerate() {
        if !p.is_empty() {
            last = i + 1;
        }
    }
    psi.truncate(last.max(1));
    let deg = psi.len() - 1;
    let m = pn::mul;
    let p = if deg == 3 {
        let (a, b, c, e) = (&psi[3], &psi[2], &psi[1], &psi[0]);
        let mut p = pn::add(&pn::scale(&m(a, &m(b, &m(c, e))), &zi(18)), &pn::scale(&m(b, &m(b, &m(b, e))), &zi(-4)));
        p = pn::add(&p, &m(&m(b, b), &m(c, c)));
        p = pn::add(&p, &pn::scale(&m(a, &m(c, &m(c, c))), &zi(-4)));
        p = pn::add(&p, &pn::scale(&m(&m(a, a), &m(e, e)), &zi(-27)));
        p
    } else if deg == 2 {
        let (b, c, e) = (&psi[2], &psi[1], &psi[0]);
        pn::add(&m(c, c), &pn::scale(&m(b, e), &zi(-4)))
    } else {
        bail!("bridgeFaceAlg:straightTie", "Psi has degree {} in the multiplier, which happens when the tie is STRAIGHT.", deg);
    };
    let mut p = pn::content_free(&p);
    if pn::degrees(&p, 3)[2] < 1 {
        bail!("bridgeFaceAlg:noZ", "the eliminated polynomial has degree 0 in z.");
    }
    let cz = pn::slices(&p, 2);
    let nz = *cz.keys().max().unwrap();
    let g = pn::content_free(&cz[&nz]);
    let is_unit = g.len() == 1 && g.get(&vec![0u32, 0]).map_or(false, |c| c.abs().is_one());
    if !is_unit {
        let mut q: BTreeMap<u32, PolyN> = BTreeMap::new();
        let mut ok = true;
        for (k, c) in &cz {
            match pn::div_exact(c, &g) {
                Ok(v) => {
                    q.insert(*k, v);
                }
                Err(me) => {
                    if !me.id.starts_with("ratQ:") {
                        return Err(me);
                    }
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            p = pn::content_free(&pn::from_slices(&q, 2));
        }
    }
    Ok(p)
}

pub fn image_conic(tie: &[Z], nk: &[Vec<Z>], mk: &[Z], d: &Z) -> R<Vec<Z>> {
    if det_exact(nk)?.is_zero() {
        return Ok(vec![]);
    }
    let adj = vec![vec![nk[1][1].clone(), -&nk[0][1]], vec![-&nk[1][0], nk[0][0].clone()]];
    let mn: Vec<Vec<Z>> = adj.iter().map(|r| r.iter().map(|x| d * x).collect()).collect();
    let tn = vec![-(&adj[0][0] * &mk[0] + &adj[0][1] * &mk[1]), -(&adj[1][0] * &mk[0] + &adj[1][1] * &mk[1])];
    let md = det_exact(nk)?;
    let q = tie;
    let (n2, _) = subst_affine(&[zi(0), zi(0), zi(0), zi(0), zi(2) * &q[0], q[1].clone(), zi(2) * &q[2], q[3].clone(), q[4].clone(), q[5].clone()], &zi(1), &mn, &tn, &md)?;
    conic(&face_conic(&n2))
}

/// FP restricted to the vertical line over x = xn/xd: integer coefficients in z, highest first.
pub fn fibre(fp: &PolyN, xn: &[Z], xd: &Z) -> Vec<Z> {
    let deg = fp.keys().map(|e| e[0] + e[1]).max().unwrap_or(0);
    let nz = fp.keys().map(|e| e[2]).max().unwrap_or(0) as usize;
    let mut acc = vec![Z::zero(); nz + 1];
    for (e, c) in fp {
        acc[e[2] as usize] += c * num_traits::pow(xn[0].clone(), e[0] as usize) * num_traits::pow(xn[1].clone(), e[1] as usize) * num_traits::pow(xd.clone(), (deg - e[0] - e[1]) as usize);
    }
    let mut pz: Vec<Z> = acc.into_iter().rev().collect();
    let g = gcd_all(&pz).abs();
    if g > Z::one() {
        pz = pz.iter().map(|v| v / &g).collect();
    }
    pz
}

fn dyadic(x: f64, floor: bool) -> Q {
    let s = 2f64.powi(20);
    let v = if floor { (x * s).floor() } else { (x * s).ceil() };
    q_from_f64(v) / Q::from_integer(Z::from(1i64 << 20))
}

fn psi_value(mu: f64, pr: &AlgPair, q1: f64, q2: f64, x: [f64; 2]) -> f64 {
    let n1: Vec<Vec<f64>> = pr.n1.iter().map(|r| r.iter().map(to_f64).collect()).collect();
    let n2: Vec<Vec<f64>> = pr.n2.iter().map(|r| r.iter().map(to_f64).collect()).collect();
    let m = [[(1.0 + mu) * n1[0][0] - mu * n2[0][0], (1.0 + mu) * n1[0][1] - mu * n2[0][1]], [(1.0 + mu) * n1[1][0] - mu * n2[1][0], (1.0 + mu) * n1[1][1] - mu * n2[1][1]]];
    let c = [(1.0 + mu) * to_f64(&pr.m1[0]) - mu * to_f64(&pr.m2[0]), (1.0 + mu) * to_f64(&pr.m1[1]) - mu * to_f64(&pr.m2[1])];
    let k = (1.0 + mu) * q1 - mu * q2;
    let dt = m[0][0] * m[1][1] - m[0][1] * m[1][0];
    if dt == 0.0 {
        return f64::INFINITY;
    }
    let adj = [[m[1][1], -m[0][1]], [-m[1][0], m[0][0]]];
    let d = to_f64(&pr.d);
    let w = [d * x[0] - c[0], d * x[1] - c[1]];
    let waw = w[0] * (adj[0][0] * w[0] + adj[0][1] * w[1]) + w[1] * (adj[1][0] * w[0] + adj[1][1] * w[1]);
    (waw / (2.0 * dt) - k) / d
}

fn branch_index(fp: &PolyN, pr: &AlgPair, q1: &Z, q2: &Z, xn: &[Z], xd: &Z) -> R<i64> {
    let ph = bridge_phi(pr, xn, xd);
    let phi = poly_trim(&rev(&ph.phi));
    let x = [q_to_f64(&Q::new(xn[0].clone(), xd.clone())), q_to_f64(&Q::new(xn[1].clone(), xd.clone()))];
    let phif: Vec<f64> = phi.iter().map(to_f64).collect();
    let mut mus: Vec<f64> = roots(&phif).into_iter().filter(|z| z.im.abs() <= 1e-10 * 1f64.max(cabs(*z))).map(|z| z.re).filter(|m| *m >= -1.0 - 1e-12 && *m <= 1e-12).collect();
    if !mus.is_empty() {
        let dscale = ph.dt.iter().map(|v| v.abs()).max().unwrap();
        let dtr: Vec<f64> = ph.dt.iter().rev().map(to_f64).collect();
        let lim = 1e-9 * 1f64.max(to_f64(&dscale));
        mus.retain(|m| polyval_f(&dtr, *m).abs() > lim);
    }
    if mus.is_empty() {
        bail!("bridgeFaceAlg:noStationary", "the witness was VERIFIED to be in the region, yet no root of Phi survived the root finder.");
    }
    let (q1f, q2f) = (to_f64(q1), to_f64(q2));
    let zs: Vec<f64> = mus.iter().map(|m| psi_value(*m, pr, q1f, q2f, x)).filter(|z| z.is_finite()).collect();
    if zs.is_empty() {
        bail!("bridgeFaceAlg:noStationary", "every stationary multiplier at the witness gave a non-finite value.");
    }
    let zstar = zs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let pz = poly_trim(&fibre(fp, xn, xd));
    if pz.len() <= 1 {
        bail!("bridgeFaceAlg:emptyFibre", "FP degenerates to a constant at the witness.");
    }
    let pzf: Vec<f64> = pz.iter().map(to_f64).collect();
    let mut rr: Vec<f64> = roots(&pzf).into_iter().filter(|z| z.im.abs() <= 1e-10 * 1f64.max(cabs(*z))).map(|z| z.re).collect();
    rr.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let scale = 1f64.max(zstar.abs());
    let gap = rr.iter().map(|r| (r - zstar).abs()).fold(f64::INFINITY, f64::min);
    if gap > 1e-6 * scale {
        bail!("bridgeFaceAlg:valueNotARoot", "the value {} computed from the multiplier is not a root of FP at the witness.", zstar);
    }
    let w = 1e-7f64.max(1e-7 * scale);
    let lo = dyadic(zstar - w, true);
    let hi = dyadic(zstar + w, false);
    let mx = pz[1..].iter().map(|v| v.abs()).max().unwrap();
    use num_integer::Integer;
    let m = zi(1) + (-((-mx).div_floor(&pz[0].abs()))) + zi(1);
    let in_bracket = count_roots_in_q(&pz, &lo, &hi);
    let above = count_roots_in_q(&pz, &hi, &Q::from_integer(m));
    if in_bracket != 1 || above != 0 {
        bail!("bridgeFaceAlg:branchUndecided", "the top branch could not be CERTIFIED at the witness: {} root(s) in the bracket and {} above it.", in_bracket, above);
    }
    Ok(-1)
}

fn side_of(c: &[Z], xn: &[Z], xd: &Z) -> i64 {
    if c.is_empty() {
        return 0;
    }
    sign_at_rational(c, xn, xd) as i64
}

pub fn bridge_face_alg(n1: &[Z], d1: &Z, n2: &[Z], d2: &Z, p: &[Vec<Z>]) -> R<BFace> {
    if !all_zero(&n1[0..4]) || !all_zero(&n2[0..4]) {
        bail!("bridgeFaceAlg:cubicFace", "a face has a CUBIC part.");
    }
    if !d1.is_positive() || !d2.is_positive() {
        bail!("bridgeFaceAlg:denominator", "denominators must be positive integers.");
    }
    let d = lcm(d1, d2);
    let k1 = &d / d1;
    let k2 = &d / d2;
    let nn1 = vec![vec![&n1[4] * &k1, &n1[5] * &k1], vec![&n1[5] * &k1, &n1[6] * &k1]];
    let m1 = vec![&n1[7] * &k1, &n1[8] * &k1];
    let q1 = &n1[9] * &k1;
    let nn2 = vec![vec![&n2[4] * &k2, &n2[5] * &k2], vec![&n2[5] * &k2, &n2[6] * &k2]];
    let m2 = vec![&n2[7] * &k2, &n2[8] * &k2];
    let q2 = &n2[9] * &k2;
    let zero = |n: &Vec<Vec<Z>>| n.iter().flatten().all(|x| x.is_zero());
    let r1 = det_exact(&nn1)?.is_zero() && !zero(&nn1);
    let r2 = det_exact(&nn2)?.is_zero() && !zero(&nn2);
    if zero(&nn1) && r2 {
        return bridge_face_rat(n1, d1, n2, d2, p);
    }
    if zero(&nn2) && r1 {
        return bridge_face_rat(n2, d2, n1, d1, p);
    }
    let dlm: Vec<Vec<Z>> = (0..2).map(|i| (0..2).map(|j| &nn1[i][j] - &nn2[i][j]).collect()).collect();
    let dl = vec![&m1[0] - &m2[0], &m1[1] - &m2[1]];
    let kp = &q1 - &q2;
    let tie = diff_conic(n1, d1, n2, d2)?;
    let fp = implicit_poly(&d, &nn1, &m1, &q1, &nn2, &m2, &q2)?;
    let c1 = image_conic(&tie, &nn1, &m1, &d)?;
    let c2 = image_conic(&tie, &nn2, &m2, &d)?;
    let pr = AlgPair { d: d.clone(), n1: nn1.clone(), m1: m1.clone(), n2: nn2.clone(), m2: m2.clone(), dlm, dl, kp, p: p.to_vec(), fcache: Default::default() };
    let w = bridge_witness(&pr, &c1, &c2)?;
    let mut fb: Option<i64> = None;
    let mut chosen: Option<(Vec<Z>, Z, i64, i64)> = None;
    let mut last_msg = String::new();
    for wrow in &w {
        let xn = wrow[0..2].to_vec();
        let xd = wrow[2].clone();
        let s1 = side_of(&c1, &xn, &xd);
        let s2 = side_of(&c2, &xn, &xd);
        match branch_index(&fp, &pr, &q1, &q2, &xn, &xd) {
            Ok(b) => {
                fb = Some(b);
                chosen = Some((xn, xd, s1, s2));
                break;
            }
            Err(me) => {
                if !["ratQ:overflow", "bridgeFaceAlg:branchUndecided", "ratQ:sturmEndpointIsRoot"].contains(&me.id.as_str()) {
                    return Err(me);
                }
                last_msg = me.msg;
            }
        }
    }
    let fb = match fb {
        None => bail!("bridgeFaceAlg:branchUndecided", "none of the {} verified witnesses could certify the branch index; the last one reported \"{}\".", w.len(), last_msg),
        Some(b) => b,
    };
    let (xn, xd, s1, s2) = chosen.unwrap();
    let mut base = vec![];
    if !c1.is_empty() {
        let mut r = c1.clone();
        r.push(zi(s1));
        r.push(zi(0));
        base.push(r);
    }
    if !c2.is_empty() {
        let mut r = c2.clone();
        r.push(zi(s2));
        r.push(zi(0));
        base.push(r);
    }
    let wit = [xn[0].clone(), xn[1].clone(), xd];
    let mut cells = vec![base.clone()];
    if !p.is_empty() {
        match bridge_arc_cut(&tie, p, &pr, &wit) {
            Ok(blocks) => {
                cells = blocks
                    .into_iter()
                    .map(|b| {
                        let mut c = base.clone();
                        c.extend(b);
                        c
                    })
                    .collect();
            }
            Err(me) => {
                if me.id != "bridgeArcCut:conicTie" {
                    return Err(me);
                }
            }
        }
    }
    let mut face = BFace::blank(FaceKind::Algebraic);
    face.fp = Arc::new(fp);
    face.fb = fb;
    face.cells = cells;
    face.pair = Some(Arc::new(PairData::Alg(pr)));
    face.wit = Some(wit);
    Ok(face)
}

// ==========================================================================================
// bridge_face_rat
// ==========================================================================================
fn poly_mul_aff(f: &[Z], g: &[Z]) -> Vec<Z> {
    vec![
        &f[0] * &g[0],
        &f[0] * &g[1] + &f[1] * &g[0],
        &f[1] * &g[1],
        &f[0] * &g[2] + &f[2] * &g[0],
        &f[1] * &g[2] + &f[2] * &g[1],
        &f[2] * &g[2],
    ]
}

pub fn bridge_face_rat(n1: &[Z], d1: &Z, n2: &[Z], d2: &Z, p: &[Vec<Z>]) -> R<BFace> {
    if !all_zero(&n1[4..7]) {
        bail!("bridgeFaceRat:notAffine", "n1 must be the AFFINE face.");
    }
    let nn2 = vec![vec![n2[4].clone(), n2[5].clone()], vec![n2[5].clone(), n2[6].clone()]];
    if nn2.iter().flatten().all(|x| x.is_zero()) || !det_exact(&nn2)?.is_zero() {
        bail!("bridgeFaceRat:notRankOne", "n2 must have a RANK-1 Hessian.");
    }
    let d = lcm(d1, d2);
    let b1 = vec![&n1[7] * (&d / d1), &n1[8] * (&d / d1)];
    let k1 = &n1[9] * (&d / d1);
    let mm2: Vec<Vec<Z>> = nn2.iter().map(|r| r.iter().map(|x| x * (&d / d2)).collect()).collect();
    let b2 = vec![&n2[7] * (&d / d2), &n2[8] * (&d / d2)];
    let k2 = &n2[9] * (&d / d2);
    let dl = vec![&b1[0] - &b2[0], &b1[1] - &b2[1]];
    let kp = &k1 - &k2;
    let (u, a) = if !mm2[0][0].is_zero() { (vec![mm2[0][0].clone(), mm2[0][1].clone()], mm2[0][0].clone()) } else { (vec![mm2[0][1].clone(), mm2[1][1].clone()], mm2[1][1].clone()) };
    if a.is_zero() {
        bail!("bridgeFaceRat:internal", "a rank-1 symmetric matrix has a nonzero diagonal entry; this one appears not to.");
    }
    let w = vec![-&u[1], u[0].clone()];
    let nu = &u[0] * &u[0] + &u[1] * &u[1];
    let dw = &dl[0] * &w[0] + &dl[1] * &w[1];
    if dw.is_zero() {
        bail!("bridgeFaceRat:noSeparateFace", "the tie is a pair of PARALLEL LINES and this face is AFFINE.");
    }
    let rq = nu.clone();
    let rl = -zi(2) * &a * (&dl[0] * &u[0] + &dl[1] * &u[1]);
    let rc = -zi(2) * &a * &kp * &nu;
    let rd = zi(2) * &a * &dw;
    let uy = vec![&d * &u[0], &d * &u[1], -(&u[0] * &b1[0] + &u[1] * &b1[1])];
    let wy = vec![&d * &w[0], &d * &w[1], -(&w[0] * &b1[0] + &w[1] * &b1[1])];
    let cc: Vec<Z> = wy.iter().map(|x| &rq * x).collect();
    let bb: Vec<Z> = uy.iter().zip(wy.iter()).map(|(x, y)| &rd * x + &rl * y).collect();
    let mut aa: Vec<Z> = wy.iter().map(|x| &rc * x).collect();
    aa[2] -= &nu * &rd * &k1;
    let scale = &nu * &d * &rd;
    let num: Vec<Z> = poly_mul_aff(&aa, &cc).iter().zip(poly_mul_aff(&bb, &bb).iter()).map(|(x, y)| zi(4) * x - y).collect();
    let den: Vec<Z> = cc.iter().map(|x| zi(4) * &scale * x).collect();
    let rfn = vec![zi(0), zi(0), zi(0), zi(0), zi(2) * &num[0], num[1].clone(), zi(2) * &num[2], num[3].clone(), num[4].clone(), num[5].clone()];
    let pr = RatPair { u, w, nu, rq, rl, rc, rd, b: bb, c: cc, scale, p: p.to_vec() };
    let wit = witness_rat(&pr)?;
    let cells = vec![region_rows(&pr)?];
    let fp = implicit_value(&rfn, &den);
    let mut face = BFace::blank(FaceKind::Rational);
    face.rfn = rfn;
    face.rfd = zi(1);
    face.rden = den;
    face.wit = Some(wit);
    face.cells = cells;
    face.fp = Arc::new(fp);
    face.fb = -1;
    face.pair = Some(Arc::new(PairData::Rat(pr)));
    Ok(face)
}

fn para_point(pr: &RatPair, t: f64) -> [f64; 2] {
    let r = (to_f64(&pr.rq) * t * t + to_f64(&pr.rl) * t + to_f64(&pr.rc)) / to_f64(&pr.rd);
    let nu = to_f64(&pr.nu);
    [(t * to_f64(&pr.u[0]) + r * to_f64(&pr.w[0])) / nu, (t * to_f64(&pr.u[1]) + r * to_f64(&pr.w[1])) / nu]
}

fn candidate_params(pr: &RatPair) -> Vec<f64> {
    let mut rts = vec![];
    for r in &pr.p {
        let aw = &r[0] * &pr.w[0] + &r[1] * &pr.w[1];
        let au = &r[0] * &pr.u[0] + &r[1] * &pr.u[1];
        let c2 = idiv(&(&aw * &pr.rq), &(&pr.rd * &pr.nu));
        let c1 = (to_f64(&au) + idiv(&(&aw * &pr.rl), &pr.rd)) / to_f64(&pr.nu);
        let c0 = idiv(&(&aw * &pr.rc), &(&pr.rd * &pr.nu)) + to_f64(&r[2]);
        for z in roots(&[c2, c1, c0]) {
            if z.im.abs() < 1e-12 * 1f64.max(cabs(z)) {
                rts.push(z.re);
            }
        }
    }
    rts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut t = if rts.is_empty() {
        vec![0.0, 1.0, -1.0]
    } else {
        let span = 1f64.max(rts[rts.len() - 1] - rts[0]);
        let mut t: Vec<f64> = rts.windows(2).map(|w| (w[0] + w[1]) / 2.0).collect();
        t.push(rts[0] - span);
        t.push(rts[rts.len() - 1] + span);
        t
    };
    t.extend([0.0, 1.0, -1.0, 2.0, -2.0, 0.5, -0.5, 4.0, -4.0, 10.0, -10.0]);
    t
}

fn member_rat_at(pr: &RatPair, xn: &[Z], xd: &Z) -> bool {
    let bv = idiv(&(&pr.b[0] * &xn[0] + &pr.b[1] * &xn[1]), xd) + to_f64(&pr.b[2]);
    let cv = idiv(&(&pr.c[0] * &xn[0] + &pr.c[1] * &xn[1]), xd) + to_f64(&pr.c[2]);
    if cv * to_f64(&pr.scale) >= 0.0 {
        return false;
    }
    let t = -bv / (2.0 * cv);
    let s = para_point(pr, t);
    pr.p.is_empty() || pr.p.iter().all(|r| s[0] * to_f64(&r[0]) + s[1] * to_f64(&r[1]) + to_f64(&r[2]) >= 0.0)
}

fn witness_rat(pr: &RatPair) -> R<[Z; 3]> {
    let bf: Vec<f64> = pr.b.iter().map(to_f64).collect();
    let cf: Vec<f64> = pr.c.iter().map(to_f64).collect();
    for t in candidate_params(pr) {
        let s = para_point(pr, t);
        if !(s[0].is_finite() && s[1].is_finite()) {
            continue;
        }
        if !pr.p.is_empty() && pr.p.iter().any(|r| s[0] * to_f64(&r[0]) + s[1] * to_f64(&r[1]) + to_f64(&r[2]) <= 0.0) {
            continue;
        }
        let l = [bf[0] + 2.0 * t * cf[0], bf[1] + 2.0 * t * cf[1], bf[2] + 2.0 * t * cf[2]];
        let g = [l[0], l[1]];
        if g[0] == 0.0 && g[1] == 0.0 {
            continue;
        }
        let gg = g[0] * g[0] + g[1] * g[1];
        let x0 = [-l[2] * g[0] / gg, -l[2] * g[1] / gg];
        let dr = [-g[1], g[0]];
        let c0 = cf[0] * x0[0] + cf[1] * x0[1] + cf[2];
        let cd = cf[0] * dr[0] + cf[1] * dr[1];
        let want = -(sign(&pr.scale) as f64) * 1f64.max(c0.abs());
        let x = if cd != 0.0 {
            let k = (want - c0) / cd;
            [x0[0] + k * dr[0], x0[1] + k * dr[1]]
        } else if c0 * to_f64(&pr.scale) < 0.0 {
            x0
        } else {
            continue;
        };
        for e in 0..31 {
            let xd = 2f64.powi(e);
            // the rounding in doubles (exact below 2^53), compared BEFORE any integer conversion, so
            // a huge or non-finite candidate ends the search as Python's bigint test does
            let xr = [crate::rat_q::round_half_away(x[0] * xd), crate::rat_q::round_half_away(x[1] * xd)];
            if !(xr[0].abs() < 2f64.powi(52)) || !(xr[1].abs() < 2f64.powi(52)) {
                break;
            }
            let xnz = [Z::from(xr[0] as i64), Z::from(xr[1] as i64)];
            let xdz = Z::from(1i64 << e);
            if !member_rat_at(pr, &xnz, &xdz) {
                continue;
            }
            return Ok([xnz[0].clone(), xnz[1].clone(), xdz]);
        }
    }
    bail!("bridgeFaceRat:noWitness", "no rational point of this face was found on the search grid.")
}

fn region_rows(pr: &RatPair) -> R<Vec<Vec<Z>>> {
    let mut r = vec![];
    let s = zi(-(sign(&pr.scale) as i64));
    let mut raw = vec![zi(0), zi(0), zi(0)];
    raw.extend(pr.c.iter().map(|x| &s * x));
    add_row(&mut r, &raw)?;
    let rs = zi(sign(&pr.rd) as i64);
    let bb = poly_mul_aff(&pr.b, &pr.b);
    let bc = poly_mul_aff(&pr.b, &pr.c);
    let cc = poly_mul_aff(&pr.c, &pr.c);
    for row in &pr.p {
        let a = &row[0..2];
        let c = &row[2];
        let aw = &a[0] * &pr.w[0] + &a[1] * &pr.w[1];
        let au = &a[0] * &pr.u[0] + &a[1] * &pr.u[1];
        let al = &aw * &pr.rq;
        let be = &au * &pr.rd + &aw * &pr.rl;
        let ga = &aw * &pr.rc + c * &pr.rd * &pr.nu;
        let raw: Vec<Z> = (0..6).map(|i| &rs * (&al * &bb[i] - zi(2) * &be * &bc[i] + zi(4) * &ga * &cc[i])).collect();
        add_row(&mut r, &raw)?;
    }
    Ok(r)
}

fn add_row(r: &mut Vec<Vec<Z>>, raw: &[Z]) -> R<()> {
    if all_zero(&raw[0..5]) {
        return Ok(());
    }
    let nz = raw.iter().find(|x| !x.is_zero()).unwrap();
    let mut row = conic(raw)?;
    row.push(zi(sign(nz) as i64));
    row.push(zi(0));
    r.push(row);
    Ok(())
}

fn implicit_value(fnn: &[Z], den: &[Z]) -> PolyN {
    let q = [Q::new(fnn[4].clone(), zi(2)), qz(&fnn[5]), Q::new(fnn[6].clone(), zi(2)), qz(&fnn[7]), qz(&fnn[8]), qz(&fnn[9])];
    let (q, den): (Vec<Z>, Vec<Z>) = if q.iter().any(|x| !x.denom().is_one()) {
        (q.iter().map(|x| (x * qi(2)).to_integer()).collect(), den.iter().map(|x| zi(2) * x).collect())
    } else {
        (q.iter().map(|x| x.to_integer()).collect(), den.to_vec())
    };
    let mut fp = PolyN::new();
    let mut put = |e: [u32; 3], c: Z| {
        if !c.is_zero() {
            fp.insert(e.to_vec(), c);
        }
    };
    put([2, 0, 0], -&q[0]);
    put([1, 1, 0], -&q[1]);
    put([0, 2, 0], -&q[2]);
    put([1, 0, 0], -&q[3]);
    put([0, 1, 0], -&q[4]);
    put([0, 0, 0], -&q[5]);
    put([1, 0, 1], den[0].clone());
    put([0, 1, 1], den[1].clone());
    put([0, 0, 1], den[2].clone());
    fp
}

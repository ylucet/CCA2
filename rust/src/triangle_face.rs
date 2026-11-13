//! Contact triangles (ports of triangle_face.py, dual_point.py, cut_triangles.py).

use crate::alg_con::{SubFace, SubKind};
use crate::bail;
use crate::biv::{self, M2};
use crate::cells::{cell_interior, side_i};
use crate::conic_meet::solve2f;
use crate::errors::R;
use crate::isolated_root_sign::isolated_root_sign;
use crate::numeric::polyval_f;
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::rc::Rc;

/// The triangle as biconjAtoms hands it over (doubles plus the forms' exact conics).
#[derive(Clone, Debug)]
pub struct AtomTri {
    pub t: [[f64; 2]; 3],
    pub s: [f64; 2],
    pub c: f64,
    pub g: Vec<[f64; 2]>,
    pub forms: Vec<usize>,
    pub gn: Vec<Vec<Z>>,
    pub gd: Vec<Z>,
}

#[derive(Clone, Debug)]
pub struct DualPoint {
    pub k: i64,
    pub r: Vec<Z>,
    pub iv: (Q, Q),
    pub s1: Vec<Z>,
    pub s0: Vec<Z>,
    pub s: [f64; 2],
}

#[derive(Clone, Debug)]
pub struct TriFace {
    pub rational: bool,
    pub s0n: Vec<Z>,
    pub s0d: Z,
    pub tn: Vec<Vec<Z>>,
    pub td: Vec<Z>,
    pub dp: Option<DualPoint>,
    pub gn: Vec<Vec<Z>>,
    pub gd: Vec<Z>,
    pub s: [f64; 2],
    pub c: f64,
    pub t: [[f64; 2]; 3],
    pub orient: i32,
}

pub fn triangle_face(tri: &AtomTri) -> R<TriFace> {
    if !tri.g.is_empty() {
        bail!("triangleFace:unbounded", "the face holds on T + cone(G) (erratum Case E); only the bounded triangle is carried so far.");
    }
    let gn = tri.gn.clone();
    let gd = tri.gd.clone();
    if gn.len() != 3 || gn.iter().any(|r| r.len() != 6) || gd.len() != 3 {
        bail!("triangleFace:noExactForms", "the forms have no exact conic, so s0 cannot be named.");
    }
    let t12 = conic(&gn[0].iter().zip(gn[1].iter()).map(|(a, b)| a * &gd[1] - b * &gd[0]).collect::<Vec<_>>())?;
    let t13 = conic(&gn[0].iter().zip(gn[2].iter()).map(|(a, b)| a * &gd[2] - b * &gd[0]).collect::<Vec<_>>())?;
    let mut tf = TriFace { rational: false, s0n: vec![], s0d: Z::one(), tn: vec![], td: vec![], dp: None, gn: gn.clone(), gd: gd.clone(), s: [0.0, 0.0], c: 0.0, t: [[0.0; 2]; 3], orient: 0 };
    if all_zero(&t12[0..3]) && all_zero(&t13[0..3]) {
        let (mut xn, mut xd) = solve2(&[t12[3..5].to_vec(), t13[3..5].to_vec()], &[-&t12[5], -&t13[5]])?;
        if xd.is_negative() {
            xn = neg(&xn);
            xd = -xd;
        }
        tf.rational = true;
        tf.s = [q_to_f64(&Q::new(xn[0].clone(), xd.clone())), q_to_f64(&Q::new(xn[1].clone(), xd.clone()))];
        for i in 0..3 {
            let g = &gn[i];
            let h = [[zi(2) * &g[0], g[1].clone()], [g[1].clone(), zi(2) * &g[2]]];
            let tn = vec![&h[0][0] * &xn[0] + &h[0][1] * &xn[1] + &g[3] * &xd, &h[1][0] * &xn[0] + &h[1][1] * &xn[1] + &g[4] * &xd];
            let (tn, td) = canon(&tn, &(&xd * &gd[i]))?;
            tf.tn.push(tn);
            tf.td.push(td);
        }
        tf.s0n = xn;
        tf.s0d = xd;
    } else {
        let dp = dual_point(&t12, &t13, tri.s)?;
        tf.s = dp.s;
        tf.dp = Some(dp);
    }
    let s = tf.s;
    let g: Vec<f64> = gn[0].iter().map(|v| q_to_f64(&Q::new(v.clone(), gd[0].clone()))).collect();
    tf.c = g[0] * s[0].powi(2) + g[1] * s[0] * s[1] + g[2] * s[1].powi(2) + g[3] * s[0] + g[4] * s[1] + g[5];
    for i in 0..3 {
        let g: Vec<f64> = gn[i].iter().map(|v| q_to_f64(&Q::new(v.clone(), gd[i].clone()))).collect();
        tf.t[i] = [2.0 * g[0] * s[0] + g[1] * s[1] + g[3], g[1] * s[0] + 2.0 * g[2] * s[1] + g[4]];
    }
    tf.orient = triangle_face_sign(&tf, &area_poly(&gn, &gd))?;
    if tf.orient == 0 {
        bail!("triangleFace:flat", "the contact triangle has zero area: its corners are collinear.");
    }
    Ok(tf)
}

fn contact_poly(g: &[Z]) -> [M2; 2] {
    [vec![vec![g[3].clone(), g[1].clone()], vec![zi(2) * &g[0], zi(0)]], vec![vec![g[4].clone(), zi(2) * &g[2]], vec![g[1].clone(), zi(0)]]]
}

fn area_poly(gn: &[Vec<Z>], gd: &[Z]) -> M2 {
    let m: Vec<[M2; 2]> = gn.iter().map(|g| contact_poly(g)).collect();
    let mut a = vec![vec![zi(0)]];
    for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
        let o = 3 - i - j;
        let cr = biv::pad_add(&biv::conv2(&m[i][0], &m[j][1]), &biv::scali(&biv::conv2(&m[i][1], &m[j][0]), -1));
        a = biv::pad_add(&a, &biv::scal(&cr, &gd[o]));
    }
    a
}

pub fn triangle_face_sign(tf: &TriFace, p: &M2) -> R<i32> {
    if !tf.rational {
        return dual_point_sign(tf.dp.as_ref().unwrap(), p);
    }
    let deg = (p.len() - 1) + (p[0].len() - 1);
    let mut acc = Z::zero();
    for (i, r) in p.iter().enumerate() {
        for (j, c) in r.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            acc += c * num_traits::pow(tf.s0n[0].clone(), i) * num_traits::pow(tf.s0n[1].clone(), j) * num_traits::pow(tf.s0d.clone(), deg - i - j);
        }
    }
    Ok(sign(&acc))
}

// ---- dual_point -------------------------------------------------------------------------------
fn padd(a: &[Z], b: &[Z]) -> Vec<Z> {
    poly_add1(a, b)
}
fn scal(a: &[Z], k: &Z) -> Vec<Z> {
    a.iter().map(|x| x * k).collect()
}

fn in_v(c: &[Z], k: i64) -> (Z, Vec<Z>, Vec<Z>) {
    let k = zi(k);
    let (a, b, cc, d, e, f) = (&c[0], &c[1], &c[2], &c[3], &c[4], &c[5]);
    (a * &k * &k + b * &k + cc, vec![zi(2) * a * &k + b, d * &k + e], vec![a.clone(), d.clone(), f.clone()])
}

pub fn dual_point(c1: &[Z], c2: &[Z], s_approx: [f64; 2]) -> R<DualPoint> {
    let a = conic(c1)?;
    let b = conic(c2)?;
    let s = s_approx;
    if !s.iter().all(|v| v.is_finite()) {
        bail!("dualPoint:approx", "sApprox must be a finite 2-vector.");
    }
    if all_zero(&a[0..3]) && all_zero(&b[0..3]) {
        bail!("dualPoint:twoLines", "both curves are lines, so their crossing is RATIONAL.");
    }
    for k in 0..13i64 {
        let (a2, a1, a0) = in_v(&a, k);
        let (b2, b1, b0) = in_v(&b, k);
        if a2.is_zero() && b2.is_zero() {
            continue;
        }
        let x = padd(&scal(&b0, &a2), &scal(&a0, &-&b2));
        let y = padd(&scal(&b1, &a2), &scal(&a1, &-&b2));
        let z = padd(&poly_mul1(&a1, &b0), &scal(&poly_mul1(&b1, &a0), &zi(-1)));
        let r = poly_trim(&padd(&poly_mul1(&x, &x), &scal(&poly_mul1(&y, &z), &zi(-1))));
        if r.is_empty() {
            bail!("dualPoint:commonComponent", "the conics share a component, so they meet in a curve and not in a point.");
        }
        let r = poly_primitive(&poly_squarefree(&r));
        let s1 = poly_trim(&padd(&scal(&a1, &b2), &scal(&b1, &-&a2)));
        let s0 = poly_trim(&padd(&scal(&a0, &b2), &scal(&b0, &-&a2)));
        if s1.is_empty() {
            continue;
        }
        let u = s[0] - (k as f64) * s[1];
        let iv = match isolate_near(&r, u) {
            None => bail!("dualPoint:noMatch", "no real root of the eliminant lies near u = {}.", u),
            Some(iv) => iv,
        };
        if isolated_root_sign(&s1, &r, &iv)? == 0 {
            continue;
        }
        let u0 = narrow(&r, &iv);
        let s0f: Vec<f64> = s0.iter().map(to_f64).collect();
        let s1f: Vec<f64> = s1.iter().map(to_f64).collect();
        let v0 = -polyval_f(&s0f, u0) / polyval_f(&s1f, u0);
        let sp = [u0 + (k as f64) * v0, v0];
        if (sp[0] - s[0]).hypot(sp[1] - s[1]) > 1e-6 * 1f64.max(s[0].hypot(s[1])) {
            bail!("dualPoint:noMatch", "no common root lies near the approximation.");
        }
        return Ok(DualPoint { k, r, iv, s1, s0, s: sp });
    }
    bail!("dualPoint:noShear", "no shear k = 0..12 gives the conics a second coordinate determined at the root.")
}

fn isolate_near(r: &[Z], u: f64) -> Option<(Q, Q)> {
    let sc = 1f64.max(u.abs());
    for e in -48..5 {
        let d = sc * 2f64.powi(e);
        let lo = q_from_f64(u - d);
        let hi = q_from_f64(u + d);
        if polyval_q(r, &lo).is_zero() || polyval_q(r, &hi).is_zero() {
            continue;
        }
        let n = count_roots_in_q(r, &lo, &hi);
        if n == 1 {
            return Some((lo, hi));
        }
        if n > 1 {
            break;
        }
    }
    None
}

fn narrow(r: &[Z], iv: &(Q, Q)) -> f64 {
    let mut lo = q_to_f64(&iv.0);
    let mut hi = q_to_f64(&iv.1);
    for _ in 0..200 {
        let mid = (lo + hi) / 2.0;
        if mid <= lo || mid >= hi {
            break;
        }
        if polyval_q(r, &q_from_f64(mid)).is_zero() {
            return mid;
        }
        if count_roots_in_q(r, &q_from_f64(lo), &q_from_f64(mid)) == 1 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (lo + hi) / 2.0
}

pub fn dual_point_sign(dp: &DualPoint, q: &M2) -> R<i32> {
    if biv::is_zero(q) {
        return Ok(0);
    }
    let k = zi(dp.k);
    let ni = q.len();
    let nj = q[0].len();
    let m0 = (ni - 1) + (nj - 1);
    let mut quv = vec![vec![Z::zero(); m0 + 1]; ni];
    for i in 0..ni {
        for j in 0..nj {
            if q[i][j].is_zero() {
                continue;
            }
            for t in 0..=i {
                quv[i - t][t + j] += &q[i][j] * Z::from(crate::poly_sign_on::comb(i, t)) * num_traits::pow(k.clone(), t);
            }
        }
    }
    while quv[0].len() > 1 && quv.iter().all(|r| r.last().unwrap().is_zero()) {
        for r in quv.iter_mut() {
            r.pop();
        }
    }
    let m = quv[0].len() - 1;
    let ppow = |a: &[Z], n: usize| -> Vec<Z> {
        let mut p = vec![Z::one()];
        for _ in 0..n {
            p = poly_mul1(&p, a);
        }
        p
    };
    let neg_s0 = neg(&dp.s0);
    let mut n: Vec<Z> = vec![];
    for j in 0..=m {
        let qj = poly_trim(&(0..ni).rev().map(|a| quv[a][j].clone()).collect::<Vec<_>>());
        if qj.is_empty() {
            continue;
        }
        n = padd(&n, &poly_mul1(&poly_mul1(&qj, &ppow(&neg_s0, j)), &ppow(&dp.s1, m - j)));
    }
    let sn = isolated_root_sign(&n, &dp.r, &dp.iv)?;
    if sn == 0 || m % 2 == 0 {
        return Ok(sn);
    }
    Ok(sn * isolated_root_sign(&dp.s1, &dp.r, &dp.iv)?)
}

// ---- cut_triangles ------------------------------------------------------------------------------
pub fn cut_triangles(sub: Vec<SubFace>, tfs: Vec<Rc<TriFace>>) -> R<Vec<SubFace>> {
    let mut sub: Vec<SubFace> = sub
        .into_iter()
        .map(|mut f| {
            f.tri = vec![];
            f.tf = None;
            f.tfs = tfs.clone();
            f
        })
        .collect();
    let mut tri_f = vec![];
    for (t0, tf) in tfs.iter().enumerate() {
        let t = (t0 + 1) as i64;
        let mut keep = vec![true; sub.len()];
        let mut out = vec![];
        let l = triangle_lines(tf)?;
        let mut keys = vec![];
        for x in &l {
            keys.push(conic_lin(x)?);
        }
        let flip: Vec<i64> = l.iter().map(|x| sign(x.iter().find(|v| !v.is_zero()).unwrap()) as i64).collect();
        for (k, f) in sub.iter().enumerate() {
            let mut sd = vec![];
            for i in 0..l.len() {
                let v = side_of(&f.con, &keys[i], flip[i], &l[i])?;
                sd.push(v);
                if v < 0 {
                    break;
                }
            }
            if sd.iter().any(|v| *v < 0) {
                continue;
            }
            if tf.rational {
                keep[k] = false;
                if sd.iter().all(|v| *v > 0) {
                    continue;
                }
                for e in 0..3 {
                    if sd[e] > 0 {
                        continue;
                    }
                    let mut p = f.clone();
                    for i in 0..=e {
                        let mut r = keys[i].clone();
                        r.push(zi(if i < e { flip[i] } else { -flip[e] }));
                        r.push(zi(0));
                        p.con.push(r);
                    }
                    out.push(p);
                }
                continue;
            }
            keep[k] = false;
            for e in 1..4i64 {
                let mut p = f.clone();
                for i in 1..e {
                    p.tri.push([t, i, 1]);
                }
                p.tri.push([t, e, -1]);
                out.push(p);
            }
        }
        tri_f.push(triangle_itself(tf, t, &keys, &flip, &tfs)?);
        sub = sub.into_iter().zip(keep).filter(|(_, k)| *k).map(|(f, _)| f).chain(out).collect();
    }
    sub.extend(tri_f);
    Ok(sub)
}

fn side_of(con: &[Vec<Z>], key: &[Z], flip: i64, line: &[Z]) -> R<i64> {
    for r in con {
        if &r[0..6] == key && (r.len() < 8 || r[7].is_zero()) {
            return Ok(side_i(&r[6]) * flip);
        }
    }
    let v = polygon_of(con);
    if !v.is_empty() {
        let lf: Vec<f64> = line.iter().map(to_f64).collect();
        let vals: Vec<f64> = v.iter().map(|p| p[0] * lf[0] + p[1] * lf[1] + lf[2]).collect();
        let m = 1e-9 * 1f64.max(lf.iter().map(|x| x * x).sum::<f64>().sqrt()) * 1f64.max(v.iter().flatten().fold(0.0f64, |a, c| a.max(c.abs())));
        if vals.iter().all(|x| *x > m) {
            return Ok(1);
        }
        if vals.iter().all(|x| *x < -m) {
            return Ok(-1);
        }
    }
    let kept: Vec<&Vec<Z>> = con.iter().filter(|r| !all_zero(&r[0..5])).collect();
    let mut c: Vec<Vec<Z>> = kept.iter().map(|r| {
        let mut v = r[0..6].to_vec();
        v.push(if r.len() >= 8 { r[7].clone() } else { zi(0) });
        v
    }).collect();
    let rows: Vec<(usize, i64)> = kept.iter().enumerate().map(|(i, r)| (i + 1, side_i(&r[6]))).collect();
    let mut kk = key.to_vec();
    kk.push(zi(0));
    c.push(kk);
    let j = c.len();
    let mut r1 = rows.clone();
    r1.push((j, -flip));
    if !cell_interior(&c, &r1)? {
        return Ok(1);
    }
    let mut r2 = rows;
    r2.push((j, flip));
    if !cell_interior(&c, &r2)? {
        return Ok(-1);
    }
    Ok(0)
}

fn polygon_of(con: &[Vec<Z>]) -> Vec<[f64; 2]> {
    let mut h: Vec<[f64; 3]> = con
        .iter()
        .filter(|r| all_zero(&r[0..3]) && (r.len() < 8 || r[7].is_zero()) && (!r[3].is_zero() || !r[4].is_zero()))
        .map(|r| [to_f64(&(&r[6] * &r[3])), to_f64(&(&r[6] * &r[4])), to_f64(&(&r[6] * &r[5]))])
        .collect();
    h.extend([[1.0, 0.0, 1e6], [-1.0, 0.0, 1e6], [0.0, 1.0, 1e6], [0.0, -1.0, 1e6]]);
    let mh = h.iter().flatten().fold(0.0f64, |a, x| a.max(x.abs()));
    let mut v = vec![];
    for i in 0..h.len() {
        for j in i + 1..h.len() {
            let m = [[h[i][0], h[i][1]], [h[j][0], h[j][1]]];
            let nm = m.iter().flatten().map(|x| x * x).sum::<f64>().sqrt();
            if (m[0][0] * m[1][1] - m[0][1] * m[1][0]).abs() < 1e-12 * 1f64.max(nm).powi(2) {
                continue;
            }
            let x = solve2f(&m, [-h[i][2], -h[j][2]]);
            let lim = -1e-9 * 1f64.max(x[0].hypot(x[1])) * 1f64.max(mh);
            if h.iter().all(|hh| hh[0] * x[0] + hh[1] * x[1] + hh[2] >= lim) {
                v.push(x);
            }
        }
    }
    v
}

fn triangle_lines(tf: &TriFace) -> R<Vec<Vec<Z>>> {
    if !tf.rational {
        return box_around(tf);
    }
    let mut l = vec![];
    for e in 0..3 {
        let i = e;
        let j = (e + 1) % 3;
        let o = 3 - i - j;
        let d = lcm(&tf.td[i], &tf.td[j]);
        let a: Vec<Z> = tf.tn[i].iter().map(|x| x * (&d / &tf.td[i])).collect();
        let b: Vec<Z> = tf.tn[j].iter().map(|x| x * (&d / &tf.td[j])).collect();
        let n = [-(&b[1] - &a[1]), &b[0] - &a[0]];
        let r = vec![&n[0] * &d, &n[1] * &d, -(&n[0] * &a[0] + &n[1] * &a[1])];
        let g = r[0].gcd(&r[1]).gcd(&r[2]);
        let mut r: Vec<Z> = r.iter().map(|x| x.div_floor(&g)).collect();
        if (&r[0] * &tf.tn[o][0] + &r[1] * &tf.tn[o][1] + &r[2] * &tf.td[o]).is_negative() {
            r = neg(&r);
        }
        l.push(r);
    }
    Ok(l)
}

fn box_around(tf: &TriFace) -> R<Vec<Vec<Z>>> {
    let t = &tf.t;
    let mn = [t.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min), t.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min)];
    let mx = [t.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max), t.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max)];
    let w = 1f64.max((mx[0] - mn[0]).max(mx[1] - mn[1]));
    let mut e = -8;
    while e <= 0 {
        let dl = 2f64.powi(w.log2().floor() as i32 + e);
        let lo = [((mn[0] - dl) / dl).floor() * dl, ((mn[1] - dl) / dl).floor() * dl];
        let hi = [((mx[0] + dl) / dl).ceil() * dl, ((mx[1] + dl) / dl).ceil() * dl];
        if all_inside(tf, lo, hi)? {
            let mut l = vec![];
            for k in 0..2 {
                let fl = q_from_f64(lo[k]);
                let fh = q_from_f64(hi[k]);
                let mut r1 = vec![zi(0), zi(0), -fl.numer().clone()];
                r1[k] = fl.denom().clone();
                let mut r2 = vec![zi(0), zi(0), fh.numer().clone()];
                r2[k] = -fh.denom().clone();
                l.push(r1);
                l.push(r2);
            }
            return Ok(l);
        }
        e += 2;
    }
    bail!("cutTriangles:box", "no dyadic box up to the triangle's own size provably contains it.")
}

fn all_inside(tf: &TriFace, lo: [f64; 2], hi: [f64; 2]) -> R<bool> {
    for i in 0..3 {
        let g = &tf.gn[i];
        let m = contact_poly(g);
        for k in 0..2 {
            let fl = q_from_f64(lo[k]);
            let fh = q_from_f64(hi[k]);
            let mut a = biv::scal(&m[k], fl.denom());
            a[0][0] -= fl.numer() * &tf.gd[i];
            let mut b = biv::scal(&m[k], &-fh.denom());
            b[0][0] += fh.numer() * &tf.gd[i];
            if triangle_face_sign(tf, &a)? < 0 || triangle_face_sign(tf, &b)? < 0 {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn triangle_itself(tf: &Rc<TriFace>, t: i64, keys: &[Vec<Z>], flip: &[i64], tfs: &[Rc<TriFace>]) -> R<SubFace> {
    let con: Vec<Vec<Z>> = keys.iter().zip(flip.iter()).map(|(k, f)| {
        let mut r = k.clone();
        r.push(zi(*f));
        r.push(zi(0));
        r
    }).collect();
    if tf.rational {
        let (num, den) = plane_q(tf)?;
        return Ok(SubFace { kind: SubKind::Quad, num, den, fp: None, fb: 0, con, alg: vec![], ii: vec![], tri: vec![], tf: None, tfs: tfs.to_vec() });
    }
    Ok(SubFace { kind: SubKind::Triangle, num: vec![Z::zero(); 10], den: Z::one(), fp: None, fb: 0, con, alg: vec![], ii: vec![], tri: vec![[t, 1, 1], [t, 2, 1], [t, 3, 1]], tf: Some(tf.clone()), tfs: tfs.to_vec() })
}

fn plane_q(tf: &TriFace) -> R<(Vec<Z>, Z)> {
    let sn = &tf.s0n;
    let sd = &tf.s0d;
    let g = &tf.gn[0];
    let gd = &tf.gd[0];
    let c = &g[0] * &sn[0] * &sn[0] + &g[1] * &sn[0] * &sn[1] + &g[2] * &sn[1] * &sn[1] + (&g[3] * &sn[0] + &g[4] * &sn[1]) * sd + &g[5] * sd * sd;
    let mut num = vec![Z::zero(); 10];
    num[7] = &sn[0] * sd * gd;
    num[8] = &sn[1] * sd * gd;
    num[9] = -c;
    canon(&num, &(sd * sd * gd))
}

//! alg_branch -- evaluate a NAMED real root of an implicit face polynomial (port of
//! alg_branch.py). b > 0 counts up from the smallest real root, b < 0 down from the largest.

use crate::bail;
use crate::errors::R;
use crate::numeric::{roots, Cx};
use crate::poly_n::PolyN;
use crate::rat_q::*;
use num_traits::{One, Signed, Zero};

const REALMIN: f64 = 2.2250738585072014e-308;

/// The polynomial prepared for evaluation: per z-power, its (i, j, c) terms.
pub struct Prep {
    lay: Vec<Vec<(usize, usize, Z)>>,
    lay_f: Vec<Vec<(usize, usize, f64)>>,
    nx: usize,
    ny: usize,
    nz: usize,
    big: bool,
}

pub fn prep(p: &PolyN) -> R<Prep> {
    if p.is_empty() {
        bail!("algBranch:polynomial", "P must be a non-empty polynomial.");
    }
    let nx = p.keys().map(|e| e[0]).max().unwrap() as usize + 1;
    let ny = p.keys().map(|e| e[1]).max().unwrap() as usize + 1;
    let nz = p.keys().map(|e| e[2]).max().unwrap() as usize + 1;
    let mut lay: Vec<Vec<(usize, usize, Z)>> = vec![vec![]; nz];
    for (e, c) in p {
        lay[e[2] as usize].push((e[0] as usize, e[1] as usize, c.clone()));
    }
    let lay_f = lay.iter().map(|t| t.iter().map(|(i, j, c)| (*i, *j, to_f64(c))).collect()).collect();
    let big = p.values().map(|c| c.abs()).max().unwrap() >= (Z::one() << 53usize);
    Ok(Prep { lay, lay_f, nx, ny, nz, big })
}

pub fn alg_branch(p: &PolyN, x: &[[f64; 2]], b: i64) -> R<Vec<f64>> {
    let pr = prep(p)?;
    alg_branch_prep(&pr, x, b)
}

pub fn alg_branch_prep(pr: &Prep, x: &[[f64; 2]], b: i64) -> R<Vec<f64>> {
    if b == 0 {
        bail!("algBranch:branch", "the branch index must be a nonzero integer.");
    }
    if pr.nz < 2 {
        bail!("algBranch:noZ", "P has degree 0 in z, so {{p = 0}} does not define a value.");
    }
    if pr.nz <= 3 {
        return Ok(x.iter().map(|p| fibre_low(pr, p[0], p[1], b)).collect());
    }
    let (z, _) = general(pr, x, b);
    Ok(z)
}

/// (z, allReal): the wantInfo form.
pub fn alg_branch_info(p: &PolyN, x: &[[f64; 2]], b: i64) -> R<(Vec<f64>, Vec<Vec<f64>>)> {
    let pr = prep(p)?;
    if pr.nz < 2 {
        bail!("algBranch:noZ", "P has degree 0 in z.");
    }
    Ok(general(&pr, x, b))
}

fn fpow(x: f64, i: usize) -> f64 {
    x.powf(i as f64)
}

fn fibre_double(pr: &Prep, x: f64, y: f64, deg_tol: f64) -> (Vec<f64>, Vec<bool>) {
    let mut c = vec![];
    let mut is_zero = vec![];
    for terms in &pr.lay {
        let mut v = 0.0;
        let mut a = 0.0;
        for (i, j, co) in terms {
            let t = to_f64(co) * fpow(x, *i) * fpow(y, *j);
            v += t;
            a += t.abs();
        }
        c.push(v);
        is_zero.push(v.abs() <= deg_tol * a);
    }
    (c, is_zero)
}

fn fibre_exact(pr: &Prep, x: f64, y: f64) -> (Vec<f64>, Vec<bool>) {
    let qx = q_from_f64(x);
    let qy = q_from_f64(y);
    let (a, dx) = (qx.numer().clone(), qx.denom().clone());
    let (cc, dy) = (qy.numer().clone(), qy.denom().clone());
    let nx = 1 + pr.lay.iter().flatten().map(|t| t.0).max().unwrap_or(0);
    let ny = 1 + pr.lay.iter().flatten().map(|t| t.1).max().unwrap_or(0);
    let av: Vec<Z> = (0..nx).map(|i| num_traits::pow(a.clone(), i) * num_traits::pow(dx.clone(), nx - 1 - i)).collect();
    let cv: Vec<Z> = (0..ny).map(|j| num_traits::pow(cc.clone(), j) * num_traits::pow(dy.clone(), ny - 1 - j)).collect();
    let d = num_traits::pow(dx, nx - 1) * num_traits::pow(dy, ny - 1);
    let mut c = vec![];
    let mut is_zero = vec![];
    for terms in &pr.lay {
        let mut v = Z::zero();
        for (i, j, co) in terms {
            v += co * &av[*i] * &cv[*j];
        }
        is_zero.push(v.is_zero());
        c.push(q_to_f64(&Q::new(v, d.clone())));
    }
    (c, is_zero)
}

fn abs_fibre(pr: &Prep, upto: usize, x: f64, y: f64) -> Vec<f64> {
    pr.lay_f[..upto.min(pr.lay_f.len())]
        .iter()
        .map(|terms| terms.iter().map(|(i, j, co)| co.abs() * fpow(x.abs(), *i) * fpow(y.abs(), *j)).sum())
        .collect()
}

fn exact_quad(pr: &Prep, x: f64, y: f64) -> ([f64; 3], f64) {
    let xq = q_from_f64(x);
    let yq = q_from_f64(y);
    let mut acc = vec![];
    for k in 0..3 {
        let mut v = Q::zero();
        if k < pr.lay.len() {
            for (i, j, co) in &pr.lay[k] {
                v += qz(co) * pow_q(&xq, *i) * pow_q(&yq, *j);
            }
        }
        acc.push(v);
    }
    let dx = &acc[1] * &acc[1] - qi(4) * &acc[2] * &acc[0];
    let c = [q_to_f64(&acc[0]), q_to_f64(&acc[1]), q_to_f64(&acc[2])];
    let mut d = q_to_f64(&dx);
    if d == 0.0 && !dx.is_zero() {
        d = if dx.is_positive() { REALMIN } else { -REALMIN };
    }
    (c, d)
}

fn pow_q(x: &Q, i: usize) -> Q {
    let mut r = qi(1);
    for _ in 0..i {
        r *= x;
    }
    r
}

fn double_root_noise(cz: &[f64], v: Cx) -> bool {
    let t = v.re;
    let cabs: f64 = cz.iter().enumerate().map(|(k, c)| c.abs() * fpow(t.abs(), k)).sum();
    let mut p2 = 0.0;
    for (k, c) in cz.iter().enumerate() {
        if k >= 2 {
            p2 += (k * (k - 1)) as f64 / 2.0 * c * fpow(t, k - 2);
        }
    }
    if p2 == 0.0 || cabs == 0.0 {
        return false;
    }
    v.im.abs() <= 8.0 * (2.220446049250313e-16 * cabs / p2.abs()).sqrt()
}

const PROBES: [(i64, i64); 4] = [(1, 0), (0, 1), (1, 1), (1, -1)];

/// Decide EXACTLY the near-double roots the doubles cannot (#108 residual, D-153); see
/// alg_branch.py `_nearDoubles`. A pair within the double-root noise of each other (two close
/// reals, or a complex pair `double_root_noise` accepted) is dropped when it is real only on a
/// CURVE: the fibre is evaluated exactly at x +- h e for four directions e and a tiny dyadic h, and
/// a pair with no real root near it on EITHER side in at least two directions is not the value of
/// any branch (f** is closed: its value is the limit from around). Such a pair is where a bridge's
/// dual point runs off to infinity, on a boundary line of dom f**.
fn near_doubles(pr: &Prep, x: f64, y: f64, cz: &[f64], r: &[Cx], rr: Vec<f64>, imag_tol: f64) -> Vec<f64> {
    let mut cands = vec![];
    for v in r {
        if v.im > 0.0 && v.im.abs() > imag_tol * 1f64.max(cabs(*v)) && double_root_noise(cz, *v) {
            cands.push(v.re);
        }
    }
    for w in rr.windows(2) {
        let (a, b) = (w[0], w[1]);
        if b > a && double_root_noise(cz, Cx { re: (a + b) / 2.0, im: (b - a) / 2.0 }) {
            cands.push((a + b) / 2.0);
        }
    }
    if cands.is_empty() {
        return rr;
    }
    let mut out = rr;
    for z0 in cands {
        let mut d: Vec<f64> = r.iter().map(|v| cabs(Cx { re: v.re - z0, im: v.im })).collect();
        d.sort_by(|a, b| a.total_cmp(b));
        let rad = if d.len() > 2 && d[2] > 0.0 { d[2] / 2.0 } else { 1.0 + z0.abs() };
        if isolated_pair(pr, x, y, z0, rad) {
            for _ in 0..2 {
                let k = (0..out.len()).min_by(|&i, &j| (out[i] - z0).abs().total_cmp(&(out[j] - z0).abs()));
                if let Some(k) = k {
                    if (out[k] - z0).abs() <= rad {
                        out.remove(k);
                    }
                }
            }
        }
    }
    snap_multiple_roots(pr, x, y, cz, out)
}

/// A GENUINE multiple root of the fibre, to full precision (#108 residual: the jump vertices);
/// see alg_branch.py `_snapMultipleRoots`. Rounding returns a k-fold root as k roots ~eps^(1/k)
/// apart (D-153's vertex (-3,-2): 7.5000017 for 7.5; the kink's (1,2): 1.0000000256 for 1). The
/// fibre at the exact dyadic point is split by multiplicity (Yun); each real root of a factor of
/// multiplicity k >= 2 is isolated by Sturm and bisected exactly, and whatever rounding left
/// within the noise of a k-fold root becomes exactly k copies of it (nothing there: stays out).
fn snap_multiple_roots(pr: &Prep, x: f64, y: f64, cz: &[f64], rr: Vec<f64>) -> Vec<f64> {
    if rr.is_empty() {
        return rr;
    }
    let p = poly_primitive(&exact_fibre_at(pr, &q_from_f64(x), &q_from_f64(y)));
    if p.len() <= 2 {
        return rr;
    }
    let mut out = rr;
    let fac = match yun(&p) {
        Some(f) => f,
        None => return out,
    };
    for (k, a) in fac {
        if k < 2 {
            continue;
        }
        let ivs = match isolate_roots(&a) {
            Ok(v) => v,
            Err(_) => continue,
        };
        for (lo, hi) in ivs {
            let m = bisect_root(&a, lo, hi);
            let w = multiple_root_noise(cz, m, k);
            let keep: Vec<f64> = out.iter().cloned().filter(|v| (v - m).abs() > w).collect();
            if keep.len() < out.len() {
                // rounding left something there: exactly k copies of m, whatever the LAPACK
                out = keep;
                out.extend(std::iter::repeat(m).take(k));
            }
        }
    }
    out.sort_by(|a, b| a.total_cmp(b));
    out
}

/// 8 (eps sum|c_j m^j| / |p^(k)(m)/k!|)^(1/k): how far rounding moves a k-fold root m.
fn multiple_root_noise(cz: &[f64], m: f64, k: usize) -> f64 {
    let cabs: f64 = cz.iter().enumerate().map(|(j, c)| c.abs() * fpow(m.abs(), j)).sum();
    let mut pk = 0.0;
    for (j, c) in cz.iter().enumerate() {
        if j >= k {
            pk += binom(j, k) * c * fpow(m, j - k);
        }
    }
    if pk == 0.0 || cabs == 0.0 {
        return 0.0;
    }
    8.0 * (2.220446049250313e-16 * cabs / pk.abs()).powf(1.0 / k as f64)
}

fn binom(n: usize, k: usize) -> f64 {
    let mut r = 1.0;
    for i in 0..k {
        r = r * (n - i) as f64 / (i + 1) as f64;
    }
    r
}

/// Yun's square-free decomposition: (k, a_k) with p = c prod a_k^k, a_k square-free, primitive,
/// of positive degree. None if an exact division unexpectedly fails.
fn yun(p: &[Z]) -> Option<Vec<(usize, Vec<Z>)>> {
    let g = poly_gcd(p, &polyder(p)).ok()?;
    if g.len() <= 1 {
        return Some(vec![(1, p.to_vec())]);
    }
    let mut out = vec![];
    let mut b = poly_div_exact1(&poly_primitive(p), &g)?;
    let mut k = 1;
    let mut rest = g;
    while b.len() > 1 {
        let a = poly_gcd(&b, &rest).ok()?;
        let f = if a.len() > 1 { poly_div_exact1(&b, &a)? } else { b.clone() };
        if f.len() > 1 {
            out.push((k, poly_primitive(&f)));
        }
        if a.len() <= 1 {
            break;
        }
        rest = poly_div_exact1(&rest, &a)?;
        b = a;
        k += 1;
    }
    Some(out)
}

/// The simple root of the square-free a in (lo, hi], bisected exactly until the interval is below
/// 2^-60 of its size (or 400 halvings); the double nearest the midpoint.
fn bisect_root(a: &[Z], lo: Q, hi: Q) -> f64 {
    let (mut lo, mut hi) = (lo, hi);
    if a.last().map_or(false, |c| c.is_zero()) && lo < Q::zero() && hi >= Q::zero() {
        return 0.0; // z = 0 is the root (no relative size to reach)
    }
    let tiny = Q::new(Z::one(), Z::one() << 60usize);
    for _ in 0..400 {
        if &hi - &lo <= (&hi + &lo).abs() * &tiny {
            break;
        }
        let mid = (&lo + &hi) / qi(2);
        if count_roots_in_q(a, &lo, &mid) == 1 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    q_to_f64(&((&lo + &hi) / qi(2)))
}

/// The fibre at the exact rational point (xq, yq) as an integer polynomial in z, highest power
/// first (scaled by a positive common denominator: the same roots).
fn exact_fibre_at(pr: &Prep, xq: &Q, yq: &Q) -> Vec<Z> {
    let mut c: Vec<Q> = vec![];
    for terms in &pr.lay {
        let mut v = Q::zero();
        for (i, j, co) in terms {
            v += qz(co) * pow_q(xq, *i) * pow_q(yq, *j);
        }
        c.push(v);
    }
    let mut den = Z::one();
    for v in &c {
        den = num_integer::Integer::lcm(&den, v.denom());
    }
    c.iter().rev().map(|v| (v * qz(&den)).to_integer()).collect()
}

fn isolated_pair(pr: &Prep, x: f64, y: f64, z0: f64, rad: f64) -> bool {
    let xq = q_from_f64(x);
    let yq = q_from_f64(y);
    let e = libm_frexp_exp(1f64.max(x.abs()).max(y.abs()));
    let h = if e >= 40 { qz(&(Z::one() << ((e - 40) as usize))) } else { Q::new(Z::one(), Z::one() << ((40 - e) as usize)) };
    let lo = q_from_f64(z0 - rad);
    let hi = q_from_f64(z0 + rad);
    let mut both = 0;
    for (ex, ey) in PROBES {
        let mut none = true;
        for sg in [1i64, -1] {
            let px = &xq + &h * qi(sg * ex);
            let py = &yq + &h * qi(sg * ey);
            if count_roots_in_q(&exact_fibre_at(pr, &px, &py), &lo, &hi) != 0 {
                none = false;
                break;
            }
        }
        if none {
            both += 1;
        }
    }
    both >= 2
}

/// The exponent e of math.frexp(v) = (m, e), v = m 2^e with 0.5 <= |m| < 1, for v >= 1 finite.
fn libm_frexp_exp(v: f64) -> i64 {
    let bits = v.to_bits();
    let exp = ((bits >> 52) & 0x7ff) as i64;
    exp - 1022
}

/// Is p(x, y, .) identically zero at each point -- a BASE POINT, where every z is a root and a
/// row comparing this value is undefined? Every coefficient zero by alg_branch's own zero test.
pub fn fibre_vanishes(p: &PolyN, x: &[[f64; 2]]) -> R<Vec<bool>> {
    let pr = prep(p)?;
    Ok(x.iter()
        .map(|q| {
            let (_, z) = if pr.big { fibre_exact(&pr, q[0], q[1]) } else { fibre_double(&pr, q[0], q[1], 1e-12) };
            z.iter().all(|b| *b)
        })
        .collect())
}

fn sgnz(v: f64) -> f64 {
    if v >= 0.0 {
        1.0
    } else {
        -1.0
    }
}

fn cabs(z: Cx) -> f64 {
    z.re.hypot(z.im)
}

fn general(pr: &Prep, x: &[[f64; 2]], b: i64) -> (Vec<f64>, Vec<Vec<f64>>) {
    let imag_tol = 1e-10;
    let deg_tol = 1e-12;
    let mut zv = vec![];
    let mut all = vec![];
    for p in x {
        let (xx, yy) = (p[0], p[1]);
        let (c, is_zero) = if pr.big { fibre_exact(pr, xx, yy) } else { fibre_double(pr, xx, yy, deg_tol) };
        let mut zt = f64::NAN;
        let mut rr: Vec<f64> = vec![];
        let scale = c.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let keep = (0..pr.nz).filter(|k| !is_zero[*k]).max();
        let cz: Vec<f64> = match keep {
            None => vec![],
            Some(k) => c[..=k].to_vec(),
        };
        if scale != 0.0 && cz.len() >= 2 {
            let rev: Vec<f64> = cz.iter().rev().cloned().collect();
            let r = roots(&rev);
            rr = r.iter().filter(|v| v.im.abs() <= imag_tol * 1f64.max(cabs(**v)) || double_root_noise(&cz, **v)).map(|v| v.re).collect();
            rr.sort_by(|a, b| a.total_cmp(b)); // NaN-safe (review #6)
            rr = near_doubles(pr, xx, yy, &cz, &r, rr, imag_tol);
            if cz.len() == 3 {
                let ca = abs_fibre(pr, 3, xx, yy);
                let dq = cz[1].powi(2) - 4.0 * cz[2] * cz[0];
                let t = ca[1].powi(2) + 4.0 * ca[2] * ca[0];
                if dq.abs() <= 1e-12 * t {
                    let (cx, dx) = exact_quad(pr, xx, yy);
                    if dx > 0.0 {
                        let qq = -(cx[1] + sgnz(cx[1]) * dx.sqrt()) / 2.0;
                        let mut v = vec![qq / cx[2], cx[0] / qq];
                        v.sort_by(|a, b| a.total_cmp(b)); // NaN-safe (review #6)
                        rr = v;
                    } else {
                        rr = vec![-cz[1] / (2.0 * cz[2]); 2];
                    }
                }
            }
            let idx = if b > 0 { b } else { rr.len() as i64 + 1 + b };
            if idx >= 1 && idx as usize <= rr.len() {
                zt = rr[idx as usize - 1];
            }
        }
        zv.push(zt);
        all.push(rr);
    }
    (zv, all)
}

fn fibre_low(pr: &Prep, x: f64, y: f64, b: i64) -> f64 {
    let imag_tol = 1e-10;
    let deg_tol = 1e-12;
    let xp: Vec<f64> = (0..pr.nx).map(|i| fpow(x, i)).collect();
    let yp: Vec<f64> = (0..pr.ny).map(|j| fpow(y, j)).collect();
    let axp: Vec<f64> = (0..pr.nx).map(|i| fpow(x.abs(), i)).collect();
    let ayp: Vec<f64> = (0..pr.ny).map(|j| fpow(y.abs(), j)).collect();
    let mut ca: Vec<f64> = pr.lay_f.iter().map(|t| t.iter().map(|(i, j, co)| co.abs() * axp[*i] * ayp[*j]).sum()).collect();
    ca.resize(3, 0.0);
    let (mut c, mut is_zero) = if pr.big {
        fibre_exact(pr, x, y)
    } else {
        let mut c = vec![];
        let mut z = vec![];
        for terms in &pr.lay_f {
            let mut v = 0.0;
            let mut a = 0.0;
            for (i, j, co) in terms {
                let t = co * xp[*i] * yp[*j];
                v += t;
                a += t.abs();
            }
            c.push(v);
            z.push(v.abs() <= deg_tol * a);
        }
        (c, z)
    };
    c.resize(3, 0.0);
    is_zero.resize(3, true);
    let lead = if !is_zero[2] { 3 } else if !is_zero[1] { 2 } else { 0 };
    if lead == 2 {
        if c.iter().fold(0.0f64, |m, v| m.max(v.abs())) > 0.0 && (b == 1 || b == -1) {
            return -c[0] / c[1];
        }
        return f64::NAN;
    }
    if lead != 3 {
        return f64::NAN;
    }
    let (mut a, mut bb, mut c0) = (c[2], c[1], c[0]);
    let mut dq = bb * bb - 4.0 * a * c0;
    let t = ca[1].powi(2) + 4.0 * ca[2] * ca[0];
    let win = dq.abs() <= 1e-12 * t;
    if win {
        let (cx, dx) = exact_quad(pr, x, y);
        if dx > 0.0 {
            a = cx[2];
            bb = cx[1];
            c0 = cx[0];
            dq = dx;
        } else {
            dq = 0.0;
        }
    }
    let mean2 = -bb / (2.0 * a);
    let (mut lo, mut hi) = (f64::NAN, f64::NAN);
    if win && dq <= 0.0 {
        lo = mean2;
        hi = mean2;
    } else if dq > 0.0 {
        let qq = -(bb + sgnz(bb) * dq.sqrt()) / 2.0;
        let r1 = qq / a;
        let r2 = if qq != 0.0 { c0 / qq } else { f64::NAN };
        // Python min/max with a NaN second argument return the first one
        if r2.is_nan() {
            lo = r1;
            hi = r1;
        } else {
            lo = r1.min(r2);
            hi = r1.max(r2);
        }
    } else if dq < 0.0 {
        let im = (-dq).sqrt() / (2.0 * a.abs());
        if im <= imag_tol * 1f64.max(mean2.hypot(im)) {
            lo = mean2;
            hi = mean2;
        }
    }
    match b {
        1 | -2 => lo,
        2 | -1 => hi,
        _ => f64::NAN,
    }
}

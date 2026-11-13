//! alg_row_fails -- a CERTIFICATE that a set of algebraic rows cannot all hold anywhere in the
//! part of a box that lies in a region (port of alg_row_fails.py). True is a proof; false means
//! "not proved". Square-free parts by Yun, through FLINT when available (src/flint.rs).

use crate::biv::{self, M2};
use crate::errors::R;
use crate::numeric::{roots, Cx, EPS};
use crate::poly_n::PolyN;
use crate::poly_sign_on::{box_sign, comb, FMat};
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::cell::RefCell;
use std::collections::HashMap;
use std::rc::Rc;

const REALMIN: f64 = 2.2250738585072014e-308;

/// A row: "side * (q - z) >= 0" (Quad) or "z_k <= z_w" (Face).
#[derive(Clone, Debug)]
pub enum RowSpec {
    Quad { fp: Rc<PolyN>, fb: i64, num: Vec<Z>, den: Z, side: i64 },
    Face { fp: Rc<PolyN>, fb: i64, wfp: Rc<PolyN>, wfb: i64 },
}

/// The region of a cell: its box, its row polynomials (also the factors), lines and corners.
pub struct Region {
    pub bx: [f64; 4],
    row_f: Vec<FMat>,
    pub factors: Vec<M2>,
    lines: Vec<(PairEval, f64)>,
    pub corners: Vec<[f64; 2]>,
}

impl Region {
    pub fn new(bx: [f64; 4], row_polys: Vec<M2>, lines: Vec<(Vec<Z>, i64, i64)>, corners: Vec<[f64; 2]>) -> R<Region> {
        let row_f = row_polys.iter().map(|m| fmat(m)).collect();
        let mut ls = vec![];
        for (c, side, br) in lines {
            ls.push((pair_eval_prep(&c, br)?, side as f64));
        }
        Ok(Region { bx, row_f, factors: row_polys, lines: ls, corners })
    }

    pub fn in_box(&self, b: &[f64; 4]) -> bool {
        for p in &self.row_f {
            if box_sign(p, b) == -1 {
                return false;
            }
        }
        let x = [[b[0], b[2]], [b[1], b[2]], [b[0], b[3]], [b[1], b[3]]];
        for (pe, s) in &self.lines {
            let v: Vec<f64> = x.iter().map(|q| s * pe.eval_m(q[0], q[1]).0).collect();
            let mx = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let sc = v.iter().fold(1f64, |m, x| m.max(x.abs()));
            if mx < -1e-7 * sc {
                return false;
            }
        }
        true
    }
}

pub fn fmat(m: &M2) -> FMat {
    m.iter().map(|r| r.iter().map(to_f64).collect()).collect()
}

// ---- sign polynomials ----------------------------------------------------------------------------
#[derive(Clone, Debug)]
#[allow(dead_code)] // the integer forms are kept as in Python's dict, for inspection
struct Surd {
    u: M2,
    v: M2,
    e: Option<M2>,
    us: Option<M2>,
    vs: Option<M2>,
    lo: FMat,
    hi: FMat,
    uz: Option<M2>,
    vz: Option<M2>,
    ez: Option<M2>,
    usz: Option<M2>,
    vsz: Option<M2>,
    k: Z,
    // float caches
    ef: Option<FMat>,
    usf: Option<FMat>,
    vsf: Option<FMat>,
    uf: FMat,
    vf: FMat,
}

#[derive(Clone, Debug)]
enum SP {
    Num { f: FMat, z: Option<M2> },
    Surd(Box<Surd>),
}

fn num_z(m: M2) -> SP {
    SP::Num { f: fmat(&m), z: Some(m) }
}

#[derive(Clone, Debug)]
enum Test {
    Top { want: i32 },
    Quad { hi: bool, want: i32 },
    Lin { want: i32 },
    Fq { hik: bool, hiw: bool },
    Ql { hik: bool },
    Lq { hiw: bool },
    Surd,
    TopFace { odd: bool },
}

struct Row {
    p: Vec<SP>,
    test: Test,
    corner: bool,
}

fn run_test(t: &Test, s: &[i32]) -> bool {
    match t {
        Test::Top { want } => s[0] != 0 && s[1..].iter().all(|v| s[0] * v == *want),
        Test::Quad { hi, want } => quad_ok(s, *hi, *want),
        Test::Lin { want } => s[0] != 0 && s[0] * s[1] == *want,
        Test::Fq { hik, hiw } => face_ok(s, *hik, *hiw),
        Test::Ql { hik } => s[1] != 0 && quad_ok(&[s[0], s[1] * s[2], s[3], s[4]], *hik, 1),
        Test::Lq { hiw } => s[1] != 0 && quad_ok(&[s[0], s[1] * s[2], s[3], s[4]], *hiw, -1),
        Test::Surd => surd_ok(s),
        Test::TopFace { odd } => top_ok(s, *odd),
    }
}

// ---- memos ----------------------------------------------------------------------------------------
thread_local! {
    static MEMO_SF: RefCell<HashMap<(PolyN, Vec<Z>, Z, i32, bool), Option<Vec<M2>>>> = RefCell::new(HashMap::new());
    static MEMO_FACE: RefCell<HashMap<(u8, PolyN, PolyN), Option<Vec<M2>>>> = RefCell::new(HashMap::new());
    static MEMO_SURD: RefCell<HashMap<(PolyN, PolyN, bool, bool), Option<(Vec<M2>, Vec<Surd>)>>> = RefCell::new(HashMap::new());
    static MEMO_TOP: RefCell<HashMap<(PolyN, PolyN, bool), Option<(Vec<M2>, Vec<Surd>, bool)>>> = RefCell::new(HashMap::new());
    static MEMO_SURDFACE: RefCell<HashMap<(PolyN, bool), Option<SurdFace>>> = RefCell::new(HashMap::new());
    static MEMO_SURDOF: RefCell<HashMap<(u64, u64), Option<[i64; 6]>>> = RefCell::new(HashMap::new());
    static MEMO_FACTOR_INTS: RefCell<HashMap<M2, Option<M2>>> = RefCell::new(HashMap::new());
    static MEMO_FACTOR_BU: RefCell<HashMap<(M2, u64, u64, i64, u8), Option<SP>>> = RefCell::new(HashMap::new());
    static MEMO_CURVE_PTS: RefCell<HashMap<M2, Rc<(bool, Option<Vec<(f64, Vec<Cx>)>>)>>> = RefCell::new(HashMap::new());
}

pub fn clear_memos() {
    MEMO_SURDOF.with(|m| m.borrow_mut().clear());
    MEMO_FACTOR_INTS.with(|m| m.borrow_mut().clear());
    MEMO_CURVE_PTS.with(|m| m.borrow_mut().clear());
}

// ---- the entry point ------------------------------------------------------------------------------
pub fn alg_row_fails(a: &[RowSpec], reg: &Region, depth: usize) -> R<bool> {
    alg_row_fails_box(a, reg.bx, &|b: &[f64; 4]| reg.in_box(b), depth, &reg.factors, &reg.corners)
}

pub fn alg_row_fails_box(a: &[RowSpec], bx: [f64; 4], in_box: &dyn Fn(&[f64; 4]) -> bool, depth: usize, factors: &[M2], corners: &[[f64; 2]]) -> R<bool> {
    let mut rows = vec![];
    for at in a {
        if let Some(r) = row_polys(at, factors)? {
            rows.push(r);
        }
    }
    if rows.is_empty() {
        return Ok(false);
    }
    let tf = box_rec(&rows, bx, in_box, depth, &[]);
    if tf || corners.is_empty() {
        return Ok(tf);
    }
    let w = (bx[1] - bx[0]).max(bx[3] - bx[2]);
    let mut sq: Vec<[f64; 3]> = vec![];
    let mut bu_memo: HashMap<(usize, usize, u64, u64, i64, u8), Option<SP>> = HashMap::new();
    for p0 in corners {
        let mut p = *p0;
        let q = [round(p[0] * 256.0) / 256.0, round(p[1] * 256.0) / 256.0];
        if q.iter().zip(p.iter()).all(|(qi, pi)| (qi - pi).abs() <= 1e-9 * 1f64.max(pi.abs())) {
            p = q;
        }
        if p[0] < bx[0] || p[0] > bx[1] || p[1] < bx[2] || p[1] > bx[3] {
            continue;
        }
        if !any_vanishes(&rows, p) {
            continue;
        }
        for k in 2..9 {
            let mut e = w / 2f64.powi(k);
            if corner_ok(&rows, p, e, factors, &mut bu_memo)? {
                if !is_dyadic(p) {
                    e *= 1.0 - 1e-9;
                }
                sq.push([p[0], p[1], e]);
                break;
            }
        }
    }
    if sq.is_empty() {
        return Ok(false);
    }
    let mut tf = box_rec(&rows, bx, in_box, depth, &sq);
    if !tf {
        tf = box_rec(&rows, bx, in_box, depth + 2, &sq);
    }
    Ok(tf)
}

fn round(x: f64) -> f64 {
    (x.abs() + 0.5).floor() * if x >= 0.0 { 1.0 } else { -1.0 }
}

fn slices(fp: &PolyN) -> Vec<M2> {
    let nx = fp.keys().map(|e| e[0]).max().unwrap() as usize + 1;
    let ny = fp.keys().map(|e| e[1]).max().unwrap() as usize + 1;
    let nz = fp.keys().map(|e| e[2]).max().unwrap() as usize + 1;
    let mut s = vec![biv::zeros(nx, ny); nz];
    for (e, c) in fp {
        s[e[2] as usize][e[0] as usize][e[1] as usize] = c.clone();
    }
    s
}

fn nz_of(fp: &PolyN) -> usize {
    fp.keys().map(|e| e[2]).max().unwrap() as usize + 1
}

fn row_polys(a: &RowSpec, factors: &[M2]) -> R<Option<Row>> {
    let (fp, fb, num, den, side) = match a {
        RowSpec::Face { fp, fb, wfp, wfb } => return face_row_polys(fp, *fb, wfp, *wfb, factors),
        RowSpec::Quad { fp, fb, num, den, side } => (fp, *fb, num, den, *side),
    };
    if num.len() != 10 || !all_zero(&num[0..4]) || !den.is_positive() {
        return Ok(None);
    }
    let q = vec![vec![zi(2) * &num[9], zi(2) * &num[8], num[6].clone()], vec![zi(2) * &num[7], zi(2) * &num[5], zi(0)], vec![num[4].clone(), zi(0), zi(0)]];
    let dn = zi(2) * den;
    let want = -side as i32;
    let nz = nz_of(fp);
    let (test, top) = if nz - 1 >= 3 {
        if fb != -1 {
            return Ok(None);
        }
        (Test::Top { want }, true)
    } else if nz - 1 == 2 {
        let hi = fb == -1 || fb == 2;
        if !hi && fb != 1 && fb != -2 {
            return Ok(None);
        }
        (Test::Quad { hi, want }, false)
    } else if nz - 1 == 1 {
        if fb != 1 && fb != -1 {
            return Ok(None);
        }
        (Test::Lin { want }, false)
    } else {
        return Ok(None);
    };
    let key = ((**fp).clone(), num.clone(), den.clone(), if top { want } else { 0 }, top);
    let p = match MEMO_SF.with(|m| m.borrow().get(&key).cloned()) {
        Some(p) => p,
        None => {
            let p = if top {
                let v = top_row_polys(fp, &q, &dn, nz - 1, want);
                if v.is_empty() {
                    None
                } else {
                    Some(v)
                }
            } else {
                Some(sign_polys(&slices(fp), &q, &dn, nz - 1).iter().map(square_free_sign).collect())
            };
            MEMO_SF.with(|m| {
                let mut m = m.borrow_mut();
                if m.len() > 5000 {
                    m.clear();
                }
                m.insert(key, p.clone());
            });
            p
        }
    };
    let p = match p {
        None => return Ok(None),
        Some(p) => p,
    };
    let p: Vec<SP> = p.iter().map(|x| num_z(divide_out(x, factors))).collect();
    Ok(Some(Row { p, test, corner: true }))
}

fn sign_polys(s: &[M2], q: &M2, dn: &Z, deg: usize) -> Vec<M2> {
    let c2 = biv::conv2;
    let pa = biv::pad_add;
    if deg == 2 {
        let (a, b, c0) = (&s[2], &s[1], &s[0]);
        return vec![
            a.clone(),
            pa(&biv::scali(&c2(a, q), 2), &biv::scal(b, dn)),
            pa(&pa(&c2(a, &c2(q, q)), &biv::scal(&c2(b, q), dn)), &biv::scal(c0, &(dn * dn))),
            pa(&c2(b, b), &biv::scali(&c2(a, c0), -4)),
        ];
    }
    let (a, c0) = (&s[1], &s[0]);
    vec![a.clone(), pa(&c2(a, q), &biv::scal(c0, dn))]
}

fn top_row_polys(fp: &PolyN, q: &M2, dn: &Z, n: usize, want: i32) -> Vec<M2> {
    let r: R<Vec<M2>> = (|| {
        let g = slices(fp);
        let nj = if want == 1 { n - 1 } else { 0 };
        let mut p = vec![biv::content_free(&g[n])];
        let fact = |i: usize| -> Z { (1..=i).fold(Z::one(), |a, k| a * zi(k as i64)) };
        for j in 0..=nj {
            let c = |i: usize| -> Z { fact(i) / fact(i - j) };
            let mut rr = biv::scal(&g[n], &c(n));
            let mut i = n as i64 - 1;
            while i >= j as i64 {
                let iu = i as usize;
                rr = biv::pad_add(&biv::conv2(&rr, q), &biv::scal(&g[iu], &(c(iu) * num_traits::pow(dn.clone(), n - iu))));
                i -= 1;
            }
            p.push(square_free_sign(&biv::content_free(&biv::trim_n(&rr))));
        }
        Ok(p)
    })();
    r.unwrap_or_default()
}

// ---- the box recursion ----------------------------------------------------------------------------
fn box_rec(rows: &[Row], b: [f64; 4], in_box: &dyn Fn(&[f64; 4]) -> bool, depth: usize, sq: &[[f64; 3]]) -> bool {
    for c in sq {
        if b[0] >= c[0] - c[2] && b[1] <= c[0] + c[2] && b[2] >= c[1] - c[2] && b[3] <= c[1] + c[2] {
            return true;
        }
    }
    for row in rows {
        let s: Vec<i32> = row.p.iter().map(|p| poly_sign(p, &b)).collect();
        if run_test(&row.test, &s) {
            return true;
        }
    }
    if !in_box(&b) {
        return true;
    }
    if depth == 0 {
        return false;
    }
    let xm = (b[0] + b[1]) / 2.0;
    let ym = (b[2] + b[3]) / 2.0;
    box_rec(rows, [b[0], xm, b[2], ym], in_box, depth - 1, sq)
        && box_rec(rows, [xm, b[1], b[2], ym], in_box, depth - 1, sq)
        && box_rec(rows, [b[0], xm, ym, b[3]], in_box, depth - 1, sq)
        && box_rec(rows, [xm, b[1], ym, b[3]], in_box, depth - 1, sq)
}

fn poly_sign(p: &SP, b: &[f64; 4]) -> i32 {
    match p {
        SP::Num { f, .. } => box_sign(f, b),
        SP::Surd(s) => {
            let lo = box_sign(&s.lo, b);
            if lo != 0 && box_sign(&s.hi, b) == lo {
                return lo;
            }
            let ef = match &s.ef {
                None => return 0,
                Some(e) => e,
            };
            let (su, sv) = if let (Some(u), Some(v)) = (&s.usf, &s.vsf) { (box_sign(u, b), box_sign(v, b)) } else { (box_sign(&s.uf, b), box_sign(&s.vf, b)) };
            let se = box_sign(ef, b);
            sq_sign(su, sv, se)
        }
    }
}

fn eval_terms_f(m: &FMat, p: [f64; 2]) -> (f64, f64) {
    let mut v = 0.0;
    let mut mm = 0.0;
    for (i, r) in m.iter().enumerate() {
        let xi = p[0].powf(i as f64);
        for (j, c) in r.iter().enumerate() {
            if *c == 0.0 {
                continue;
            }
            let t = c * xi * p[1].powf(j as f64);
            v += t;
            mm += t.abs();
        }
    }
    (v, mm)
}

fn cmul(a: Cx, b: Cx) -> Cx {
    Cx { re: a.re * b.re - a.im * b.im, im: a.re * b.im + a.im * b.re }
}

fn cpow(z: Cx, n: usize) -> Cx {
    let mut r = Cx { re: 1.0, im: 0.0 };
    for _ in 0..n {
        r = cmul(r, z);
    }
    r
}

fn eval_terms_c(m: &M2, x: f64, y: Cx) -> (f64, f64) {
    let mut v = Cx { re: 0.0, im: 0.0 };
    let mut mm = 0.0;
    for (i, r) in m.iter().enumerate() {
        let xi = x.powf(i as f64);
        for (j, c) in r.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let yj = cpow(y, j);
            let s = to_f64(c) * xi;
            let t = Cx { re: s * yj.re, im: s * yj.im };
            v.re += t.re;
            v.im += t.im;
            mm += t.re.hypot(t.im);
        }
    }
    (v.re.hypot(v.im), mm)
}

fn vanishes_at(p: &SP, q: [f64; 2]) -> bool {
    match p {
        SP::Num { f, .. } => {
            let (v, m) = eval_terms_f(f, q);
            v.abs() <= 1e-9 * m.max(REALMIN)
        }
        SP::Surd(s) => {
            let a = eval_terms_f(&s.uf, q);
            let b = eval_terms_f(&s.vf, q);
            a.0.abs() <= 1e-9 * a.1.max(REALMIN) && b.0.abs() <= 1e-9 * b.1.max(REALMIN)
        }
    }
}

fn any_vanishes(rows: &[Row], p: [f64; 2]) -> bool {
    rows.iter().any(|r| r.p.iter().any(|q| vanishes_at(q, p)))
}

// ---- corners --------------------------------------------------------------------------------------
fn corner_ok(rows: &[Row], p: [f64; 2], e: f64, factors: &[M2], memo: &mut HashMap<(usize, usize, u64, u64, i64, u8), Option<SP>>) -> R<bool> {
    if !(e > 0.0) {
        return Ok(false);
    }
    for major in [1u8, 2] {
        for sg in [1i64, -1] {
            let mut rows_c = vec![];
            for (ri, row) in rows.iter().enumerate() {
                if !row.corner {
                    continue;
                }
                let mut ps = vec![];
                let mut ok = true;
                for (pi, pp) in row.p.iter().enumerate() {
                    let key = (ri, pi, p[0].to_bits(), p[1].to_bits(), sg, major);
                    let g = match memo.get(&key) {
                        Some(g) => g.clone(),
                        None => {
                            let g = blow_up_raw(pp, p, sg, major)?;
                            memo.insert(key, g.clone());
                            g
                        }
                    };
                    match g {
                        None => {
                            ok = false;
                            break;
                        }
                        Some(g) => ps.push(g),
                    }
                }
                if ok {
                    rows_c.push(Row { p: ps, test: row.test.clone(), corner: true });
                }
            }
            if rows_c.is_empty() {
                return Ok(false);
            }
            let mut fac_c: Vec<SP> = vec![];
            for g in factors {
                let key = (g.clone(), p[0].to_bits(), p[1].to_bits(), sg, major);
                let hit = MEMO_FACTOR_BU.with(|m| m.borrow().get(&key).cloned());
                let gg = match hit {
                    Some(x) => x,
                    None => {
                        let gn = num_z(g.clone());
                        let x = if vanishes_at(&gn, p) { blow_up_raw(&gn, p, sg, major)? } else { None };
                        MEMO_FACTOR_BU.with(|m| {
                            let mut m = m.borrow_mut();
                            if m.len() > 20000 {
                                m.clear();
                            }
                            m.insert(key, x.clone());
                        });
                        x
                    }
                };
                if let Some(gg) = gg {
                    fac_c.push(gg);
                }
            }
            let in_c = |cb: &[f64; 4]| !fac_c.iter().any(|g| poly_sign(g, cb) == -1);
            if !box_rec(&rows_c, [0.0, e, -1.0, 1.0], &in_c, 6, &[]) {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

fn is_dyadic(p: [f64; 2]) -> bool {
    p.iter().all(|v| v * 256.0 == round(v * 256.0))
}

fn to_integers_f(m: &FMat) -> Option<M2> {
    let mut m = m.clone();
    for _ in 0..30 {
        if m.iter().flatten().all(|x| *x == x.trunc()) {
            break;
        }
        m = m.iter().map(|r| r.iter().map(|x| 2.0 * x).collect()).collect();
    }
    if m.iter().flatten().any(|x| *x != x.trunc() || !x.is_finite()) {
        return None;
    }
    Some(m.iter().map(|r| r.iter().map(|x| q_from_f64(*x).to_integer()).collect()).collect())
}

fn blow_up_raw(p: &SP, q: [f64; 2], sg: i64, major: u8) -> R<Option<SP>> {
    let (f, z) = match p {
        SP::Surd(s) => return Ok(blow_up_pair(s, q, sg, major)),
        SP::Num { f, z } => (f, z),
    };
    let m = match z {
        Some(z) => z.clone(),
        None => match to_integers_f(f) {
            None => return Ok(None),
            Some(m) => m,
        },
    };
    if !is_dyadic(q) {
        return Ok(blow_up_surd(&m, q, sg, major));
    }
    let mut sc = 1.0f64;
    while q.iter().any(|v| v * sc != round(v * sc)) && sc < 256.0 {
        sc *= 2.0;
    }
    let (n1, n2) = biv::shape(&m);
    let deg = n1 + n2 - 2;
    let sum: f64 = m.iter().flatten().map(|x| to_f64(&x.abs())).sum();
    let bound = sum * (q[0].abs() + 1.0).powf((n1 - 1) as f64) * (q[1].abs() + 1.0).powf((n2 - 1) as f64) * 2f64.powf(deg as f64) * sc.powf(deg as f64);
    if bound >= 2f64.powi(52) {
        return Ok(None);
    }
    let px = q_from_f64(q[0]);
    let py = q_from_f64(q[1]);
    let mq: Vec<Vec<Q>> = m.iter().map(|r| r.iter().map(qz).collect()).collect();
    let ps = transpose_q(&shift_to(&transpose_q(&shift_to(&mq, &px)), &py));
    let mut g = vec![vec![Q::zero(); n1.max(n2)]; n1 + n2 - 1];
    for i in 0..n1 {
        for j in 0..n2 {
            let mut c = ps[i][j].clone();
            if sg < 0 && (i + j) % 2 == 1 {
                c = -c;
            }
            if major == 1 {
                g[i + j][j] += c;
            } else {
                g[i + j][i] += c;
            }
        }
    }
    let k = match g.iter().position(|r| r.iter().any(|x| !x.is_zero())) {
        None => return Ok(Some(SP::Num { f: vec![vec![0.0]], z: None })),
        Some(k) => k,
    };
    let f: FMat = g[k..].iter().map(|r| r.iter().map(q_to_f64).collect()).collect();
    Ok(Some(SP::Num { f, z: None }))
}

fn transpose_q(a: &[Vec<Q>]) -> Vec<Vec<Q>> {
    let n1 = a.len();
    let n2 = a[0].len();
    (0..n2).map(|j| (0..n1).map(|i| a[i][j].clone()).collect()).collect()
}

fn shift_to(p: &[Vec<Q>], x0: &Q) -> Vec<Vec<Q>> {
    let n1 = p.len();
    let mut s = vec![vec![Q::zero(); p[0].len()]; n1];
    for i in 0..n1 {
        for k in 0..=i {
            let mut f = Q::from_integer(Z::from(comb(i, k)));
            for _ in 0..(i - k) {
                f = f * x0;
            }
            for (j, v) in p[i].iter().enumerate() {
                if !v.is_zero() {
                    s[k][j] += &f * v;
                }
            }
        }
    }
    s
}

fn chart_z(zm: &M2, p: [f64; 2], sg: i64, major: u8) -> M2 {
    let mut c = zm.clone();
    let mut sc = 1.0f64;
    while p.iter().any(|v| v * sc != round(v * sc)) && sc < 256.0 {
        sc *= 2.0;
    }
    let scz = Z::from(sc as i64);
    let a = [Z::from(round(p[0] * sc) as i64), Z::from(round(p[1] * sc) as i64)];
    let (n1, n2) = biv::shape(&c);
    let deg = n1 + n2 - 2;
    for i in 0..n1 {
        for j in 0..n2 {
            let e = deg as i64 - i as i64 - j as i64;
            if e > 0 {
                c[i][j] = &c[i][j] * num_traits::pow(scz.clone(), e as usize);
            }
        }
    }
    for j in 0..n2 {
        for i in 0..n1.saturating_sub(1) {
            let mut k = n1 as i64 - 2;
            while k >= i as i64 {
                let ku = k as usize;
                let add = &a[0] * &c[ku + 1][j];
                c[ku][j] += add;
                k -= 1;
            }
        }
    }
    for i in 0..n1 {
        for j in 0..n2.saturating_sub(1) {
            let mut k = n2 as i64 - 2;
            while k >= j as i64 {
                let ku = k as usize;
                let add = &a[1] * &c[i][ku + 1];
                c[i][ku] += add;
                k -= 1;
            }
        }
    }
    let mut g = biv::zeros(n1 + n2 - 1, n1.max(n2));
    for k in 0..n1 {
        for l in 0..n2 {
            let mut v = c[k][l].clone();
            if v.is_zero() {
                continue;
            }
            if k + l > 0 {
                v = v * num_traits::pow(scz.clone(), k + l);
            }
            if sg < 0 && (k + l) % 2 == 1 {
                v = -v;
            }
            if major == 1 {
                g[k + l][l] += v;
            } else {
                g[k + l][k] += v;
            }
        }
    }
    g
}

fn trim_s(z: &M2) -> M2 {
    for i in 0..z.len() {
        if z[i].iter().any(|x| !x.is_zero()) {
            return biv::content_free(&z[i..].to_vec());
        }
    }
    vec![vec![Z::zero()]]
}

fn blow_up_pair(p: &Surd, q: [f64; 2], sg: i64, major: u8) -> Option<SP> {
    if p.ez.is_none() || !is_dyadic(q) {
        return None;
    }
    let gu = chart_z(p.uz.as_ref()?, q, sg, major);
    let gv = chart_z(p.vz.as_ref()?, q, sg, major);
    let n1 = gu.len().max(gv.len());
    let n2 = gu[0].len().max(gv[0].len());
    let gu = biv::pad(&gu, n1, n2);
    let gv = biv::pad(&gv, n1, n2);
    let k = (0..n1).find(|&i| gu[i].iter().any(|x| !x.is_zero()) || gv[i].iter().any(|x| !x.is_zero()));
    let k = match k {
        None => return Some(SP::Num { f: vec![vec![0.0]], z: Some(vec![vec![Z::zero()]]) }),
        Some(k) => k,
    };
    let ge = trim_s(&chart_z(p.ez.as_ref()?, q, sg, major));
    let gus = trim_s(&chart_z(p.usz.as_ref()?, q, sg, major));
    let gvs = trim_s(&chart_z(p.vsz.as_ref()?, q, sg, major));
    Some(SP::Surd(Box::new(make_surd(gu[k..].to_vec(), gv[k..].to_vec(), p.k.clone(), ge, gus, gvs))))
}

// ---- surds: U + sqrt(k) V -------------------------------------------------------------------------
fn lohi(u: &FMat, v: &FMat, k: &Z) -> (FMat, FMat) {
    let rd = to_f64(k).sqrt();
    let lo = rd * (1.0 - 4.0 * EPS);
    let hi = rd * (1.0 + 4.0 * EPS);
    let l = u.iter().zip(v.iter()).map(|(ru, rv)| ru.iter().zip(rv.iter()).map(|(a, b)| a + lo * b).collect()).collect();
    let h = u.iter().zip(v.iter()).map(|(ru, rv)| ru.iter().zip(rv.iter()).map(|(a, b)| a + hi * b).collect()).collect();
    (l, h)
}

fn make_surd(uz: M2, vz: M2, k: Z, ez: M2, usz: M2, vsz: M2) -> Surd {
    let e = biv::content_free(&ez);
    let e = if biv::is_zero(&e) { None } else { Some(e) };
    let n1 = uz.len().max(vz.len());
    let n2 = uz[0].len().max(vz[0].len());
    let u = biv::pad(&uz, n1, n2);
    let v = biv::pad(&vz, n1, n2);
    let uf = fmat(&u);
    let vf = fmat(&v);
    let (lo, hi) = lohi(&uf, &vf, &k);
    let us = biv::content_free(&usz);
    let vs = biv::content_free(&vsz);
    Surd {
        ef: e.as_ref().map(fmat),
        usf: Some(fmat(&us)),
        vsf: Some(fmat(&vs)),
        uf: fmat(&uz),
        vf: fmat(&vz),
        u: uz.clone(),
        v: vz.clone(),
        e,
        us: Some(us),
        vs: Some(vs),
        lo,
        hi,
        uz: Some(uz),
        vz: Some(vz),
        ez: Some(ez),
        usz: Some(usz),
        vsz: Some(vsz),
        k,
    }
}

fn surd_poly(u: &M2, v: &M2, k: &Z) -> Surd {
    let n1 = u.len().max(v.len());
    let n2 = u[0].len().max(v[0].len());
    let mut u = biv::pad(u, n1, n2);
    let mut v = biv::pad(v, n1, n2);
    let mut g = Z::zero();
    for r in u.iter().chain(v.iter()) {
        for x in r {
            g = g.gcd(x);
        }
    }
    if g > Z::one() {
        u = u.iter().map(|r| r.iter().map(|x| x.div_floor(&g)).collect()).collect();
        v = v.iter().map(|r| r.iter().map(|x| x.div_floor(&g)).collect()).collect();
    }
    let e = biv::pad_add(&biv::conv2(&u, &u), &biv::scal(&biv::conv2(&v, &v), &-k));
    let (oe, ou, ov) = (odd_part_z(&e), odd_part_z(&u), odd_part_z(&v));
    make_surd(u, v, k.clone(), oe, ou, ov)
}

fn odd_part_z(zm: &M2) -> M2 {
    let d = biv::trim_n(zm);
    if biv::is_const(&d) {
        return vec![vec![zi(sign(&d[0][0]) as i64)]];
    }
    if let Ok(o) = odd_part(&d, 1, 0) {
        let r = biv::content_free(&biv::trim_n(&o));
        let sg = rel_sign(&d, &r);
        if sg != 0 {
            return biv::scali(&r, sg as i64);
        }
    }
    zm.clone()
}

const SURD_D: [i64; 30] = [2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31, 33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47];

fn surd_of(p: [f64; 2]) -> Option<[i64; 6]> {
    let key = (p[0].to_bits(), p[1].to_bits());
    if let Some(r) = MEMO_SURDOF.with(|m| m.borrow().get(&key).cloned()) {
        return r;
    }
    let r = surd_of_search(p);
    MEMO_SURDOF.with(|m| {
        let mut m = m.borrow_mut();
        if m.len() > 20000 {
            m.clear();
        }
        m.insert(key, r);
    });
    r
}

fn surd_of_search(p: [f64; 2]) -> Option<[i64; 6]> {
    let sd: Vec<f64> = SURD_D.iter().map(|d| (*d as f64).sqrt()).collect();
    for q in 1..65i64 {
        let qf = q as f64;
        for (k, &d) in SURD_D.iter().enumerate() {
            let mut ab: Vec<(i64, i64)> = vec![];
            for c in 0..2 {
                let mut best: Option<(i64, i64)> = None;
                for b in -64i64..=64 {
                    let v = qf * p[c] - (b as f64) * sd[k];
                    let a = round(v);
                    if (v - a).abs() < 1e-11 * qf * 1f64.max(p[c].abs()) {
                        if best.map_or(true, |bb| b.abs() < bb.1.abs()) {
                            best = Some((a as i64, b));
                        }
                    }
                }
                match best {
                    None => break,
                    Some(bb) => ab.push(bb),
                }
            }
            if ab.len() == 2 && ab.iter().any(|x| x.1 != 0) {
                return Some([ab[0].0, ab[0].1, ab[1].0, ab[1].1, q, d]);
            }
        }
    }
    None
}

fn surd_mul(x: &[Z; 2], y: &[Z; 2], d: &Z) -> [Z; 2] {
    [&x[0] * &y[0] + d * &x[1] * &y[1], &x[0] * &y[1] + &x[1] * &y[0]]
}

fn surd_pow(a: &[Z; 2], n: usize, d: &Z) -> Vec<[Z; 2]> {
    let mut m = vec![[Z::one(), Z::zero()]];
    for _ in 0..n {
        let nx = surd_mul(m.last().unwrap(), a, d);
        m.push(nx);
    }
    m
}

fn blow_up_surd(pm: &M2, p: [f64; 2], sg: i64, major: u8) -> Option<SP> {
    let r = surd_of(p)?;
    let a = [zi(r[0]), zi(r[1])];
    let b = [zi(r[2]), zi(r[3])];
    let q = r[4];
    let d = r[5];
    let (n1, n2) = biv::shape(pm);
    let deg = n1 + n2 - 2;
    let sdd = (d as f64).sqrt();
    let na = (r[0].abs() as f64) + sdd * (r[1].abs() as f64);
    let nb = (r[2].abs() as f64) + sdd * (r[3].abs() as f64);
    let qf = q as f64;
    let mut bound = 0.0;
    for i in 0..n1 {
        for j in 0..n2 {
            bound += to_f64(&pm[i][j].abs()) * (na + qf).powf(i as f64) * (nb + qf).powf(j as f64) * qf.powf((deg - i - j) as f64);
        }
    }
    if bound >= 2f64.powi(52) {
        return None;
    }
    let dz = zi(d);
    let qz_ = zi(q);
    let pa = surd_pow(&a, n1 - 1, &dz);
    let pb = surd_pow(&b, n2 - 1, &dz);
    let mut su = biv::zeros(n1, n2);
    let mut sv = biv::zeros(n1, n2);
    for k in 0..n1 {
        for l in 0..n2 {
            let mut acc = [Z::zero(), Z::zero()];
            for i in k..n1 {
                for j in l..n2 {
                    if pm[i][j].is_zero() {
                        continue;
                    }
                    let m = surd_mul(&pa[i - k], &pb[j - l], &dz);
                    let f = &pm[i][j] * Z::from(comb(i, k)) * Z::from(comb(j, l)) * num_traits::pow(qz_.clone(), k + l + deg - i - j);
                    acc[0] += &f * &m[0];
                    acc[1] += &f * &m[1];
                }
            }
            su[k][l] = acc[0].clone();
            sv[k][l] = acc[1].clone();
        }
    }
    let mut u = biv::zeros(n1 + n2 - 1, n1.max(n2));
    let mut v = biv::zeros(n1 + n2 - 1, n1.max(n2));
    for k in 0..n1 {
        for l in 0..n2 {
            let neg = sg < 0 && (k + l) % 2 == 1;
            let col = if major == 1 { l } else { k };
            if neg {
                u[k + l][col] -= &su[k][l];
                v[k + l][col] -= &sv[k][l];
            } else {
                u[k + l][col] += &su[k][l];
                v[k + l][col] += &sv[k][l];
            }
        }
    }
    let k0 = (0..u.len()).find(|&i| u[i].iter().any(|x| !x.is_zero()) || v[i].iter().any(|x| !x.is_zero()));
    let k0 = match k0 {
        None => return Some(SP::Num { f: vec![vec![0.0]], z: Some(vec![vec![Z::zero()]]) }),
        Some(k) => k,
    };
    let u = u[k0..].to_vec();
    let v = v[k0..].to_vec();
    let au: M2 = u.iter().map(|r| r.iter().map(|x| x.abs()).collect()).collect();
    let av: M2 = v.iter().map(|r| r.iter().map(|x| x.abs()).collect()).collect();
    let bnd = biv::pad_add(&biv::conv2(&au, &au), &biv::scal(&biv::conv2(&av, &av), &dz));
    let e = if bnd.iter().flatten().max().unwrap() < &(Z::one() << 53usize) { Some(biv::pad_add(&biv::conv2(&u, &u), &biv::scal(&biv::conv2(&v, &v), &-&dz))) } else { None };
    let uf = fmat(&u);
    let vf = fmat(&v);
    let (lo, hi) = lohi(&uf, &vf, &dz);
    Some(SP::Surd(Box::new(Surd { ef: e.as_ref().map(fmat), usf: None, vsf: None, uf, vf, u, v, e, us: None, vs: None, lo, hi, uz: None, vz: None, ez: None, usz: None, vsz: None, k: dz })))
}

// ---- rows comparing two faces ---------------------------------------------------------------------
fn face_row_polys(fk: &Rc<PolyN>, fbk: i64, fw: &Rc<PolyN>, fbw: i64, factors: &[M2]) -> R<Option<Row>> {
    let nk = nz_of(fk);
    let nw = nz_of(fw);
    if nw >= 4 {
        if !(nk == 2 || nk == 3) || fbw != -1 {
            return Ok(None);
        }
        let mut hik = false;
        if nk == 3 {
            hik = fbk == -1 || fbk == 2;
            if !hik && fbk != 1 && fbk != -2 {
                return Ok(None);
            }
        } else if fbk != 1 && fbk != -1 {
            return Ok(None);
        }
        return top_face_row(fk, fw, nk, hik, factors);
    }
    let mode = match (nk, nw) {
        (3, 2) => 1u8,
        (2, 3) => 2,
        (3, 3) => 0,
        _ => return Ok(None),
    };
    let (mut hik, mut hiw) = (false, false);
    if nk == 3 {
        hik = fbk == -1 || fbk == 2;
        if !hik && fbk != 1 && fbk != -2 {
            return Ok(None);
        }
    } else if fbk != 1 && fbk != -1 {
        return Ok(None);
    }
    if nw == 3 {
        hiw = fbw == -1 || fbw == 2;
        if !hiw && fbw != 1 && fbw != -2 {
            return Ok(None);
        }
    } else if fbw != 1 && fbw != -1 {
        return Ok(None);
    }
    if let Some(r) = surd_face_row(fk, fw, nk, nw, hik, hiw, factors)? {
        return Ok(Some(r));
    }
    let key = (mode, (**fk).clone(), (**fw).clone());
    let p = match MEMO_FACE.with(|m| m.borrow().get(&key).cloned()) {
        Some(p) => p,
        None => {
            let p = face_raw(&slices(fk), &slices(fw), -1, mode).map(|v| v.iter().map(square_free_sign).collect::<Vec<_>>()).ok();
            MEMO_FACE.with(|m| {
                let mut m = m.borrow_mut();
                if m.len() > 2000 {
                    m.clear();
                }
                m.insert(key, p.clone());
            });
            p
        }
    };
    let p = match p {
        None => return Ok(None),
        Some(p) => p,
    };
    let p: Vec<SP> = p.iter().map(|x| num_z(divide_out(x, factors))).collect();
    let test = match mode {
        0 => Test::Fq { hik, hiw },
        1 => Test::Ql { hik },
        _ => Test::Lq { hiw },
    };
    Ok(Some(Row { p, test, corner: true }))
}

fn face_raw(fk: &[M2], fw: &[M2], sg: i64, mode: u8) -> R<Vec<M2>> {
    let c2 = biv::conv2;
    let pa = biv::pad_add;
    let sc = biv::scali;
    if mode != 0 {
        let (q, l) = if mode == 1 { (fk, fw) } else { (fw, fk) };
        let (a, b, c) = (&q[2], &q[1], &q[0]);
        let (bl, cl) = (&l[1], &l[0]);
        let g = pa(&c2(b, bl), &sc(&c2(a, cl), 2 * sg));
        let h = pa(&pa(&c2(a, &c2(cl, cl)), &sc(&c2(b, &c2(cl, bl)), sg)), &c2(c, &c2(bl, bl)));
        let d = pa(&c2(b, b), &sc(&c2(a, c), 4 * sg));
        return Ok(vec![a.clone(), bl.clone(), g, h, d]);
    }
    let (ak, bk, ck) = (&fk[2], &fk[1], &fk[0]);
    let (aw, bw, cw) = (&fw[2], &fw[1], &fw[0]);
    let dw = pa(&c2(bw, bw), &sc(&c2(aw, cw), 4 * sg));
    let dk = pa(&c2(bk, bk), &sc(&c2(ak, ck), 4 * sg));
    let p1 = pa(&c2(aw, bk), &sc(&c2(ak, bw), sg));
    let e1 = pa(&c2(&p1, &p1), &sc(&c2(&c2(ak, ak), &dw), sg));
    let p2 = pa(&sc(&c2(&p1, bw), sg), &sc(&c2(aw, &pa(&c2(aw, ck), &sc(&c2(ak, cw), sg))), 2));
    let e2 = pa(&c2(&p2, &p2), &sc(&c2(&c2(&p1, &p1), &dw), sg));
    Ok(vec![aw.clone(), ak.clone(), dw, dk, p1, e1, p2, e2])
}

#[derive(Clone, Debug)]
struct SurdFace {
    lin: bool,
    den: M2,
    a0: M2,
    s: M2,
    c: Z,
    eta: i64,
}

fn surd_face(f: &PolyN, n: usize, hi: bool) -> Option<SurdFace> {
    let key = (f.clone(), hi);
    if let Some(r) = MEMO_SURDFACE.with(|m| m.borrow().get(&key).cloned()) {
        return r;
    }
    let s = slices(f);
    let mut out = None;
    if n == 2 {
        out = Some(SurdFace { lin: true, den: s[1].clone(), a0: biv::scali(&s[0], -1), s: vec![vec![zi(1)]], c: zi(0), eta: 1 });
    } else {
        let (a, b, c) = (&s[2], &s[1], &s[0]);
        let d = biv::pad_add(&biv::conv2(b, b), &biv::scali(&biv::conv2(a, c), -4));
        let g = biv::content(&d);
        if g.is_positive() {
            let dp: M2 = d.iter().map(|r| r.iter().map(|x| x.div_floor(&g)).collect()).collect();
            if let Some(sq) = z_sqrt(&dp) {
                out = Some(SurdFace { lin: false, den: biv::scali(a, 2), a0: biv::scali(b, -1), s: sq, c: g, eta: 2 * (hi as i64) - 1 });
            }
        }
    }
    MEMO_SURDFACE.with(|m| {
        let mut m = m.borrow_mut();
        if m.len() > 2000 {
            m.clear();
        }
        m.insert(key, out.clone());
    });
    out
}

fn field_mu(f: &SurdFace, kap: &Z) -> Option<Z> {
    if f.lin {
        return Some(Z::zero());
    }
    let n = &f.c * kap;
    let s = n.sqrt();
    if &s * &s == n {
        Some(s)
    } else {
        None
    }
}

fn primes() -> &'static [i64] {
    static P: std::sync::OnceLock<Vec<i64>> = std::sync::OnceLock::new();
    P.get_or_init(|| {
        let mut sieve = vec![true; 5001];
        let mut out = vec![];
        for i in 2..5001usize {
            if sieve[i] {
                out.push(i as i64);
                let mut j = i * i;
                while j < 5001 {
                    sieve[j] = false;
                    j += i;
                }
            }
        }
        out
    })
}

fn square_free_int(c: &Z) -> Z {
    if c.is_zero() {
        return Z::zero(); // callers pass a positive content; 0 would never leave the loop below
    }
    let mut k = Z::one();
    let mut r = c.clone();
    for &i in primes() {
        let p = zi(i);
        let mut e = 0;
        while (&r % &p).is_zero() {
            r = &r / &p;
            e += 1;
        }
        if e % 2 == 1 {
            k = k * &p;
        }
        if r.is_one() {
            return k;
        }
    }
    let s = r.sqrt();
    if &s * &s != r {
        k = k * r;
    }
    k
}

fn surd_face_row(fk: &Rc<PolyN>, fw: &Rc<PolyN>, nk: usize, nw: usize, hik: bool, hiw: bool, factors: &[M2]) -> R<Option<Row>> {
    let key = ((**fk).clone(), (**fw).clone(), hik, hiw);
    let p = match MEMO_SURD.with(|m| m.borrow().get(&key).cloned()) {
        Some(p) => p,
        None => {
            let mut p = None;
            let fkk = surd_face(fk, nk, hik);
            let fww = surd_face(fw, nw, hiw);
            if let (Some(fk2), Some(fw2)) = (&fkk, &fww) {
                if !(fk2.lin && fw2.lin) {
                    let kap = square_free_int(if !fw2.lin { &fw2.c } else { &fk2.c });
                    if kap.bits() <= 52 {
                        if let (Some(muk), Some(muw)) = (field_mu(fk2, &kap), field_mu(fw2, &kap)) {
                            let k = kap.clone();
                            let c2 = biv::conv2;
                            let u = biv::scal(&biv::pad_add(&c2(&fk2.den, &fw2.a0), &biv::scali(&c2(&fw2.den, &fk2.a0), -1)), &k);
                            let wk = c2(&fw2.den, &fk2.s);
                            let ww = c2(&fk2.den, &fw2.s);
                            let mut t: Vec<Option<Surd>> = vec![None, None, None, None];
                            for tk in [1i64, -1] {
                                for tw in [1i64, -1] {
                                    let v = biv::pad_add(&biv::scal(&ww, &(&muw * zi(tw * fw2.eta))), &biv::scal(&wk, &(-&muk * zi(tk * fk2.eta))));
                                    t[(if tk < 0 { 2 } else { 0 }) + (if tw < 0 { 1 } else { 0 })] = Some(surd_poly(&u, &v, &k));
                                }
                            }
                            let cf = biv::content_free;
                            p = Some((vec![cf(&fk2.den), cf(&fw2.den), cf(&fk2.s), cf(&fw2.s)], t.into_iter().map(|x| x.unwrap()).collect()));
                        }
                    }
                }
            }
            MEMO_SURD.with(|m| {
                let mut m = m.borrow_mut();
                if m.len() > 2000 {
                    m.clear();
                }
                m.insert(key, p.clone());
            });
            p
        }
    };
    let (pp, tt) = match p {
        None => return Ok(None),
        Some(x) => x,
    };
    let mut pn: Vec<SP> = pp.iter().map(|x| num_z(divide_out(x, factors))).collect();
    for t in &tt {
        pn.push(SP::Surd(Box::new(divide_out_surd(t, factors))));
    }
    Ok(Some(Row { p: pn, test: Test::Surd, corner: true }))
}

fn top_face_row(fk: &Rc<PolyN>, fw: &Rc<PolyN>, nk: usize, hik: bool, factors: &[M2]) -> R<Option<Row>> {
    let key = ((**fk).clone(), (**fw).clone(), hik);
    let p = match MEMO_TOP.with(|m| m.borrow().get(&key).cloned()) {
        Some(p) => p,
        None => {
            let mut p = None;
            if let Some(fk2) = surd_face(fk, nk, hik) {
                let kap = if fk2.lin { Z::one() } else { square_free_int(&fk2.c) };
                let mu = field_mu(&fk2, &kap);
                if kap.bits() <= 52 && mu.is_some() {
                    let mu = mu.unwrap();
                    let k = kap.clone();
                    let g = slices(fw);
                    let n = g.len() - 1;
                    let c2 = biv::conv2;
                    let x = biv::scal(&fk2.a0, &k);
                    let m = biv::scal(&fk2.den, &k);
                    let mut t: Vec<Surd> = vec![];
                    for tk in [1i64, -1] {
                        let y = biv::scal(&fk2.s, &(&mu * zi(tk * fk2.eta)));
                        let mut u = g[n].clone();
                        let mut v = vec![vec![Z::zero()]];
                        let mut mp = vec![vec![Z::one()]];
                        for i in (0..n).rev() {
                            mp = c2(&mp, &m);
                            let un = biv::pad_add(&c2(&u, &x), &biv::scal(&c2(&v, &y), &k));
                            v = biv::pad_add(&c2(&u, &y), &c2(&v, &x));
                            u = biv::pad_add(&un, &c2(&g[i], &mp));
                        }
                        t.push(surd_poly(&biv::trim_n(&u), &biv::trim_n(&v), &k));
                    }
                    let cf = biv::content_free;
                    p = Some((vec![cf(&fk2.den), cf(&fk2.s), cf(&g[n])], t, n % 2 == 1));
                }
            }
            MEMO_TOP.with(|m| {
                let mut m = m.borrow_mut();
                if m.len() > 2000 {
                    m.clear();
                }
                m.insert(key, p.clone());
            });
            p
        }
    };
    let (pp, tt, odd) = match p {
        None => return Ok(None),
        Some(x) => x,
    };
    let mut pn: Vec<SP> = pp.iter().map(|x| num_z(divide_out(x, factors))).collect();
    for t in &tt {
        pn.push(SP::Surd(Box::new(divide_out_surd(t, factors))));
    }
    Ok(Some(Row { p: pn, test: Test::TopFace { odd }, corner: false }))
}

fn top_ok(s: &[i32], odd: bool) -> bool {
    if s[0..3].iter().any(|v| *v == 0) {
        return false;
    }
    let t = s[3 + if s[0] * s[1] < 0 { 1 } else { 0 }];
    let dn = if odd { s[0] } else { 1 };
    t * s[2] * dn == -1
}

fn surd_ok(s: &[i32]) -> bool {
    if s[0..4].iter().any(|v| *v == 0) {
        return false;
    }
    let tk = s[0] * s[2];
    let tw = s[1] * s[3];
    let t = s[4 + if tk < 0 { 2 } else { 0 } + if tw < 0 { 1 } else { 0 }];
    t != 0 && t * s[0] * s[1] == 1
}

fn z_sqrt(dp: &M2) -> Option<M2> {
    let p = biv::trim_n(dp);
    if biv::is_const(&p) {
        let v = &p[0][0];
        if v.is_positive() {
            let s = v.sqrt();
            if &s * &s == *v {
                return Some(vec![vec![s]]);
            }
        }
        return None;
    }
    let l = yun_parts(&p, 1, 0).ok()?;
    let mut q = vec![vec![Z::one()]];
    for (fac, mult) in l {
        if mult % 2 == 1 {
            return None;
        }
        for _ in 0..mult / 2 {
            q = biv::trim_n(&biv::conv2(&q, &fac));
        }
    }
    let q = biv::content_free(&biv::trim_n(&q));
    if biv::trim_n(&biv::conv2(&q, &q)) == p {
        Some(q)
    } else {
        None
    }
}

fn yun_parts(s: &M2, d: usize, lvl: usize) -> R<Vec<(M2, usize)>> {
    let mut l = vec![];
    let s = biv::trim_n(s);
    if biv::is_const(&s) || lvl > 2 {
        return Ok(l);
    }
    let sd = biv::diff(&s, d);
    if biv::is_zero(&sd) {
        if d == 1 {
            l = yun_parts(&s, 2, lvl + 1)?;
        }
        return Ok(l);
    }
    let mut c = biv::trim_n(&biv::poly_gcd_biv(&s, &sd)?);
    let mut w = biv::trim_n(&biv::poly_div_exact(&s, &c)?);
    let mut i = 1;
    while !biv::is_const(&w) {
        if i > 60 {
            return Err(crate::errors::cerr("algRowFails:yun", "square-free decomposition did not end"));
        }
        let y = biv::trim_n(&biv::poly_gcd_biv(&w, &c)?);
        let zz = biv::trim_n(&biv::poly_div_exact(&w, &y)?);
        if !biv::is_const(&zz) {
            l.push((zz, i));
        }
        c = biv::trim_n(&biv::poly_div_exact(&c, &y)?);
        w = y;
        i += 1;
    }
    if d == 1 && !biv::is_const(&c) {
        l.extend(yun_parts(&c, 2, lvl + 1)?);
    }
    Ok(l)
}

// ---- the sign logic --------------------------------------------------------------------------------
fn face_ok(s: &[i32], hik: bool, hiw: bool) -> bool {
    if s[2] == -1 || s[3] == -1 {
        return true;
    }
    if s[0] == 0 || s[1] == 0 {
        return false;
    }
    let e = if hiw { s[0] } else { -s[0] };
    let u = s[0] * s[1] * sq_sign(s[4], e * s[1], s[5]);
    let v = s[1] * sq_sign(s[6], e * s[4], s[7]);
    if hik {
        return u == 1 && v == 1;
    }
    u == 1 || v == -1
}

fn sq_sign(sp: i32, sq: i32, se: i32) -> i32 {
    if sp != 0 && sp == sq {
        return sp;
    }
    if se == 1 && sp != 0 {
        return sp;
    }
    if se == -1 && sq != 0 {
        return sq;
    }
    0
}

fn quad_ok(s: &[i32], hi: bool, want: i32) -> bool {
    if s[0] == 0 {
        return false;
    }
    let u = s[1] * s[0];
    let v = s[2] * s[0];
    if hi && want == -1 {
        return u == -1 || v == -1;
    }
    if hi {
        return u == 1 && v == 1;
    }
    if want == -1 {
        return u == -1 && v == 1;
    }
    u == 1 || v == -1
}

// ---- square-free parts -----------------------------------------------------------------------------
fn square_free_sign(r: &M2) -> M2 {
    let first: R<Option<M2>> = (|| {
        let s0 = biv::content_free(&biv::trim_n(r));
        if biv::is_const(&s0) {
            return Ok(Some(s0));
        }
        let o = biv::content_free(&biv::trim_n(&odd_part(&s0, 1, 0)?));
        let sg = rel_sign(&s0, &o);
        if sg != 0 {
            return Ok(Some(biv::scali(&o, sg as i64)));
        }
        Ok(None)
    })();
    if let Ok(Some(x)) = first {
        return x;
    }
    let mut s = r.clone();
    let second: R<M2> = (|| {
        s = biv::content_free(&biv::trim_n(r));
        for _ in 0..4 {
            let mut g = vec![vec![Z::one()]];
            for d in [1usize, 2] {
                let sd = biv::diff(&s, d);
                if !biv::is_zero(&sd) {
                    g = biv::trim_n(&biv::poly_gcd_biv(&s, &sd)?);
                }
                if !biv::is_const(&g) {
                    break;
                }
            }
            if biv::is_const(&g) {
                return Ok(s.clone());
            }
            s = biv::content_free(&biv::trim_n(&biv::poly_div_exact(&s, &biv::conv2(&g, &g))?));
        }
        Ok(s.clone())
    })();
    match second {
        Ok(x) => x,
        Err(_) => s,
    }
}

/// The product of the square-free factors of ODD multiplicity (Yun).
fn odd_part(s: &M2, d: usize, lvl: usize) -> R<M2> {
    if lvl == 0 && d == 1 && crate::flint::available() {
        if let Some(o) = crate::flint::odd_part(&biv::trim_n(s)) {
            return Ok(o);
        }
    }
    let mut o = vec![vec![Z::one()]];
    let s = biv::trim_n(s);
    if biv::is_const(&s) || lvl > 2 {
        return Ok(o);
    }
    let sd = biv::diff(&s, d);
    if biv::is_zero(&sd) {
        if d == 1 {
            o = odd_part(&s, 2, lvl + 1)?;
        }
        return Ok(o);
    }
    let mut c = biv::trim_n(&biv::poly_gcd_biv(&s, &sd)?);
    let mut w = biv::trim_n(&biv::poly_div_exact(&s, &c)?);
    let mut i = 1;
    while !biv::is_const(&w) {
        if i > 60 {
            return Err(crate::errors::cerr("algRowFails:yun", "square-free decomposition did not end"));
        }
        let y = biv::trim_n(&biv::poly_gcd_biv(&w, &c)?);
        let zz = biv::trim_n(&biv::poly_div_exact(&w, &y)?);
        if i % 2 == 1 {
            o = biv::trim_n(&biv::conv2(&o, &zz));
        }
        c = biv::trim_n(&biv::poly_div_exact(&c, &y)?);
        w = y;
        i += 1;
    }
    if d == 1 && !biv::is_const(&c) {
        o = biv::trim_n(&biv::conv2(&o, &odd_part(&c, 2, lvl + 1)?));
    }
    Ok(o)
}

const RELPTS: [[f64; 2]; 6] = [[0.6180339887, 0.3247179572], [-0.7182818285, 0.4142135624], [0.2360679775, -0.8660254038], [1.7320508076, 1.4142135624], [-1.2599210499, -0.5772156649], [3.1415926536, -2.7182818285]];

fn rel_sign(s: &M2, o: &M2) -> i32 {
    let sf = fmat(s);
    let of = fmat(o);
    for p in RELPTS {
        let (vs, ms) = eval_terms_f(&sf, p);
        let (vo, mo) = eval_terms_f(&of, p);
        if vs.abs() > 1e-6 * ms.max(REALMIN) && vo.abs() > 1e-6 * mo.max(REALMIN) {
            return signf(vs) * signf(vo);
        }
    }
    0
}

// ---- dividing out the region's own boundary ----------------------------------------------------------
fn factor_ints(g: &M2) -> Option<M2> {
    if let Some(h) = MEMO_FACTOR_INTS.with(|m| m.borrow().get(g).cloned()) {
        return h;
    }
    let h = if biv::is_const(g) {
        None
    } else {
        let g2 = biv::content_free(&biv::trim_n(g));
        if biv::is_const(&g2) {
            None
        } else {
            Some(g2)
        }
    };
    MEMO_FACTOR_INTS.with(|m| {
        let mut m = m.borrow_mut();
        if m.len() > 20000 {
            m.clear();
        }
        m.insert(g.clone(), h.clone());
    });
    h
}

fn divide_out(s: &M2, factors: &[M2]) -> M2 {
    let mut s = s.clone();
    for g0 in factors {
        let g = match factor_ints(g0) {
            None => continue,
            Some(g) => g,
        };
        if !may_divide(&s, &g) {
            continue;
        }
        for _ in 0..4 {
            match biv::poly_div_exact(&s, &g) {
                Ok(q) => s = biv::trim_n(&q),
                Err(_) => break,
            }
        }
    }
    s
}

fn divide_out_surd(t: &Surd, factors: &[M2]) -> Surd {
    let mut u = t.uz.clone().unwrap();
    let mut v = t.vz.clone().unwrap();
    let mut ch = false;
    for g0 in factors {
        let g = match factor_ints(g0) {
            None => continue,
            Some(g) => g,
        };
        for _ in 0..4 {
            if !may_divide(&u, &g) || (!biv::is_zero(&v) && !may_divide(&v, &g)) {
                break;
            }
            let u2 = match biv::poly_div_exact(&u, &g) {
                Ok(x) => biv::trim_n(&x),
                Err(_) => break,
            };
            let v2 = if !biv::is_zero(&v) {
                match biv::poly_div_exact(&v, &g) {
                    Ok(x) => biv::trim_n(&x),
                    Err(_) => break,
                }
            } else {
                vec![vec![Z::zero()]]
            };
            u = u2;
            v = v2;
            ch = true;
        }
    }
    let e = t.ez.clone().unwrap();
    let e2 = divide_out(&e, factors);
    ch = ch || e2 != e;
    let us = t.usz.clone().unwrap();
    let u3 = divide_out(&us, factors);
    ch = ch || u3 != us;
    let mut vs = t.vsz.clone().unwrap();
    if !biv::is_zero(&vs) {
        let v3 = divide_out(&vs, factors);
        ch = ch || v3 != vs;
        vs = v3;
    }
    if !ch {
        return t.clone();
    }
    let n1 = u.len().max(v.len());
    let n2 = u[0].len().max(v[0].len());
    make_surd(biv::pad(&u, n1, n2), biv::pad(&v, n1, n2), t.k.clone(), e2, u3, vs)
}

fn may_divide(s: &M2, g: &M2) -> bool {
    let nonconst = s.iter().enumerate().any(|(i, r)| r.iter().enumerate().any(|(j, x)| (i, j) != (0, 0) && !x.is_zero()));
    if !nonconst {
        return s[0][0].is_zero();
    }
    let cp = curve_points(g);
    let (tr, pts) = (&cp.0, &cp.1);
    let pts = match pts {
        None => return true,
        Some(p) => p,
    };
    let st;
    let s = if *tr {
        st = biv::transpose(s);
        &st
    } else {
        s
    };
    for (x, ys) in pts {
        for y in ys {
            let (v, m) = eval_terms_c(s, *x, *y);
            if v > 1e-6 * m.max(REALMIN) {
                return false;
            }
        }
    }
    true
}

fn curve_points(g: &M2) -> Rc<(bool, Option<Vec<(f64, Vec<Cx>)>>)> {
    if let Some(h) = MEMO_CURVE_PTS.with(|m| m.borrow().get(g).cloned()) {
        return h;
    }
    let mut tr = false;
    let mut gg = g.clone();
    if gg[0].len() < 2 {
        gg = biv::transpose(&gg);
        tr = true;
    }
    let hit = if gg[0].len() < 2 {
        Rc::new((tr, None))
    } else {
        let mut pts = vec![];
        for x in [0.6180339887f64, -1.3247179572, 2.7182818285] {
            let mut cy: Vec<f64> = (0..gg[0].len())
                .map(|j| {
                    let mut s = 0.0;
                    for i in 0..gg.len() {
                        s += to_f64(&gg[i][j]) * x.powf(i as f64);
                    }
                    s
                })
                .collect();
            while cy.last() == Some(&0.0) {
                cy.pop();
            }
            if cy.len() < 2 {
                continue;
            }
            let rev: Vec<f64> = cy.iter().rev().cloned().collect();
            pts.push((x, roots(&rev)));
        }
        Rc::new((tr, Some(pts)))
    };
    MEMO_CURVE_PTS.with(|m| {
        let mut m = m.borrow_mut();
        if m.len() > 20000 {
            m.clear();
        }
        m.insert(g.clone(), hit.clone());
    });
    hit
}



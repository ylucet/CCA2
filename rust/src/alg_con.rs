//! The result lattice of biconj (ports of alg_crv.py, assemble_alg_crv.py, envelope_min.py):
//! an AlgCon is a QuaCon, an AlgCrv or an EnvelopeMin.

use crate::alg_branch::{alg_branch, alg_branch_prep, fibre_vanishes, prep};
use crate::envelope::{envelope_eval, Env};
use crate::errors::R;
use crate::poly_n::PolyN;
use crate::qua_con::QuaCon;
use crate::rat_q::*;
use crate::triangle_face::TriFace;
use num_traits::Zero;
use std::sync::Arc;

#[derive(Debug)]
pub enum AlgCon {
    QuaCon(QuaCon),
    AlgCrv(AlgCrv),
    EnvelopeMin(EnvelopeMin),
}

impl AlgCon {
    pub fn kind(&self) -> &'static str {
        match self {
            AlgCon::QuaCon(_) => "QuaCon",
            AlgCon::AlgCrv(_) => "AlgCrv",
            AlgCon::EnvelopeMin(_) => "EnvelopeMin",
        }
    }
    pub fn nf(&self) -> Option<usize> {
        match self {
            AlgCon::QuaCon(g) => Some(g.nf),
            AlgCon::AlgCrv(g) => Some(g.nf),
            AlgCon::EnvelopeMin(_) => None,
        }
    }
    pub fn eval(&self, x: &[[f64; 2]]) -> R<Vec<f64>> {
        match self {
            AlgCon::QuaCon(g) => Ok(g.eval(x)?.0),
            AlgCon::AlgCrv(g) => Ok(g.eval(x)?.0),
            AlgCon::EnvelopeMin(g) => g.eval(x),
        }
    }
}

// ==========================================================================================
// subdivision faces (envelopeSubdivisionAlg / cutTriangles output)
// ==========================================================================================
#[derive(Clone, Debug, PartialEq)]
pub struct AlgRow {
    pub fp: Arc<PolyN>,
    pub fb: i64,
    pub num: Vec<Z>,
    pub den: Z,
    pub side: i64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IiRow {
    pub fp1: Arc<PolyN>,
    pub fb1: i64,
    pub fp2: Arc<PolyN>,
    pub fb2: i64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SubKind {
    Quad,
    Bridge,
    Triangle,
}

#[derive(Clone, Debug)]
pub struct SubFace {
    pub kind: SubKind,
    pub num: Vec<Z>,
    pub den: Z,
    pub fp: Option<Arc<PolyN>>,
    pub fb: i64,
    pub con: Vec<Vec<Z>>,
    pub alg: Vec<AlgRow>,
    pub ii: Vec<IiRow>,
    pub tri: Vec<[i64; 3]>,
    pub tf: Option<Arc<TriFace>>,
    pub tfs: Vec<Arc<TriFace>>,
}

// ==========================================================================================
// AlgCrv
// ==========================================================================================
#[derive(Clone, Debug)]
pub struct EcA {
    pub fp: Arc<PolyN>,
    pub fb: i64,
    pub num: Vec<Z>,
    pub den: Z,
}

#[derive(Clone, Debug, Default)]
pub struct AlgCrv {
    pub nf: usize,
    pub ecq: Vec<Vec<Z>>,
    pub fn_: Vec<Vec<Z>>,
    pub fd: Vec<Z>,
    pub fc: Vec<Vec<[i64; 4]>>,
    pub fp: Vec<Option<Arc<PolyN>>>,
    pub fb: Vec<i64>,
    pub eca: Vec<EcA>,
    pub eci: Vec<IiRow>,
    pub tri: Vec<Arc<TriFace>>,
    pub ft: Vec<usize>,
}

impl AlgCrv {
    pub fn eval(&self, x: &[[f64; 2]]) -> R<(Vec<f64>, Vec<usize>)> {
        let n = x.len();
        let mut val = vec![f64::INFINITY; n];
        let mut idx = vec![0usize; n];
        // A CONTACT TRIANGLE decides where it holds: its plane <s,x> - f*(s) is an affine
        // minorant of f** equal to f** on the whole closed triangle. The closed remnants that
        // cutting it out of a face leaves can meet it on a ZERO-AREA set (a segment of its edge)
        // carrying their pair's value, which for a discontinuous f is larger there (#108
        // residual: D-132, D-111) and would otherwise win the max.
        let mut tri = vec![false; n];
        for k in 0..self.nf {
            let m = self.face_closure_points(k, x)?;
            let pts: Vec<usize> = (0..n).filter(|&i| m[i]).collect();
            if pts.is_empty() {
                continue;
            }
            let xs: Vec<[f64; 2]> = pts.iter().map(|&i| x[i]).collect();
            let v = self.face_values(k, &xs)?;
            let is_tri = self.ft[k] > 0;
            for (&i, vi) in pts.iter().zip(v.iter()) {
                if vi.is_nan() {
                    continue;
                }
                if tri[i] && !is_tri {
                    continue; // a triangle already decides this point
                }
                if (is_tri && !tri[i]) || val[i].is_infinite() || *vi > val[i] {
                    val[i] = *vi;
                    idx[i] = k + 1;
                    tri[i] = tri[i] || is_tri;
                }
            }
        }
        Ok((val, idx))
    }

    /// face_points, with a row that is UNDEFINED at a point decided by the face's closure: an
    /// algebraic row whose whole fibre vanishes there (a BASE POINT, every z a root -- a 0/0 of a
    /// bridge between two atoms sharing a vertex) says nothing, and reading it as "fails" left a
    /// hole: D-153's domain corner (0,0), f** = -2, was +inf (#108 residual). There the face holds
    /// iff it holds at x nudged by 1e-7 in one of eight directions (envelope_eval's nudge); the
    /// value is still the face's own at x. See alg_crv.py `_faceClosurePoints`.
    fn face_closure_points(&self, k: usize, x: &[[f64; 2]]) -> R<Vec<bool>> {
        let (mut m, u) = self.face_points_und(k, x)?;
        let r2 = std::f64::consts::FRAC_1_SQRT_2;
        let dirs = [[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0], [r2, r2], [r2, -r2], [-r2, r2], [-r2, -r2]];
        for i in 0..x.len() {
            if m[i] || !u[i] {
                continue;
            }
            let h = 1e-7 * 1f64.max(x[i][0].abs()).max(x[i][1].abs());
            let y: Vec<[f64; 2]> = dirs.iter().map(|d| [x[i][0] + h * d[0], x[i][1] + h * d[1]]).collect();
            if self.face_points(k, &y)?.iter().any(|b| *b) {
                m[i] = true;
            }
        }
        Ok(m)
    }

    fn face_points(&self, k: usize, x: &[[f64; 2]]) -> R<Vec<bool>> {
        Ok(self.face_points_und(k, x)?.0)
    }

    /// (held, undecided): undecided = held but for rows undefined at a base point.
    fn face_points_und(&self, k: usize, x: &[[f64; 2]]) -> R<(Vec<bool>, Vec<bool>)> {
        let mut m = vec![true; x.len()];
        let mut u = vec![false; x.len()];
        for row in &self.fc[k] {
            let live: Vec<usize> = (0..x.len()).filter(|&i| m[i] || u[i]).collect();
            if live.is_empty() {
                return Ok((m, u));
            }
            let p: Vec<[f64; 2]> = live.iter().map(|&i| x[i]).collect();
            let tbl = row[3];
            let mut und = vec![false; p.len()];
            let ok: Vec<bool> = if tbl == 3 {
                let t = &self.tri[row[0] as usize - 1].t;
                let i0 = row[1] as usize;
                let j0 = i0 % 3 + 1;
                let o0 = 6 - i0 - j0;
                let (ti, tj, to) = (t[i0 - 1], t[j0 - 1], t[o0 - 1]);
                let e = [tj[0] - ti[0], tj[1] - ti[1]];
                let so = signf(e[0] * (to[1] - ti[1]) - e[1] * (to[0] - ti[0])) as f64;
                let ne = 1f64.max(e[0].hypot(e[1]));
                p.iter()
                    .map(|q| {
                        let cr = (e[0] * (q[1] - ti[1]) - e[1] * (q[0] - ti[0])) * so;
                        !((row[2] as f64) * cr < -1e-9 * ne * 1f64.max(q[0].hypot(q[1])))
                    })
                    .collect()
            } else if tbl == 2 {
                let r = &self.eci[row[0] as usize - 1];
                let v1 = alg_branch(&r.fp1, &p, r.fb1)?;
                let v2 = alg_branch(&r.fp2, &p, r.fb2)?;
                if v1.iter().chain(v2.iter()).any(|v| !v.is_finite()) {
                    let z1 = fibre_vanishes(&r.fp1, &p)?;
                    let z2 = fibre_vanishes(&r.fp2, &p)?;
                    for i in 0..p.len() {
                        und[i] = (!v1[i].is_finite() && z1[i]) || (!v2[i].is_finite() && z2[i]);
                    }
                }
                v1.iter().zip(v2.iter()).map(|(a, b)| a.is_finite() && b.is_finite() && !((row[1] as f64) * (a - b) < -1e-9)).collect()
            } else if tbl == 1 {
                let a = &self.eca[row[0] as usize - 1];
                let v = alg_branch(&a.fp, &p, a.fb)?;
                if v.iter().any(|vi| !vi.is_finite()) {
                    let z = fibre_vanishes(&a.fp, &p)?;
                    for i in 0..p.len() {
                        und[i] = !v[i].is_finite() && z[i];
                    }
                }
                p.iter().zip(v.iter()).map(|(q, vi)| {
                    let qi = eval_face(&a.num, &a.den, q[0], q[1]);
                    vi.is_finite() && !((row[1] as f64) * (qi - vi) < -1e-9)
                }).collect()
            } else {
                let c = &self.ecq[row[0] as usize - 1];
                let br = row[2];
                let v: Vec<f64> = if br != 0 {
                    let pe = pair_eval_prep(c, br)?;
                    p.iter().map(|q| pe.eval_m(q[0], q[1]).0).collect()
                } else {
                    p.iter().map(|q| eval_conic(c, q[0], q[1])).collect()
                };
                if row[1] == 0 {
                    v.iter().map(|vi| !(vi.abs() > 1e-9)).collect()
                } else {
                    v.iter().map(|vi| !((row[1] as f64) * vi < -1e-9)).collect()
                }
            };
            for ((&i, o), un) in live.iter().zip(ok).zip(und) {
                if m[i] {
                    m[i] = o;
                    u[i] = !o && un;
                } else {
                    u[i] = o || un; // an undecided point stays undecided, or fails
                }
            }
        }
        Ok((m, u))
    }

    fn face_values(&self, k: usize, x: &[[f64; 2]]) -> R<Vec<f64>> {
        if self.ft[k] > 0 {
            let t = &self.tri[self.ft[k] - 1];
            return Ok(x.iter().map(|q| q[0] * t.s[0] + q[1] * t.s[1] - t.c).collect());
        }
        if let Some(fp) = &self.fp[k] {
            let pr = prep(fp)?;
            return alg_branch_prep(&pr, x, self.fb[k]);
        }
        Ok(x.iter().map(|q| eval_face(&self.fn_[k], &self.fd[k], q[0], q[1])).collect())
    }
}

fn rc_eq(a: &Arc<PolyN>, b: &Arc<PolyN>) -> bool {
    Arc::ptr_eq(a, b) || **a == **b
}

pub fn assemble_alg_crv(sub: &[SubFace]) -> R<AlgCrv> {
    let mut g = AlgCrv::default();
    if let Some(f0) = sub.first() {
        g.tri = f0.tfs.clone();
    }
    let mut fc = vec![];
    let mut fnv = vec![];
    let mut fdv = vec![];
    let mut fpv = vec![];
    let mut fbv = vec![];
    let mut ftv = vec![];
    for f in sub {
        let mut rows: Vec<[i64; 4]> = vec![];
        for r in &f.con {
            let c = r[0..6].to_vec();
            if all_zero(&c[0..5]) {
                continue;
            }
            let hit = match g.ecq.iter().position(|e| *e == c) {
                Some(h) => h + 1,
                None => {
                    g.ecq.push(c);
                    g.ecq.len()
                }
            };
            let br = if r.len() >= 8 { crate::cells::side_i(&r[7]) } else { 0 };
            rows.push([hit as i64, crate::cells::side_i(&r[6]), br, 0]);
        }
        for a in &f.alg {
            let mut hit = 0;
            for (j, e) in g.eca.iter().enumerate() {
                if e.fb == a.fb && rc_eq(&e.fp, &a.fp) && e.num == a.num && e.den == a.den {
                    hit = j + 1;
                    break;
                }
            }
            if hit == 0 {
                g.eca.push(EcA { fp: a.fp.clone(), fb: a.fb, num: a.num.clone(), den: a.den.clone() });
                hit = g.eca.len();
            }
            rows.push([hit as i64, a.side, 0, 1]);
        }
        for r in &f.ii {
            let mut hit = 0;
            for (j, e) in g.eci.iter().enumerate() {
                if e.fb1 == r.fb1 && e.fb2 == r.fb2 && rc_eq(&e.fp1, &r.fp1) && rc_eq(&e.fp2, &r.fp2) {
                    hit = j + 1;
                    break;
                }
            }
            if hit == 0 {
                g.eci.push(r.clone());
                hit = g.eci.len();
            }
            rows.push([hit as i64, 1, 0, 2]);
        }
        for t in &f.tri {
            rows.push([t[0], t[1], t[2], 3]);
        }
        fc.push(rows);
        match f.kind {
            SubKind::Triangle => {
                fpv.push(None);
                fbv.push(0);
                fnv.push(vec![Z::zero(); 10]);
                fdv.push(zi(1));
                let tf = f.tf.as_ref().unwrap();
                ftv.push(g.tri.iter().position(|u| u.t == tf.t && u.s == tf.s).unwrap() + 1);
            }
            SubKind::Bridge => {
                fpv.push(f.fp.clone());
                fbv.push(f.fb);
                fnv.push(vec![Z::zero(); 10]);
                fdv.push(zi(1));
                ftv.push(0);
            }
            SubKind::Quad => {
                fpv.push(None);
                fbv.push(0);
                fnv.push(f.num.clone());
                fdv.push(f.den.clone());
                ftv.push(0);
            }
        }
    }
    let keys: Vec<Vec<[i64; 4]>> = fc.iter().map(|r: &Vec<[i64; 4]>| r.iter().cloned().collect::<std::collections::BTreeSet<_>>().into_iter().collect()).collect();
    let mut keep = vec![true; fc.len()];
    for k in 1..fc.len() {
        for j in 0..k {
            let same_fp = match (&fpv[j], &fpv[k]) {
                (None, None) => true,
                (Some(a), Some(b)) => rc_eq(a, b),
                _ => false,
            };
            if keep[j] && keys[j] == keys[k] && ftv[j] == ftv[k] && fbv[j] == fbv[k] && same_fp && fnv[j] == fnv[k] && fdv[j] == fdv[k] {
                keep[k] = false;
                break;
            }
        }
    }
    let sel = |i: usize| keep[i];
    g.fc = fc.into_iter().enumerate().filter(|(i, _)| sel(*i)).map(|(_, v)| v).collect();
    g.fn_ = fnv.into_iter().enumerate().filter(|(i, _)| sel(*i)).map(|(_, v)| v).collect();
    g.fd = fdv.into_iter().enumerate().filter(|(i, _)| sel(*i)).map(|(_, v)| v).collect();
    g.fp = fpv.into_iter().enumerate().filter(|(i, _)| sel(*i)).map(|(_, v)| v).collect();
    g.fb = fbv.into_iter().enumerate().filter(|(i, _)| sel(*i)).map(|(_, v)| v).collect();
    g.ft = ftv.into_iter().enumerate().filter(|(i, _)| sel(*i)).map(|(_, v)| v).collect();
    g.nf = g.fc.len();
    Ok(g)
}

// ==========================================================================================
// EnvelopeMin
// ==========================================================================================
#[derive(Clone, Debug)]
pub struct EnvTri {
    pub t: [[f64; 2]; 3],
    pub s: [f64; 2],
    pub c: f64,
    pub g: Vec<[f64; 2]>,
}

#[derive(Clone, Debug, Default)]
pub struct EnvelopeMin {
    pub pairs: Vec<Env>,
    pub tris: Vec<EnvTri>,
}

impl EnvelopeMin {
    pub fn eval(&self, x: &[[f64; 2]]) -> R<Vec<f64>> {
        let n = x.len();
        let mut val = vec![f64::INFINITY; n];
        for pr in &self.pairs {
            let v = envelope_eval(pr, x)?;
            for i in 0..n {
                val[i] = val[i].min(v[i]);
            }
        }
        for t in &self.tris {
            for i in 0..n {
                let inside = if !t.g.is_empty() { in_triangle_plus_cone(x[i], &t.t, &t.g) } else { in_triangle(x[i], &t.t) };
                if inside {
                    val[i] = val[i].min(x[i][0] * t.s[0] + x[i][1] * t.s[1] - t.c);
                }
            }
        }
        Ok(val)
    }
}

fn tri_tol(t: &[[f64; 2]; 3]) -> f64 {
    1e-9 * 1f64.max(t.iter().flatten().fold(0.0f64, |m, v| m.max(v.abs()))).powi(2)
}

fn in_triangle(x: [f64; 2], t: &[[f64; 2]; 3]) -> bool {
    let mut c = [0.0; 3];
    for i in 0..3 {
        let a = t[i];
        let b = t[(i + 1) % 3];
        c[i] = (b[0] - a[0]) * (x[1] - a[1]) - (b[1] - a[1]) * (x[0] - a[0]);
    }
    let tol = tri_tol(t);
    c.iter().all(|v| *v >= -tol) || c.iter().all(|v| *v <= tol)
}

fn in_triangle_plus_cone(x: [f64; 2], t: &[[f64; 2]; 3], g: &[[f64; 2]]) -> bool {
    let k = g.len();
    let tol = tri_tol(t);
    let mut a: Vec<Vec<f64>> = vec![];
    let mut b: Vec<f64> = vec![];
    for e in 0..3 {
        let p = t[e];
        let c = t[(e + 1) % 3];
        let o = t[(e + 2) % 3];
        let nrm = [-(c[1] - p[1]), c[0] - p[0]];
        let mut sg = signf(nrm[0] * (o[0] - p[0]) + nrm[1] * (o[1] - p[1])) as f64;
        if sg == 0.0 {
            sg = 1.0;
        }
        a.push(g.iter().map(|gg| sg * (nrm[0] * gg[0] + nrm[1] * gg[1])).collect());
        b.push(sg * (nrm[0] * (x[0] - p[0]) + nrm[1] * (x[1] - p[1])));
    }
    for j in 0..k {
        a.push((0..k).map(|i| if i == j { -1.0 } else { 0.0 }).collect());
        b.push(0.0);
    }
    let mut cand: Vec<Vec<f64>> = vec![];
    if k == 1 {
        for (ai, bi) in a.iter().zip(b.iter()) {
            if ai[0] != 0.0 {
                cand.push(vec![bi / ai[0]]);
            }
        }
        cand.push(vec![0.0]);
    } else {
        for p in 0..a.len() {
            for q in p + 1..a.len() {
                if k != 2 {
                    continue;
                }
                let m = [[a[p][0], a[p][1]], [a[q][0], a[q][1]]];
                let det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
                if det.abs() < 1e-14 {
                    continue;
                }
                cand.push(crate::conic_meet::solve2f(&m, [b[p], b[q]]).to_vec());
            }
        }
    }
    for c in cand {
        if a.iter().zip(b.iter()).all(|(ai, bi)| ai.iter().zip(c.iter()).map(|(x, y)| x * y).sum::<f64>() <= bi + tol) {
            return true;
        }
    }
    false
}

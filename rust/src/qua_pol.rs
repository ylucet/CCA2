//! QuaPol -- a piecewise QUAdratic function on a POLyhedral subdivision: the INPUT type
//! (port of qua_pol.py). Vertex / edge / face indices are 1-based in the data, 0 = none.

use crate::bail;
use crate::errors::R;
use crate::rat_q::*;
use num_traits::Zero;

#[derive(Clone, Debug)]
pub struct Dom {
    pub dim: usize,
    pub is_convex: bool,
}

#[derive(Clone, Debug)]
pub struct QuaPol {
    pub nv: usize,
    pub ne: usize,
    pub nf: usize,
    pub vf: Vec<[f64; 2]>,      // the caller's vertices as doubles (legacy dom tests)
    pub ff: Vec<Vec<f64>>,      // the caller's function rows as doubles (degree())
    pub e: Vec<[usize; 3]>,
    pub f: Vec<[usize; 2]>,
    pub fn_: Vec<Vec<Z>>,
    pub fd: Vec<Z>,
    pub vn: Vec<Vec<Z>>,
    pub vd: Vec<Z>,
    pub exact: bool,
    pub f_is_convex: Option<bool>,
    pub dom: Dom,
}

fn num_f(x: &Num) -> f64 {
    match x {
        Num::F(f) => *f,
        Num::Q(q) => q_to_f64(q),
    }
}

impl QuaPol {
    /// QuaPol(f): one quadratic on the whole plane.
    pub fn whole(f: Vec<Num>) -> R<QuaPol> {
        Self::build(vec![], vec![], vec![f], vec![], true)
    }

    /// QuaPol(f) with several rows (one-argument form: nf is 1 whatever f holds).
    pub fn whole_rows(f: Vec<Vec<Num>>) -> R<QuaPol> {
        Self::build(vec![], vec![], f, vec![], true)
    }

    /// QuaPol(V, E, f, F).
    pub fn new(v: Vec<Vec<Num>>, e: Vec<[usize; 3]>, f: Vec<Vec<Num>>, ff: Vec<[usize; 2]>) -> R<QuaPol> {
        Self::build(v, e, f, ff, false)
    }

    /// Convenience: exact rational data.
    pub fn new_q(v: Vec<[Q; 2]>, e: Vec<[usize; 3]>, f: Vec<Vec<Q>>, ff: Vec<[usize; 2]>) -> R<QuaPol> {
        let v = v.into_iter().map(|r| r.iter().map(|x| Num::Q(x.clone())).collect()).collect();
        let f = f.into_iter().map(|r| r.into_iter().map(Num::Q).collect()).collect();
        Self::new(v, e, f, ff)
    }

    fn build(v: Vec<Vec<Num>>, e: Vec<[usize; 3]>, f: Vec<Vec<Num>>, fe: Vec<[usize; 2]>, single: bool) -> R<QuaPol> {
        let f: Vec<Vec<Num>> = f
            .into_iter()
            .map(|r| {
                let mut out = vec![Num::Q(qi(0)); 10usize.saturating_sub(r.len())];
                out.extend(r);
                out
            })
            .collect();
        let (nv, ne, mut nf);
        if single {
            nv = 0;
            ne = 0;
            nf = 1;
        } else {
            nv = v.len();
            ne = e.len();
            nf = f.len();
            let ok = [v.iter().all(|r| r.len() == 2), ne == fe.len(), nf > 0, ne == fe.len()];
            if !ok.iter().all(|x| *x) {
                bail!("QuaPol:size", "Invalid matrix size in QuaPol constructor.");
            }
            if (nv, ne) == (1, 0) || (nv, ne) == (2, 1) {
                nf = 0;
            }
            if !fe.is_empty() && fe.iter().all(|r| r[0] == 0 && r[1] == 0) {
                nf = 0;
            }
            if nf > 0 {
                let mx = fe.iter().map(|r| r[0].max(r[1])).max().unwrap_or(0);
                if ne > 1 && mx != nf {
                    bail!("QuaPol:badF", "max index in F must be equal to the number of faces");
                }
                if nv > 1 && e.iter().map(|r| r[0].max(r[1])).max().unwrap_or(0) != nv {
                    bail!("QuaPol:badE", "max index in E must be equal to size(V,1)");
                }
            }
        }
        let vf: Vec<[f64; 2]> = v.iter().map(|r| [num_f(&r[0]), num_f(&r[1])]).collect();
        let ffl: Vec<Vec<f64>> = f.iter().map(|r| r.iter().map(num_f).collect()).collect();
        // the exact reading
        let mut exact = true;
        let mut fn_ = vec![];
        let mut fd = vec![];
        let mut vn = vec![];
        let mut vdd = vec![];
        let res: R<()> = (|| {
            for row in &f {
                let (n, d) = from_double(row)?;
                fn_.push(n);
                fd.push(d);
            }
            for row in &v {
                let (n, d) = from_double(row)?;
                vn.push(n);
                vdd.push(d);
            }
            Ok(())
        })();
        if let Err(e) = res {
            if !["ratQ:notExact", "ratQ:overflow", "ratQ:notFinite"].contains(&e.id.as_str()) {
                return Err(e);
            }
            exact = false;
        }
        let mut q = QuaPol {
            nv,
            ne,
            nf,
            vf,
            ff: ffl,
            e,
            f: fe,
            fn_,
            fd,
            vn,
            vd: vdd,
            exact,
            f_is_convex: None,
            dom: Dom { dim: 0, is_convex: false },
        };
        q.dom = q.create_dom()?;
        Ok(q)
    }

    pub fn is_exact(&self) -> bool {
        self.exact
    }

    pub fn assert_exact(&self) -> R<()> {
        if !self.exact {
            bail!("PLQ:QuaPol:notExact", "this QuaPol does not carry exact coefficients, so it cannot be used on the exact path.");
        }
        Ok(())
    }

    /// Face k (1-based), exactly.
    pub fn face_q(&self, k: usize) -> (Vec<Z>, Z) {
        (self.fn_[k - 1].clone(), self.fd[k - 1].clone())
    }

    pub fn vertex_q(&self, i: usize) -> (Vec<Z>, Z) {
        (self.vn[i - 1].clone(), self.vd[i - 1].clone())
    }

    /// Vertex i (1-based) as exact rationals.
    pub fn vx(&self, i: usize) -> [Q; 2] {
        [Q::new(self.vn[i - 1][0].clone(), self.vd[i - 1].clone()), Q::new(self.vn[i - 1][1].clone(), self.vd[i - 1].clone())]
    }

    /// Function row k (1-based) as exact rationals.
    pub fn fx(&self, k: usize) -> Vec<Q> {
        self.fn_[k - 1].iter().map(|x| Q::new(x.clone(), self.fd[k - 1].clone())).collect()
    }

    fn create_dom(&self) -> R<Dom> {
        if self.nv == 0 {
            if self.nf != 1 {
                bail!("QuaPol:dom", "a function with no vertex must have exactly 1 face");
            }
            return Ok(Dom { dim: 2, is_convex: true });
        }
        if self.nv == 1 {
            return Ok(Dom { dim: 0, is_convex: true });
        }
        if self.nf == 0 {
            let conv = (self.nv == 2 && self.ne == 1) || is_collinear(&self.vf, (2f64.powi(-52)).sqrt());
            return Ok(Dom { dim: 1, is_convex: conv });
        }
        if !self.f.iter().any(|r| r[0] == 0 || r[1] == 0) {
            return Ok(Dom { dim: 2, is_convex: true });
        }
        let c = domain_is_convex_d(&self.vf, &self.e, &self.f)?;
        Ok(Dom { dim: 2, is_convex: c })
    }

    pub fn is_dom_bounded(&self) -> bool {
        if self.ne == 0 {
            return self.nv == 1;
        }
        self.e.iter().all(|r| r[2] != 0)
    }

    pub fn degree(&self) -> usize {
        if self.ff.is_empty() {
            return 0;
        }
        let tol = (2f64.powi(-52)).sqrt();
        let nrm = |cols: std::ops::Range<usize>| -> f64 {
            let mut s = 0.0;
            for r in &self.ff {
                for c in cols.clone() {
                    s += r[c] * r[c];
                }
            }
            s.sqrt()
        };
        if nrm(0..4) > tol {
            3
        } else if nrm(4..7) > tol {
            2
        } else if nrm(7..9) > tol {
            1
        } else {
            0
        }
    }

    pub fn assert_operable(&self) -> R<()> {
        if self.degree() > 2 {
            bail!("PLQ:op:unsupportedType", "Convex-analysis operators require a quadratic numerator (degree<=2); got degree {}.", self.degree());
        }
        Ok(())
    }

    /// f is CONTINUOUS across every shared edge -- exact, three points per edge.
    ///
    /// d = f_L - f_R restricted to the edge LINE is a univariate quadratic, and it vanishes on
    /// the whole edge iff it vanishes at three points of it: both endpoints and the midpoint
    /// (the same condition (b) three-point test as MATLAB isConvexQ.m). A degenerate edge
    /// (both endpoints at one vertex) leaves one point, and the midpoint formula below
    /// evaluates to that point, so the single test that remains decides d on it.
    pub fn assert_continuous(&self) -> R<()> {
        self.assert_exact()?;
        if self.nf < 2 {
            return Ok(());
        }
        for (j, row) in self.f.iter().enumerate() {
            let (il, ir) = (row[0], row[1]);
            if il == 0 || ir == 0 || il == ir {
                continue;
            }
            // d = f_il - f_ir as exact rationals, in the weighted basis: the value is
            // (d5 x^2 + 2 d6 xy + d7 y^2)/2 + d8 x + d9 y + d10 (QuaPol's header; isConvexQ.m
            // evalNum multiplies this by the positive 2*D^2).
            let d: Vec<Q> = (0..10)
                .map(|k| Q::new(self.fn_[il - 1][k].clone(), self.fd[il - 1].clone())
                    - Q::new(self.fn_[ir - 1][k].clone(), self.fd[ir - 1].clone()))
                .collect();
            let half = Q::new(zi(1), zi(2));
            let eval = |p: [Q; 2]| -> Q {
                d[4].clone() * p[0].clone() * p[0].clone() * half.clone()
                    + d[5].clone() * p[0].clone() * p[1].clone()
                    + d[6].clone() * p[1].clone() * p[1].clone() * half.clone()
                    + d[7].clone() * p[0].clone()
                    + d[8].clone() * p[1].clone()
                    + d[9].clone()
            };
            let e = self.e[j];
            let p1 = self.vx(e[0]);
            let p2 = self.vx(e[1]);
            let pm = [p1[0].clone() * half.clone() + p2[0].clone() * half.clone(),
                      p1[1].clone() * half.clone() + p2[1].clone() * half.clone()];
            if !eval(p1).is_zero() || !eval(p2).is_zero() || !eval(pm).is_zero() {
                bail!("PLQ:notContinuous",
                    "f is discontinuous across edge {}: faces {} and {} do not agree on the edge, and the conjugate and biconjugate are implemented only for continuous inputs.", j + 1, il, ir);
            }
        }
        Ok(())
    }

    pub fn n_rows(&self) -> usize {
        self.fn_.len()
    }
}

pub fn is_collinear(v: &[[f64; 2]], tol: f64) -> bool {
    if v.len() < 3 {
        return true;
    }
    let p0 = v[0];
    let mut d = None;
    for p in &v[1..] {
        let u = [p[0] - p0[0], p[1] - p0[1]];
        if u[0].hypot(u[1]) > tol {
            d = Some(u);
            break;
        }
    }
    let d = match d {
        None => return true,
        Some(d) => d,
    };
    let n = d[0].hypot(d[1]);
    for p in &v[1..] {
        let u = [p[0] - p0[0], p[1] - p0[1]];
        if (d[0] * u[1] - d[1] * u[0]).abs() / n > tol * 1f64.max(u[0].hypot(u[1])) {
            return false;
        }
    }
    true
}

/// Is the DOMAIN of a polyhedral subdivision convex? (port of domain_is_convex_d.py; doubles)
pub fn domain_is_convex_d(v: &[[f64; 2]], e: &[[usize; 3]], f: &[[usize; 2]]) -> R<bool> {
    let tol = (2f64.powi(-52)).sqrt();
    if f.is_empty() || e.is_empty() || v.is_empty() {
        return Ok(true);
    }
    if f.iter().all(|r| r[0] == 0 && r[1] == 0) {
        bail!("PLQ:domainIsConvexD:thinSubdivision", "every edge has +infinity on both sides.");
    }
    let mut dirs = vec![];
    for ed in e {
        if ed[2] == 0 {
            let u = [v[ed[1] - 1][0] - v[ed[0] - 1][0], v[ed[1] - 1][1] - v[ed[0] - 1][1]];
            let n = u[0].hypot(u[1]);
            if n > 0.0 {
                dirs.push([u[0] / n, u[1] / n]);
            }
        }
    }
    for j in 0..f.len() {
        let (il, ir) = (f[j][0], f[j][1]);
        if (il != 0 && ir != 0) || (il == 0 && ir == 0) {
            continue;
        }
        let v1 = v[e[j][0] - 1];
        let v2 = v[e[j][1] - 1];
        let a = [-(v2[1] - v1[1]), v2[0] - v1[0]];
        let na = a[0].hypot(a[1]);
        if na == 0.0 {
            continue;
        }
        let c = v1[0] * v2[1] - v2[0] * v1[1];
        let s = if il == 0 { -1.0 } else { 1.0 };
        let a = [s * a[0] / na, s * a[1] / na];
        let c = s * c / na;
        for p in v {
            if a[0] * p[0] + a[1] * p[1] + c < -tol * 1f64.max(p[0].hypot(p[1])) {
                return Ok(false);
            }
        }
        for d in &dirs {
            if a[0] * d[0] + a[1] * d[1] < -tol {
                return Ok(false);
            }
        }
    }
    Ok(true)
}

/// Helper: all-zero test for exact rows.
pub fn row_is_zero(r: &[Z]) -> bool {
    r.iter().all(|x| x.is_zero())
}

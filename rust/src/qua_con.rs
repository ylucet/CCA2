//! QuaCon -- a quadratic function on a CONIC subdivision, stored EXACTLY (port of qua_con.py).
//! FC[k] rows are [edgeIndex (1-based), side (-1/0/+1), branch (0/1/2)].

use crate::bail;
use crate::conic_meet::conic_meet;
use crate::errors::R;
use crate::rat_q::*;

#[derive(Clone, Debug, Default)]
pub struct QuaCon {
    pub nv: usize,
    pub ne: usize,
    pub nf: usize,
    pub vname: Vec<[usize; 3]>,
    pub e: Vec<[usize; 3]>,
    pub ecq: Vec<Vec<Z>>,
    pub fn_: Vec<Vec<Z>>,
    pub fd: Vec<Z>,
    pub f: Vec<[usize; 2]>,
    pub fc: Vec<Vec<[i64; 3]>>,
}

impl QuaCon {
    pub fn empty() -> QuaCon {
        QuaCon::default()
    }

    pub fn new(
        vname: Vec<[usize; 3]>,
        e: Vec<[usize; 3]>,
        ecq: Vec<Vec<Z>>,
        fn_: Vec<Vec<Z>>,
        fd: Vec<Z>,
        f: Vec<[usize; 2]>,
        fc: Vec<Vec<[i64; 3]>>,
    ) -> R<QuaCon> {
        let nv = vname.len();
        let ne = e.len();
        let nf = fn_.len();
        if ecq.len() != ne || f.len() != ne {
            bail!("QuaCon:size", "E must be ne x 3 and EcQ / F must have one row per edge.");
        }
        if fd.len() != nf || fc.len() != nf {
            bail!("QuaCon:size", "fD and FC need one entry per face.");
        }
        if nv > 0 && ne > 0 {
            let mx = e.iter().map(|r| r[0].max(r[1])).max().unwrap();
            if mx > nv {
                bail!("QuaCon:badEdge", "E names vertex {} but there are only {} vertices.", mx, nv);
            }
        }
        let mut fn2 = Vec::with_capacity(nf);
        let mut fd2 = Vec::with_capacity(nf);
        for k in 0..nf {
            let (n, d) = canon(&fn_[k], &fd[k])?;
            if !all_zero(&n[0..4]) {
                bail!("QuaCon:cubicFace", "face {} has a cubic numerator.", k + 1);
            }
            fn2.push(n);
            fd2.push(d);
        }
        let mut ecq2 = Vec::with_capacity(ne);
        for c in &ecq {
            ecq2.push(conic(c)?);
        }
        for (i, v) in vname.iter().enumerate() {
            let (a, b, r) = (v[0], v[1], v[2]);
            if a < 1 || b < 1 || a > ne || b > ne {
                bail!("QuaCon:badVertexName", "vertex {} names edges {} and {}, but there are {} edges.", i + 1, a, b, ne);
            }
            if a == b {
                bail!("QuaCon:badVertexName", "vertex {} names edge {} twice.", i + 1, a);
            }
            let m = conic_meet(&ecq2[a - 1], &ecq2[b - 1])?;
            if m.1.degenerate {
                bail!("QuaCon:degenerateVertex", "vertex {} names edges {} and {}, whose conics share a component.", i + 1, a, b);
            }
            if r < 1 || r > m.0.len() {
                bail!("QuaCon:badRootIndex", "vertex {} asks for real intersection {} of edges {} and {}, which meet in {} real points.", i + 1, r, a, b, m.0.len());
            }
        }
        for (k, ck) in fc.iter().enumerate() {
            for (ri, r) in ck.iter().enumerate() {
                if r[0] < 1 || r[0] as usize > ne {
                    bail!("QuaCon:badConstraint", "FC{{{}}} names an edge outside 1..{}.", k + 1, ne);
                }
                if !(-1..=1).contains(&r[1]) {
                    bail!("QuaCon:badConstraint", "FC{{{}}} carries a side that is not -1, 0 or +1.", k + 1);
                }
                if !(0..=2).contains(&r[2]) {
                    bail!("QuaCon:badConstraint", "FC{{{}}} carries a branch that is not 0, 1 or 2.", k + 1);
                }
                if r[2] != 0 {
                    if let Err(me) = pair_line_data_m(&ecq2[r[0] as usize - 1]) {
                        bail!("QuaCon:badConstraint", "FC{{{}}} row {} names branch {} of edge {}, but that conic is not two distinct crossing lines ({}).", k + 1, ri + 1, r[2], r[0], me.id);
                    }
                }
            }
        }
        Ok(QuaCon { nv, ne, nf, vname, e, ecq: ecq2, fn_: fn2, fd: fd2, f, fc })
    }

    /// Realise the named vertices as double coordinates, on demand.
    pub fn vertex_coords(&self) -> R<Vec<[f64; 2]>> {
        let mut v = vec![];
        for n in &self.vname {
            let m = conic_meet(&self.ecq[n[0] - 1], &self.ecq[n[1] - 1])?;
            v.push(m.0[n[2] - 1]);
        }
        Ok(v)
    }

    /// Value at each point of X and the 1-based face that supplied it (0 where +inf).
    pub fn eval(&self, x: &[[f64; 2]]) -> R<(Vec<f64>, Vec<usize>)> {
        let n = x.len();
        let mut val = vec![f64::INFINITY; n];
        let mut idx = vec![0usize; n];
        if n == 0 {
            return Ok((val, idx));
        }
        let mut best = vec![f64::INFINITY; n];
        let xmax: Vec<f64> = x.iter().map(|p| 1f64.max(p[0].abs()).max(p[1].abs())).collect();
        let tol = (2f64.powi(-52)).sqrt();
        let ecf: Vec<Vec<f64>> = self.ecq.iter().map(|c| c.iter().map(to_f64).collect()).collect();
        for k in 0..self.nf {
            let mut viol = vec![f64::NEG_INFINITY; n];
            for r in &self.fc[k] {
                let (e, side, br) = (r[0] as usize, r[1], r[2]);
                let c = &ecf[e - 1];
                if br == 0 {
                    let cm = c.iter().fold(0.0f64, |m, v| m.max(v.abs()));
                    for i in 0..n {
                        let q = eval_conic_f(c, x[i][0], x[i][1]);
                        let s = 1f64.max(cm * xmax[i] * xmax[i]);
                        viol[i] = if side == 0 { viol[i].max(q.abs() / s) } else { viol[i].max(-(side as f64) * q / s) };
                    }
                } else {
                    let pe = pair_eval_prep(&self.ecq[e - 1], br)?;
                    for i in 0..n {
                        let q = pe.eval_m(x[i][0], x[i][1]).0;
                        let s = xmax[i];
                        viol[i] = if side == 0 { viol[i].max(q.abs() / s) } else { viol[i].max(-(side as f64) * q / s) };
                    }
                }
            }
            let nf: Vec<f64> = self.fn_[k].iter().map(to_f64).collect();
            let df = to_f64(&self.fd[k]);
            for i in 0..n {
                if viol[i] <= tol && viol[i] < best[i] {
                    val[i] = eval_poly_f(&nf, x[i][0], x[i][1]) / df;
                    idx[i] = k + 1;
                    best[i] = viol[i];
                }
            }
        }
        Ok((val, idx))
    }
}

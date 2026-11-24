//! is_convex_q -- is an exact QuaPol convex? Decided exactly (port of is_convex_q.py).

use crate::bail;
use crate::errors::R;
use crate::qua_pol::QuaPol;
use crate::rat_q::*;
use num_traits::{Signed, Zero};

pub fn is_convex_q(obj: &QuaPol) -> R<(bool, String)> {
    obj.assert_operable()?;
    obj.assert_exact()?;
    let nf = obj.fn_.len();
    for k in 0..nf {
        let (h, _) = hess_q(&obj.fn_[k], &obj.fd[k])?;
        if !is_psd2(&h) {
            return Ok((false, format!("face {} has a Hessian that is not positive semidefinite", k + 1)));
        }
    }
    if obj.nv == 0 {
        return Ok((true, String::new()));
    }
    if obj.dom.dim < 2 {
        bail!("PLQ:isConvexQ:thinDomain", "isConvexQ is for a 2-dimensional domain; this one has dimension {}.", obj.dom.dim);
    }
    let (ok, jbad) = domain_is_convex(obj)?;
    if !ok {
        return Ok((false, format!("the domain is not convex: it leaves the half-plane of boundary edge {}", jbad)));
    }
    for j in 0..obj.f.len() {
        let (il, ir) = (obj.f[j][0], obj.f[j][1]);
        if il == 0 || ir == 0 || il == ir {
            continue;
        }
        let (dn, dd) = sub(&obj.fn_[il - 1], &obj.fd[il - 1], &obj.fn_[ir - 1], &obj.fd[ir - 1])?;
        let (p1, d1) = vertex_of(obj, obj.e[j][0])?;
        let (p2, d2) = vertex_of(obj, obj.e[j][1])?;
        let a = vec![-(&p2[1] * &d1 - &p1[1] * &d2), &p2[0] * &d1 - &p1[0] * &d2];
        if a[0].is_zero() && a[1].is_zero() {
            bail!("PLQ:isConvexQ:degenerateEdge", "edge {} has both endpoints at the same vertex, so it has no normal.", j + 1);
        }
        let pm = vec![&p1[0] * &d2 + &p2[0] * &d1, &p1[1] * &d2 + &p2[1] * &d1];
        let dm = zi(2) * &d1 * &d2;
        if !eval_num(&dn, &p1, &d1).is_zero() || !eval_num(&dn, &p2, &d2).is_zero() || !eval_num(&dn, &pm, &dm).is_zero() {
            return Ok((false, format!("f is discontinuous across edge {}", j + 1)));
        }
        let (g, _) = grad_q(&dn, &dd)?;
        let s1 = delta_num(&g, &a, &p1, &d1);
        let s2 = delta_num(&g, &a, &p2, &d2);
        if s1.is_negative() || s2.is_negative() {
            return Ok((false, format!("the subdifferential is empty on edge {}", j + 1)));
        }
        if obj.e[j][2] == 0 && &s1 * &d2 > &s2 * &d1 {
            return Ok((false, format!("the subdifferential empties out along the unbounded edge {}", j + 1)));
        }
    }
    Ok((true, String::new()))
}

fn vertex_of(obj: &QuaPol, i: usize) -> R<(Vec<Z>, Z)> {
    let p = obj.vn[i - 1].clone();
    let d = obj.vd[i - 1].clone();
    if !d.is_positive() {
        bail!("PLQ:isConvexQ:badVertex", "vertex {} has a nonpositive denominator.", i);
    }
    Ok((p, d))
}

fn eval_num(n: &[Z], p: &[Z], d: &Z) -> Z {
    &n[4] * &p[0] * &p[0] + zi(2) * &n[5] * &p[0] * &p[1] + &n[6] * &p[1] * &p[1] + zi(2) * &n[7] * &p[0] * d + zi(2) * &n[8] * &p[1] * d + zi(2) * &n[9] * d * d
}

fn delta_num(g: &[Vec<Z>], a: &[Z], p: &[Z], d: &Z) -> Z {
    let gv = mv(g, &[p[0].clone(), p[1].clone(), d.clone()]);
    &gv[0] * &a[0] + &gv[1] * &a[1]
}

fn domain_is_convex(obj: &QuaPol) -> R<(bool, usize)> {
    let mut rays = vec![];
    for e in &obj.e {
        if e[2] == 0 {
            let (q1, c1) = vertex_of(obj, e[0])?;
            let (q2, c2) = vertex_of(obj, e[1])?;
            rays.push(vec![&q2[0] * &c1 - &q1[0] * &c2, &q2[1] * &c1 - &q1[1] * &c2]);
        }
    }
    for j in 0..obj.f.len() {
        let (il, ir) = (obj.f[j][0], obj.f[j][1]);
        if il != 0 && ir != 0 {
            continue;
        }
        if il == 0 && ir == 0 {
            bail!("PLQ:isConvexQ:danglingEdge", "edge {} has +infinity on both sides.", j + 1);
        }
        let (p1, d1) = vertex_of(obj, obj.e[j][0])?;
        let (p2, d2) = vertex_of(obj, obj.e[j][1])?;
        let a = [-(&p2[1] * &d1 - &p1[1] * &d2), &p2[0] * &d1 - &p1[0] * &d2];
        let c = &p1[0] * &p2[1] - &p2[0] * &p1[1];
        let sgn = if il == 0 { zi(-1) } else { zi(1) };
        for i in 1..=obj.vn.len() {
            let (q, cc) = vertex_of(obj, i)?;
            if (&sgn * (&a[0] * &q[0] + &a[1] * &q[1] + &c * &cc)).is_negative() {
                return Ok((false, j + 1));
            }
        }
        for r in &rays {
            if (&sgn * (&a[0] * &r[0] + &a[1] * &r[1])).is_negative() {
                return Ok((false, j + 1));
            }
        }
    }
    Ok((true, 0))
}

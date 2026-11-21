//! Cell interior tests and the QuaCon assembly (ports of cell_has_interior.py, cell_interior.py,
//! cell_sample_point.py, compact_cells.py, merge_adjacent_cells.py, assemble_qua_con_cells.py).

use crate::bail;
use crate::conic_meet::{conic_meet, solve2f};
use crate::errors::R;
use crate::numeric::{det3, linspace, round_half_away, solve3};
use crate::probe_data::{PROBE_X, PROBE_Y};
use crate::qua_con::QuaCon;
use crate::rat_q::*;
use num_traits::{Signed, Zero};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::sync::Arc;

// ==========================================================================================
// cell_has_interior
// ==========================================================================================
fn conic_f(c: &[Z]) -> ([f64; 6], f64) {
    let f = [to_f64(&c[0]), to_f64(&c[1]), to_f64(&c[2]), to_f64(&c[3]), to_f64(&c[4]), to_f64(&c[5])];
    let cm = c[0..6].iter().map(|x| x.abs()).max().unwrap();
    (f, to_f64(&cm))
}

/// Do all the cell's sign conditions hold at the double point x, to a 1e-9 relative band?
pub fn satisfies_all_double(ecq: &[Vec<Z>], rows: &[(usize, i64)], x: [f64; 2]) -> bool {
    for &(e, s) in rows {
        let (c, cm) = conic_f(&ecq[e - 1]);
        let v = c[0] * x[0].powi(2) + c[1] * x[0] * x[1] + c[2] * x[1].powi(2) + c[3] * x[0] + c[4] * x[1] + c[5];
        let scale = 1f64.max(cm * 1f64.max(x[0].abs().max(x[1].abs())).powi(2));
        if s == 0 {
            if v.abs() > 1e-9 * scale {
                return false;
            }
        } else if (s as f64) * v < -1e-9 * scale {
            return false;
        }
    }
    true
}

fn neg_semidefinite(c: &[Z]) -> bool {
    let m = [
        [-zi(2) * &c[0], -c[1].clone(), -c[3].clone()],
        [-c[1].clone(), -zi(2) * &c[2], -c[4].clone()],
        [-c[3].clone(), -c[4].clone(), -zi(2) * &c[5]],
    ];
    let lim = zi(1 << 16);
    if m.iter().flatten().any(|x| x.abs() > lim) || m.iter().flatten().all(|x| x.is_zero()) {
        return false;
    }
    if (0..3).any(|i| m[i][i].is_negative()) {
        return false;
    }
    for (a, b) in [(0usize, 1usize), (0, 2), (1, 2)] {
        if (&m[a][a] * &m[b][b] - &m[a][b] * &m[a][b]).is_negative() {
            return false;
        }
    }
    let d = &m[0][0] * &m[1][1] * &m[2][2] + zi(2) * &m[0][1] * &m[1][2] * &m[0][2]
        - &m[0][0] * &m[1][2] * &m[1][2]
        - &m[1][1] * &m[0][2] * &m[0][2]
        - &m[2][2] * &m[0][1] * &m[0][1];
    !d.is_negative()
}

/// Does a cell cut out by conic sign conditions have interior? rows: (index into ecq, side).
pub fn cell_has_interior(ecq: &[Vec<Z>], rows: &[(usize, i64)]) -> R<bool> {
    if rows.is_empty() {
        return Ok(true);
    }
    let is_line: Vec<bool> = rows.iter().map(|r| all_zero(&ecq[r.0 - 1][0..3])).collect();
    let lin_rows = |only: bool| -> Vec<Vec<Z>> {
        rows.iter()
            .zip(is_line.iter())
            .filter(|(_, l)| !only || **l)
            .map(|(r, _)| ecq[r.0 - 1][3..6].iter().map(|x| x * zi(r.1)).collect())
            .collect()
    };
    if is_line.iter().all(|x| *x) {
        return Ok(feasible2(&lin_rows(false), true));
    }
    if is_line.iter().any(|x| *x) && !feasible2(&lin_rows(true), true) {
        return Ok(false);
    }
    for (r, l) in rows.iter().zip(is_line.iter()) {
        if r.1 != 0 && !l {
            let c: Vec<Z> = ecq[r.0 - 1][0..6].iter().map(|x| x * zi(r.1)).collect();
            if neg_semidefinite(&c) {
                return Ok(false);
            }
        }
    }
    let idx: Vec<usize> = rows.iter().map(|r| r.0).collect::<BTreeSet<_>>().into_iter().collect();
    for a in 0..idx.len() {
        for b in a + 1..idx.len() {
            let m = conic_meet(&ecq[idx[a] - 1], &ecq[idx[b] - 1])?;
            if m.1.degenerate || m.0.is_empty() {
                continue;
            }
            for p in &m.0 {
                if satisfies_all_double(ecq, rows, *p) {
                    return Ok(true);
                }
            }
        }
    }
    if probe_interior(ecq, rows)? {
        return Ok(true);
    }
    cad_interior(ecq, rows)
}

/// EXACT and COMPLETE decision of "does the open cell {side_i * q_i > 0} have a point?", run when the
/// witness search finds nothing. Mirrors python/cca2/cell_has_interior.py `_cadInterior`; the
/// completeness argument is in that docstring (cylindrical decomposition over Q: shear so every
/// non-line conic has a y^2 term; critical abscissae = discriminants, vertical-line roots, pairwise
/// resultants, and for conics with a zero resultant the crossing of their non-shared lines; one
/// rational sample column per interval between them, one rational point per sector; exact signs).
/// Equalities (side 0) have no open cell: false.
fn cad_interior(ecq: &[Vec<Z>], rows: &[(usize, i64)]) -> R<bool> {
    if rows.iter().any(|r| r.1 == 0) {
        return Ok(false);
    }
    let idx: Vec<usize> = rows.iter().map(|r| r.0).collect::<BTreeSet<_>>().into_iter().collect();
    let mut sg: Vec<i64> = vec![];
    for &i in &idx {
        let sides: BTreeSet<i64> = rows.iter().filter(|r| r.0 == i).map(|r| r.1).collect();
        if sides.len() > 1 {
            return Ok(false); // q > 0 and q < 0 on the same conic
        }
        sg.push(*sides.iter().next().unwrap());
    }
    let cs0: Vec<Vec<Z>> = idx.iter().map(|&i| ecq[i - 1][0..6].to_vec()).collect();
    let mut k = 0i64;
    while cs0.iter().any(|c| {
        (!c[0].is_zero() || !c[1].is_zero() || !c[2].is_zero()) && (&c[0] * zi(k * k) + &c[1] * zi(k) + &c[2]).is_zero()
    }) {
        k += 1;
    }
    let cs: Vec<Vec<Z>> = cs0
        .iter()
        .map(|c| vec![c[0].clone(), zi(2 * k) * &c[0] + &c[1], &c[0] * zi(k * k) + &c[1] * zi(k) + &c[2], c[3].clone(), &c[3] * zi(k) + &c[4], c[5].clone()])
        .collect();
    let (xpolys, special) = cad_critical_polys(&cs)?;
    let xs = cad_samples(&xpolys, &[])?;
    for x0 in &xs {
        let (ypolys, yspecial) = match cad_y_polys(&cs, &sg, x0, &special) {
            Some(v) => v,
            None => continue,
        };
        let ys = cad_samples(&yspecial, &ypolys)?;
        for y0 in &ys {
            let mut ok = true;
            for (c, s) in cs.iter().zip(sg.iter()) {
                let v = cad_eval(c, x0, y0);
                if (*s > 0 && !v.is_positive()) || (*s < 0 && !v.is_negative()) {
                    ok = false;
                    break;
                }
            }
            if ok {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn cad_eval(c: &[Z], x: &Q, y: &Q) -> Q {
    qz(&c[0]) * x * x + qz(&c[1]) * x * y + qz(&c[2]) * y * y + qz(&c[3]) * x + qz(&c[4]) * y + qz(&c[5])
}

fn padd(a: &[Z], b: &[Z], sb: i64) -> Vec<Z> {
    let nb: Vec<Z> = b.iter().map(|x| x * zi(sb)).collect();
    poly_trim(&poly_add1(a, &nb))
}

fn pmul(a: &[Z], b: &[Z]) -> Vec<Z> {
    poly_trim(&poly_mul1(a, b))
}

fn cad_critical_polys(cs: &[Vec<Z>]) -> R<(Vec<Vec<Z>>, BTreeSet<usize>)> {
    let mut polys: Vec<Vec<Z>> = vec![];
    let mut special: BTreeSet<usize> = BTreeSet::new();
    let full: Vec<usize> = (0..cs.len()).filter(|&i| !cs[i][2].is_zero()).collect();
    let slant: Vec<usize> = (0..cs.len())
        .filter(|&i| cs[i][2].is_zero() && cs[i][0].is_zero() && cs[i][1].is_zero() && !cs[i][4].is_zero())
        .collect();
    for (i, c) in cs.iter().enumerate() {
        if !c[2].is_zero() {
            let d = vec![
                &c[1] * &c[1] - zi(4) * &c[2] * &c[0],
                zi(2) * &c[1] * &c[4] - zi(4) * &c[2] * &c[3],
                &c[4] * &c[4] - zi(4) * &c[2] * &c[5],
            ];
            if poly_trim(&d).is_empty() {
                special.insert(i);
            }
            polys.push(d);
        } else if c[4].is_zero() && c[1].is_zero() && c[0].is_zero() {
            polys.push(vec![c[3].clone(), c[5].clone()]); // vertical line (or a constant)
        }
    }
    for (ii, &i) in full.iter().enumerate() {
        let p = &cs[i];
        let (al, be, ga) = (vec![p[2].clone()], vec![p[1].clone(), p[4].clone()], vec![p[0].clone(), p[3].clone(), p[5].clone()]);
        for &j in &full[ii + 1..] {
            let q = &cs[j];
            let (al2, be2, ga2) = (vec![q[2].clone()], vec![q[1].clone(), q[4].clone()], vec![q[0].clone(), q[3].clone(), q[5].clone()]);
            let u = padd(&pmul(&al2, &be), &pmul(&al, &be2), -1);
            let w = padd(&pmul(&al2, &ga), &pmul(&al, &ga2), -1);
            let v = padd(&pmul(&be, &ga2), &pmul(&be2, &ga), -1);
            let res = padd(&pmul(&w, &w), &pmul(&u, &v), 1);
            if !res.is_empty() {
                polys.push(res);
            } else {
                special.insert(i);
                special.insert(j);
                if !u.is_empty() || !w.is_empty() {
                    polys.extend(cad_shared_line_crossing(p, q, &u, &w)?);
                }
            }
        }
        for &j in &slant {
            let q = &cs[j];
            let (b, cc) = (vec![q[4].clone()], vec![q[3].clone(), q[5].clone()]);
            let t1 = pmul(&al, &pmul(&cc, &cc));
            let t2 = pmul(&be, &pmul(&b, &cc));
            let t3 = pmul(&ga, &pmul(&b, &b));
            let res = padd(&padd(&t1, &t2, -1), &t3, 1);
            if !res.is_empty() {
                polys.push(res);
            } else {
                special.insert(i);
                special.insert(j);
            }
        }
    }
    for (ii, &i) in slant.iter().enumerate() {
        let p = &cs[i];
        for &j in &slant[ii + 1..] {
            let q = &cs[j];
            let res = padd(&pmul(&[p[4].clone()], &[q[3].clone(), q[5].clone()]), &pmul(&[q[4].clone()], &[p[3].clone(), p[5].clone()]), -1);
            if !res.is_empty() {
                polys.push(res);
            } else {
                special.insert(i);
                special.insert(j);
            }
        }
    }
    Ok((polys, special))
}

/// p = L M_p and q = L M_q share the line L (resultant == 0, not proportional): T = alpha' p - alpha q
/// = u(x) y + w(x) = L * N with N vertical, so L = T / gcd(u, w). The abscissa of M_p ^ M_q.
fn cad_shared_line_crossing(p: &[Z], q: &[Z], u: &[Z], w: &[Z]) -> R<Vec<Vec<Z>>> {
    if u.is_empty() {
        bail!("cellHasInterior:internal", "shared line factor with no y term: {:?} {:?}", p, q);
    }
    let g = if w.is_empty() { poly_primitive(u) } else { poly_gcd(u, w)? };
    let lu = cad_fdiv(u, &g)?;
    let lw = cad_fdiv(w, &g)?;
    if lu.len() != 1 || lw.len() > 2 || lu[0].is_zero() {
        bail!("cellHasInterior:internal", "shared factor of two conics is not a line: {:?} {:?}", p, q);
    }
    let ell = lu[0].clone();
    let pp = if lw.len() >= 2 { lw[lw.len() - 2].clone() } else { Q::zero() };
    let rr = if !lw.is_empty() { lw[lw.len() - 1].clone() } else { Q::zero() };
    let mut m = vec![];
    for c in [p, q] {
        let mm = qz(&c[2]) / &ell;
        let kk = (qz(&c[1]) - &pp * &mm) / &ell;
        let hh = (qz(&c[4]) - &rr * &mm) / &ell;
        m.push((mm, kk, hh));
    }
    let det = &m[0].0 * &m[1].1 - &m[1].0 * &m[0].1;
    if det.is_zero() {
        return Ok(vec![]);
    }
    let xp = (&m[1].0 * &m[0].2 - &m[0].0 * &m[1].2) / det;
    Ok(vec![vec![xp.denom().clone(), -xp.numer().clone()]])
}

fn cad_fdiv(p: &[Z], g: &[Z]) -> R<Vec<Q>> {
    let mut p: Vec<Q> = poly_trim(p).iter().map(qz).collect();
    let g: Vec<Q> = poly_trim(g).iter().map(qz).collect();
    if p.is_empty() {
        return Ok(vec![]);
    }
    let mut out = vec![];
    if p.len() >= g.len() {
        for i in 0..=(p.len() - g.len()) {
            let c = &p[i] / &g[0];
            for (j, gj) in g.iter().enumerate() {
                let t = &c * gj;
                p[i + j] -= t;
            }
            out.push(c);
        }
    }
    if p.len() < g.len() || p.iter().any(|x| !x.is_zero()) {
        bail!("cellHasInterior:internal", "inexact polynomial division in a shared line factor");
    }
    Ok(out)
}

/// The y-polynomials q(x0, .) of one column scaled to integers, as (rows in `special`, the others);
/// None when the column holds no point of the cell. Exact shortcuts: a row constant in y with the
/// wrong sign, or a quadratic without real roots with the wrong sign, kills the column; a quadratic
/// without real roots with the right sign holds everywhere, so it adds no sample boundary.
fn cad_y_polys(cs: &[Vec<Z>], sg: &[i64], x0: &Q, special: &BTreeSet<usize>) -> Option<(Vec<Vec<Z>>, Vec<Vec<Z>>)> {
    let (mut sp, mut plain) = (vec![], vec![]);
    for (i, c) in cs.iter().enumerate() {
        let co = [qz(&c[2]), qz(&c[1]) * x0 + qz(&c[4]), qz(&c[0]) * x0 * x0 + qz(&c[3]) * x0 + qz(&c[5])];
        let mut den = zi(1);
        for t in &co {
            den = num_integer::Integer::lcm(&den, t.denom());
        }
        let ints: Vec<Z> = co.iter().map(|t| (t * qz(&den)).to_integer()).collect();
        let ints = poly_primitive(&ints);
        if ints.is_empty() {
            return None;
        }
        if ints.len() == 1 {
            if sg[i] as i32 * sign(&ints[0]) <= 0 {
                return None;
            }
        } else if ints.len() == 3 && (&ints[1] * &ints[1] - zi(4) * &ints[0] * &ints[2]).is_negative() {
            if sg[i] as i32 * sign(&ints[0]) < 0 {
                return None;
            }
        } else if special.contains(&i) {
            sp.push(ints);
        } else {
            plain.push(ints);
        }
    }
    Some((plain, sp))
}

/// Pairwise coprime squarefree primitive polynomials with the same real roots as the product.
fn cad_coprime_base(polys: &[Vec<Z>]) -> R<Vec<Vec<Z>>> {
    let mut work: Vec<Vec<Z>> = polys.iter().filter(|p| poly_trim(p).len() > 1).cloned().collect();
    let mut base: Vec<Vec<Z>> = vec![];
    while let Some(p0) = work.pop() {
        let p = poly_primitive(&poly_squarefree(&p0));
        if p.len() <= 1 {
            continue;
        }
        let mut hit = None;
        for (i, b) in base.iter().enumerate() {
            let g = poly_gcd(&p, b)?;
            if g.len() > 1 {
                hit = Some((i, g));
                break;
            }
        }
        match hit {
            Some((i, g)) => {
                let b = base.remove(i);
                for piece in [poly_div_exact1(&b, &g), poly_div_exact1(&p, &g), Some(g)] {
                    match piece {
                        Some(x) => work.push(x),
                        None => bail!("cellHasInterior:internal", "inexact gcd division"),
                    }
                }
            }
            None => base.push(p),
        }
    }
    Ok(base)
}

struct CadRoot {
    chain: Arc<Vec<Vec<Z>>>,
    lo: Q,
    hi: Q,
}

impl CadRoot {
    fn refine(&mut self) {
        let mid = (&self.lo + &self.hi) / qi(2);
        if sturm_v_q(&self.chain, &mid) - sturm_v_q(&self.chain, &self.hi) == 1 {
            self.lo = mid;
        } else {
            self.hi = mid;
        }
    }
}

/// r < s for roots of coprime polynomials (so r != s): refine until the intervals separate.
fn cad_before(r: &mut CadRoot, s: &mut CadRoot) -> bool {
    loop {
        if r.hi <= s.lo {
            return true;
        }
        if s.hi <= r.lo {
            return false;
        }
        if &r.hi - &r.lo >= &s.hi - &s.lo {
            r.refine()
        } else {
            s.refine()
        }
    }
}

/// Rational sample points: one in every open interval between consecutive distinct real roots of
/// the polynomials `sp` (made coprime first) and `plain` (already squarefree and coprime), one
/// below and one above all roots (0 when there is none). Exact: Sturm isolation of each factor,
/// pairwise separation by bisection.
fn cad_samples(sp: &[Vec<Z>], plain: &[Vec<Z>]) -> R<Vec<Q>> {
    let mut polys = cad_coprime_base(sp)?;
    polys.extend(plain.iter().filter(|p| p.len() > 1).cloned());
    let mut roots: Vec<CadRoot> = vec![];
    for p in &polys {
        let chain = sturm_chain_shared(p);
        for (lo, hi) in isolate_roots(p)? {
            let mut r = CadRoot { chain: chain.clone(), lo, hi };
            // binary insertion keeps `roots` sorted
            let (mut a, mut b) = (0usize, roots.len());
            while a < b {
                let m = (a + b) / 2;
                if cad_before(&mut roots[m], &mut r) {
                    a = m + 1;
                } else {
                    b = m;
                }
            }
            roots.insert(a, r);
        }
    }
    if roots.is_empty() {
        return Ok(vec![Q::zero()]);
    }
    let mut out = vec![roots[0].lo.clone()];
    for i in 0..roots.len() - 1 {
        let (left, right) = roots.split_at_mut(i + 1);
        let (r, t) = (&mut left[i], &mut right[0]);
        while r.hi >= t.lo {
            if &r.hi - &r.lo >= &t.hi - &t.lo {
                r.refine()
            } else {
                t.refine()
            }
        }
        out.push((&r.hi + &t.lo) / qi(2));
    }
    out.push(&roots[roots.len() - 1].hi + qi(1));
    Ok(out)
}

fn probe_interior(ecq: &[Vec<Z>], rows: &[(usize, i64)]) -> R<bool> {
    let idx: Vec<usize> = rows.iter().map(|r| r.0).collect::<BTreeSet<_>>().into_iter().collect();
    let mut pts: Vec<[f64; 2]> = vec![];
    for a in 0..idx.len() {
        for b in a + 1..idx.len() {
            let m = conic_meet(&ecq[idx[a] - 1], &ecq[idx[b] - 1])?;
            pts.extend(m.0.iter().cloned());
        }
    }
    let (ctr, r) = if pts.is_empty() {
        ([0.0, 0.0], 4.0)
    } else {
        let mut sx = 0.0;
        let mut sy = 0.0;
        for p in &pts {
            sx += p[0];
            sy += p[1];
        }
        let ctr = [sx / pts.len() as f64, sy / pts.len() as f64];
        let mut mx = f64::NEG_INFINITY;
        for p in &pts {
            mx = mx.max((p[0] - ctr[0]).abs().max((p[1] - ctr[1]).abs()));
        }
        (ctr, 2f64.max(2.0 * mx))
    };
    let conics: Vec<([f64; 6], f64, i64)> = rows
        .iter()
        .map(|&(e, s)| {
            let (c, cm) = conic_f(&ecq[e - 1]);
            (c, cm, s)
        })
        .collect();
    for i in 0..400 {
        let x0 = ctr[0] + (PROBE_X[i] - 0.5) * 2.0 * r;
        let x1 = ctr[1] + (PROBE_Y[i] - 0.5) * 2.0 * r;
        let m = 1f64.max(x0.abs().max(x1.abs())).powi(2);
        let mut ok = true;
        for (f, cm, s) in &conics {
            let v = f[0] * x0.powi(2) + f[1] * x0 * x1 + f[2] * x1.powi(2) + f[3] * x0 + f[4] * x1 + f[5];
            let scale = 1f64.max(cm * m);
            if *s == 0 {
                if v.abs() > 1e-9 * scale {
                    ok = false;
                    break;
                }
            } else if (*s as f64) * v < -1e-9 * scale {
                ok = false;
                break;
            }
        }
        if ok {
            return Ok(true);
        }
    }
    Ok(false)
}

// ==========================================================================================
// cell_interior: C rows are [conic(6), branch]; rows are (index into C, side)
// ==========================================================================================
pub fn cell_interior(c: &[Vec<Z>], rows: &[(usize, i64)]) -> R<bool> {
    Ok(cell_interior_info(c, rows)?.0)
}

/// (tf, unsure): unsure is MATLAB's info.path == 'unsureTie'.
pub fn cell_interior_info(c: &[Vec<Z>], rows: &[(usize, i64)]) -> R<(bool, bool)> {
    let r = cell_interior_impl(c, rows)?;
    Ok(r)
}

fn cell_interior_impl(c: &[Vec<Z>], rows: &[(usize, i64)]) -> R<(bool, bool)> {
    if c.iter().any(|r| r.len() != 7) {
        bail!("cellInterior:shape", "C must be k x 7 -- [conic(1:6), branch].");
    }
    if rows.is_empty() {
        return Ok((true, false));
    }
    let mut lin: Vec<Vec<Z>> = vec![];
    let mut curved: Vec<Vec<Z>> = vec![];
    let mut sides: Vec<i64> = vec![];
    let mut h: Vec<[f64; 3]> = vec![];
    let mut src: Vec<usize> = vec![];
    for (t, row) in rows.iter().enumerate() {
        let cc = &c[row.0 - 1][0..6];
        let br = &c[row.0 - 1][6];
        let sd = row.1;
        if all_zero(&cc[0..5]) {
            bail!("cellInterior:constant", "row {} is a CONSTANT, not a curve.", t + 1);
        }
        if !br.is_zero() {
            let l = branch_line(cc, br_i(br))?;
            let s = sd as f64;
            h.push([s * l[0], s * l[1], s * l[2]]);
            src.push(t);
        } else if all_zero(&cc[0..3]) {
            let l: Vec<Z> = cc[3..6].iter().map(|x| x * zi(sd)).collect();
            h.push([to_f64(&l[0]), to_f64(&l[1]), to_f64(&l[2])]);
            lin.push(l);
            src.push(t);
        } else {
            curved.push(cc.to_vec());
            sides.push(sd);
        }
    }
    if !lin.is_empty() && !feasible2(&lin, true) {
        return Ok((false, false));
    }
    if !curved.is_empty() {
        let mut tab = curved.clone();
        for l in &lin {
            let mut r = vec![zi(0), zi(0), zi(0)];
            r.extend(l.iter().cloned());
            tab.push(r);
        }
        let mut rr: Vec<(usize, i64)> = sides.iter().enumerate().map(|(i, s)| (i + 1, *s)).collect();
        for i in 0..lin.len() {
            rr.push((curved.len() + i + 1, 1));
        }
        if !cell_has_interior(&tab, &rr)? {
            return Ok((false, false));
        }
    }
    if h.len() < 2 {
        return Ok((true, false));
    }
    let m = margin_lp(&h);
    let hm = h.iter().flatten().fold(1f64, |a, x| a.max(x.abs()));
    let tol = 1e-9 * hm;
    if m > tol {
        return Ok((true, false));
    }
    if m < -tol {
        return Ok((false, false));
    }
    let (empty, skipped) = certificate(c, rows, &src)?;
    if empty {
        return Ok((false, false));
    }
    Ok((true, skipped))
}

fn br_i(b: &Z) -> i64 {
    if b == &zi(1) {
        1
    } else {
        2
    }
}

fn certificate(c: &[Vec<Z>], rows: &[(usize, i64)], src: &[usize]) -> R<(bool, bool)> {
    let mut empty = false;
    let mut skipped = false;
    let m = src.len();
    if m < 2 {
        return Ok((empty, skipped));
    }
    for i in 0..m - 1 {
        for j in i + 1..m {
            let (e, sk) = try_one(c, rows, &[src[i], src[j]])?;
            empty = empty || e;
            skipped = skipped || sk;
            if empty {
                return Ok((empty, skipped));
            }
            for k in j + 1..m {
                let (e, sk) = try_one(c, rows, &[src[i], src[j], src[k]])?;
                empty = empty || e;
                skipped = skipped || sk;
                if empty {
                    return Ok((empty, skipped));
                }
            }
        }
    }
    Ok((empty, skipped))
}

fn unique_rows(m: &[Vec<Z>]) -> Vec<Vec<Z>> {
    m.iter().cloned().collect::<BTreeSet<_>>().into_iter().collect()
}

fn qs(a: Z, b: Z, d: Z) -> Qs {
    [a, b, d]
}
fn q1(a: Z) -> Qs {
    [a, zi(0), zi(1)]
}

fn try_one(c: &[Vec<Z>], rows: &[(usize, i64)], idx: &[usize]) -> R<(bool, bool)> {
    let conics: Vec<Vec<Z>> = idx.iter().filter(|&&t| !c[rows[t].0 - 1][6].is_zero()).map(|&t| c[rows[t].0 - 1][0..6].to_vec()).collect();
    let mut quadp: Vec<Z> = vec![];
    let mut iv = (qi(-1), qi(1));
    let mut rf = 0i64;
    if !conics.is_empty() {
        let nu = unique_rows(&conics).len();
        if nu == 2 {
            return try_two_fields(c, rows, idx);
        }
        if nu > 2 {
            return Ok((false, true));
        }
        let pd = pair_line_data_m(&conics[0])?;
        quadp = pd.quad.clone();
        for &t in idx {
            if !c[rows[t].0 - 1][6].is_zero() {
                rf = br_i(&c[rows[t].0 - 1][6]);
                break;
            }
        }
        if quadp.len() >= 3 {
            iv = pd.iv[(rf - 1) as usize].clone();
        }
    }
    let mut a1 = vec![];
    let mut a2 = vec![];
    let mut cc = vec![];
    for &t in idx {
        match exact_half_plane(c, rows, t, rf)? {
            None => return Ok((false, true)),
            Some((a, k)) => {
                a1.push(a[0].clone());
                a2.push(a[1].clone());
                cc.push(k);
            }
        }
    }
    if quadp.is_empty() || quadp.len() == 2 {
        quadp = vec![zi(1), zi(0)];
    }
    let lam: Vec<Qs>;
    if a1.len() == 2 {
        let cr = cross2(&a1[0], &a2[0], &a1[1], &a2[1], &quadp)?;
        if qs_sign(&cr, &quadp, &iv)? != 0 {
            return Ok((false, false));
        }
        let dt = qs_add(&qs_mul(&a1[0], &a1[1], &quadp)?, &qs_mul(&a2[0], &a2[1], &quadp)?)?;
        if qs_sign(&dt, &quadp, &iv)? >= 0 {
            return Ok((false, false));
        }
        let n2 = qs_add(&qs_mul(&a1[1], &a1[1], &quadp)?, &qs_mul(&a2[1], &a2[1], &quadp)?)?;
        lam = vec![n2, qneg(&dt)];
    } else {
        let mut l = vec![
            cross2(&a1[1], &a2[1], &a1[2], &a2[2], &quadp)?,
            cross2(&a1[2], &a2[2], &a1[0], &a2[0], &quadp)?,
            cross2(&a1[0], &a2[0], &a1[1], &a2[1], &quadp)?,
        ];
        let mut sg = vec![];
        for x in &l {
            sg.push(qs_sign(x, &quadp, &iv)?);
        }
        if sg.iter().all(|s| *s <= 0) {
            l = l.iter().map(qneg).collect();
            sg = sg.iter().map(|s| -s).collect();
        }
        if !sg.iter().all(|s| *s >= 0) {
            return Ok((false, false));
        }
        if sg.iter().all(|s| *s == 0) {
            return Ok((false, false));
        }
        lam = l;
    }
    let mut tot = qs(zi(0), zi(0), zi(1));
    for k in 0..cc.len() {
        tot = qs_add(&tot, &qs_mul(&lam[k], &cc[k], &quadp)?)?;
    }
    Ok((qs_sign(&tot, &quadp, &iv)? <= 0, false))
}

fn exact_half_plane(c: &[Vec<Z>], rows: &[(usize, i64)], t: usize, rf: i64) -> R<Option<([Qs; 2], Qs)>> {
    let con = &c[rows[t].0 - 1][0..6];
    let br = &c[rows[t].0 - 1][6];
    let sd = zi(rows[t].1);
    if br.is_zero() {
        if !all_zero(&con[0..3]) {
            return Ok(None);
        }
        return Ok(Some(([q1(&sd * &con[3]), q1(&sd * &con[4])], q1(&sd * &con[5]))));
    }
    let b = br_i(br);
    let pd = pair_line_data_m(con)?;
    if !pd.parallel && pd.vertical && b == 2 {
        let p = &pd.apex;
        return Ok(Some(([q1(&sd * &p[2]), q1(zi(0))], q1(-&sd * &p[0]))));
    }
    let tt = match root_element(&pd.quad, b, rf)? {
        None => return Ok(None),
        Some(t) => t,
    };
    if pd.parallel {
        let a = [q1(&sd * &pd.form[0]), q1(&sd * &pd.form[1])];
        let k = qs_mul(&q1(-&sd), &tt, &pd.quad)?;
        Ok(Some((a, k)))
    } else {
        let p = &pd.apex;
        let a = [qs_mul(&q1(-&sd * &p[2]), &tt, &pd.quad)?, q1(&sd * &p[2])];
        let k = qs_add(&q1(-&sd * &p[1]), &qs_mul(&q1(&sd * &p[0]), &tt, &pd.quad)?)?;
        Ok(Some((a, k)))
    }
}

fn try_two_fields(c: &[Vec<Z>], rows: &[(usize, i64)], idx: &[usize]) -> R<(bool, bool)> {
    let cs: Vec<Vec<Z>> = idx.iter().filter(|&&t| !c[rows[t].0 - 1][6].is_zero()).map(|&t| c[rows[t].0 - 1][0..6].to_vec()).collect();
    let u = unique_rows(&cs);
    let pd1 = pair_line_data_m(&u[0])?;
    let pd2 = pair_line_data_m(&u[1])?;
    let qq1 = pd1.quad.clone();
    let qq2 = pd2.quad.clone();
    if qq1.len() != 3 || qq2.len() != 3 {
        return Ok((false, true));
    }
    let mut ref1 = 0;
    let mut ref2 = 0;
    for &t in idx {
        let b = &c[rows[t].0 - 1][6];
        if b.is_zero() {
            continue;
        }
        let cc = &c[rows[t].0 - 1][0..6];
        if cc == &u[0][..] && ref1 == 0 {
            ref1 = br_i(b);
        }
        if cc == &u[1][..] && ref2 == 0 {
            ref2 = br_i(b);
        }
    }
    let iv1 = pd1.iv[(ref1 - 1) as usize].clone();
    let iv2 = pd2.iv[(ref2 - 1) as usize].clone();
    let mut a1 = vec![];
    let mut a2 = vec![];
    let mut cc = vec![];
    for &t in idx {
        let b = &c[rows[t].0 - 1][6];
        let (res, slot) = if b.is_zero() {
            (exact_half_plane(c, rows, t, 0)?, 0)
        } else if &c[rows[t].0 - 1][0..6] == &u[0][..] {
            (exact_half_plane(c, rows, t, ref1)?, 1)
        } else {
            (exact_half_plane(c, rows, t, ref2)?, 2)
        };
        let (a, k) = match res {
            None => return Ok((false, true)),
            Some(x) => x,
        };
        a1.push(embed(&a[0], slot)?);
        a2.push(embed(&a[1], slot)?);
        cc.push(embed(&k, slot)?);
    }
    let n = a1.len();
    let lam: Vec<Qs2>;
    if n == 2 {
        let cr = cross2two(&a1[0], &a2[0], &a1[1], &a2[1], &qq1, &qq2)?;
        let (sg, ok) = qs2_sign(&cr, &qq1, &iv1, &qq2, &iv2)?;
        if !ok {
            return Ok((false, true));
        }
        if sg != 0 {
            return Ok((false, false));
        }
        let dt = qs2_add(&qs2_mul(&a1[0], &a1[1], &qq1, &qq2)?, &qs2_mul(&a2[0], &a2[1], &qq1, &qq2)?)?;
        let (sg, ok) = qs2_sign(&dt, &qq1, &iv1, &qq2, &iv2)?;
        if !ok {
            return Ok((false, true));
        }
        if sg >= 0 {
            return Ok((false, false));
        }
        let n2 = qs2_add(&qs2_mul(&a1[1], &a1[1], &qq1, &qq2)?, &qs2_mul(&a2[1], &a2[1], &qq1, &qq2)?)?;
        lam = vec![n2, qneg2(&dt)];
    } else {
        let mut l = vec![
            cross2two(&a1[1], &a2[1], &a1[2], &a2[2], &qq1, &qq2)?,
            cross2two(&a1[2], &a2[2], &a1[0], &a2[0], &qq1, &qq2)?,
            cross2two(&a1[0], &a2[0], &a1[1], &a2[1], &qq1, &qq2)?,
        ];
        let mut sg = vec![0; 3];
        for k in 0..3 {
            let (s, ok) = qs2_sign(&l[k], &qq1, &iv1, &qq2, &iv2)?;
            if !ok {
                return Ok((false, true));
            }
            sg[k] = s;
        }
        if sg.iter().all(|s| *s <= 0) {
            l = l.iter().map(qneg2).collect();
            sg = sg.iter().map(|s| -s).collect();
        }
        if !sg.iter().all(|s| *s >= 0) {
            return Ok((false, false));
        }
        if sg.iter().all(|s| *s == 0) {
            return Ok((false, false));
        }
        lam = l;
    }
    let mut tot: Qs2 = [zi(0), zi(0), zi(0), zi(0), zi(1)];
    for k in 0..cc.len() {
        tot = qs2_add(&tot, &qs2_mul(&lam[k], &cc[k], &qq1, &qq2)?)?;
    }
    let (sg, ok) = qs2_sign(&tot, &qq1, &iv1, &qq2, &iv2)?;
    if !ok {
        return Ok((false, true));
    }
    Ok((sg <= 0, false))
}

fn embed(u: &Qs, slot: i32) -> R<Qs2> {
    let w = match slot {
        1 => [u[0].clone(), u[1].clone(), zi(0), zi(0), u[2].clone()],
        2 => [u[0].clone(), zi(0), u[1].clone(), zi(0), u[2].clone()],
        _ => {
            if !u[1].is_zero() {
                bail!("cellInterior:embed", "a rational row cannot carry a t term.");
            }
            [u[0].clone(), zi(0), zi(0), zi(0), u[2].clone()]
        }
    };
    qs2_red(&w)
}

fn cross2two(u1: &Qs2, u2: &Qs2, v1: &Qs2, v2: &Qs2, q1: &[Z], q2: &[Z]) -> R<Qs2> {
    qs2_add(&qs2_mul(u1, v2, q1, q2)?, &qneg2(&qs2_mul(u2, v1, q1, q2)?))
}

fn qneg2(u: &Qs2) -> Qs2 {
    [-&u[0], -&u[1], -&u[2], -&u[3], u[4].clone()]
}

fn root_element(quad: &[Z], br: i64, rf: i64) -> R<Option<Qs>> {
    if quad.len() == 2 {
        return Ok(Some([-&quad[1], zi(0), quad[0].clone()]));
    }
    if quad.len() != 3 {
        return Ok(None);
    }
    let t = if br == rf { [zi(0), zi(1), zi(1)] } else { [-&quad[1], -&quad[0], quad[0].clone()] };
    Ok(Some(qs_red(&t)?))
}

fn cross2(u1: &Qs, u2: &Qs, v1: &Qs, v2: &Qs, q: &[Z]) -> R<Qs> {
    qs_add(&qs_mul(u1, v2, q)?, &qneg(&qs_mul(u2, v1, q)?))
}

fn qneg(u: &Qs) -> Qs {
    [-&u[0], -&u[1], u[2].clone()]
}

fn branch_line(c: &[Z], b: i64) -> R<[f64; 3]> {
    let pe = pair_eval_prep(c, b)?;
    let v0 = pe.eval_m(0.0, 0.0).0;
    let v1 = pe.eval_m(1.0, 0.0).0;
    let v2 = pe.eval_m(0.0, 1.0).0;
    let mut l = [v1 - v0, v2 - v0, v0];
    let n = l[0].hypot(l[1]);
    if n > 0.0 {
        l = [l[0] / n, l[1] / n, l[2] / n];
    }
    Ok(l)
}

/// Largest t with H*[x;1] >= t inside a big box, by vertex enumeration (doubles).
fn margin_lp(h: &[[f64; 3]]) -> f64 {
    let r = 1e4;
    let mut m: Vec<[f64; 3]> = vec![];
    let mut b: Vec<f64> = vec![];
    for row in h {
        let n = row[0].hypot(row[1]);
        let hn = if n > 0.0 { [row[0] / n, row[1] / n, row[2] / n] } else { *row };
        m.push([hn[0], hn[1], -1.0]);
        b.push(-hn[2]);
    }
    for (row, bb) in [([1.0, 0.0, 0.0], -r), ([-1.0, 0.0, 0.0], -r), ([0.0, 1.0, 0.0], -r), ([0.0, -1.0, 0.0], -r)] {
        m.push(row);
        b.push(bb);
    }
    let k = m.len();
    let mut best = f64::NEG_INFINITY;
    let mut any = false;
    let mut anyfeas = false;
    for i in 0..k - 2 {
        for j in i + 1..k - 1 {
            for l in j + 1..k {
                let s = [m[i], m[j], m[l]];
                if det3(&s).abs() < 1e-12 {
                    continue;
                }
                any = true;
                let z = match solve3(&s, [b[i], b[j], b[l]]) {
                    None => continue,
                    Some(z) => z,
                };
                let feas = (0..k).all(|t| m[t][0] * z[0] + m[t][1] * z[1] + m[t][2] * z[2] >= b[t] - 1e-9);
                if feas {
                    anyfeas = true;
                    best = best.max(z[2]);
                }
            }
        }
    }
    if !any || !anyfeas {
        return -1.0;
    }
    best
}

// ==========================================================================================
// cell_sample_point
// ==========================================================================================
/// A strictly interior point of a cell with curved rows, by sampling the polygon of its straight
/// rows. Returns (pt, ok, P).
pub fn cell_sample_point(c: &[Vec<Z>], rows: &[(usize, i64)], h: &[[f64; 3]], r: f64, dense: bool) -> R<([f64; 2], bool, Vec<[f64; 2]>)> {
    let mut pt = [0.0, 0.0];
    let mut ok = false;
    let mut p: Vec<[f64; 2]> = vec![];
    let mut v = polygon_vertices(h, r);
    if v.len() < 3 {
        return Ok((pt, ok, p));
    }
    let mut sx = 0.0;
    let mut sy = 0.0;
    for q in &v {
        sx += q[0];
        sy += q[1];
    }
    let c0 = [sx / v.len() as f64, sy / v.len() as f64];
    v.sort_by(|a, b| {
        let ta = (a[1] - c0[1]).atan2(a[0] - c0[0]);
        let tb = (b[1] - c0[1]).atan2(b[0] - c0[0]);
        ta.total_cmp(&tb)                       // total order: NaN-safe, never panics (review #6)
    });
    let prep = SlackPrep::new(c, rows)?;
    for n in [12usize, 48, 192] {
        let (bi, bj) = barycentric(n);
        let mut q: Vec<[f64; 2]> = Vec::with_capacity(bi.len() * v.len());
        for e in 0..v.len() {
            let a = v[e];
            let b = v[(e + 1) % v.len()];
            for t in 0..bi.len() {
                q.push([c0[0] + bi[t] * (a[0] - c0[0]) + bj[t] * (b[0] - c0[0]), c0[1] + bi[t] * (a[1] - c0[1]) + bj[t] * (b[1] - c0[1])]);
            }
        }
        let s: Vec<f64> = q.iter().map(|x| prep.slack(x[0], x[1])).collect();
        let any_pos = s.iter().any(|x| *x > 0.0);
        let argmax = || {
            let mut bi = 0;
            for i in 1..s.len() {
                if s[i] > s[bi] {
                    bi = i;
                }
            }
            bi
        };
        if dense && n < 192 && (n < 48 || !any_pos) {
            if any_pos && !ok {
                pt = q[argmax()];
                ok = true;
            }
            p.extend(q.iter().zip(s.iter()).filter(|(_, x)| **x > 0.0).map(|(a, _)| *a));
            continue;
        }
        if any_pos || ok {
            if !ok {
                pt = q[argmax()];
                ok = true;
            }
            p.extend(q.iter().zip(s.iter()).filter(|(_, x)| **x > 0.0).map(|(a, _)| *a));
            let cap = if dense { 4000 } else { 400 };
            if p.len() > cap {
                let idx: Vec<usize> = linspace(1.0, p.len() as f64, cap).iter().map(|t| round_half_away(*t) as usize).collect();
                p = idx.iter().map(|i| p[i - 1]).collect();
            }
            return Ok((pt, ok, p));
        }
    }
    Ok((pt, ok, p))
}

fn barycentric(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut bi = vec![];
    let mut bj = vec![];
    for i in 0..=n {
        for j in 0..=n {
            if i + j < n && i > 0 && j > 0 {
                bi.push(i as f64 / n as f64);
                bj.push(j as f64 / n as f64);
            }
        }
    }
    (bi, bj)
}

/// The per-row data of _slackArr, prepared once.
pub struct SlackPrep {
    rows: Vec<SlackRow>,
}
enum SlackRow {
    Branch(PairEval, f64),
    Conic([f64; 6], [f64; 6], f64, f64),
}

impl SlackPrep {
    pub fn new(c: &[Vec<Z>], rows: &[(usize, i64)]) -> R<SlackPrep> {
        let mut out = vec![];
        for &(i, sd) in rows {
            let cc = &c[i - 1][0..6];
            let br = &c[i - 1][6];
            if !br.is_zero() {
                out.push(SlackRow::Branch(pair_eval_prep(cc, br_i(br))?, sd as f64));
            } else {
                let f: Vec<f64> = cc.iter().map(to_f64).collect();
                let a: Vec<f64> = cc.iter().map(|x| to_f64(&x.abs())).collect();
                let cm = to_f64(&cc.iter().map(|x| x.abs()).max().unwrap());
                out.push(SlackRow::Conic([f[0], f[1], f[2], f[3], f[4], f[5]], [a[0], a[1], a[2], a[3], a[4], a[5]], cm, sd as f64));
            }
        }
        Ok(SlackPrep { rows: out })
    }

    /// min over the rows of side * value / scale at (x, y), with the 32-eps rounding band.
    pub fn slack(&self, x: f64, y: f64) -> f64 {
        let xx = x * x;
        let xy = x * y;
        let yy = y * y;
        let nx = 1f64.max((xx + yy).sqrt());
        let nx2 = nx * nx;
        let eps32 = 32.0 * 2f64.powi(-52);
        let mut s = f64::INFINITY;
        for r in &self.rows {
            let (v, m, sc, sd) = match r {
                SlackRow::Branch(pe, sd) => {
                    let (v, m) = pe.eval_m(x, y);
                    (v, m, nx, *sd)
                }
                SlackRow::Conic(c, a, cm, sd) => {
                    let v = c[0] * xx + c[1] * xy + c[2] * yy + c[3] * x + c[4] * y + c[5];
                    let m = a[0] * xx + (c[1] * xy).abs() + a[2] * yy + (c[3] * x).abs() + (c[4] * y).abs() + a[5];
                    (v, m, cm * nx2, *sd)
                }
            };
            let mut w = sd * v;
            if w > 0.0 && w <= eps32 * m {
                w = 0.0;
            }
            s = s.min(w / sc);
        }
        s
    }
}

fn polygon_vertices(h: &[[f64; 3]], r: f64) -> Vec<[f64; 2]> {
    let mut a: Vec<[f64; 3]> = h.to_vec();
    a.extend([[1.0, 0.0, r], [-1.0, 0.0, r], [0.0, 1.0, r], [0.0, -1.0, r]]);
    let mut v = vec![];
    for i in 0..a.len() {
        for j in i + 1..a.len() {
            let m = [[a[i][0], a[i][1]], [a[j][0], a[j][1]]];
            if (m[0][0] * m[1][1] - m[0][1] * m[1][0]).abs() < 1e-12 {
                continue;
            }
            let x = solve2f(&m, [-a[i][2], -a[j][2]]);
            let lim = -1e-9 * 1f64.max((x[0] * x[0] + x[1] * x[1]).sqrt());
            if a.iter().all(|aa| aa[0] * x[0] + aa[1] * x[1] + aa[2] >= lim) {
                v.push(x);
            }
        }
    }
    if v.is_empty() {
        return v;
    }
    v.sort_by(|p, q| p[0].total_cmp(&q[0]).then(p[1].total_cmp(&q[1])));   // NaN-safe (review #6)
    let mut out: Vec<[f64; 2]> = vec![];
    for p in v {
        if out.iter().any(|u| (p[0] - u[0]).abs() <= 1e-12 && (p[1] - u[1]).abs() <= 1e-12) {
            continue;
        }
        out.push(p);
    }
    out
}

// ==========================================================================================
// compact_cells / merge_adjacent_cells
// ==========================================================================================
pub type FcRow = [i64; 3];

pub fn compact_cells(fn_: Vec<Vec<Z>>, fd: Vec<Z>, fc: Vec<Vec<FcRow>>, ecq: &[Vec<Z>]) -> R<(Vec<Vec<Z>>, Vec<Z>, Vec<Vec<FcRow>>)> {
    if fc.is_empty() {
        return Ok((fn_, fd, fc));
    }
    let (mut fn_, mut fd, mut fc) = (fn_, fd, fc);
    for _ in 0..8 {
        let mut drops = 0;
        // every cell sheds its implied rows independently
        let res = crate::par::map_range(fc.len(), |k| drop_implied(fc[k].clone(), ecq))?;
        for (k, (r, d)) in res.into_iter().enumerate() {
            fc[k] = r;
            drops += d;
        }
        let n0 = fc.len();
        let m = merge_adjacent_cells(fn_, fd, fc);
        fn_ = m.0;
        fd = m.1;
        fc = m.2;
        if drops == 0 && fc.len() == n0 {
            break;
        }
    }
    Ok((fn_, fd, fc))
}

fn drop_implied(rows: Vec<FcRow>, ecq: &[Vec<Z>]) -> R<(Vec<FcRow>, usize)> {
    let mut dropped = 0;
    if rows.len() < 2 {
        return Ok((rows, dropped));
    }
    let (c0, l0) = local_table(&rows, ecq);
    match cell_interior(&c0, &l0) {
        Ok(true) => {}
        Ok(false) => return Ok((rows, dropped)),
        Err(_) => return Ok((rows, dropped)),
    }
    let mut rows = rows;
    let mut t = 0;
    while t < rows.len() {
        if rows.len() < 2 || rows[t][1] == 0 {
            t += 1;
            continue;
        }
        let mut flipped = rows.clone();
        flipped[t][1] = -flipped[t][1];
        let (c, l) = local_table(&flipped, ecq);
        match cell_interior(&c, &l) {
            Ok(false) => {
                rows.remove(t);
                dropped += 1;
                continue;
            }
            Ok(true) => {}
            Err(me) => {
                if me.id != "ratQ:overflow" && me.id != "cellInterior:constant" {
                    return Err(me);
                }
            }
        }
        t += 1;
    }
    Ok((rows, dropped))
}

fn local_table(rows: &[FcRow], ecq: &[Vec<Z>]) -> (Vec<Vec<Z>>, Vec<(usize, i64)>) {
    let mut c: Vec<Vec<Z>> = vec![];
    let mut local = vec![];
    for r in rows {
        let mut key = ecq[r[0] as usize - 1].clone();
        key.push(zi(r[2]));
        let hit = match c.iter().position(|x| *x == key) {
            Some(h) => h,
            None => {
                c.push(key);
                c.len() - 1
            }
        };
        local.push((hit + 1, r[1]));
    }
    (c, local)
}

fn normalize_rows(r: &[FcRow]) -> Vec<FcRow> {
    r.iter().cloned().collect::<BTreeSet<_>>().into_iter().collect()
}

fn mergeable_on(a: &[FcRow], b: &[FcRow]) -> Option<i64> {
    if a.len() != b.len() {
        return None;
    }
    let d: Vec<usize> = (0..a.len()).filter(|&i| a[i] != b[i]).collect();
    if d.len() != 1 {
        return None;
    }
    let d = d[0];
    if a[d][0] != b[d][0] || a[d][1] == b[d][1] {
        return None;
    }
    Some(a[d][0])
}

pub fn merge_adjacent_cells(fn_: Vec<Vec<Z>>, fd: Vec<Z>, fc: Vec<Vec<FcRow>>) -> (Vec<Vec<Z>>, Vec<Z>, Vec<Vec<FcRow>>) {
    let mut fn_ = fn_;
    let mut fd = fd;
    let mut fc: Vec<Vec<FcRow>> = fc.iter().map(|r| normalize_rows(r)).collect();
    let mut changed = true;
    while changed {
        changed = false;
        let mut groups: BTreeMap<(Vec<Z>, Z), Vec<usize>> = BTreeMap::new();
        for k in 0..fc.len() {
            groups.entry((fn_[k].clone(), fd[k].clone())).or_default().push(k);
        }
        let mut dead: BTreeSet<usize> = BTreeSet::new();
        for (_, idx) in groups.iter() {
            if idx.len() < 2 {
                continue;
            }
            let mut done = false;
            for a in 0..idx.len() {
                for b in a + 1..idx.len() {
                    let (ka, kb) = (idx[a], idx[b]);
                    let c = match mergeable_on(&fc[ka], &fc[kb]) {
                        None => continue,
                        Some(c) => c,
                    };
                    let kept: Vec<FcRow> = fc[ka].iter().filter(|r| r[0] != c).cloned().collect();
                    fc[ka] = normalize_rows(&kept);
                    fc[kb] = vec![];
                    dead.insert(kb);
                    changed = true;
                    done = true;
                    break;
                }
                if done {
                    break;
                }
            }
        }
        if !dead.is_empty() {
            let alive: Vec<usize> = (0..fc.len()).filter(|k| !dead.contains(k)).collect();
            fn_ = alive.iter().map(|&k| fn_[k].clone()).collect();
            fd = alive.iter().map(|&k| fd[k].clone()).collect();
            fc = alive.iter().map(|&k| fc[k].clone()).collect();
        }
    }
    (fn_, fd, fc)
}

// ==========================================================================================
// assemble_qua_con_cells
// ==========================================================================================
/// A cell: face function num/den and H-form rows [conic(6), side, branch].
#[derive(Clone, Debug)]
pub struct Cell {
    pub num: Vec<Z>,
    pub den: Z,
    pub con: Vec<Vec<Z>>,
}

/// conic([0 0 0 row]) + [side, branch] as an 8-wide con row.
pub fn con_row(conic6: Vec<Z>, side: i64, branch: i64) -> Vec<Z> {
    let mut r = conic6;
    r.push(zi(side));
    r.push(zi(branch));
    r
}

fn both_sides(rows: &[Vec<Z>]) -> bool {
    let mut sides: HashMap<Vec<Z>, [bool; 2]> = HashMap::new();
    for r in rows {
        let mut key = r[0..6].to_vec();
        key.push(r[7].clone());
        let s = sides.entry(key).or_insert([false, false]);
        if r[6].is_positive() {
            s[0] = true;
        } else if r[6].is_negative() {
            s[1] = true;
        }
    }
    sides.values().any(|s| s[0] && s[1])
}

fn conic_at(c: &[Z], xn: &[Z], xd: &Z) -> Z {
    &c[0] * &xn[0] * &xn[0] + &c[1] * &xn[0] * &xn[1] + &c[2] * &xn[1] * &xn[1] + &c[3] * &xn[0] * xd + &c[4] * &xn[1] * xd + &c[5] * xd * xd
}

fn cell_holds_at(ecq: &[Vec<Z>], rows: &[FcRow], xn: &[Z], xd: &Z) -> bool {
    for r in rows {
        let v = conic_at(&ecq[r[0] as usize - 1], xn, xd);
        if r[1] == 0 {
            if !v.is_zero() {
                return false;
            }
        } else if sign(&v) * (r[1] as i32) < 0 {
            return false;
        }
    }
    true
}

fn cell_holds_at_point(rows: &[Vec<Z>], xn: &[Z], xd: &Z) -> bool {
    for r in rows {
        let v = conic_at(&r[0..6], xn, xd);
        if r[6].is_zero() {
            if !v.is_zero() {
                return false;
            }
        } else if sign(&v) * sign(&r[6]) < 0 {
            return false;
        }
    }
    true
}

fn single_point_domain(cells: &[Cell]) -> R<Option<(Vec<Z>, Z)>> {
    let mut seen: Option<(Z, Z, Z)> = None;
    for c in cells {
        let lin: Vec<&Vec<Z>> = c.con.iter().filter(|r| all_zero(&r[0..3])).collect();
        for a in 0..lin.len() {
            for b in a + 1..lin.len() {
                let am = vec![vec![lin[a][3].clone(), lin[a][4].clone()], vec![lin[b][3].clone(), lin[b][4].clone()]];
                if det_exact(&am)?.is_zero() {
                    continue;
                }
                let (x, xd) = solve2(&am, &[-&lin[a][5], -&lin[b][5]])?;
                if !cell_holds_at_point(&c.con, &x, &xd) {
                    continue;
                }
                match &seen {
                    None => seen = Some((x[0].clone(), x[1].clone(), xd)),
                    Some(s) => {
                        if &s.0 * &xd != &x[0] * &s.2 || &s.1 * &xd != &x[1] * &s.2 {
                            return Ok(None);
                        }
                    }
                }
            }
        }
    }
    Ok(seen.map(|s| (vec![s.0, s.1], s.2)))
}

pub fn assemble_qua_con_cells(cells: Vec<Cell>) -> R<QuaCon> {
    if cells.is_empty() {
        return Ok(QuaCon::empty());
    }
    let mut cells = cells;
    for (k, c) in cells.iter_mut().enumerate() {
        for r in c.con.iter_mut() {
            if r.len() == 7 {
                r.push(zi(0));
            }
            if r.len() != 8 {
                bail!("assembleQuaConCells:conShape", "cell {} has a `con` that is not k x 7 or k x 8.", k + 1);
            }
        }
    }
    // every cell is tested for liveness on its own (a cell's rows are simplified in the process)
    let verdicts = crate::par::map_range(cells.len(), |k| -> R<(bool, Option<Vec<Vec<Z>>>)> {
        let rows = cells[k].con.clone();
        let has_branch = rows.iter().any(|r| !r[7].is_zero());
        if rows.iter().any(|r| r[6].is_zero()) {
            return Ok((!both_sides(&rows), None));
        }
        let mut keep = vec![true; rows.len()];
        for (i, r) in rows.iter().enumerate() {
            if !r[7].is_zero() {
                continue;
            }
            let cs = conic_sign(&r[0..6]);
            if cs == 0 {
                continue;
            }
            if zi(cs as i64) == r[6] {
                keep[i] = false;
            } else {
                return Ok((false, None));
            }
        }
        let rows: Vec<Vec<Z>> = rows.into_iter().zip(keep.iter()).filter(|(_, kp)| **kp).map(|(r, _)| r).collect();
        if rows.is_empty() {
            return Ok((true, Some(rows)));
        }
        if both_sides(&rows) {
            return Ok((false, Some(rows)));
        }
        let is_line: Vec<bool> = rows.iter().map(|r| all_zero(&r[0..3]) && r[7].is_zero()).collect();
        let lin: Vec<Vec<Z>> = rows.iter().zip(is_line.iter()).filter(|(_, l)| **l).map(|(r, _)| r[3..6].iter().map(|x| &r[6] * x).collect()).collect();
        let mut lv = feasible2(&lin, true);
        if lv && has_branch {
            return Ok((lv, Some(rows)));
        }
        if lv && !is_line.iter().all(|x| *x) {
            let tab: Vec<Vec<Z>> = rows.iter().map(|r| r[0..6].to_vec()).collect();
            let rr: Vec<(usize, i64)> = rows.iter().enumerate().map(|(i, r)| (i + 1, side_i(&r[6]))).collect();
            lv = cell_has_interior(&tab, &rr)?;
        }
        Ok((lv, Some(rows)))
    })?;
    let mut live = Vec::with_capacity(cells.len());
    for (k, (lv, rows)) in verdicts.into_iter().enumerate() {
        if let Some(r) = rows {
            cells[k].con = r;
        }
        live.push(lv);
    }

    if !live.iter().any(|x| *x) {
        if let Some((pt, den)) = single_point_domain(&cells)? {
            let mut best: Option<(Vec<Z>, Z)> = None;
            let xf = [q_to_f64(&Q::new(pt[0].clone(), den.clone())), q_to_f64(&Q::new(pt[1].clone(), den.clone()))];
            for c in &cells {
                if !cell_holds_at_point(&c.con, &pt, &den) {
                    continue;
                }
                let better = match &best {
                    None => true,
                    Some((bn, bd)) => eval_face(&c.num, &c.den, xf[0], xf[1]) > eval_face(bn, bd, xf[0], xf[1]),
                };
                if better {
                    best = Some((c.num.clone(), c.den.clone()));
                }
            }
            if let Some((bn, bd)) = best {
                let r1 = vec![den.clone(), zi(0), -&pt[0]];
                let r2 = vec![zi(0), den.clone(), -&pt[1]];
                let con = vec![con_row(conic_lin(&r1)?, 0, 0), con_row(conic_lin(&r2)?, 0, 0)];
                cells = vec![Cell { num: bn, den: bd, con }];
                live = vec![true];
            }
        }
    }
    let cells: Vec<Cell> = cells.into_iter().zip(live.iter()).filter(|(_, l)| **l).map(|(c, _)| c).collect();
    if cells.is_empty() {
        return Ok(QuaCon::empty());
    }

    // ---- tabulate the conics
    let mut ecq: Vec<Vec<Z>> = vec![];
    let mut index: HashMap<Vec<Z>, usize> = HashMap::new();
    let mut fc: Vec<Vec<FcRow>> = vec![];
    for c in &cells {
        let mut f = vec![];
        for r in &c.con {
            let key = r[0..6].to_vec();
            let hit = match index.get(&key) {
                Some(h) => *h,
                None => {
                    ecq.push(key.clone());
                    index.insert(key, ecq.len());
                    ecq.len()
                }
            };
            f.push([hit as i64, side_i(&r[6]), side_i(&r[7])]);
        }
        fc.push(f);
    }
    let mut fn_ = vec![];
    let mut fd = vec![];
    for c in &cells {
        let (n, d) = canon(&c.num, &c.den)?;
        fn_.push(n);
        fd.push(d);
    }
    let (fn_, fd, mut fc) = compact_cells(fn_, fd, fc, &ecq)?;

    let mut used = vec![false; ecq.len()];
    for f in &fc {
        for r in f {
            used[r[0] as usize - 1] = true;
        }
    }
    let mut renum = vec![];
    let mut acc = 0;
    for u in &used {
        if *u {
            acc += 1;
        }
        renum.push(acc);
    }
    let ecq: Vec<Vec<Z>> = ecq.into_iter().zip(used.iter()).filter(|(_, u)| **u).map(|(c, _)| c).collect();
    for f in fc.iter_mut() {
        for r in f.iter_mut() {
            r[0] = renum[r[0] as usize - 1];
        }
    }

    // ---- name the corners
    let mut seen: BTreeSet<(i64, i64)> = BTreeSet::new();
    let mut vname: Vec<[usize; 3]> = vec![];
    for rows in &fc {
        for p1 in 0..rows.len() {
            for p2 in p1 + 1..rows.len() {
                let i = rows[p1][0];
                let j = rows[p2][0];
                if i == j {
                    continue;
                }
                let (lo, hi) = (i.min(j), i.max(j));
                if seen.contains(&(lo, hi)) {
                    continue;
                }
                let ci = &ecq[i as usize - 1];
                let cj = &ecq[j as usize - 1];
                if all_zero(&ci[0..3]) && all_zero(&cj[0..3]) {
                    let a = vec![vec![ci[3].clone(), ci[4].clone()], vec![cj[3].clone(), cj[4].clone()]];
                    if det_exact(&a)?.is_zero() {
                        continue;
                    }
                    let (xn, xd) = solve2(&a, &[-&ci[5], -&cj[5]])?;
                    if !cell_holds_at(&ecq, rows, &xn, &xd) {
                        continue;
                    }
                    seen.insert((lo, hi));
                    vname.push([lo as usize, hi as usize, 1]);
                } else {
                    let m = conic_meet(ci, cj)?;
                    if m.1.degenerate {
                        continue;
                    }
                    let rr: Vec<(usize, i64)> = rows.iter().map(|r| (r[0] as usize, r[1])).collect();
                    for (rt, p) in m.0.iter().enumerate() {
                        if !satisfies_all_double(&ecq, &rr, *p) {
                            continue;
                        }
                        seen.insert((lo, hi));
                        vname.push([lo as usize, hi as usize, rt + 1]);
                        break;
                    }
                }
            }
        }
    }
    let (e, f) = build_incidence(&ecq, &fc, &vname);
    QuaCon::new(vname, e, ecq, fn_, fd, f, fc)
}

pub fn side_i(x: &Z) -> i64 {
    num_traits::ToPrimitive::to_i64(x).expect("a side / branch fits in i64")
}

fn build_incidence(ecq: &[Vec<Z>], fc: &[Vec<FcRow>], vname: &[[usize; 3]]) -> (Vec<[usize; 3]>, Vec<[usize; 2]>) {
    let ne = ecq.len();
    let mut e = vec![[0usize; 3]; ne];
    let mut f = vec![[0usize; 2]; ne];
    for j in 1..=ne {
        let vs: Vec<usize> = vname.iter().enumerate().filter(|(_, v)| v[0] == j || v[1] == j).map(|(i, _)| i + 1).collect();
        if vs.len() >= 2 {
            e[j - 1] = [vs[0], vs[1], 1];
        } else if vs.len() == 1 {
            e[j - 1] = [vs[0], 0, 0];
        }
        for (k, row) in fc.iter().enumerate() {
            let r = match row.iter().find(|x| x[0] as usize == j) {
                None => continue,
                Some(r) => r,
            };
            if r[1] > 0 {
                if f[j - 1][0] == 0 {
                    f[j - 1][0] = k + 1;
                }
            } else if r[1] < 0 && f[j - 1][1] == 0 {
                f[j - 1][1] = k + 1;
            }
        }
    }
    (e, f)
}


//! envelope_subdivision -- subdivide an assembled envelope (ONE atom pair) on a CONIC
//! subdivision (port of envelope_subdivision.py), and biconj_q's subdividePair around it.

use crate::alg_con::{assemble_alg_crv, AlgCon, EnvelopeMin, SubKind};
use crate::bail;
use crate::bridge::bridge_member;
use crate::cells::{assemble_qua_con_cells, cell_interior, cell_sample_point, side_i, Cell};
use crate::envelope::{envelope_eval, envelope_rational_corners, BFace, Env, FaceKind, QuadC};
use crate::envelope_subdivision_alg::envelope_subdivision_alg;
use crate::errors::R;
use crate::numeric::{det3, solve3};
use crate::qua_con::QuaCon;
use crate::rat_q::*;
use num_traits::Zero;

pub fn has_irrational_face(env: &Env) -> R<bool> {
    Ok(!envelope_rational_corners(env)?.bridge.is_empty())
}

pub fn unsubdivided(env: Env) -> AlgCon {
    AlgCon::EnvelopeMin(EnvelopeMin { pairs: vec![env], tris: vec![] })
}

pub fn subdivide_pair(env: Env, dom: &[Vec<Z>]) -> R<AlgCon> {
    if has_irrational_face(&env)? {
        match envelope_subdivision_alg(&[env.clone()], dom, &[], &Default::default()) {
            Ok(m) => {
                if m.iter().any(|x| x.kind == SubKind::Bridge) {
                    return Ok(AlgCon::AlgCrv(assemble_alg_crv(&m)?));
                }
            }
            Err(me) => {
                if me.id != "envelopeSubdivisionAlg:faceKind" && me.id != "envelopeSubdivisionAlg:arrangementSize" {
                    return Err(me);
                }
                return Ok(unsubdivided(env));
            }
        }
    }
    match envelope_subdivision(&env) {
        Ok(g) => Ok(AlgCon::QuaCon(g)),
        Err(me) => {
            if me.id == "envelopeSubdivision:bridgeCompetes" {
                return match envelope_subdivision_alg(&[env.clone()], &[], &[], &Default::default()) {
                    Ok(m) => Ok(AlgCon::AlgCrv(assemble_alg_crv(&m)?)),
                    Err(me2) => {
                        if me2.id != "envelopeSubdivisionAlg:faceKind" {
                            return Err(me2);
                        }
                        Ok(unsubdivided(env))
                    }
                };
            }
            if me.id == "envelopeSubdivision:bridgeFaces" || me.id == "envelopeSubdivision:arrangementSize" {
                return Ok(unsubdivided(env));
            }
            Err(me)
        }
    }
}

pub fn envelope_subdivision(env0: &Env) -> R<QuaCon> {
    let max_curves = 40;
    let env = envelope_rational_corners(env0)?;
    let mut corners: Vec<&BFace> = vec![];
    let mut bridges: Vec<&BFace> = vec![];
    for f in &env.bridge {
        if f.kind == FaceKind::Corner && !f.con.is_empty() {
            corners.push(f);
        } else if !f.cells.is_empty() && f.pair.is_some() {
            bridges.push(f);
        } else {
            bail!("envelopeSubdivision:bridgeFaces", "this envelope has a face whose region is a PREDICATE with no rows and no outer bound.");
        }
    }
    for f in &bridges {
        let w = match &f.wit {
            None => continue,
            Some(w) => w,
        };
        let x = [q_to_f64(&Q::new(w[0].clone(), w[2].clone())), q_to_f64(&Q::new(w[1].clone(), w[2].clone()))];
        if bridge_member(f.pair.as_ref().unwrap(), &[x])?[0] && envelope_eval(&env, &[x])?[0].is_finite() {
            bail!("envelopeSubdivision:bridgeCompetes", "a bridge face owns part of the envelope's domain, and its value competes with the quadratic faces there.");
        }
    }
    if env.quad.is_empty() {
        return Ok(QuaCon::empty());
    }
    let (mut c, req, pair_of) = curve_table(&env.quad, &corners);
    let n_quad = env.quad.len();
    add_bridge_curves(&mut c, &bridges);
    let mut all_pairs: Vec<usize> = pair_of.iter().flatten().cloned().chain(bridges.iter().map(|f| f.pair_index)).collect();
    all_pairs.sort();
    all_pairs.dedup();
    if c.len() > max_curves {
        bail!("envelopeSubdivision:arrangementSize", "the arrangement has {} curves (budget {}).", c.len(), max_curves);
    }
    let mut cells: Vec<Vec<(usize, i64)>> = vec![vec![]];
    for j in 1..=c.len() {
        let mut nxt = vec![];
        for cell in &cells {
            for sd in [1i64, -1] {
                let mut cand = cell.clone();
                cand.push((j, sd));
                if cell_interior(&c, &cand)? {
                    nxt.push(cand);
                }
            }
        }
        cells = nxt;
    }
    let mut out: Vec<Cell> = vec![];
    for cell in &cells {
        let app = applicable(&req, cell);
        let (pt, ok_pt) = interior_point(&c, cell)?;
        let mut here: Vec<usize> = app.iter().flat_map(|&i| pair_of[i - 1].iter().cloned()).collect();
        let mut bridge_here = false;
        if ok_pt {
            for f in &bridges {
                if bridge_member(f.pair.as_ref().unwrap(), &[pt])?[0] {
                    here.push(f.pair_index);
                    bridge_here = true;
                }
            }
        }
        if app.is_empty() && !bridge_here {
            continue;
        }
        if !all_pairs.iter().all(|p| here.contains(p)) {
            continue;
        }
        if bridge_here {
            bail!("envelopeSubdivision:bridgeCompetes", "a bridge face owns part of a cell, and its value competes with the quadratic faces there.");
        }
        let quad_app: Vec<usize> = app.iter().cloned().filter(|&i| i <= n_quad).collect();
        if quad_app.is_empty() {
            bail!("envelopeSubdivision:algebraicFace", "a cell is covered only by CORNER faces.");
        }
        let sub = winner_cells(&mut c, cell, &env.quad, &quad_app)?;
        out.extend(sub);
    }
    assemble_qua_con_cells(out)
}

pub fn curve_table(faces: &[QuadC], corners: &[&BFace]) -> (Vec<Vec<Z>>, Vec<Option<Vec<(usize, i64)>>>, Vec<Vec<usize>>) {
    let mut c: Vec<Vec<Z>> = vec![];
    let n = faces.len() + corners.len();
    let mut req = vec![None; n];
    let mut pair_of = vec![vec![]; n];
    for i in 0..n {
        let con = if i < faces.len() {
            pair_of[i] = faces[i].pair.clone();
            &faces[i].con
        } else {
            pair_of[i] = vec![corners[i - faces.len()].pair_index];
            &corners[i - faces.len()].con
        };
        let mut r = vec![];
        let mut dead = false;
        for row in con {
            let cc = &row[0..6];
            if all_zero(&cc[0..5]) {
                if !const_holds(&cc[5], side_i(&row[6])) {
                    dead = true;
                    break;
                }
                continue;
            }
            let br = if row.len() >= 8 { row[7].clone() } else { zi(0) };
            let mut key = cc.to_vec();
            key.push(br);
            let hit = add_key(&mut c, key);
            r.push((hit, side_i(&row[6])));
        }
        req[i] = if dead { None } else { Some(r) };
    }
    (c, req, pair_of)
}

pub fn const_holds(c6: &Z, side: i64) -> bool {
    (side > 0 && !c6.is_negative_z()) || (side < 0 && !c6.is_positive_z()) || (side == 0 && c6.is_zero())
}

trait SignZ {
    fn is_negative_z(&self) -> bool;
    fn is_positive_z(&self) -> bool;
}
impl SignZ for Z {
    fn is_negative_z(&self) -> bool {
        sign(self) < 0
    }
    fn is_positive_z(&self) -> bool {
        sign(self) > 0
    }
}

pub fn add_key(c: &mut Vec<Vec<Z>>, key: Vec<Z>) -> usize {
    match c.iter().position(|x| *x == key) {
        Some(h) => h + 1,
        None => {
            c.push(key);
            c.len()
        }
    }
}

fn add_bridge_curves(c: &mut Vec<Vec<Z>>, bridges: &[&BFace]) {
    for f in bridges {
        for con in &f.cells {
            for row in con {
                let cc = &row[0..6];
                if all_zero(&cc[0..5]) {
                    continue;
                }
                let mut key = cc.to_vec();
                key.push(if row.len() >= 8 { row[7].clone() } else { zi(0) });
                if !c.contains(&key) {
                    c.push(key);
                }
            }
        }
    }
}

/// Chebyshev centre of the straight rows (doubles), as numpy does it.
pub fn chebyshev(h: &[[f64; 3]], r: f64) -> (f64, [f64; 2]) {
    let mut m: Vec<[f64; 3]> = h.iter().map(|x| [x[0], x[1], -1.0]).collect();
    let mut b: Vec<f64> = h.iter().map(|x| -x[2]).collect();
    for (row, bb) in [([1.0, 0.0, 0.0], -r), ([-1.0, 0.0, 0.0], -r), ([0.0, 1.0, 0.0], -r), ([0.0, -1.0, 0.0], -r)] {
        m.push(row);
        b.push(bb);
    }
    let k = m.len();
    let mut best = f64::NEG_INFINITY;
    let mut cand = [0.0, 0.0];
    for i in 0..k - 2 {
        for j in i + 1..k - 1 {
            for l in j + 1..k {
                let s = [m[i], m[j], m[l]];
                if det3(&s).abs() < 1e-12 {
                    continue;
                }
                let z = match solve3(&s, [b[i], b[j], b[l]]) {
                    None => continue,
                    Some(z) => z,
                };
                if (0..k).all(|t| m[t][0] * z[0] + m[t][1] * z[1] + m[t][2] * z[2] >= b[t] - 1e-9) && z[2] > best {
                    best = z[2];
                    cand = [z[0], z[1]];
                }
            }
        }
    }
    (best, cand)
}

pub fn straight_h(c: &[Vec<Z>], rows: &[(usize, i64)]) -> Vec<[f64; 3]> {
    let mut h = vec![];
    for &(j, s) in rows {
        let cc = &c[j - 1];
        if cc[6].is_zero() && all_zero(&cc[0..3]) {
            let v = [to_f64(&(&cc[3] * zi(s))), to_f64(&(&cc[4] * zi(s))), to_f64(&(&cc[5] * zi(s)))];
            let n = v[0].hypot(v[1]);
            h.push(if n > 0.0 { [v[0] / n, v[1] / n, v[2] / n] } else { v });
        }
    }
    h
}

pub fn eval_curve(c: &[Vec<Z>], j: usize, p: [f64; 2]) -> R<f64> {
    let cc = &c[j - 1];
    if !cc[6].is_zero() {
        return pair_line_eval(&cc[0..6], side_i(&cc[6]), p[0], p[1]);
    }
    Ok(eval_conic(&cc[0..6], p[0], p[1]))
}

fn interior_point(c: &[Vec<Z>], rows: &[(usize, i64)]) -> R<([f64; 2], bool)> {
    let h = straight_h(c, rows);
    let r = 1e3;
    let cand = if h.is_empty() {
        [0.0, 0.0]
    } else {
        let (best, cand) = chebyshev(&h, r);
        if best <= 0.0 {
            return Ok(([0.0, 0.0], false));
        }
        cand
    };
    for &(j, s) in rows {
        let v = eval_curve(c, j, cand)?;
        if (s as f64) * v <= 0.0 {
            let (pt, ok, _) = cell_sample_point(c, rows, &h, r, false)?;
            return Ok((pt, ok));
        }
    }
    Ok((cand, true))
}

fn applicable(req: &[Option<Vec<(usize, i64)>>], rows: &[(usize, i64)]) -> Vec<usize> {
    let mut side = std::collections::HashMap::new();
    for r in rows {
        side.entry(r.0).or_insert(r.1);
    }
    let mut app = vec![];
    for (i, rq) in req.iter().enumerate() {
        if let Some(rq) = rq {
            if rq.iter().all(|t| side.get(&t.0) == Some(&t.1)) {
                app.push(i + 1);
            }
        }
    }
    app
}

fn winner_cells(c: &mut Vec<Vec<Z>>, rows: &[(usize, i64)], faces: &[QuadC], app: &[usize]) -> R<Vec<Cell>> {
    let mut open: Vec<(Vec<(usize, i64)>, usize)> = vec![(rows.to_vec(), app[0])];
    for &j in &app[1..] {
        let mut nxt = vec![];
        for cur in &open {
            let (d, sgn) = match diff_conic_of(&faces[j - 1], &faces[cur.1 - 1])? {
                None => {
                    nxt.push(cur.clone());
                    continue;
                }
                Some(x) => x,
            };
            let mut key = d;
            key.push(zi(0));
            let idx = add_key(c, key);
            for sd in [1i64, -1] {
                let mut rr = cur.0.clone();
                rr.push((idx, sd * sgn));
                let win = if sd > 0 { j } else { cur.1 };
                if cell_interior(c, &rr)? {
                    nxt.push((rr, win));
                }
            }
        }
        open = nxt;
    }
    let mut out = vec![];
    for (rr, win) in open {
        let f = &faces[win - 1];
        let con = rr.iter().map(|&(j, s)| {
            let mut r = c[j - 1][0..6].to_vec();
            r.push(zi(s));
            r.push(c[j - 1][6].clone());
            r
        }).collect();
        out.push(Cell { num: f.num.clone(), den: f.den.clone(), con });
    }
    Ok(out)
}

pub fn diff_conic_of(fj: &QuadC, fi: &QuadC) -> R<Option<(Vec<Z>, i64)>> {
    let num: Vec<Z> = fj.num.iter().zip(fi.num.iter()).map(|(a, b)| a * &fi.den - b * &fj.den).collect();
    let raw = face_conic(&num);
    if all_zero(&raw) || !all_zero(&num[0..4]) {
        return Ok(None);
    }
    let nz = raw.iter().find(|x| !x.is_zero()).unwrap();
    Ok(Some((conic(&raw)?, sign(nz) as i64)))
}

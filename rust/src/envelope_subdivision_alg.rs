//! envelope_subdivision_alg -- subdivide an assembled envelope when a BRIDGE (or irrational
//! corner) face competes with the quadratic faces (port of envelope_subdivision_alg.py).

use crate::memo::Memo;
use std::sync::LazyLock;
use crate::alg_branch::{alg_branch, alg_branch_info};
use crate::alg_con::{AlgRow, IiRow, SubFace, SubKind};
use crate::alg_face_below::{alg_face_below, BelowFace};
use crate::alg_row_fails::{alg_row_fails, Region, RowCache, RowSpec};
use crate::bail;
use crate::biv::M2;
use crate::bridge::bridge_member;
use crate::cells::{cell_sample_point, side_i};
use crate::conic_meet::solve2f;
use crate::envelope::{envelope_rational_corners, BFace, Env, FaceKind, PairData, QuadC};
use crate::envelope_subdivision::{add_key, chebyshev, eval_curve, straight_h};
use crate::errors::R;
use crate::numeric::{linspace, round_half_away};
use crate::poly_n::PolyN;
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::{BTreeSet, HashMap};
use std::sync::Arc;

pub const PAIR_STRIDE: usize = 1_000_000;

#[derive(Clone, Debug)]
pub struct AlgOpts {
    pub overlay_budget: f64,
    pub max_cells: usize,
}

impl Default for AlgOpts {
    fn default() -> Self {
        AlgOpts { overlay_budget: 3000.0, max_cells: 3000 }
    }
}

pub struct SubInfo {
    pub cell_con: Vec<Vec<Vec<Z>>>,
    pub face_cell: Vec<usize>,
}

type Rows = Vec<(usize, i64)>;

#[derive(Clone)]
struct Irr {
    fp: Arc<PolyN>,
    fb: i64,
}

pub fn envelope_subdivision_alg(envs: &[Env], dom: &[Vec<Z>], cuts: &[Vec<Z>], opts: &AlgOpts) -> R<Vec<SubFace>> {
    Ok(esa(envs, dom, cuts, opts)?.0)
}

pub fn esa(envs: &[Env], dom: &[Vec<Z>], cuts: &[Vec<Z>], opts: &AlgOpts) -> R<(Vec<SubFace>, SubInfo)> {
    if envs.len() > 1 {
        // one subdivision per envelope: independent
        let res = crate::par::map(envs, |e| esa(&[e.clone()], dom, cuts, opts))?;
        let mut subs = vec![];
        let mut infos = vec![];
        for (sb, inf) in res {
            subs.push(sb);
            infos.push(inf);
        }
        return Ok((min_overlay(&subs, &infos, opts.overlay_budget)?, SubInfo { cell_con: vec![], face_cell: vec![] }));
    }
    let env = combine_pairs(envs)?;
    let (corners, bridges) = split_bridge(&env.bridge)?;
    let (mut c, req, pair_of) = conic_table(&env.quad, &corners);
    let n_quad = env.quad.len();
    let mut all_pairs: Vec<usize> = pair_of.iter().flatten().cloned().chain(bridges.iter().map(|f| f.pair_index)).collect();
    all_pairs.sort();
    all_pairs.dedup();
    let c_ends = add_corner_crossings(&mut c, &corners)?;
    let b_ends = add_bridge_curves(&mut c, &bridges)?;
    let l_ends = split_line_curves(&mut c)?;
    for cu in cuts {
        let mut k = conic_lin(cu)?;
        k.push(zi(0));
        add_key(&mut c, k);
    }
    let mut owners = curve_owners(&c, &req, n_quad, &bridges, &c_ends, &b_ends, &l_ends);
    let mut cells: Vec<Rows> = vec![vec![]];
    if !dom.is_empty() {
        let mut base = vec![];
        for drow in dom {
            let mut raw = vec![zi(0), zi(0), zi(0)];
            raw.extend(drow.iter().cloned());
            let mut key = conic(&raw)?;
            key.push(zi(0));
            let hit = add_key(&mut c, key);
            base.push((hit, sign(raw.iter().find(|x| !x.is_zero()).unwrap()) as i64));
        }
        cells = vec![base];
        while owners.len() < c.len() {
            owners.push(vec![vec![]]);
        }
    }
    let mut skips: Vec<Vec<usize>> = vec![vec![]];
    for j in 1..=c.len() {
        // every cell is tested against curve j independently; the outcomes are assembled in cell order
        let owners_j = &owners[j - 1];
        let outcomes = crate::par::map_range(cells.len(), |ci| -> R<Vec<(Rows, Vec<usize>)>> {
            let (cell, sk) = (&cells[ci], &skips[ci]);
            if cell.iter().any(|r| r.0 == j) {
                return Ok(vec![(cell.clone(), sk.clone())]);
            }
            let mut cp = cell.clone();
            cp.push((j, 1));
            let mut cm = cell.clone();
            cm.push((j, -1));
            let live_p = live_cell(&c, &cp)?;
            let live_m = live_cell(&c, &cm)?;
            if live_p && live_m && !relevant_here(&c, owners_j, cell)? {
                let mut s2 = sk.clone();
                s2.push(j);
                Ok(vec![(cell.clone(), s2)])
            } else if live_p && live_m {
                Ok(vec![(cp, sk.clone()), (cm, sk.clone())])
            } else if live_p || live_m {
                Ok(vec![(cell.clone(), sk.clone())])
            } else {
                Ok(vec![])
            }
        })?;
        let mut nxt = vec![];
        let mut nxt_sk = vec![];
        for o in outcomes {
            for (r, k) in o {
                nxt.push(r);
                nxt_sk.push(k);
            }
        }
        cells = nxt;
        skips = nxt_sk;
        if cells.len() > opts.max_cells {
            bail!("envelopeSubdivisionAlg:arrangementSize", "the arrangement passed {} cells at curve {} of {}.", opts.max_cells, j, c.len());
        }
    }
    // the final loop: every cell gets its faces independently (this is where `drop_below` spends
    // its time); a cell with no interior point is None, a cell that gets no faces is Some with none
    let per_cell = crate::par::map_range(cells.len(), |ci| -> R<Option<(Vec<Vec<Z>>, Vec<SubFace>)>> {
        let (rows, sk) = (&cells[ci], &skips[ci]);
        let (pt, ok_pt) = interior_point(&c, rows)?;
        if !ok_pt {
            return Ok(None);
        }
        let con = rows_to_con(&c, rows);
        let app = applicable(&req, rows, &c, pt)?;
        let mut here: Vec<usize> = app.iter().flat_map(|&i| pair_of[i - 1].iter().cloned()).collect();
        let mut b_here = vec![];
        let (q, pw) = vote_points(&c, rows, pt)?;
        let rows_out = full_rows(&c, rows, &q, sk)?;
        for (b, f) in bridges.iter().enumerate() {
            let m = bridge_member(f.pair.as_ref().unwrap(), &q)?;
            if m.iter().filter(|x| **x).count() * 2 > q.len() {
                here.push(f.pair_index);
                b_here.push(b);
            }
        }
        if !all_pairs.iter().all(|p| here.contains(p)) {
            return Ok(Some((con, vec![])));
        }
        let quad_app: Vec<usize> = app.iter().cloned().filter(|&i| i <= n_quad).collect();
        let mut irr: Vec<Irr> = b_here.iter().map(|&b| Irr { fp: bridges[b].fp.clone(), fb: bridges[b].fb }).collect();
        for &j in app.iter().filter(|&&i| i > n_quad) {
            if let Some(fc) = corner_value_face(corners[j - n_quad - 1], pt)? {
                irr.push(fc);
            }
        }
        if quad_app.is_empty() && irr.is_empty() {
            return Ok(Some((con, vec![])));
        }
        if irr.len() > 1 {
            irr = drop_below(irr, &c, rows, &pw, &env.quad, &quad_app)?;
        }
        let sub = winner_faces(&env.quad, &quad_app, &irr, rows, &c, &pw, &rows_out)?;
        Ok(Some((con, sub)))
    })?;
    let mut subdivision = vec![];
    let mut cell_con = vec![];
    let mut face_cell = vec![];
    for r in per_cell.into_iter().flatten() {
        cell_con.push(r.0);
        for sf in r.1 {
            subdivision.push(sf);
            face_cell.push(cell_con.len());
        }
    }
    Ok((subdivision, SubInfo { cell_con, face_cell }))
}

fn quad_face(q: &QuadC, con: Vec<Vec<Z>>, alg: Vec<AlgRow>) -> SubFace {
    SubFace { kind: SubKind::Quad, num: q.num.clone(), den: q.den.clone(), fp: None, fb: 0, con, alg, ii: vec![], tri: vec![], tf: None, tfs: vec![] }
}

fn winner_faces(quad: &[QuadC], quad_app: &[usize], irr: &[Irr], rows: &Rows, c: &[Vec<Z>], q: &[[f64; 2]], rows_out: &Rows) -> R<Vec<SubFace>> {
    let mut out = vec![];
    let base = rows_to_con(c, rows_out);
    if irr.is_empty() {
        // every applicable quadratic is tested independently
        let faces = crate::par::map(quad_app, |&j| -> R<Option<SubFace>> {
            let extra = quad_over(quad, quad_app, j)?;
            if !piece_live(c, rows, &extra, q)? {
                return Ok(None);
            }
            let mut con = base.clone();
            con.extend(extra);
            Ok(Some(quad_face(&quad[j - 1], con, vec![])))
        })?;
        out.extend(faces.into_iter().flatten());
        return Ok(out);
    }
    for f in irr {
        let alg: Vec<AlgRow> = quad_app.iter().map(|&j| AlgRow { fp: f.fp.clone(), fb: f.fb, num: quad[j - 1].num.clone(), den: quad[j - 1].den.clone(), side: -1 }).collect();
        out.push(SubFace { kind: SubKind::Bridge, num: vec![], den: zi(0), fp: Some(f.fp.clone()), fb: f.fb, con: base.clone(), alg, ii: vec![], tri: vec![], tf: None, tfs: vec![] });
    }
    let faces = crate::par::map(quad_app, |&j| -> R<Option<SubFace>> {
        let extra = quad_over(quad, quad_app, j)?;
        if !piece_live(c, rows, &extra, q)? {
            return Ok(None);
        }
        let alg: Vec<AlgRow> = irr.iter().map(|f| AlgRow { fp: f.fp.clone(), fb: f.fb, num: quad[j - 1].num.clone(), den: quad[j - 1].den.clone(), side: 1 }).collect();
        if alg_empty(c, rows_out, &extra, &alg, q)? {
            return Ok(None);
        }
        let mut con = base.clone();
        con.extend(extra);
        Ok(Some(quad_face(&quad[j - 1], con, alg)))
    })?;
    out.extend(faces.into_iter().flatten());
    Ok(out)
}

fn full_rows(c: &[Vec<Z>], rows: &Rows, p: &[[f64; 2]], skipped: &[usize]) -> R<Rows> {
    let mut rows = rows.clone();
    let have: BTreeSet<usize> = rows.iter().map(|r| r.0).chain(skipped.iter().cloned()).collect();
    for j in 1..=c.len() {
        if have.contains(&j) {
            continue;
        }
        let mut sum = 0i64;
        let mut any = false;
        for q in p {
            let s = signf(eval_curve(c, j, *q)?);
            if s != 0 {
                any = true;
                sum += s as i64;
            }
        }
        if !any {
            continue;
        }
        rows.push((j, sum.signum()));
    }
    Ok(rows)
}

fn vote_points(c: &[Vec<Z>], rows: &Rows, pt: [f64; 2]) -> R<(Vec<[f64; 2]>, Vec<[f64; 2]>)> {
    let mut q = vec![pt];
    let (_, ok, p) = cell_sample_point(c, rows, &straight_h(c, rows), 1e3, false)?;
    if ok && p.len() >= 4 {
        let idx: Vec<usize> = linspace(1.0, p.len() as f64, 4).iter().map(|t| round_half_away(*t) as usize).collect();
        q = vec![pt];
        q.extend(idx.iter().map(|i| p[i - 1]));
    }
    let mut pw = vec![pt];
    pw.extend(p);
    Ok((q, pw))
}

pub fn live_cell(c: &[Vec<Z>], rows: &Rows) -> R<bool> {
    let (tf, unsure) = crate::cells::cell_interior_info(c, rows)?;
    if !tf {
        return Ok(false);
    }
    let has_branch = rows.iter().any(|r| !c[r.0 - 1][6].is_zero());
    let has_curved = rows.iter().any(|r| c[r.0 - 1][6].is_zero() && !all_zero(&c[r.0 - 1][0..3]));
    if unsure || (has_branch && has_curved) {
        let (_, tf2, _) = cell_sample_point(c, rows, &straight_h(c, rows), 1e3, false)?;
        return Ok(tf2);
    }
    Ok(tf)
}

fn conic_holds_f(r: &[Z], q: [f64; 2]) -> f64 {
    to_f64(&r[6]) * eval_conic(&r[0..6], q[0], q[1])
}

fn alg_empty(c: &[Vec<Z>], rows: &Rows, con_rows: &[Vec<Z>], alg_rows: &[AlgRow], q: &[[f64; 2]]) -> R<bool> {
    if alg_rows.is_empty() {
        return Ok(false);
    }
    let qk: Vec<[f64; 2]> = q.iter().cloned().filter(|p| con_rows.iter().all(|r| conic_holds_f(r, *p) > 0.0)).collect();
    let mut holds = vec![true; qk.len()];
    let mut n_fail = vec![];
    for a in alg_rows {
        let z = alg_branch(&a.fp, &qk, a.fb)?;
        let ok: Vec<bool> = qk.iter().zip(z.iter()).map(|(p, zi)| {
            let qi = eval_face(&a.num, &a.den, p[0], p[1]);
            zi.is_finite() && (a.side as f64) * (qi - zi) > 1e-9 * 1f64.max(zi.abs())
        }).collect();
        for (h, o) in holds.iter_mut().zip(ok.iter()) {
            *h = *h && *o;
        }
        n_fail.push(ok.iter().filter(|o| !**o).count());
    }
    if holds.iter().any(|h| *h) {
        return Ok(false);
    }
    for (t, a) in alg_rows.iter().enumerate() {
        if qk.is_empty() || n_fail[t] < qk.len() {
            continue;
        }
        if row_fails(c, rows, con_rows, &[a.clone()])? {
            return Ok(true);
        }
    }
    if alg_rows.len() > 1 || qk.is_empty() {
        let mut o: Vec<usize> = (0..alg_rows.len()).collect();
        o.sort_by_key(|&t| std::cmp::Reverse(n_fail[t]));
        return row_fails(c, rows, con_rows, &o.iter().map(|&t| alg_rows[t].clone()).collect::<Vec<_>>());
    }
    Ok(false)
}

type RowFailsKey = (Vec<(Vec<Z>, i64)>, Vec<Vec<Z>>, Vec<(Vec<(Vec<u32>, Z)>, Vec<Z>, Z, i64, i64)>);

static ROW_FAILS_CACHE: LazyLock<Memo<RowFailsKey, bool>> = LazyLock::new(|| Memo::new(1000000));

pub fn clear_memos() {
    ROW_FAILS_CACHE.clear();
}

fn row_fails(c: &[Vec<Z>], rows: &Rows, con_rows: &[Vec<Z>], a: &[AlgRow]) -> R<bool> {
    let mut k1: Vec<(Vec<Z>, i64)> = rows.iter().map(|r| (c[r.0 - 1].clone(), r.1)).collect();
    k1.sort();
    let mut k2: Vec<Vec<Z>> = con_rows.to_vec();
    k2.sort();
    let k3: Vec<(Vec<(Vec<u32>, Z)>, Vec<Z>, Z, i64, i64)> = a.iter().map(|x| (x.fp.iter().map(|(e, c)| (e.clone(), c.clone())).collect(), x.num.clone(), x.den.clone(), x.fb, x.side)).collect();
    let key: RowFailsKey = (k1, k2, k3);
    if let Some(v) = ROW_FAILS_CACHE.get(&key) {
        return Ok(v);
    }
    let mut ok = false;
    if let Some(reg) = region_box(c, rows, con_rows)? {
        let specs: Vec<RowSpec> = a.iter().map(|x| RowSpec::Quad { fp: x.fp.clone(), fb: x.fb, num: x.num.clone(), den: x.den.clone(), side: x.side }).collect();
        ok = alg_row_fails(&specs, &reg, 6).unwrap_or(false);
    }
    ROW_FAILS_CACHE.insert(key, ok);
    Ok(ok)
}

/// The cell's bounding box, its region predicate pieces, its row polynomials and corners.
pub fn region_box(c: &[Vec<Z>], rows: &Rows, k: &[Vec<Z>]) -> R<Option<Region>> {
    let h = straight_h(c, rows);
    let mut a: Vec<[f64; 3]> = h.clone();
    a.extend([[1.0, 0.0, 1e3], [-1.0, 0.0, 1e3], [0.0, 1.0, 1e3], [0.0, -1.0, 1e3]]);
    a.extend(branch_lines(c, rows, k)?);
    let mut v = vec![];
    for i in 0..a.len() {
        for j in i + 1..a.len() {
            let mx = [[a[i][0], a[i][1]], [a[j][0], a[j][1]]];
            if (mx[0][0] * mx[1][1] - mx[0][1] * mx[1][0]).abs() < 1e-12 {
                continue;
            }
            let x = solve2f(&mx, [-a[i][2], -a[j][2]]);
            let lim = -1e-9 * 1f64.max((x[0] * x[0] + x[1] * x[1]).sqrt());
            if a.iter().all(|aa| aa[0] * x[0] + aa[1] * x[1] + aa[2] >= lim) {
                v.push(x);
            }
        }
    }
    if v.len() < 3 {
        return Ok(None);
    }
    let bx = [
        v.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min),
        v.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max),
        v.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min),
        v.iter().map(|p| p[1]).fold(f64::NEG_INFINITY, f64::max),
    ];
    let mut row_polys: Vec<M2> = vec![];
    for r in rows {
        if !c[r.0 - 1][6].is_zero() {
            continue;
        }
        row_polys.push(conic_poly(&c[r.0 - 1][0..6], r.1));
    }
    for r in k {
        if r.len() >= 8 && !r[7].is_zero() {
            continue;
        }
        row_polys.push(conic_poly(&r[0..6], side_i(&r[6])));
    }
    let mut lines: Vec<(Vec<Z>, i64, i64)> = vec![];
    for r in rows {
        if !c[r.0 - 1][6].is_zero() {
            lines.push((c[r.0 - 1][0..6].to_vec(), r.1, side_i(&c[r.0 - 1][6])));
        }
    }
    for r in k {
        if r.len() >= 8 && !r[7].is_zero() {
            lines.push((r[0..6].to_vec(), side_i(&r[6]), side_i(&r[7])));
        }
    }
    Ok(Some(Region::new(bx, row_polys, lines, v)?))
}

fn branch_lines(c: &[Vec<Z>], rows: &Rows, k: &[Vec<Z>]) -> R<Vec<[f64; 3]>> {
    let mut l: Vec<[f64; 3]> = vec![];
    let line = |cc: &[Z], br: i64, sd: i64| -> R<[f64; 3]> {
        let pe = pair_eval_prep(cc, br)?;
        let s = sd as f64;
        let v = [s * pe.eval_m(0.0, 0.0).0, s * pe.eval_m(1.0, 0.0).0, s * pe.eval_m(0.0, 1.0).0];
        Ok([v[1] - v[0], v[2] - v[0], v[0]])
    };
    for r in rows {
        let br = side_i(&c[r.0 - 1][6]);
        if br == 0 {
            continue;
        }
        l.push(line(&c[r.0 - 1][0..6], br, r.1)?);
    }
    for r in k {
        if r.len() >= 8 && !r[7].is_zero() {
            l.push(line(&r[0..6], side_i(&r[7]), side_i(&r[6]))?);
        }
    }
    let mut out = vec![];
    for x in l {
        if !x.iter().all(|v| v.is_finite()) {
            continue;
        }
        let nr = x[0].hypot(x[1]);
        if nr > 0.0 {
            out.push([x[0] / nr, x[1] / nr, x[2] / nr]);
        }
    }
    Ok(out)
}

/// The conic as a 3x3 coefficient matrix P[i][j] of x^i y^j, times s.
pub fn conic_poly(c: &[Z], s: i64) -> M2 {
    let s = zi(s);
    vec![
        vec![&s * &c[5], &s * &c[4], &s * &c[2]],
        vec![&s * &c[3], &s * &c[1], zi(0)],
        vec![&s * &c[0], zi(0), zi(0)],
    ]
}

fn piece_live(c: &[Vec<Z>], rows: &Rows, extra: &[Vec<Z>], q: &[[f64; 2]]) -> R<bool> {
    let k0: Vec<bool> = extra.iter().map(|r| all_zero(&r[0..5])).collect();
    let mut extra: Vec<Vec<Z>> = extra.to_vec();
    if k0.iter().any(|x| *x) {
        if extra.iter().zip(k0.iter()).any(|(r, z)| *z && (&r[6] * &r[5]).is_negative()) {
            return Ok(false);
        }
        extra = extra.into_iter().zip(k0.iter()).filter(|(_, z)| !**z).map(|(r, _)| r).collect();
    }
    if extra.is_empty() {
        return Ok(true);
    }
    for p in q {
        if extra.iter().all(|r| {
            let m = to_f64(&r[0..6].iter().map(|x| x.abs()).max().unwrap());
            conic_holds_f(r, *p) > 1e-9 * m
        }) {
            return Ok(true);
        }
    }
    let fails: Vec<usize> = extra.iter().map(|r| q.iter().filter(|p| conic_holds_f(r, **p) <= 0.0).count()).collect();
    let mut o: Vec<usize> = (0..extra.len()).collect();
    o.sort_by_key(|&t| std::cmp::Reverse(fails[t]));
    for &t in o.iter().take(3) {
        let r = &extra[t];
        let mut cc = c.to_vec();
        let mut key = r[0..6].to_vec();
        key.push(r[7].clone());
        cc.push(key);
        let mut rr = rows.clone();
        rr.push((c.len() + 1, side_i(&r[6])));
        if !live_cell(&cc, &rr)? {
            return Ok(false);
        }
    }
    let mut cc = c.to_vec();
    let mut rr = rows.clone();
    for (t, r) in extra.iter().enumerate() {
        let mut key = r[0..6].to_vec();
        key.push(r[7].clone());
        cc.push(key);
        rr.push((c.len() + 1 + t, side_i(&r[6])));
    }
    live_cell(&cc, &rr)
}

fn quad_over(quad: &[QuadC], quad_app: &[usize], j: usize) -> R<Vec<Vec<Z>>> {
    let mut extra = vec![];
    for &i in quad_app {
        if i == j {
            continue;
        }
        if let Some((c, sgn)) = diff_conic(&quad[j - 1].num, &quad[j - 1].den, &quad[i - 1].num, &quad[i - 1].den)? {
            let mut r = c;
            r.push(zi(sgn));
            r.push(zi(0));
            extra.push(r);
        }
    }
    Ok(extra)
}

/// (conic, sign) of fw - fj, or None when identical / cubic.
fn diff_conic(wn: &[Z], wd: &Z, jn: &[Z], jd: &Z) -> R<Option<(Vec<Z>, i64)>> {
    let num: Vec<Z> = wn.iter().zip(jn.iter()).map(|(a, b)| a * jd - b * wd).collect();
    let raw = face_conic(&num);
    if all_zero(&raw) || !all_zero(&num[0..4]) {
        return Ok(None);
    }
    let nz = raw.iter().find(|x| !x.is_zero()).unwrap();
    Ok(Some((conic(&raw)?, sign(nz) as i64)))
}

fn add_corner_crossings(c: &mut Vec<Vec<Z>>, corners: &[&BFace]) -> R<Vec<(Vec<Z>, usize)>> {
    let mut ends = vec![];
    for (j, cf) in corners.iter().enumerate() {
        let sg = match &cf.sigma {
            None => continue,
            Some(s) => s,
        };
        if sg.root == 0 {
            continue;
        }
        let src = cf.src.as_ref().unwrap();
        let n = &src.n;
        let d = &src.d;
        let nn = vec![vec![n[4].clone(), n[5].clone()], vec![n[5].clone(), n[6].clone()]];
        let m = vec![n[7].clone(), n[8].clone()];
        let (al, be) = (&sg.quad[0], &sg.quad[1]);
        let t1 = zi(2) * al * quad(&sg.s0n, &nn, &sg.w) - be * quad(&sg.w, &nn, &sg.w) + zi(2) * &sg.na * al * dot(&m, &sg.w);
        let e1 = vec![zi(2) * al * &sg.na * d * &sg.w[0], zi(2) * al * &sg.na * d * &sg.w[1], -t1];
        if e1[0].is_zero() && e1[1].is_zero() {
            continue;
        }
        let mut key = conic_lin(&e1)?;
        key.push(zi(0));
        add_key(c, key.clone());
        ends.push((key, j));
    }
    Ok(ends)
}

fn corner_value_face(cf: &BFace, pt: [f64; 2]) -> R<Option<Irr>> {
    let src = cf.src.as_ref().unwrap();
    let n: Vec<f64> = src.n.iter().map(to_f64).collect();
    let s = cf.s;
    let z = pt[0] * s[0] + pt[1] * s[1] - ((s[0] * (n[4] * s[0] + n[5] * s[1]) + s[1] * (n[5] * s[0] + n[6] * s[1])) / 2.0 + n[7] * s[0] + n[8] * s[1] + n[9]) / to_f64(&src.d);
    let (_, all) = alg_branch_info(&cf.fp, &[pt], -1)?;
    let rr = &all[0];
    if rr.is_empty() {
        return Ok(None);
    }
    let gaps: Vec<f64> = rr.iter().map(|r| (r - z).abs()).collect();
    let gap = gaps.iter().cloned().fold(f64::INFINITY, f64::min);
    let j = gaps.iter().position(|g| *g == gap).unwrap();
    if gap > 1e-6 * 1f64.max(z.abs()) {
        return Ok(None);
    }
    Ok(Some(Irr { fp: cf.fp.clone(), fb: (j + 1) as i64 - rr.len() as i64 - 1 }))
}

fn same_irr(a: &Irr, b: &Irr) -> bool {
    a.fb == b.fb && (Arc::ptr_eq(&a.fp, &b.fp) || *a.fp == *b.fp)
}

fn drop_below(irr: Vec<Irr>, c: &[Vec<Z>], rows: &Rows, p: &[[f64; 2]], quad: &[QuadC], quad_app: &[usize]) -> R<Vec<Irr>> {
    let mut keep = vec![true; irr.len()];
    for i in 1..irr.len() {
        for j in 0..i {
            if keep[j] && same_irr(&irr[i], &irr[j]) {
                keep[i] = false;
                break;
            }
        }
    }
    let irr: Vec<Irr> = irr.into_iter().zip(keep).filter(|(_, k)| *k).map(|(f, _)| f).collect();
    let n = irr.len();
    if n < 2 {
        return Ok(irr);
    }
    let mut v: Vec<Vec<f64>> = vec![];
    for f in &irr {
        let z = alg_branch(&f.fp, p, f.fb)?;
        v.push(z.into_iter().map(|x| if x.is_finite() { x } else { f64::NEG_INFINITY }).collect());
    }
    let mut top = vec![];
    for t in 0..p.len() {
        let col: Vec<f64> = (0..n).map(|i| v[i][t]).collect();
        let mx = col.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        top.push(col.iter().position(|x| *x == mx).unwrap_or(0));
    }
    let wins: Vec<usize> = (0..n).map(|i| top.iter().filter(|t| **t == i).count()).collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| std::cmp::Reverse(wins[i]));
    let above: Vec<usize> = order[..n.min(2)].to_vec();
    let reg = match region_box(c, rows, &[])? {
        None => return Ok(irr),
        Some(r) => r,
    };
    let own: Vec<Vec<AlgRow>> = irr.iter().map(|f| quad_app.iter().map(|&j| AlgRow { fp: f.fp.clone(), fb: f.fb, num: quad[j - 1].num.clone(), den: quad[j - 1].den.clone(), side: -1 }).collect()).collect();
    // The drop decision of each face is taken against ALL the other candidates, not against the
    // ones that were still kept at that moment as in the old ordered loops. "lo is below hi" is a
    // certificate of a pointwise order that is monotone in hi (adding faces only raises the max)
    // and transitive (lo <= hi and hi <= hi' imply lo <= hi'), so for distinct faces the kept set
    // is the same as the old ordered tests gave (checked on the 200-fixture digest); and the
    // decisions are now independent per face, so the two loops run in parallel and their dependent
    // chain is gone.
    let cache = RowCache::new();
    let drop_top: Vec<bool> = crate::par::map_range(n, |i| {
        for &j in &above {
            if j == i {
                continue;
            }
            if v[i].iter().zip(v[j].iter()).all(|(a, b)| *a <= b + 1e-9 * 1f64.max(b.abs())) {
                let lo = BelowFace { fp: irr[i].fp.clone(), fb: irr[i].fb };
                let hi = vec![BelowFace { fp: irr[j].fp.clone(), fb: irr[j].fb }];
                if alg_face_below(&lo, &hi, &reg, 5, &own[i], Some(&cache))? {
                    return Ok(true);
                }
            }
        }
        Ok(false)
    })?;
    let kept: Vec<usize> = (0..n).filter(|&i| !drop_top[i]).collect();
    if kept.len() < 2 {
        return Ok(kept.into_iter().map(|i| irr[i].clone()).collect());
    }
    let drop_rest: Vec<bool> = crate::par::map_range(n, |i| {
        if drop_top[i] || wins[i] != 0 {
            return Ok(false);
        }
        let hi: Vec<BelowFace> = kept.iter().filter(|&&j| j != i).map(|&j| BelowFace { fp: irr[j].fp.clone(), fb: irr[j].fb }).collect();
        alg_face_below(&BelowFace { fp: irr[i].fp.clone(), fb: irr[i].fb }, &hi, &reg, 5, &own[i], Some(&cache))
    })?;
    Ok(kept.into_iter().filter(|&i| !drop_rest[i]).map(|i| irr[i].clone()).collect())
}

fn combine_pairs(envs: &[Env]) -> R<Env> {
    let mut env = Env::default();
    for (g, e0) in envs.iter().enumerate() {
        let e = envelope_rational_corners(e0)?;
        let off = g * PAIR_STRIDE;
        for q in &e.quad {
            env.quad.push(QuadC { num: q.num.clone(), den: q.den.clone(), con: q.con.clone(), pair: q.pair.iter().map(|p| p + off).collect() });
        }
        for b in &e.bridge {
            let mut b = b.clone();
            b.pair_index += off;
            env.bridge.push(b);
        }
    }
    Ok(env)
}

fn split_bridge(bl: &[BFace]) -> R<(Vec<&BFace>, Vec<&BFace>)> {
    let mut corners = vec![];
    let mut bridges = vec![];
    for f in bl {
        if f.kind == FaceKind::Corner && !f.con.is_empty() {
            corners.push(f);
        } else if f.kind != FaceKind::Corner && f.pair.is_some() {
            bridges.push(f);
        } else {
            bail!("envelopeSubdivisionAlg:faceKind", "a face of kind {:?} has neither rows nor an implicit value.", f.kind);
        }
    }
    Ok((corners, bridges))
}

type Req = Vec<Option<Rows>>;

fn conic_table(faces: &[QuadC], corners: &[&BFace]) -> (Vec<Vec<Z>>, Req, Vec<Vec<usize>>) {
    crate::envelope_subdivision::curve_table(faces, corners)
}

fn applicable(req: &Req, rows: &Rows, c: &[Vec<Z>], pt: [f64; 2]) -> R<Vec<usize>> {
    let mut side: HashMap<usize, i64> = HashMap::new();
    for r in rows {
        side.entry(r.0).or_insert(r.1);
    }
    let mut app = vec![];
    for (i, rq) in req.iter().enumerate() {
        let rq = match rq {
            None => continue,
            Some(r) => r,
        };
        let mut ok = true;
        for t in rq {
            let s = match side.get(&t.0) {
                Some(s) => *s,
                None => signf(eval_curve(c, t.0, pt)?) as i64,
            };
            if s != t.1 {
                ok = false;
                break;
            }
        }
        if ok {
            app.push(i + 1);
        }
    }
    Ok(app)
}

fn add_bridge_curves(c: &mut Vec<Vec<Z>>, bridges: &[&BFace]) -> R<Vec<(Vec<Z>, usize)>> {
    let mut ends = vec![];
    for (b, f) in bridges.iter().enumerate() {
        for con in &f.cells {
            for row in con {
                if all_zero(&row[0..5]) {
                    continue;
                }
                let mut key = row[0..6].to_vec();
                key.push(if row.len() >= 8 { row[7].clone() } else { zi(0) });
                add_key(c, key);
            }
        }
        let pr = match &**f.pair.as_ref().unwrap() {
            PairData::Alg(a) => a,
            PairData::Rat(_) => continue,
        };
        for k in [1, 2] {
            let (nk, mk) = if k == 1 { (&pr.n1, &pr.m1) } else { (&pr.n2, &pr.m2) };
            if nk.iter().flatten().all(|x| x.is_zero()) || !det_exact(nk)?.is_zero() {
                continue;
            }
            let col = (0..2).find(|&j| !nk[0][j].is_zero() || !nk[1][j].is_zero()).unwrap();
            let u = [nk[0][col].clone(), nk[1][col].clone()];
            let w = [-&u[1], u[0].clone()];
            let raw = vec![zi(0), zi(0), zi(0), &pr.d * &w[0], &pr.d * &w[1], -(&w[0] * &mk[0] + &w[1] * &mk[1])];
            let mut key = conic(&raw)?;
            key.push(zi(0));
            add_key(c, key.clone());
            ends.push((key, b));
        }
    }
    Ok(ends)
}

fn split_line_curves(c: &mut Vec<Vec<Z>>) -> R<Vec<(Vec<Z>, usize)>> {
    let mut ends = vec![];
    let n = c.len();
    for j in 0..n {
        let cc = c[j][0..6].to_vec();
        if !c[j][6].is_zero() || all_zero(&cc[0..3]) {
            continue;
        }
        if let Some(l) = square_root(&cc) {
            let mut key = conic_lin(&l)?;
            key.push(zi(0));
            add_key(c, key.clone());
            ends.push((key, j + 1));
        }
    }
    Ok(ends)
}

fn curve_owners(c: &[Vec<Z>], req: &Req, n_quad: usize, bridges: &[&BFace], c_ends: &[(Vec<Z>, usize)], b_ends: &[(Vec<Z>, usize)], l_ends: &[(Vec<Z>, usize)]) -> Vec<Vec<Rows>> {
    let n = c.len();
    let mut owners: Vec<Vec<Rows>> = vec![vec![]; n];
    for r in req.iter().flatten() {
        if r.is_empty() {
            continue;
        }
        for t in r {
            owners[t.0 - 1].push(r.iter().cloned().filter(|x| x.0 != t.0).collect());
        }
    }
    for f in bridges {
        for r in bridge_blocks(c, f) {
            for u in &r {
                owners[u.0 - 1].push(r.iter().cloned().filter(|x| x.0 != u.0).collect());
            }
        }
    }
    let index_of = |key: &Vec<Z>| c.iter().position(|x| x == key).unwrap();
    for (key, cj) in c_ends {
        let j = index_of(key);
        if let Some(r) = &req[n_quad + cj] {
            owners[j].push(r.iter().cloned().filter(|x| x.0 != j + 1).collect());
        }
    }
    for (key, b) in b_ends {
        let j = index_of(key);
        for blk in bridge_blocks(c, bridges[*b]) {
            owners[j].push(blk.iter().cloned().filter(|x| x.0 != j + 1).collect());
        }
    }
    for (key, src) in l_ends {
        let j = index_of(key);
        let s = if owners[src - 1].is_empty() { vec![vec![]] } else { owners[src - 1].clone() };
        owners[j].extend(s);
    }
    for j in 0..n {
        if owners[j].is_empty() {
            owners[j] = vec![vec![]];
            continue;
        }
        let mut seen: BTreeSet<Vec<(usize, i64)>> = BTreeSet::new();
        let mut uniq = vec![];
        for r in &owners[j] {
            let mut k = r.clone();
            k.sort();
            if seen.insert(k) {
                uniq.push(r.clone());
            }
        }
        owners[j] = uniq;
    }
    owners
}

fn bridge_blocks(c: &[Vec<Z>], f: &BFace) -> Vec<Rows> {
    let mut blocks = vec![];
    for con in &f.cells {
        let mut r = vec![];
        for row in con {
            if all_zero(&row[0..5]) {
                continue;
            }
            let mut key = row[0..6].to_vec();
            key.push(if row.len() >= 8 { row[7].clone() } else { zi(0) });
            r.push((c.iter().position(|x| *x == key).unwrap() + 1, side_i(&row[6])));
        }
        blocks.push(r);
    }
    blocks
}

fn relevant_here(c: &[Vec<Z>], owners: &[Rows], rows: &Rows) -> R<bool> {
    let mut side: HashMap<usize, i64> = HashMap::new();
    for r in rows {
        side.entry(r.0).or_insert(r.1);
    }
    for r in owners {
        if r.is_empty() {
            return Ok(true);
        }
        let mut clash = false;
        let mut rest = vec![];
        for t in r {
            match side.get(&t.0) {
                None => rest.push(*t),
                Some(s) => {
                    if *s != t.1 {
                        clash = true;
                        break;
                    }
                }
            }
        }
        if clash {
            continue;
        }
        if rest.is_empty() {
            return Ok(true);
        }
        let mut rr = rows.clone();
        rr.extend(rest);
        if live_cell(c, &rr)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// L = [p, q, r] with c = k * (p x + q y + r)^2 exactly, or None.
fn square_root(c: &[Z]) -> Option<Vec<Z>> {
    let (a, b, cc) = (&c[0], &c[1], &c[2]);
    if b * b != zi(4) * a * cc {
        return None;
    }
    let g = a.gcd(cc).gcd(b);
    if g.is_zero() {
        return None;
    }
    let (mut aa, mut bb, mut cq) = (a / &g, b / &g, cc / &g);
    let mut sg = 1i64;
    if aa.is_negative() || cq.is_negative() {
        aa = -aa;
        bb = -bb;
        cq = -cq;
        sg = -1;
    }
    let p = aa.sqrt();
    let mut q = cq.sqrt();
    if &p * &p != aa || &q * &q != cq {
        return None;
    }
    if bb.is_negative() {
        q = -q;
    }
    if zi(2) * &p * &q != bb {
        return None;
    }
    let k = zi(sg) * &g;
    let r = if !p.is_zero() { Q::new(c[3].clone(), zi(2) * &k * &p) } else { Q::new(c[4].clone(), zi(2) * &k * &q) };
    let (lc, kk) = if !r.denom().is_one() {
        let r2 = &r * qi(2);
        if !r2.denom().is_one() {
            return None;
        }
        (vec![zi(2) * &p, zi(2) * &q, r2.to_integer()], Q::new(k.clone(), zi(4)))
    } else {
        (vec![p.clone(), q.clone(), r.to_integer()], qz(&k))
    };
    let want = [
        &kk * qz(&(&lc[0] * &lc[0])),
        &kk * qz(&(zi(2) * &lc[0] * &lc[1])),
        &kk * qz(&(&lc[1] * &lc[1])),
        &kk * qz(&(zi(2) * &lc[0] * &lc[2])),
        &kk * qz(&(zi(2) * &lc[1] * &lc[2])),
        &kk * qz(&(&lc[2] * &lc[2])),
    ];
    if want.iter().zip(c.iter()).all(|(w, x)| *w == qz(x)) {
        Some(lc)
    } else {
        None
    }
}

fn rows_to_con(c: &[Vec<Z>], rows: &Rows) -> Vec<Vec<Z>> {
    rows.iter()
        .map(|&(j, s)| {
            let mut r = c[j - 1][0..6].to_vec();
            r.push(zi(s));
            r.push(c[j - 1][6].clone());
            r
        })
        .collect()
}

fn holds_all(c: &[Vec<Z>], rows: &Rows, p: [f64; 2]) -> R<bool> {
    for &(j, s) in rows {
        if (s as f64) * eval_curve(c, j, p)? <= 0.0 {
            return Ok(false);
        }
    }
    Ok(true)
}

fn interior_point(c: &[Vec<Z>], rows: &Rows) -> R<([f64; 2], bool)> {
    let h = straight_h(c, rows);
    let r = 1e3;
    let mut cand = if h.is_empty() {
        [0.0, 0.0]
    } else {
        let (best, cand) = chebyshev(&h, r);
        if best <= 0.0 {
            return Ok(([0.0, 0.0], false));
        }
        cand
    };
    for pas in [1, 2] {
        for g in [8.0, 16.0, 64.0, 256.0, 1024.0] {
            let p = [round_half_away(cand[0] * g) as f64 / g, round_half_away(cand[1] * g) as f64 / g];
            if holds_all(c, rows, p)? {
                return Ok((p, true));
            }
        }
        if holds_all(c, rows, cand)? {
            return Ok((cand, true));
        }
        if pas == 1 {
            let (cd, found, _) = cell_sample_point(c, rows, &h, r, false)?;
            if !found {
                return Ok(([0.0, 0.0], false));
            }
            cand = cd;
        }
    }
    Ok(([0.0, 0.0], false))
}

// ---- the overlay of several pairs: a MIN of maxes ------------------------------------------------
fn key_rows(keys: &mut Vec<Vec<Z>>, con: &[Vec<Z>]) -> Rows {
    con.iter()
        .map(|r| {
            let mut k = r[0..6].to_vec();
            k.push(r[7].clone());
            (add_key(keys, k), side_i(&r[6]))
        })
        .collect()
}

fn merge_rows(a: &Rows, b: &Rows) -> Option<Rows> {
    let mut side: HashMap<usize, i64> = HashMap::new();
    for r in a {
        side.entry(r.0).or_insert(r.1);
    }
    let mut out = a.clone();
    for t in b {
        match side.get(&t.0) {
            None => out.push(*t),
            Some(s) => {
                if *s != t.1 {
                    return None;
                }
            }
        }
    }
    Some(out)
}

fn unique_stable(rows: Vec<Vec<Z>>) -> Vec<Vec<Z>> {
    let mut seen: BTreeSet<Vec<Z>> = BTreeSet::new();
    let mut out = vec![];
    for r in rows {
        if seen.insert(r.clone()) {
            out.push(r);
        }
    }
    out
}

struct Samples {
    x: Vec<[f64; 2]>,
    z: HashMap<(*const PolyN, i64), HashMap<usize, f64>>,
    keep: Vec<Arc<PolyN>>,
}

impl Samples {
    fn branch(&mut self, fp: &Arc<PolyN>, fb: i64, idx: &[usize]) -> R<Vec<f64>> {
        let key = (Arc::as_ptr(fp), fb);
        if !self.keep.iter().any(|k| Arc::ptr_eq(k, fp)) {
            self.keep.push(fp.clone());
        }
        let e = self.z.entry(key).or_default();
        let need: Vec<usize> = idx.iter().cloned().filter(|i| !e.contains_key(i)).collect();
        if !need.is_empty() {
            let pts: Vec<[f64; 2]> = need.iter().map(|&i| self.x[i]).collect();
            let v = alg_branch(fp, &pts, fb)?;
            for (i, vi) in need.iter().zip(v) {
                e.insert(*i, vi);
            }
        }
        Ok(idx.iter().map(|i| e[i]).collect())
    }
}

fn con_val(r: &[Z], x: [f64; 2]) -> R<f64> {
    if !r[7].is_zero() {
        return pair_line_eval(&r[0..6], side_i(&r[7]), x[0], x[1]);
    }
    Ok(eval_conic(&r[0..6], x[0], x[1]))
}

fn witnessed(con: &[Vec<Z>], alg: &[AlgRow], ii: &[IiRow], pp: &mut Samples) -> R<bool> {
    if pp.x.is_empty() {
        return Ok(false);
    }
    let mut idx: Vec<usize> = (0..pp.x.len()).collect();
    for r in con {
        let lim = 1e-9 * 1f64.max(to_f64(&r[0..6].iter().map(|x| x.abs()).max().unwrap()));
        let s = to_f64(&r[6]);
        let mut keep = vec![];
        for &i in &idx {
            if s * con_val(r, pp.x[i])? > lim {
                keep.push(i);
            }
        }
        idx = keep;
        if idx.is_empty() {
            return Ok(false);
        }
    }
    let k = 96;
    if idx.len() > k {
        let pick: Vec<usize> = linspace(1.0, idx.len() as f64, k).iter().map(|t| round_half_away(*t) as usize).collect();
        idx = pick.iter().map(|p| idx[p - 1]).collect::<BTreeSet<_>>().into_iter().collect();
    }
    for a in alg {
        let z = pp.branch(&a.fp, a.fb, &idx)?;
        idx = idx.iter().zip(z.iter()).filter(|(i, zi)| {
            let x = pp.x[**i];
            let qi = eval_face(&a.num, &a.den, x[0], x[1]);
            zi.is_finite() && (a.side as f64) * (qi - **zi) > 1e-9
        }).map(|(i, _)| *i).collect();
        if idx.is_empty() {
            return Ok(false);
        }
    }
    for t in ii {
        let z1 = pp.branch(&t.fp1, t.fb1, &idx)?;
        let z2 = pp.branch(&t.fp2, t.fb2, &idx)?;
        idx = idx.iter().zip(z1.iter().zip(z2.iter())).filter(|(_, (a1, a2))| a1.is_finite() && a2.is_finite() && *a1 - *a2 > 1e-9).map(|(i, _)| *i).collect();
        if idx.is_empty() {
            return Ok(false);
        }
    }
    Ok(!idx.is_empty())
}

fn overlay_live(keys: &[Vec<Z>], rows: &Rows, con: &[Vec<Z>], alg: &[AlgRow], ii: &[IiRow], pp: &mut Samples) -> R<bool> {
    if witnessed(con, alg, ii, pp)? {
        return Ok(true);
    }
    let k0: Vec<bool> = con.iter().map(|r| all_zero(&r[0..5])).collect();
    if con.iter().zip(k0.iter()).any(|(r, z)| *z && (&r[6] * &r[5]).is_negative()) {
        return Ok(false);
    }
    let mut kk = keys.to_vec();
    let nz: Vec<Vec<Z>> = con.iter().zip(k0.iter()).filter(|(_, z)| !**z).map(|(c, _)| c.clone()).collect();
    let r = key_rows(&mut kk, &nz);
    let all: Rows = rows.iter().cloned().chain(r).collect::<BTreeSet<_>>().into_iter().collect();
    if !live_cell(&kk, &all)? {
        return Ok(false);
    }
    if alg.is_empty() {
        return Ok(true);
    }
    let plain: Vec<Vec<Z>> = con.iter().zip(k0.iter()).filter(|(c, z)| !**z && c[7].is_zero()).map(|(c, _)| {
        let mut r = c[0..7].to_vec();
        r.push(zi(0));
        r
    }).collect();
    Ok(!alg_empty(keys, rows, &plain, alg, &pp.x.clone())?)
}

fn tie_empty(keys: &[Vec<Z>], rows: &Rows, con: &[Vec<Z>], ii: &[IiRow], pp: &mut Samples) -> R<bool> {
    if ii.is_empty() {
        return Ok(false);
    }
    let k: Vec<Vec<Z>> = con.iter().filter(|r| !all_zero(&r[0..5])).cloned().collect();
    let mut idx: Vec<usize> = (0..pp.x.len()).collect();
    for r in &k {
        let lim = -1e-9 * 1f64.max(to_f64(&r[0..6].iter().map(|x| x.abs()).max().unwrap()));
        let s = to_f64(&r[6]);
        let mut keep = vec![];
        for &i in &idx {
            if s * con_val(r, pp.x[i])? >= lim {
                keep.push(i);
            }
        }
        idx = keep;
    }
    let mut reg: Option<Option<Region>> = None;
    for t in ii {
        if !idx.is_empty() {
            let z1 = pp.branch(&t.fp1, t.fb1, &idx)?;
            let z2 = pp.branch(&t.fp2, t.fb2, &idx)?;
            if z1.iter().zip(z2.iter()).any(|(a, b)| a.is_finite() && b.is_finite() && a > b) {
                continue;
            }
        }
        if reg.is_none() {
            reg = Some(region_box(keys, rows, &k)?);
        }
        let rg = match reg.as_ref().unwrap() {
            None => return Ok(false),
            Some(r) => r,
        };
        let lo = BelowFace { fp: t.fp1.clone(), fb: t.fb1 };
        let hi = vec![BelowFace { fp: t.fp2.clone(), fb: t.fb2 }];
        let tf = alg_face_below(&lo, &hi, rg, 5, &[], None).unwrap_or(false);
        if tf {
            return Ok(true);
        }
    }
    Ok(false)
}

fn le_rows(f: &SubFace, o: &SubFace, con: &mut Vec<Vec<Z>>, alg: &mut Vec<AlgRow>, ii: &mut Vec<IiRow>) -> R<()> {
    if f.kind == SubKind::Quad && o.kind == SubKind::Quad {
        if let Some((c, sgn)) = diff_conic(&o.num, &o.den, &f.num, &f.den)? {
            let mut r = c;
            r.push(zi(sgn));
            r.push(zi(0));
            con.push(r);
        }
    } else if f.kind == SubKind::Quad {
        alg.push(AlgRow { fp: o.fp.clone().unwrap(), fb: o.fb, num: f.num.clone(), den: f.den.clone(), side: -1 });
    } else if o.kind == SubKind::Quad {
        alg.push(AlgRow { fp: f.fp.clone().unwrap(), fb: f.fb, num: o.num.clone(), den: o.den.clone(), side: 1 });
    } else {
        let (ff, of) = (f.fp.clone().unwrap(), o.fp.clone().unwrap());
        if !((Arc::ptr_eq(&ff, &of) || *ff == *of) && f.fb == o.fb) {
            ii.push(IiRow { fp1: of, fb1: o.fb, fp2: ff, fb2: f.fb });
        }
    }
    Ok(())
}

fn min_overlay(subs: &[Vec<SubFace>], infos: &[SubInfo], budget: f64) -> R<Vec<SubFace>> {
    let np = subs.len();
    let mut keys: Vec<Vec<Z>> = vec![];
    let mut cell_rows: Vec<Vec<Rows>> = vec![];
    let mut covered: Vec<Vec<bool>> = vec![];
    let mut faces_of: Vec<Vec<Vec<usize>>> = vec![];
    for g in 0..np {
        let cc = &infos[g].cell_con;
        cell_rows.push(cc.iter().map(|c| key_rows(&mut keys, c)).collect());
        let fo: Vec<Vec<usize>> = (0..cc.len()).map(|k| infos[g].face_cell.iter().enumerate().filter(|(_, fc)| **fc == k + 1).map(|(i, _)| i).collect()).collect();
        covered.push(fo.iter().map(|x| !x.is_empty()).collect());
        faces_of.push(fo);
    }
    let mut tup: Vec<Vec<usize>> = (0..cell_rows[0].len()).map(|k| vec![k]).collect();
    let mut trow: Vec<Rows> = cell_rows[0].clone();
    for g in 1..np {
        // every current tuple is merged with every cell of envelope g independently; results in order
        let merged = crate::par::map_range(trow.len(), |a| -> R<Vec<(Vec<usize>, Rows)>> {
            let mut o = vec![];
            for k in 0..cell_rows[g].len() {
                let r = match merge_rows(&trow[a], &cell_rows[g][k]) {
                    None => continue,
                    Some(r) => r,
                };
                if !live_cell(&keys, &r)? {
                    continue;
                }
                let mut t = tup[a].clone();
                t.push(k);
                o.push((t, r));
            }
            Ok(o)
        })?;
        let mut nt = vec![];
        let mut nr = vec![];
        for (t, r) in merged.into_iter().flatten() {
            nt.push(t);
            nr.push(r);
        }
        tup = nt;
        trow = nr;
    }
    let mut bound = 0f64;
    for t in &tup {
        let gs: Vec<usize> = (0..np).filter(|&g| covered[g][t[g]]).collect();
        let mut pr = 1f64;
        for &g in &gs {
            pr *= faces_of[g][t[g]].len() as f64;
        }
        bound += pr * gs.len() as f64;
    }
    if bound > budget {
        bail!("envelopeSubdivisionAlg:overlaySize", "the overlay of {} pairs would visit up to {} pieces over {} tuples, past the budget of {}.", np, bound, tup.len(), budget);
    }
    // every tuple of cells is resolved on its own; the faces are concatenated in tuple order
    let per_tuple = crate::par::map_range(tup.len(), |ti| -> R<Vec<SubFace>> {
        let (t, rows) = (&tup[ti], &trow[ti]);
        let mut out = vec![];
        let gs: Vec<usize> = (0..np).filter(|&g| covered[g][t[g]]).collect();
        if gs.is_empty() {
            return Ok(out);
        }
        let mut base: Vec<Vec<Z>> = vec![];
        for g in 0..np {
            base.extend(infos[g].cell_con[t[g]].iter().cloned());
        }
        let (_, _, px) = cell_sample_point(&keys, rows, &straight_h(&keys, rows), 1e3, true)?;
        let mut pp = Samples { x: px, z: HashMap::new(), keep: vec![] };
        let mut combos: Vec<(Vec<usize>, Vec<usize>, Vec<Vec<Z>>, Vec<AlgRow>)> = vec![(vec![], vec![], base.clone(), vec![])];
        for &g in &gs {
            let mut nc = vec![];
            for cb in &combos {
                for &k in &faces_of[g][t[g]] {
                    let fk = &subs[g][k];
                    let mut con = cb.2.clone();
                    con.extend(fk.con.iter().map(|r| r[0..8].to_vec()));
                    let con = unique_stable(con);
                    let mut alg = cb.3.clone();
                    alg.extend(fk.alg.iter().cloned());
                    if !overlay_live(&keys, rows, &con, &alg, &[], &mut pp)? {
                        continue;
                    }
                    let mut pg = cb.0.clone();
                    pg.push(g);
                    let mut pk = cb.1.clone();
                    pk.push(k);
                    nc.push((pg, pk, con, alg));
                }
            }
            combos = nc;
        }
        for cb in &combos {
            let pick: Vec<&SubFace> = cb.0.iter().zip(cb.1.iter()).map(|(&g, &k)| &subs[g][k]).collect();
            for w in 0..gs.len() {
                let f = pick[w];
                let mut con = cb.2.clone();
                let mut alg = cb.3.clone();
                let mut ii = vec![];
                for i in 0..gs.len() {
                    if i == w {
                        continue;
                    }
                    le_rows(f, pick[i], &mut con, &mut alg, &mut ii)?;
                }
                let con = unique_stable(con);
                if !overlay_live(&keys, rows, &con, &alg, &ii, &mut pp)? {
                    continue;
                }
                if !witnessed(&con, &alg, &ii, &mut pp)? && tie_empty(&keys, rows, &con, &ii, &mut pp)? {
                    continue;
                }
                out.push(SubFace { kind: f.kind.clone(), num: f.num.clone(), den: f.den.clone(), fp: f.fp.clone(), fb: f.fb, con, alg, ii, tri: vec![], tf: None, tfs: vec![] });
            }
        }
        Ok(out)
    })?;
    let out: Vec<SubFace> = per_tuple.into_iter().flatten().collect();
    Ok(out)
}


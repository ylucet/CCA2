//! The assembled envelope of one atom pair: its data types and the modules that build, prune
//! and evaluate it (ports of envelope_cell_pairs.py, envelope_cell_regions.py,
//! envelope_assemble.py, envelope_prune.py, envelope_rational_corners.py, envelope_eval.py).

use crate::alg_branch::alg_branch;
use crate::bail;
use crate::bridge::{bridge_face_alg, bridge_member};
use crate::conic_meet::solve2f;
use crate::corner::corner_face_alg;
use crate::errors::R;
use crate::poly_n::PolyN;
use crate::qua_con::QuaCon;
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeSet;
use std::sync::Arc;

// ==========================================================================================
// types
// ==========================================================================================
#[derive(Clone, Debug)]
pub struct QuadC {
    pub num: Vec<Z>,
    pub den: Z,
    pub con: Vec<Vec<Z>>,
    pub pair: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum FaceKind {
    Algebraic,
    Rational,
    Corner,
}

#[derive(Clone, Debug)]
pub struct AlgPair {
    pub d: Z,
    pub n1: Vec<Vec<Z>>,
    pub m1: Vec<Z>,
    pub n2: Vec<Vec<Z>>,
    pub m2: Vec<Z>,
    pub dlm: Vec<Vec<Z>>,
    pub dl: Vec<Z>,
    pub kp: Z,
    pub p: Vec<Vec<Z>>,
    pub fcache: std::sync::OnceLock<(Arc<crate::bridge::AlgPairF>, Arc<crate::bridge::AlgPairF>)>,
}

#[derive(Clone, Debug)]
pub struct RatPair {
    pub u: Vec<Z>,
    pub w: Vec<Z>,
    pub nu: Z,
    pub rq: Z,
    pub rl: Z,
    pub rc: Z,
    pub rd: Z,
    pub b: Vec<Z>,
    pub c: Vec<Z>,
    pub scale: Z,
    pub p: Vec<Vec<Z>>,
}

#[derive(Clone, Debug)]
pub enum PairData {
    Alg(AlgPair),
    Rat(RatPair),
}

#[derive(Clone, Debug)]
pub struct Sigma {
    pub quad: [Z; 3],
    pub root: i64,
    pub facet: usize,
    pub s0n: Vec<Z>,
    pub w: Vec<Z>,
    pub na: Z,
}

#[derive(Clone, Debug)]
pub struct Src {
    pub n: Vec<Z>,
    pub d: Z,
    pub tie: Vec<Z>,
    pub side: i64,
}

#[derive(Clone, Debug)]
pub struct BFace {
    pub kind: FaceKind,
    pub fp: Arc<PolyN>,
    pub fb: i64,
    // algebraic / rational
    pub cells: Vec<Vec<Vec<Z>>>,
    pub pair: Option<Arc<PairData>>,
    pub wit: Option<[Z; 3]>,
    pub rfn: Vec<Z>,
    pub rfd: Z,
    pub rden: Vec<Z>,
    // corner
    pub s: [f64; 2],
    pub sigma: Option<Sigma>,
    pub mm: [[f64; 2]; 2],
    pub p: [f64; 2],
    pub src: Option<Src>,
    pub con: Vec<Vec<Z>>,
    pub pair_index: usize,
}

impl BFace {
    pub fn blank(kind: FaceKind) -> BFace {
        BFace {
            kind,
            fp: Arc::new(PolyN::new()),
            fb: 0,
            cells: vec![],
            pair: None,
            wit: None,
            rfn: vec![],
            rfd: Z::one(),
            rden: vec![],
            s: [0.0, 0.0],
            sigma: None,
            mm: [[0.0; 2]; 2],
            p: [0.0, 0.0],
            src: None,
            con: vec![],
            pair_index: 0,
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct Env {
    pub quad: Vec<QuadC>,
    pub bridge: Vec<BFace>,
}

fn h2(n: &[Z]) -> Vec<Vec<Z>> {
    vec![vec![n[4].clone(), n[5].clone()], vec![n[5].clone(), n[6].clone()]]
}

pub fn first_nonzero_sign(c: &[Z]) -> i64 {
    match c.iter().find(|x| !x.is_zero()) {
        None => 1,
        Some(v) => sign(v) as i64,
    }
}

// ==========================================================================================
// envelope_cell_pairs
// ==========================================================================================
#[derive(Clone, Debug)]
pub struct CellPair {
    pub poly: bool,
    pub ij: [usize; 2],
    pub p: Vec<Vec<Z>>,
    pub q: Vec<Z>,
    pub win: char,
}

fn cell_half_planes(g: &QuaCon, k: usize) -> R<(Vec<Vec<Z>>, bool)> {
    let mut rows = vec![];
    let mut thin = false;
    for r in &g.fc[k - 1] {
        let cc = &g.ecq[r[0] as usize - 1];
        if !all_zero(&cc[0..3]) {
            bail!("envelopeCellPairs:curvedCell", "cell {} has a genuinely CONIC edge (edge {}).", k, r[0]);
        }
        if r[1] == 0 {
            thin = true;
            continue;
        }
        rows.push(cc[3..6].iter().map(|x| x * zi(r[1])).collect());
    }
    Ok((rows, thin))
}

fn tie_conic(ga: &QuaCon, i: usize, gb: &QuaCon, j: usize) -> Vec<Z> {
    let d: Vec<Z> = ga.fn_[i - 1].iter().zip(gb.fn_[j - 1].iter()).map(|(a, b)| a * &gb.fd[j - 1] - b * &ga.fd[i - 1]).collect();
    face_conic(&d)
}

pub fn envelope_cell_pairs(ga: &QuaCon, gb: &QuaCon, d: &[Vec<Z>]) -> R<Vec<CellPair>> {
    let mut pairs = vec![];
    for i in 1..=ga.nf {
        let (ra, thin_a) = cell_half_planes(ga, i)?;
        if thin_a {
            continue;
        }
        for j in 1..=gb.nf {
            let (rb, thin_b) = cell_half_planes(gb, j)?;
            if thin_b {
                continue;
            }
            let mut p = ra.clone();
            p.extend(rb.iter().cloned());
            p.extend(d.iter().cloned());
            if !feasible2(&p, true) {
                continue;
            }
            let q = tie_conic(ga, i, gb, j);
            if all_zero(&q) {
                pairs.push(CellPair { poly: true, ij: [i, j], p, q, win: 'A' });
                continue;
            }
            let (sgn, strict) = quad_sign_on_poly(&q, &p)?;
            if sgn != 0 && strict {
                pairs.push(CellPair { poly: true, ij: [i, j], p, q, win: if sgn > 0 { 'A' } else { 'B' } });
            } else {
                pairs.push(CellPair { poly: false, ij: [i, j], p, q, win: ' ' });
            }
        }
    }
    Ok(pairs)
}

// ==========================================================================================
// envelope_cell_regions
// ==========================================================================================
#[derive(Clone, Debug)]
pub struct Smap {
    pub mn: Vec<Vec<Z>>,
    pub tn: Vec<Z>,
    pub d: Z,
}

#[derive(Clone, Debug)]
pub struct RegionFace {
    pub r: Vec<Vec<Z>>,
    pub fn_: Vec<Z>,
    pub fd: Z,
    pub s: Smap,
}

fn mk_map(mn: Vec<Vec<Z>>, tn: Vec<Z>, d: Z) -> Smap {
    if d.is_negative() {
        Smap { mn: mn.iter().map(|r| neg(r)).collect(), tn: neg(&tn), d: -d }
    } else {
        Smap { mn, tn, d }
    }
}

fn mk_face(r: Vec<Vec<Z>>, fnn: Vec<Z>, fd: Z, s: Smap) -> R<RegionFace> {
    let (fnn, fd) = if fd.is_negative() { (neg(&fnn), -fd) } else { (fnn, fd) };
    let (fnn, fd) = canon(&fnn, &fd)?;
    Ok(RegionFace { r, fn_: fnn, fd, s })
}

pub fn envelope_cell_regions(n: &[Vec<Z>], m: &[Z], q: &Z, dd: &Z, p: &[Vec<Z>]) -> R<Vec<RegionFace>> {
    if n[0][1] != n[1][0] {
        bail!("envelopeCellRegions:shape", "N must be a symmetric 2x2.");
    }
    if !dd.is_positive() {
        bail!("envelopeCellRegions:denominator", "dd must be positive; got {}.", dd);
    }
    let mut faces = vec![];
    let k = p.len();
    if k == 0 {
        return Ok(faces);
    }
    let det_n = &n[0][0] * &n[1][1] - &n[0][1] * &n[1][0];
    let kk = vec![vec![n[1][1].clone(), -&n[0][1]], vec![-&n[0][1], n[0][0].clone()]];
    if det_n.is_positive() {
        let mut r = vec![];
        for row in p {
            let a = [&row[0] * &kk[0][0] + &row[1] * &kk[1][0], &row[0] * &kk[0][1] + &row[1] * &kk[1][1]];
            r.push(vec![dd * &a[0], dd * &a[1], -(&a[0] * &m[0] + &a[1] * &m[1]) + &row[2] * &det_n]);
        }
        if feasible2(&r, true) {
            let km = mv(&kk, m);
            let dd2 = dd * dd;
            let mut fnn = vec![Z::zero(); 10];
            fnn[4] = zi(2) * &dd2 * &kk[0][0];
            fnn[5] = zi(2) * &dd2 * &kk[0][1];
            fnn[6] = zi(2) * &dd2 * &kk[1][1];
            fnn[7] = -zi(2) * dd * &km[0];
            fnn[8] = -zi(2) * dd * &km[1];
            fnn[9] = dot(m, &km) - zi(2) * q * &det_n;
            let s = mk_map(kk.iter().map(|row| row.iter().map(|x| dd * x).collect()).collect(), vec![-&km[0], -&km[1]], det_n.clone());
            faces.push(mk_face(r, fnn, zi(2) * &det_n * dd, s)?);
        }
    }
    for i in 0..k {
        let a = vec![p[i][0].clone(), p[i][1].clone()];
        let c = &p[i][2];
        let na = &a[0] * &a[0] + &a[1] * &a[1];
        if na.is_zero() {
            continue;
        }
        let w = vec![a[1].clone(), -&a[0]];
        let dnd = quad(&w, n, &w);
        if !dnd.is_positive() {
            continue;
        }
        let p0 = vec![-c * &a[0], -c * &a[1]];
        let alpha = vec![&na * dd * &w[0], &na * dd * &w[1]];
        let beta = -&na * dot(&w, m) - quad(&w, n, &p0);
        let mut r = vec![];
        let mut dead = false;
        for j in 0..k {
            if j == i {
                continue;
            }
            let aj = &p[j][0..2];
            let cj = &p[j][2];
            let ajw = dot(aj, &w);
            let bj = dot(aj, &p0) + cj * &na;
            if ajw.is_zero() {
                if bj.is_negative() {
                    dead = true;
                    break;
                }
                continue;
            }
            r.push(vec![&ajw * &alpha[0], &ajw * &alpha[1], &ajw * &beta + &bj * &dnd]);
        }
        if dead {
            continue;
        }
        let nav = mv(n, &a);
        let lhs_x = [&dnd * &na * dd * &a[0], &dnd * &na * dd * &a[1]];
        let lhs_c = -&dnd * &na * dot(&a, m) - &dnd * dot(&a, &mv(n, &p0));
        let wn = dot(&w, &nav);
        r.push(vec![-(&lhs_x[0] - &wn * &alpha[0]), -(&lhs_x[1] - &wn * &alpha[1]), -(&lhs_c - &wn * &beta)]);
        if !feasible2(&r, true) {
            continue;
        }
        let mut fnn = vec![Z::zero(); 10];
        fnn[4] = zi(2) * &alpha[0] * &alpha[0];
        fnn[5] = zi(2) * &alpha[0] * &alpha[1];
        fnn[6] = zi(2) * &alpha[1] * &alpha[1];
        fnn[7] = zi(2) * &na * &dnd * dd * &p0[0] + zi(2) * &beta * &alpha[0];
        fnn[8] = zi(2) * &na * &dnd * dd * &p0[1] + zi(2) * &beta * &alpha[1];
        fnn[9] = -zi(2) * &na * &dnd * dot(&p0, m) - &dnd * quad(&p0, n, &p0) - zi(2) * &na * &na * &dnd * q + &beta * &beta;
        let s = mk_map(
            vec![vec![&w[0] * &alpha[0], &w[0] * &alpha[1]], vec![&w[1] * &alpha[0], &w[1] * &alpha[1]]],
            vec![&p0[0] * &dnd + &beta * &w[0], &p0[1] * &dnd + &beta * &w[1]],
            &na * &dnd,
        );
        faces.push(mk_face(r, fnn, zi(2) * &na * &na * &dnd * dd, s)?);
    }
    for i in 0..k {
        for j in i + 1..k {
            let ai = &p[i][0..2];
            let aj = &p[j][0..2];
            let d2 = &ai[0] * &aj[1] - &ai[1] * &aj[0];
            if d2.is_zero() {
                continue;
            }
            let rhs = [-&p[i][2], -&p[j][2]];
            let mut v0 = vec![&aj[1] * &rhs[0] - &ai[1] * &rhs[1], -&aj[0] * &rhs[0] + &ai[0] * &rhs[1]];
            let mut vd = d2.clone();
            if vd.is_negative() {
                v0 = neg(&v0);
                vd = -vd;
            }
            if !p.iter().all(|r| !(&r[0] * &v0[0] + &r[1] * &v0[1] + &r[2] * &vd).is_negative()) {
                continue;
            }
            let adj = [[aj[1].clone(), -&aj[0]], [-&ai[1], ai[0].clone()]];
            let sg = zi(sign(&d2) as i64);
            let mc = [-&vd * &m[0] - (&n[0][0] * &v0[0] + &n[0][1] * &v0[1]), -&vd * &m[1] - (&n[1][0] * &v0[0] + &n[1][1] * &v0[1])];
            let mut r = vec![];
            for row in &adj {
                r.push(vec![-&sg * &row[0] * &vd * dd, -&sg * &row[1] * &vd * dd, -&sg * (&row[0] * &mc[0] + &row[1] * &mc[1])]);
            }
            if !feasible2(&r, true) {
                continue;
            }
            let mut fnn = vec![Z::zero(); 10];
            fnn[7] = zi(2) * &vd * dd * &v0[0];
            fnn[8] = zi(2) * &vd * dd * &v0[1];
            fnn[9] = -zi(2) * &vd * dot(&v0, m) - quad(&v0, n, &v0) - zi(2) * q * &vd * &vd;
            let s = mk_map(vec![vec![Z::zero(), Z::zero()], vec![Z::zero(), Z::zero()]], v0.clone(), vd.clone());
            faces.push(mk_face(r, fnn, zi(2) * &vd * &vd * dd, s)?);
        }
    }
    Ok(faces)
}

// ==========================================================================================
// envelope_assemble
// ==========================================================================================
fn rows_to_conics(r: &[Vec<Z>]) -> R<Vec<Vec<Z>>> {
    let mut out = vec![];
    for row in r {
        let mut c = conic_lin(row)?;
        c.push(zi(first_nonzero_sign(row)));
        c.push(zi(0));
        out.push(c);
    }
    Ok(out)
}

struct CC {
    num: Vec<Z>,
    den: Z,
    con: Vec<Vec<Z>>,
}

fn cells_for(n: &[Z], d: &Z, p: &[Vec<Z>]) -> R<Vec<CC>> {
    let nn = h2(n);
    let m = vec![n[7].clone(), n[8].clone()];
    let faces = envelope_cell_regions(&nn, &m, &n[9], d, p)?;
    let mut out = vec![];
    for f in faces {
        out.push(CC { num: f.fn_, den: f.fd, con: rows_to_conics(&f.r)? });
    }
    Ok(out)
}

fn pullback(q: &[Z], mn: &[Vec<Z>], tn: &[Z], dd: &Z) -> Vec<Z> {
    let l1 = [mn[0][0].clone(), mn[0][1].clone(), tn[0].clone()];
    let l2 = [mn[1][0].clone(), mn[1][1].clone(), tn[1].clone()];
    let mul_aff = |f: &[Z; 3], g: &[Z; 3]| -> Vec<Z> {
        vec![
            &f[0] * &g[0],
            &f[0] * &g[1] + &f[1] * &g[0],
            &f[1] * &g[1],
            &f[0] * &g[2] + &f[2] * &g[0],
            &f[1] * &g[2] + &f[2] * &g[1],
            &f[2] * &g[2],
        ]
    };
    let a = mul_aff(&l1, &l1);
    let b = mul_aff(&l1, &l2);
    let c = mul_aff(&l2, &l2);
    let mut out: Vec<Z> = (0..6).map(|i| &q[0] * &a[i] + &q[1] * &b[i] + &q[2] * &c[i]).collect();
    let lin: Vec<Z> = (0..3).map(|i| dd * (&q[3] * &l1[i] + &q[4] * &l2[i])).collect();
    out[3] += &lin[0];
    out[4] += &lin[1];
    out[5] += &lin[2] + &q[5] * dd * dd;
    out
}

fn side_cells(n: &[Z], d: &Z, p: &[Vec<Z>], tie: &[Z], side: i64) -> R<Vec<CC>> {
    let nn = h2(n);
    let m = vec![n[7].clone(), n[8].clone()];
    let faces = envelope_cell_regions(&nn, &m, &n[9], d, p)?;
    let mut cells = vec![];
    for f in faces {
        let c = pullback(tie, &f.s.mn, &f.s.tn, &f.s.d);
        let mut con = rows_to_conics(&f.r)?;
        if !all_zero(&c) {
            let mut row = conic(&c)?;
            row.push(zi(side * first_nonzero_sign(&c)));
            row.push(zi(0));
            con.push(row);
        }
        cells.push(CC { num: f.fn_, den: f.fd, con });
    }
    Ok(cells)
}

fn side_winner(na: &[Z], da: &Z, nb: &[Z], db: &Z, p: &[Vec<Z>]) -> R<i32> {
    let c = diff_conic(na, da, nb, db)?;
    let (s, _) = quad_sign_on_poly(&c, p)?;
    Ok(if s != 0 { s } else { 1 })
}

fn straight_tie_guard(na: &[Z], nb: &[Z], q: &[Z]) -> R<()> {
    let (a, _b, c, d, e, f) = (&q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
    let (disc, w) = if !a.is_zero() {
        (zi(4) * d * d - zi(16) * a * f, vec![-&q[1], zi(2) * a])
    } else {
        (zi(4) * e * e - zi(16) * c * f, vec![zi(1), zi(0)])
    };
    if !disc.is_positive() {
        return Ok(());
    }
    let r = disc.sqrt();
    if &r * &r == disc {
        return Ok(());
    }
    for n in [na, nb] {
        if !quad(&w, &h2(n), &w).is_zero() {
            bail!("envelopeAssemble:irrationalStraightFacet", "a cell pair ties on irrational parallel lines and a face CURVES along them.");
        }
    }
    Ok(())
}

pub fn envelope_assemble(ga: &QuaCon, gb: &QuaCon, d: &[Vec<Z>]) -> R<(Env, ())> {
    let pairs = envelope_cell_pairs(ga, gb, d)?;
    let mut env = Env::default();
    for (k0, pr) in pairs.iter().enumerate() {
        let k = k0 + 1;
        let (na, da) = (&ga.fn_[pr.ij[0] - 1], &ga.fd[pr.ij[0] - 1]);
        let (nb, db) = (&gb.fn_[pr.ij[1] - 1], &gb.fd[pr.ij[1] - 1]);
        let collect = |cells: Vec<CC>, env: &mut Env| {
            for c in cells {
                env.quad.push(QuadC { num: c.num, den: c.den, con: c.con, pair: vec![k] });
            }
        };
        if pr.poly {
            let cells = if pr.win == 'A' { cells_for(na, da, &pr.p)? } else { cells_for(nb, db, &pr.p)? };
            collect(cells, &mut env);
            continue;
        }
        let (sgn, strict) = quad_sign_on_poly(&pr.q, &pr.p)?;
        if sgn != 0 && !strict {
            let cells = if sgn > 0 { cells_for(na, da, &pr.p)? } else { cells_for(nb, db, &pr.p)? };
            collect(cells, &mut env);
            continue;
        }
        if all_zero(&pr.q[0..3]) {
            let mut pa = pr.p.clone();
            pa.push(pr.q[3..6].to_vec());
            collect(cells_for(na, da, &pa)?, &mut env);
            let mut pb = pr.p.clone();
            pb.push(neg(&pr.q[3..6]));
            collect(cells_for(nb, db, &pb)?, &mut env);
            continue;
        }
        let corners = |env: &mut Env| -> R<()> {
            let mut cf = corner_face_alg(na, da, &pr.p, &pr.q, 1)?;
            cf.extend(corner_face_alg(nb, db, &pr.p, &pr.q, -1)?);
            for mut c in cf {
                c.pair_index = k;
                env.bridge.push(c);
            }
            Ok(())
        };
        match bridge_face_alg(na, da, nb, db, &pr.p) {
            Err(me) => {
                if me.id == "bridgeFaceAlg:tangentialArc" {
                    let s2 = side_winner(na, da, nb, db, &pr.p)?;
                    let cells = if s2 > 0 { cells_for(na, da, &pr.p)? } else { cells_for(nb, db, &pr.p)? };
                    collect(cells, &mut env);
                    continue;
                }
                if me.id == "bridgeFaceAlg:straightTie" || me.id == "bridgeFaceRat:noSeparateFace" {
                    collect(side_cells(na, da, &pr.p, &pr.q, 1)?, &mut env);
                    collect(side_cells(nb, db, &pr.p, &pr.q, -1)?, &mut env);
                    straight_tie_guard(na, nb, &pr.q)?;
                    corners(&mut env)?;
                    continue;
                }
                return Err(me);
            }
            Ok(mut face) => {
                face.pair_index = k;
                env.bridge.push(face);
                collect(side_cells(na, da, &pr.p, &pr.q, 1)?, &mut env);
                collect(side_cells(nb, db, &pr.p, &pr.q, -1)?, &mut env);
                corners(&mut env)?;
            }
        }
    }
    Ok((env, ()))
}

// ==========================================================================================
// envelope_prune
// ==========================================================================================
fn union(a: &[usize], b: &[usize]) -> Vec<usize> {
    a.iter().chain(b.iter()).cloned().collect::<BTreeSet<_>>().into_iter().collect()
}

fn linear_rows(con: &[Vec<Z>]) -> Vec<Vec<Z>> {
    let mut l = vec![];
    for r in con {
        if !all_zero(&r[0..3]) {
            continue;
        }
        if r.len() >= 8 && !r[7].is_zero() {
            continue;
        }
        l.push(r[3..6].iter().map(|x| &r[6] * x).collect());
    }
    l
}

fn contains_region(li: &[Vec<Z>], conj: &[Vec<Z>]) -> bool {
    for r in conj {
        if !all_zero(&r[0..3]) {
            return false;
        }
        if r.len() >= 8 && !r[7].is_zero() {
            return false;
        }
        let row: Vec<Z> = r[3..6].iter().map(|x| &r[6] * x).collect();
        let mut rows = li.to_vec();
        rows.push(neg(&row));
        if feasible2(&rows, true) {
            return false;
        }
    }
    true
}

fn at_least_as_large(fj: &QuadC, fi: &QuadC, li: &[Vec<Z>]) -> R<bool> {
    let num: Vec<Z> = fj.num.iter().zip(fi.num.iter()).map(|(a, b)| a * &fi.den - b * &fj.den).collect();
    if !all_zero(&num[0..4]) {
        return Ok(false);
    }
    let c = vec![num[4].clone(), zi(2) * &num[5], num[6].clone(), zi(2) * &num[7], zi(2) * &num[8], zi(2) * &num[9]];
    if all_zero(&c) {
        return Ok(true);
    }
    let nz = c.iter().find(|x| !x.is_zero()).unwrap();
    let (s, _) = quad_sign_on_poly(&conic(&c)?, li)?;
    Ok(s * sign(nz) >= 1)
}

pub fn envelope_prune(env: &Env) -> R<Env> {
    let mut quad: Vec<QuadC> = env.quad.clone();
    let n0 = quad.len();
    let mut keep = vec![true; n0];
    for i in 0..n0 {
        if !keep[i] {
            continue;
        }
        for j in i + 1..n0 {
            if !keep[j] {
                continue;
            }
            if quad[i].num == quad[j].num && quad[i].den == quad[j].den && quad[i].con == quad[j].con {
                quad[i].pair = union(&quad[i].pair, &quad[j].pair);
                keep[j] = false;
            }
        }
    }
    let idx: Vec<usize> = (0..n0).filter(|&i| keep[i]).collect();
    for &i in &idx {
        if !keep[i] {
            continue;
        }
        let li = linear_rows(&quad[i].con);
        if li.is_empty() {
            continue;
        }
        for &j in &idx {
            if j == i || !keep[j] || !keep[i] {
                continue;
            }
            if !contains_region(&li, &quad[j].con) {
                continue;
            }
            if !at_least_as_large(&quad[j], &quad[i], &li)? {
                continue;
            }
            quad[j].pair = union(&quad[j].pair, &quad[i].pair);
            keep[i] = false;
        }
    }
    Ok(Env { quad: quad.into_iter().zip(keep).filter(|(_, k)| *k).map(|(q, _)| q).collect(), bridge: env.bridge.clone() })
}

// ==========================================================================================
// envelope_rational_corners
// ==========================================================================================
fn rational_corner(sg: &Sigma) -> Option<(Vec<Z>, Z)> {
    let (al, be, ga) = (&sg.quad[0], &sg.quad[1], &sg.quad[2]);
    let (p, q);
    if sg.root == 0 || al.is_zero() {
        if be.is_zero() {
            return None;
        }
        p = -ga;
        q = be.clone();
    } else {
        let d = be * be - zi(4) * al * ga;
        if d.is_negative() {
            return None;
        }
        let r = d.sqrt();
        if &r * &r != d {
            return None;
        }
        p = -be + zi(sg.root) * r;
        q = zi(2) * al;
    }
    let mut sn = vec![&sg.s0n[0] * &q + &p * &sg.w[0], &sg.s0n[1] * &q + &p * &sg.w[1]];
    let mut sd = &sg.na * &q;
    if sd.is_negative() {
        sn = neg(&sn);
        sd = -sd;
    }
    let gs = sn[0].gcd(&sn[1]).gcd(&sd);
    Some((vec![sn[0].div_floor(&gs), sn[1].div_floor(&gs)], sd.div_floor(&gs)))
}

pub fn envelope_rational_corners(env: &Env) -> R<Env> {
    let mut out = Env { quad: env.quad.clone(), bridge: vec![] };
    for c in &env.bridge {
        if c.kind != FaceKind::Corner || c.con.is_empty() || c.sigma.is_none() {
            out.bridge.push(c.clone());
            continue;
        }
        let (sn, sd) = match rational_corner(c.sigma.as_ref().unwrap()) {
            None => {
                out.bridge.push(c.clone());
                continue;
            }
            Some(x) => x,
        };
        let src = c.src.as_ref().unwrap();
        let n = &src.n;
        let d = &src.d;
        let nn = h2(n);
        let m = vec![n[7].clone(), n[8].clone()];
        let g0 = quad(&sn, &nn, &sn) + zi(2) * &sd * dot(&m, &sn) + zi(2) * &sd * &sd * &n[9];
        let mut num = vec![Z::zero(); 7];
        num.extend([zi(2) * d * &sd * &sn[0], zi(2) * d * &sd * &sn[1], -g0]);
        let den = zi(2) * d * &sd * &sd;
        let gg = den.gcd(&num[7].gcd(&num[8])).gcd(&num[9]);
        out.quad.push(QuadC { num: num.iter().map(|x| x.div_floor(&gg)).collect(), den: den.div_floor(&gg), con: c.con.clone(), pair: vec![c.pair_index] });
    }
    Ok(out)
}

// ==========================================================================================
// envelope_eval
// ==========================================================================================
const SQEPS: f64 = 1.4901161193847656e-08;

pub fn envelope_eval(env: &Env, x: &[[f64; 2]]) -> R<Vec<f64>> {
    let n = x.len();
    let mut val = vec![f64::INFINITY; n];
    if env.quad.is_empty() && env.bridge.is_empty() {
        return Ok(val);
    }
    let mut pairs: BTreeSet<usize> = BTreeSet::new();
    for q in &env.quad {
        pairs.extend(q.pair.iter().cloned());
    }
    for f in &env.bridge {
        pairs.insert(f.pair_index);
    }
    let r2 = 2f64.sqrt();
    let dirs = [[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0], [1.0 / r2, 1.0 / r2], [1.0 / r2, -1.0 / r2], [-1.0 / r2, 1.0 / r2], [-1.0 / r2, -1.0 / r2]];
    for i in 0..n {
        let xi = x[i];
        let mut best = f64::NEG_INFINITY;
        let mut diverges = false;
        for &p in &pairs {
            let (mut pv, _, _, mut cover) = pair_value(env, p, xi)?;
            if !cover {
                let h = 1e-7 * 1f64.max(xi[0].abs()).max(xi[1].abs());
                for d in &dirs {
                    let (_, nk, ni, nc) = pair_value(env, p, [xi[0] + h * d[0], xi[1] + h * d[1]])?;
                    if nc {
                        let v = face_value_at(env, nk, ni, xi)?;
                        if v.is_finite() {
                            pv = v;
                            cover = true;
                        }
                        break;
                    }
                }
            }
            if !cover {
                diverges = true;
                break;
            }
            if pv > best {
                best = pv;
            }
        }
        if !diverges && best.is_finite() {
            val[i] = best;
        }
    }
    Ok(val)
}

/// (value, kind 0 none / 1 quad / 2 bridge, 1-based index, covered)
fn pair_value(env: &Env, p: usize, x: [f64; 2]) -> R<(f64, u8, usize, bool)> {
    let mut cover = false;
    let mut pv = f64::NEG_INFINITY;
    let mut pk = 0u8;
    let mut pi = 0usize;
    for (k, f) in env.quad.iter().enumerate() {
        if !f.pair.contains(&p) {
            continue;
        }
        if !in_region(&f.con, x)? {
            continue;
        }
        let v = eval_face(&f.num, &f.den, x[0], x[1]);
        if v > pv {
            pv = v;
            pk = 1;
            pi = k + 1;
        }
        cover = true;
    }
    for (k, face) in env.bridge.iter().enumerate() {
        if face.pair_index != p {
            continue;
        }
        if !face_owns(face, x)? {
            continue;
        }
        let v = bridge_value(face, x)?;
        if !v.is_finite() {
            continue;
        }
        cover = true;
        if v > pv {
            pv = v;
            pk = 2;
            pi = k + 1;
        }
    }
    Ok((pv, pk, pi, cover))
}

fn face_value_at(env: &Env, kind: u8, k: usize, x: [f64; 2]) -> R<f64> {
    if kind == 1 {
        let f = &env.quad[k - 1];
        return Ok(eval_face(&f.num, &f.den, x[0], x[1]));
    }
    bridge_value(&env.bridge[k - 1], x)
}

fn face_owns(face: &BFace, x: [f64; 2]) -> R<bool> {
    if face.kind == FaceKind::Corner {
        let mu = solve2f(&face.mm, [face.p[0] - x[0], face.p[1] - x[1]]);
        if !(mu[0].is_finite() && mu[1].is_finite()) {
            return Ok(false);
        }
        let lim = -1e-9 * 1f64.max(mu[0].abs()).max(mu[1].abs());
        return Ok(mu[0] >= lim && mu[1] >= lim);
    }
    Ok(bridge_member(face.pair.as_ref().unwrap(), &[x])?[0])
}

pub fn bridge_value(face: &BFace, x: [f64; 2]) -> R<f64> {
    match face.kind {
        FaceKind::Corner => {
            let src = face.src.as_ref().unwrap();
            let n: Vec<f64> = src.n.iter().map(to_f64).collect();
            let s = face.s;
            let g = ((s[0] * (n[4] * s[0] + n[5] * s[1]) + s[1] * (n[5] * s[0] + n[6] * s[1])) / 2.0 + n[7] * s[0] + n[8] * s[1] + n[9]) / to_f64(&src.d);
            Ok(x[0] * s[0] + x[1] * s[1] - g)
        }
        FaceKind::Algebraic => Ok(alg_branch(&face.fp, &[x], face.fb)?[0]),
        FaceKind::Rational => {
            let q: Vec<f64> = face.rfn.iter().map(to_f64).collect();
            let num = q[4] / 2.0 * x[0].powi(2) + q[5] * x[0] * x[1] + q[6] / 2.0 * x[1].powi(2) + q[7] * x[0] + q[8] * x[1] + q[9];
            let dn: Vec<f64> = face.rden.iter().map(to_f64).collect();
            let den = dn[0] * x[0] + dn[1] * x[1] + dn[2];
            if den == 0.0 {
                return Ok(f64::NAN);
            }
            Ok(num / den / to_f64(&face.rfd))
        }
    }
}

pub fn in_region(con: &[Vec<Z>], x: [f64; 2]) -> R<bool> {
    let xmax = 1f64.max(x[0].abs()).max(x[1].abs());
    for r in con {
        let c = &r[0..6];
        let (q, scale) = if r.len() >= 8 && !r[7].is_zero() {
            (pair_line_eval(c, crate::cells::side_i(&r[7]), x[0], x[1])?, xmax)
        } else {
            let cf: Vec<f64> = c.iter().map(to_f64).collect();
            let cm = cf.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            (eval_conic_f(&cf, x[0], x[1]), 1f64.max(cm * xmax * xmax))
        };
        if !(to_f64(&r[6]) * q >= -SQEPS * scale) {
            return Ok(false);
        }
    }
    Ok(true)
}

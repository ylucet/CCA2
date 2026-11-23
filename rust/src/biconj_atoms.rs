//! biconj_atoms -- the envelope of a QuaPol as a MIN over atom pairs (port of biconj_atoms.py),
//! and biconj_q's atoms route.

use crate::alg_con::{assemble_alg_crv, AlgCon, EnvTri, EnvelopeMin};
use crate::bail;
use crate::biconj_q::{hull_rows, is_declined, subdivide_pair, Opts};
use crate::conic_meet::solve2f;
use crate::conj_q::conj_q;
use crate::envelope::{envelope_assemble, envelope_prune, Env};
use crate::envelope_subdivision_alg::{envelope_subdivision_alg, AlgOpts};
use crate::errors::R;
use crate::numeric::{convf, polyval_f, rcond2, roots, EPS};
use crate::qua_con::QuaCon;
use crate::qua_pol::QuaPol;
use crate::rat_q::*;
use crate::triangle_face::{cut_triangles, triangle_face, AtomTri};
use num_traits::{One, Zero};
use std::sync::OnceLock;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

// ==========================================================================================
// the atoms route of biconjQ
// ==========================================================================================
pub fn atoms_route(obj: &QuaPol, opts: &Opts) -> R<AlgCon> {
    let (pairs, tris) = biconj_atoms(obj)?;
    let aopts = AlgOpts { overlay_budget: opts.overlay_budget.unwrap_or(3000.0), ..Default::default() };
    let fallback = || AlgCon::EnvelopeMin(EnvelopeMin { pairs: pairs.clone(), tris: tris.iter().map(env_tri).collect() });
    if pairs.len() == 1 && tris.is_empty() {
        let p = envelope_prune(&pairs[0])?;
        return match subdivide_pair(p, &hull_rows(obj)?) {
            Ok(h) => Ok(h),
            Err(me) => {
                if !is_declined(&me) {
                    return Err(me);
                }
                Ok(fallback())
            }
        };
    }
    let run = || -> R<AlgCon> {
        let envs = crate::par::map(&pairs, envelope_prune)?;
        if tris.is_empty() {
            let sub = envelope_subdivision_alg(&envs, &hull_rows(obj)?, &[], &aopts)?;
            return Ok(AlgCon::AlgCrv(assemble_alg_crv(&sub)?));
        }
        let tfs = crate::par::map(&tris, |t| Ok(Arc::new(triangle_face(t)?)))?;
        let sub = envelope_subdivision_alg(&envs, &hull_rows(obj)?, &[], &aopts)?;
        Ok(AlgCon::AlgCrv(assemble_alg_crv(&cut_triangles(sub, tfs)?)?))
    };
    match run() {
        Ok(h) => Ok(h),
        Err(me) => {
            if !is_declined(&me) {
                return Err(me);
            }
            Ok(fallback())
        }
    }
}

fn env_tri(t: &AtomTri) -> EnvTri {
    EnvTri { t: t.t, s: t.s, c: t.c, g: t.g.clone() }
}

// ==========================================================================================
// types
// ==========================================================================================
#[derive(Clone, Debug)]
enum PieceRegion {
    Poly(Vec<[f64; 2]>),
    Rows(Vec<[f64; 3]>),
}

impl PartialEq for PieceRegion {
    fn eq(&self, o: &Self) -> bool {
        match (self, o) {
            (PieceRegion::Poly(a), PieceRegion::Poly(b)) => a == b,
            (PieceRegion::Rows(a), PieceRegion::Rows(b)) => a == b,
            _ => false,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
enum AKind {
    Piece,
    Edge,
    Ray,
    Vertex,
}

struct Atom {
    kind: AKind,
    q: QuaPol,
    region: Option<PieceRegion>,
    u: [f64; 2],
    w: [f64; 2],
    dr: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
enum FKind {
    Piece,
    Edge,
    Vertex,
}

struct Form {
    kind: FKind,
    a: [[f64; 2]; 2],
    b: [f64; 2],
    c: f64,
    u: [f64; 2],
    d: [f64; 2],
    al: f64,
    be: f64,
    h: [[f64; 2]; 2],
    l: [f64; 2],
    poly: Option<PieceRegion>,
    tmax: f64,
    gn: Vec<Z>,
    gd: Z,
    dom: OnceLock<Vec<[f64; 3]>>,
}

fn norm(v: [f64; 2]) -> f64 {
    v[0].hypot(v[1])
}

// ==========================================================================================
// biconjAtoms
// ==========================================================================================
pub fn biconj_atoms(f: &QuaPol) -> R<(Vec<Env>, Vec<AtomTri>)> {
    let mut d: Vec<Vec<Z>> = vec![];
    if f.e.iter().any(|e| e[2] == 0) {
        let (_, dd) = conj_q(f)?;
        d = dd.into_iter().collect::<BTreeSet<_>>().into_iter().collect();
        if d.iter().any(|r| r[0].is_zero() && r[1].is_zero() && sign(&r[2]) < 0) {
            bail!("biconjAtoms:emptyDomain", "dom f* is empty: f has no affine minorant, so co f is -infinity everywhere.");
        }
    }
    let df: Vec<[f64; 3]> = d.iter().map(|r| [to_f64(&r[0]), to_f64(&r[1]), to_f64(&r[2])]).collect();
    let (atoms, forms) = collect_atoms(f)?;
    if atoms.len() == 1 && !d.is_empty() {
        let (g1, _) = conj_q(&atoms[0].q)?;
        let tris = triple_faces(&forms, &df);
        return Ok((vec![envelope_assemble(&g1, &g1, &d)?.0], tris));
    }
    if atoms.len() < 2 {
        bail!("biconjAtoms:oneAtom", "the input reduces to {} atom(s).", atoms.len());
    }
    let tris = triple_faces(&forms, &df);
    let (ab, tris) = active_atom_pairs(&atoms, &forms, tris, &df)?;
    // the atoms' conjugates are independent, and so is the envelope of each pair of them
    let ids: Vec<usize> = ab.iter().flatten().cloned().collect::<BTreeSet<_>>().into_iter().collect();
    let gs = crate::par::map(&ids, |&a| Ok(conj_q(&atoms[a - 1].q)?.0))?;
    let g: BTreeMap<usize, QuaCon> = ids.into_iter().zip(gs).collect();
    let pairs = crate::par::map(&ab, |r| Ok(envelope_assemble(&g[&r[0]], &g[&r[1]], &d)?.0))?;
    Ok((pairs, tris))
}

// ------------------------------------------------------------------------------------------
fn active_atom_pairs(atoms: &[Atom], forms: &[Form], tris: Vec<AtomTri>, d: &[[f64; 3]]) -> R<(Vec<[usize; 2]>, Vec<AtomTri>)> {
    let na = atoms.len();
    let holds: Vec<Vec<usize>> = forms.iter().map(|fm| (0..na).filter(|&a| atom_holds(&atoms[a], fm)).map(|a| a + 1).collect()).collect();
    let mut req: Vec<(Vec<[usize; 2]>, Vec<usize>)> = vec![];
    let nf = forms.len();
    for i in 0..nf {
        for j in i + 1..nf {
            if sweeps_no_area(&forms[i], &forms[j]) {
                continue;
            }
            if !tie_arc(i, j, forms, d) {
                continue;
            }
            req.push(requirement(&[holds[i].clone(), holds[j].clone()], na));
        }
    }
    for i in 0..nf {
        if forms[i].kind == FKind::Piece || (!d.is_empty() && forms[i].kind == FKind::Edge) {
            req.push((vec![], holds[i].clone()));
        }
    }
    let mut keep_t = vec![true; tris.len()];
    for (t, tr) in tris.iter().enumerate() {
        if !triangle_active(tr, forms, d) {
            keep_t[t] = false;
            continue;
        }
        let h: Vec<Vec<usize>> = tr.forms.iter().map(|c| holds[c - 1].clone()).collect();
        let r = requirement(&h, na);
        if !r.0.is_empty() || !r.1.is_empty() {
            keep_t[t] = false;
            req.push(r);
        }
    }
    let tris: Vec<AtomTri> = tris.into_iter().zip(keep_t).filter(|(_, k)| *k).map(|(t, _)| t).collect();
    if req.is_empty() {
        let mut ab = vec![];
        for a in 1..=na {
            for b in a + 1..=na {
                ab.push([a, b]);
            }
        }
        return Ok((ab, tris));
    }
    Ok((cover_requirements(&req, na)?, tris))
}

fn requirement(h: &[Vec<usize>], na: usize) -> (Vec<[usize; 2]>, Vec<usize>) {
    let mut all1: BTreeSet<usize> = h[0].iter().cloned().collect();
    for x in &h[1..] {
        let s: BTreeSet<usize> = x.iter().cloned().collect();
        all1 = all1.intersection(&s).cloned().collect();
    }
    let mut p = vec![];
    for a in 1..=na {
        for b in a + 1..=na {
            if h.iter().all(|x| x.contains(&a) || x.contains(&b)) {
                p.push([a, b]);
            }
        }
    }
    (p, all1.into_iter().collect())
}

fn cover_requirements(req: &[(Vec<[usize; 2]>, Vec<usize>)], na: usize) -> R<Vec<[usize; 2]>> {
    let mut cand = vec![];
    for a in 1..=na {
        for b in a + 1..=na {
            cand.push([a, b]);
        }
    }
    let sat: Vec<Vec<bool>> = req.iter().map(|r| cand.iter().map(|c| r.0.contains(c) || r.1.contains(&c[0]) || r.1.contains(&c[1])).collect()).collect();
    if sat.iter().any(|row| !row.iter().any(|x| *x)) {
        bail!("biconjAtoms:uncovered", "a tie no atom pair holds -- the form holders are inconsistent.");
    }
    let mut ab = vec![];
    let mut open = vec![true; req.len()];
    while open.iter().any(|x| *x) {
        let score: Vec<usize> = (0..cand.len()).map(|c| (0..req.len()).filter(|&r| open[r] && sat[r][c]).count()).collect();
        let mx = *score.iter().max().unwrap();
        let c = score.iter().position(|s| *s == mx).unwrap();
        ab.push(cand[c]);
        for r in 0..req.len() {
            open[r] = open[r] && !sat[r][c];
        }
    }
    ab.sort();
    Ok(ab)
}

fn atom_holds(at: &Atom, fm: &Form) -> bool {
    let tol = 1e-12;
    let (pts, dirs): (Vec<[f64; 2]>, Vec<[f64; 2]>) = match fm.kind {
        FKind::Vertex => (vec![fm.u], vec![]),
        FKind::Edge => {
            if fm.tmax.is_infinite() {
                (vec![fm.u], vec![fm.d])
            } else {
                (vec![fm.u, [fm.u[0] + fm.d[0], fm.u[1] + fm.d[1]]], vec![])
            }
        }
        FKind::Piece => return at.kind == AKind::Piece && at.region == fm.poly,
    };
    match at.kind {
        AKind::Vertex => dirs.is_empty() && pts.iter().all(|p| norm([p[0] - at.u[0], p[1] - at.u[1]]) < tol),
        AKind::Edge => dirs.is_empty() && pts.iter().all(|p| on_segment(*p, at.u, at.w)),
        AKind::Ray => pts.iter().all(|p| on_ray(*p, at.u, at.dr)) && dirs.iter().all(|d| parallel_same(*d, at.dr)),
        AKind::Piece => {
            let p = at.region.as_ref().unwrap();
            let tf = pts.iter().all(|q| in_piece(*q, p));
            if tf && !dirs.is_empty() {
                let r = match p {
                    PieceRegion::Poly(_) => return false,
                    PieceRegion::Rows(r) => r,
                };
                let m = r.iter().flatten().fold(1f64, |a, x| a.max(x.abs()));
                for d in &dirs {
                    let nd = 1f64.max(norm(*d));
                    if !r.iter().all(|row| row[0] * d[0] + row[1] * d[1] >= -1e-9 * m * nd) {
                        return false;
                    }
                }
            }
            tf
        }
    }
}

fn on_segment(p: [f64; 2], u: [f64; 2], w: [f64; 2]) -> bool {
    let d = [w[0] - u[0], w[1] - u[1]];
    let l2 = d[0] * d[0] + d[1] * d[1];
    let t = ((p[0] - u[0]) * d[0] + (p[1] - u[1]) * d[1]) / l2;
    (-1e-12..=1.0 + 1e-12).contains(&t) && norm([u[0] + t * d[0] - p[0], u[1] + t * d[1] - p[1]]) <= 1e-12 * 1f64.max(norm(p))
}

fn on_ray(p: [f64; 2], u: [f64; 2], d: [f64; 2]) -> bool {
    let t = ((p[0] - u[0]) * d[0] + (p[1] - u[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1]);
    t >= -1e-12 && norm([u[0] + t * d[0] - p[0], u[1] + t * d[1] - p[1]]) <= 1e-12 * 1f64.max(norm(p))
}

fn parallel_same(a: [f64; 2], b: [f64; 2]) -> bool {
    (a[0] * b[1] - a[1] * b[0]).abs() <= 1e-12 * norm(a) * norm(b) && a[0] * b[0] + a[1] * b[1] > 0.0
}

fn sweeps_no_area(f1: &Form, f2: &Form) -> bool {
    if f1.kind == FKind::Piece || f2.kind == FKind::Piece {
        return false;
    }
    if f1.kind == FKind::Vertex && f2.kind == FKind::Vertex {
        return true;
    }
    let (f1, f2) = if f1.kind == FKind::Vertex { (f2, f1) } else { (f1, f2) };
    let on_line = |p: [f64; 2]| (f1.d[0] * (p[1] - f1.u[1]) - f1.d[1] * (p[0] - f1.u[0])).abs() <= 1e-12 * norm(f1.d) * 1f64.max(norm([p[0] - f1.u[0], p[1] - f1.u[1]]));
    if f2.kind == FKind::Vertex {
        return on_line(f2.u);
    }
    on_line(f2.u) && on_line([f2.u[0] + f2.d[0], f2.u[1] + f2.d[1]])
}

#[derive(Clone)]
enum Par {
    Line { p: [f64; 2], w: [f64; 2] },
    Rat { s1: Vec<f64>, s2: Vec<f64>, dn: Vec<f64> },
}

fn tie_arc(i: usize, j: usize, forms: &[Form], d: &[[f64; 3]]) -> bool {
    let q = conic_of(&forms[i], &forms[j]);
    for par in conic_params(&q) {
        let mut bp: Vec<f64> = vec![];
        for c in 0..forms.len() {
            for r in dom_rows(&forms[c]) {
                bp.extend(real_roots(&compose(&par, &[0.0, 0.0, 0.0, r[0], r[1], r[2]])));
            }
            if c != i && c != j {
                bp.extend(real_roots(&compose(&par, &conic_of(&forms[c], &forms[i]))));
            }
        }
        for r in d {
            bp.extend(real_roots(&compose(&par, &[0.0, 0.0, 0.0, r[0], r[1], r[2]])));
        }
        if let Par::Rat { dn, .. } = &par {
            bp.extend(real_roots(dn));
        }
        bp.sort_by(|a, b| a.total_cmp(b));   // total order: NaN-safe, never panics (review #6)
        bp.dedup();
        let samp: Vec<f64> = if bp.is_empty() {
            vec![0.0]
        } else {
            let mut s = vec![bp[0] - 1.0 - bp[0].abs()];
            for w in bp.windows(2) {
                s.push((w[0] + w[1]) / 2.0);
            }
            s.push(bp[bp.len() - 1] + 1.0 + bp[bp.len() - 1].abs());
            s
        };
        for t in samp {
            let (s, ok) = point_at(&par, t);
            if ok && both_active(s, i, j, forms, d) {
                return true;
            }
        }
    }
    false
}

fn both_active(s: [f64; 2], i: usize, j: usize, forms: &[Form], d: &[[f64; 3]]) -> bool {
    if !(in_dom(s, dom_rows(&forms[i]), -1.0) && in_dom(s, dom_rows(&forms[j]), -1.0) && in_dom(s, d, -1.0)) {
        return false;
    }
    let vi = gval(&forms[i], s);
    for (c, fc) in forms.iter().enumerate() {
        if c == i || c == j || !in_dom(s, dom_rows(fc), 1.0) {
            continue;
        }
        if gval(fc, s) > vi + 1e-9 * 1f64.max(vi.abs()) {
            return false;
        }
    }
    true
}

fn triangle_active(tr: &AtomTri, forms: &[Form], d: &[[f64; 3]]) -> bool {
    let s = tr.s;
    let mut tf = in_dom(s, d, -1.0);
    let v = tr.c;
    for (c, fc) in forms.iter().enumerate() {
        if tr.forms.contains(&(c + 1)) {
            continue;
        }
        if !tf {
            return tf;
        }
        if in_dom(s, dom_rows(fc), 1.0) && gval(fc, s) > v + 1e-9 * 1f64.max(v.abs()) {
            tf = false;
        }
    }
    tf
}

fn in_dom(s: [f64; 2], r: &[[f64; 3]], strictness: f64) -> bool {
    if r.is_empty() {
        return true;
    }
    let tol = 1e-9 * r.iter().flatten().fold(1f64, |a, x| a.max(x.abs())) * 1f64.max(norm(s));
    r.iter().all(|x| x[0] * s[0] + x[1] * s[1] + x[2] >= -strictness * tol)
}

fn inv2(m: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
    // numpy.linalg.inv via LU with partial pivoting (gesv against the identity)
    let c0 = solve2f(&m, [1.0, 0.0]);
    let c1 = solve2f(&m, [0.0, 1.0]);
    [[c0[0], c1[0]], [c0[1], c1[1]]]
}

fn dom_rows(fm: &Form) -> &[[f64; 3]] {
    fm.dom.get_or_init(|| match fm.kind {
        FKind::Vertex => vec![],
        FKind::Edge => {
            let al2 = 2.0 * fm.al;
            let r = [fm.d[0] / al2, fm.d[1] / al2, -fm.be / al2];
            let mut v = vec![r];
            if fm.tmax.is_finite() {
                v.push([-r[0], -r[1], fm.tmax - r[2]]);
            }
            v
        }
        FKind::Piece => {
            let hi = inv2(fm.h);
            let rm: Vec<[f64; 3]> = match fm.poly.as_ref().unwrap() {
                PieceRegion::Rows(r) => r.clone(),
                PieceRegion::Poly(p) => {
                    let n = p.len();
                    (0..n)
                        .map(|e| {
                            let a = p[e];
                            let b = p[(e + 1) % n];
                            let nr = [-(b[1] - a[1]), b[0] - a[0]];
                            [nr[0], nr[1], -(nr[0] * a[0] + nr[1] * a[1])]
                        })
                        .collect()
                }
            };
            rm.iter()
                .map(|r| {
                    let a = [r[0] * hi[0][0] + r[1] * hi[1][0], r[0] * hi[0][1] + r[1] * hi[1][1]];
                    [a[0], a[1], r[2] - (a[0] * fm.l[0] + a[1] * fm.l[1])]
                })
                .collect()
        }
    })
}

/// Eigen-decomposition of a symmetric 2x2: ascending eigenvalues, eigenvectors as columns.
fn eigh2(a: [[f64; 2]; 2]) -> ([f64; 2], [[f64; 2]; 2]) {
    let (p, q, r) = (a[0][0], a[0][1], a[1][1]);
    if q == 0.0 {
        if p <= r {
            return ([p, r], [[1.0, 0.0], [0.0, 1.0]]);
        }
        return ([r, p], [[0.0, 1.0], [1.0, 0.0]]);
    }
    let t = (p + r) / 2.0;
    let dd = ((p - r) / 2.0).hypot(q);
    let l1 = t - dd;
    let l2 = t + dd;
    // eigenvector for l2: (q, l2 - p) or (l2 - r, q), whichever is better conditioned
    let v2 = if (l2 - p).abs() > (l2 - r).abs() { [q, l2 - p] } else { [l2 - r, q] };
    let n2 = norm(v2);
    let v2 = [v2[0] / n2, v2[1] / n2];
    let v1 = [-v2[1], v2[0]];
    ([l1, l2], [[v1[0], v2[0]], [v1[1], v2[1]]])
}

fn conic_params(q: &[f64; 6]) -> Vec<Par> {
    let mut p = vec![];
    let sc = q.iter().fold(0.0f64, |a, x| a.max(x.abs()));
    if sc == 0.0 {
        return p;
    }
    let q: Vec<f64> = q.iter().map(|x| x / sc).collect();
    let a = [[q[0], q[1] / 2.0], [q[1] / 2.0, q[2]]];
    let g = [q[3], q[4]];
    let f0 = q[5];
    let (la, _) = eigh2(a);
    let norm_a = la[0].abs().max(la[1].abs());
    if norm_a <= 1e-12 {
        if norm(g) <= 1e-12 {
            return p;
        }
        let gg = g[0] * g[0] + g[1] * g[1];
        p.push(Par::Line { p: [-f0 * g[0] / gg, -f0 * g[1] / gg], w: [-g[1], g[0]] });
        return p;
    }
    let m = [[a[0][0], a[0][1], g[0] / 2.0], [a[1][0], a[1][1], g[1] / 2.0], [g[0] / 2.0, g[1] / 2.0, f0]];
    let det_m = crate::numeric::det3(&m);
    let norm_m = sym3_norm(&m);
    if det_m.abs() > 1e-10 * 1f64.max(norm_m).powi(3) {
        let p0 = match real_point(a, g, f0) {
            None => return p,
            Some(x) => x,
        };
        let gv = [2.0 * (a[0][0] * p0[0] + a[0][1] * p0[1]) + g[0], 2.0 * (a[1][0] * p0[0] + a[1][1] * p0[1]) + g[1]];
        let dn = vec![a[1][1], 2.0 * a[0][1], a[0][0]];
        let nn = [gv[1], gv[0]];
        let s1 = padd(&dn.iter().map(|x| p0[0] * x).collect::<Vec<_>>(), &[-nn[0], -nn[1]]);
        let s2 = padd(&dn.iter().map(|x| p0[1] * x).collect::<Vec<_>>(), &[-nn[0], -nn[1], 0.0]);
        p.push(Par::Rat { s1, s2, dn });
        return p;
    }
    let da = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if da > 1e-12 {
        return p;
    }
    if da < -1e-12 {
        let sol = solve2f(&a, g);
        let c0 = [-sol[0] / 2.0, -sol[1] / 2.0];
        let (l, v) = eigh2(a);
        // descending order
        let (l0, l1) = (l[1], l[0]);
        let v0 = [v[0][1], v[1][1]];
        let v1 = [v[0][0], v[1][0]];
        for sg in [1.0, -1.0] {
            let w = [v0[0] / l0.sqrt() + sg * v1[0] / (-l1).sqrt(), v0[1] / l0.sqrt() + sg * v1[1] / (-l1).sqrt()];
            p.push(Par::Line { p: c0, w });
        }
        return p;
    }
    let (l, v) = eigh2(a);
    let k = if l[0].abs() >= l[1].abs() { 0 } else { 1 };
    let lk = l[k];
    let vv = [v[0][k], v[1][k]];
    let vp = [-vv[1], vv[0]];
    for r in roots(&[lk, g[0] * vv[0] + g[1] * vv[1], f0]) {
        if r.im.abs() > 1e-9 * 1f64.max(r.re.hypot(r.im)) {
            continue;
        }
        p.push(Par::Line { p: [r.re * vv[0], r.re * vv[1]], w: vp });
    }
    p
}

// Jacobi eigenvalues: the indices address a[p][q]-style PAIRS, so they are load-bearing
#[allow(clippy::needless_range_loop)]
fn sym3_norm(m: &[[f64; 3]; 3]) -> f64 {
    // largest |eigenvalue| of a symmetric 3x3 (Jacobi)
    let mut a = *m;
    for _ in 0..50 {
        let mut off = 0.0;
        for p in 0..3 {
            for q in p + 1..3 {
                off += a[p][q].abs();
            }
        }
        if off < 1e-300 {
            break;
        }
        for p in 0..3 {
            for q in p + 1..3 {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..3 {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..3 {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    (0..3).map(|i| a[i][i].abs()).fold(0.0, f64::max)
}

fn real_point(a: [[f64; 2]; 2], g: [f64; 2], f0: f64) -> Option<[f64; 2]> {
    let (l, v) = eigh2(a);
    let mut base = vec![[0.0, 0.0]];
    if (a[0][0] * a[1][1] - a[0][1] * a[1][0]).abs() > 1e-12 {
        let s = solve2f(&a, g);
        base.push([-s[0] / 2.0, -s[1] / 2.0]);
    }
    let mut dirs: Vec<[f64; 2]> = vec![[v[0][0], v[1][0]], [v[0][1], v[1][1]], [1.0, 0.0], [0.0, 1.0]];
    // null(A): eigenvectors whose singular value is below 2 * smax * eps
    let smax = l[0].abs().max(l[1].abs());
    let tol = 2.0 * smax * EPS;
    let mut sv: Vec<(f64, [f64; 2])> = vec![(l[0].abs(), [v[0][0], v[1][0]]), (l[1].abs(), [v[0][1], v[1][1]])];
    sv.sort_by(|x, y| y.0.total_cmp(&x.0));
    for (s, vv) in sv {
        if !(s > tol) {
            dirs.push(vv);
        }
    }
    for p in &base {
        for w in &dirs {
            let aw = [a[0][0] * w[0] + a[0][1] * w[1], a[1][0] * w[0] + a[1][1] * w[1]];
            let ap = [a[0][0] * p[0] + a[0][1] * p[1], a[1][0] * p[0] + a[1][1] * p[1]];
            let c = [w[0] * aw[0] + w[1] * aw[1], (2.0 * ap[0] + g[0]) * w[0] + (2.0 * ap[1] + g[1]) * w[1], p[0] * ap[0] + p[1] * ap[1] + g[0] * p[0] + g[1] * p[1] + f0];
            let t: Vec<f64> = roots(&c).into_iter().filter(|z| z.im.abs() <= 1e-9 * 1f64.max(z.re.hypot(z.im))).map(|z| z.re).collect();
            // the root of SMALLEST magnitude: along a parabola's axis w'Aw is round-off and its
            // other root a spurious ~1e34 point that wrecks the parametrisation (review #9c)
            if let Some(t0) = t.iter().copied().min_by(|a, b| a.abs().total_cmp(&b.abs())) {
                return Some([p[0] + t0 * w[0], p[1] + t0 * w[1]]);
            }
        }
    }
    None
}

fn padd(a: &[f64], b: &[f64]) -> Vec<f64> {
    let n = a.len().max(b.len());
    let mut out = vec![0.0; n];
    for (i, x) in a.iter().enumerate() {
        out[n - a.len() + i] += x;
    }
    for (i, x) in b.iter().enumerate() {
        out[n - b.len() + i] += x;
    }
    out
}

fn compose(par: &Par, q: &[f64]) -> Vec<f64> {
    match par {
        Par::Line { p, w } => {
            let a = [[q[0], q[1] / 2.0], [q[1] / 2.0, q[2]]];
            let aw = [a[0][0] * w[0] + a[0][1] * w[1], a[1][0] * w[0] + a[1][1] * w[1]];
            let ap = [a[0][0] * p[0] + a[0][1] * p[1], a[1][0] * p[0] + a[1][1] * p[1]];
            vec![w[0] * aw[0] + w[1] * aw[1], (2.0 * ap[0] + q[3]) * w[0] + (2.0 * ap[1] + q[4]) * w[1], p[0] * ap[0] + p[1] * ap[1] + q[3] * p[0] + q[4] * p[1] + q[5]]
        }
        Par::Rat { s1, s2, dn } => {
            let sc = |v: &[f64], k: f64| -> Vec<f64> { v.iter().map(|x| k * x).collect() };
            padd(
                &padd(&padd(&sc(&convf(s1, s1), q[0]), &sc(&convf(s1, s2), q[1])), &sc(&convf(s2, s2), q[2])),
                &padd(&convf(&padd(&sc(s1, q[3]), &sc(s2, q[4])), dn), &sc(&convf(dn, dn), q[5])),
            )
        }
    }
}

fn trim_lead(p: &[f64]) -> Vec<f64> {
    if p.is_empty() {
        return vec![];
    }
    let m = p.iter().fold(1f64, |a, x| a.max(x.abs()));
    match p.iter().position(|x| x.abs() > 1e-12 * m) {
        None => vec![],
        Some(i) => p[i..].to_vec(),
    }
}

fn real_roots(c: &[f64]) -> Vec<f64> {
    let c = trim_lead(c);
    if c.len() < 2 {
        return vec![];
    }
    // Non-finite coefficients make the companion-matrix eigenvalue run produce NaN roots;
    // they are not breakpoints (review 2026-11-19 #6), and sorting them would panic.
    roots(&c).into_iter()
        .filter(|z| z.re.is_finite() && z.im.abs() <= 1e-6 * 1f64.max(z.re.hypot(z.im)))
        .map(|z| z.re).collect()
}

fn point_at(par: &Par, t: f64) -> ([f64; 2], bool) {
    match par {
        Par::Line { p, w } => ([p[0] + t * w[0], p[1] + t * w[1]], true),
        Par::Rat { s1, s2, dn } => {
            let d = polyval_f(dn, t);
            let ok = d.abs() > 1e-12 * dn.iter().fold(1f64, |a, x| a.max(x.abs())) * 1f64.max(t * t);
            if d == 0.0 {
                return ([f64::NAN, f64::NAN], false);
            }
            ([polyval_f(s1, t) / d, polyval_f(s2, t) / d], ok)
        }
    }
}

// ------------------------------------------------------------------------------------------
fn collect_atoms(f: &QuaPol) -> R<(Vec<Atom>, Vec<Form>)> {
    let mut atoms: Vec<Atom> = vec![];
    let mut forms: Vec<Form> = vec![];
    let mut cvx: Vec<PieceRegion> = vec![];
    let mut edge_at: Vec<[f64; 4]> = vec![];
    let mut ray_at: Vec<[f64; 4]> = vec![];
    let mut vtx: Vec<([f64; 2], f64, usize, usize)> = vec![];
    let vf: Vec<[f64; 2]> = (1..=f.vn.len()).map(|i| {
        let v = f.vx(i);
        [q_to_f64(&v[0]), q_to_f64(&v[1])]
    }).collect();
    for k in 1..=f.nf {
        let (n, d) = f.face_q(k);
        let h = [[idiv(&n[4], &d), idiv(&n[5], &d)], [idiv(&n[5], &d), idiv(&n[6], &d)]];
        let l = [idiv(&n[7], &d), idiv(&n[8], &d)];
        let k0 = idiv(&n[9], &d);
        let is_cvx = sign(&n[4]) >= 0 && sign(&n[6]) >= 0 && sign(&(&n[4] * &n[6] - &n[5] * &n[5])) >= 0;
        let q = |x: [f64; 2]| 0.5 * ((x[0] * h[0][0] + x[1] * h[1][0]) * x[0] + (x[0] * h[0][1] + x[1] * h[1][1]) * x[1]) + x[0] * l[0] + x[1] * l[1] + k0;
        let own: Vec<usize> = f.f.iter().enumerate().filter(|(_, r)| r[0] == k || r[1] == k).map(|(j, _)| j + 1).collect();
        struct Ed {
            u: [f64; 2],
            dv: [f64; 2],
            ray: bool,
            iu: usize,
            iw: usize,
        }
        let (poly, eidx, ed): (PieceRegion, Vec<usize>, Vec<Ed>) = if own.iter().any(|&j| f.e[j - 1][2] == 0) {
            let mut ed = vec![];
            for &j in &own {
                let (iu, iw) = (f.e[j - 1][0], f.e[j - 1][1]);
                let u = vf[iu - 1];
                let w = vf[iw - 1];
                ed.push(Ed { u, dv: [w[0] - u[0], w[1] - u[1]], ray: f.e[j - 1][2] == 0, iu, iw });
            }
            let rows = piece_rows(&ed.iter().map(|e| (e.u, e.dv, e.ray)).collect::<Vec<_>>());
            (PieceRegion::Rows(rows), own.clone(), ed)
        } else {
            let (order, eidx) = piece_polygon(f, k)?;
            let poly: Vec<[f64; 2]> = order.iter().map(|&i| vf[i - 1]).collect();
            let m = order.len();
            let mut ed = vec![];
            for i in 0..m {
                let u = poly[i];
                let w = poly[(i + 1) % m];
                ed.push(Ed { u, dv: [w[0] - u[0], w[1] - u[1]], ray: false, iu: order[i], iw: order[(i + 1) % m] });
            }
            (PieceRegion::Poly(poly), eidx, ed)
        };
        for e in &ed {
            add_vertex(&mut vtx, e.u, q(e.u), k, e.iu);
            if !e.ray {
                let w = [e.u[0] + e.dv[0], e.u[1] + e.dv[1]];
                add_vertex(&mut vtx, w, q(w), k, e.iw);
            }
        }
        if is_cvx {
            cvx.push(poly.clone());
            atoms.push(Atom { kind: AKind::Piece, q: single_piece(f, k, &eidx)?, region: Some(poly.clone()), u: [0.0; 2], w: [0.0; 2], dr: [0.0; 2] });
            if (&n[4] * &n[6] - &n[5] * &n[5]) > Z::zero() {
                let hi = inv2(h);
                let hil = [hi[0][0] * l[0] + hi[0][1] * l[1], hi[1][0] * l[0] + hi[1][1] * l[1]];
                let lhl = l[0] * (hi[0][0] * l[0] + hi[0][1] * l[1]) + l[1] * (hi[1][0] * l[0] + hi[1][1] * l[1]);
                let (gn, gd) = form_q(FKind::Piece, f, k, 0, 0)?;
                forms.push(Form { kind: FKind::Piece, a: hi, b: [-hil[0], -hil[1]], c: 0.5 * lhl - k0, u: [0.0; 2], d: [0.0; 2], al: 0.0, be: 0.0, h, l, poly: Some(poly.clone()), tmax: 1.0, gn, gd, dom: OnceLock::new() });
            }
        }
        for e in &ed {
            let (u, dv) = (e.u, e.dv);
            let al = 0.5 * ((dv[0] * h[0][0] + dv[1] * h[1][0]) * dv[0] + (dv[0] * h[0][1] + dv[1] * h[1][1]) * dv[1]);
            if al <= 0.0 {
                continue;
            }
            let be = (u[0] * h[0][0] + u[1] * h[1][0] + l[0]) * dv[0] + (u[0] * h[0][1] + u[1] * h[1][1] + l[1]) * dv[1];
            let ga = q(u);
            let (gn, gd) = form_q(FKind::Edge, f, k, e.iu, e.iw)?;
            let a2 = 2.0 * al;
            let fm = Form {
                kind: FKind::Edge,
                a: [[dv[0] * dv[0] / a2, dv[0] * dv[1] / a2], [dv[1] * dv[0] / a2, dv[1] * dv[1] / a2]],
                b: [u[0] - be * dv[0] / a2, u[1] - be * dv[1] / a2],
                c: -ga + be.powi(2) / (4.0 * al),
                u,
                d: dv,
                al,
                be,
                h: [[0.0; 2]; 2],
                l: [0.0; 2],
                poly: None,
                tmax: if e.ray { f64::INFINITY } else { 1.0 },
                gn,
                gd,
                dom: OnceLock::new(),
            };
            forms.push(fm);
            if is_cvx {
                continue;
            }
            if e.ray {
                let nd = norm(dv);
                let key = [u[0], u[1], dv[0] / nd, dv[1] / nd];
                if ray_at.iter().any(|r| r.iter().zip(key.iter()).all(|(a, b)| (a - b).abs() < 1e-12)) {
                    continue;
                }
                ray_at.push(key);
                let qp = QuaPol::new_q(vec![f.vx(e.iu), f.vx(e.iw)], vec![[1, 2, 0]], vec![f.fx(k)], vec![[0, 0]])?;
                atoms.push(Atom { kind: AKind::Ray, q: qp, region: None, u, w: [0.0; 2], dr: dv });
                continue;
            }
            let w = [u[0] + dv[0], u[1] + dv[1]];
            let mut pr = [u, w];
            pr.sort_by(|a, b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));   // NaN-safe (review #6)
            let key = [pr[0][0], pr[0][1], pr[1][0], pr[1][1]];
            if edge_at.iter().any(|r| r.iter().zip(key.iter()).all(|(a, b)| (a - b).abs() < 1e-12)) {
                continue;
            }
            edge_at.push(key);
            let qp = QuaPol::new_q(vec![f.vx(e.iu), f.vx(e.iw)], vec![[1, 2, 1]], vec![f.fx(k)], vec![[0, 0]])?;
            atoms.push(Atom { kind: AKind::Edge, q: qp, region: None, u, w, dr: [0.0; 2] });
        }
    }
    for (p, val, face, iv) in &vtx {
        let (gn, gd) = form_q(FKind::Vertex, f, *face, *iv, 0)?;
        forms.push(Form { kind: FKind::Vertex, a: [[0.0; 2]; 2], b: *p, c: -val, u: *p, d: [0.0; 2], al: 0.0, be: 0.0, h: [[0.0; 2]; 2], l: [0.0; 2], poly: None, tmax: 1.0, gn, gd, dom: OnceLock::new() });
        let mut covered = cvx.iter().any(|pp| in_piece(*p, pp));
        for r in &edge_at {
            covered = covered || norm([p[0] - r[0], p[1] - r[1]]) < 1e-12 || norm([p[0] - r[2], p[1] - r[3]]) < 1e-12;
        }
        for r in &ray_at {
            covered = covered || norm([p[0] - r[0], p[1] - r[1]]) < 1e-12;
        }
        if !covered {
            let (nq, dq) = f.face_q(*face);
            let vx = f.vx(*iv);
            let vex = eval_face_exact(&nq, &dq, &vx[0], &vx[1]);
            let mut row = vec![qi(0); 9];
            row.push(vex);
            let qp = QuaPol::new_q(vec![vx], vec![], vec![row], vec![])?;
            atoms.push(Atom { kind: AKind::Vertex, q: qp, region: None, u: *p, w: [0.0; 2], dr: [0.0; 2] });
        }
    }
    let atoms: Vec<Atom> = atoms.into_iter().filter(|at| !(at.kind == AKind::Edge && cvx.iter().any(|pp| in_piece(at.u, pp) && in_piece(at.w, pp)))).collect();
    Ok((atoms, forms))
}

fn add_vertex(vtx: &mut Vec<([f64; 2], f64, usize, usize)>, p: [f64; 2], v: f64, k: usize, iv: usize) {
    if vtx.iter().any(|r| (r.0[0] - p[0]).abs() < 1e-12 && (r.0[1] - p[1]).abs() < 1e-12) {
        return;
    }
    vtx.push((p, v, k, iv));
}

fn form_q(kind: FKind, f: &QuaPol, k: usize, iu: usize, iw: usize) -> R<(Vec<Z>, Z)> {
    let (n, d) = f.face_q(k);
    let fr = |x: &Z| Q::new(x.clone(), d.clone());
    let h = [[fr(&n[4]), fr(&n[5])], [fr(&n[5]), fr(&n[6])]];
    let l = [fr(&n[7]), fr(&n[8])];
    let k0 = fr(&n[9]);
    let half = Q::new(zi(1), zi(2));
    let quadf = |a: &[Q; 2], b: &[Q; 2]| -> Q {
        let mut s = qi(0);
        for i in 0..2 {
            for j in 0..2 {
                s += &a[i] * &h[i][j] * &b[j];
            }
        }
        s
    };
    let qv = |u: &[Q; 2]| -> Q { &half * quadf(u, u) + &l[0] * &u[0] + &l[1] * &u[1] + &k0 };
    let c: Vec<Q> = match kind {
        FKind::Vertex => {
            let u = f.vx(iu);
            vec![qi(0), qi(0), qi(0), u[0].clone(), u[1].clone(), -qv(&u)]
        }
        FKind::Edge => {
            let u = f.vx(iu);
            let w = f.vx(iw);
            let dv = [&w[0] - &u[0], &w[1] - &u[1]];
            let al = &half * quadf(&dv, &dv);
            let hu = [&h[0][0] * &u[0] + &h[0][1] * &u[1] + &l[0], &h[1][0] * &u[0] + &h[1][1] * &u[1] + &l[1]];
            let be = &hu[0] * &dv[0] + &hu[1] * &dv[1];
            let i4 = qi(1) / (qi(4) * &al);
            let i2 = qi(1) / (qi(2) * &al);
            vec![&dv[0] * &dv[0] * &i4, qi(2) * &dv[0] * &dv[1] * &i4, &dv[1] * &dv[1] * &i4, &u[0] - &be * &dv[0] * &i2, &u[1] - &be * &dv[1] * &i2, &be * &be * &i4 - qv(&u)]
        }
        FKind::Piece => {
            let dt = &h[0][0] * &h[1][1] - &h[0][1] * &h[0][1];
            let hi = [[&h[1][1] / &dt, -&h[0][1] / &dt], [-&h[0][1] / &dt, &h[0][0] / &dt]];
            let hil = [&hi[0][0] * &l[0] + &hi[0][1] * &l[1], &hi[1][0] * &l[0] + &hi[1][1] * &l[1]];
            vec![&half * &hi[0][0], hi[0][1].clone(), &half * &hi[1][1], -hil[0].clone(), -hil[1].clone(), &half * (&l[0] * &hil[0] + &l[1] * &hil[1]) - &k0]
        }
    };
    let mut dd = Z::one();
    for x in &c {
        dd = lcm(&dd, x.denom());
    }
    let nn: Vec<Z> = c.iter().map(|x| (x * qz(&dd)).to_integer()).collect();
    canon(&nn, &dd)
}

fn piece_polygon(f: &QuaPol, k: usize) -> R<(Vec<usize>, Vec<usize>)> {
    let eidx: Vec<usize> = f.f.iter().enumerate().filter(|(_, r)| r[0] == k || r[1] == k).map(|(j, _)| j + 1).collect();
    let mut e = vec![];
    for &j in &eidx {
        let (mut a, mut b) = (f.e[j - 1][0], f.e[j - 1][1]);
        if f.f[j - 1][0] != k {
            std::mem::swap(&mut a, &mut b);
        }
        e.push([a, b]);
    }
    let mut order = vec![e[0][0]];
    let mut nxt = e[0][1];
    while nxt != order[0] {
        order.push(nxt);
        let cand: Vec<usize> = e.iter().filter(|x| x[0] == nxt).map(|x| x[1]).collect();
        if cand.len() != 1 {
            bail!("biconjAtoms:subdivision", "face {} is not a simple polygon.", k);
        }
        nxt = cand[0];
    }
    Ok((order, eidx))
}

fn single_piece(f: &QuaPol, k: usize, eidx: &[usize]) -> R<QuaPol> {
    let vs: Vec<usize> = eidx.iter().flat_map(|&j| [f.e[j - 1][0], f.e[j - 1][1]]).collect::<BTreeSet<_>>().into_iter().collect();
    let idx: BTreeMap<usize, usize> = vs.iter().enumerate().map(|(i, v)| (*v, i + 1)).collect();
    let e: Vec<[usize; 3]> = eidx.iter().map(|&j| [idx[&f.e[j - 1][0]], idx[&f.e[j - 1][1]], f.e[j - 1][2]]).collect();
    let ff: Vec<[usize; 2]> = eidx.iter().map(|&j| [(f.f[j - 1][0] == k) as usize, (f.f[j - 1][1] == k) as usize]).collect();
    QuaPol::new_q(vs.iter().map(|&v| f.vx(v)).collect(), e, vec![f.fx(k)], ff)
}

fn in_piece(p: [f64; 2], pr: &PieceRegion) -> bool {
    match pr {
        PieceRegion::Rows(r) => {
            let m = r.iter().flatten().fold(1f64, |a, x| a.max(x.abs()));
            r.iter().all(|x| x[0] * p[0] + x[1] * p[1] + x[2] >= -1e-9 * m * 1f64.max(norm(p)))
        }
        PieceRegion::Poly(pp) => in_poly(p, pp),
    }
}

fn piece_rows(ed: &[([f64; 2], [f64; 2], bool)]) -> Vec<[f64; 3]> {
    let n = ed.len() as f64;
    let mut ip = [ed.iter().map(|e| e.0[0]).sum::<f64>() / n, ed.iter().map(|e| e.0[1]).sum::<f64>() / n];
    for e in ed.iter().filter(|e| e.2) {
        let nr = norm(e.1);
        ip[0] += e.1[0] / nr;
        ip[1] += e.1[1] / nr;
    }
    let mut r = vec![];
    for e in ed {
        let nrm = [-e.1[1], e.1[0]];
        let mut row = [nrm[0], nrm[1], -(nrm[0] * e.0[0] + nrm[1] * e.0[1])];
        if row[0] * ip[0] + row[1] * ip[1] + row[2] < 0.0 {
            row = [-row[0], -row[1], -row[2]];
        }
        r.push(row);
    }
    r
}

fn in_poly(p: [f64; 2], pp: &[[f64; 2]]) -> bool {
    let n = pp.len();
    let m = pp.iter().flatten().fold(1f64, |a, x| a.max(x.abs())).powi(2);
    for i in 0..n {
        let a = pp[i];
        let b = pp[(i + 1) % n];
        if (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]) < -1e-9 * m {
            return false;
        }
    }
    true
}

// ------------------------------------------------------------------------------------------
fn triple_faces(forms: &[Form], d: &[[f64; 3]]) -> Vec<AtomTri> {
    let mut tris = vec![];
    let n = forms.len();
    for a in 0..n {
        for b in a + 1..n {
            for c in b + 1..n {
                let f3 = [&forms[a], &forms[b], &forms[c]];
                for s0 in common_ties(&f3) {
                    let mut m = [[0.0; 2]; 3];
                    let mut ok = true;
                    for t in 0..3 {
                        let (mm, ok1) = contact(f3[t], s0);
                        m[t] = mm;
                        ok = ok && ok1;
                    }
                    if !ok {
                        continue;
                    }
                    let ar = (m[1][0] - m[0][0]) * (m[2][1] - m[0][1]) - (m[1][1] - m[0][1]) * (m[2][0] - m[0][0]);
                    if ar.abs() <= 1e-9 * m.iter().flatten().fold(1f64, |acc, x| acc.max(x.abs())).powi(2) {
                        continue;
                    }
                    tris.push(AtomTri {
                        t: m,
                        s: s0,
                        c: gval(f3[0], s0),
                        g: normal_cone(d, s0),
                        forms: vec![a + 1, b + 1, c + 1],
                        gn: f3.iter().filter(|f| !f.gn.is_empty()).map(|f| f.gn.clone()).collect(),
                        gd: f3.iter().filter(|f| !f.gn.is_empty()).map(|f| f.gd.clone()).collect(),
                    });
                }
            }
        }
    }
    tris
}

fn normal_cone(d: &[[f64; 3]], s0: [f64; 2]) -> Vec<[f64; 2]> {
    if d.is_empty() {
        return vec![];
    }
    let r: Vec<f64> = d.iter().map(|x| x[0] * s0[0] + x[1] * s0[1] + x[2]).collect();
    let tol = 1e-9 * d.iter().flatten().fold(1f64, |a, x| a.max(x.abs()));
    if r.iter().any(|x| *x < -tol) {
        return vec![];
    }
    d.iter().zip(r.iter()).filter(|(_, x)| x.abs() <= tol * 1f64.max(norm(s0))).map(|(dd, _)| [-dd[0], -dd[1]]).collect()
}

fn gval(fm: &Form, s: [f64; 2]) -> f64 {
    let a = fm.a;
    0.5 * (s[0] * (a[0][0] * s[0] + a[0][1] * s[1]) + s[1] * (a[1][0] * s[0] + a[1][1] * s[1])) + fm.b[0] * s[0] + fm.b[1] * s[1] + fm.c
}

fn contact(fm: &Form, s: [f64; 2]) -> ([f64; 2], bool) {
    match fm.kind {
        FKind::Vertex => (fm.u, true),
        FKind::Edge => {
            let t = (s[0] * fm.d[0] + s[1] * fm.d[1] - fm.be) / (2.0 * fm.al);
            ([fm.u[0] + t * fm.d[0], fm.u[1] + t * fm.d[1]], t >= -1e-9 && t <= fm.tmax + 1e-9)
        }
        FKind::Piece => {
            let m = solve2f(&fm.h, [s[0] - fm.l[0], s[1] - fm.l[1]]);
            (m, in_piece(m, fm.poly.as_ref().unwrap()))
        }
    }
}

fn conic_of(fa: &Form, fb: &Form) -> [f64; 6] {
    let a = [[fa.a[0][0] - fb.a[0][0], fa.a[0][1] - fb.a[0][1]], [fa.a[1][0] - fb.a[1][0], fa.a[1][1] - fb.a[1][1]]];
    let b = [fa.b[0] - fb.b[0], fa.b[1] - fb.b[1]];
    [a[0][0] / 2.0, a[0][1], a[1][1] / 2.0, b[0], b[1], fa.c - fb.c]
}

fn conic_at(c: &[f64; 6], s: [f64; 2]) -> f64 {
    c[0] * s[0].powi(2) + c[1] * s[0] * s[1] + c[2] * s[1].powi(2) + c[3] * s[0] + c[4] * s[1] + c[5]
}

fn conic_grad(c: &[f64; 6], s: [f64; 2]) -> [f64; 2] {
    [2.0 * c[0] * s[0] + c[1] * s[1] + c[3], c[1] * s[0] + 2.0 * c[2] * s[1] + c[4]]
}

fn common_ties(f3: &[&Form; 3]) -> Vec<[f64; 2]> {
    let c1 = conic_of(f3[0], f3[1]);
    let c2 = conic_of(f3[0], f3[2]);
    let mut s: Vec<[f64; 2]> = vec![];
    if c1.iter().all(|x| x.abs() < 1e-14) || c2.iter().all(|x| x.abs() < 1e-14) {
        return s;
    }
    for swap in [false, true] {
        let (p1, p2) = if swap {
            ([c1[2], c1[1], c1[0], c1[4], c1[3], c1[5]], [c2[2], c2[1], c2[0], c2[4], c2[3], c2[5]])
        } else {
            (c1, c2)
        };
        let mut cand = conic_meet_double(&p1, &p2);
        if swap {
            cand = cand.into_iter().map(|q| [q[1], q[0]]).collect();
        }
        let any = !cand.is_empty();
        s.extend(cand);
        if any {
            break;
        }
    }
    let mut out: Vec<[f64; 2]> = vec![];
    let mut kept: Vec<bool> = vec![];
    for s0 in s {
        let mut s = s0;
        for _ in 0..30 {
            let r = [conic_at(&c1, s), conic_at(&c2, s)];
            let j = [conic_grad(&c1, s), conic_grad(&c2, s)];
            if rcond2(&j) < 1e-14 {
                break;
            }
            let ds = solve2f(&j, [-r[0], -r[1]]);
            s = [s[0] + ds[0], s[1] + ds[1]];
            if norm(ds) <= 1e-15 * 1f64.max(norm(s)) {
                break;
            }
        }
        let sc = 1f64.max(norm(s)).powi(2) * c1.iter().chain(c2.iter()).fold(1f64, |a, x| a.max(x.abs()));
        let mut keep = s[0].is_finite() && s[1].is_finite() && conic_at(&c1, s).abs() < 1e-8 * sc && conic_at(&c2, s).abs() < 1e-8 * sc;
        for (sj, kj) in out.iter().zip(kept.iter()) {
            if *kj && norm([sj[0] - s[0], sj[1] - s[1]]) < 1e-9 * 1f64.max(norm(s)) {
                keep = false;
            }
        }
        out.push(s);
        kept.push(keep);
    }
    out.into_iter().zip(kept).filter(|(_, k)| *k).map(|(s, _)| s).collect()
}

fn conic_meet_double(c1: &[f64; 6], c2: &[f64; 6]) -> Vec<[f64; 2]> {
    let a1 = c1[2];
    let b1 = [c1[1], c1[4]];
    let cc1 = [c1[0], c1[3], c1[5]];
    let a2 = c2[2];
    let b2 = [c2[1], c2[4]];
    let cc2 = [c2[0], c2[3], c2[5]];
    let rr = if a1 == 0.0 && a2 == 0.0 {
        padd(&convf(&b1, &cc2), &convf(&b2, &cc1).iter().map(|x| -x).collect::<Vec<_>>())
    } else {
        let x = padd(&cc2.iter().map(|v| a1 * v).collect::<Vec<_>>(), &cc1.iter().map(|v| -a2 * v).collect::<Vec<_>>());
        let y = padd(&b2.iter().map(|v| a1 * v).collect::<Vec<_>>(), &b1.iter().map(|v| -a2 * v).collect::<Vec<_>>());
        let z = padd(&convf(&b1, &cc2), &convf(&b2, &cc1).iter().map(|v| -v).collect::<Vec<_>>());
        padd(&convf(&x, &x), &convf(&y, &z).iter().map(|v| -v).collect::<Vec<_>>())
    };
    let rr = trim_lead(&rr);
    let mut s = vec![];
    if rr.is_empty() || rr.iter().all(|x| *x == 0.0) {
        return s;
    }
    let r: Vec<f64> = if rr.len() < 2 { vec![] } else { roots(&rr).into_iter().filter(|z| z.im.abs() <= 1e-7 * 1f64.max(z.re.hypot(z.im))).map(|z| z.re).collect() };
    for s1 in r {
        let a = a1;
        let bb = polyval_f(&b1, s1);
        let cc = polyval_f(&cc1, s1);
        let mut s2c: Vec<crate::numeric::Cx> = if a.abs() > 0.0 {
            roots(&[a, bb, cc])
        } else if bb != 0.0 {
            vec![crate::numeric::Cx { re: -cc / bb, im: 0.0 }]
        } else {
            vec![]
        };
        let b2v = polyval_f(&b2, s1);
        let cc2v = polyval_f(&cc2, s1);
        if s2c.is_empty() || (a.abs() == 0.0 && bb == 0.0) {
            if a2.abs() > 0.0 {
                s2c = roots(&[a2, b2v, cc2v]);
            } else if b2v != 0.0 {
                s2c = vec![crate::numeric::Cx { re: -cc2v / b2v, im: 0.0 }];
            }
        }
        for z in s2c {
            if z.im.abs() <= 1e-7 * 1f64.max(z.re.hypot(z.im)) {
                s.push([s1, z.re]);
            }
        }
    }
    s
}


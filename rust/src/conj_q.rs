//! conj_q -- the Legendre-Fenchel conjugate computed EXACTLY over the rationals (port of
//! conj_q.py / conjQ.m). Returns (g, dom): g a QuaCon and dom integer half-planes [a b c].

use crate::bail;
use crate::cells::{assemble_qua_con_cells, con_row, Cell};
use crate::errors::R;
use crate::max_q::max_q;
use crate::qua_con::QuaCon;
use crate::qua_pol::QuaPol;
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeSet;

pub type Dom = Vec<Vec<Z>>;

pub fn assert_no_closed_edges(e: &[[usize; 3]], who: &str) -> R<()> {
    let bad: Vec<usize> = e.iter().enumerate().filter(|(_, r)| r[2] == 2).map(|(j, _)| j + 1).collect();
    if !bad.is_empty() {
        bail!("RatCon:closedEdge", "{} cannot yet handle a CLOSED edge (E(:,3) = 2). Edge(s): {:?}.", who, bad);
    }
    Ok(())
}

pub fn conj_q(obj: &QuaPol) -> R<(QuaCon, Dom)> {
    assert_no_closed_edges(&obj.e, "conjQ")?;
    obj.assert_operable()?;
    obj.assert_exact()?;
    let nf = obj.fn_.len();
    if obj.dom.dim < 2 && nf > 1 {
        if nf != obj.ne {
            bail!("PLQ:conjQ:ambiguousChain", "a domain of dimension < 2 carrying {} function rows on {} edges.", nf, obj.ne);
        }
        let mut g: Option<QuaCon> = None;
        let mut dom = vec![];
        for j in 1..=obj.ne {
            let e = obj.e[j - 1];
            let ej = QuaPol::new_q(vec![obj.vx(e[0]), obj.vx(e[1])], vec![[1, 2, e[2]]], vec![obj.fx(j)], vec![[0, 0]])?;
            let (gj, dj) = conj_piece_q(&ej, 1)?;
            g = Some(match g {
                None => gj,
                Some(g0) => max_q(&g0, &gj)?,
            });
            dom.extend(dj);
        }
        return Ok((g.unwrap(), dom));
    }
    // the pieces are independent (map), their conjugates are combined by a maximum (reduction)
    let parts = crate::par::map_range(nf, |k| conj_piece_q(obj, k + 1))?;
    let mut dom = vec![];
    let mut gs = Vec::with_capacity(nf);
    for (g, d) in parts {
        gs.push(g);
        dom.extend(d);
    }
    let g = crate::par::reduce(gs, |a, b| max_q(&a, &b))?;
    Ok((g, dom))
}

#[derive(Clone, Debug)]
pub struct Edge {
    pub a: usize,
    pub b: usize,
    pub d: Vec<Z>,
    pub is_ray: bool,
    pub j: usize,
    pub n_in: Option<Vec<Z>>,
}

#[derive(Clone, Debug)]
pub struct Shape {
    pub vs: Vec<usize>,
    pub ed: Vec<Edge>,
    pub rays: Vec<Vec<Z>>,
    pub hp: Vec<Vec<Z>>,
    pub bounded: bool,
    pub convex: bool,
}

fn hmat(f: &[Z]) -> Vec<Vec<Z>> {
    vec![vec![f[4].clone(), f[5].clone()], vec![f[5].clone(), f[6].clone()]]
}

pub fn conj_piece_q(obj: &QuaPol, k: usize) -> R<(QuaCon, Dom)> {
    let (fnk, fdk) = obj.face_q(k);
    let hn = hmat(&fnk);
    let d = det_exact(&hn)?;
    let dom: Dom = vec![];
    if obj.nv == 0 {
        return Ok((case_a_full_domain(&fnk, &fdk)?, dom));
    }
    let (vi, vd) = combine_den(&obj.vn, &obj.vd)?;
    let ln = vec![fnk[7].clone(), fnk[8].clone()];
    if obj.dom.dim < 2 {
        if obj.fn_.len() > 1 {
            bail!("PLQ:conjQ:ambiguousChain", "a domain of dimension < 2 carrying {} different functions.", obj.fn_.len());
        }
        let sh = degenerate_shape(obj, &vi);
        let (ok, dom, _) = recession_conditions(&sh, &hn, &ln, &fdk, &vi, &vd)?;
        if !ok {
            return Ok((QuaCon::empty(), vec![zv(&[0, 0, -1])]));
        }
        let g = case_d_boundary_max(&fnk, &fdk, &vi, &vd, &sh, &dom)?;
        return Ok((g, dom));
    }
    let sh = piece_shape(obj, k, &vi, &vd)?;
    if sh.vs.is_empty() {
        bail!("PLQ:conjQ:noVertex", "piece {} has no vertex.", k);
    }
    if !sh.convex {
        if !sh.bounded {
            return conj_unbounded_non_convex(obj, k);
        }
        let tris = ear_clip(&vi, boundary_cycle(&sh, &vi)?)?;
        // the triangles are independent (map), their conjugates are combined by a maximum (reduction)
        let gts = crate::par::map(&tris, |tri| {
            let sht = triangle_shape(&vi, &vd, tri)?;
            if hn[0][0].is_positive() && d.is_positive() {
                case_b_convex_on_piece(&fnk, &fdk, &vi, &vd, &sht, &[])
            } else {
                case_d_boundary_max(&fnk, &fdk, &vi, &vd, &sht, &[])
            }
        })?;
        let g = crate::par::reduce(gts, |a, b| max_q(&a, &b))?;
        return Ok((g, dom));
    }
    let (ok, dom, _why) = recession_conditions(&sh, &hn, &ln, &fdk, &vi, &vd)?;
    if !ok {
        return Ok((QuaCon::empty(), vec![zv(&[0, 0, -1])]));
    }
    if let Some(gt) = lineality_conj(obj, k, &sh, &hn, &fnk, &fdk, &vi, &vd)? {
        return Ok((gt, dom));
    }
    let g = if hn[0][0].is_positive() && d.is_positive() {
        case_b_convex_on_piece(&fnk, &fdk, &vi, &vd, &sh, &dom)?
    } else {
        case_d_boundary_max(&fnk, &fdk, &vi, &vd, &sh, &dom)?
    };
    Ok((g, dom))
}

fn degenerate_shape(obj: &QuaPol, vi: &[Vec<Z>]) -> Shape {
    let mut sh = Shape { vs: vec![], ed: vec![], rays: vec![], hp: vec![], bounded: true, convex: true };
    for e in &obj.e {
        let (a, b) = (e[0], e[1]);
        let d = vec![&vi[b - 1][0] - &vi[a - 1][0], &vi[b - 1][1] - &vi[a - 1][1]];
        if e[2] != 0 {
            sh.vs.extend([a, b]);
            sh.ed.push(Edge { a, b, d, is_ray: false, j: 0, n_in: None });
        } else {
            sh.vs.push(a);
            sh.ed.push(Edge { a, b: 0, d: d.clone(), is_ray: true, j: 0, n_in: None });
            sh.rays.push(d);
            sh.bounded = false;
        }
    }
    if sh.ed.is_empty() {
        sh.vs = (1..=obj.vn.len()).collect();
    }
    sh.vs = sh.vs.iter().cloned().collect::<BTreeSet<_>>().into_iter().collect();
    sh
}

pub fn case_a_full_domain(fnk: &[Z], fd: &Z) -> R<QuaCon> {
    let hn = hmat(fnk);
    let ln = vec![fnk[7].clone(), fnk[8].clone()];
    let d = det_exact(&hn)?;
    if !(hn[0][0].is_positive() && d.is_positive()) {
        if !is_psd2(&hn) {
            return Ok(QuaCon::empty());
        }
        if hn.iter().flatten().all(|x| x.is_zero()) {
            let r1 = vec![fd.clone(), zi(0), -&ln[0]];
            let r2 = vec![zi(0), fd.clone(), -&ln[1]];
            let con = vec![con_row(conic_lin(&r1)?, sgn_of(&r1, 0)?, 0), con_row(conic_lin(&r2)?, sgn_of(&r2, 0)?, 0)];
            let mut num = vec![zi(0); 9];
            num.push(-&fnk[9]);
            return assemble_qua_con_cells(vec![Cell { num, den: fd.clone(), con }]);
        }
        let nl = if !hn[0][0].is_zero() { vec![-&hn[0][1], hn[0][0].clone()] } else { zv(&[1, 0]) };
        let mut m = vec![hn[0][0].clone(), hn[0][1].clone()];
        if all_zero(&m) {
            m = vec![hn[0][1].clone(), hn[1][1].clone()];
        }
        let mhm = quad(&m, &hn, &m);
        let p = [fd * &m[0], fd * &m[1], -dot(&ln, &m)];
        let den = zi(2) * &mhm * fd;
        let num = vec![
            zi(0),
            zi(0),
            zi(0),
            zi(0),
            zi(2) * &p[0] * &p[0],
            zi(2) * &p[0] * &p[1],
            zi(2) * &p[1] * &p[1],
            zi(2) * &p[0] * &p[2],
            zi(2) * &p[1] * &p[2],
            &p[2] * &p[2] - zi(2) * &mhm * &fnk[9],
        ];
        let row = vec![fd * &nl[0], fd * &nl[1], -dot(&ln, &nl)];
        return assemble_qua_con_cells(vec![Cell { num, den, con: vec![con_row(conic_lin(&row)?, sgn_of(&row, 0)?, 0)] }]);
    }
    let a = adj(&hn);
    let aln = mv(&a, &ln);
    let den = zi(2) * fd * &d;
    let fd2 = fd * fd;
    let num = vec![
        zi(0),
        zi(0),
        zi(0),
        zi(0),
        zi(2) * &fd2 * &a[0][0],
        zi(2) * &fd2 * &a[0][1],
        zi(2) * &fd2 * &a[1][1],
        -zi(2) * fd * &aln[0],
        -zi(2) * fd * &aln[1],
        dot(&ln, &aln) - zi(2) * &d * &fnk[9],
    ];
    let (gn, gd) = canon(&num, &den)?;
    QuaCon::new(vec![], vec![], vec![], vec![gn], vec![gd], vec![], vec![vec![]])
}

fn adj(hn: &[Vec<Z>]) -> Vec<Vec<Z>> {
    vec![vec![hn[1][1].clone(), -&hn[0][1]], vec![-&hn[0][1], hn[0][0].clone()]]
}

fn vplus(a: &[Z], b: &[Z], s: &Z) -> Vec<Z> {
    a.iter().zip(b.iter()).map(|(x, y)| x + s * y).collect()
}

fn qval(v: &[Z], hn: &[Vec<Z>], ln: &[Z], kn: &Z, vd: &Z) -> Z {
    quad(v, hn, v) + zi(2) * vd * dot(ln, v) + zi(2) * vd * vd * kn
}

pub fn case_b_convex_on_piece(fnk: &[Z], fd: &Z, vi: &[Vec<Z>], vd: &Z, sh: &Shape, dom: &[Vec<Z>]) -> R<QuaCon> {
    let hn = hmat(fnk);
    let ln = vec![fnk[7].clone(), fnk[8].clone()];
    let kn = &fnk[9];
    let d = det_exact(&hn)?;
    if !(hn[0][0].is_positive() && d.is_positive()) {
        bail!("PLQ:conjQ:notStrictlyConvex", "caseB needs a strictly convex Q.");
    }
    let den_v = zi(2) * fd * vd * vd;
    let fdvd = fd * vd;
    let mut cells = vec![];
    for &vix in &sh.vs {
        let v = &vi[vix - 1];
        let gv = vplus(&mv(&hn, v), &ln, vd);
        let qv = qval(v, &hn, &ln, kn, vd);
        let mut num = vec![zi(0); 7];
        num.extend([zi(2) * &fdvd * &v[0], zi(2) * &fdvd * &v[1], -qv]);
        let mut rows: Vec<Vec<Z>> = dom.to_vec();
        for ed in dirs_at_vertex(sh, vix) {
            rows.push(vec![-&fdvd * &ed[0], -&fdvd * &ed[1], dot(&gv, &ed)]);
        }
        if let Some(con) = to_con(&rows)? {
            cells.push(Cell { num, den: den_v.clone(), con });
        }
    }
    for (j, e) in sh.ed.iter().enumerate() {
        let a = &vi[e.a - 1];
        let dd = &e.d;
        let alpha = quad(dd, &hn, dd);
        if !alpha.is_positive() {
            bail!("PLQ:conjQ:degenerateEdge", "edge {} has d'Qd = {} <= 0 under a strictly convex Q.", j + 1, alpha);
        }
        let ga = vplus(&mv(&hn, a), &ln, vd);
        let qa = qval(a, &hn, &ln, kn, vd);
        let td = vec![&fdvd * &dd[0], &fdvd * &dd[1], -dot(&ga, dd)];
        let den = zi(2) * &alpha * fd * vd * vd;
        let lin = [zi(2) * &alpha * &fdvd * &a[0], zi(2) * &alpha * &fdvd * &a[1]];
        let num = vec![
            zi(0),
            zi(0),
            zi(0),
            zi(0),
            zi(2) * &td[0] * &td[0],
            zi(2) * &td[0] * &td[1],
            zi(2) * &td[1] * &td[1],
            &lin[0] + zi(2) * &td[0] * &td[2],
            &lin[1] + zi(2) * &td[1] * &td[2],
            &td[2] * &td[2] - &alpha * &qa,
        ];
        let mut rows: Vec<Vec<Z>> = dom.to_vec();
        rows.push(td.clone());
        if !e.is_ray {
            rows.push(vec![-&td[0], -&td[1], &alpha - &td[2]]);
        }
        let n = outward_normal_of(sh, j);
        let hdn = quad(dd, &hn, &n);
        rows.push(vec![
            &alpha * &fdvd * &n[0] - &td[0] * &hdn,
            &alpha * &fdvd * &n[1] - &td[1] * &hdn,
            -&alpha * dot(&ga, &n) - &td[2] * &hdn,
        ]);
        if let Some(con) = to_con(&rows)? {
            cells.push(Cell { num, den, con });
        }
    }
    let a = adj(&hn);
    let aln = mv(&a, &ln);
    let fd2 = fd * fd;
    let num_i = vec![
        zi(0),
        zi(0),
        zi(0),
        zi(0),
        zi(2) * &fd2 * &a[0][0],
        zi(2) * &fd2 * &a[0][1],
        zi(2) * &fd2 * &a[1][1],
        -zi(2) * fd * &aln[0],
        -zi(2) * fd * &aln[1],
        dot(&ln, &aln) - zi(2) * &d * kn,
    ];
    let mut rows: Vec<Vec<Z>> = dom.to_vec();
    for hp in &sh.hp {
        let p = &hp[0..2];
        let ap = vec![&a[0][0] * &p[0] + &a[1][0] * &p[1], &a[0][1] * &p[0] + &a[1][1] * &p[1]];
        rows.push(vec![fd * &ap[0], fd * &ap[1], -dot(&ap, &ln) + &d * &hp[2]]);
    }
    if let Some(con) = to_con(&rows)? {
        cells.push(Cell { num: num_i, den: zi(2) * fd * &d, con });
    }
    assemble_qua_con_cells(cells)
}

fn outward_normal_of(sh: &Shape, j: usize) -> Vec<Z> {
    let e = &sh.ed[j];
    if let Some(n) = &e.n_in {
        return neg(n);
    }
    vec![e.d[1].clone(), -&e.d[0]]
}

/// [a b c] half-planes -> H-form rows [conic6, side]; None when a row is unsatisfiable.
pub fn to_con(rows: &[Vec<Z>]) -> R<Option<Vec<Vec<Z>>>> {
    let mut con = vec![];
    for r in rows {
        if r[0].is_zero() && r[1].is_zero() {
            if r[2].is_negative() {
                return Ok(None);
            }
            continue;
        }
        con.push(con_row(conic_lin(r)?, sgn_of(r, 1)?, 0));
    }
    Ok(Some(con))
}

/// The H-form sign for `row >= 0` (want = +1) or `row <= 0` (want = -1) after conic().
pub fn sgn_of(row: &[Z], want: i64) -> R<i64> {
    if want == 0 {
        return Ok(0);
    }
    let r: Vec<Z> = if row.len() == 3 {
        let mut v = vec![zi(0), zi(0), zi(0)];
        v.extend(row.iter().cloned());
        v
    } else {
        row.to_vec()
    };
    match r.iter().find(|x| !x.is_zero()) {
        None => bail!("ratQ:zeroConic", "a constraint with all-zero coefficients names no half-plane."),
        Some(nz) => Ok(if nz.is_negative() { -want } else { want }),
    }
}

pub fn piece_shape(obj: &QuaPol, k: usize, vi: &[Vec<Z>], vd: &Z) -> R<Shape> {
    let own: Vec<usize> = obj.f.iter().enumerate().filter(|(_, r)| r[0] == k || r[1] == k).map(|(j, _)| j + 1).collect();
    if own.is_empty() {
        bail!("PLQ:conjQ:noFace", "face {} has no edges.", k);
    }
    let mut sh = Shape { vs: vec![], ed: vec![], rays: vec![], hp: vec![], bounded: true, convex: true };
    for &jg in &own {
        let e = obj.e[jg - 1];
        let (a, b) = (e[0], e[1]);
        let d = vec![&vi[b - 1][0] - &vi[a - 1][0], &vi[b - 1][1] - &vi[a - 1][1]];
        if e[2] != 0 {
            sh.vs.extend([a, b]);
            sh.ed.push(Edge { a, b, d, is_ray: false, j: jg, n_in: None });
        } else {
            sh.vs.push(a);
            sh.ed.push(Edge { a, b: 0, d: d.clone(), is_ray: true, j: jg, n_in: None });
            sh.rays.push(d);
        }
    }
    sh.vs = sh.vs.iter().cloned().collect::<BTreeSet<_>>().into_iter().collect();
    sh.bounded = sh.rays.is_empty();
    for ei in 0..sh.ed.len() {
        let ed = sh.ed[ei].clone();
        let p = &vi[ed.a - 1];
        let n = vec![ed.d[1].clone(), -&ed.d[0]];
        let c0 = -dot(&n, p);
        let mut sgn = 0;
        for &v in &sh.vs {
            let v = &vi[v - 1];
            let t = &n[0] * &v[0] + &n[1] * &v[1] + &c0;
            if !t.is_zero() {
                sgn = sign(&t);
                break;
            }
        }
        if sgn == 0 {
            for r in &sh.rays {
                let t = dot(&n, r);
                if !t.is_zero() {
                    sgn = sign(&t);
                    break;
                }
            }
        }
        if sgn == 0 {
            let jg = ed.j;
            if obj.f[jg - 1][1] == k {
                sgn = 1;
            } else if obj.f[jg - 1][0] == k {
                sgn = -1;
            } else {
                bail!("PLQ:conjQ:badSubdivision", "edge {} lists faces [{} {}], neither of which is face {}.", jg, obj.f[jg - 1][0], obj.f[jg - 1][1], k);
            }
        }
        let s = zi(sgn as i64);
        for &v in &sh.vs {
            let v = &vi[v - 1];
            if (&s * (&n[0] * &v[0] + &n[1] * &v[1] + &c0)).is_negative() {
                sh.convex = false;
            }
        }
        for r in &sh.rays {
            if (&s * dot(&n, r)).is_negative() {
                sh.convex = false;
            }
        }
        sh.hp.push(vec![&s * vd * &n[0], &s * vd * &n[1], &s * &c0]);
        sh.ed[ei].n_in = Some(vec![&s * &n[0], &s * &n[1]]);
    }
    Ok(sh)
}

fn dirs_at_vertex(sh: &Shape, i: usize) -> Vec<Vec<Z>> {
    let mut dirs = vec![];
    for e in &sh.ed {
        if e.a == i {
            dirs.push(e.d.clone());
        } else if !e.is_ray && e.b == i {
            dirs.push(neg(&e.d));
        }
    }
    dirs
}

fn sign2(n: &[Vec<Z>]) -> Vec<Vec<Z>> {
    n.iter()
        .map(|r| {
            let g = r[0].gcd(&r[1]);
            if g.is_zero() {
                r.clone()
            } else {
                vec![&r[0] / &g, &r[1] / &g]
            }
        })
        .collect()
}

fn recession_conditions(sh: &Shape, hn: &[Vec<Z>], ln: &[Z], fd: &Z, vi: &[Vec<Z>], vd: &Z) -> R<(bool, Dom, String)> {
    let mut dom = vec![];
    let r = &sh.rays;
    if r.is_empty() {
        return Ok((true, dom, String::new()));
    }
    let fdvd = fd * vd;
    if !sh.hp.is_empty() {
        let np: BTreeSet<Vec<Z>> = sign2(&sh.hp.iter().map(|h| h[0..2].to_vec()).collect::<Vec<_>>()).into_iter().collect();
        if np.len() == 1 {
            let p = np.into_iter().next().unwrap();
            if !is_psd2(hn) {
                return Ok((false, dom, "half-plane of directions with negative curvature".into()));
            }
            if det_exact(hn)?.is_zero() {
                let cand: Vec<Vec<Z>> = if hn.iter().flatten().any(|x| !x.is_zero()) {
                    let d0 = if !hn[0][0].is_zero() { vec![-&hn[0][1], hn[0][0].clone()] } else { zv(&[1, 0]) };
                    vec![d0.clone(), neg(&d0)]
                } else {
                    vec![zv(&[1, 0]), zv(&[-1, 0]), zv(&[0, 1]), zv(&[0, -1])]
                };
                for d in cand {
                    if dot(&p, &d).is_negative() {
                        continue;
                    }
                    let row = vec![-&fdvd * &d[0], -&fdvd * &d[1], vd * dot(ln, &d)];
                    if all_zero(&row) {
                        continue;
                    }
                    dom.push(row);
                }
            }
            return Ok((true, dom, String::new()));
        }
    }
    for d in r {
        if quad(d, hn, d).is_negative() {
            return Ok((false, dom, "recedes where the objective rises".into()));
        }
    }
    if r.len() >= 2 {
        let a = quad(&r[0], hn, &r[0]);
        let b = quad(&r[0], hn, &r[1]);
        let c = quad(&r[1], hn, &r[1]);
        if b.is_negative() && (&b * &b - &a * &c).is_positive() {
            return Ok((false, dom, "negative curvature inside the cone".into()));
        }
    }
    for d in r {
        if !quad(d, hn, d).is_zero() {
            continue;
        }
        let hd = mv(hn, d);
        for q in r {
            if dot(q, &hd).is_negative() {
                return Ok((false, dom, "unbounded slope".into()));
            }
        }
        let mut mu: Option<Z> = None;
        for &v in &sh.vs {
            let t = dot(&vi[v - 1], &hd);
            if mu.as_ref().map_or(true, |m| &t < m) {
                mu = Some(t);
            }
        }
        let mu = mu.unwrap_or_else(Z::zero);
        let row = vec![-&fdvd * &d[0], -&fdvd * &d[1], fd * &mu + vd * dot(ln, d)];
        if all_zero(&row) {
            continue;
        }
        dom.push(row);
    }
    // A NULL DIRECTION STRICTLY INSIDE the cone: H PSD singular with null(H) off the edge rays.
    // The loop above never sees it, and without its row dom f* is the whole plane (review
    // 2026-11-19, finding 1). For such a d, Hd = 0, so the slope test is void and mu = 0, and
    // the row is exactly the null edge ray's. The closed recession cone is the non-negative span
    // of the rays (the half-plane case already returned above, and a proper cone holds at most
    // one of +/-d0; if a null edge ray exists it is that one, already carrying the row).
    if is_psd2(hn) && det_exact(hn)?.is_zero() && hn.iter().flatten().any(|x| !x.is_zero()) {
        let d0 = if !hn[0][0].is_zero() { vec![-&hn[0][1], hn[0][0].clone()] } else { zv(&[1, 0]) };
        for d in [d0.clone(), neg(&d0)] {
            if !in_ray_cone(r, &d) {
                continue;
            }
            if r.iter().any(|rr| quad(rr, hn, rr).is_zero()) {
                break;
            }
            let row = vec![-&fdvd * &d[0], -&fdvd * &d[1], vd * dot(ln, &d)];
            if all_zero(&row) {
                continue;
            }
            dom.push(row);
        }
    }
    Ok((true, dom, String::new()))
}

/// Is `d` in the NON-NEGATIVE SPAN of the rays? In the plane one pair of rays suffices
/// (Caratheodory), so test every pair; a collinear pair is a ray (same direction) or a line
/// (opposite directions).
fn in_ray_cone(r: &[Vec<Z>], d: &[Z]) -> bool {
    let pair = |r1: &Vec<Z>, r2: &Vec<Z>| -> bool {
        let cr = &r1[0] * &r2[1] - &r1[1] * &r2[0];
        if cr.is_zero() {
            if !(&r1[0] * &d[1] - &r1[1] * &d[0]).is_zero() {
                return false;
            }
            if (&r1[0] * &r2[0] + &r1[1] * &r2[1]).is_negative() {
                return true;                            // opposite: the line spans both ways
            }
            return !(&r1[0] * &d[0] + &r1[1] * &d[1]).is_negative();
        }
        let a = &cr * (&r1[0] * &d[1] - &r1[1] * &d[0]);
        let b = cr * (&d[0] * &r2[1] - &d[1] * &r2[0]);
        !a.is_negative() && !b.is_negative()
    };
    for i in 0..r.len() {
        for j in 0..r.len() {
            if pair(&r[i], &r[j]) {
                return true;
            }
        }
    }
    false
}

struct EdC {
    td: Vec<Z>,
    alpha: Z,
    num: Vec<Z>,
    den: Z,
    is_ray: bool,
}

pub fn case_d_boundary_max(fnk: &[Z], fd: &Z, vi: &[Vec<Z>], vd: &Z, sh: &Shape, dom: &[Vec<Z>]) -> R<QuaCon> {
    let hn = hmat(fnk);
    let ln = vec![fnk[7].clone(), fnk[8].clone()];
    let kn = &fnk[9];
    let m = sh.vs.len();
    let fdvd = fd * vd;
    let den_v = zi(2) * fd * vd * vd;
    let mut qv = vec![];
    let mut num_v = vec![];
    for &v in &sh.vs {
        let v = &vi[v - 1];
        let q = qval(v, &hn, &ln, kn, vd);
        let mut num = vec![zi(0); 7];
        num.extend([zi(2) * &fdvd * &v[0], zi(2) * &fdvd * &v[1], -&q]);
        qv.push(q);
        num_v.push(num);
    }
    let mut ed_c: Vec<EdC> = vec![];
    for e in &sh.ed {
        let d = &e.d;
        let alpha = quad(d, &hn, d);
        if !alpha.is_positive() {
            continue;
        }
        let a = &vi[e.a - 1];
        let ga = vplus(&mv(&hn, a), &ln, vd);
        let qa = qval(a, &hn, &ln, kn, vd);
        let td = vec![&fdvd * &d[0], &fdvd * &d[1], -dot(&ga, d)];
        let lin = [zi(2) * &alpha * &fdvd * &a[0], zi(2) * &alpha * &fdvd * &a[1]];
        let num = vec![
            zi(0),
            zi(0),
            zi(0),
            zi(0),
            zi(2) * &td[0] * &td[0],
            zi(2) * &td[0] * &td[1],
            zi(2) * &td[1] * &td[1],
            &lin[0] + zi(2) * &td[0] * &td[2],
            &lin[1] + zi(2) * &td[1] * &td[2],
            &td[2] * &td[2] - &alpha * &qa,
        ];
        let den = zi(2) * &alpha * fd * vd * vd;
        ed_c.push(EdC { td, alpha, num, den, is_ray: e.is_ray });
    }
    struct C0 {
        lin: Vec<Vec<Z>>,
        vi: usize,
        ei: Vec<usize>,
    }
    let mut cells0: Vec<C0> = vec![];
    for i in 0..m {
        let mut rows: Vec<Vec<Z>> = dom.to_vec();
        let v_i = &vi[sh.vs[i] - 1];
        for j in 0..m {
            if j == i {
                continue;
            }
            let v_j = &vi[sh.vs[j] - 1];
            let r = vec![zi(2) * &fdvd * (&v_i[0] - &v_j[0]), zi(2) * &fdvd * (&v_i[1] - &v_j[1]), -(&qv[i] - &qv[j])];
            if all_zero(&r) {
                bail!("PLQ:conjQ:repeatedVertex", "piece vertices {} and {} coincide.", sh.vs[i], sh.vs[j]);
            }
            rows.push(r);
        }
        if feasible2(&rows, true) {
            cells0.push(C0 { lin: rows, vi: i, ei: vec![] });
        }
    }
    for (e, ed) in ed_c.iter().enumerate() {
        let td = &ed.td;
        let ng = neg(td);
        let (regimes, avail): (Vec<Vec<Vec<Z>>>, Vec<bool>) = if ed.is_ray {
            (vec![vec![ng], vec![td.clone()]], vec![false, true])
        } else {
            (
                vec![
                    vec![ng],
                    vec![td.clone(), vec![-&td[0], -&td[1], &ed.alpha - &td[2]]],
                    vec![vec![td[0].clone(), td[1].clone(), &td[2] - &ed.alpha]],
                ],
                vec![false, true, false],
            )
        };
        let mut nxt = vec![];
        for c in &cells0 {
            for r in 0..regimes.len() {
                let mut lin = c.lin.clone();
                lin.extend(regimes[r].iter().cloned());
                if !feasible2(&lin, true) {
                    continue;
                }
                let mut ei = c.ei.clone();
                if avail[r] {
                    ei.push(e);
                }
                nxt.push(C0 { lin, vi: c.vi, ei });
            }
        }
        cells0 = nxt;
    }
    let mut out = vec![];
    for c in &cells0 {
        let base = match to_con(&c.lin)? {
            None => continue,
            Some(b) => b,
        };
        let mut cand: Vec<(&Vec<Z>, &Z)> = vec![(&num_v[c.vi], &den_v)];
        for &e in &c.ei {
            cand.push((&ed_c[e].num, &ed_c[e].den));
        }
        for i in 0..cand.len() {
            let mut con = base.clone();
            for j in 0..cand.len() {
                if j == i {
                    continue;
                }
                let (dn, _) = sub(cand[i].0, cand[i].1, cand[j].0, cand[j].1)?;
                if all_zero(&dn) {
                    continue;
                }
                let rr = face_conic(&dn);
                con.push(con_row(conic(&rr)?, sgn_of(&rr, 1)?, 0));
            }
            out.push(Cell { num: cand[i].0.clone(), den: cand[i].1.clone(), con });
        }
    }
    assemble_qua_con_cells(out)
}

fn boundary_cycle(sh: &Shape, vi: &[Vec<Z>]) -> R<Vec<usize>> {
    let m = sh.ed.len();
    let mut adj: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for e in &sh.ed {
        adj.entry(e.a).or_default().push(e.b);
        adj.entry(e.b).or_default().push(e.a);
    }
    let mut cyc = vec![sh.ed[0].a, sh.ed[0].b];
    for t in 2..m {
        let nxt: Vec<usize> = adj[&cyc[t - 1]].iter().cloned().filter(|v| *v != cyc[t - 2]).collect();
        if nxt.is_empty() {
            bail!("PLQ:conjQ:notSimple", "the face boundary is not a single closed cycle.");
        }
        cyc.push(nxt[0]);
    }
    if cyc.iter().collect::<BTreeSet<_>>().len() != m {
        bail!("PLQ:conjQ:notSimple", "the face boundary revisits a vertex, so it is not a simple polygon.");
    }
    if twice_area(vi, &cyc).is_negative() {
        cyc.reverse();
    }
    Ok(cyc)
}

fn twice_area(vi: &[Vec<Z>], cyc: &[usize]) -> Z {
    let mut a2 = Z::zero();
    let m = cyc.len();
    for i in 0..m {
        let p = &vi[cyc[i] - 1];
        let q = &vi[cyc[(i + 1) % m] - 1];
        a2 += &p[0] * &q[1] - &p[1] * &q[0];
    }
    a2
}

fn ear_clip(vi: &[Vec<Z>], cyc: Vec<usize>) -> R<Vec<Vec<usize>>> {
    let mut cyc = drop_collinear(vi, cyc);
    let m = cyc.len();
    if m < 3 {
        bail!("PLQ:conjQ:degenerate", "the face collapses to fewer than three corners.");
    }
    let mut tris = vec![];
    let mut guard = 0;
    while cyc.len() > 3 {
        guard += 1;
        if guard > m + 2 {
            bail!("PLQ:conjQ:noEar", "no ear was found on a polygon of {} corners.", cyc.len());
        }
        let n = cyc.len();
        let mut clipped = false;
        for i in 0..n {
            let a = cyc[(i + n - 1) % n];
            let b = cyc[i];
            let c = cyc[(i + 1) % n];
            if !is_ear(vi, &cyc, a, b, c) {
                continue;
            }
            tris.push(vec![a, b, c]);
            cyc.remove(i);
            clipped = true;
            break;
        }
        if !clipped {
            bail!("PLQ:conjQ:noEar", "no ear was found, so the boundary is self-intersecting or not counter-clockwise.");
        }
    }
    tris.push(cyc);
    Ok(tris)
}

fn is_ear(vi: &[Vec<Z>], cyc: &[usize], a: usize, b: usize, c: usize) -> bool {
    if !cross2(vi, a, b, c).is_positive() {
        return false;
    }
    for &v in cyc {
        if v == a || v == b || v == c {
            continue;
        }
        if in_triangle(vi, a, b, c, v) {
            return false;
        }
    }
    true
}

fn cross2(vi: &[Vec<Z>], a: usize, b: usize, c: usize) -> Z {
    let (p, q, r) = (&vi[a - 1], &vi[b - 1], &vi[c - 1]);
    (&q[0] - &p[0]) * (&r[1] - &q[1]) - (&q[1] - &p[1]) * (&r[0] - &q[0])
}

fn in_triangle(vi: &[Vec<Z>], a: usize, b: usize, c: usize, v: usize) -> bool {
    !cross2v(vi, a, b, v).is_negative() && !cross2v(vi, b, c, v).is_negative() && !cross2v(vi, c, a, v).is_negative()
}

fn cross2v(vi: &[Vec<Z>], a: usize, b: usize, v: usize) -> Z {
    let (p, q, r) = (&vi[a - 1], &vi[b - 1], &vi[v - 1]);
    (&q[0] - &p[0]) * (&r[1] - &p[1]) - (&q[1] - &p[1]) * (&r[0] - &p[0])
}

fn drop_collinear(vi: &[Vec<Z>], mut cyc: Vec<usize>) -> Vec<usize> {
    let mut changed = true;
    while changed && cyc.len() > 3 {
        changed = false;
        let n = cyc.len();
        for i in 0..n {
            let a = cyc[(i + n - 1) % n];
            let b = cyc[i];
            let c = cyc[(i + 1) % n];
            if cross2(vi, a, b, c).is_zero() {
                cyc.remove(i);
                changed = true;
                break;
            }
        }
    }
    cyc
}

fn triangle_shape(vi: &[Vec<Z>], vd: &Z, tri: &[usize]) -> R<Shape> {
    let mut vs = tri.to_vec();
    vs.sort();
    let mut sh = Shape { vs, ed: vec![], rays: vec![], hp: vec![], bounded: true, convex: true };
    for i in 0..3 {
        let a = tri[i];
        let b = tri[(i + 1) % 3];
        let o = tri[(i + 2) % 3];
        let d = vec![&vi[b - 1][0] - &vi[a - 1][0], &vi[b - 1][1] - &vi[a - 1][1]];
        sh.ed.push(Edge { a, b, d: d.clone(), is_ray: false, j: 0, n_in: None });
        let n = vec![d[1].clone(), -&d[0]];
        let c0 = -dot(&n, &vi[a - 1]);
        let t = &n[0] * &vi[o - 1][0] + &n[1] * &vi[o - 1][1] + &c0;
        if t.is_zero() {
            bail!("PLQ:conjQ:degenerate", "a split triangle has zero area.");
        }
        let s = zi(sign(&t) as i64);
        sh.hp.push(vec![&s * vd * &n[0], &s * vd * &n[1], &s * &c0]);
    }
    Ok(sh)
}

fn conj_unbounded_non_convex(obj: &QuaPol, k: usize) -> R<(QuaCon, Dom)> {
    let (vi, _) = combine_den(&obj.vn, &obj.vd)?;
    let own: Vec<usize> = obj.f.iter().enumerate().filter(|(_, r)| r[0] == k || r[1] == k).map(|(j, _)| j + 1).collect();
    let mut segs = vec![];
    let mut rays = vec![];
    let mut n_in = vec![];
    for &e in &own {
        let (a, b) = (obj.e[e - 1][0], obj.e[e - 1][1]);
        let d = [&vi[b - 1][0] - &vi[a - 1][0], &vi[b - 1][1] - &vi[a - 1][1]];
        let n = if obj.f[e - 1][0] == k { vec![-&d[1], d[0].clone()] } else { vec![d[1].clone(), -&d[0]] };
        n_in.push(n);
        if obj.e[e - 1][2] == 0 {
            rays.push([a, b]);
        } else {
            segs.push([a, b]);
        }
    }
    if rays.len() != 2 {
        bail!("PLQ:conjQ:nonConvexPiece", "face {} is non-convex and unbounded with {} boundary rays; the cover needs exactly two.", k, rays.len());
    }
    let da: Vec<Z> = (0..2).map(|i| &vi[rays[0][1] - 1][i] - &vi[rays[0][0] - 1][i]).collect();
    let db: Vec<Z> = (0..2).map(|i| &vi[rays[1][1] - 1][i] - &vi[rays[1][0] - 1][i]).collect();
    let mut r: Option<Vec<Z>> = None;
    for c in [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2)] {
        let cand = vec![zi(c.0) * &da[0] + zi(c.1) * &db[0], zi(c.0) * &da[1] + zi(c.1) * &db[1]];
        if n_in.iter().all(|n| dot(n, &cand).is_positive()) {
            let g = cand[0].gcd(&cand[1]).max(Z::one());
            r = Some(vec![&cand[0] / &g, &cand[1] / &g]);
            break;
        }
    }
    let r = match r {
        None => bail!("PLQ:conjQ:nonConvexPiece", "face {} is non-convex and unbounded, and no direction between its rays points into it from every edge.", k),
        Some(r) => r,
    };
    let rq = [qz(&r[0]), qz(&r[1])];
    let q = obj.fx(k);
    let half = Q::new(zi(1), zi(2));
    let mut parts = vec![];
    for s in &segs {
        let p = obj.vx(s[0]);
        let w = obj.vx(s[1]);
        let pr = [&p[0] + &rq[0], &p[1] + &rq[1]];
        let wr = [&w[0] + &rq[0], &w[1] + &rq[1]];
        let ip = [(&p[0] + &w[0]) * &half + &rq[0], (&p[1] + &w[1]) * &half + &rq[1]];
        parts.push(one_part(vec![p, w, pr, wr], vec![[1, 2, 1], [1, 3, 0], [2, 4, 0]], ip, &q)?);
    }
    for i in 0..2 {
        let a = obj.vx(rays[i][0]);
        let b = obj.vx(rays[i][1]);
        let ar = [&a[0] + &rq[0], &a[1] + &rq[1]];
        let ip = [(&a[0] + &b[0]) * &half + &rq[0], (&a[1] + &b[1]) * &half + &rq[1]];
        parts.push(one_part(vec![a, b, ar], vec![[1, 2, 0], [1, 3, 0]], ip, &q)?);
    }
    let cs = crate::par::map(&parts, |part| conj_piece_q(part, 1))?;
    let mut dom = vec![];
    let mut gs = Vec::with_capacity(cs.len());
    for (g, d) in cs {
        gs.push(g);
        dom.extend(d);
    }
    let g = crate::par::reduce(gs, |a, b| max_q(&a, &b))?;
    Ok((g, dom))
}

fn one_part(v: Vec<[Q; 2]>, e: Vec<[usize; 3]>, ip: [Q; 2], q: &[Q]) -> R<QuaPol> {
    let mut f = vec![];
    for ed in &e {
        let a = &v[ed[0] - 1];
        let d = [&v[ed[1] - 1][0] - &a[0], &v[ed[1] - 1][1] - &a[1]];
        if (&d[0] * (&ip[1] - &a[1]) - &d[1] * (&ip[0] - &a[0])).is_positive() {
            f.push([1, 0]);
        } else {
            f.push([0, 1]);
        }
    }
    QuaPol::new_q(v, e, vec![q.to_vec()], f)
}

#[allow(clippy::too_many_arguments)]
fn lineality_conj(_obj: &QuaPol, _k: usize, sh: &Shape, hn: &[Vec<Z>], fnk: &[Z], fd: &Z, vi: &[Vec<Z>], vd: &Z) -> R<Option<QuaCon>> {
    if sh.hp.is_empty() || (!sh.rays.is_empty() && sh.rays.len() != 2) {
        return Ok(None);
    }
    let n: Vec<&[Z]> = sh.hp.iter().map(|h| &h[0..2]).collect();
    let g0 = n[0][0].gcd(&n[0][1]).max(Z::one());
    let p = vec![&n[0][0] / &g0, &n[0][1] / &g0];
    for i in 1..n.len() {
        if !(&n[i][0] * &p[1] - &n[i][1] * &p[0]).is_zero() {
            return Ok(None);
        }
    }
    let d = vec![-&p[1], p[0].clone()];
    if mv(hn, &d).iter().any(|x| !x.is_zero()) {
        return Ok(None);
    }
    let v = [Q::new(vi[sh.vs[0] - 1][0].clone(), vd.clone()), Q::new(vi[sh.vs[0] - 1][1].clone(), vd.clone())];
    let mut lo: Option<Q> = None;
    let mut hi: Option<Q> = None;
    for hp in &sh.hp {
        let coef = dot(&hp[0..2], &p);
        let rhs = -(qz(&hp[0]) * &v[0] + qz(&hp[1]) * &v[1] + qz(&hp[2]));
        let b = rhs / qz(&coef);
        if coef.is_positive() {
            lo = Some(match lo {
                None => b,
                Some(l) => l.max(b),
            });
        } else {
            hi = Some(match hi {
                None => b,
                Some(h) => h.min(b),
            });
        }
    }
    let row: Vec<Q> = fnk.iter().map(|x| Q::new(x.clone(), fd.clone())).collect();
    let at = |t: &Q| -> [Q; 2] { [&v[0] + t * qz(&p[0]), &v[1] + t * qz(&p[1])] };
    let pp = match (&lo, &hi) {
        (Some(l), Some(h)) => {
            if l == h {
                QuaPol::new_q(vec![at(l)], vec![], vec![row], vec![])?
            } else {
                QuaPol::new_q(vec![at(l), at(h)], vec![[1, 2, 1]], vec![row], vec![[0, 0]])?
            }
        }
        (Some(l), None) => QuaPol::new_q(vec![at(l), at(&(l + qi(1)))], vec![[1, 2, 0]], vec![row], vec![[0, 0]])?,
        (None, Some(h)) => QuaPol::new_q(vec![at(h), at(&(h - qi(1)))], vec![[1, 2, 0]], vec![row], vec![[0, 0]])?,
        (None, None) => return Ok(None),
    };
    let (gp, _) = conj_piece_q(&pp, 1)?;
    let ln = vec![fnk[7].clone(), fnk[8].clone()];
    let eq = vec![fd * &d[0], fd * &d[1], -dot(&ln, &d)];
    let eq_row = con_row(conic_lin(&eq)?, sgn_of(&eq, 0)?, 0);
    let mut cells = vec![];
    for f in 0..gp.nf {
        let mut con: Vec<Vec<Z>> = gp.fc[f]
            .iter()
            .map(|r| {
                let mut v = gp.ecq[r[0] as usize - 1].clone();
                v.push(zi(r[1]));
                v.push(zi(r[2]));
                v
            })
            .collect();
        con.push(eq_row.clone());
        cells.push(Cell { num: gp.fn_[f].clone(), den: gp.fd[f].clone(), con });
    }
    Ok(Some(assemble_qua_con_cells(cells)?))
}

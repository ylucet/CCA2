//! biconj_q -- the biconjugate f** = cl co f (the closed convex envelope), computed EXACTLY over
//! the rationals (port of biconj_q.py / biconjQ.m). Returns an AlgCon: a QuaCon, an AlgCrv or an
//! EnvelopeMin.

use crate::alg_con::AlgCon;
use crate::bail;
use crate::cells::{assemble_qua_con_cells, con_row, Cell};
use crate::conj_q::{assert_no_closed_edges, conj_q};
use crate::envelope::{envelope_assemble, envelope_prune, Env};
use crate::errors::{CErr, R};
use crate::is_convex_q::is_convex_q;
use crate::qua_pol::QuaPol;
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Debug, Default)]
pub struct Opts {
    pub overlay_budget: Option<f64>,
}

pub fn biconj_q(obj: &QuaPol, opts: &Opts) -> R<AlgCon> {
    assert_no_closed_edges(&obj.e, "biconjQ")?;
    obj.assert_operable()?;
    obj.assert_exact()?;
    obj.assert_continuous()?;
    let nf = obj.fn_.len();

    if obj.dom.dim < 2 {
        if is_chain_subdivision(obj) {
            assert_chain_edge_concave(obj)?;
            if is_chain(obj)? {
                return Ok(AlgCon::QuaCon(concave_envelope(obj)?));
            }
            return Ok(AlgCon::QuaCon(collinear_chain_envelope(obj)?));
        }
        return Ok(AlgCon::QuaCon(thin_envelope(obj)?));
    }

    let hof = |k: usize| -> Vec<Vec<Z>> {
        let r = &obj.fn_[k];
        vec![vec![r[4].clone(), r[5].clone()], vec![r[5].clone(), r[6].clone()]]
    };
    let all_psd = (0..nf).all(|k| is_psd2(&hof(k)));
    let flagged = obj.f_is_convex == Some(true);
    if flagged {
        assert_flag_not_refuted(obj)?;
    }
    if all_psd && (nf == 1 || flagged || is_convex_q(obj)?.0) {
        if obj.nv == 0 {
            return Ok(AlgCon::QuaCon(assemble_qua_con_cells(vec![Cell { num: obj.fn_[0].clone(), den: obj.fd[0].clone(), con: vec![] }])?));
        }
        return Ok(AlgCon::QuaCon(qua_pol_as_qua_con(obj)?));
    }
    if obj.f_is_convex == Some(true) && !all_psd {
        bail!("PLQ:biconjQ:notConvex", "fIsConvex is asserted true but some piece has a Hessian that is not positive semidefinite.");
    }

    let mut assembly_why = String::new();
    if all_psd && nf == 2 {
        match assembled_envelope(obj) {
            Ok(h) => return Ok(h),
            Err(me) => {
                if !is_declined(&me) {
                    return Err(me);
                }
                assembly_why = format!(" The assembly route declined first: {} -- {}", me.id, me.msg.chars().take(200).collect::<String>());
            }
        }
    }

    let (vi, vd) = combine_den(&obj.vn, &obj.vd)?;
    for k in 1..=nf {
        let hk = hof(k - 1);
        let (vsk, rays, bnd) = piece_geometry_b(obj, k)?;
        if !bnd {
            let lk = vec![obj.fn_[k - 1][7].clone(), obj.fn_[k - 1][8].clone()];
            let mut diverges = false;
            'outer: for d in &rays {
                let a = quad(d, &hk, d);
                if a.is_negative() {
                    diverges = true;
                    break;
                }
                if a.is_zero() {
                    for &v in &vsk {
                        let v = &vi[v - 1];
                        let g: Vec<Z> = mv(&hk, v).iter().zip(lk.iter()).map(|(x, y)| x + &vd * y).collect();
                        if dot(&g, d).is_negative() {
                            diverges = true;
                            break 'outer;
                        }
                    }
                }
            }
            if !diverges && rays.len() >= 2 {
                let a = quad(&rays[0], &hk, &rays[0]);
                let b = quad(&rays[0], &hk, &rays[1]);
                let c = quad(&rays[1], &hk, &rays[1]);
                if b.is_negative() && (&b * &b - &a * &c).is_positive() {
                    diverges = true;
                }
            }
            if diverges {
                bail!("PLQ:biconjQ:minusInfinity", "piece {} is unbounded and q runs DOWNWARD along one of its recession directions, so co f is -infinity there.", k);
            }
        }
        let own: Vec<usize> = obj.f.iter().enumerate().filter(|(_, r)| r[0] == k || r[1] == k).map(|(j, _)| j + 1).collect();
        for j in own {
            let (a, b) = (obj.e[j - 1][0], obj.e[j - 1][1]);
            let va = obj.vx(a);
            let vb = obj.vx(b);
            let (dn, _) = from_double(&[Num::Q(&vb[0] - &va[0]), Num::Q(&vb[1] - &va[1])])?;
            if quad(&dn, &hk, &dn).is_positive() {
                if all_psd {
                    bail!("PLQ:biconjQ:convexPieceCoupling", "every piece is convex but f is not, so co f must BRIDGE the gradient jump across a shared edge (piece {} curves upward along one of its own edges).{}", k, assembly_why);
                }
                return atoms_route(obj, opts);
            }
        }
        if bnd {
            if vsk.len() < 3 {
                bail!("PLQ:biconjQ:degenerate", "piece {} has {} vertices; a 2-D bounded face needs at least three.", k, vsk.len());
            }
        } else if vsk.len() + rays.len() < 3 {
            bail!("PLQ:biconjQ:degenerate", "piece {} has {} vertices and {} rays; a 2-D face needs three generators.", k, vsk.len(), rays.len());
        }
    }
    Ok(AlgCon::QuaCon(concave_envelope(obj)?))
}

fn atoms_route(obj: &QuaPol, opts: &Opts) -> R<AlgCon> {
    crate::biconj_atoms::atoms_route(obj, opts)
}

pub fn assembled_envelope(obj: &QuaPol) -> R<AlgCon> {
    let pieces = split_pieces(obj)?;
    let (ga, _) = conj_q(&pieces[0])?;
    let (gb, _) = conj_q(&pieces[1])?;
    let env = envelope_prune(&envelope_assemble(&ga, &gb, &[])?.0)?;
    subdivide_pair(env, &hull_rows(obj)?)
}

pub fn subdivide_pair(env: Env, dom: &[Vec<Z>]) -> R<AlgCon> {
    crate::envelope_subdivision::subdivide_pair(env, dom)
}

/// MATLAB convhull(x, y): indices (1-based) of the hull, counter-clockwise, CLOSED.
pub fn convhull2(pts: &[[Q; 2]]) -> Vec<usize> {
    let n = pts.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| (&pts[i][0], &pts[i][1]).cmp(&(&pts[j][0], &pts[j][1])));
    let cross = |o: usize, a: usize, b: usize| -> Q { (&pts[a][0] - &pts[o][0]) * (&pts[b][1] - &pts[o][1]) - (&pts[a][1] - &pts[o][1]) * (&pts[b][0] - &pts[o][0]) };
    let mut lower: Vec<usize> = vec![];
    for &i in &order {
        while lower.len() >= 2 && !cross(lower[lower.len() - 2], lower[lower.len() - 1], i).is_positive() {
            lower.pop();
        }
        lower.push(i);
    }
    let mut upper: Vec<usize> = vec![];
    for &i in order.iter().rev() {
        while upper.len() >= 2 && !cross(upper[upper.len() - 2], upper[upper.len() - 1], i).is_positive() {
            upper.pop();
        }
        upper.push(i);
    }
    let mut hull: Vec<usize> = lower[..lower.len() - 1].to_vec();
    hull.extend_from_slice(&upper[..upper.len() - 1]);
    let mut out: Vec<usize> = hull.iter().map(|i| i + 1).collect();
    out.push(hull[0] + 1);
    out
}

/// The half-planes of conv(dom f) for a BOUNDED domain; empty when the domain is unbounded.
pub fn hull_rows(obj: &QuaPol) -> R<Vec<Vec<Z>>> {
    let mut r = vec![];
    if obj.e.iter().any(|e| e[2] == 0) {
        return Ok(r);
    }
    let v: Vec<[Q; 2]> = (1..=obj.vn.len()).map(|i| obj.vx(i)).collect();
    let k = convhull2(&v);
    for i in 0..k.len() - 1 {
        let (a, b) = (k[i], k[i + 1]);
        let (pn, pd) = (&obj.vn[a - 1], &obj.vd[a - 1]);
        let (qn, qd) = (&obj.vn[b - 1], &obj.vd[b - 1]);
        let n = vec![-(&qn[1] * pd - &pn[1] * qd), &qn[0] * pd - &pn[0] * qd];
        let row = vec![&n[0] * pd, &n[1] * pd, -dot(&n, pn)];
        let g = row[0].gcd(&row[1]).gcd(&row[2]);
        r.push(row.iter().map(|x| x.div_floor(&g)).collect());
    }
    Ok(r)
}

pub fn is_declined(me: &CErr) -> bool {
    ["envelope", "bridgeFace", "bridgeArcCut", "cornerFace", "triangleFace", "cutTriangles"].iter().any(|p| me.id.starts_with(p))
}

/// The two pieces of a two-face QuaPol as one-face QuaPols.
pub fn split_pieces(obj: &QuaPol) -> R<Vec<QuaPol>> {
    let mut out = vec![];
    for k in 1..=2usize {
        let own: Vec<usize> = obj.f.iter().enumerate().filter(|(_, r)| r[0] == k || r[1] == k).map(|(j, _)| j).collect();
        let vs: Vec<usize> = own.iter().flat_map(|&j| [obj.e[j][0], obj.e[j][1]]).filter(|v| *v > 0).collect::<BTreeSet<_>>().into_iter().collect();
        let idx: BTreeMap<usize, usize> = vs.iter().enumerate().map(|(i, v)| (*v, i + 1)).collect();
        let v: Vec<[Q; 2]> = vs.iter().map(|&x| obj.vx(x)).collect();
        let e: Vec<[usize; 3]> = own.iter().map(|&j| [idx[&obj.e[j][0]], idx[&obj.e[j][1]], obj.e[j][2]]).collect();
        let f: Vec<[usize; 2]> = own.iter().map(|&j| [(obj.f[j][0] == k) as usize, (obj.f[j][1] == k) as usize]).collect();
        out.push(QuaPol::new_q(v, e, vec![obj.fx(k)], f)?);
    }
    Ok(out)
}

fn assert_flag_not_refuted(obj: &QuaPol) -> R<()> {
    let lvl: i64 = std::env::var("BICONJ_ASSERT").ok().and_then(|s| s.parse().ok()).unwrap_or(0);
    if lvl < 1 {
        return Ok(());
    }
    let (tf, why) = is_convex_q(obj)?;
    if !tf {
        bail!("PLQ:biconjQ:notConvex", "fIsConvex is asserted true but the exact test refutes it: {}.", why);
    }
    Ok(())
}

fn hk_of(fnk: &[Z]) -> (Vec<Vec<Z>>, Vec<Z>, Z) {
    (vec![vec![fnk[4].clone(), fnk[5].clone()], vec![fnk[5].clone(), fnk[6].clone()]], vec![fnk[7].clone(), fnk[8].clone()], fnk[9].clone())
}

fn qv_of(v: &[Z], hn: &[Vec<Z>], ln: &[Z], kn: &Z, vd: &Z) -> Z {
    quad(v, hn, v) + zi(2) * vd * dot(ln, v) + zi(2) * vd * vd * kn
}

/// co f for EDGE-CONCAVE pieces: the lower convex hull of the lifted vertices and rays.
pub fn concave_envelope(obj: &QuaPol) -> R<crate::qua_con::QuaCon> {
    let (vi, vd) = combine_den(&obj.vn, &obj.vd)?;
    let nf = obj.fn_.len();
    let chain = is_chain(obj)?;
    let mut lden = Z::one();
    for k in 0..nf {
        lden = lcm(&lden, &obj.fd[k]);
    }
    let dl = zi(2) * &lden * &vd * &vd;
    let mut lifted: BTreeMap<usize, Z> = BTreeMap::new();
    let mut r3: Vec<Vec<Z>> = vec![];
    for k in 1..=nf {
        let (fnk, fdk) = obj.face_q(k);
        let (hn, ln, kn) = hk_of(&fnk);
        let scale = &lden / &fdk;
        let (vsk, raysk): (Vec<usize>, Vec<Vec<Z>>) = if chain {
            let e = obj.e[k - 1];
            if e[2] != 0 {
                (vec![e[0], e[1]], vec![])
            } else {
                (vec![e[0]], vec![vec![&vi[e[1] - 1][0] - &vi[e[0] - 1][0], &vi[e[1] - 1][1] - &vi[e[0] - 1][1]]])
            }
        } else {
            let (a, b, _) = piece_geometry_b(obj, k)?;
            (a, b)
        };
        for d in &raysk {
            let mut rho: Option<Z> = None;
            for &v in &vsk {
                let v = &vi[v - 1];
                let g: Vec<Z> = mv(&hn, v).iter().zip(ln.iter()).map(|(x, y)| x + &vd * y).collect();
                let ri = dot(&g, d);
                if rho.as_ref().map_or(true, |r| &ri < r) {
                    rho = Some(ri);
                }
            }
            let rho = rho.unwrap_or_else(Z::zero);
            r3.push(vec![zi(2) * &lden * &vd * &d[0], zi(2) * &lden * &vd * &d[1], zi(2) * &scale * rho]);
        }
        for &v in &vsk {
            let qv = &scale * qv_of(&vi[v - 1], &hn, &ln, &kn, &vd);
            let e = lifted.entry(v).or_insert_with(|| qv.clone());
            if qv < *e {
                *e = qv;
            }
        }
    }
    let vs: Vec<usize> = lifted.keys().cloned().collect();
    let m = vs.len();
    let r3: Vec<Vec<Z>> = r3.into_iter().collect::<BTreeSet<_>>().into_iter().collect();
    let nr = r3.len();
    if m + nr < 3 || m < 1 {
        bail!("PLQ:biconjQ:degenerate", "the domain has {} vertices and {} recession rays; a 2-D one needs three generators.", m, nr);
    }
    let q3: Vec<Vec<Z>> = vs.iter().map(|v| vec![zi(2) * &lden * &vd * &vi[v - 1][0], zi(2) * &lden * &vd * &vi[v - 1][1], lifted[v].clone()]).collect();
    let mut g: Vec<Vec<Z>> = q3.clone();
    g.extend(r3.iter().cloned());
    let is_p: Vec<bool> = (0..m + nr).map(|i| i < m).collect();
    let ng = m + nr;
    let mut facets: BTreeSet<[Z; 4]> = BTreeSet::new();
    for i in 0..ng {
        for j in i + 1..ng {
            for k in j + 1..ng {
                let idx = [i, j, k];
                let pts: Vec<usize> = idx.iter().cloned().filter(|t| is_p[*t]).collect();
                if pts.is_empty() {
                    continue;
                }
                let base = pts[0];
                let mut dirs = vec![];
                for &t in &idx {
                    if t == base {
                        continue;
                    }
                    if is_p[t] {
                        dirs.push((0..3).map(|c| &g[t][c] - &g[base][c]).collect::<Vec<Z>>());
                    } else {
                        dirs.push(g[t].clone());
                    }
                }
                let (u, w) = (&dirs[0], &dirs[1]);
                let mut n = vec![&u[1] * &w[2] - &u[2] * &w[1], &u[2] * &w[0] - &u[0] * &w[2], &u[0] * &w[1] - &u[1] * &w[0]];
                if n[2].is_zero() {
                    continue;
                }
                if n[2].is_negative() {
                    n = neg(&n);
                }
                let d = dot(&n, &g[base]);
                let mut below = q3.iter().any(|q| (dot(&n, q) - &d).is_negative());
                if !below {
                    below = r3.iter().any(|r| dot(&n, r).is_negative());
                }
                if below {
                    continue;
                }
                let mut gg = n[0].gcd(&n[1]).gcd(&n[2].gcd(&d));
                if gg.is_zero() {
                    gg = Z::one();
                }
                facets.insert([n[0].div_floor(&gg), n[1].div_floor(&gg), n[2].div_floor(&gg), d.div_floor(&gg)]);
            }
        }
    }
    let facets: Vec<[Z; 4]> = facets.into_iter().collect();
    if facets.is_empty() {
        bail!("PLQ:biconjQ:noFacet", "no lower-hull facet was found, which cannot happen for three non-collinear vertices.");
    }
    let nfac = facets.len();
    let r3d: Vec<Vec<Z>> = r3.iter().map(|r| r[0..2].to_vec()).collect();
    let dom = hull_half_planes(&vi, &vd, &vs, &r3d);
    let face = |t: usize| -> (Vec<Z>, Z) {
        let [a, b, c, d] = &facets[t];
        let mut num = vec![zi(0); 7];
        num.extend([-a * &dl, -b * &dl, d.clone()]);
        (num, &dl * c)
    };
    let mut cells = vec![];
    for t in 0..nfac {
        let (num, den) = face(t);
        let mut rows: Vec<Vec<Z>> = dom.clone();
        for s in 0..nfac {
            if s == t {
                continue;
            }
            let (n2, d2) = face(s);
            let (dn, _) = sub(&num, &den, &n2, &d2)?;
            if all_zero(&dn) {
                continue;
            }
            rows.push(vec![zi(2) * &dn[7], zi(2) * &dn[8], zi(2) * &dn[9]]);
        }
        if !feasible2(&rows, true) {
            continue;
        }
        let mut con = vec![];
        for r in &rows {
            con.push(con_row(conic_lin(r)?, sgn_of_b(r, 1)?, 0));
        }
        cells.push(Cell { num, den, con });
    }
    assemble_qua_con_cells(cells)
}

/// The half-planes of conv(vertices) + cone(dirs), as primitive integer rows.
pub fn hull_half_planes(vi: &[Vec<Z>], vd: &Z, vs: &[usize], dirs: &[Vec<Z>]) -> Vec<Vec<Z>> {
    let mut rows: Vec<Vec<Z>> = vec![];
    let m = vs.len();
    let r = dirs.len();
    for a in 0..m {
        for b in 0..m + r {
            let p = &vi[vs[a] - 1];
            let n: Vec<Z> = if b < m {
                if a == b {
                    continue;
                }
                let q = &vi[vs[b] - 1];
                vec![&q[1] - &p[1], -(&q[0] - &p[0])]
            } else {
                let e = &dirs[b - m];
                vec![e[1].clone(), -&e[0]]
            };
            if n[0].is_zero() && n[1].is_zero() {
                continue;
            }
            let c0 = -dot(&n, p);
            let mut sgn = 0;
            let mut ok = true;
            for i in 0..m {
                let v = &vi[vs[i] - 1];
                let t = &n[0] * &v[0] + &n[1] * &v[1] + &c0;
                if t.is_zero() {
                    continue;
                }
                if sgn == 0 {
                    sgn = sign(&t);
                } else if sign(&t) != sgn {
                    ok = false;
                    break;
                }
            }
            if ok {
                for i in 0..r {
                    let t = dot(&n, &dirs[i]);
                    if t.is_zero() {
                        continue;
                    }
                    if sgn == 0 {
                        sgn = sign(&t);
                    } else if sign(&t) != sgn {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok || sgn == 0 {
                continue;
            }
            let s = zi(sgn as i64);
            let row = vec![&s * vd * &n[0], &s * vd * &n[1], &s * &c0];
            let g = gcd_row(&row);
            let row: Vec<Z> = row.iter().map(|x| x.div_floor(&g)).collect();
            if !rows.contains(&row) {
                rows.push(row);
            }
        }
    }
    rows
}

pub fn gcd_row(row: &[Z]) -> Z {
    let g = gcd_all(row);
    if g.is_zero() {
        Z::one()
    } else {
        g.abs()
    }
}

fn polygon_half_planes(obj: &QuaPol, vi: &[Vec<Z>], vd: &Z, vs: &[usize], k: usize, rays: &[Vec<Z>]) -> Vec<Vec<Z>> {
    let mut rows = vec![];
    for j in (0..obj.f.len()).filter(|&j| obj.f[j][0] == k || obj.f[j][1] == k) {
        let p = &vi[obj.e[j][0] - 1];
        let q = &vi[obj.e[j][1] - 1];
        let n = vec![&q[1] - &p[1], -(&q[0] - &p[0])];
        let c0 = -dot(&n, p);
        let mut sgn = 0;
        for &v in vs {
            let v = &vi[v - 1];
            let t = &n[0] * &v[0] + &n[1] * &v[1] + &c0;
            if !t.is_zero() {
                sgn = sign(&t);
                break;
            }
        }
        if sgn == 0 {
            for ray in rays {
                let t = dot(&n, ray);
                if !t.is_zero() {
                    sgn = sign(&t);
                    break;
                }
            }
        }
        if sgn == 0 {
            if obj.f[j][1] == k {
                sgn = 1;
            } else if obj.f[j][0] == k {
                sgn = -1;
            } else {
                continue;
            }
        }
        let s = zi(sgn as i64);
        rows.push(vec![&s * vd * &n[0], &s * vd * &n[1], &s * &c0]);
    }
    rows
}

pub fn sgn_of_b(row: &[Z], want: i64) -> R<i64> {
    let nz = row.iter().find(|x| !x.is_zero());
    match nz {
        None => bail!("ratQ:zeroConic", "a constraint with all-zero coefficients names no half-plane."),
        Some(x) => Ok(if x.is_negative() { -want } else { want }),
    }
}

/// f is convex, so co f = f: the same faces, as a QuaCon. A NON-CONVEX piece is refused.
fn qua_pol_as_qua_con(obj: &QuaPol) -> R<crate::qua_con::QuaCon> {
    let (vi, vd) = combine_den(&obj.vn, &obj.vd)?;
    let mut cells = vec![];
    for k in 1..=obj.fn_.len() {
        let (vs, rays, _) = piece_geometry_b(obj, k)?;
        let rows = polygon_half_planes(obj, &vi, &vd, &vs, k, &rays);
        for r in &rows {
            for &v in &vs {
                let v = &vi[v - 1];
                if (&r[0] * &v[0] + &r[1] * &v[1] + &vd * &r[2]).is_negative() {
                    bail!("PLQ:biconjQ:nonConvexPiece", "piece {} is NOT CONVEX, and co f = f here, so the answer must carry that same non-convex region.", k);
                }
            }
            for q in &rays {
                if (&r[0] * &q[0] + &r[1] * &q[1]).is_negative() {
                    bail!("PLQ:biconjQ:nonConvexPiece", "piece {} recedes along ({},{}), leaving one of its own half-planes.", k, q[0], q[1]);
                }
            }
        }
        let mut con = vec![];
        for r in &rows {
            con.push(con_row(conic_lin(r)?, sgn_of_b(r, 1)?, 0));
        }
        cells.push(Cell { num: obj.fn_[k - 1].clone(), den: obj.fd[k - 1].clone(), con });
    }
    assemble_qua_con_cells(cells)
}

fn all_f_zero(obj: &QuaPol) -> bool {
    !obj.f.is_empty() && obj.f.iter().all(|r| r[0] == 0 && r[1] == 0)
}

pub fn is_chain_subdivision(obj: &QuaPol) -> bool {
    all_f_zero(obj) && obj.e.len() >= 2 && obj.nv >= 3
}

fn collinear_exact(obj: &QuaPol) -> R<bool> {
    let (vi, _) = combine_den(&obj.vn, &obj.vd)?;
    let a = &vi[0];
    let mut d = None;
    for p in &vi[1..] {
        let c = [&p[0] - &a[0], &p[1] - &a[1]];
        if !c[0].is_zero() || !c[1].is_zero() {
            d = Some(c);
            break;
        }
    }
    let d = match d {
        None => return Ok(true),
        Some(d) => d,
    };
    Ok(vi.iter().all(|p| (&d[0] * (&p[1] - &a[1]) - &d[1] * (&p[0] - &a[0])).is_zero()))
}

pub fn is_chain(obj: &QuaPol) -> R<bool> {
    if !all_f_zero(obj) || obj.e.len() < 2 || obj.nv < 3 {
        return Ok(false);
    }
    Ok(!collinear_exact(obj)?)
}

fn assert_chain_edge_concave(obj: &QuaPol) -> R<()> {
    let (vi, _) = combine_den(&obj.vn, &obj.vd)?;
    for j in 1..=obj.e.len() {
        let (fnk, _) = obj.face_q(j);
        let (hn, _, _) = hk_of(&fnk);
        let e = obj.e[j - 1];
        let d = vec![&vi[e[1] - 1][0] - &vi[e[0] - 1][0], &vi[e[1] - 1][1] - &vi[e[0] - 1][1]];
        if quad(&d, &hn, &d).is_positive() {
            bail!("PLQ:biconjQ:notImplemented", "chain edge {} curves UPWARD along itself, so its lifted arc dips below its own chord.", j);
        }
    }
    Ok(())
}

/// (vs, rays, bounded): the vertices (1-based, sorted) and integer ray directions of face k.
pub fn piece_geometry_b(obj: &QuaPol, k: usize) -> R<(Vec<usize>, Vec<Vec<Z>>, bool)> {
    let own: Vec<usize> = (0..obj.f.len()).filter(|&j| obj.f[j][0] == k || obj.f[j][1] == k).collect();
    if own.is_empty() {
        bail!("PLQ:biconjQ:noFace", "face {} has no edges.", k);
    }
    let (vi, _) = combine_den(&obj.vn, &obj.vd)?;
    let mut vs = vec![];
    let mut rays = vec![];
    for j in own {
        let e = obj.e[j];
        if e[2] != 0 {
            vs.push(e[0]);
            vs.push(e[1]);
        } else {
            vs.push(e[0]);
            rays.push(vec![&vi[e[1] - 1][0] - &vi[e[0] - 1][0], &vi[e[1] - 1][1] - &vi[e[0] - 1][1]]);
        }
    }
    let vs: Vec<usize> = vs.into_iter().collect::<BTreeSet<_>>().into_iter().collect();
    let b = rays.is_empty();
    Ok((vs, rays, b))
}

fn collinear_chain_envelope(obj: &QuaPol) -> R<crate::qua_con::QuaCon> {
    let (vi, vd) = combine_den(&obj.vn, &obj.vd)?;
    let ne = obj.e.len();
    let a = vi[0].clone();
    let mut d = None;
    for p in &vi[1..] {
        let c = vec![&p[0] - &a[0], &p[1] - &a[1]];
        if !c[0].is_zero() || !c[1].is_zero() {
            d = Some(c);
            break;
        }
    }
    let mut d = match d {
        None => bail!("PLQ:biconjQ:degenerate", "every vertex of the chain is the same point."),
        Some(d) => d,
    };
    let g = d[0].gcd(&d[1]);
    if g > Z::one() {
        d = vec![&d[0] / &g, &d[1] / &g];
    }
    let mut lden = Z::one();
    for j in 0..ne {
        lden = lcm(&lden, &obj.fd[j]);
    }
    let dl = zi(2) * &lden * &vd * &vd;
    let mut lifted: BTreeMap<usize, Z> = BTreeMap::new();
    for j in 1..=ne {
        let (fnk, fdk) = obj.face_q(j);
        let (hn, ln, kn) = hk_of(&fnk);
        let scale = &lden / &fdk;
        for e in 0..2 {
            let iv = obj.e[j - 1][e];
            let qv = &scale * qv_of(&vi[iv - 1], &hn, &ln, &kn, &vd);
            let en = lifted.entry(iv).or_insert_with(|| qv.clone());
            if qv < *en {
                *en = qv;
            }
        }
    }
    let vs: Vec<usize> = lifted.keys().cloned().collect();
    let mut tz: Vec<(Z, Z)> = vs.iter().map(|v| (dot(&[&vi[v - 1][0] - &a[0], &vi[v - 1][1] - &a[1]], &d), lifted[v].clone())).collect();
    tz.sort_by(|x, y| x.0.cmp(&y.0));
    let t: Vec<Z> = tz.iter().map(|x| x.0.clone()).collect();
    let z: Vec<Z> = tz.iter().map(|x| x.1.clone()).collect();
    let mut st: Vec<usize> = vec![];
    for i in 0..t.len() {
        while st.len() >= 2 {
            let o = st[st.len() - 2];
            let a1 = st[st.len() - 1];
            let cr = (&t[a1] - &t[o]) * (&z[i] - &z[o]) - (&z[a1] - &z[o]) * (&t[i] - &t[o]);
            if !cr.is_positive() {
                st.pop();
            } else {
                break;
            }
        }
        st.push(i);
    }
    if st.len() < 2 {
        bail!("PLQ:biconjQ:degenerate", "the collinear chain produced fewer than two hull points.");
    }
    let n0 = vec![d[1].clone(), -&d[0]];
    let c0 = -dot(&n0, &a);
    let on_line = vec![&vd * &n0[0], &vd * &n0[1], c0];
    let ad = dot(&a, &d);
    let mut cells = vec![];
    for k in 0..st.len() - 1 {
        let (p, q) = (st[k], st[k + 1]);
        let dt = &t[q] - &t[p];
        let dz = &z[q] - &z[p];
        let mut num = vec![zi(0); 7];
        num.extend([&dz * &vd * &d[0], &dz * &vd * &d[1], &z[p] * &dt - &dz * (&ad + &t[p])]);
        let den = &dl * &dt;
        let lo = vec![&vd * &d[0], &vd * &d[1], -(&ad + &t[p])];
        let hi = vec![-&vd * &d[0], -&vd * &d[1], &ad + &t[q]];
        let con = vec![con_row(conic_lin(&on_line)?, 0, 0), con_row(conic_lin(&lo)?, sgn_of_b(&lo, 1)?, 0), con_row(conic_lin(&hi)?, sgn_of_b(&hi, 1)?, 0)];
        cells.push(Cell { num, den, con });
    }
    assemble_qua_con_cells(cells)
}

fn thin_envelope(obj: &QuaPol) -> R<crate::qua_con::QuaCon> {
    let (fnk, fdk) = obj.face_q(1);
    let (hn, ln, kn) = hk_of(&fnk);
    let (vi, vd) = combine_den(&obj.vn, &obj.vd)?;
    if obj.e.is_empty() {
        let p = &vi[0];
        let qp = qv_of(p, &hn, &ln, &kn, &vd);
        let r1 = vec![vd.clone(), zi(0), -&p[0]];
        let r2 = vec![zi(0), vd.clone(), -&p[1]];
        let con = vec![con_row(conic_lin(&r1)?, 0, 0), con_row(conic_lin(&r2)?, 0, 0)];
        let mut num = vec![zi(0); 9];
        num.push(qp);
        return assemble_qua_con_cells(vec![Cell { num, den: zi(2) * &fdk * &vd * &vd, con }]);
    }
    if obj.e.len() != 1 {
        bail!("PLQ:biconjQ:chainNotImplemented", "a chain of {} edges reached thinEnvelope.", obj.e.len());
    }
    if obj.e[0][2] == 0 {
        bail!("PLQ:biconjQ:unbounded", "the domain is a RAY.");
    }
    let a = &vi[obj.e[0][0] - 1];
    let b = &vi[obj.e[0][1] - 1];
    let d = vec![&b[0] - &a[0], &b[1] - &a[1]];
    let alpha = quad(&d, &hn, &d);
    let n = vec![d[1].clone(), -&d[0]];
    let c0 = -dot(&n, a);
    let rows = vec![&vd * &n[0], &vd * &n[1], c0];
    let lo = vec![&vd * &d[0], &vd * &d[1], -dot(&d, a)];
    let hi = vec![-&vd * &d[0], -&vd * &d[1], dot(&d, b)];
    let con = vec![con_row(conic_lin(&rows)?, 0, 0), con_row(conic_lin(&lo)?, sgn_of_b(&lo, 1)?, 0), con_row(conic_lin(&hi)?, sgn_of_b(&hi, 1)?, 0)];
    let (num, den) = if !alpha.is_negative() {
        (fnk.clone(), fdk.clone())
    } else {
        let qa = qv_of(a, &hn, &ln, &kn, &vd);
        let qb = qv_of(b, &hn, &ln, &kn, &vd);
        let dd = dot(&d, &d);
        let mut num = vec![zi(0); 7];
        num.extend([(&qb - &qa) * &vd * &d[0], (&qb - &qa) * &vd * &d[1], &qa * &dd - (&qb - &qa) * dot(a, &d)]);
        (num, zi(2) * &fdk * &vd * &vd * &dd)
    };
    assemble_qua_con_cells(vec![Cell { num, den, con }])
}

//! Corner faces (ports of corner_face_alg.py and corner_wedge_rows.py).

use crate::alg_branch::alg_branch_info;
use crate::bail;
use crate::envelope::{BFace, FaceKind, Sigma, Src};
use crate::errors::R;
use crate::isolated_root_sign::isolated_root_sign;
use crate::poly_n::{self as pn, PolyN};
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::rc::Rc;

pub fn corner_face_alg(n: &[Z], d: &Z, p: &[Vec<Z>], tie: &[Z], side: i64) -> R<Vec<BFace>> {
    if n.len() != 10 || !all_zero(&n[0..4]) {
        bail!("cornerFaceAlg:shape", "n must be a 1x10 QUADRATIC face numerator.");
    }
    let nn = vec![vec![n[4].clone(), n[5].clone()], vec![n[5].clone(), n[6].clone()]];
    let m = vec![n[7].clone(), n[8].clone()];
    let qq = &n[9];
    let nf: Vec<Vec<f64>> = nn.iter().map(|r| r.iter().map(to_f64).collect()).collect();
    let mf = [to_f64(&m[0]), to_f64(&m[1])];
    let qqf = to_f64(qq);
    let df = to_f64(d);
    let tf: Vec<f64> = tie.iter().map(to_f64).collect();
    let mut faces = vec![];
    for (i, row) in p.iter().enumerate() {
        let a = vec![row[0].clone(), row[1].clone()];
        let c = &row[2];
        let na = &a[0] * &a[0] + &a[1] * &a[1];
        if na.is_zero() {
            continue;
        }
        let w = vec![a[1].clone(), -&a[0]];
        let s0n = vec![-c * &a[0], -c * &a[1]];
        let (al, be, ga) = restrict_tie(tie, &s0n, &w, &na);
        if al.is_zero() {
            if be.is_zero() {
                continue;
            }
            faces.extend(rational_corner(n, d, &nn, &m, qq, p, &a, &s0n, &w, &na, &be, &ga, tie, side, i + 1)?);
            continue;
        }
        let g = al.gcd(&be).gcd(&ga);
        let (al, be, ga) = (al.div_floor(&g), be.div_floor(&g), ga.div_floor(&g));
        let dd = &be * &be - zi(4) * &al * &ga;
        if !dd.is_positive() {
            continue;
        }
        let (e0, e1, k) = value_form(&nn, &m, qq, d, &s0n, &w, &na, &al, &be, &ga);
        let fp = Rc::new(implicit_value(&e0, &e1, &k, &al, &be, &dd));
        let (alf, bef, ddf, naf) = (to_f64(&al), to_f64(&be), to_f64(&dd), to_f64(&na));
        let s0f = [to_f64(&s0n[0]), to_f64(&s0n[1])];
        let wf = [to_f64(&w[0]), to_f64(&w[1])];
        for root in [1i64, -1] {
            let sg = (-bef + root as f64 * ddf.sqrt()) / (2.0 * alf);
            let s = [(s0f[0] + sg * wf[0]) / naf, (s0f[1] + sg * wf[1]) / naf];
            let lim = -1e-9 * 1f64.max(s[0].abs()).max(s[1].abs());
            if !p.iter().all(|r| to_f64(&r[0]) * s[0] + to_f64(&r[1]) * s[1] + to_f64(&r[2]) >= lim) {
                continue;
            }
            let gq = [2.0 * tf[0] * s[0] + tf[1] * s[1] + tf[3], tf[1] * s[0] + 2.0 * tf[2] * s[1] + tf[4]];
            let sd = side as f64;
            let mm = [[to_f64(&a[0]), sd * gq[0]], [to_f64(&a[1]), sd * gq[1]]];
            let det_m = mm[0][0] * mm[1][1] - mm[0][1] * mm[1][0];
            let mx = mm.iter().flatten().fold(1f64, |acc, x| acc.max(x.abs()));
            if det_m.abs() < 1e-12 * mx {
                continue;
            }
            let pp = [(nf[0][0] * s[0] + nf[0][1] * s[1] + mf[0]) / df, (nf[1][0] * s[0] + nf[1][1] * s[1] + mf[1]) / df];
            let (fb, ok) = branch_of(&fp, s, pp, df, &nf, mf, qqf)?;
            if !ok {
                continue;
            }
            let mut face = BFace::blank(FaceKind::Corner);
            face.fp = fp.clone();
            face.fb = fb;
            face.s = s;
            face.sigma = Some(Sigma { quad: [al.clone(), be.clone(), ga.clone()], root, facet: i + 1, s0n: s0n.clone(), w: w.clone(), na: na.clone() });
            face.mm = mm;
            face.p = pp;
            face.src = Some(Src { n: n.to_vec(), d: d.clone(), tie: tie.to_vec(), side });
            face.con = corner_wedge_rows(&face, n, d, tie, side)?;
            faces.push(face);
        }
    }
    Ok(faces)
}

#[allow(clippy::too_many_arguments)]
fn rational_corner(n: &[Z], d: &Z, nn: &[Vec<Z>], m: &[Z], qq: &Z, p: &[Vec<Z>], a: &[Z], s0n: &[Z], w: &[Z], na: &Z, be: &Z, ga: &Z, tie: &[Z], side: i64, i: usize) -> R<Vec<BFace>> {
    let mut sn = vec![&s0n[0] * be - ga * &w[0], &s0n[1] * be - ga * &w[1]];
    let mut sd = na * be;
    if sd.is_negative() {
        sn = neg(&sn);
        sd = -sd;
    }
    let gs = sn[0].gcd(&sn[1]).gcd(&sd);
    let sn = vec![sn[0].div_floor(&gs), sn[1].div_floor(&gs)];
    let sd = sd.div_floor(&gs);
    if p.iter().any(|r| (&r[0] * &sn[0] + &r[1] * &sn[1] + &r[2] * &sd).is_negative()) {
        return Ok(vec![]);
    }
    let mut fp = PolyN::new();
    for (e, c) in [
        ([0u32, 0, 1], zi(2) * d * &sd * &sd),
        ([1, 0, 0], -zi(2) * d * &sd * &sn[0]),
        ([0, 1, 0], -zi(2) * d * &sd * &sn[1]),
        ([0, 0, 0], quad(&sn, nn, &sn) + zi(2) * &sd * dot(m, &sn) + zi(2) * &sd * &sd * qq),
    ] {
        if !c.is_zero() {
            fp.insert(e.to_vec(), c);
        }
    }
    let fp = pn::content_free(&fp);
    let t = vec![vec![zi(2) * &tie[0], tie[1].clone(), tie[3].clone()], vec![tie[1].clone(), zi(2) * &tie[2], tie[4].clone()]];
    let gqn = mv(&t, &[sn[0].clone(), sn[1].clone(), sd.clone()]);
    let s = zi(side);
    let mn = [[a[0].clone(), &s * &gqn[0]], [a[1].clone(), &s * &gqn[1]]];
    let dm = &mn[0][0] * &mn[1][1] - &mn[0][1] * &mn[1][0];
    if dm.is_zero() {
        return Ok(vec![]);
    }
    let pnum = [&nn[0][0] * &sn[0] + &nn[0][1] * &sn[1] + &m[0] * &sd, &nn[1][0] * &sn[0] + &nn[1][1] * &sn[1] + &m[1] * &sd];
    let pd = d * &sd;
    let adj = [[mn[1][1].clone(), -&mn[0][1]], [-&mn[1][0], mn[0][0].clone()]];
    let sdm = zi(sign(&dm) as i64);
    let mut con = vec![];
    for r in 0..2 {
        let rr = vec![-&sdm * &pd * &adj[r][0], -&sdm * &pd * &adj[r][1], &sdm * (&adj[r][0] * &pnum[0] + &adj[r][1] * &pnum[1])];
        let first = rr.iter().find(|x| !x.is_zero()).unwrap();
        let mut row = conic_lin(&rr)?;
        row.push(zi(sign(first) as i64));
        row.push(zi(0));
        con.push(row);
    }
    let sdf = to_f64(&sd);
    let pdf = to_f64(&pd);
    let mut face = BFace::blank(FaceKind::Corner);
    face.fp = Rc::new(fp);
    face.fb = -1;
    face.s = [q_to_f64(&Q::new(sn[0].clone(), sd.clone())), q_to_f64(&Q::new(sn[1].clone(), sd.clone()))];
    face.sigma = Some(Sigma { quad: [zi(0), be.clone(), ga.clone()], root: 0, facet: i, s0n: s0n.to_vec(), w: w.to_vec(), na: na.clone() });
    face.mm = [[to_f64(&a[0]), q_to_f64(&Q::new(&s * &gqn[0], sd.clone()))], [to_f64(&a[1]), q_to_f64(&Q::new(&s * &gqn[1], sd.clone()))]];
    face.p = [q_to_f64(&Q::new(pnum[0].clone(), pd.clone())), q_to_f64(&Q::new(pnum[1].clone(), pd.clone()))];
    let _ = (sdf, pdf);
    face.src = Some(Src { n: n.to_vec(), d: d.clone(), tie: tie.to_vec(), side });
    face.con = con;
    let _ = (q_to_f64(&Q::new(pnum[0].clone(), pd.clone())),);
    Ok(vec![face])
}

fn restrict_tie(tie: &[Z], s0n: &[Z], w: &[Z], na: &Z) -> (Z, Z, Z) {
    let (a, b, c, dc, e, f) = (&tie[0], &tie[1], &tie[2], &tie[3], &tie[4], &tie[5]);
    let al = a * &w[0] * &w[0] + b * &w[0] * &w[1] + c * &w[1] * &w[1];
    let be = zi(2) * a * &s0n[0] * &w[0] + b * (&s0n[0] * &w[1] + &s0n[1] * &w[0]) + zi(2) * c * &s0n[1] * &w[1] + na * (dc * &w[0] + e * &w[1]);
    let ga = a * &s0n[0] * &s0n[0] + b * &s0n[0] * &s0n[1] + c * &s0n[1] * &s0n[1] + na * (dc * &s0n[0] + e * &s0n[1]) + f * na * na;
    (al, be, ga)
}

#[allow(clippy::too_many_arguments)]
fn value_form(nn: &[Vec<Z>], m: &[Z], qq: &Z, d: &Z, s0n: &[Z], w: &[Z], na: &Z, al: &Z, be: &Z, ga: &Z) -> (Vec<Z>, Vec<Z>, Z) {
    let t0 = al * quad(s0n, nn, s0n) - ga * quad(w, nn, w) + zi(2) * na * al * dot(m, s0n) + zi(2) * al * na * na * qq;
    let t1 = zi(2) * al * quad(s0n, nn, w) - be * quad(w, nn, w) + zi(2) * na * al * dot(m, w);
    let e0 = vec![zi(2) * al * na * d * &s0n[0], zi(2) * al * na * d * &s0n[1], -t0];
    let e1 = vec![zi(2) * al * na * d * &w[0], zi(2) * al * na * d * &w[1], -t1];
    let k = zi(2) * al * na * na * d;
    (e0, e1, k)
}

fn aff_plus_z(e: &[Z], kz: &Z) -> PolyN {
    let mut out = PolyN::new();
    for (ex, c) in [([1u32, 0, 0], &e[0]), ([0, 1, 0], &e[1]), ([0, 0, 0], &e[2]), ([0, 0, 1], kz)] {
        if !c.is_zero() {
            out.insert(ex.to_vec(), c.clone());
        }
    }
    out
}

fn implicit_value(e0: &[Z], e1: &[Z], k: &Z, al: &Z, be: &Z, dd: &Z) -> PolyN {
    let lin: Vec<Z> = e0.iter().zip(e1.iter()).map(|(x, y)| -zi(2) * al * x + be * y).collect();
    let l = aff_plus_z(&lin, &(zi(2) * al * k));
    let e = aff_plus_z(e1, &zi(0));
    let fp = pn::subp(&pn::mul(&l, &l), &pn::scale(&pn::mul(&e, &e), dd));
    pn::content_free(&fp)
}

fn branch_of(fp: &PolyN, s: [f64; 2], p: [f64; 2], d: f64, n: &[Vec<f64>], m: [f64; 2], qq: f64) -> R<(i64, bool)> {
    let z = s[0] * p[0] + s[1] * p[1] - ((s[0] * (n[0][0] * s[0] + n[0][1] * s[1]) + s[1] * (n[1][0] * s[0] + n[1][1] * s[1])) / 2.0 + m[0] * s[0] + m[1] * s[1] + qq) / d;
    let (_, all) = alg_branch_info(fp, &[p], -1)?;
    let rr = &all[0];
    if rr.is_empty() {
        return Ok((0, false));
    }
    let gaps: Vec<f64> = rr.iter().map(|r| (r - z).abs()).collect();
    let gap = gaps.iter().cloned().fold(f64::INFINITY, f64::min);
    let j = gaps.iter().position(|g| *g == gap).unwrap();
    if gap > 1e-6 * 1f64.max(z.abs()) {
        return Ok((0, false));
    }
    Ok(((j + 1) as i64 - rr.len() as i64 - 1, true))
}

// ==========================================================================================
// corner_wedge_rows
// ==========================================================================================
fn q1(a: Z) -> Qs {
    [a, zi(0), zi(1)]
}

pub fn corner_wedge_rows(f: &BFace, n: &[Z], d: &Z, tie: &[Z], side: i64) -> R<Vec<Vec<Z>>> {
    let nn = vec![vec![n[4].clone(), n[5].clone()], vec![n[5].clone(), n[6].clone()]];
    let m = vec![n[7].clone(), n[8].clone()];
    let sg = f.sigma.as_ref().unwrap();
    let qq: Vec<Z> = sg.quad.to_vec();
    let na = &sg.na;
    let qmul = |u: &Qs, v: &Qs| qs_mul(u, v, &qq);
    let qneg = |u: &Qs| -> Qs { [-&u[0], -&u[1], u[2].clone()] };
    let s1 = qs_red(&[sg.s0n[0].clone(), sg.w[0].clone(), na.clone()])?;
    let s2 = qs_red(&[sg.s0n[1].clone(), sg.w[1].clone(), na.clone()])?;
    let g1 = qs_add(&qs_add(&qmul(&q1(zi(2) * &tie[0]), &s1)?, &qmul(&q1(tie[1].clone()), &s2)?)?, &q1(tie[3].clone()))?;
    let g2 = qs_add(&qs_add(&qmul(&q1(tie[1].clone()), &s1)?, &qmul(&q1(zi(2) * &tie[2]), &s2)?)?, &q1(tie[4].clone()))?;
    let a = [Z::from(f.mm[0][0] as i64), Z::from(f.mm[1][0] as i64)];
    let mq: [[Qs; 2]; 2] = [[q1(a[0].clone()), qmul(&q1(zi(side)), &g1)?], [q1(a[1].clone()), qmul(&q1(zi(side)), &g2)?]];
    let p1 = qs_red(&qs_add(&qs_add(&qmul(&q1(nn[0][0].clone()), &s1)?, &qmul(&q1(nn[0][1].clone()), &s2)?)?, &q1(m[0].clone()))?)?;
    let p2 = qs_red(&qs_add(&qs_add(&qmul(&q1(nn[1][0].clone()), &s1)?, &qmul(&q1(nn[1][1].clone()), &s2)?)?, &q1(m[1].clone()))?)?;
    let p1 = qs_red(&[p1[0].clone(), p1[1].clone(), &p1[2] * d])?;
    let p2 = qs_red(&[p2[0].clone(), p2[1].clone(), &p2[2] * d])?;
    let adj = [[mq[1][1].clone(), qneg(&mq[0][1])], [qneg(&mq[1][0]), mq[0][0].clone()]];
    let iv = isolate_roots(&qq)?;
    let ivs = iv[if sg.root * (sign(&qq[0]) as i64) > 0 { iv.len() - 1 } else { 0 }].clone();
    let mf = f.mm;
    let xw = [f.p[0] - (mf[0][0] + mf[0][1]), f.p[1] - (mf[1][0] + mf[1][1])];
    let mut rows = vec![];
    for i in 0..2 {
        let (aa, bb) = edge_form(&adj[i], &[p1.clone(), p2.clone()], &qq)?;
        match row_for(&aa, &bb, &qq, xw, &ivs)? {
            None => return Ok(vec![]),
            Some(r) => rows.push(r),
        }
    }
    Ok(rows)
}

fn lcm2(a: &Z, b: &Z) -> Z {
    (a * b).abs() / a.gcd(b)
}

fn edge_form(adj_row: &[Qs; 2], p: &[Qs; 2], q: &[Z]) -> R<(Vec<Z>, Vec<Z>)> {
    let (c1, c2) = (&adj_row[0], &adj_row[1]);
    let k = qs_red(&qs_add(&qs_mul(c1, &p[0], q)?, &qs_mul(c2, &p[1], q)?)?)?;
    let d = lcm2(&lcm2(&c1[2], &c2[2]), &k[2]);
    let mut a = vec![-&c1[0] * (&d / &c1[2]), -&c2[0] * (&d / &c2[2]), &k[0] * (&d / &k[2])];
    let mut b = vec![-&c1[1] * (&d / &c1[2]), -&c2[1] * (&d / &c2[2]), &k[1] * (&d / &k[2])];
    let mut all = a.clone();
    all.extend(b.iter().cloned());
    let g = gcd_all(&all).abs();
    if g.is_positive() {
        a = a.iter().map(|x| x.div_floor(&g)).collect();
        b = b.iter().map(|x| x.div_floor(&g)).collect();
    }
    Ok((a, b))
}

fn lin_mul(a: &[Z], b: &[Z]) -> Vec<Z> {
    vec![
        &a[0] * &b[0],
        &a[0] * &b[1] + &a[1] * &b[0],
        &a[1] * &b[1],
        &a[0] * &b[2] + &a[2] * &b[0],
        &a[1] * &b[2] + &a[2] * &b[1],
        &a[2] * &b[2],
    ]
}

fn maxabs_i(l: &[Z]) -> f64 {
    let m = l.iter().map(|x| x.abs()).max().unwrap();
    to_f64(&m).max(1.0)
}

fn row_for(a: &[Z], b: &[Z], q: &[Z], xw: [f64; 2], ivs: &(Q, Q)) -> R<Option<Vec<Z>>> {
    if all_zero(b) {
        let g = gcd_all(a).abs();
        if g.is_zero() {
            return Ok(None);
        }
        let l: Vec<Z> = a.iter().map(|x| x.div_floor(&g)).collect();
        let v = to_f64(&l[0]) * xw[0] + to_f64(&l[1]) * xw[1] + to_f64(&l[2]);
        if v.abs() <= 1e-9 * maxabs_i(&l) {
            return Ok(None);
        }
        let mut r = vec![zi(0), zi(0), zi(0)];
        r.extend(l);
        r.push(zi(if v > 0.0 { 1 } else { -1 }));
        r.push(zi(0));
        return Ok(Some(r));
    }
    if (&a[0] * &b[1] - &a[1] * &b[0]).is_zero() && (&a[0] * &b[2] - &a[2] * &b[0]).is_zero() && (&a[1] * &b[2] - &a[2] * &b[1]).is_zero() {
        let g = gcd_all(b).abs();
        let l: Vec<Z> = b.iter().map(|x| x.div_floor(&g)).collect();
        let v = to_f64(&l[0]) * xw[0] + to_f64(&l[1]) * xw[1] + to_f64(&l[2]);
        if v.abs() <= 1e-9 * maxabs_i(&l) {
            return Ok(None);
        }
        let mut r = vec![zi(0), zi(0), zi(0)];
        r.extend(l);
        r.push(zi(if v > 0.0 { 1 } else { -1 }));
        r.push(zi(0));
        return Ok(Some(r));
    }
    let (al, be, ga) = (&q[0], &q[1], &q[2]);
    let aa = lin_mul(a, a);
    let ab = lin_mul(a, b);
    let bb = lin_mul(b, b);
    let c: Vec<Z> = (0..6).map(|i| al * &aa[i] - be * &ab[i] + ga * &bb[i]).collect();
    if all_zero(&c) {
        return Ok(None);
    }
    let c = conic(&c)?;
    let pd = match pair_line_data_m(&c) {
        Err(_) => return Ok(None),
        Ok(p) => p,
    };
    let (br, ok) = branch_of_pair(&pd, a, b, q, ivs)?;
    if !ok {
        return Ok(None);
    }
    let v = pair_line_eval(&c, br, xw[0], xw[1])?;
    if v.abs() <= 1e-9 * maxabs_i(&c) {
        return Ok(None);
    }
    let mut r = c;
    r.push(zi(if v > 0.0 { 1 } else { -1 }));
    r.push(zi(br));
    Ok(Some(r))
}

fn qs_sign_big(u: &[Z; 3], q: &[Z], iv: &(Q, Q)) -> R<i32> {
    if u[1].is_zero() {
        return Ok(sign(&u[0]));
    }
    Ok(isolated_root_sign(&[u[1].clone(), u[0].clone()], q, iv)? * sign(&u[2]))
}

fn branch_of_pair(pd: &PairLine, a: &[Z], b: &[Z], q: &[Z], ivs: &(Q, Q)) -> R<(i64, bool)> {
    let (den, num): ([Z; 3], Box<dyn Fn(&Z, &Z) -> [Z; 3]>);
    if pd.parallel {
        let k = pd.form.iter().position(|x| !x.is_zero()).unwrap();
        let fk = pd.form[k].clone();
        let (a2, b2, ak, bk) = (a[2].clone(), b[2].clone(), a[k].clone(), b[k].clone());
        num = Box::new(move |c0: &Z, c1: &Z| [-&a2 * &fk * c1 - c0 * &ak, -&b2 * &fk * c1 - c0 * &bk, zi(1)]);
        den = [a[k].clone(), b[k].clone(), zi(1)];
    } else {
        let sv = qs_sign_big(&[a[1].clone(), b[1].clone(), zi(1)], q, ivs)?;
        if pd.vertical {
            return Ok((1 + if sv == 0 { 1 } else { 0 }, true));
        }
        if sv == 0 {
            return Ok((0, false));
        }
        let (a0, b0, a1, b1) = (a[0].clone(), b[0].clone(), a[1].clone(), b[1].clone());
        num = Box::new(move |c0: &Z, c1: &Z| [-&a0 * c1 - c0 * &a1, -&b0 * c1 - c0 * &b1, zi(1)]);
        den = [a[1].clone(), b[1].clone(), zi(1)];
    }
    let sd = qs_sign_big(&den, q, ivs)?;
    if sd == 0 {
        return Ok((0, false));
    }
    let mm = short_between(&pd.iv[0].1, &pd.iv[1].0);
    let st = sd * qs_sign_big(&num(mm.numer(), mm.denom()), q, ivs)?;
    Ok((1 + if st > 0 { 1 } else { 0 }, true))
}

/// A dyadic rational with few bits strictly between a and b (or a itself when a == b).
fn short_between(a: &Q, b: &Q) -> Q {
    if a == b {
        return a.clone();
    }
    let mx = if a.abs() > b.abs() { a.abs() } else { b.abs() };
    let k0 = (q_to_f64(&mx) + 1.0).log2().ceil() as i64;
    let mut k = k0;
    while k > -1075 {
        let step = if k >= 0 { Q::from_integer(Z::one() << k as usize) } else { Q::new(Z::one(), Z::one() << (-k) as usize) };
        let mut mm = (a / &step).ceil() * &step;
        if &mm == a {
            mm = mm + &step;
        }
        if &mm < b {
            return mm;
        }
        k -= 1;
    }
    (a + b) / qi(2)
}

pub fn _unused() -> Option<i64> {
    zi(0).to_i64()
}

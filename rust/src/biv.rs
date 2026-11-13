//! biv -- exact BIVARIATE integer polynomials as dense coefficient matrices (port of biv.py).
//! M[i][j] multiplies x^i y^j; a matrix always has at least one row and one column, and its
//! SHAPE is kept as MATLAB keeps it.

use crate::bail;
use crate::errors::R;
use crate::poly_n;
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

pub type M2 = Vec<Vec<Z>>;

pub fn shape(a: &M2) -> (usize, usize) {
    (a.len(), a[0].len())
}

pub fn zeros(n1: usize, n2: usize) -> M2 {
    vec![vec![Z::zero(); n2]; n1]
}

pub fn is_const(a: &M2) -> bool {
    a.len() == 1 && a[0].len() == 1
}

pub fn is_zero(a: &M2) -> bool {
    a.iter().all(|r| r.iter().all(|x| x.is_zero()))
}

pub fn trim_n(a: &M2) -> M2 {
    let mut a = a.clone();
    while a.len() > 1 && a.last().unwrap().iter().all(|x| x.is_zero()) {
        a.pop();
    }
    let mut n2 = a[0].len();
    while n2 > 1 && a.iter().all(|r| r[n2 - 1].is_zero()) {
        n2 -= 1;
    }
    a.iter().map(|r| r[..n2].to_vec()).collect()
}

pub fn content(a: &M2) -> Z {
    let mut g = Z::zero();
    for r in a {
        for x in r {
            if !x.is_zero() {
                g = g.gcd(x);
                if g.is_one() {
                    return g;
                }
            }
        }
    }
    g
}

pub fn content_free(a: &M2) -> M2 {
    let g = content(a);
    if g > Z::one() {
        a.iter().map(|r| r.iter().map(|x| x / &g).collect()).collect()
    } else {
        a.clone()
    }
}

pub fn scal(a: &M2, k: &Z) -> M2 {
    a.iter().map(|r| r.iter().map(|x| x * k).collect()).collect()
}

pub fn scali(a: &M2, k: i64) -> M2 {
    scal(a, &zi(k))
}

pub fn conv2(a: &M2, b: &M2) -> M2 {
    let (a1, a2) = shape(a);
    let (b1, b2) = shape(b);
    let mut c = zeros(a1 + b1 - 1, a2 + b2 - 1);
    for i in 0..a1 {
        for j in 0..a2 {
            let v = &a[i][j];
            if v.is_zero() {
                continue;
            }
            for k in 0..b1 {
                for l in 0..b2 {
                    if !b[k][l].is_zero() {
                        c[i + k][j + l] += v * &b[k][l];
                    }
                }
            }
        }
    }
    c
}

pub fn pad_add(a: &M2, b: &M2) -> M2 {
    let (a1, a2) = shape(a);
    let (b1, b2) = shape(b);
    let mut s = zeros(a1.max(b1), a2.max(b2));
    for i in 0..a1 {
        for j in 0..a2 {
            s[i][j] = a[i][j].clone();
        }
    }
    for i in 0..b1 {
        for j in 0..b2 {
            s[i][j] += &b[i][j];
        }
    }
    s
}

pub fn pad(a: &M2, n1: usize, n2: usize) -> M2 {
    let (a1, a2) = shape(a);
    if (a1, a2) == (n1, n2) {
        return a.clone();
    }
    let mut s = zeros(n1, n2);
    for i in 0..a1 {
        for j in 0..a2 {
            s[i][j] = a[i][j].clone();
        }
    }
    s
}

pub fn transpose(a: &M2) -> M2 {
    let (n1, n2) = shape(a);
    (0..n2).map(|j| (0..n1).map(|i| a[i][j].clone()).collect()).collect()
}

/// Derivative in x (d = 1) or y (d = 2).
pub fn diff(s: &M2, d: usize) -> M2 {
    if d == 2 {
        return transpose(&diff(&transpose(s), 1));
    }
    if s.len() < 2 {
        return vec![vec![Z::zero(); s[0].len()]];
    }
    (1..s.len()).map(|i| s[i].iter().map(|x| x * zi(i as i64)).collect()).collect()
}

pub fn equal(a: &M2, b: &M2) -> bool {
    a == b
}

pub fn to_dict(a: &M2) -> poly_n::PolyN {
    let mut d = poly_n::PolyN::new();
    for (i, r) in a.iter().enumerate() {
        for (j, x) in r.iter().enumerate() {
            if !x.is_zero() {
                d.insert(vec![i as u32, j as u32], x.clone());
            }
        }
    }
    d
}

pub fn from_dict(d: &poly_n::PolyN) -> M2 {
    if d.is_empty() {
        return vec![vec![Z::zero()]];
    }
    let n1 = d.keys().map(|e| e[0]).max().unwrap() as usize + 1;
    let n2 = d.keys().map(|e| e[1]).max().unwrap() as usize + 1;
    let mut m = zeros(n1, n2);
    for (e, c) in d {
        m[e[0] as usize][e[1] as usize] = c.clone();
    }
    m
}

pub fn poly_div_exact(p: &M2, d: &M2) -> R<M2> {
    Ok(from_dict(&poly_n::div_exact(&to_dict(p), &to_dict(d))?))
}

// ---- the gcd -----------------------------------------------------------------------------------
fn trim1(a: &[Z]) -> Vec<Z> {
    let mut a = a.to_vec();
    while a.len() > 1 && a.last().unwrap().is_zero() {
        a.pop();
    }
    a
}

fn rev(a: &[Z]) -> Vec<Z> {
    a.iter().rev().cloned().collect()
}

fn y_gcd(a: &[Z], b: &[Z]) -> R<Vec<Z>> {
    let a = trim1(a);
    let b = trim1(b);
    if a.iter().all(|x| x.is_zero()) {
        return Ok(b);
    }
    if b.iter().all(|x| x.is_zero()) {
        return Ok(a);
    }
    let g = poly_gcd(&rev(&a), &rev(&b))?;
    Ok(trim1(&rev(&g)))
}

fn y_div_exact(a: &[Z], b: &[Z]) -> R<Vec<Z>> {
    let a = trim1(a);
    let b = trim1(b);
    if b.len() == 1 {
        if b[0].is_zero() {
            bail!("ratQ:divideByZero", "y-division by zero.");
        }
        if a.iter().any(|x| !x.mod_floor(&b[0]).is_zero()) {
            bail!("ratQ:notDivisible", "the y-division is not exact.");
        }
        return Ok(a.iter().map(|x| x.div_floor(&b[0])).collect());
    }
    match poly_div_exact1(&rev(&a), &rev(&b)) {
        None => bail!("ratQ:notDivisible", "the y-division is not exact."),
        Some(q) => {
            if q.is_empty() {
                Ok(vec![Z::zero()])
            } else {
                Ok(trim1(&rev(&q)))
            }
        }
    }
}

type Cellp = Vec<Vec<Z>>;

fn to_cell(a: &M2) -> Cellp {
    let a = trim_n(a);
    let mut c: Cellp = a.iter().map(|r| trim1(r)).collect();
    while c.len() > 1 && c.last().unwrap().iter().all(|x| x.is_zero()) {
        c.pop();
    }
    if c.len() == 1 && c[0].iter().all(|x| x.is_zero()) {
        return vec![];
    }
    c
}

fn from_cell(c: &Cellp) -> M2 {
    if c.is_empty() {
        return vec![vec![Z::zero()]];
    }
    let w = c.iter().map(|r| r.len()).max().unwrap();
    trim_n(&c.iter().map(|r| {
        let mut v = r.clone();
        v.resize(w, Z::zero());
        v
    }).collect())
}

fn primitive(c: &Cellp) -> R<(Cellp, Vec<Z>)> {
    let mut cont = vec![Z::zero()];
    for r in c {
        cont = y_gcd(&cont, r)?;
    }
    if cont.len() == 1 && cont[0].abs().is_one() {
        return Ok((c.clone(), vec![Z::one()]));
    }
    let mut out = vec![];
    for r in c {
        out.push(y_div_exact(r, &cont)?);
    }
    Ok((out, cont))
}

fn conv1(a: &[Z], b: &[Z]) -> Vec<Z> {
    let mut out = vec![Z::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

fn is_zero_cell(r: &Cellp) -> bool {
    r.is_empty() || (r.len() == 1 && r[0].iter().all(|x| x.is_zero()))
}

fn pseudo_rem(a: &Cellp, b: &Cellp) -> R<Cellp> {
    let db = b.len() - 1;
    let mut r: Cellp = a.clone();
    let mut guard = 0;
    loop {
        while r.len() > 1 && r.last().unwrap().iter().all(|x| x.is_zero()) {
            r.pop();
        }
        if is_zero_cell(&r) {
            return Ok(vec![]);
        }
        let dr = r.len() - 1;
        if dr < db {
            return Ok(r);
        }
        guard += 1;
        if guard > 10000 {
            bail!("ratQ:premNoConverge", "the bivariate pseudo-remainder did not terminate.");
        }
        let lb = b.last().unwrap().clone();
        let lr = r.last().unwrap().clone();
        r = r.iter().map(|x| conv1(x, &lb)).collect();
        let sh = dr - db;
        for i in 0..=db {
            let t = conv1(&lr, &b[i]);
            let n = r[sh + i].len().max(t.len());
            let mut u = r[sh + i].clone();
            u.resize(n, Z::zero());
            for (k, v) in t.iter().enumerate() {
                u[k] -= v;
            }
            r[sh + i] = trim1(&u);
        }
        let last = r.len() - 1;
        r[last] = vec![Z::zero()];
        while r.len() > 1 && r.last().unwrap().iter().all(|x| x.is_zero()) {
            r.pop();
        }
        if is_zero_cell(&r) {
            return Ok(vec![]);
        }
        r = primitive(&r)?.0;
    }
}

/// The gcd of two bivariate integer polynomials, content-free (sign not normalised).
pub fn poly_gcd_biv(a: &M2, b: &M2) -> R<M2> {
    let a0 = content_free(a);
    let b0 = content_free(b);
    let a = to_cell(&a0);
    let b = to_cell(&b0);
    if a.is_empty() {
        return Ok(b0);
    }
    if b.is_empty() {
        return Ok(a0);
    }
    let (mut a, ca) = primitive(&a)?;
    let (mut b, cb) = primitive(&b)?;
    let cg = y_gcd(&ca, &cb)?;
    let mut guard = 0;
    loop {
        guard += 1;
        if guard > 200 {
            bail!("ratQ:gcdNoConverge", "the bivariate gcd did not terminate.");
        }
        if b.len() == 1 {
            return Ok(content_free(&from_cell(&vec![cg])));
        }
        let r = pseudo_rem(&a, &b)?;
        if r.is_empty() {
            let bb: Cellp = b.iter().map(|x| trim1(&conv1(x, &cg))).collect();
            return Ok(content_free(&from_cell(&bb)));
        }
        a = b;
        b = primitive(&r)?.0;
    }
}

pub fn from_i(v: &[&[i64]]) -> M2 {
    v.iter().map(|r| r.iter().map(|x| zi(*x)).collect()).collect()
}

pub fn konst(c: Z) -> M2 {
    vec![vec![c]]
}

pub fn neg(a: &M2) -> M2 {
    a.iter().map(|r| r.iter().map(|x| -x).collect()).collect()
}

pub fn is_neg_one_or_one(x: &Z) -> bool {
    x.abs().is_one()
}

pub fn sgn(x: &Z) -> i64 {
    if x.is_positive() {
        1
    } else if x.is_negative() {
        -1
    } else {
        0
    }
}

//! rat_q -- EXACT rational arithmetic on coefficient vectors (port of rat_q.py / ratQ.m).
//!
//! Conventions (as in Python): a univariate polynomial is a Vec of integers, HIGHEST degree
//! first; the empty Vec is the zero polynomial. A conic is [a b c d e f] for
//! a x^2 + b xy + c y^2 + d x + e y + f. A face row is 10 wide in the weighted cubic basis.
//! There is ONE path on unbounded integers (no 2^53 wall), as in the Python port.

use crate::memo::Memo;
use std::sync::LazyLock;
use crate::bail;
use crate::errors::R;
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::sync::Arc;

pub type Z = BigInt;
pub type Q = BigRational;

#[inline]
pub fn zi(x: i64) -> Z {
    BigInt::from(x)
}
#[inline]
pub fn qz(x: &Z) -> Q {
    Q::from_integer(x.clone())
}
#[inline]
pub fn qi(x: i64) -> Q {
    Q::from_integer(zi(x))
}
pub fn zv(v: &[i64]) -> Vec<Z> {
    v.iter().map(|&x| zi(x)).collect()
}

// ------------------------------------------------------------------------------------------
// small helpers
// ------------------------------------------------------------------------------------------
#[inline]
pub fn sign(x: &Z) -> i32 {
    match x.sign() {
        num_bigint::Sign::Plus => 1,
        num_bigint::Sign::Minus => -1,
        num_bigint::Sign::NoSign => 0,
    }
}
#[inline]
pub fn signq(x: &Q) -> i32 {
    if x.is_positive() {
        1
    } else if x.is_negative() {
        -1
    } else {
        0
    }
}
#[inline]
pub fn signf(x: f64) -> i32 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

pub fn gcd(a: &Z, b: &Z) -> Z {
    a.gcd(b)
}

pub fn lcm(a: &Z, b: &Z) -> Z {
    let a = a.abs();
    let b = b.abs();
    if a.is_zero() || b.is_zero() {
        return Z::zero();
    }
    &a / a.gcd(&b) * b
}

pub fn dot(u: &[Z], v: &[Z]) -> Z {
    let mut s = Z::zero();
    for (a, b) in u.iter().zip(v.iter()) {
        s += a * b;
    }
    s
}

pub fn mv(m: &[Vec<Z>], v: &[Z]) -> Vec<Z> {
    m.iter().map(|r| dot(r, v)).collect()
}

pub fn quad(u: &[Z], m: &[Vec<Z>], v: &[Z]) -> Z {
    dot(u, &mv(m, v))
}

pub fn neg(v: &[Z]) -> Vec<Z> {
    v.iter().map(|x| -x).collect()
}

pub fn all_zero(v: &[Z]) -> bool {
    v.iter().all(|x| x.is_zero())
}

pub fn gcd_all(v: &[Z]) -> Z {
    let mut g = Z::zero();
    for x in v {
        if !x.is_zero() {
            g = g.gcd(x);
            if g.is_one() {
                break;
            }
        }
    }
    g
}

/// Python `math.gcd` over several values: nonnegative, 0 only when all are 0.
pub fn gcd_list(v: &[&Z]) -> Z {
    let mut g = Z::zero();
    for x in v {
        g = g.gcd(x);
    }
    g
}

pub fn to_f64(x: &Z) -> f64 {
    x.to_f64().unwrap_or(if x.is_negative() { f64::NEG_INFINITY } else { f64::INFINITY })
}
/// a / b as the correctly rounded double (Python's int / int).
pub fn idiv(a: &Z, b: &Z) -> f64 {
    q_to_f64(&Q::new(a.clone(), b.clone()))
}

pub fn q_to_f64(x: &Q) -> f64 {
    // Python float(Fraction): correctly rounded n/d.
    x.to_f64().unwrap_or_else(|| to_f64(x.numer()) / to_f64(x.denom()))
}

// ------------------------------------------------------------------------------------------
// canonical forms
// ------------------------------------------------------------------------------------------
/// Canonical form of the rational vector n/d: gcd removed, denominator positive.
pub fn canon(n: &[Z], d: &Z) -> R<(Vec<Z>, Z)> {
    if d.is_zero() {
        bail!("ratQ:badDenominator", "the denominator must be a nonzero scalar.");
    }
    let mut flat: Vec<Z> = n.to_vec();
    let mut d = d.clone();
    if d.is_negative() {
        for x in flat.iter_mut() {
            *x = -&*x;
        }
        d = -d;
    }
    let mut g = d.clone();
    for x in &flat {
        if !x.is_zero() {
            g = g.gcd(x);
        }
        if g.is_one() {
            break;
        }
    }
    if g > Z::one() {
        for x in flat.iter_mut() {
            *x = &*x / &g;
        }
        d = &d / &g;
    }
    Ok((flat, d))
}

/// canon of a matrix (list of rows), returned in the same shape.
pub fn canon_rows(n: &[Vec<Z>], d: &Z) -> R<(Vec<Vec<Z>>, Z)> {
    if n.is_empty() {
        let (_, d2) = canon(&[], d)?;
        return Ok((vec![], d2));
    }
    let w = n[0].len();
    let flat: Vec<Z> = n.iter().flat_map(|r| r.iter().cloned()).collect();
    let (f, d2) = canon(&flat, d)?;
    Ok((f.chunks(w).map(|c| c.to_vec()).collect(), d2))
}

/// Canonical form of a conic: primitive, first nonzero entry positive.
pub fn conic(c: &[Z]) -> R<Vec<Z>> {
    let nz = c.iter().position(|x| !x.is_zero());
    let nz = match nz {
        Some(i) => i,
        None => bail!(
            "ratQ:zeroConic",
            "an all-zero conic names no curve. In a QuaCon every edge carries its own equation, and a straight edge is [0 0 0 d e f], not all-zero."
        ),
    };
    let g = gcd_all(c);
    let mut out: Vec<Z> = c.iter().map(|x| x / &g).collect();
    if out[nz].is_negative() {
        for x in out.iter_mut() {
            *x = -&*x;
        }
    }
    Ok(out)
}

/// conic([0,0,0] ++ row) for a 3-wide affine row.
pub fn conic_lin(row: &[Z]) -> R<Vec<Z>> {
    let mut c = vec![Z::zero(), Z::zero(), Z::zero()];
    c.extend(row.iter().cloned());
    conic(&c)
}

pub fn same_conic(c1: &[Z], c2: &[Z]) -> R<bool> {
    Ok(conic(c1)? == conic(c2)?)
}

pub fn add(n1: &[Z], d1: &Z, n2: &[Z], d2: &Z) -> R<(Vec<Z>, Z)> {
    let l = lcm(d1, d2);
    let f1 = &l / d1;
    let f2 = &l / d2;
    let n: Vec<Z> = n1.iter().zip(n2.iter()).map(|(a, b)| a * &f1 + b * &f2).collect();
    canon(&n, &l)
}

pub fn sub(n1: &[Z], d1: &Z, n2: &[Z], d2: &Z) -> R<(Vec<Z>, Z)> {
    add(n1, d1, &neg(n2), d2)
}

pub fn scale(n1: &[Z], d1: &Z, kn: &Z, kd: &Z) -> R<(Vec<Z>, Z)> {
    let n: Vec<Z> = n1.iter().map(|x| x * kn).collect();
    canon(&n, &(d1 * kd))
}

pub fn eq_rat(n1: &[Z], d1: &Z, n2: &[Z], d2: &Z) -> R<bool> {
    let (a, p) = canon(n1, d1)?;
    let (b, q) = canon(n2, d2)?;
    Ok(a == b && p == q)
}

/// The conic on which two QUADRATIC faces agree, canonically.
pub fn diff_conic(n1: &[Z], d1: &Z, n2: &[Z], d2: &Z) -> R<Vec<Z>> {
    if !all_zero(&n1[0..4]) || !all_zero(&n2[0..4]) {
        bail!("ratQ:cubicFace", "diffConic is for QUADRATIC faces; a cubic numerator has no conic level set.");
    }
    let (dn, _) = sub(n1, d1, n2, d2)?;
    conic(&face_conic(&dn))
}

/// [dn5, 2 dn6, dn7, 2 dn8, 2 dn9, 2 dn10]: the conic 2*(face) in plain coefficients.
pub fn face_conic(dn: &[Z]) -> Vec<Z> {
    vec![
        dn[4].clone(),
        zi(2) * &dn[5],
        dn[6].clone(),
        zi(2) * &dn[7],
        zi(2) * &dn[8],
        zi(2) * &dn[9],
    ]
}

pub fn discriminant(c: &[Z]) -> Z {
    &c[1] * &c[1] - zi(4) * &c[0] * &c[2]
}

pub fn conic_kind(c: &[Z]) -> &'static str {
    if all_zero(&c[0..3]) {
        return "line";
    }
    let d = discriminant(c);
    match sign(&d) {
        0 => "parabola",
        1 => "hyperbolic",
        _ => "elliptic",
    }
}

pub fn div(n1: &[Z], d1: &Z, n2: &Z, d2: &Z) -> R<(Vec<Z>, Z)> {
    if n2.is_zero() {
        bail!("ratQ:divideByZero", "division by the exact zero: a defect in the caller, not a case to regularise.");
    }
    scale(n1, d1, d2, n2)
}

/// Rows Ns[i]/ds[i] rewritten over ONE common denominator. Returns (N, d).
pub fn combine_den(ns: &[Vec<Z>], ds: &[Z]) -> R<(Vec<Vec<Z>>, Z)> {
    if ns.len() != ds.len() {
        bail!("ratQ:sizeMismatch", "combineDen needs one denominator per row: {} rows, {} denominators.", ns.len(), ds.len());
    }
    if ds.iter().any(|x| x.is_zero()) {
        bail!("ratQ:badDenominator", "every denominator must be nonzero.");
    }
    let mut d = Z::one();
    for x in ds {
        d = lcm(&d, &x.abs());
    }
    let n = ns
        .iter()
        .zip(ds.iter())
        .map(|(row, di)| {
            let f = &d / di.abs();
            let s = if di.is_positive() { zi(1) } else { zi(-1) };
            row.iter().map(|x| x * &f * &s).collect()
        })
        .collect();
    Ok((n, d))
}

pub fn det_exact(a: &[Vec<Z>]) -> R<Z> {
    match a.len() {
        1 => Ok(a[0][0].clone()),
        2 => Ok(&a[0][0] * &a[1][1] - &a[0][1] * &a[1][0]),
        3 => Ok(&a[0][0] * (&a[1][1] * &a[2][2] - &a[1][2] * &a[2][1])
            - &a[0][1] * (&a[1][0] * &a[2][2] - &a[1][2] * &a[2][0])
            + &a[0][2] * (&a[1][0] * &a[2][1] - &a[1][1] * &a[2][0])),
        n => bail!("ratQ:badSystem", "detExact handles n <= 3; got {}.", n),
    }
}

pub fn det2(a: &Z, b: &Z, c: &Z, d: &Z) -> Z {
    a * d - b * c
}

/// Solve A x = b exactly (2x2 or 3x3) by Cramer. Returns (x, xd).
pub fn solve2(a: &[Vec<Z>], b: &[Z]) -> R<(Vec<Z>, Z)> {
    let n = a.len();
    if !(n == 2 || n == 3) || a.iter().any(|r| r.len() != n) || b.len() != n {
        bail!("ratQ:badSystem", "solve2 takes a 2x2 or 3x3 matrix with a matching n x 1 column.");
    }
    let d = det_exact(a)?;
    if d.is_zero() {
        bail!("ratQ:singular", "the system is exactly singular (determinant 0), so it has no unique solution.");
    }
    let mut x = Vec::with_capacity(n);
    for j in 0..n {
        let mut aj: Vec<Vec<Z>> = a.to_vec();
        for i in 0..n {
            aj[i][j] = b[i].clone();
        }
        x.push(det_exact(&aj)?);
    }
    canon(&x, &d)
}

pub fn hess_q(n: &[Z], dn: &Z) -> R<(Vec<Vec<Z>>, Z)> {
    if !all_zero(&n[0..4]) {
        bail!("ratQ:cubicFace", "hessQ is for QUADRATIC faces.");
    }
    canon_rows(&[vec![n[4].clone(), n[5].clone()], vec![n[5].clone(), n[6].clone()]], dn)
}

pub fn grad_q(n: &[Z], dn: &Z) -> R<(Vec<Vec<Z>>, Z)> {
    if !all_zero(&n[0..4]) {
        bail!("ratQ:cubicFace", "gradQ is for QUADRATIC faces.");
    }
    canon_rows(
        &[
            vec![n[4].clone(), n[5].clone(), n[7].clone()],
            vec![n[5].clone(), n[6].clone(), n[8].clone()],
        ],
        dn,
    )
}

// ------------------------------------------------------------------------------------------
// linear feasibility in the plane (Fourier-Motzkin), exact
// ------------------------------------------------------------------------------------------
/// Is { s : a*s1 + b*s2 + c >= 0 for every row [a b c] } nonempty (strict: has interior)?
pub fn feasible2<T: AsRef<[Z]>>(p: &[T], strict: bool) -> bool {
    if p.is_empty() {
        return true;
    }
    let mut qv: Vec<(Z, Z)> = Vec::new();
    let mut lo: Vec<&[Z]> = Vec::new();
    let mut up: Vec<&[Z]> = Vec::new();
    for r in p {
        let r = r.as_ref();
        if r[1].is_positive() {
            lo.push(r);
        } else if r[1].is_negative() {
            up.push(r);
        }
    }
    for r in p {
        let r = r.as_ref();
        if r[1].is_zero() {
            qv.push((r[0].clone(), r[2].clone()));
        }
    }
    for ri in &lo {
        for rj in &up {
            qv.push((
                -&rj[1] * &ri[0] + &ri[1] * &rj[0],
                -&rj[1] * &ri[2] + &ri[1] * &rj[2],
            ));
        }
    }
    let mut have_lo = false;
    let mut have_hi = false;
    let (mut lon, mut lod, mut hin, mut hid) = (Z::zero(), Z::one(), Z::zero(), Z::one());
    for (a, b) in qv {
        if a.is_zero() {
            if (strict && !b.is_positive()) || (!strict && b.is_negative()) {
                return false;
            }
        } else if a.is_positive() {
            if !have_lo || (-&b * &lod - &lon * &a).is_positive() {
                lon = -b;
                lod = a;
                have_lo = true;
            }
        } else {
            let nn = b;
            let nd = -a;
            if !have_hi || (&nn * &hid - &hin * &nd).is_negative() {
                hin = nn;
                hid = nd;
                have_hi = true;
            }
        }
    }
    if have_lo && have_hi {
        let c = &lon * &hid - &hin * &lod;
        if (strict && !c.is_negative()) || (!strict && c.is_positive()) {
            return false;
        }
    }
    true
}

// ------------------------------------------------------------------------------------------
// definiteness
// ------------------------------------------------------------------------------------------
fn conic_matrix(c: &[Z]) -> Vec<Vec<Z>> {
    vec![
        vec![zi(2) * &c[0], c[1].clone(), c[3].clone()],
        vec![c[1].clone(), zi(2) * &c[2], c[4].clone()],
        vec![c[3].clone(), c[4].clone(), zi(2) * &c[5]],
    ]
}

pub fn is_psd2(h: &[Vec<Z>]) -> bool {
    !h[0][0].is_negative() && !h[1][1].is_negative() && !(&h[0][0] * &h[1][1] - &h[0][1] * &h[1][0]).is_negative()
}

pub fn is_psd3(m: &[Vec<Z>]) -> bool {
    if (0..3).any(|i| m[i][i].is_negative()) {
        return false;
    }
    for (a, b) in [(0usize, 1usize), (0, 2), (1, 2)] {
        if (&m[a][a] * &m[b][b] - &m[a][b] * &m[b][a]).is_negative() {
            return false;
        }
    }
    !det_exact(m).unwrap().is_negative()
}

/// +1 if the conic's quadratic form is >= 0 on the whole plane, -1 if <= 0, else 0.
pub fn conic_sign(c: &[Z]) -> i32 {
    let m = conic_matrix(c);
    if is_psd3(&m) {
        return 1;
    }
    let mn: Vec<Vec<Z>> = m.iter().map(|r| neg(r)).collect();
    if is_psd3(&mn) {
        return -1;
    }
    0
}

fn leading_positive(m: &[Vec<Z>]) -> bool {
    m[0][0].is_positive()
        && (&m[0][0] * &m[1][1] - &m[0][1] * &m[1][0]).is_positive()
        && det_exact(m).unwrap().is_positive()
}

pub fn conic_definite(c: &[Z]) -> bool {
    let m = conic_matrix(c);
    let mn: Vec<Vec<Z>> = m.iter().map(|r| neg(r)).collect();
    leading_positive(&m) || leading_positive(&mn)
}

// ------------------------------------------------------------------------------------------
// univariate integer polynomials, highest degree first
// ------------------------------------------------------------------------------------------
pub fn poly_trim(p: &[Z]) -> Vec<Z> {
    match p.iter().position(|x| !x.is_zero()) {
        None => vec![],
        Some(k) => p[k..].to_vec(),
    }
}

pub fn poly_primitive(p: &[Z]) -> Vec<Z> {
    let p = poly_trim(p);
    if p.is_empty() {
        return p;
    }
    let g = gcd_all(&p);
    if g > Z::one() {
        p.iter().map(|x| x / &g).collect()
    } else {
        p
    }
}

/// MATLAB polyder: [0] for a constant (never empty).
pub fn polyder(p: &[Z]) -> Vec<Z> {
    let n = p.len() as i64 - 1;
    if n <= 0 {
        return vec![Z::zero()];
    }
    p[..p.len() - 1].iter().enumerate().map(|(i, c)| c * zi(n - i as i64)).collect()
}

pub fn polyval(p: &[Z], x: &Z) -> Z {
    let mut v = Z::zero();
    for c in p {
        v = v * x + c;
    }
    v
}

pub fn polyval_q(p: &[Z], x: &Q) -> Q {
    let mut v = Q::zero();
    for c in p {
        v = v * x + qz(c);
    }
    v
}

pub fn poly_mul1(a: &[Z], b: &[Z]) -> Vec<Z> {
    if a.is_empty() || b.is_empty() {
        return vec![];
    }
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

pub fn poly_add1(a: &[Z], b: &[Z]) -> Vec<Z> {
    let n = a.len().max(b.len());
    let mut out = vec![Z::zero(); n];
    for (i, x) in a.iter().enumerate() {
        out[n - a.len() + i] += x;
    }
    for (i, x) in b.iter().enumerate() {
        out[n - b.len() + i] += x;
    }
    out
}

/// p / g for integer polynomials when the division is exact; None when it is not.
pub fn poly_div_exact1(p: &[Z], g: &[Z]) -> Option<Vec<Z>> {
    let p = poly_trim(p);
    let g = poly_trim(g);
    if g.is_empty() {
        return None;
    }
    if p.len() < g.len() {
        return if p.is_empty() { Some(vec![]) } else { None };
    }
    let mut r = p.clone();
    let mut q = Vec::with_capacity(p.len() - g.len() + 1);
    for i in 0..(p.len() - g.len() + 1) {
        let (c, rem) = r[i].div_mod_floor(&g[0]);
        if !rem.is_zero() {
            return None;
        }
        if !c.is_zero() {
            for (j, gj) in g.iter().enumerate() {
                r[i + j] -= &c * gj;
            }
        }
        q.push(c);
    }
    if r.iter().any(|x| !x.is_zero()) {
        return None;
    }
    Some(q)
}

pub fn poly_prem_primitive(a: &[Z], b: &[Z]) -> R<Vec<Z>> {
    let a = poly_primitive(a);
    let b = poly_primitive(b);
    if b.is_empty() {
        bail!("ratQ:divideByZero", "pseudo-remainder by the zero polynomial.");
    }
    let db = b.len() - 1;
    let mut r = a;
    let mut guard = 0;
    loop {
        r = poly_trim(&r);
        if r.is_empty() || r.len() - 1 < db {
            return Ok(r);
        }
        guard += 1;
        if guard > 10000 {
            bail!("ratQ:premNoConverge", "the reduced pseudo-remainder did not terminate.");
        }
        let lcb = b[0].clone();
        let lcr = r[0].clone();
        for x in r.iter_mut() {
            *x = &*x * &lcb;
        }
        for i in 0..b.len() {
            r[i] -= &lcr * &b[i];
        }
        r[0] = Z::zero();
        r = poly_primitive(&poly_trim(&r));
    }
}

pub fn poly_pseudo_rem(a: &[Z], b: &[Z]) -> R<Vec<Z>> {
    let a = poly_primitive(a);
    let b = poly_primitive(b);
    if b.is_empty() {
        bail!("ratQ:divideByZero", "pseudo-remainder by the zero polynomial.");
    }
    if a.is_empty() || a.len() < b.len() {
        return Ok(a);
    }
    let da = a.len() - 1;
    let db = b.len() - 1;
    let f = num_traits::pow(b[0].clone(), da - db + 1);
    let mut r: Vec<Z> = a.iter().map(|x| x * &f).collect();
    loop {
        r = poly_trim(&r);
        if r.is_empty() || r.len() - 1 < db {
            break;
        }
        let (c, rem) = r[0].div_mod_floor(&b[0]);
        if !rem.is_zero() {
            bail!("ratQ:pseudoRemNotIntegral", "the pseudo-remainder step did not divide evenly.");
        }
        for k in 0..b.len() {
            r[k] -= &c * &b[k];
        }
        r[0] = Z::zero();
    }
    Ok(poly_primitive(&poly_trim(&r)))
}

/// gcd of two integer polynomials, primitive; [1] when coprime.
pub fn poly_gcd(a: &[Z], b: &[Z]) -> R<Vec<Z>> {
    let mut a = poly_primitive(a);
    let mut b = poly_primitive(b);
    if a.is_empty() {
        return Ok(b);
    }
    if b.is_empty() {
        return Ok(a);
    }
    loop {
        let r = poly_trim(&poly_prem_primitive(&a, &b)?);
        if r.is_empty() {
            return Ok(poly_primitive(&b));
        }
        a = b;
        b = poly_primitive(&r);
        if b.len() <= 1 {
            return Ok(vec![Z::one()]);
        }
    }
}

fn poly_squarefree_build(p: &[Z]) -> Vec<Z> {
    let p = poly_trim(p);
    if p.len() <= 2 {
        return p;
    }
    let g = match poly_gcd(&p, &polyder(&p)) {
        Ok(g) => g,
        Err(_) => return p,
    };
    if g.len() <= 1 {
        return p;
    }
    match poly_div_exact1(&p, &g) {
        None => p,
        Some(q) => poly_trim(&q),
    }
}

static SQF_MEMO: LazyLock<Memo<Vec<Z>, Arc<Vec<Z>>>> = LazyLock::new(|| Memo::new(20000));
static CHAIN_MEMO: LazyLock<Memo<Vec<Z>, Arc<Vec<Vec<Z>>>>> = LazyLock::new(|| Memo::new(20000));

pub fn clear_memos() {
    SQF_MEMO.clear();
    CHAIN_MEMO.clear();
}

/// The squarefree part of p (same roots, each simple); memoised.
pub fn poly_squarefree(p: &[Z]) -> Vec<Z> {
    if let Some(r) = SQF_MEMO.get(p) {
        return (*r).clone();
    }
    let r = poly_squarefree_build(p);
    SQF_MEMO.insert(p.to_vec(), Arc::new(r.clone()));
    r
}

fn sturm_chain_build(p: &[Z]) -> Vec<Vec<Z>> {
    let p = poly_trim(p);
    let p = poly_squarefree(&p);
    let mut chain = vec![p.clone()];
    if p.len() <= 1 {
        return chain;
    }
    chain.push(poly_trim(&polyder(&p)));
    while chain.last().unwrap().len() > 1 {
        let a = poly_primitive(&chain[chain.len() - 2]);
        let b = poly_primitive(&chain[chain.len() - 1]);
        let mut r = match poly_pseudo_rem(&a, &b) {
            Ok(r) => poly_trim(&r),
            Err(_) => break,
        };
        if r.is_empty() {
            break;
        }
        if b[0].is_negative() && (a.len() as i64 - b.len() as i64 + 1).rem_euclid(2) == 1 {
            r = neg(&r);
        }
        chain.push(poly_primitive(&neg(&r)));
    }
    chain
}

/// sturmChain(p), shared from a memo (callers must not mutate it).
pub fn sturm_chain_shared(p: &[Z]) -> Arc<Vec<Vec<Z>>> {
    if let Some(c) = CHAIN_MEMO.get(p) {
        return c;
    }
    let c = Arc::new(sturm_chain_build(p));
    CHAIN_MEMO.insert(p.to_vec(), c.clone());
    c
}

pub fn sturm_chain(p: &[Z]) -> Vec<Vec<Z>> {
    (*sturm_chain_shared(p)).clone()
}

/// An extended point for Sturm counts: a rational or +-infinity.
#[derive(Clone, Debug)]
pub enum XPt {
    Fin(Q),
    PosInf,
    NegInf,
}

fn sign_at_inf(c: &[Z], pos: bool) -> i32 {
    if c.is_empty() {
        return 0;
    }
    let mut s = sign(&c[0]);
    if !pos && (c.len() - 1) % 2 == 1 {
        s = -s;
    }
    s
}

fn variations<I: IntoIterator<Item = i32>>(signs: I) -> usize {
    let mut v = 0;
    let mut last = 0;
    for s in signs {
        if s == 0 {
            continue;
        }
        if last != 0 && s != last {
            v += 1;
        }
        last = s;
    }
    v
}

/// sign of c(n/d) for d > 0, dp[i] = d^i: the homogeneous Horner sum in integers.
fn sign_hom(c: &[Z], n: &Z, dp: &[Z]) -> i32 {
    if c.is_empty() {
        return 0;
    }
    let mut v = c[0].clone();
    for i in 1..c.len() {
        v = v * n + &c[i] * &dp[i];
    }
    sign(&v)
}

pub fn sturm_v_q(chain: &[Vec<Z>], x: &Q) -> usize {
    let n = x.numer();
    let d = x.denom();
    if d.is_one() {
        return variations(chain.iter().map(|c| sign(&polyval(c, n))));
    }
    let k = chain.iter().map(|c| c.len()).max().unwrap_or(0);
    let mut dp = vec![Z::one(); k.max(1)];
    for i in 1..k {
        dp[i] = &dp[i - 1] * d;
    }
    variations(chain.iter().map(|c| sign_hom(c, n, &dp)))
}

pub fn sturm_v(chain: &[Vec<Z>], x: &XPt) -> usize {
    match x {
        XPt::PosInf => variations(chain.iter().map(|c| sign_at_inf(c, true))),
        XPt::NegInf => variations(chain.iter().map(|c| sign_at_inf(c, false))),
        XPt::Fin(q) => sturm_v_q(chain, q),
    }
}

/// Number of distinct real roots of p in the half-open interval (a, b]. Exact.
pub fn count_roots_in(p: &[Z], a: &XPt, b: &XPt) -> usize {
    let chain = sturm_chain_shared(&poly_trim(p));
    if chain[0].len() <= 1 {
        return 0;
    }
    let va = sturm_v(&chain, a);
    let vb = sturm_v(&chain, b);
    va.saturating_sub(vb)
}

pub fn count_roots_in_q(p: &[Z], a: &Q, b: &Q) -> usize {
    count_roots_in(p, &XPt::Fin(a.clone()), &XPt::Fin(b.clone()))
}

fn bisect(chain: &[Vec<Z>], lo: Q, hi: Q, n: usize, depth: usize, out: &mut Vec<(Q, Q)>) -> R<()> {
    if n == 0 {
        return Ok(());
    }
    if n == 1 {
        out.push((lo, hi));
        return Ok(());
    }
    if depth > 200 {
        bail!("ratQ:isolateNoConverge", "bisection did not separate {} roots in 200 levels.", n);
    }
    let mid = (&lo + &hi) / qi(2);
    let vl = sturm_v_q(chain, &lo);
    let vm = sturm_v_q(chain, &mid);
    let vh = sturm_v_q(chain, &hi);
    bisect(chain, lo, mid.clone(), vl - vm, depth + 1, out)?;
    bisect(chain, mid, hi, vm - vh, depth + 1, out)
}

/// Isolating intervals (lo, hi] for the real roots of p, sorted. Endpoints are rationals.
pub fn isolate_roots(p: &[Z]) -> R<Vec<(Q, Q)>> {
    let p = poly_squarefree(&poly_trim(p));
    if p.len() <= 1 {
        return Ok(vec![]);
    }
    let mx = p[1..].iter().map(|x| x.abs()).max().unwrap();
    let m = zi(1) + (-((-mx).div_floor(&p[0].abs())));
    let chain = sturm_chain_shared(&p);
    let lo = Q::from_integer(-&m);
    let hi = Q::from_integer(m);
    let n = sturm_v_q(&chain, &lo) - sturm_v_q(&chain, &hi);
    let mut out = Vec::new();
    bisect(&chain, lo, hi, n, 0, &mut out)?;
    out.sort();
    Ok(out)
}

fn is_dyadic(x: &Q) -> bool {
    let d = x.denom();
    let dm1: Z = d - Z::one();
    (d & dm1).is_zero()
}

fn bitlen(x: &Z) -> u64 {
    x.bits()
}

/// sign of c(N / 2^e): the homogeneous Horner sum, in integers.
fn sign_dy(c: &[Z], n: &Z, e: u64) -> i32 {
    if c.is_empty() {
        return 0;
    }
    let mut v = c[0].clone();
    for i in 1..c.len() {
        v = v * n + (&c[i] << (e as usize * i));
    }
    sign(&v)
}

fn var_dy(chain: &[Vec<Z>], n: &Z, e: u64) -> usize {
    variations(chain.iter().map(|c| sign_dy(c, n, e)))
}

fn sign_at_dyadic(q: &[Z], chain_p: &[Vec<Z>], chain_q: &[Vec<Z>], q_const: bool, lo: &Q, hi: &Q) -> R<i32> {
    let bl_lo = bitlen(lo.denom()) - 1;
    let bl_hi = bitlen(hi.denom()) - 1;
    let mut e = bl_lo.max(bl_hi);
    let mut l: Z = lo.numer() << ((e - bl_lo) as usize);
    let mut h: Z = hi.numer() << ((e - bl_hi) as usize);
    let mut vqlo = if q_const { 0 } else { var_dy(chain_q, &l, e) };
    let mut vqhi = if q_const { 0 } else { var_dy(chain_q, &h, e) };
    let mut vplo = var_dy(chain_p, &l, e);
    for _ in 0..400 {
        if q_const || vqlo == vqhi {
            let s = sign_dy(q, &(&l + &h), e + 1);
            if s != 0 {
                return Ok(s);
            }
            let s = sign_dy(q, &h, e);
            if s != 0 {
                return Ok(s);
            }
            return Ok(sign_dy(q, &l, e));
        }
        let m = &l + &h;
        e += 1;
        l = l << 1usize;
        h = h << 1usize;
        let vpmid = var_dy(chain_p, &m, e);
        let vqmid = var_dy(chain_q, &m, e);
        if vplo as i64 - vpmid as i64 == 1 {
            h = m;
            vqhi = vqmid;
        } else {
            l = m;
            vplo = vpmid;
            vqlo = vqmid;
        }
    }
    bail!("ratQ:signAtNoConverge", "bisection did not separate the root from q's zeros.")
}

/// The sign of q at the root of p isolated by iv = (lo, hi]. Exact.
pub fn sign_at(q: &[Z], p: &[Z], iv: &(Q, Q)) -> R<i32> {
    let p = poly_squarefree(&poly_trim(p));
    let q = poly_trim(q);
    if q.is_empty() {
        return Ok(0);
    }
    let lo = iv.0.clone();
    let hi = iv.1.clone();
    let g = poly_gcd(&p, &q)?;
    if g.len() > 1 && count_roots_in_q(&g, &lo, &hi) > 0 {
        return Ok(0);
    }
    let chain_p = sturm_chain_shared(&p);
    let chain_q = sturm_chain_shared(&q);
    let q_const = chain_q[0].len() <= 1;
    if is_dyadic(&lo) && is_dyadic(&hi) {
        return sign_at_dyadic(&q, &chain_p, &chain_q, q_const, &lo, &hi);
    }
    let mut lo = lo;
    let mut hi = hi;
    let mut vqlo = if q_const { 0 } else { sturm_v_q(&chain_q, &lo) };
    let mut vqhi = if q_const { 0 } else { sturm_v_q(&chain_q, &hi) };
    let mut vplo = sturm_v_q(&chain_p, &lo);
    let two = qi(2);
    for _ in 0..400 {
        if q_const || vqlo == vqhi {
            let v = polyval_q(&q, &((&lo + &hi) / &two));
            if !v.is_zero() {
                return Ok(signq(&v));
            }
            let v = polyval_q(&q, &hi);
            if !v.is_zero() {
                return Ok(signq(&v));
            }
            return Ok(signq(&polyval_q(&q, &lo)));
        }
        let mid = (&lo + &hi) / &two;
        let vpmid = sturm_v_q(&chain_p, &mid);
        let vqmid = sturm_v_q(&chain_q, &mid);
        if vplo as i64 - vpmid as i64 == 1 {
            hi = mid;
            vqhi = vqmid;
        } else {
            lo = mid;
            vplo = vpmid;
            vqlo = vqmid;
        }
    }
    bail!("ratQ:signAtNoConverge", "bisection did not separate the root from q's zeros.")
}

// ------------------------------------------------------------------------------------------
// reading the caller's doubles exactly
// ------------------------------------------------------------------------------------------
pub fn round_half_away(x: f64) -> f64 {
    (x.abs() + 0.5).floor() * if x >= 0.0 { 1.0 } else { -1.0 }
}

/// The denominator of the rational a double was MEANT to be (continued fractions).
pub fn denom_of(x: f64, max_den: f64) -> R<Z> {
    if !x.is_finite() {
        bail!("ratQ:notFinite", "value must be real and finite.");
    }
    if x == x.round() {
        return Ok(Z::one());
    }
    let ax = x.abs();
    let (mut p0, mut q0, mut p1, mut q1) = (1.0f64, 0.0f64, ax.floor(), 1.0f64);
    let mut r = ax - p1;
    for _ in 0..64 {
        if r == 0.0 {
            break;
        }
        let a = (1.0 / r).floor();
        r = 1.0 / r - a;
        let p2 = a * p1 + p0;
        let q2 = a * q1 + q0;
        if q2 > max_den {
            break;
        }
        p0 = p1;
        q0 = q1;
        p1 = p2;
        q1 = q2;
        if p1 / q1 == ax {
            return Ok(Z::from(q1 as i64));
        }
    }
    bail!("ratQ:notExact", "{:.17e} is not a rational with denominator <= {} to within a rounding.", x, max_den)
}

/// A number as the caller wrote it: an exact rational or a double (MATLAB's reading).
#[derive(Clone, Debug)]
pub enum Num {
    Q(Q),
    F(f64),
}

/// The exact rational vector a caller's numbers stand for: (n, d), canonical.
pub fn from_double(xs: &[Num]) -> R<(Vec<Z>, Z)> {
    let max_den = 1e6;
    let limit = Z::one() << 53usize;
    let mut d = Z::one();
    for v in xs {
        let dv = match v {
            Num::Q(q) => q.denom().clone(),
            Num::F(f) => denom_of(*f, max_den)?,
        };
        d = lcm(&d, &dv);
        if d > limit {
            bail!("ratQ:overflow", "the common denominator exceeded 2^53.");
        }
    }
    let mut n = Vec::with_capacity(xs.len());
    let df = to_f64(&d);
    for v in xs {
        match v {
            Num::Q(q) => n.push((q * qz(&d)).to_integer()),
            Num::F(f) => {
                let ni = round_half_away(f * df);
                if ni / df != *f {
                    bail!("ratQ:notExact", "{:.17e} does not round-trip through the common denominator {}.", f, d);
                }
                n.push(Z::from(ni as i128));
            }
        }
    }
    canon(&n, &d)
}

// ------------------------------------------------------------------------------------------
// evaluation (doubles on purpose -- see QuaCon.eval)
// ------------------------------------------------------------------------------------------
const CW: [f64; 10] = [1.0 / 6.0, 0.5, 0.5, 1.0 / 6.0, 0.5, 1.0, 0.5, 1.0, 1.0, 1.0];

/// Value of the weighted cubic row c (6 or 10 wide) at the point (x, y).
pub fn eval_poly_f(c: &[f64], x: f64, y: f64) -> f64 {
    let off = 10 - c.len();
    let mon = [x * x * x, x * x * y, x * y * y, y * y * y, x * x, x * y, y * y, x, y, 1.0];
    let mut s = 0.0;
    for i in 0..10 {
        let ci = if i < off { 0.0 } else { c[i - off] };
        s += mon[i] * CW[i] * ci;
    }
    s
}

pub fn eval_poly(c: &[Z], x: f64, y: f64) -> f64 {
    let cf: Vec<f64> = c.iter().map(to_f64).collect();
    eval_poly_f(&cf, x, y)
}

pub fn eval_face(n: &[Z], d: &Z, x: f64, y: f64) -> f64 {
    eval_poly(n, x, y) / to_f64(d)
}

pub fn eval_conic(c: &[Z], x: f64, y: f64) -> f64 {
    let c: Vec<f64> = c.iter().map(to_f64).collect();
    eval_conic_f(&c, x, y)
}

#[inline]
pub fn eval_conic_f(c: &[f64], x: f64, y: f64) -> f64 {
    c[0] * x.powi(2) + c[1] * x * y + c[2] * y.powi(2) + c[3] * x + c[4] * y + c[5]
}

/// Exact value of the face n/d at the rational point (x1, x2).
pub fn eval_face_exact(n: &[Z], d: &Z, x1: &Q, x2: &Q) -> Q {
    let f = |k: usize, den: i64| Q::new(n[k].clone(), zi(den));
    let v = f(0, 6) * x1 * x1 * x1
        + f(1, 2) * x1 * x1 * x2
        + f(2, 2) * x1 * x2 * x2
        + f(3, 6) * x2 * x2 * x2
        + f(4, 2) * x1 * x1
        + qz(&n[5]) * x1 * x2
        + f(6, 2) * x2 * x2
        + qz(&n[7]) * x1
        + qz(&n[8]) * x2
        + qz(&n[9]);
    v / qz(d)
}

// ------------------------------------------------------------------------------------------
// the sign of a quadratic on a polyhedron
// ------------------------------------------------------------------------------------------
pub fn rec_cone2(p: &[Vec<Z>]) -> (Vec<Vec<Z>>, bool) {
    let n: Vec<[Z; 2]> = p
        .iter()
        .filter(|r| !r[0].is_zero() || !r[1].is_zero())
        .map(|r| [r[0].clone(), r[1].clone()])
        .collect();
    if n.is_empty() {
        return (vec![], false);
    }
    let mut c: Vec<Vec<Z>> = Vec::new();
    for v in &n {
        c.push(vec![v[1].clone(), -&v[0]]);
        c.push(vec![-&v[1], v[0].clone()]);
    }
    let c: Vec<Vec<Z>> = c
        .into_iter()
        .filter(|cc| n.iter().all(|v| !(&v[0] * &cc[0] + &v[1] * &cc[1]).is_negative()))
        .collect();
    let mut u: Vec<Vec<Z>> = Vec::new();
    for cc in c {
        let mut new = true;
        for uu in &u {
            if (&uu[0] * &cc[1] - &uu[1] * &cc[0]).is_zero() && (&uu[0] * &cc[0] + &uu[1] * &cc[1]).is_positive() {
                new = false;
                break;
            }
        }
        if new {
            u.push(cc);
        }
    }
    if u.len() == 2 && (&u[0][0] * &u[1][1] - &u[0][1] * &u[1][0]).is_zero() {
        return (vec![], false);
    }
    if u.len() > 2 {
        return (vec![], false);
    }
    (u, true)
}

fn cone_form_min_sign(h: &[Vec<Z>], u1: &[Z], u2: &[Z]) -> i32 {
    let a = quad(u1, h, u1);
    let b = quad(u1, h, u2);
    let c = quad(u2, h, u2);
    let a2 = &a - zi(2) * &b + &c;
    let a1 = zi(2) * (&b - &c);
    let mut s = sign(&a).min(sign(&c));
    if a2.is_positive() && a1.is_negative() && -&a1 < zi(2) * &a2 {
        s = s.min(sign(&(zi(4) * &a2 * &c - &a1 * &a1)));
    }
    s
}

fn flat_ray_sign(h: &[Vec<Z>], g0: &[Z], p: &[Vec<Z>], u: &[Z]) -> i32 {
    let a = mv(h, u);
    let b = dot(g0, u);
    if all_zero(&a) {
        return if b.is_negative() { -1 } else { 0 };
    }
    let mut rows: Vec<Vec<Z>> = p.to_vec();
    rows.push(vec![-&a[0], -&a[1], -b]);
    if feasible2(&rows, true) {
        -1
    } else {
        0
    }
}

fn quad_sign_at_solve(c: &[Z], p: &[Vec<Z>], a: &[Vec<Z>], b: &[Z]) -> i32 {
    if det_exact(a).unwrap().is_zero() {
        return 1;
    }
    let (xn, xd) = solve2(a, b).unwrap();
    if p.iter().any(|r| (&r[0] * &xn[0] + &r[1] * &xn[1] + &r[2] * &xd).is_negative()) {
        return 1;
    }
    sign(
        &(&c[0] * &xn[0] * &xn[0] + &c[1] * &xn[0] * &xn[1] + &c[2] * &xn[1] * &xn[1]
            + &c[3] * &xn[0] * &xd
            + &c[4] * &xn[1] * &xd
            + &c[5] * &xd * &xd),
    )
}

/// Sign of the minimum of the quadratic c over {P >= 0}: Some(-1/0/1), or None (declined).
pub fn quad_min_sign(c: &[Z], p: &[Vec<Z>]) -> Option<i32> {
    let h = vec![vec![zi(2) * &c[0], c[1].clone()], vec![c[1].clone(), zi(2) * &c[2]]];
    let g0 = vec![c[3].clone(), c[4].clone()];
    let (rec, ok) = rec_cone2(p);
    if !ok {
        return None;
    }
    if rec.len() == 1 {
        let a = quad(&rec[0], &h, &rec[0]);
        if a.is_negative() {
            return Some(-1);
        }
        if a.is_zero() {
            let r = flat_ray_sign(&h, &g0, p, &rec[0]);
            if r < 0 {
                return Some(r);
            }
        }
    } else if rec.len() == 2 {
        let (u1, u2) = (&rec[0], &rec[1]);
        let t = cone_form_min_sign(&h, u1, u2);
        if t < 0 {
            return Some(-1);
        }
        if t == 0 {
            let f1 = quad(u1, &h, u1).is_zero();
            let f2 = quad(u2, &h, u2).is_zero();
            if !f1 && !f2 {
                return None;
            }
            if f1 {
                let r = flat_ray_sign(&h, &g0, p, u1);
                if r < 0 {
                    return Some(r);
                }
            }
            if f2 {
                let r = flat_ray_sign(&h, &g0, p, u2);
                if r < 0 {
                    return Some(r);
                }
            }
        }
    }
    let mut s = 1;
    let k = p.len();
    s = s.min(quad_sign_at_solve(c, p, &h, &[-&g0[0], -&g0[1]]));
    for i in 0..k {
        let n = [p[i][0].clone(), p[i][1].clone()];
        if n[0].is_zero() && n[1].is_zero() {
            continue;
        }
        let u = [n[1].clone(), -&n[0]];
        let uh = vec![&u[0] * &h[0][0] + &u[1] * &h[1][0], &u[0] * &h[0][1] + &u[1] * &h[1][1]];
        s = s.min(quad_sign_at_solve(c, p, &[uh, n.to_vec()], &[-dot(&u, &g0), -&p[i][2]]));
    }
    for i in 0..k.saturating_sub(1) {
        for j in (i + 1)..k {
            s = s.min(quad_sign_at_solve(
                c,
                p,
                &[vec![p[i][0].clone(), p[i][1].clone()], vec![p[j][0].clone(), p[j][1].clone()]],
                &[-&p[i][2], -&p[j][2]],
            ));
        }
    }
    Some(s)
}

/// (s, strict): s = +1 when q >= 0 on {P >= 0}, -1 when q <= 0, 0 undecided.
pub fn quad_sign_on_poly(c: &[Z], p: &[Vec<Z>]) -> R<(i32, bool)> {
    if c.len() != 6 || p.iter().any(|r| r.len() != 3) {
        bail!("ratQ:badQuadSignInput", "quadSignOnPoly takes a 1 x 6 conic and k x 3 half-plane rows.");
    }
    if let Some(mn) = quad_min_sign(c, p) {
        if mn >= 0 {
            return Ok((1, mn > 0));
        }
    }
    if let Some(mx) = quad_min_sign(&neg(c), p) {
        if mx >= 0 {
            return Ok((-1, mx > 0));
        }
    }
    if !all_zero(c) {
        let s = conic_sign(c);
        return Ok((s, s != 0 && conic_definite(c)));
    }
    Ok((0, false))
}

// ------------------------------------------------------------------------------------------
// a degenerate conic as a NAMED pair of lines
// ------------------------------------------------------------------------------------------
#[derive(Clone, Debug)]
pub struct PairLine {
    pub parallel: bool,
    pub apex: [Z; 3],
    pub quad: Vec<Z>,
    pub iv: Vec<(Q, Q)>,
    pub vertical: bool,
    pub form: [Z; 2],
}

fn parallel_pair_data(c: &[Z]) -> R<PairLine> {
    let (a, b, cc, d, e, f) = (&c[0], &c[1], &c[2], &c[3], &c[4], &c[5]);
    let (u, v, q);
    if !a.is_zero() {
        u = zi(2) * a;
        v = b.clone();
        if !(zi(2) * a * e - d * b).is_zero() {
            bail!("ratQ:notALinePair", "the linear part is not a multiple of the repeated factor, so this conic is a PARABOLA rather than two parallel lines.");
        }
        q = vec![zi(1), zi(2) * d, zi(4) * a * f];
    } else if !cc.is_zero() {
        u = zi(0);
        v = zi(1);
        if !d.is_zero() {
            bail!("ratQ:notALinePair", "the linear part is not a multiple of the repeated factor, so this conic is a PARABOLA rather than two parallel lines.");
        }
        q = vec![cc.clone(), e.clone(), f.clone()];
    } else {
        bail!("ratQ:notALinePair", "the quadratic part vanishes, so this conic is a single line, not a pair.");
    }
    let q = poly_trim(&q);
    let iv = isolate_roots(&q)?;
    if iv.len() != 2 {
        bail!("ratQ:notALinePair", "the level quadratic has {} real roots, not two distinct lines.", iv.len());
    }
    Ok(PairLine { parallel: true, apex: [zi(0), zi(0), zi(1)], quad: q, iv, vertical: false, form: [u, v] })
}

static PAIR_MEMO: LazyLock<Memo<Vec<Z>, Arc<PairLine>>> = LazyLock::new(|| Memo::new(20000));

/// pairLineData, memoised on the raw conic (a pure function of it).
pub fn pair_line_data_m(c: &[Z]) -> R<Arc<PairLine>> {
    if let Some(p) = PAIR_MEMO.get(c) {
        return Ok(p);
    }
    let p = Arc::new(pair_line_data(c)?);
    PAIR_MEMO.insert(c.to_vec(), p.clone());
    Ok(p)
}

/// The two lines of a degenerate conic, NAMED.
pub fn pair_line_data(c: &[Z]) -> R<PairLine> {
    let c = conic(c)?;
    let (a, b, cc, d, e, f) = (&c[0], &c[1], &c[2], &c[3], &c[4], &c[5]);
    let disc = b * b - zi(4) * a * cc;
    if disc.is_negative() {
        bail!("ratQ:notALinePair", "this conic's quadratic part is DEFINITE, so its zero set is a point or empty -- not two lines, and there is no branch to name.");
    }
    if disc.is_zero() {
        return parallel_pair_data(&c);
    }
    let m3 = vec![
        vec![zi(2) * a, b.clone(), d.clone()],
        vec![b.clone(), zi(2) * cc, e.clone()],
        vec![d.clone(), e.clone(), zi(2) * f],
    ];
    if !det_exact(&m3)?.is_zero() {
        bail!("ratQ:notDegenerate", "this conic is NOT degenerate, so it is one curve and not two lines. A branch index names nothing on it.");
    }
    let m2 = vec![vec![zi(2) * a, b.clone()], vec![b.clone(), zi(2) * cc]];
    if det_exact(&m2)?.is_zero() {
        bail!("ratQ:parallelPair", "the two lines are PARALLEL.");
    }
    let (mut xn, mut xd) = solve2(&m2, &[-d, -e])?;
    if xd.is_negative() {
        xn = neg(&xn);
        xd = -xd;
    }
    let q = poly_trim(&[cc.clone(), b.clone(), a.clone()]);
    let vertical = cc.is_zero();
    let iv;
    if vertical {
        if b.is_zero() {
            bail!("ratQ:notALinePair", "both lines are vertical, so they are parallel and cannot cross.");
        }
        let t = Q::new(-a, b.clone());
        iv = vec![(&t - qi(1), &t + qi(1)), (qi(0), qi(0))];
    } else {
        let i = isolate_roots(&q)?;
        if i.len() != 2 {
            bail!("ratQ:notALinePair", "the slope quadratic has {} real roots, not two.", i.len());
        }
        iv = i;
    }
    Ok(PairLine {
        parallel: false,
        apex: [xn[0].clone(), xn[1].clone(), xd],
        quad: q,
        iv,
        vertical,
        form: [zi(0), zi(0)],
    })
}

/// Exact sign of branch b of the line pair c at the rational point xn/xd.
pub fn pair_line_sign(c: &[Z], b: i64, xn: &[Z], xd: &Z) -> R<i32> {
    let info = pair_line_data_m(c)?;
    if b != 1 && b != 2 {
        bail!("ratQ:branch", "the branch must be 1 or 2; got {}.", b);
    }
    if info.parallel {
        let tv = &info.form[0] * &xn[0] + &info.form[1] * &xn[1];
        return sign_at(&poly_trim(&[-xd, tv]), &info.quad, &info.iv[(b - 1) as usize]);
    }
    let p = &info.apex;
    let ux = &xn[0] * &p[2] - &p[0] * xd;
    let uy = &xn[1] * &p[2] - &p[1] * xd;
    if info.vertical && b == 2 {
        return Ok(sign(&ux));
    }
    if info.vertical {
        let td = &info.quad[0];
        let tn = -&info.quad[1];
        return Ok(sign(&(&uy * td - &tn * &ux)) * sign(td));
    }
    let q = poly_trim(&[-ux, uy]);
    if q.is_empty() {
        return Ok(0);
    }
    sign_at(&q, &info.quad, &info.iv[(b - 1) as usize])
}

/// sort(roots(quad)) for a quadratic with two real roots, in doubles.
pub fn sorted_real_roots2(q: &[Z]) -> [f64; 2] {
    let a = to_f64(&q[0]);
    let b = to_f64(&q[1]);
    let c = to_f64(&q[2]);
    let disc = b * b - 4.0 * a * c;
    let sq = disc.max(0.0).sqrt();
    let (r1, r2);
    if b != 0.0 {
        let qq = -0.5 * (b + sq.copysign(b));
        r1 = qq / a;
        r2 = if qq != 0.0 { c / qq } else { r1 };
    } else {
        r1 = sq / (2.0 * a);
        r2 = -r1;
    }
    if r1 <= r2 {
        [r1, r2]
    } else {
        [r2, r1]
    }
}

/// The double-precision branch value and term magnitude at one point (pairLineEvalM).
pub struct PairEval {
    pub kind: u8, // 0 parallel, 1 vertical b==2, 2 general (slope t through p0)
    pub f1: f64,
    pub f2: f64,
    pub r: f64,
    pub p0: [f64; 2],
    pub t: f64,
}

pub fn pair_eval_prep(c: &[Z], b: i64) -> R<PairEval> {
    let info = pair_line_data_m(&c[0..6])?;
    if info.parallel {
        let r = sorted_real_roots2(&info.quad)[(b - 1) as usize];
        return Ok(PairEval { kind: 0, f1: to_f64(&info.form[0]), f2: to_f64(&info.form[1]), r, p0: [0.0, 0.0], t: 0.0 });
    }
    let p0 = [
        q_to_f64(&Q::new(info.apex[0].clone(), info.apex[2].clone())),
        q_to_f64(&Q::new(info.apex[1].clone(), info.apex[2].clone())),
    ];
    // Python: info['apex'][0] / info['apex'][2] is int/int true division -> correctly rounded.
    if info.vertical && b == 2 {
        return Ok(PairEval { kind: 1, f1: 0.0, f2: 0.0, r: 0.0, p0, t: 0.0 });
    }
    let t = if info.vertical {
        q_to_f64(&Q::new(-&info.quad[1], info.quad[0].clone()))
    } else {
        sorted_real_roots2(&info.quad)[(b - 1) as usize]
    };
    Ok(PairEval { kind: 2, f1: 0.0, f2: 0.0, r: 0.0, p0, t })
}

impl PairEval {
    #[inline]
    pub fn eval_m(&self, x: f64, y: f64) -> (f64, f64) {
        match self.kind {
            0 => (self.f1 * x + self.f2 * y - self.r, (self.f1 * x).abs() + (self.f2 * y).abs() + self.r.abs()),
            1 => (x - self.p0[0], x.abs() + self.p0[0].abs()),
            _ => (
                (y - self.p0[1]) - self.t * (x - self.p0[0]),
                y.abs() + self.p0[1].abs() + self.t.abs() * (x.abs() + self.p0[0].abs()),
            ),
        }
    }
}

pub fn pair_line_eval(c: &[Z], b: i64, x: f64, y: f64) -> R<f64> {
    Ok(pair_eval_prep(c, b)?.eval_m(x, y).0)
}

// ------------------------------------------------------------------------------------------
// Q(t): elements [n1 n2 D] = (n1 + n2 t)/D with t a root of Q = [al be ga]
// ------------------------------------------------------------------------------------------
pub type Qs = [Z; 3];

pub fn qs_red(u: &[Z]) -> R<Qs> {
    if u.len() != 3 {
        bail!("ratQ:qsShape", "a Q(t) element is [n1 n2 D]; got length {}.", u.len());
    }
    if u[2].is_zero() {
        bail!("ratQ:qsZeroDen", "a Q(t) element cannot have a zero denominator.");
    }
    let mut u = [u[0].clone(), u[1].clone(), u[2].clone()];
    if u[2].is_negative() {
        for x in u.iter_mut() {
            *x = -&*x;
        }
    }
    let g = u[0].gcd(&u[1]).gcd(&u[2]);
    if g.is_positive() {
        for x in u.iter_mut() {
            *x = &*x / &g;
        }
    }
    Ok(u)
}

pub fn qs_add(u: &Qs, v: &Qs) -> R<Qs> {
    qs_red(&[&u[0] * &v[2] + &v[0] * &u[2], &u[1] * &v[2] + &v[1] * &u[2], &u[2] * &v[2]])
}

pub fn qs_mul(u: &Qs, v: &Qs, q: &[Z]) -> R<Qs> {
    if q.len() == 2 {
        if !u[1].is_zero() || !v[1].is_zero() {
            bail!("ratQ:qsField", "a linear Q means t does not occur, but an operand has a t term.");
        }
        return qs_red(&[&u[0] * &v[0], zi(0), &u[2] * &v[2]]);
    }
    let (al, be, ga) = (&q[0], &q[1], &q[2]);
    qs_red(&[
        al * &u[0] * &v[0] - ga * &u[1] * &v[1],
        al * (&u[0] * &v[1] + &u[1] * &v[0]) - be * &u[1] * &v[1],
        al * &u[2] * &v[2],
    ])
}

pub fn qs_sign(u: &Qs, q: &[Z], iv: &(Q, Q)) -> R<i32> {
    let u = qs_red(u)?;
    if u[1].is_zero() {
        return Ok(sign(&u[0]));
    }
    if q.len() == 2 {
        bail!("ratQ:qsField", "a linear Q means t does not occur, but u has a t term.");
    }
    sign_at(&poly_trim(&[u[1].clone(), u[0].clone()]), q, iv)
}

pub type Qs2 = [Z; 5];

pub fn qs2_red(u: &[Z]) -> R<Qs2> {
    if u.len() != 5 {
        bail!("ratQ:qs2Shape", "a Q(t1,t2) element is [n1 n2 n3 n4 D]; got length {}.", u.len());
    }
    if u[4].is_zero() {
        bail!("ratQ:qs2ZeroDen", "a Q(t1,t2) element cannot have a zero denominator.");
    }
    let mut w: Qs2 = [u[0].clone(), u[1].clone(), u[2].clone(), u[3].clone(), u[4].clone()];
    if w[4].is_negative() {
        for x in w.iter_mut() {
            *x = -&*x;
        }
    }
    let g = gcd_all(&w);
    if g.is_positive() {
        for x in w.iter_mut() {
            *x = &*x / &g;
        }
    }
    Ok(w)
}

pub fn qs2_add(u: &Qs2, v: &Qs2) -> R<Qs2> {
    qs2_red(&[
        &u[0] * &v[4] + &v[0] * &u[4],
        &u[1] * &v[4] + &v[1] * &u[4],
        &u[2] * &v[4] + &v[2] * &u[4],
        &u[3] * &v[4] + &v[3] * &u[4],
        &u[4] * &v[4],
    ])
}

pub fn qs2_mul(u: &Qs2, v: &Qs2, q1: &[Z], q2: &[Z]) -> R<Qs2> {
    let (al1, be1, ga1) = (&q1[0], &q1[1], &q1[2]);
    let (al2, be2, ga2) = (&q2[0], &q2[1], &q2[2]);
    let a = &u[0..4];
    let b = &v[0..4];
    let c00 = &a[0] * &b[0];
    let c10 = &a[0] * &b[1] + &a[1] * &b[0];
    let c01 = &a[0] * &b[2] + &a[2] * &b[0];
    let c11 = &a[0] * &b[3] + &a[3] * &b[0] + &a[1] * &b[2] + &a[2] * &b[1];
    let c20 = &a[1] * &b[1];
    let c02 = &a[2] * &b[2];
    let c21 = &a[1] * &b[3] + &a[3] * &b[1];
    let c12 = &a[2] * &b[3] + &a[3] * &b[2];
    let c22 = &a[3] * &b[3];
    let n1 = al1 * al2 * &c00 - al2 * ga1 * &c20 - al1 * ga2 * &c02 + ga1 * ga2 * &c22;
    let n2 = al1 * al2 * &c10 - al2 * be1 * &c20 - al1 * ga2 * &c12 + be1 * ga2 * &c22;
    let n3 = al1 * al2 * &c01 - al1 * be2 * &c02 - al2 * ga1 * &c21 + ga1 * be2 * &c22;
    let n4 = al1 * al2 * &c11 - al2 * be1 * &c21 - al1 * be2 * &c12 + be1 * be2 * &c22;
    qs2_red(&[n1, n2, n3, n4, al1 * al2 * &u[4] * &v[4]])
}

fn qs2_interval(u: &Qs2, i1: &(Q, Q), i2: &(Q, Q)) -> (i32, bool) {
    let c: Vec<f64> = u[0..4].iter().map(to_f64).collect();
    let xs = [q_to_f64(&i1.0), q_to_f64(&i1.1)];
    let ys = [q_to_f64(&i2.0), q_to_f64(&i2.1)];
    let mut v = Vec::with_capacity(4);
    for a in xs {
        for b in ys {
            v.push(c[0] + c[1] * a + c[2] * b + c[3] * a * b);
        }
    }
    let eps = 2f64.powi(-52);
    let cm = c.iter().fold(1.0f64, |m, x| m.max(x.abs()));
    let xm = xs.iter().chain(ys.iter()).fold(1.0f64, |m, x| m.max(x.abs()));
    let tol = 8.0 * eps * cm * xm * xm;
    let mn = v.iter().cloned().fold(f64::INFINITY, f64::min);
    let mx = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if mn > tol {
        return (1, true);
    }
    if mx < -tol {
        return (-1, true);
    }
    (0, false)
}

fn qs2_bisect(q: &[Z], lo: &Q, hi: &Q) -> (Q, Q) {
    let mid = (lo + hi) / qi(2);
    let sm = signq(&polyval_q(q, &mid));
    let sl = signq(&polyval_q(q, lo));
    if sm == 0 {
        return (mid.clone(), mid);
    }
    if sm == sl {
        return (mid, hi.clone());
    }
    (lo.clone(), mid)
}

/// (s, ok): the sign of an element of Q(t1, t2); ok false when it could not be decided.
pub fn qs2_sign(u: &Qs2, q1: &[Z], iv1: &(Q, Q), q2: &[Z], iv2: &(Q, Q)) -> R<(i32, bool)> {
    let u = qs2_red(u)?;
    let d1 = &q1[1] * &q1[1] - zi(4) * &q1[0] * &q1[2];
    let d2 = &q2[1] * &q2[1] - zi(4) * &q2[0] * &q2[2];
    if d1.is_positive() && d2.is_positive() {
        let p = &d1 * &d2;
        let r = p.sqrt();
        if &r * &r == p {
            return Ok((0, false));
        }
    }
    if all_zero(&u[1..4]) {
        return Ok((sign(&u[0]), true));
    }
    if all_zero(&u[0..4]) {
        return Ok((0, true));
    }
    let (mut lo1, mut hi1) = iv1.clone();
    let (mut lo2, mut hi2) = iv2.clone();
    for _ in 0..60 {
        let (s, decided) = qs2_interval(&u, &(lo1.clone(), hi1.clone()), &(lo2.clone(), hi2.clone()));
        if decided {
            return Ok((s, true));
        }
        let a = qs2_bisect(q1, &lo1, &hi1);
        lo1 = a.0;
        hi1 = a.1;
        let b = qs2_bisect(q2, &lo2, &hi2);
        lo2 = b.0;
        hi2 = b.1;
    }
    Ok((0, false))
}

/// Exact rational of a double (Python Fraction(x)).
pub fn q_from_f64(x: f64) -> Q {
    Q::from_float(x).expect("finite double")
}

// ------------------------------------------------------------------------------------------
// change of variables in a face
// ------------------------------------------------------------------------------------------
/// The face n/d composed with x = (Mn*u + tn)/md, exactly: returns (n2, d2), canonical.
pub fn subst_affine(n: &[Z], d: &Z, mn: &[Vec<Z>], tn: &[Z], md: &Z) -> R<(Vec<Z>, Z)> {
    if md.is_zero() {
        bail!("ratQ:badDenominator", "the map denominator must be a nonzero scalar.");
    }
    if det_exact(mn)?.is_zero() {
        bail!("ratQ:singular", "the change of variables is singular.");
    }
    let (mut mn, mut tn, mut md) = (mn.to_vec(), tn.to_vec(), md.clone());
    if md.is_negative() {
        mn = mn.iter().map(|r| neg(r)).collect();
        tn = neg(&tn);
        md = -md;
    }
    use std::collections::BTreeMap;
    type P2 = BTreeMap<(usize, usize), Z>;
    let mut p: Vec<((usize, usize), Z)> = vec![
        ((3, 0), n[0].clone()),
        ((2, 1), zi(3) * &n[1]),
        ((1, 2), zi(3) * &n[2]),
        ((0, 3), n[3].clone()),
        ((2, 0), zi(3) * &n[4]),
        ((1, 1), zi(6) * &n[5]),
        ((0, 2), zi(3) * &n[6]),
        ((1, 0), zi(6) * &n[7]),
        ((0, 1), zi(6) * &n[8]),
        ((0, 0), zi(6) * &n[9]),
    ];
    p.retain(|(_, c)| !c.is_zero());
    let xm: P2 = [((1, 0), mn[0][0].clone()), ((0, 1), mn[0][1].clone()), ((0, 0), tn[0].clone())].into_iter().collect();
    let ym: P2 = [((1, 0), mn[1][0].clone()), ((0, 1), mn[1][1].clone()), ((0, 0), tn[1].clone())].into_iter().collect();
    fn pmul(a: &P2, b: &P2) -> P2 {
        let mut out = P2::new();
        for (ea, ca) in a {
            if ca.is_zero() {
                continue;
            }
            for (eb, cb) in b {
                if cb.is_zero() {
                    continue;
                }
                let e = (ea.0 + eb.0, ea.1 + eb.1);
                *out.entry(e).or_insert_with(Z::zero) += ca * cb;
            }
        }
        out
    }
    fn ppow(a: &P2, k: usize) -> P2 {
        let mut r: P2 = [((0, 0), Z::one())].into_iter().collect();
        for _ in 0..k {
            r = pmul(&r, a);
        }
        r
    }
    let mut res = P2::new();
    for ((i, j), cij) in &p {
        let term = pmul(&ppow(&xm, *i), &ppow(&ym, *j));
        let lift = num_traits::pow(md.clone(), 3 - i - j);
        for (e, c) in term {
            *res.entry(e).or_insert_with(Z::zero) += cij * &lift * c;
        }
    }
    let r = |i: usize, j: usize| res.get(&(i, j)).cloned().unwrap_or_else(Z::zero);
    let den = zi(6) * d * num_traits::pow(md.clone(), 3);
    let n2 = vec![
        zi(6) * r(3, 0),
        zi(2) * r(2, 1),
        zi(2) * r(1, 2),
        zi(6) * r(0, 3),
        zi(2) * r(2, 0),
        r(1, 1),
        zi(2) * r(0, 2),
        r(1, 0),
        r(0, 1),
        r(0, 0),
    ];
    canon(&n2, &den)
}

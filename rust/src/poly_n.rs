//! poly_n -- exact multivariate integer polynomials, sparse: {exponent tuple: coefficient}
//! (port of poly_n.py). A BTreeMap keeps the terms in lex order, so the lex-greatest term is the
//! last one.

use crate::bail;
use crate::errors::R;
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeMap;

pub type PolyN = BTreeMap<Vec<u32>, Z>;

pub fn mono(e: &[u32], c: Z) -> PolyN {
    let mut p = PolyN::new();
    if !c.is_zero() {
        p.insert(e.to_vec(), c);
    }
    p
}

pub fn add(a: &PolyN, b: &PolyN) -> PolyN {
    let mut out = a.clone();
    for (e, c) in b {
        let v = out.get(e).cloned().unwrap_or_else(Z::zero) + c;
        if v.is_zero() {
            out.remove(e);
        } else {
            out.insert(e.clone(), v);
        }
    }
    out
}

pub fn negp(a: &PolyN) -> PolyN {
    a.iter().map(|(e, c)| (e.clone(), -c)).collect()
}

pub fn subp(a: &PolyN, b: &PolyN) -> PolyN {
    add(a, &negp(b))
}

pub fn scale(a: &PolyN, k: &Z) -> PolyN {
    if k.is_zero() {
        return PolyN::new();
    }
    a.iter().map(|(e, c)| (e.clone(), c * k)).collect()
}

pub fn mul(a: &PolyN, b: &PolyN) -> PolyN {
    let mut out = PolyN::new();
    for (ea, ca) in a {
        for (eb, cb) in b {
            let e: Vec<u32> = ea.iter().zip(eb.iter()).map(|(x, y)| x + y).collect();
            *out.entry(e).or_insert_with(Z::zero) += ca * cb;
        }
    }
    out.retain(|_, c| !c.is_zero());
    out
}

pub fn content(a: &PolyN) -> Z {
    let mut g = Z::zero();
    for c in a.values() {
        g = g.gcd(c);
        if g.is_one() {
            return g;
        }
    }
    g
}

/// Divide out the (positive) gcd of the coefficients; the SIGN is left alone.
pub fn content_free(a: &PolyN) -> PolyN {
    let g = content(a);
    if g > Z::one() {
        a.iter().map(|(e, c)| (e.clone(), c / &g)).collect()
    } else {
        a.clone()
    }
}

pub fn degrees(a: &PolyN, nv: usize) -> Vec<u32> {
    let mut d = vec![0u32; nv];
    for e in a.keys() {
        for i in 0..nv {
            d[i] = d[i].max(e[i]);
        }
    }
    d
}

pub fn max_coeff(a: &PolyN) -> Z {
    a.values().map(|c| c.abs()).max().unwrap_or_else(Z::zero)
}

/// P / D when the division is exact; ratQ:notDivisible otherwise.
pub fn div_exact(p: &PolyN, d: &PolyN) -> R<PolyN> {
    let (ed, cd) = match d.iter().next_back() {
        None => bail!("ratQ:divideByZero", "the divisor is the zero polynomial."),
        Some((e, c)) => (e.clone(), c.clone()),
    };
    let mut r = p.clone();
    let mut q = PolyN::new();
    let mut guard = 0;
    while let Some((er, cr)) = r.iter().next_back().map(|(e, c)| (e.clone(), c.clone())) {
        guard += 1;
        if guard > 100000 {
            bail!("ratQ:divideNoConverge", "exact division did not terminate.");
        }
        if er.iter().zip(ed.iter()).any(|(x, y)| x < y) || !cr.mod_floor(&cd).is_zero() {
            bail!("ratQ:notDivisible", "the division is not exact.");
        }
        let eq: Vec<u32> = er.iter().zip(ed.iter()).map(|(x, y)| x - y).collect();
        let t = mono(&eq, cr.div_floor(&cd));
        q = add(&q, &t);
        r = subp(&r, &mul(&t, d));
    }
    Ok(q)
}

/// {k: polynomial in the other variables} -- the coefficients of variable d.
pub fn slices(a: &PolyN, d: usize) -> BTreeMap<u32, PolyN> {
    let mut out: BTreeMap<u32, PolyN> = BTreeMap::new();
    for (e, c) in a {
        let mut rest = e.clone();
        let k = rest.remove(d);
        out.entry(k).or_default().insert(rest, c.clone());
    }
    out
}

pub fn from_slices(sl: &BTreeMap<u32, PolyN>, d: usize) -> PolyN {
    let mut out = PolyN::new();
    for (k, p) in sl {
        for (e, c) in p {
            let mut e2 = e.clone();
            e2.insert(d, *k);
            out.insert(e2, c.clone());
        }
    }
    out
}

//! isolated_root_sign -- the exact sign of N(alpha) where alpha is the root of R isolated by
//! iv = (lo, hi] (port of isolated_root_sign.py). R need not be squarefree; rational roots are
//! met exactly.

use crate::bail;
use crate::errors::R;
use crate::rat_q::*;
use num_integer::Integer;
use num_traits::Signed;

fn strip(p: &[Z]) -> Vec<Z> {
    poly_trim(p)
}

/// MATLAB ratQ.bigPrem: a * |lc(b)|^(da-db+1) mod b (no content removal).
fn prem(a: &[Z], b: &[Z]) -> Vec<Z> {
    if a.len() < b.len() {
        return a.to_vec();
    }
    let da = a.len() - 1;
    let db = b.len() - 1;
    let lc = &b[0];
    let f = num_traits::pow(lc.abs(), da - db + 1);
    let mut r: Vec<Z> = a.iter().map(|c| c * &f).collect();
    while !r.is_empty() && r.len() > db {
        let q = r[0].div_floor(lc);
        for j in 0..b.len() {
            r[j] -= &q * &b[j];
        }
        r.remove(0);
        while !r.is_empty() && r[0].is_zero() {
            r.remove(0);
        }
        if r.is_empty() {
            return r;
        }
    }
    r
}

fn bgcd(a: &[Z], b: &[Z]) -> Vec<Z> {
    let mut a = poly_primitive(a);
    let mut b = poly_primitive(b);
    loop {
        let r = strip(&prem(&a, &b));
        if r.is_empty() {
            return b;
        }
        if r.len() == 1 {
            return vec![zi(1)];
        }
        a = b;
        b = poly_primitive(&r);
    }
}

use num_traits::Zero;

pub fn isolated_root_sign(n: &[Z], r: &[Z], iv: &(Q, Q)) -> R<i32> {
    let n = strip(n);
    if n.is_empty() {
        return Ok(0);
    }
    let rb = strip(r);
    let mut lo = iv.0.clone();
    let mut hi = iv.1.clone();
    let two = qi(2);
    if signq(&polyval_q(&rb, &hi)) == 0 {
        return Ok(signq(&polyval_q(&n, &hi)));
    }
    while signq(&polyval_q(&rb, &lo)) == 0 {
        let mid = (&lo + &hi) / &two;
        if signq(&polyval_q(&rb, &mid)) == 0 {
            return Ok(signq(&polyval_q(&n, &mid)));
        }
        if count_roots_in_q(&rb, &mid, &hi) == 1 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let rr = strip(&prem(&n, &rb));
    if rr.is_empty() {
        return Ok(0);
    }
    if rr.len() == 1 {
        return Ok(sign(&rr[0]));
    }
    let g = bgcd(&rb, &rr);
    if g.len() > 1 && count_roots_in_q(&g, &lo, &hi) > 0 {
        return Ok(0);
    }
    for _ in 0..400 {
        if count_roots_in_q(&rr, &lo, &hi) == 0 && signq(&polyval_q(&rr, &lo)) != 0 {
            return Ok(signq(&polyval_q(&rr, &((&lo + &hi) / &two))));
        }
        let mid = (&lo + &hi) / &two;
        if signq(&polyval_q(&rb, &mid)) == 0 {
            return Ok(signq(&polyval_q(&rr, &mid)));
        }
        if count_roots_in_q(&rb, &lo, &mid) == 1 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    bail!("isolatedRootSign:noConverge", "400 bisections did not separate the root from the zeros of the reduced polynomial.")
}

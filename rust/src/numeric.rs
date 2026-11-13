//! numeric -- the few double-precision primitives the float screens rely on (port of numeric.py):
//! polynomial roots (numpy.roots: eigenvalues of the companion matrix), linspace, round, and the
//! small dense LU solves numpy/LAPACK does for 3x3 systems.

pub const EPS: f64 = 2.220446049250313e-16;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cx {
    pub re: f64,
    pub im: f64,
}

/// numpy.roots(p): complex roots of the polynomial p (highest degree first).
pub fn roots(p: &[f64]) -> Vec<Cx> {
    let nz: Vec<usize> = p.iter().enumerate().filter(|(_, x)| **x != 0.0).map(|(i, _)| i).collect();
    if nz.is_empty() {
        return vec![];
    }
    let first = nz[0];
    let last = *nz.last().unwrap();
    let trailing = p.len() - last - 1;
    let q = &p[first..=last];
    let n = q.len() - 1;
    let mut out = Vec::new();
    if n >= 1 {
        // companion matrix, numpy layout: first row -q[1:]/q[0], ones on the subdiagonal
        let mut a = vec![vec![0.0f64; n]; n];
        for j in 0..n {
            a[0][j] = -q[j + 1] / q[0];
        }
        for i in 1..n {
            a[i][i - 1] = 1.0;
        }
        out = eig_hessenberg(a);
    }
    for _ in 0..trailing {
        out.push(Cx { re: 0.0, im: 0.0 });
    }
    out
}

/// Eigenvalues of an upper Hessenberg matrix (balanced first), by the shifted QR algorithm
/// (EISPACK balanc + hqr). The companion matrix is already Hessenberg.
fn eig_hessenberg(mut a: Vec<Vec<f64>>) -> Vec<Cx> {
    let n = a.len();
    balance(&mut a);
    let mut wr = vec![0.0; n];
    let mut wi = vec![0.0; n];
    hqr(&mut a, &mut wr, &mut wi);
    (0..n).map(|i| Cx { re: wr[i], im: wi[i] }).collect()
}

fn balance(a: &mut [Vec<f64>]) {
    let n = a.len();
    let radix = 2.0f64;
    let sqrdx = radix * radix;
    let mut last = false;
    while !last {
        last = true;
        for i in 0..n {
            let mut r = 0.0;
            let mut c = 0.0;
            for j in 0..n {
                if j != i {
                    c += a[j][i].abs();
                    r += a[i][j].abs();
                }
            }
            if c != 0.0 && r != 0.0 {
                let mut g = r / radix;
                let mut f = 1.0;
                let s = c + r;
                while c < g {
                    f *= radix;
                    c *= sqrdx;
                }
                g = r * radix;
                while c > g {
                    f /= radix;
                    c /= sqrdx;
                }
                if (c + r) / f < 0.95 * s {
                    last = false;
                    g = 1.0 / f;
                    for j in 0..n {
                        a[i][j] *= g;
                    }
                    for j in 0..n {
                        a[j][i] *= f;
                    }
                }
            }
        }
    }
}

fn hqr(a: &mut [Vec<f64>], wr: &mut [f64], wi: &mut [f64]) {
    let n = a.len() as i64;
    let mut anorm = 0.0;
    for i in 0..n as usize {
        for j in (if i > 0 { i - 1 } else { 0 })..n as usize {
            anorm += a[i][j].abs();
        }
    }
    let mut nn = n - 1;
    let mut t = 0.0;
    while nn >= 0 {
        let mut its = 0;
        let mut l;
        loop {
            l = nn;
            while l >= 1 {
                let s = a[(l - 1) as usize][(l - 1) as usize].abs() + a[l as usize][l as usize].abs();
                let s = if s == 0.0 { anorm } else { s };
                if a[l as usize][(l - 1) as usize].abs() + s == s {
                    a[l as usize][(l - 1) as usize] = 0.0;
                    break;
                }
                l -= 1;
            }
            let nnu = nn as usize;
            let x = a[nnu][nnu];
            if l == nn {
                wr[nnu] = x + t;
                wi[nnu] = 0.0;
                nn -= 1;
                break;
            }
            let y = a[nnu - 1][nnu - 1];
            let w = a[nnu][nnu - 1] * a[nnu - 1][nnu];
            if l == nn - 1 {
                let p = 0.5 * (y - x);
                let q = p * p + w;
                let z = q.abs().sqrt();
                let xx = x + t;
                if q >= 0.0 {
                    let z = p + if p >= 0.0 { z.abs() } else { -z.abs() };
                    wr[nnu - 1] = xx + z;
                    wr[nnu] = wr[nnu - 1];
                    if z != 0.0 {
                        wr[nnu] = xx - w / z;
                    }
                    wi[nnu - 1] = 0.0;
                    wi[nnu] = 0.0;
                } else {
                    wr[nnu - 1] = xx + p;
                    wr[nnu] = xx + p;
                    wi[nnu - 1] = -z;
                    wi[nnu] = z;
                }
                nn -= 2;
                break;
            }
            if its == 60 {
                // no convergence: report what we have (NaN marks the rest)
                for k in 0..=nnu {
                    wr[k] = f64::NAN;
                    wi[k] = 0.0;
                }
                return;
            }
            let (mut xx, mut yy, mut ww) = (x, y, w);
            if its == 10 || its == 20 {
                t += xx;
                for i in 0..=nnu {
                    a[i][i] -= xx;
                }
                let s = a[nnu][nnu - 1].abs() + a[nnu - 1][nnu - 2].abs();
                xx = 0.75 * s;
                yy = xx;
                ww = -0.4375 * s * s;
            }
            its += 1;
            let mut m = nn - 2;
            let (mut p, mut q, mut r);
            loop {
                let mu = m as usize;
                let z = a[mu][mu];
                let rr = xx - z;
                let ss = yy - z;
                p = (rr * ss - ww) / a[mu + 1][mu] + a[mu][mu + 1];
                q = a[mu + 1][mu + 1] - z - rr - ss;
                r = a[mu + 2][mu + 1];
                let s = p.abs() + q.abs() + r.abs();
                p /= s;
                q /= s;
                r /= s;
                if m == l {
                    break;
                }
                let u = a[mu][mu - 1].abs() * (q.abs() + r.abs());
                let v = p.abs() * (a[mu - 1][mu - 1].abs() + z.abs() + a[mu + 1][mu + 1].abs());
                if u + v == v {
                    break;
                }
                m -= 1;
            }
            for i in (m + 2) as usize..=nnu {
                a[i][i - 2] = 0.0;
                if i != (m + 2) as usize {
                    a[i][i - 3] = 0.0;
                }
            }
            let mut k = m;
            while k <= nn - 1 {
                let ku = k as usize;
                if k != m {
                    p = a[ku][ku - 1];
                    q = a[ku + 1][ku - 1];
                    r = 0.0;
                    if k != nn - 1 {
                        r = a[ku + 2][ku - 1];
                    }
                    let xa = p.abs() + q.abs() + r.abs();
                    if xa != 0.0 {
                        p /= xa;
                        q /= xa;
                        r /= xa;
                    }
                    xx = xa;
                }
                let s0 = (p * p + q * q + r * r).sqrt();
                let s = if p >= 0.0 { s0 } else { -s0 };
                if s != 0.0 {
                    if k == m {
                        if l != m {
                            a[ku][ku - 1] = -a[ku][ku - 1];
                        }
                    } else {
                        a[ku][ku - 1] = -s * xx;
                    }
                    p += s;
                    xx = p / s;
                    yy = q / s;
                    let z = r / s;
                    q /= p;
                    r /= p;
                    for j in ku..=nnu {
                        let mut pp = a[ku][j] + q * a[ku + 1][j];
                        if k != nn - 1 {
                            pp += r * a[ku + 2][j];
                            a[ku + 2][j] -= pp * z;
                        }
                        a[ku + 1][j] -= pp * yy;
                        a[ku][j] -= pp * xx;
                    }
                    let mmin = if nn < k + 3 { nn } else { k + 3 } as usize;
                    for i in l as usize..=mmin {
                        let mut pp = xx * a[i][ku] + yy * a[i][ku + 1];
                        if k != nn - 1 {
                            pp += z * a[i][ku + 2];
                            a[i][ku + 2] -= pp * r;
                        }
                        a[i][ku + 1] -= pp * q;
                        a[i][ku] -= pp;
                    }
                }
                k += 1;
            }
            if l >= nn - 1 {
                continue;
            }
        }
    }
}

/// MATLAB/Python round(x): halves away from zero.
pub fn round_half_away(x: f64) -> i64 {
    ((x.abs() + 0.5).floor() * if x >= 0.0 { 1.0 } else { -1.0 }) as i64
}

/// MATLAB linspace: d1 + (0:n-1).*(d2-d1)./(n-1), endpoints exact.
pub fn linspace(d1: f64, d2: f64, n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![d2];
    }
    let n1 = (n - 1) as f64;
    let mut y: Vec<f64> = (0..n).map(|i| d1 + (i as f64 * (d2 - d1)) / n1).collect();
    if d1 == d2 {
        y = vec![d1; n];
    } else {
        y[0] = d1;
        y[n - 1] = d2;
    }
    y
}

pub fn norm2(v: [f64; 2]) -> f64 {
    v[0].hypot(v[1])
}

/// LU with partial pivoting of a 3x3 (LAPACK dgetrf); returns (LU, piv, sign) or None if singular.
fn lu3(m: &[[f64; 3]; 3]) -> ([[f64; 3]; 3], [usize; 3], f64, bool) {
    let mut a = *m;
    let mut piv = [0usize, 1, 2];
    let mut sgn = 1.0;
    let mut singular = false;
    for k in 0..3 {
        let mut p = k;
        let mut mx = a[k][k].abs();
        for i in (k + 1)..3 {
            if a[i][k].abs() > mx {
                mx = a[i][k].abs();
                p = i;
            }
        }
        if p != k {
            a.swap(p, k);
            piv.swap(p, k);
            sgn = -sgn;
        }
        if a[k][k] == 0.0 {
            singular = true;
            continue;
        }
        for i in (k + 1)..3 {
            a[i][k] /= a[k][k];
            let f = a[i][k];
            for j in (k + 1)..3 {
                a[i][j] -= f * a[k][j];
            }
        }
    }
    (a, piv, sgn, singular)
}

/// numpy.linalg.det of a 3x3.
pub fn det3(m: &[[f64; 3]; 3]) -> f64 {
    let (a, _, s, _) = lu3(m);
    s * a[0][0] * a[1][1] * a[2][2]
}

/// numpy.linalg.solve of a 3x3 (None if exactly singular).
pub fn solve3(m: &[[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let (a, piv, _, singular) = lu3(m);
    if singular {
        return None;
    }
    let mut y = [b[piv[0]], b[piv[1]], b[piv[2]]];
    for i in 0..3 {
        for j in 0..i {
            y[i] -= a[i][j] * y[j];
        }
    }
    for i in (0..3).rev() {
        for j in (i + 1)..3 {
            y[i] -= a[i][j] * y[j];
        }
        y[i] /= a[i][i];
    }
    Some(y)
}

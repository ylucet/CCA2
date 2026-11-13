//! flint -- the ONE place the engine may call FLINT (LGPL), as an optional accelerator, exactly
//! as python/cca2/_poly.py does: the odd-multiplicity part of a bivariate integer polynomial
//! (Yun's square-free decomposition) for algRowFails. Without FLINT (or with CCA2_FLINT=0) the
//! pure-Rust path is used and the answers are the same.

use crate::biv::M2;
use crate::rat_q::Z;

pub fn available() -> bool {
    imp::available()
}

/// The product of the square-free factors of ODD multiplicity of S (dense, S[i][j] x^i y^j),
/// primitive; sign left to the caller. None when FLINT is not in use.
pub fn odd_part(s: &M2) -> Option<M2> {
    imp::odd_part(s)
}

#[cfg(have_flint)]
mod imp {
    use super::*;
    use num_traits::{One, Zero};
    use std::ffi::{CStr, CString};
    use std::os::raw::{c_char, c_int, c_long, c_ulong, c_void};

    #[repr(C, align(16))]
    struct Buf([u8; 4096]);

    extern "C" {
        fn fmpz_mpoly_ctx_init(ctx: *mut c_void, nvars: c_long, ord: c_int);
        fn fmpz_mpoly_ctx_clear(ctx: *mut c_void);
        fn fmpz_mpoly_init(a: *mut c_void, ctx: *const c_void);
        fn fmpz_mpoly_clear(a: *mut c_void, ctx: *const c_void);
        fn fmpz_mpoly_set_coeff_fmpz_ui(a: *mut c_void, c: *const c_long, exp: *const c_ulong, ctx: *const c_void);
        fn fmpz_mpoly_factor_init(f: *mut c_void, ctx: *const c_void);
        fn fmpz_mpoly_factor_clear(f: *mut c_void, ctx: *const c_void);
        fn fmpz_mpoly_factor_squarefree(f: *mut c_void, a: *const c_void, ctx: *const c_void) -> c_int;
        fn fmpz_mpoly_factor_length(f: *const c_void, ctx: *const c_void) -> c_long;
        fn fmpz_mpoly_factor_get_base(b: *mut c_void, f: *const c_void, i: c_long, ctx: *const c_void);
        fn fmpz_mpoly_factor_get_exp_si(f: *mut c_void, i: c_long, ctx: *const c_void) -> c_long;
        fn fmpz_mpoly_length(a: *const c_void, ctx: *const c_void) -> c_long;
        fn fmpz_mpoly_get_term_coeff_fmpz(c: *mut c_long, a: *const c_void, i: c_long, ctx: *const c_void);
        fn fmpz_mpoly_get_term_exp_ui(exp: *mut c_ulong, a: *const c_void, i: c_long, ctx: *const c_void);
        fn fmpz_mpoly_mul(a: *mut c_void, b: *const c_void, c: *const c_void, ctx: *const c_void);
        fn fmpz_mpoly_one(a: *mut c_void, ctx: *const c_void);
        fn fmpz_mpoly_set(a: *mut c_void, b: *const c_void, ctx: *const c_void);
        fn fmpz_init(f: *mut c_long);
        fn fmpz_clear(f: *mut c_long);
        fn fmpz_set_str(f: *mut c_long, s: *const c_char, b: c_int) -> c_int;
        fn fmpz_get_str(s: *mut c_char, b: c_int, f: *const c_long) -> *mut c_char;
        fn flint_free(p: *mut c_void);
    }

    pub fn available() -> bool {
        static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
        *ON.get_or_init(|| std::env::var("CCA2_FLINT").map(|v| v != "0").unwrap_or(true))
    }

    fn boxed() -> Box<Buf> {
        Box::new(Buf([0u8; 4096]))
    }

    pub fn odd_part(s: &M2) -> Option<M2> {
        if !available() {
            return None;
        }
        if s.iter().all(|r| r.iter().all(|x| x.is_zero())) {
            return None;
        }
        unsafe {
            let mut ctx = boxed();
            let c = ctx.0.as_mut_ptr() as *mut c_void;
            fmpz_mpoly_ctx_init(c, 2, 0);
            let mut pa = boxed();
            let a = pa.0.as_mut_ptr() as *mut c_void;
            fmpz_mpoly_init(a, c);
            let mut z: c_long = 0;
            fmpz_init(&mut z);
            for (i, r) in s.iter().enumerate() {
                for (j, v) in r.iter().enumerate() {
                    if v.is_zero() {
                        continue;
                    }
                    let cs = CString::new(v.to_str_radix(16)).unwrap();
                    fmpz_set_str(&mut z, cs.as_ptr(), 16);
                    let e: [c_ulong; 2] = [i as c_ulong, j as c_ulong];
                    fmpz_mpoly_set_coeff_fmpz_ui(a, &z, e.as_ptr(), c);
                }
            }
            let mut pf = boxed();
            let f = pf.0.as_mut_ptr() as *mut c_void;
            fmpz_mpoly_factor_init(f, c);
            let ok = fmpz_mpoly_factor_squarefree(f, a, c);
            let mut po = boxed();
            let o = po.0.as_mut_ptr() as *mut c_void;
            fmpz_mpoly_init(o, c);
            fmpz_mpoly_one(o, c);
            let mut pb = boxed();
            let b = pb.0.as_mut_ptr() as *mut c_void;
            fmpz_mpoly_init(b, c);
            let mut pt = boxed();
            let t = pt.0.as_mut_ptr() as *mut c_void;
            fmpz_mpoly_init(t, c);
            let mut out = None;
            if ok != 0 {
                let n = fmpz_mpoly_factor_length(f, c);
                for i in 0..n {
                    let m = fmpz_mpoly_factor_get_exp_si(f, i, c);
                    if m % 2 == 1 {
                        fmpz_mpoly_factor_get_base(b, f, i, c);
                        fmpz_mpoly_mul(t, o, b, c);
                        fmpz_mpoly_set(o, t, c);
                    }
                }
                let len = fmpz_mpoly_length(o, c);
                let mut terms = vec![];
                for i in 0..len {
                    let mut e: [c_ulong; 2] = [0, 0];
                    fmpz_mpoly_get_term_exp_ui(e.as_mut_ptr(), o, i, c);
                    fmpz_mpoly_get_term_coeff_fmpz(&mut z, o, i, c);
                    let p = fmpz_get_str(std::ptr::null_mut(), 16, &z);
                    let sv = CStr::from_ptr(p).to_str().unwrap().to_string();
                    flint_free(p as *mut c_void);
                    terms.push((e[0] as usize, e[1] as usize, Z::parse_bytes(sv.as_bytes(), 16).unwrap()));
                }
                let n1 = terms.iter().map(|t| t.0).max().unwrap_or(0) + 1;
                let n2 = terms.iter().map(|t| t.1).max().unwrap_or(0) + 1;
                let mut mm = vec![vec![Z::zero(); n2]; n1];
                for (i, j, v) in terms {
                    mm[i][j] = v;
                }
                if mm.iter().all(|r| r.iter().all(|x| x.is_zero())) {
                    mm[0][0] = Z::one();
                }
                out = Some(mm);
            }
            fmpz_mpoly_clear(t, c);
            fmpz_mpoly_clear(b, c);
            fmpz_mpoly_clear(o, c);
            fmpz_mpoly_factor_clear(f, c);
            fmpz_clear(&mut z);
            fmpz_mpoly_clear(a, c);
            fmpz_mpoly_ctx_clear(c);
            out
        }
    }
}

#[cfg(not(have_flint))]
mod imp {
    use super::*;
    pub fn available() -> bool {
        false
    }
    pub fn odd_part(_s: &M2) -> Option<M2> {
        None
    }
}

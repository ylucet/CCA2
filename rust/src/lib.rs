//! cca2 -- exact conjugate and biconjugate of bivariate PLQ functions (Rust port of the Python
//! port in CCA2/python, itself a port of the MATLAB reference). Module names follow the Python
//! modules one for one; index values stay 1-based in the data, 0 = none, exactly as there.

pub mod errors;
pub mod rat_q;
pub mod numeric;
pub mod probe_data;
pub mod conic_meet;
pub mod qua_pol;
pub mod qua_con;
pub mod cells;
pub mod max_q;
pub mod conj_q;
pub mod io;
pub mod is_convex_q;
pub mod biconj_q;
pub mod poly_n;
pub mod biv;
pub mod alg_branch;
pub mod isolated_root_sign;
pub mod alg_con;
pub mod envelope;
pub mod bridge;
pub mod corner;
pub mod envelope_subdivision;
pub mod envelope_subdivision_alg;
pub mod poly_sign_on;
pub mod alg_row_fails;
pub mod alg_face_below;
pub mod triangle_face;
pub mod biconj_atoms;
pub mod flint;

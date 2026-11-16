//! digest -- run every fixture once and print a hash of the COMPLETE result (Debug text of the
//! returned object: every face, row, coefficient, in order), so two builds can be compared for
//! bit-identical output, not just for kind and face count.
//!
//!     digest fixtures.json > out.tsv         columns: id, kind, nf, err, hash, bytes
//!     CCA2_DIGEST_ONLY=191,196 digest ...    only these ids
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use cca2::io::num;
use cca2::qua_pol::QuaPol;
use serde_json::Value;

fn build(x: &Value) -> QuaPol {
    let rows = |v: &Value| -> Vec<Vec<cca2::rat_q::Num>> { v.as_array().map(|a| a.iter().map(|r| r.as_array().unwrap().iter().map(num).collect()).collect()).unwrap_or_default() };
    let f = rows(&x["f"]);
    let v = rows(&x["V"]);
    if v.is_empty() {
        return QuaPol::whole_rows(f).unwrap();
    }
    let idx = |v: &Value| v.as_f64().unwrap() as usize;
    let e: Vec<[usize; 3]> = x["E"].as_array().map(|a| a.iter().map(|r| [idx(&r[0]), idx(&r[1]), idx(&r[2])]).collect()).unwrap_or_default();
    let ff: Vec<[usize; 2]> = x["F"].as_array().map(|a| a.iter().map(|r| [idx(&r[0]), idx(&r[1])]).collect()).unwrap_or_default();
    QuaPol::new(v, e, f, ff).unwrap()
}


fn main() {
    let a: Vec<String> = std::env::args().collect();
    eprintln!("{}", cca2::flint::status());
    let fix: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let only: Option<Vec<i64>> = std::env::var("CCA2_DIGEST_ONLY").ok().map(|s| s.split(',').map(|x| x.parse().unwrap()).collect());
    let h = std::thread::Builder::new()
        .stack_size(256 << 20)
        .spawn(move || {
            for x in &fix {
                let id = x["id"].as_i64().unwrap();
                if let Some(o) = &only {
                    if !o.contains(&id) {
                        continue;
                    }
                }
                let f = build(x);
                let (kind, nf, err, txt) = if x["fn"].as_str() == Some("conjQ") {
                    match cca2::conj_q::conj_q(&f) {
                        Ok((g, d)) => ("QuaCon".to_string(), g.nf.to_string(), String::new(), format!("{:?}{:?}", g, d)),
                        Err(e) => (String::new(), String::new(), e.id.clone(), format!("{:?}", e.id)),
                    }
                } else {
                    match cca2::biconj_q::biconj_q(&f, &Default::default()) {
                        Ok(r) => (r.kind().to_string(), r.nf().map(|n| n.to_string()).unwrap_or_default(), String::new(), format!("{:?}", r)),
                        Err(e) => (String::new(), String::new(), e.id.clone(), format!("{:?}", e.id)),
                    }
                };
                // OnceLock and OnceCell (before 2026-11-15) print under different names
                let txt = txt.replace("OnceLock", "OnceCell");
                let mut hs = DefaultHasher::new();
                txt.hash(&mut hs);
                println!("{}\t{}\t{}\t{}\t{:016x}\t{}", id, kind, nf, err, hs.finish(), txt.len());
            }
        })
        .unwrap();
    h.join().unwrap();
}

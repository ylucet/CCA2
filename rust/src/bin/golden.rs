//! Replay a golden trace (JSON lines) through the Rust engine and print one JSON result per
//! record: {"err": id} for a refusal, else the QuaCon fields (QuaCon) or {"kind"} plus "evalV"
//! sampled at the record's evalX. Compared against the Python port by tools/compare.py.
#[cfg(feature = "mimalloc")]
#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

use cca2::io::*;
use serde_json::{json, Value};
use std::io::BufRead;

fn main() {
    let path = std::env::args().nth(1).expect("trace path");
    let only: Option<String> = std::env::args().nth(2);
    let f = std::fs::File::open(&path).unwrap();
    for line in std::io::BufReader::new(f).lines() {
        let line = line.unwrap();
        if line.trim().is_empty() {
            continue;
        }
        let rec: Value = serde_json::from_str(&line).unwrap();
        let fname = rec["fn"].as_str().unwrap_or("").to_string();
        if let Some(o) = &only {
            if &fname != o {
                println!("null");
                continue;
            }
        }
        let t0 = std::time::Instant::now();
        let out = run(&fname, &rec);
        let mut out = out;
        out["t"] = json!(t0.elapsed().as_secs_f64());
        println!("{}", out);
    }
}

fn run(fname: &str, rec: &Value) -> Value {
    let p = match qua_pol_from(&rec["in"]) {
        Ok(p) => p,
        Err(e) => return json!({"err": e.id}),
    };
    let x = eval_x(&rec["out"]);
    match fname {
        "conjQ" => match cca2::conj_q::conj_q(&p) {
            Err(e) => json!({"err": e.id, "msg": e.msg}),
            Ok((g, dom)) => {
                let mut j = qua_con_json(&g);
                j["dom"] = dom_json(&dom);
                match g.eval(&x) {
                    Ok((v, _)) => j["evalV"] = Value::Array(v.into_iter().map(fj).collect()),
                    Err(e) => j["evalErr"] = json!(e.id),
                }
                j
            }
        },
        "biconjQ" => match cca2::biconj_q::biconj_q(&p, &Default::default()) {
            Err(e) => json!({"err": e.id, "msg": e.msg}),
            Ok(h) => {
                let mut j = match &h {
                    cca2::alg_con::AlgCon::QuaCon(g) => qua_con_json(g),
                    _ => json!({"kind": h.kind()}),
                };
                if let Some(nf) = h.nf() {
                    j["nf"] = json!(nf);
                }
                match h.eval(&x) {
                    Ok(v) => j["evalV"] = Value::Array(v.into_iter().map(fj).collect()),
                    Err(e) => j["evalErr"] = json!(e.id),
                }
                j
            }
        },
        _ => json!({"err": "rust:notPorted"}),
    }
}

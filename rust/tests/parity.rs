//! Parity with the Python port (the oracle): every conjQ / biconjQ call of the MATLAB golden
//! traces is replayed and compared with what the Python port returned for it (ref/*.py.jsonl,
//! written by tools/py_ref.py): the refusal identifier, or the QuaCon field for field (nf, ne, nv,
//! fN, fD, EcQ, FC, Vname, E, F, dom), the result kind and face count, and the value at every
//! sampled point to 1e-9 relative.
//!
//!     cargo test --release                 # conj + biconj traces (a few seconds)
//!     cargo test --release -- --ignored    # also the slow MATLAB suites (kink, Example 2, ...)
use cca2::alg_con::AlgCon;
use cca2::io::*;
use serde_json::{json, Value};
use std::path::PathBuf;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn lines(p: &PathBuf) -> Vec<Value> {
    std::fs::read_to_string(p).unwrap().lines().filter(|l| !l.trim().is_empty()).map(|l| serde_json::from_str(l).unwrap()).collect()
}

fn as_f(v: &Value) -> f64 {
    match v {
        Value::String(s) => match s.as_str() {
            "inf" => f64::INFINITY,
            "-inf" => f64::NEG_INFINITY,
            _ => f64::NAN,
        },
        _ => numf(v),
    }
}

/// Integer-valued JSON as canonical strings, so 3 and 3.0 compare equal.
fn norm(v: &Value) -> Value {
    match v {
        Value::Array(a) => Value::Array(a.iter().map(norm).collect()),
        Value::Number(_) => json!(zint(v).to_string()),
        Value::String(s) => json!(s.clone()),
        x => x.clone(),
    }
}

fn run(fname: &str, rec: &Value) -> Value {
    let p = match qua_pol_from(&rec["in"]) {
        Ok(p) => p,
        Err(e) => return json!({"err": e.id}),
    };
    let x = eval_x(&rec["out"]);
    let (mut j, h): (Value, Option<AlgCon>) = if fname == "conjQ" {
        match cca2::conj_q::conj_q(&p) {
            Err(e) => return json!({"err": e.id}),
            Ok((g, dom)) => {
                let mut j = qua_con_json(&g);
                j["dom"] = dom_json(&dom);
                (j, Some(AlgCon::QuaCon(g)))
            }
        }
    } else {
        match cca2::biconj_q::biconj_q(&p, &Default::default()) {
            Err(e) => return json!({"err": e.id}),
            Ok(h) => {
                let mut j = match &h {
                    AlgCon::QuaCon(g) => qua_con_json(g),
                    _ => json!({"kind": h.kind()}),
                };
                if let Some(nf) = h.nf() {
                    j["nf"] = json!(nf);
                }
                (j, Some(h))
            }
        }
    };
    if !x.is_empty() {
        let v = h.unwrap().eval(&x).unwrap();
        j["evalV"] = Value::Array(v.into_iter().map(fj).collect());
    }
    j
}

fn check(trace: &str, refname: &str) {
    let recs = lines(&root().join("../python/golden").join(trace));
    let refs = lines(&root().join("ref").join(refname));
    assert_eq!(recs.len(), refs.len(), "{}: trace and reference differ in length", trace);
    let mut bad = vec![];
    for (i, (rec, py)) in recs.iter().zip(refs.iter()).enumerate() {
        if py.is_null() {
            continue;
        }
        let r = run(rec["fn"].as_str().unwrap(), rec);
        let why = if r.get("err").is_some() || py.get("err").is_some() {
            if r.get("err") != py.get("err") {
                format!("refusal {:?} vs python {:?}", r.get("err"), py.get("err"))
            } else {
                String::new()
            }
        } else if r["kind"] != py["kind"] {
            format!("kind {} vs python {}", r["kind"], py["kind"])
        } else {
            let mut w = String::new();
            for k in ["nf", "ne", "nv", "fN", "fD", "EcQ", "FC", "Vname", "E", "F", "dom"] {
                if py.get(k).is_some() && norm(&r[k]) != norm(&py[k]) {
                    w = format!("{} differs", k);
                    break;
                }
            }
            if w.is_empty() {
                if let (Some(a), Some(b)) = (r["evalV"].as_array(), py["evalV"].as_array()) {
                    for (t, (x, y)) in a.iter().zip(b.iter()).enumerate() {
                        let (x, y) = (as_f(x), as_f(y));
                        let same = if x.is_finite() && y.is_finite() { (x - y).abs() <= 1e-9 * 1f64.max(y.abs()) } else { x == y || (x.is_nan() && y.is_nan()) };
                        if !same {
                            w = format!("eval[{}] {} vs python {}", t, x, y);
                            break;
                        }
                    }
                }
            }
            w
        };
        if !why.is_empty() {
            bad.push(format!("{} #{}: {}", trace, i, why));
        }
    }
    assert!(bad.is_empty(), "{} record(s) differ from the Python port:\n{}", bad.len(), bad.join("\n"));
}

#[test]
fn conj_trace_matches_python() {
    check("conj.jsonl", "conj.py.jsonl");
}

#[test]
fn biconj_trace_matches_python() {
    check("biconj.jsonl", "biconj.py.jsonl");
}

#[test]
#[ignore]
fn slow_suites_match_python() {
    for s in ["slow_algCrvSubdivisionTest", "slow_envelopeSubdivisionSlowTest", "slow_minSubdivisionTest", "slow_kinkSubdivisionTest", "slow_erratumBiconjEx2Test"] {
        check(&format!("{}.jsonl", s), &format!("{}.py.jsonl", s));
    }
}

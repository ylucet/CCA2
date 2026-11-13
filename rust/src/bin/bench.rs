//! Time conjQ / biconjQ on the benchmark fixtures (CCA2/python/bench/fixtures.json), with the
//! protocol of bench_python.py / bench_matlab.m: the first call, then `reps` more unless the first
//! took longer than slowCap seconds; wall clock inside one process. The memos that MATLAB does NOT
//! have are emptied before every timed call (the same set bench_python.py empties). One JSON line
//! per fixture, in bench_python.py's format plus the result's kind and face count.
//!
//!     cargo run --release --bin bench -- fixtures.json out.jsonl [reps] [slowCap]
use cca2::alg_con::AlgCon;
use cca2::io::num;
use cca2::qua_pol::QuaPol;
use serde_json::{json, Value};
use std::io::Write;
use std::time::Instant;

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

fn cold_rust_only_memos() {
    cca2::conic_meet::clear_memo();
    cca2::alg_row_fails::clear_memos();
    cca2::rat_q::clear_memos();
}

enum Out {
    Con(usize),
    Alg(&'static str, Option<usize>),
}

fn call(fname: &str, p: &QuaPol) -> Result<Out, String> {
    if fname == "conjQ" {
        cca2::conj_q::conj_q(p).map(|(g, _)| Out::Con(g.nf)).map_err(|e| e.id)
    } else {
        cca2::biconj_q::biconj_q(p, &Default::default()).map(|h: AlgCon| Out::Alg(h.kind(), h.nf())).map_err(|e| e.id)
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let fix: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let mut out = std::fs::File::create(&a[2]).unwrap();
    let reps: usize = a.get(3).map(|s| s.parse().unwrap()).unwrap_or(3);
    let slow_cap: f64 = a.get(4).map(|s| s.parse().unwrap()).unwrap_or(20.0);
    let only: Option<String> = a.get(5).cloned();
    let w = QuaPol::whole(vec![1, 0, 1, 0, 0, 0].into_iter().map(|v| cca2::rat_q::Num::Q(cca2::rat_q::qi(v))).collect()).unwrap();
    let _ = call("conjQ", &w);
    let _ = call("biconjQ", &w);
    for x in &fix {
        if let Some(o) = &only {
            if x["bucket"].as_str() != Some(o.as_str()) {
                continue;
            }
        }
        let f = build(x);
        let fname = x["fn"].as_str().unwrap();
        let mut t: Vec<f64> = vec![];
        let mut err = String::new();
        let mut kind = String::new();
        let mut nf: Option<usize> = None;
        cold_rust_only_memos();
        let t0 = Instant::now();
        match call(fname, &f) {
            Err(e) => err = e,
            Ok(o) => {
                t.push(t0.elapsed().as_secs_f64());
                match o {
                    Out::Con(n) => {
                        kind = "QuaCon".into();
                        nf = Some(n);
                    }
                    Out::Alg(k, n) => {
                        kind = k.into();
                        nf = n;
                    }
                }
                if t[0] <= slow_cap {
                    for _ in 0..reps {
                        cold_rust_only_memos();
                        let t0 = Instant::now();
                        let _ = call(fname, &f);
                        t.push(t0.elapsed().as_secs_f64());
                    }
                }
            }
        }
        let best = t.iter().cloned().fold(f64::INFINITY, f64::min);
        let worst = t.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let rec = json!({"id": x["id"], "first": t.first().cloned().unwrap_or(-1.0), "best": if t.is_empty() { -1.0 } else { best },
                         "worst": if t.is_empty() { -1.0 } else { worst }, "n": t.len(), "err": err, "kind": kind, "nf": nf});
        writeln!(out, "{}", rec).unwrap();
        out.flush().unwrap();
        println!("{:4} {:<8} {:<11} first {:9.4}  best {:9.4}  {} {} {:?}", x["id"], fname, x["kind"].as_str().unwrap_or(""), t.first().cloned().unwrap_or(-1.0), best, err, kind, nf);
    }
}

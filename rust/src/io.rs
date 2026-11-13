//! io -- read the golden traces (JSON lines) and write results in the same shape (port of the
//! parts of io_matlab.py the parity tests use).

use crate::errors::R;
use crate::qua_con::QuaCon;
use crate::qua_pol::QuaPol;
use crate::rat_q::*;
use serde_json::{json, Value};

/// A JSON number or string as the caller's number: an int stays exact, a float is a double.
pub fn num(v: &Value) -> Num {
    match v {
        Value::Number(n) => {
            let s = n.to_string();
            if s.contains('.') || s.contains('e') || s.contains('E') {
                Num::F(s.parse::<f64>().unwrap())
            } else {
                Num::Q(Q::from_integer(s.parse::<Z>().unwrap()))
            }
        }
        Value::String(s) => Num::F(match s.as_str() {
            "inf" => f64::INFINITY,
            "-inf" => f64::NEG_INFINITY,
            _ => f64::NAN,
        }),
        Value::Bool(b) => Num::Q(qi(*b as i64)),
        _ => panic!("not a number: {}", v),
    }
}

pub fn numf(v: &Value) -> f64 {
    match num(v) {
        Num::F(f) => f,
        Num::Q(q) => q_to_f64(&q),
    }
}

pub fn zint(v: &Value) -> Z {
    match num(v) {
        Num::Q(q) => q.to_integer(),
        Num::F(f) => Z::from(f as i128),
    }
}

/// MATLAB _mat: a vector becomes one row.
pub fn mat(a: &Value) -> Vec<Vec<Value>> {
    let arr = match a.as_array() {
        None => return vec![],
        Some(x) => x,
    };
    if arr.is_empty() {
        return vec![];
    }
    if !arr[0].is_array() {
        return vec![arr.clone()];
    }
    arr.iter().map(|r| r.as_array().unwrap().clone()).collect()
}

fn is_empty(v: &Value) -> bool {
    match v {
        Value::Array(a) => a.is_empty(),
        Value::Null => true,
        _ => false,
    }
}

/// Build the QuaPol from a record's input exactly as the MATLAB caller did (quaPolFrom).
pub fn qua_pol_from(j: &Value) -> R<QuaPol> {
    let f: Vec<Vec<Num>> = mat(&j["f"]).iter().map(|r| r.iter().map(num).collect()).collect();
    let nv = j["nv"].as_u64().unwrap_or(0);
    let mut p = if nv == 0 && is_empty(&j["V"]) {
        // QuaPol(f[0] if len(f) == 1 else f): a list of rows has len != 1 -> treated as rows
        QuaPol::whole_rows(f)?
    } else {
        let v: Vec<Vec<Num>> = mat(&j["V"]).iter().map(|r| r.iter().map(num).collect()).collect();
        let e: Vec<[usize; 3]> = mat(&j["E"]).iter().map(|r| [idx(&r[0]), idx(&r[1]), idx(&r[2])]).collect();
        let ff: Vec<[usize; 2]> = mat(&j["F"]).iter().map(|r| [idx(&r[0]), idx(&r[1])]).collect();
        QuaPol::new(v, e, f, ff)?
    };
    if !j["exact"].as_bool().unwrap_or(true) {
        p.exact = false;
    }
    if let Some(b) = j.get("fIsConvex") {
        p.f_is_convex = Some(b.as_bool().unwrap_or_else(|| numf(b) != 0.0));
    }
    Ok(p)
}

fn idx(v: &Value) -> usize {
    numf(v) as usize
}

fn zs(v: &[Z]) -> Value {
    Value::Array(v.iter().map(|x| serde_json::from_str::<Value>(&x.to_string()).unwrap()).collect())
}

pub fn qua_con_json(g: &QuaCon) -> Value {
    json!({
        "kind": "QuaCon",
        "nv": g.nv, "ne": g.ne, "nf": g.nf,
        "Vname": g.vname, "E": g.e,
        "EcQ": g.ecq.iter().map(|c| zs(c)).collect::<Vec<_>>(),
        "fN": g.fn_.iter().map(|c| zs(c)).collect::<Vec<_>>(),
        "fD": zs(&g.fd),
        "F": g.f, "FC": g.fc,
    })
}

pub fn dom_json(d: &[Vec<Z>]) -> Value {
    Value::Array(d.iter().map(|r| zs(r)).collect())
}

pub fn eval_x(out: &Value) -> Vec<[f64; 2]> {
    match out.get("evalX") {
        None => vec![],
        Some(x) => mat(x).iter().map(|r| [numf(&r[0]), numf(&r[1])]).collect(),
    }
}

/// A double as JSON (inf / nan as strings, as MATLAB wrote them).
pub fn fj(x: f64) -> Value {
    if x.is_nan() {
        json!("nan")
    } else if x.is_infinite() {
        json!(if x > 0.0 { "inf" } else { "-inf" })
    } else {
        json!(x)
    }
}

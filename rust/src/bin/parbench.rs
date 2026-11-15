//! parbench -- how conjQ / biconjQ scale when SEVERAL CALLS run at once, one call per thread.
//!
//! The engine itself is single-threaded, but a call owns everything it touches (its data are built
//! inside the calling thread and every memo is `thread_local!`), so independent calls can run on
//! different threads without any change to the engine. This driver measures that:
//!
//!     parbench fixtures.json out.jsonl THREADS            # every fixture once, THREADS workers
//!     parbench fixtures.json out.jsonl THREADS copies ID  # THREADS copies of fixture ID at once
//!     parbench fixtures.json out.jsonl 0 [copies ID]      # on the main thread, no thread created
//!
//! `batch` (the default) hands the fixtures to the workers from a shared queue, heaviest first when
//! a previous out.jsonl-style file is given in CCA2_PARBENCH_ORDER (else in file order), and writes
//! one JSON line per fixture (id, seconds, kind, nf, err) so a THREADS = 1 run can be compared with
//! a parallel one. `copies` is the calibration run for PARALLEL.md: N identical calls on N threads
//! take as long as one call only if the machine (cores, memory allocator, caches) scales.
//! THREADS = 0 matters because glibc's allocator is faster while the process has only ever had one
//! thread: the heaviest fixture takes about 20% longer once a second thread exists.
//!
//! As in bench.rs, the memos MATLAB does not have are emptied before each call.
use cca2::io::num;
use cca2::qua_pol::QuaPol;
use serde_json::{json, Value};
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
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

/// One call, timed: (seconds, kind, nf, err).
fn run(x: &Value) -> (f64, String, Option<usize>, String) {
    let f = build(x);
    cca2::conic_meet::clear_memo();
    cca2::alg_row_fails::clear_memos();
    cca2::rat_q::clear_memos();
    let t0 = Instant::now();
    let r = if x["fn"].as_str() == Some("conjQ") {
        cca2::conj_q::conj_q(&f).map(|(g, _)| ("QuaCon".to_string(), Some(g.nf))).map_err(|e| e.id)
    } else {
        cca2::biconj_q::biconj_q(&f, &Default::default()).map(|h| (h.kind().to_string(), h.nf())).map_err(|e| e.id)
    };
    let t = t0.elapsed().as_secs_f64();
    match r {
        Ok((k, n)) => (t, k, n, String::new()),
        Err(e) => (t, String::new(), None, e),
    }
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let fix: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let threads: usize = a[3].parse().unwrap();
    let copies = a.get(4).map(|s| s == "copies").unwrap_or(false);
    let jobs: Vec<&Value> = if copies {
        let id: i64 = a[5].parse().unwrap();
        let x = fix.iter().find(|x| x["id"].as_i64() == Some(id)).expect("no such fixture id");
        vec![x; threads.max(1)]
    } else {
        let mut j: Vec<&Value> = fix.iter().collect();
        if let Ok(p) = std::env::var("CCA2_PARBENCH_ORDER") {
            let mut t = std::collections::HashMap::new();
            for l in std::fs::read_to_string(p).unwrap().lines() {
                let v: Value = serde_json::from_str(l).unwrap();
                t.insert(v["id"].as_i64().unwrap(), v["t"].as_f64().or(v["best"].as_f64()).unwrap_or(0.0));
            }
            j.sort_by(|x, y| t.get(&y["id"].as_i64().unwrap()).partial_cmp(&t.get(&x["id"].as_i64().unwrap())).unwrap());
        }
        j
    };
    let next = AtomicUsize::new(0);
    let out: Mutex<Vec<Value>> = Mutex::new(vec![]);
    if threads == 0 {
        // no thread is ever created: the single-threaded reference (see PARALLEL.md, "thread tax")
        let t0 = Instant::now();
        let recs: Vec<(f64, String, Option<usize>, String)> = jobs.iter().map(|x| run(x)).collect();
        let sum: f64 = recs.iter().map(|r| r.0).sum();
        println!("main thread only, calls {} wall {:.3} s  sum of calls {:.3} s", recs.len(), t0.elapsed().as_secs_f64(), sum);
        return;
    }
    let t0 = Instant::now();
    std::thread::scope(|s| {
        for _ in 0..threads {
            // the main thread has an 8 MB stack; give the workers the same room for the recursions
            std::thread::Builder::new()
                .stack_size(64 << 20)
                .spawn_scoped(s, || loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= jobs.len() {
                        break;
                    }
                    let (t, kind, nf, err) = run(jobs[i]);
                    out.lock().unwrap().push(json!({"id": jobs[i]["id"], "t": t, "kind": kind, "nf": nf, "err": err}));
                })
                .unwrap();
        }
    });
    let wall = t0.elapsed().as_secs_f64();
    let mut recs = out.into_inner().unwrap();
    recs.sort_by_key(|r| r["id"].as_i64().unwrap());
    let sum: f64 = recs.iter().map(|r| r["t"].as_f64().unwrap()).sum();
    let mx = recs.iter().map(|r| r["t"].as_f64().unwrap()).fold(0.0, f64::max);
    let mut f = std::fs::File::create(&a[2]).unwrap();
    for r in &recs {
        writeln!(f, "{}", r).unwrap();
    }
    println!("threads {} calls {} wall {:.3} s  sum of calls {:.3} s  longest call {:.3} s  sum/wall {:.2}", threads, recs.len(), wall, sum, mx, sum / wall);
}

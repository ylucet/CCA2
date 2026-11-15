//! prof -- a sampling profile of ONE conjQ / biconjQ call, for machines where `perf` is not
//! permitted (kernel.perf_event_paranoid = 4 on Frances). SIGPROF every 2 ms of CPU time, glibc
//! backtrace() in the handler, raw return addresses (relative to the load base) written to a file;
//! tools/prof_report.py resolves them with addr2line and prints inclusive / self time.
//!
//!     cargo run --release --bin prof -- ../python/bench/fixtures.json <id> out.stacks
//!
//! A development tool: the engine does not depend on it. Linux + glibc only.
use cca2::io::num;
use cca2::qua_pol::QuaPol;
use serde_json::Value;
use std::io::Write;
use std::os::raw::{c_int, c_long, c_void};
use std::sync::atomic::{AtomicUsize, Ordering};

const DEPTH: usize = 96;
const CAP: usize = 400_000;
static mut BUF: *mut [usize; DEPTH + 1] = std::ptr::null_mut();
static N: AtomicUsize = AtomicUsize::new(0);

#[repr(C)]
struct Timeval {
    sec: c_long,
    usec: c_long,
}
#[repr(C)]
struct Itimerval {
    interval: Timeval,
    value: Timeval,
}
extern "C" {
    fn backtrace(buf: *mut *mut c_void, size: c_int) -> c_int;
    fn signal(sig: c_int, handler: extern "C" fn(c_int)) -> usize;
    fn setitimer(which: c_int, new: *const Itimerval, old: *mut Itimerval) -> c_int;
}

extern "C" fn on_prof(_: c_int) {
    let i = N.fetch_add(1, Ordering::Relaxed);
    if i >= CAP {
        return;
    }
    unsafe {
        let row = &mut *BUF.add(i);
        let n = backtrace(row.as_mut_ptr().add(1) as *mut *mut c_void, DEPTH as c_int);
        row[0] = n as usize;
    }
}

fn timer(usec: c_long) {
    let t = Itimerval { interval: Timeval { sec: 0, usec }, value: Timeval { sec: 0, usec } };
    unsafe {
        setitimer(2, &t, std::ptr::null_mut()); // ITIMER_PROF
    }
}

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
    let fix: Vec<Value> = serde_json::from_str(&std::fs::read_to_string(&a[1]).unwrap()).unwrap();
    let id: i64 = a[2].parse().unwrap();
    let x = fix.iter().find(|x| x["id"].as_i64() == Some(id)).expect("no such fixture id");
    let f = build(x);
    let fname = x["fn"].as_str().unwrap();
    // the load base of the executable, so the addresses can be given to addr2line
    let maps = std::fs::read_to_string("/proc/self/maps").unwrap();
    let exe = std::env::current_exe().unwrap();
    let base = maps
        .lines()
        .find(|l| l.ends_with(exe.to_str().unwrap()))
        .map(|l| usize::from_str_radix(l.split('-').next().unwrap(), 16).unwrap())
        .unwrap();
    unsafe {
        BUF = Box::leak(vec![[0usize; DEPTH + 1]; CAP].into_boxed_slice()).as_mut_ptr();
        let mut warm = [std::ptr::null_mut::<c_void>(); 4];
        backtrace(warm.as_mut_ptr(), 4); // loads the unwinder outside the handler
        signal(27, on_prof); // SIGPROF
    }
    let t0 = std::time::Instant::now();
    timer(2000);
    let ok = if fname == "conjQ" { cca2::conj_q::conj_q(&f).is_ok() } else { cca2::biconj_q::biconj_q(&f, &Default::default()).is_ok() };
    timer(0);
    let wall = t0.elapsed().as_secs_f64();
    let n = N.load(Ordering::Relaxed).min(CAP);
    let mut out = std::io::BufWriter::new(std::fs::File::create(&a[3]).unwrap());
    writeln!(out, "# exe {} id {} fn {} ok {} wall {:.3} samples {}", exe.display(), id, fname, ok, wall, n).unwrap();
    for i in 0..n {
        let row = unsafe { &*BUF.add(i) };
        let k = row[0];
        // frames 0..1 are the handler and the signal trampoline
        let s: Vec<String> = row[1..=k].iter().skip(2).map(|p| format!("{:x}", p.wrapping_sub(base))).collect();
        writeln!(out, "{}", s.join(" ")).unwrap();
    }
    eprintln!("id {} {} ok={} wall {:.3} s, {} samples", id, fname, ok, wall, n);
}

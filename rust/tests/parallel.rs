//! The parallel engine returns EXACTLY what the serial engine returned: `digest` hashes the complete
//! result (every face, row and coefficient, in order) of every one of the 200 benchmark fixtures, and
//! the hashes recorded from the serial engine (bench/parallel/digest_serial.tsv, taken before any
//! parallelism was added) must be reproduced with 1 thread (CCA2_PAR=0), 2, 8 and 32 threads.
//!
//!     cargo test --release --test parallel     # ~30 s
use std::path::PathBuf;
use std::process::Command;

fn run(env: &[(&str, &str)]) -> String {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut c = Command::new(env!("CARGO_BIN_EXE_digest"));
    c.arg(root.join("../python/bench/fixtures.json"));
    for (k, v) in env {
        c.env(k, v);
    }
    let out = c.output().expect("digest did not start");
    assert!(out.status.success(), "digest failed: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

fn check(env: &[(&str, &str)]) {
    let want = std::fs::read_to_string(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("bench/parallel/digest_serial.tsv")).unwrap();
    let got = run(env);
    let (w, g): (Vec<&str>, Vec<&str>) = (want.lines().collect(), got.lines().collect());
    assert_eq!(w.len(), g.len(), "different number of fixtures");
    let bad: Vec<String> = w.iter().zip(g.iter()).filter(|(a, b)| a != b).map(|(a, b)| format!("serial {} / now {}", a, b)).collect();
    assert!(bad.is_empty(), "{} fixtures differ with {:?}:\n{}", bad.len(), env, bad.join("\n"));
}

#[test]
fn one_thread_is_the_serial_engine() {
    check(&[("CCA2_PAR", "0")]);
}

#[test]
fn two_threads() {
    check(&[("CCA2_THREADS", "2")]);
}

#[test]
fn eight_threads() {
    check(&[("CCA2_THREADS", "8")]);
}

#[test]
fn thirty_two_threads() {
    check(&[("CCA2_THREADS", "32")]);
}

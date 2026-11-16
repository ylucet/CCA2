//! Optional accelerator: link the FLINT library that ships with python-flint (the same one the
//! Python port's benchmark used) when it can be found. Without it the pure-Rust Yun
//! decomposition is used; the answers are the same (python/cca2/_poly.py has the rationale) but the
//! heavy fixtures (erratum Example 2, kink) take minutes instead of seconds.
//!
//! Search order for a directory holding libflint*.so*:
//!   1. CCA2_FLINT_LIB_DIR
//!   2. ../python/.pydeps/python_flint.libs                    (this checkout)
//!   3. <main checkout of this repository>/python/.pydeps/python_flint.libs
//!      (a git worktree or clone has no .pydeps of its own: it is untracked)
//!   4. the system library directories (a distribution's libflint)
//!
//! CCA2_FLINT_BUILD=0        build without FLINT, on purpose (no warning)
//! CCA2_FLINT_BUILD=require  fail the build when FLINT is not found (benchmarks, suite runs)
//! Anything else: link FLINT if found, else WARN LOUDLY and fall back.
use std::path::{Path, PathBuf};
use std::process::Command;

fn libflint_in(dir: &Path) -> Option<String> {
    std::fs::read_dir(dir).ok().and_then(|rd| {
        rd.filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().to_string())
            .find(|n| n.starts_with("libflint") && n.contains(".so"))
    })
}

/// The working tree of the main checkout, from the shared git directory (differs from this one in
/// a worktree). Best effort: any failure is just "not found".
fn main_checkout(manifest: &Path) -> Option<PathBuf> {
    let out = Command::new("git").args(["rev-parse", "--path-format=absolute", "--git-common-dir"]).current_dir(manifest).output().ok()?;
    if !out.status.success() {
        return None;
    }
    let gd = PathBuf::from(String::from_utf8_lossy(&out.stdout).trim());
    gd.parent().map(|p| p.to_path_buf())
}

fn main() {
    println!("cargo:rustc-check-cfg=cfg(have_flint)");
    println!("cargo:rerun-if-env-changed=CCA2_FLINT_LIB_DIR");
    println!("cargo:rerun-if-env-changed=CCA2_FLINT_BUILD");
    println!("cargo:rerun-if-changed=build.rs");
    let mode = std::env::var("CCA2_FLINT_BUILD").unwrap_or_default();
    if mode == "0" {
        return;
    }
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
    let mut cands: Vec<PathBuf> = vec![];
    if let Ok(d) = std::env::var("CCA2_FLINT_LIB_DIR") {
        cands.push(d.into());
    }
    cands.push(manifest.join("../python/.pydeps/python_flint.libs"));
    if let Some(m) = main_checkout(&manifest) {
        cands.push(m.join("python/.pydeps/python_flint.libs"));
    }
    for d in ["/usr/lib/x86_64-linux-gnu", "/usr/lib64", "/usr/lib", "/usr/local/lib"] {
        cands.push(d.into());
    }
    for dir in &cands {
        if let Some(lib) = libflint_in(dir) {
            let dir = dir.canonicalize().unwrap();
            println!("cargo:warning=FLINT linked from {}", dir.display());
            println!("cargo:rustc-link-search=native={}", dir.display());
            println!("cargo:rustc-link-lib=dylib:+verbatim={}", lib);
            println!("cargo:rustc-link-arg=-Wl,-rpath,{}", dir.display());
            println!("cargo:rustc-cfg=have_flint");
            return;
        }
    }
    let tried: Vec<String> = cands.iter().map(|p| p.display().to_string()).collect();
    if mode == "require" {
        panic!("CCA2_FLINT_BUILD=require but no libflint was found; tried: {}. Set CCA2_FLINT_LIB_DIR.", tried.join(", "));
    }
    println!("cargo:warning=FLINT NOT FOUND: building the slow pure-Rust square-free path. Fixtures 191 and 196 take minutes instead of seconds, so do not benchmark this build. Set CCA2_FLINT_LIB_DIR (tried: {}).", tried.join(", "));
}

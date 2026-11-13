//! Optional accelerator: link the FLINT library that ships with python-flint (the same one the
//! Python port's benchmark used) when it can be found. Without it the pure-Rust Yun
//! decomposition is used; the answers are the same (python/cca2/_poly.py has the rationale).
//!
//! Set CCA2_FLINT_LIB_DIR to the directory holding libflint*.so*, or CCA2_FLINT_BUILD=0 to build
//! without it.
use std::path::PathBuf;

fn main() {
    println!("cargo:rustc-check-cfg=cfg(have_flint)");
    println!("cargo:rerun-if-env-changed=CCA2_FLINT_LIB_DIR");
    println!("cargo:rerun-if-env-changed=CCA2_FLINT_BUILD");
    if std::env::var("CCA2_FLINT_BUILD").map(|v| v == "0").unwrap_or(false) {
        return;
    }
    let dir: PathBuf = match std::env::var("CCA2_FLINT_LIB_DIR") {
        Ok(d) => d.into(),
        Err(_) => {
            let m = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
            m.join("../python/.pydeps/python_flint.libs")
        }
    };
    let lib = std::fs::read_dir(&dir).ok().and_then(|rd| {
        rd.filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().to_string())
            .find(|n| n.starts_with("libflint") && n.contains(".so"))
    });
    if let Some(lib) = lib {
        let dir = dir.canonicalize().unwrap();
        println!("cargo:rustc-link-search=native={}", dir.display());
        println!("cargo:rustc-link-lib=dylib:+verbatim={}", lib);
        println!("cargo:rustc-link-arg=-Wl,-rpath,{}", dir.display());
        println!("cargo:rustc-cfg=have_flint");
    }
}

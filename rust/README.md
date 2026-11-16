# rust/ -- Rust port of the CCA2 exact conjugate / biconjugate engine

Status (2026-11-13): **complete for conjQ and biconjQ, at parity with the Python port.**

The crate is a module-for-module port of `CCA2/python/cca2` (itself the port of the MATLAB
reference). Module names follow the Python ones; index values stay 1-based in the data, 0 = none;
exact data are `num-bigint` / `num-rational`; refusals carry the SAME identifier strings
(`Err(CErr { id, msg })` in place of `CCA2Error`).

## Verification

`cargo test --release` (add `-- --include-ignored` for the slow suites, ~20 s in all) replays every
`conjQ` / `biconjQ` call of the MATLAB golden traces (`python/golden`) and compares with what the
Python port returns for the same call (`ref/*.py.jsonl`, written by `tools/py_ref.py`):

| trace | records | match Python |
|---|---|---|
| conj.jsonl | 176 | 176 -- every QuaCon field (fN, fD, EcQ, FC, Vname, E, F), dom, sampled values |
| biconj.jsonl | 114 | 114 -- refusals by identifier, kind, face count, sampled values |
| slow suites (kink, Example 2, min-subdivision, ...) | 20 | 20 -- kink AlgCrv 317 faces, Example 2 EnvelopeMin |

Both with and without FLINT (`CCA2_FLINT=0`).

## Running

    cargo build --release
    ./target/release/golden <trace.jsonl>            # one JSON result per record
    python3 tools/compare.py <trace> <rust.jsonl> <python.jsonl>
    ./target/release/bench ../python/bench/fixtures.json out.jsonl
    python3 tools/compare3.py                        # MATLAB / Python / Rust table

FLINT: `build.rs` links the libflint shipped with python-flint (`../python/.pydeps`), the same
accelerator the Python port uses, for exactly one kernel (the odd-multiplicity part of a bivariate
polynomial, `src/flint.rs`). Without it (`CCA2_FLINT_BUILD=0` at build time, or `CCA2_FLINT=0` at
run time) the pure-Rust Yun decomposition is used; the answers are identical, but erratum Example 2
and the kink then take minutes (as the pure-Python path does).

Where `build.rs` looks, in order: `CCA2_FLINT_LIB_DIR`; `../python/.pydeps/python_flint.libs`; the
same folder in the MAIN checkout of the repository (`.pydeps` is untracked, so a git worktree or a
fresh clone has none of its own); the system library directories. Build-time switches:
`CCA2_FLINT_BUILD=0` builds without FLINT on purpose; `CCA2_FLINT_BUILD=require` fails the build if
FLINT is not found (use it for benchmark and suite builds); otherwise a missing FLINT gives a
`cargo:warning=FLINT NOT FOUND` and the slow build. `bench`, `parbench` and `prof` print one line,
`FLINT: linked, in use` or `FLINT: NOT linked`, on stderr at start; do not compare timings from a
build that says NOT linked. Check an existing binary with `ldd target/release/bench | grep flint`.

## Speed (bench/RESULTS.md)

Same 200 fixtures and protocol as `python/bench` (MATLAB and Python not re-measured: run5):

| | all 200 fixtures (best times) | vs MATLAB |
|---|---|---|
| MATLAB | 1126 s | 1x |
| Python (numba + FLINT) | 88.5 s | 12.7x |
| Rust (single thread, FLINT) | 15.7 s | 72x (5.7x over Python) |

Shared machine; the runs were on different days under different load -- the ratios are indicative.

## Parallelism (PARALLEL.md)

The engine runs on a rayon pool: independent pieces / cells / pairs / corners run in parallel and the
result is bit-identical to the serial engine (`tests/parallel.rs`). `CCA2_THREADS=n` sets the pool
size (default all hardware threads), `CCA2_THREADS=1` or `CCA2_PAR=0` is the serial reference. On
Frances (32 cores / 64 threads, shared) the 200-fixture benchmark goes from 14.8 s to 3.0 s (4.9x at
16 threads); details and the reduction study in PARALLEL.md section 11.

## Not ported (Python-only)

- `check.py` (the definition oracle `checkBiconj`) and the intermediate-function trace tests
  (`trace_*.jsonl`); parity is checked end to end instead.
- Whatever the Python port itself declines (see `CCA2/python/PORT_NOTES.md` section 5).

The earlier, partial Rust attempt this replaces stays outside the repository (`AI/cca2-rs/legacy`).

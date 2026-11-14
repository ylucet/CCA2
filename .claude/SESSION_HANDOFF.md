# Session Handoff

_Last updated: 2026-11-13T(late evening)-08:00_

## What happened this session
The Rust port of the exact conj/biconj engine (`rust/`, moved there from `~/AI/cca2-rs`) was finished
and brought to parity with the Python port: conj 176, biconj 114 and the 20 slow golden records all
match Python, with and without FLINT. It links python-flint's bundled libflint via `build.rs` for
square-free parts, which is what Python does too; `CCA2_FLINT=0` switches it off at runtime. Two
`/code-review` passes (low on ca912e6, high on five risky modules) were run and all findings fixed.
The MATLAB sources moved from the root to `MATLAB/` and every caller was updated (suite driver, stage
cache, golden export scripts, coverage, bench). The MATLAB fast suite gives 834/0 after the move.

Benchmark (`rust/bench/RESULTS.md`): all 200 fixtures give MATLAB 1126 s, Python 88.5 s (12.7x) and
Rust 15.7 s (5.7x over Python, 72x over MATLAB). The MATLAB and Python numbers are from python/bench
run5 and were not re-measured. The machine was SHARED, so read these as ratios, not absolute claims.

## Where things stand
- Branch: `main` @ `27e9464` ("chore: update session handoff"); a later bookkeeping commit added CURRENT_STATE.md
- Pushed: everything up to 27e9464 is on origin/main; the CURRENT_STATE.md commit is pending confirmation.
- Tests:
  - Rust: `cd rust && cargo test --release` covers parity. The slow suites are `#[ignore]`; run them with `-- --ignored`.
  - MATLAB: fast 834/0 (2026-11-13; ~574 s because the move invalidated the stage cache, ~90 s once it is warm).
  - Python: unchanged.
- Untracked and deliberately left alone:
  - `RUST_PORT.md` (not created by Claude) and `test.txt`.
  - The `python/dbg_*.py`, `cmp_order.py` and `classify_biconj.py` scratch scripts.
- The shared-config repo (`~/AI/.shared.git`) got a mybiblio remote, which was then removed again; it is
  back to its original state. Nothing was pushed to mybiblio, because the Frances token has no access to it.

## Next steps (not yet chosen)
1. Delete `~/AI/cca2-rs`: it only holds legacy/ and target/, and everything else is in `rust/`. Ask first.
2. If biblio/ should go to github.com/ylucet/mybiblio:
   - Grant the token access first.
   - Decide what goes there. biblio/ currently lives inside the shared-config work tree.
3. Older items:
   - The bridgeFaceAlgTest tieMax/tiePoints sweep window.
   - Release notes on the dependencies (MATLAB + SMT, Python and Rust ports).
   - Example 2 with overlayBudget=Inf, now cheap enough to time in Rust or Python.

## Relevant files
- `rust/README.md`: how to build and test, and the FLINT switch.
- `rust/bench/RESULTS.md`: the speedup table.
- `rust/tools/compare3.py`: the three-way bench table.
- `rust/tests/parity.rs`: replays the golden traces against the Python outputs in `rust/ref/`.
- `rust/tools/py_ref.py`: regenerates those Python outputs.
- `rust/src/flint.rs`, `rust/build.rs`: the FLINT FFI.
- `MATLAB/`: the MATLAB sources.
- `.claude/suite.sh`: runs MATLAB from `MATLAB/`.
- `MATLAB/plqStage.m`: its cache is in `.claude/stagecache`.

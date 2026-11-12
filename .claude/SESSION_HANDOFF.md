# Session Handoff

_Last updated: 2026-11-12T(evening)-07:00_

## What happened this session
The Python port of the exact conj/biconj engine (`python/`, overnight/2026-11-10) was merged to
main, the kink face-count mismatch was fixed (334 -> 317 faces, now identical to MATLAB), and the
Python engine was sped up in four rounds, each bit-identical to the code before it (differential
tests against the previous module plus the golden tests):
1. Optional accelerators: python-flint (square-free parts, univariate gcd) and numba (Bernstein test).
2. conjQ path: batched margin LP, conicMeet pair memo, vectorised interior probe.
3. AlgCrv path: Sturm chains built once, vectorised cell sampling, memoised factor screens.
4. EnvelopeMin path (erratum Example 2): integer exact fibres, integer dyadic bisection in signAt,
   vectorised surd-corner search.

Benchmark (python/bench/RESULTS.md, last section; best of 2 passes, MATLAB `-singleCompThread` and
Python run concurrently on a SHARED machine -- ratios, not absolute claims):
- all 200 fixtures: Python 89 s vs MATLAB 1126 s (Python's own path: 939 -> 495 -> 372 -> 113 -> 89 s)
- Example 2 (EnvelopeMin) 27 s vs 566 s; kink AlgCrv 42 s vs 413 s; min-subdivision 10.5 s vs 73 s
- conjQ: 2.7 s vs 10.0 s

## Where things stand
- Branch: `main` @ handoff commit after `8fd2b36` ("bench: rerun after the EnvelopeMin speed-up")
- Pushed: yes, at the end of this session (on request)
- Tests: Python fast 2733 passed / 1 skipped (~20 s, with and without FLINT/numba), Python slow
  23/0 (~48 s). MATLAB last run: fast 833/0, normal 35/0, 8 exact slow suites green (2026-11-11).
- Known reds: none. Blocked: nothing.
- Running Python: `cd python && PYTHONPATH=.pydeps python3 -m pytest -q -n 32` (python-flint 0.9.0
  lives in python/.pydeps; numba 0.62.1 in ~/.local). Both are optional; results are identical without.
- Untracked and deliberately left alone: `RUST_PORT.md` (not created by Claude), `test.txt`.

## Next steps (not yet chosen)
1. Python profile on Example 2 is now flat (largest item: numpy slack evaluation in cell sampling,
   ~12%); further speed-ups would be small. Stop here unless a new workload appears.
2. Test oracle: `tieMax`/`tiePoints` in bridgeFaceAlgTest.m should scale their sweep window to the
   tie centre (7 false disagreements at large data).
3. Release/paper: state the MATLAB + Symbolic Math Toolbox dependency (and now the Python port);
   cite `convexdb` by its Zenodo DOI.
4. Example 2 with `overlayBudget=Inf` (full subdivision) -- now worth re-timing in Python.

## Relevant files
- `python/cca2/rat_q.py` -- Sturm chain memo, dyadic signAt, square-free memo; `_poly.py` FLINT facade
- `python/cca2/alg_row_fails.py`, `alg_branch.py`, `cell_sample_point.py`, `conic_meet.py`,
  `cell_interior.py`, `cell_has_interior.py`, `poly_sign_on.py` -- the speed-ups
- `python/bench/` -- bench_matlab.m, bench_python.py (clears Python-only memos per call), compare.py,
  RESULTS.md, run2..run5 raw data
- `TODO.md` 'Next up -- 2026-11-10'; `MORNING.md`; `DECISIONS.md` 2026-11-11

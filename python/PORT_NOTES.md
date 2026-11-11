# Port notes — where the Python port deliberately differs from MATLAB, and what it found

_Started 2026-11-10 (overnight run). The plan is `PORT_PLAN.md`; the measured scope is `PORT_SCOPE.md`._

## 1. Conventions chosen during the port

- **Index values stay 1-based in the data** (`E`, `F`, `FC`, `Vname`, face/edge numbers, branch
  1/2), with 0 = "none", exactly as in MATLAB. Only Python list subscripts shift by one.
  The plan said "0-based internally"; this was changed on purpose: fixtures in the MATLAB tests,
  the papers and the goldens are then VERBATIM, the sentinels (`F == 0` outside, `E(j,2) == 0`
  for a ray) keep their meaning, and golden comparison is field-for-field with no translation
  layer to get wrong. A 0-based public API can be put on top later if wanted.
- **Function names keep the MATLAB spelling** inside each module (`rat_q.feasible2`), module
  names are snake_case (`conjQ.m` -> `cca2/conj_q.py`). The two code bases grep against each other.
- **Exact data are plain Python `int` / `Fraction`** in lists. No numpy integer ever touches
  exact data (numpy's int64 wraps silently). numpy is used only for double-precision linear
  algebra and polynomial roots in the documented float screens.
- **Comments.** Each module keeps the objective and the load-bearing warnings and points at the
  `.m` file for the long rationale; the MATLAB headers remain the reference text.

## 2. The 2^53 wall is gone (deliberate, DECISION.md section 10)

MATLAB stores integers in doubles, raises `ratQ:overflow` at 2^53 and keeps a second BigInteger
path behind most kernels. The port has ONE path on unbounded ints:

- `rat_q.chk` / `chkB` check integrality only; `ratQ:overflow` is never raised by arithmetic.
- Every `try ... catch ratQ:overflow` fallback in MATLAB (keep the cell, leave a corner unnamed,
  skip the squarefree step, answer "+1") is therefore never taken in Python: on an input large
  enough to overflow in MATLAB, Python computes the exact answer instead, and the two results
  may differ in exactly those places (Python's is the decided one).
- Isolating intervals and bisection midpoints are `Fraction`s (dyadic), so the Sturm bisection is
  exact at any depth instead of stalling at double precision.

## 3. MATLAB defects found while porting -- all three FIXED in MATLAB on 2026-11-11

DECISIONS.md 2026-11-11 has the measurements. The port had the correct behaviour from the start.

1. **`ratQ.sturmChain` sign** (wrong root count when a chain element has a negative leading
   coefficient and the degree drops by two). MATLAB now scales by `|lc|^(gap+1)`.
2. **`conjQ/linealityConj` dropped the face denominator** (conjugated fD*f). MATLAB now passes
   `fN / fD`.
3. **`cellSamplePoint` accepted a sample on a boundary line** (zero-width cells kept as faces).
   MATLAB and Python now both require a row's value to clear its own rounding
   (32 eps times the sum of its terms' magnitudes) -- see section 6.

## 4. Float screens kept as in MATLAB

`conicMeet` (vertex coordinates), `QuaCon.eval`, `cellHasInterior`'s witness search (including the
seeded 400-point probe: MATLAB `rng(20260905); rand(400,2)` is reproduced by
`numpy.random.RandomState(20260905)`, column-major) and `cellInterior`'s margin LP are
double-precision in both. `numpy.roots` and MATLAB `roots` agree to rounding, not bitwise; every
decision taken from them is behind a tolerance (1e-9 .. 1e-7).

## 5. What is not ported yet (and what happens instead)

| MATLAB | status in Python | behaviour today |
|---|---|---|
| `envelopeSubdivisionAlg/minOverlay` (several atom pairs) | not ported | named refusal `envelopeSubdivisionAlg:overlayNotPorted`; `biconjQ` answers with the `EnvelopeMin` (same function, evaluable, not subdivided) |
| `triangleFace`, `cutTriangles`, `dualPoint` (contact triangles) | not ported | `ModuleNotFoundError` on the third branch of `biconjQ` -- to be turned into the same decline |
| `algRowFails`: `topRowPolys` / `topFaceRow` (rows against a face of degree >= 3 in z) | not ported | "no proof" (sound: the face is kept) |
| `algCrvVertices`, `algCrvIncidence` (naming an AlgCrv's vertices) | not ported | an `AlgCrv` has no `Vname` (as MATLAB's default result) |
| legacy `'cplq'` engine, `partialConj`, `moreau`, `proxAverage`, `infConv`, `lasryLions` | out of scope | `PLQ:conj:engine` names MATLAB |

## 6. Zero-width cells (resolved 2026-11-11)

Until 2026-11-11 MATLAB kept a cell whose sample point lay exactly on a boundary line (rounding gave
a positive slack on both sides of it) and emitted bridge faces for that zero-width cell; the port,
whose rounding fell the other way, did not. `cellSamplePoint` now discards a slack that does not
clear the row's own rounding, in both code bases, and the goldens were re-recorded.
`tests/test_trace_subdivision.py` still accepts an extra MATLAB face only after proving it lies in
such a cell -- with the fixed MATLAB that branch is no longer taken.

## 7. How the port is verified

1. **Golden parity** (`tests/test_golden_*.py`): every `conjQ` / `biconjQ` call the MATLAB suites
   make, replayed (`python/tools/export_golden.sh` records them through a shadow wrapper).
2. **Trace parity** (`tests/test_trace_*.py`): every recorded call of 25 intermediate MATLAB
   functions, replayed against its Python twin (`python/tools/export_trace.sh`).
3. **Definition oracle** (`tests/test_definition.py`, `cca2/check.py`): MATLAB-free.
4. **CAS-free lint** (`tests/test_no_cas.py`): no module of `cca2` imports a CAS.

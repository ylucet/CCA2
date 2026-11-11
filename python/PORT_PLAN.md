# Plan — port CCA2's exact conjugate / biconjugate engine from MATLAB to Python

_Planner: Opus 5.5, 2026-11-10. Executor: Fable. Repo: `/home/ylucet/AI/CCA2`, branch `main`._

## Context

**Why.** `DECISION/DECISION.md` (closed 2026-11-09) chose **Python** over Julia as the primary
language for CCA2. Both were measured on CCA2's own hard inputs and gave identical answers. Python
won on coverage (14/14 CAS ops vs 11/13 for Nemo, which lacks one-call `numden`/`collect`), on
port friction, and on ecosystem. Speed was a wash on the exact core (Python slightly faster: 31 µs
vs 34 µs `signExact`). Julia's "15 000× cheaper CAS" was **sympy vs FLINT**, a backend gap that
Python closes with **one C-core binding (python-flint, or cypari2/PARI)**. Licence: CCA2 is GPL-3,
so FLINT (LGPL) and PARI (GPL-2+) are fine and Sage (AGPL) is rejected. Bignums, not int64, in the
port: this **drops MATLAB's int64 wall** (a correctness gain). A first slice already exists:
`PORT/exactq_sign_py.py` (the `exactQ` sign core, 38/38 green against a sympy oracle).

**What gets ported.** The **exact engine**, meaning `QuaPol.conj` → `conjQ` (returns a `QuaCon`)
and `QuaPol.biconj` → `biconjQ` (returns an `AlgCon`: a `QuaCon`, `AlgCrv` or `EnvelopeMin`). Both are
the defaults since 2026-10-11. They are **essentially sym-free**: they use integer and rational
arithmetic through `ratQ`, `exactQ` and `bigZ`, and every decision is an integer predicate. So the
port is a translation of exact integer code. It is not a CAS port. sympy is needed only as a test
oracle, and python-flint only where it measurably accelerates `ratQ`'s polynomial kernels.

**Out of scope (stays in MATLAB, used as an oracle only):** the legacy `'cplq'`/`'symbolic'`
engine (`region.m`, `plq_1piece.m`, `plq_1p.m`, `symbolicFunction.m`, `functionNDomain.m`,
`QuaParCPLQ`, the symbolic branches of `conjCPLQ`/`biconjCPLQ`), `proveStage*`, `verify*`,
`sweep*`, `partialConj`, and the derived operators (`moreau`, `proxAverage`, `infConv`,
`lasryLions`). These can go in a later phase once conj/biconj are green. Also out of scope: the
Rust→C SCIP core that the decision names as the shipped artifact. It is a separate project that
this Python port will serve as reference for.

**Outcome.** A Python package `cca2` that returns the same exact answers as MATLAB `conjQ`/`biconjQ`
on every fixture in the exact-engine test suites. It refuses with the same named errors where
MATLAB refuses, and it passes a MATLAB-independent definition oracle.

## Size of the job (measured statically, 2026-11-10)

The function-level closure of `conjQ` is 11 files, ~6.3k lines. Of that, `ratQ` is 3.1k and
`conjQ` is 1.3k. The closure of `biconjQ` is 61 files, ~26k lines. The largest files there are
`maxQuaPar` (4.7k), `ratQ` (3.1k), `algRowFails` (1.5k), `envelopeSubdivisionAlg` (1.3k),
`conjQ` (1.3k), `biconjQ` (1.2k) and `bridgeFaceAlg` (0.9k). The classes add ~5.8k lines.
**Static grep overstates the closure**: comments and the legacy branches name many files. Phase 0
replaces it with a measured one. Files with `sym` calls inside the closure are `algRowFails`
(lines ~578 and ~1198), `pieceRecessionRays`, `convEnvUnbounded` and `envelopeAssemble`. Phase 0
must establish whether each one is on the exact default path or a legacy/oracle branch.

## Layout and toolchain

- New tree **`CCA2/python/`** in this repo (one repo, so the MATLAB oracle and the port stay in
  lockstep). Package `python/cca2/`, tests `python/tests/`, golden data `python/golden/`.
  `pyproject.toml`, GPL-3, Python ≥ 3.10 (host has 3.10.12).
- Runtime deps: **none beyond the stdlib** for the core (`int`, `fractions.Fraction`). `numpy` for
  float evaluation paths only. **`python-flint`** is an optional accelerator behind a thin
  `cca2/_poly.py` facade. Add it only when a kernel is measured hot, using operation counts
  rather than single wall-clock runs (umbrella §3). Use cypari2 only if flint lacks a needed op.
- Test deps: `pytest` (installed 9.1.1), `pytest-xdist` (`-n auto`), `sympy` 1.13.1 (oracle only;
  the compute path must never import it, and a lint test must enforce that, mirroring
  `conjSymFreeTest`/`checkConjSymFree.m`), `hypothesis` optional.
- Neither python-flint nor cypari2 is installed. `pip install --user python-flint` when needed.

## Translation rules (apply to every file)

1. **One MATLAB file → one Python module of the same name in snake_case** (`conjQ.m` →
   `cca2/conj_q.py`). Keep the MATLAB header comments: they hold the design rationale, so translate
   them and do not drop them. Add a docstring line `Ported from <file>.m @ <git sha>`.
2. **Indexing.** Use 0-based internally. Vertex, edge and face index columns inside `V/E/F`-style
   matrices are data: convert them **only** at the I/O boundary (`cca2/io_matlab.py`), never ad hoc.
3. **Integers.** MATLAB `int64` → Python `int`. MATLAB int64 `/` **rounds to nearest**, and Python
   `//` **floors**. Translate every integer division explicitly (`idivide` default is `'fix'`).
   `mod` ↔ `%` agree in sign, and `rem` does not. int64 **saturates**: any code that relies on
   saturation or `mulChecked` overflow becomes a no-op. Keep the guard sites as comments that
   point to the dropped wall.
4. **`bigZ`** (BigInteger via Java) → plain `int`. Its MATLAB-side "past 2^53" mode switches
   (`bridgePhi` 'big' mode, `ratQ` `*Big` variants, `implicitPolyBig`) collapse into the single
   int path. Port the **logic once**, not twice.
5. **`exactQ`** → start from `PORT/exactq_sign_py.py`. Do not re-port it. Keep its known
   subtleties: `exactQ(p,1)` is the rational `p`, not `√p`, and rational addition uses
   LCM `(d1//g)*d2`.
6. **Value semantics.** MATLAB classes here are value classes. Use `@dataclass(frozen=True)` or copy
   on write, because Python aliasing is the most likely silent bug.
7. **Named refusals.** `error('PLQ:conjQ:nonConvexPiece', ...)` → `raise CCA2Error('PLQ:conjQ:nonConvexPiece', msg)`.
   This is one exception class carrying the identifier string verbatim, so `verifyError` tests port
   to `pytest.raises(...)` plus an identifier check.
8. **Empties and shapes.** Make `[]`, `0x3` and `1x0` explicit, and copy the order and output of
   `unique`/`sortrows`/`find`/`sort` exactly. Determinism of face, edge and vertex order is what
   makes the golden comparison exact.
9. **Environment toggles** (`BICONJ_ASSERT`, `CCA2_*`) → same names via `os.environ`. `persistent`
   caches → module-level caches with an explicit clear function for tests.
10. **No fallbacks.** The rule in `conjQ.m`'s header applies to the port: answer exactly or raise a
    named refusal. Never silently degrade to floats.

## Phases (bottom-up; each phase ends green and committed)

### Phase 0 — measure the real surface and freeze the oracle (**needs MATLAB → UBC VPN**)
- Run the exact-engine suites with `bash .claude/suite.sh --fast --coverage` (Cobertura output is
  `.claude/coverage/cobertura.xml`). The **lines actually executed** under `conjQTest`, `biconjQTest`,
  `conjQExamplesSweepTest`, `erratumExamplesTest`, `erratumBiconjTest`, `biconjAtomsTest`,
  `bridgeFaceAlgTest` and the unit suites of each closure file are the port's true scope. Write the
  file list with per-file covered/total lines to `python/PORT_SCOPE.md`. Files with zero coverage
  on the exact path are dropped from the port.
- Write `python/tools/export_golden.m`. For every fixture those suites build, it serialises input
  `QuaPol` (exact `fN/fD`, `V`, `E`, `F`), output `QuaCon` (exact face N/D, primitive integer edge
  conics, named vertices, `dom` rows), and `AlgCon` kind with its exact data, to
  **`python/golden/*.json`**. Integers go out as decimal strings so nothing passes through double.
  For `AlgCrv`/`EnvelopeMin`, also export exact values at a set of rational sample points.
  Also export each named refusal: the fixture plus its error identifier.
- Copy the test fixture definitions (the `examples{k}` cell arrays etc.) as data, not by hand.
- **If the VPN is down:** skip coverage, use the static closure listed above, and build goldens
  later. Phases 1 to 3 can proceed on ported unit tests, which carry literal expected values, plus
  the definition oracle (Phase 1b). Record `BLOCKED: golden export` in `TODO.md`.
- Commit: `python/` skeleton, `PORT_SCOPE.md`, the exporter, and the golden JSON.

### Phase 1 — exact arithmetic kernel
- `exact_q.py` (from `PORT/exactq_sign_py.py`), `rat_q.py` (port of `ratQ.m`: rational and
  integer polynomial ops, Sturm/root isolation, `countRootsInHalfOpen*`, `isolateRoots*`, `signAt*`,
  `bigGcdRow`, `bigDivLinear`, collapsed to one int path per rule 4), `poly_sign_on.py`,
  `isolated_root_sign.py`, `dual_point_sign.py`.
- Port `ratQTest.m`, `exactQTest.m` and `bigZTest.m` (the last for whichever semantics survive) to
  pytest. Tests that pin int64 **overflow refusals** become tests that the value is now
  **computed**: assert it against sympy exact arithmetic. Record each one in `python/PORT_NOTES.md`
  as a deliberate divergence.
- **1b. Definition oracle**: port `plqCheck.m` to `cca2/check.py`. It checks `co f ≤ f` with
  equality at vertices, `f*(s) = sup_D <s,x> − f(x)` against a sampled sup (the sampled sup is a
  lower bound, so only `f* < sup_sampled` is a defect), and that `f** ≤ f` and is convex on
  segments. This oracle does not depend on MATLAB, so it covers VPN outages.

### Phase 2 — data types
- `domain.py`, `conic.py`, `crv.py`, `alg.py`, `rat_con.py`, `qua_con.py`, `alg_con.py`,
  `alg_crv.py`, `envelope_min.py`, and `qua_pol.py`. `qua_pol.py` covers only the exact subset
  `conj`/`biconj` reach: constructor, `readExactly`, `isExact`/`assertExact`, `faceQ`, `vertexQ`,
  `eval`, `createDom`, `orderEdges`, `isDomBounded`, `isConvex`, `assertOperable`, and `degree`.
  Skip plotting and the cubic-coefficient paths beyond what `assertOperable` needs.
- Use MATLAB inheritance only where the exact path dispatches on it (`RatCon` ⊃ `QuaCon`;
  `AlgCon` ⊃ {`QuaCon`, `AlgCrv`, `EnvelopeMin`} via `kind()`).
- Port `QuaPolExactTest`, `QuaPolEvalTest`, `QuaConTest`, `AlgConTest` and `AlgCrvTest`.
- Write `io_matlab.py` to load the golden JSON into these types and dump them back. The round trip
  must be identity, and a test pins it.

### Phase 3 — conjugate (`conjQ`) — **first user-visible milestone**
- `conj_q.py`, `max_q.py`, `conic_meet.py`, `cell_interior.py`, `cell_has_interior.py`,
  `compact_cells.py`, `merge_adjacent_cells.py`, `assemble_qua_con_cells.py`, `y_intercept.py`,
  plus any per-piece closed forms Phase 0 shows on the exact path (`conjConvexPolygon`,
  `conjAffinePLQ`, `conjPieceCPLQ`, `convEnvCPLQ`; their sym-free parts only).
- Wire `QuaPol.conj(engine='exact')`. Any other engine raises `PLQ:conj:engine` and names
  MATLAB for `'cplq'`.
- **Gate:** `conjQTest`, `maxQTest`, `conicMeetTest`, `cellInteriorTest`, `compactCellsTest` and
  `conjQExamplesSweepTest` are ported and green. Every conj golden must match **exactly**, field by
  field after canonical ordering. `check.py` must pass on every fixture. Commit, then update
  `TODO.md`.

### Phase 4 — biconjugate, direct cases
- `biconj_q.py` (classification, convex short-circuit, edge-concave lower hull with recession
  rays), `is_convex_q.py`, `domain_is_convex_d.py`, `merge_same_quad_faces.py`, `cut_triangles.py`,
  `piece_recession_rays.py` (only if Phase 0 shows it on the exact path; port its `sym` logic to
  `rat_q`), `to_qua_par.py` if reached.
- **Gate:** the convex, separable and edge-concave rows of `biconjQTest` are green, and
  `isConvexQTest` and `mergeSameQuadFacesTest` are green.

### Phase 5 — biconjugate, coupled pieces (atoms → envelope)
- `biconj_atoms.py`, `envelope_assemble.py`, `envelope_prune.py`, `envelope_eval.py`,
  `envelope_cell_pairs.py`, `envelope_cell_regions.py`, `envelope_rational_corners.py`,
  `EnvelopeMin`'s eval.
- **Gate:** `biconjAtomsTest`, `envelopeEvalTest`, `envelopeAssembleTest`, `envelopeCellPairsTest`,
  `envelopeCellRegionsTest`, and the `EnvelopeMin` rows of `biconjQTest`/`erratumBiconjTest`.

### Phase 6 — subdivision and algebraic faces
- `envelope_subdivision.py`, `envelope_subdivision_alg.py` (keep the `overlayBudget` option,
  default 3000), `alg_branch.py`, `alg_crv_vertices.py`, `alg_face_below.py`, `alg_row_fails.py`
  (port its two `sym` blocks to `rat_q` polynomial ops; this is where python-flint `fmpz_mpoly`
  may pay), `assemble_alg_crv.py`, `corner_face_alg.py`, `corner_wedge_rows.py`, `triangle_face*.py`,
  `dual_point.py`, `clip_arc_by_conic.py`, `clip_arc_by_half_plane.py`, `parabola_arc_frame.py`,
  `cell_sample_point.py`.
- **Gate:** the matching unit suites (`algBranchTest`, `algRowFailsTest`, `algFaceBelowTest`,
  `triangleFaceTest`, `cornerWedgeRowsTest`, `cornerFaceAlgTest`, `envelopeSubdivisionCornersTest`,
  `envelopeSubdivisionAlgTest`, `clipArcByHalfPlaneTest`, `frameAndFanTest`, `subdivisionPredicateTest`).

### Phase 7 — bridge faces
- `bridge_face_alg.py` (one int path replaces implicitPolyDouble/implicitPolyBig, discBigRaw and
  bigTail), `bridge_face_rat.py`, `bridge_phi.py`, `bridge_witness.py` (BigInteger Sturm → `rat_q`),
  `bridge_member.py`, `bridge_arc_cut.py`, `bridge_quad_cut.py`.
- **Gate:** `bridgeFaceAlgTest` and `bridgeFaceRatTest`, including
  `LARGEDataBuildsTheFaceWherePsiItselfPasses2to53`. While porting `bridgeFaceAlgTest`, also apply
  the known oracle fix: `tieMax`/`tiePoints` must scale their sweep window to the tie centre
  (`SESSION_HANDOFF` next step 1). Port the fixed version, not the buggy one, and note the change in
  `PORT_NOTES.md`.

### Phase 8 — end-to-end parity and the remaining fixtures
- Write `tests/test_golden_parity.py`, one parametrised case per golden. Exact equality applies to
  `QuaCon` outputs. For `AlgCrv`/`EnvelopeMin`, require equality of the exact values at the
  exported sample points plus `check.py`. For refusals, require the same identifier.
- Port `guardExerciseTest` for the exact-engine guards. Every guard `SUPPORT_MATRIX.md` lists as
  **OK** or **GAP** on the exact path gets a Python counterpart. Guards that existed only because of
  int64 are listed in `PORT_NOTES.md` as removed. Guards that are "not reachable" and
  "internal invariant" stay as asserts.
- Erratum Example 2 with `overlayBudget=Inf` is the verylong case: mark it `@pytest.mark.verylong`
  and exclude it from the default run.
- Test buckets mirror `suite.sh`: `pytest -m "not slow and not verylong" -n auto` = fast.

## Verification (definition of done)

1. `cd python && pytest -n auto -m "not verylong"` passes with 0 failures and 0 errors. The test
   count is reported against the MATLAB counts for the same suites.
2. Golden parity: 100 % of exported `conjQ`/`biconjQ` fixtures match, and refusals match by
   identifier. Any mismatch is either fixed or recorded as a deliberate divergence (int64 wall) in
   `PORT_NOTES.md` with the sympy-exact proof that Python is right.
3. `check.py` passes on every fixture. This is the MATLAB-free check.
4. The sym-free lint test passes: no `sympy` import reachable from `cca2/` outside `check.py`/tests.
5. Performance is **reported, not concluded on**: operation counts (Sturm sequences built,
   `conicMeet` calls), with timings as min-of-N on a shared machine (umbrella §3). There is no
   speed claim against MATLAB without that caveat.

## Working rules for the executor (from `AI/CLAUDE.md` and memory)

- Commit **on `main`** after each green step. Trailer: `Co-Authored-By: Claude <current session
  model> <noreply@anthropic.com>`, read from the session, never copied from here. **Do not push.**
- No parallel subagents. Use parallel processes (`pytest -n auto`, background MATLAB) freely.
- Ten strikes before parking a failure in `TODO.md` under **Blocked**, then move on.
- MATLAB runs only with the UBC VPN. If it is down, say so; never reason about what a MATLAB run
  would show.
- Before changing behaviour relative to MATLAB, grep `DECISIONS.md` headings and read the entry for
  that area. `SUPPORT_MATRIX.md` says whether a guard is a bug or by design. Record port-level
  decisions as new dated `DECISIONS.md` entries (e.g. "2026-11-xx — Python port: int64 wall dropped").
- Final messages to the user use the umbrella §4 status-line format.

## Critical files

- Decision: `DECISION/DECISION.md`, plus `DECISION/bench.py` (sympy usage patterns for the oracle).
- Seed port: `PORT/exactq_sign_py.py`.
- Entry points: `QuaPol.m` (`conj` ~l.764, `biconj` ~l.814), `conjQ.m`, `biconjQ.m`.
- Kernel: `ratQ.m`, `exactQ.m`, `bigZ.m`. Oracle: `plqCheck.m`. Buckets: `.claude/suite.sh`.
- Docs to grep: `ALGORITHM.md` (order of operations), `SUPPORT_MATRIX.md`, `DECISIONS.md`,
  `RETURN_TYPE.md`, `CONJ_FIELD_PROOF.md` §8.6/8.7.

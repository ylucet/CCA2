# CCA2 manuscript — working outline

_Drafted 2026-11-17, on top of main `1f9a056`. Target venue: **Mathematical Programming
Computation** (backup: COAP; SVVA only if the paper ends up led by the mathematics — the
allocation is recorded in `../convexdb/manuscript/PLAN-DECISIONS.md` "Cross-project
allocation" and in `TODO.md` 2026-09-03). This file is a plan for the manuscript, not the
manuscript; every claim below is tied to a file or test in this repository._

---

## 1. Working title

Candidates, in preference order:

1. **"Exact Conjugates and Convex Envelopes of Bivariate Piecewise Quadratic Functions:
   Algorithms, Proofs, and a Verified Implementation"**
2. "The Fenchel Conjugate and the Convex Envelope of Bivariate PLQ Functions in Exact
   Closed Form"
3. "Beyond Parabolic Subdivisions: Exact Conjugacy for Nonconvex Bivariate Piecewise
   Quadratic Functions"

Candidate 1 is the default: it names the two operators, the input class, exactness, and the
implementation + verification story that makes this an MPC paper rather than an analysis note.
Candidate 3 is the stronger claim (the face-type result below) and the right choice only if
the paper is re-scoped toward the mathematics (which is what would move it to SVVA).

## 2. The question, and why it matters

**The question.** Given a bivariate piecewise linear-quadratic (PLQ) function `f` — possibly
nonconvex, possibly on an unbounded polyhedral subdivision — compute its Fenchel conjugate
`f*` and its biconjugate `f**` (the closed convex envelope) **exactly**, as closed-form
piecewise objects, not by grid sampling, LP relaxation or symbolic computation.

**Why it matters.**

- **Convex envelopes of nonconvex bilinear/quadratic terms are the workhorse of global
  optimization.** On QPLIB's box-domain bilinear terms the envelope is McCormick
  (Al-Khayyal–Falk), which CCA2 reproduces as a special case — see
  `SCIP_READINESS.md` §0 and the pin
  `biconjugateTest/bilinearOverABoxGivesTheMcCormickEnvelope`. CCA2's interest is the part
  McCormick does not cover: non-box domains and diagonal terms such as `x² − y²`
  (`SCIP_READINESS.md`, the reading that shapes the whole readiness plan).
- **The published algorithms for this exact problem exist, but their implementation path does
  not match the computation they describe.** The linear-time algorithms of [COAP] and [JOGO]
  (Karmarkar & Lucet, 2026, both cited in `README.md` "Algorithms and references") are
  implemented in the repository as the legacy `'cplq'` engine, which runs on the Symbolic
  Math Toolbox and whose heaviest cases take more than four hours
  (`biconjCPLQ` on the poly2 kinks runs > 4 h — `TODO.md` "Open now"). Even the exact
  MATLAB reference is far from the ports on those same fixtures — 412–566 s per single
  slow fixture (`python/bench/RESULTS.md`) against 5.4–7.6 s in Rust
  (`rust/bench/RESULTS.md`). The exact engine this paper is about exists to remove the
  symbolic dependency from the compute path entirely — the plan noted in
  `../convexdb/manuscript/PLAN-DECISIONS.md` item 3, and now landed.
- **The answer's shape is richer than the published theory says.** The conjugate and
  biconjugate of a continuous PLQ function with rational data do not in general admit a
  parabolic subdivision, and their vertices are not in any tower of quadratic extensions:
  an explicit three-piece counterexample refutes [JOGO] Theorem 6 and Proposition 7
  (`doc/QuaConExample.md`), and a five-piece continuous example carries a degree-4 vertex
  with Galois group `S4` (`CONJ_FIELD_PROOF.md` §4–5). The face type of the biconjugate is
  algebraic of degree ≤ 4, and `f**` is not `C1` (`DECISIONS.md` 2026-08-23). An exact
  implementation therefore needs a data structure (rational face functions, rational edge
  conics, degree-≤4 algebraic vertices, implicit algebraic faces) that no grid method
  expresses at all — and the paper's algorithmic core is exactly that structure plus the
  machinery that fills it.

## 3. Contributions, each tied to the files and tests that establish it

### C1 — New exact algorithms for the conjugate and biconjugate of bivariate PLQ functions

- **What it is.** A direct, closed-form, exact-arithmetic algorithm for `conj(f)` (the
  Fenchel conjugate of a `QuaPol`, convex or nonconvex) and for `biconj(f)` (the closed
  convex envelope, computed directly rather than as `conj ∘ conj`), returning typed exact
  objects: `QuaCon` (piecewise quadratic, rational data), `AlgCrv` (carrying curved
  algebraic bridge faces) and `EnvelopeMin` (min-of-maxes when several atom pairs are
  active). The order-of-operations design — classify first, split last, prefer the direct
  operator — is in `ALGORITHM.md`, with the measured effect of getting the order right
  (456 s → 4.2 s on `(x²+y²)/2` over the unit square; `x² − y²` from an error to exact in
  1.9 s).
- **The mathematics behind it.** Face functions and edge conics of `f*` are always rational
  (Theorem 1 of `CONJ_FIELD_PROOF.md`); vertices have degree ≤ 4, `S4` when the three
  meeting pieces are pairwise non-adjacent; the biconjugate's ruled-cell faces are algebraic
  of degree ≤ 4 (`DECISIONS.md` 2026-08-23); the three-piece counterexample to [JOGO]
  Theorem 6 is computed in full in `doc/QuaConExample.md` with its proof.
- **Files.** The exact engine in `MATLAB/`: `conjQ.m`, `biconjQ.m`, `maxQ.m`,
  `assembleQuaConCells.m`, `biconjAtoms.m`, `envelopeSubdivisionAlg.m`, `algRowFails.m`
  (the emptiness certificate, incl. the `psTangency` perfect-square tangency rule,
  `DECISIONS.md` 2026-11-17), `bridgeFaceAlg.m` / `bridgeFaceRat.m`, `cellInteriorTest.m`
  and companions; type classes `QuaCon.m`, `AlgCrv.m`, `AlgCon.m`; kernels `ratQ.m`,
  `bigZ.m`.
- **Tests that establish it.** MATLAB suites, measured on the shared machine (Frances):
  fast bucket **834 / 0 over 64 suites** (`DECISIONS.md` 2026-11-17), normal **35 / 0**,
  and the eight exact slow suites green (kinkSubdivisionTest 3/0, minSubdivisionTest 4/0,
  erratumBiconjEx2Test 4/0, algCrvSubdivisionTest 3/0,
  envelopeSubdivisionSlowTest 1/0, algRowFailsSurdTest 5/0, unboundedFaceTest 19/0,
  guardExerciseSymbolicTest 11/0 — `DECISIONS.md` 2026-11-11). The per-stage cache for
  test splitting is `plqStage.m`; the definition verifiers are `plqCheck.m`
  (co f ≤ f with equality at vertices; `f*(s)` against a numeric sup; `f** ≤ f` and convex
  on sampled segments).

### C2 — The exact engine: exact rational arithmetic with a verification harness that owes nothing to it

- **What it is.** All exact decisions — face liveness, row-dropping, cell interiors, bridge
  witnesses, the emptiness certificate `algRowFails` — are decided in exact arithmetic
  (rationals via `ratQ.m` with the `bigZ.m` BigInteger fallback; FLINT in the ports), with
  doubles used only for sampling and for the polynomial kernels that have a symbolic
  oracle. The Symbolic Math Toolbox is out of the compute path: in the benchmark it is
  reached only in `algRowFails` fallbacks on `Yun` overflow over 2^53, "not reached here"
  (`python/bench/RESULTS.md`, "What ran").
- **Verification harness.** Four independent oracles, none of which shares code with the
  engine it checks:
  1. **Golden replay.** Every `conjQ` / `biconjQ` call the MATLAB suites make is recorded
     as exact-rational JSON lines (`python/golden/*.jsonl`, 9 files, 2912 lines: 176 conj +
     114 biconj records plus the slow-suite and trace files) and replayed by the other two
     engines.
  2. **Definition oracle.** `plqCheck.m` (MATLAB) and its twin `python/cca2/check.py`,
     which computes `sup <s,x> − f(x)` pointwise from KKT candidates and shares nothing
     with `conjQ`; all **184 fixtures** satisfy the definition in Python
     (`python/tests/test_definition.py`, `TODO.md` P1b).
  3. **Kernel oracle.** 186 kernel tests against a SymPy oracle with a lint test that keeps
     SymPy out of the production path (`python/tests/test_kernel_oracle.py`,
     `test_no_cas.py`).
  4. **Cross-engine parity.** MATLAB / Python / Rust agree field-for-field on the goldens
     and on digests of all 200 benchmark fixtures (`rust/tests/parallel.rs`,
     `rust/README.md`).
- **Why it is a contribution in its own right.** The engines disagreeing on the star-polygon
  class — then being traced to one shared heuristic, the incomplete cell-interior witness
  search (`DECISIONS.md` 2026-11-17, `python/star_defect_diagnosis.py` on the pending
  branch `feature/17-cca2-star-polygon-28u`) — is only *visible* because the harness
  exists. The paper should present the harness as the reason the correctness claims are
  checkable, which is precisely what MPC's verification process rewards.

### C3 — Three-language implementations and their speed-ups

- **What it is.** A module-for-module Python port of the exact engine
  (`python/cca2/`, `PORT_PLAN.md`, `PORT_NOTES.md`) and a Rust port of the Python port
  (`rust/`, `rust/README.md`), both at parity with the MATLAB reference on the recorded
  corpus, plus a rayon parallelization of the Rust engine (`rust/PARALLEL.md`) with a
  thread sweep from 1 to 64 threads.
- **Established by.**
  - Python: **2730 passed / 1 skipped** fast bucket (`MORNING.md` 2026-11-10/11;
    `pytest -q -n 32`); golden parity 176/176 conj, 114/114 biconj (refusals by
    identifier), slow suites 20/20 (`rust/README.md` table, which covers the same traces);
    the three MATLAB defects the port surfaced and fixed are recorded in
    `DECISIONS.md` 2026-11-11 (Sturm-chain sign, zero-width cells, `linealityConj`
    conjugating the numerator instead of the fraction) with red-first tests.
  - Rust: `cargo test --release` replays every golden record against the Python reference
    with and without FLINT; the parallel tests assert digests of all 200 fixtures are
    equal at 1, 2, 8 and 32 threads (`rust/tests/parallel.rs`); `conjstress` (star
    polygons: fold vs tree vs the definition, commit `9355207`) is the tool that found the
    open star-polygon defect.
- **Speed-ups.** The 200-fixture protocol of `python/bench/` (fixtures in
  `python/bench/fixtures.json`, derived from the goldens; best of repetitions, memos
  emptied per call):

  | engine | 200 fixtures | vs MATLAB |
  |---|---|---|
  | MATLAB R2024a (reference) | 1126 s | 1× |
  | Python (numba + FLINT) | 88.5 s | 12.7× |
  | Rust, 1 thread (FLINT) | 15.7 s | 72× |
  | Rust, 32 threads | 2.22 s | ≈ 500× |

  Heaviest fixtures, 32 threads: kink 0.81 s, erratum Example 2 0.70 s
  (`rust/PARALLEL.md` sweep). **These numbers are from a shared machine and the runs were
  on different days under different loads — the repository's own caution
  (`rust/bench/RESULTS.md`, "read the ratios as indicative") carries into the paper; see
  Gap G6.**

## 4. Section plan

Sized for MPC's roughly 20-page body plus references. Sections that already have their
substance in the repository are marked; the rest is writing.

1. **Introduction.** Bivariate PLQ, Fenchel conjugate and convex envelope, the nonconvex
   case, why exact closed form (relaxations for global optimization; McCormick as special
   case — `SCIP_READINESS.md` §0). The one-paragraph statement of the three contributions.
2. **The problem and prior work.** [COAP] and [JOGO] (the linear-time algorithms, their
   symbolic implementation path); the convex-only predecessors and the mapping to engines
   they would supply (`DESIGN.md` "PART I" and the lineage section); the scope decision —
   `QuaPol` conjugate/biconjugate is the design target and what is out of scope by design
   (`SUPPORT_MATRIX.md` §0, the OK / GAP / N/R / INV legend).
3. **What the answer looks like: face types and number fields.** Rational face functions
   and edge conics (Theorem 1, `CONJ_FIELD_PROOF.md` §2); degree-≤4 vertices and the
   `S4` five-piece example (`CONJ_FIELD_PROOF.md` §4–5); the biconjugate's algebraic
   faces, degree ≤ 4, `f**` not `C1` (`DECISIONS.md` 2026-08-23); the three-piece
   counterexample to [JOGO] Theorem 6 (`doc/QuaConExample.md` — computed, with proof,
   exact rational arithmetic). This is the section that justifies the data structure.
4. **The exact conjugate.** The classify-first order (`ALGORITHM.md`); per-piece closed
   forms and the A.2–A.5 split by convex-edge count; `maxQ` assembly and the exact
   evaluator (every candidate decision the sign of a rational, `TODO.md` 2026-09-06
   entry); the cell-interior screen and its known incompleteness (stated here, with the
   star-polygon defect — see Gap G1).
5. **The exact biconjugate and convex envelope.** The direct algorithm (never
   `conj ∘ conj`); atoms, envelope assembly and pruning, bridge faces (algebraic and
   rational), the local-arrangement subdivision and the min-of-maxes overlay
   (`DECISIONS.md` 2026-10-18 → 2026-11-05 sequence); the emptiness certificate
   `algRowFails` incl. the `psTangency` tangency rule (`DECISIONS.md` 2026-11-17); output
   types `QuaCon` / `AlgCrv` / `EnvelopeMin` and the refusal identifiers.
6. **The exact engine and its verification.** Exact arithmetic and where doubles appear;
   the four oracles of C2; the test statistics; the three MATLAB defects the port found
   (`DECISIONS.md` 2026-11-11) as a case study of why the harness matters.
7. **Implementations and performance.** MATLAB reference, Python port, Rust port and
   parallelization; the 200-fixture protocol; the speed-up table of C3 with its shared-
   machine caveats; the FLINT kernel (one bivariate-polynomial part, `rust/src/flint.rs`);
   reproducibility of the benchmark from the repository alone.
8. **Case studies.** A worked end-to-end example (primal, `f*`, `f**`, with figures —
   the material for it exists: `doc/QuaConExample.md`, `doc/fig_primal.pdf`,
   `doc/fig_dual.pdf`, `doc/QuaConExample-facemap.png`, `doc/QuaCon.svg`); the kink
   fixture and the erratum Example 2 as the heavy cases; the star-polygon class reported
   honestly, including the known-wrong regime (Gap G1).
9. **Code and data availability.** MPC compliance, how to run each implementation, the
   data files, DOIs, the declared limitations (section 7 of this outline).
10. **Conclusion and open problems.** The exact cell-interior test (the star-polygon
    fix), `partialConj`, generalisation beyond bivariate, the SCIP/QPLIB consumer
    (`SCIP_READINESS.md`, `SUPPORT_MATRIX.md` §0.0.1).

## 5. MPC code-and-data requirements, and how CCA2 meets them

From the venue record (`TODO.md` 2026-09-03; `../convexdb/manuscript/PLAN-DECISIONS.md`):
MPC requires **mandatory code-and-data submission with independent verification of
computational results** — its results verification assumes a referee can run the code.
Free route: Springer's 12-month green (self-archive the accepted manuscript); no APC.
The MATLAB dependency is the real friction and is to be **declared in the cover letter**,
offering MPC's "detailed observation / guest access" route rather than letting a referee
discover it.

**What a referee needs, stated plainly:**

- **To run the reference implementation and the MATLAB suites: MATLAB, academic licence,
  no toolbox required for the compute path.** The exact engine (`conjQ`, `biconjQ`) and
  the whole fast/normal bucket run on base MATLAB R2024a; exact integers are
  double-precision to 2^53 then `java.math.BigInteger` (`bigZ.m`). The Symbolic Math
  Toolbox is touched only in `algRowFails` fallbacks on overflow (not reached in the
  benchmark) and by the legacy `'cplq'` symbolic suites (the slow bucket), which are
  presented as a cross-check, not as the result. (Authors' own machines additionally
  route the licence server through the UBC VPN — an author-side constraint, not a
  referee-side one, but it is why the MATLAB-free route below must be first-class.)
- **To run without MATLAB at all: Python 3 or Rust.** The Python port
  (`python/cca2/`) is complete for `conjQ` / `biconjQ` and its test suite is
  MATLAB-free: golden replay from the JSON data files, the KKT definition oracle
  (`python/cca2/check.py`), the SymPy kernel oracle (test-only) and 2730 green tests.
  The Rust port (`rust/`) replays the same goldens against the Python reference and
  carries the 200-fixture digest tests; FLINT is optional (a slower pure-Rust kernel
  gives identical answers, `rust/README.md`).
- **The data files.** The test corpus is **plain JSON data, not MATLAB code**:
  `python/golden/*.jsonl` (9 files, 2912 lines, exact rational coefficients) and
  `python/bench/fixtures.json` (the 200 distinct exact inputs). They are the artifact's
  data component: a referee can verify parity and replay the benchmark with neither
  MATLAB nor the paper's engines.

**How CCA2 meets the requirements, item by item:**

| MPC requirement | How CCA2 meets it | State |
|---|---|---|
| Code available, runnable by a referee | Three implementations, two of them MATLAB-free; base MATLAB for the third; declared in cover letter | **met** |
| Data available | Exact-rational JSON corpus + 200-fixture bench protocol, plain text in the archive | **met** |
| Independent verification of results | Four oracles sharing no code with the engine; three engines agreeing at parity; refusal identifiers stable across engines | **met** |
| Reproducible computational results | Benchmark protocol + scripts committed (`python/bench/`, `rust/bench/`); timings caveated (shared machine) | **met, caveated** |
| Citable artifact | — | **gap: G5** (no release tag, no code DOI yet) |

## 6. Gaps between what exists and what the paper needs

Ordered by how much they block a submission.

- **G1 — The star-polygon conjugate defect (OPEN, diagnosed 2026-11-17).** `conjQ` of a
  star polygon (2n vertices, `f = (x²+y²)/2`) returns a value that is **too low** at some
  dual points: Rust wrong at n=6 (111/625 grid points), n=7 (173/625), n=8 (3/625);
  Python wrong at n=7 (25/625); the engines disagree with each other and with the
  definition on this class. Root cause localised: the cell-interior witness search
  (`cellHasInterior`) is an incomplete heuristic (conic/line intersections + 400 random
  probes) that misses thin cell interiors, so `_dropImplied` drops a bounding row and the
  evaluator picks the lower face. A class-based refusal is not viable — it breaks the
  correct L-shape (golden conj record 36). The proper fix — a complete exact feasibility
  test for the conic sign conditions, in all three engines — is specified in
  `DECISIONS.md` 2026-11-17 (recorded on the pending branch
  `feature/17-cca2-star-polygon-28u`, reproduction script
  `python/star_defect_diagnosis.py`). **The paper either fixes this before submission or
  states the affected input class as a declared limitation; it cannot claim exact
  conjugacy for nonconvex bounded polygons while this is open.**
- **G2 — `getInterior` on a singular quadratic (OPEN, `SUPPORT_MATRIX.md` §7, LIVE since
  2026-08-16).** On a piece where `f` is constant along an edge, the gradient map is not
  invertible and two edge cells come out on the same region: the parallelogram's `f**`
  is about 4% LOW at 2 of 10 interior probe points (`0.986` for `1.031`, `1.913` for
  `1.950`). Reproducible in about a minute (`functionNDomainTest.parallelogramPiece9`).
  Smaller than G1 but it is a wrong value in the same operators the paper is about.
- **G3 — Missing comparisons.** There is no comparison against an independent external
  implementation or against grid/LP-based envelopes in the repository. In-repo, the
  legacy `'cplq'` symbolic engine is the [COAP]/[JOGO] implementation and is available as
  a same-problem comparison, but it is not benchmarked against the exact engine on
  accuracy or speed (only the exact-vs-port speed-ups exist). And the `convexdb`
  ground-truth registry — the intended external oracle set — has **no Zenodo DOI yet**
  (`TODO.md` "Release / paper"; the DOI must be minted before any of the four group
  manuscripts is submitted, `../convexdb/manuscript/PLAN-DECISIONS.md` sequencing item 1).
- **G4 — Missing examples and figures.** The paper needs one end-to-end worked example
  (primal + `f*` + `f**`, figures, narrative) and a publication figure set (face-type
  taxonomy; the three- and five-piece counterexamples; the kink and Example 2; the
  speed-up chart). Fragments exist — `doc/QuaConExample.md` + its proof,
  `doc/fig_primal.pdf`, `doc/fig_dual.pdf`, `doc/fig_hyperbola_insets.pdf`,
  `doc/QuaConExample-facemap.png`, `doc/QuaCon.svg`, `doc/full_pipeline.svg` — but
  nothing is assembled as the paper's figure set, and the worked example is not written
  for a reader who does not know the code.
- **G5 — The citable artifact.** CCA2 is alpha with no tagged release (`README.md`
  "Status"); MPC wants a code archive with a DOI for the submitted version, and the
  `convexdb` citation it carries needs the DOI from G3. Neither exists. The API is also
  still changing (`RatPol.conj` / `biconj` / `add` missing, `partialConj`
  unimplemented — `README.md` operators table; `SUPPORT_MATRIX.md` §2, §8).
- **G6 — The timings are indicative, not final.** All quoted speed-ups come from a
  shared machine, on different days, under loads 0.5–16 (`rust/bench/RESULTS.md`,
  `rust/PARALLEL.md`). The umbrella rule (`~/AI/CLAUDE.md` §3) forbids putting a single
  timed run into a manuscript: the 200-fixture benchmark needs a controlled re-run
  (min-of-repetitions, load recorded) before its table appears in §7.
- **G7 — Residual engine items that must at least be disclosed.** `Example 2` returns an
  `EnvelopeMin` by default (overlay budget 3000; lifting it takes ~hours — `TODO.md`
  "Open now"); two `BRIDGE` faces owning one cell are refused by name, never observed;
  the legacy `'cplq'` engine carries its own open items (the > 4 h poly2 kinks, the `eig`
  tolerance that can return `+Inf` on a strictly convex quadratic, G1/G4/G10/G11 and
  companions — `TODO.md` "Open now"). None blocks the exact-engine claims, but the
  paper's availability section must say what the legacy path is and is not.
- **G8 — The independent geometric oracle is planned, not built.** The planned
  `CCA2-oracle` sibling (a CGAL-based cross-check generator, deliberately outside this
  repository so the submitted artifact stays "MATLAB plus archivable
  exact-rational data" — `../convexdb/manuscript/PLAN-DECISIONS.md`) does not exist as a
  project yet. Today's oracles are definition-based, kernel-based and cross-engine — all
  sharing the repository's own definitions. The paper should say so rather than imply
  geometric independence.

**The three biggest gaps, in one line each:** G1 (a known wrong-value class in the exact
conjugate — star polygons — diagnosed but unfixed in all three engines), G3 (no external
comparison or DOI'd ground-truth set, so the correctness claims rest on in-repo oracles
alone), G5 (no tagged release or code DOI — the citable artifact does not exist yet).

---

## Notes for the writing session

- Dates and measurements carry their machine caveats wherever they appear
  (`DECISIONS.md` entries 2026-11-11, 2026-11-17; `python/bench/RESULTS.md`;
  `rust/PARALLEL.md`).
- Citations: [COAP] / [JOGO] per `README.md`; `convexdb` by its Zenodo DOI once minted
  (not a repository URL); the master bibliography is `../biblio/biblio.bib`.
- The star-polygon diagnosis lives on the pending branch
  `feature/17-cca2-star-polygon-28u` (commit `a6da6de`); merge it before drafting §4/§8
  so `DECISIONS.md` and `TODO.md` in the submitted tree already carry the record.

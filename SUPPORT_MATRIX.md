# CCA2 — Supported / unsupported matrix

**Generated from the error guards in the source**, not from `DESIGN.md` (which is a design
*proposal* and describes intended as well as built code). Every "not supported" row below
corresponds to a specific `error(...)` call, cited by file and identifier, so this file can be
re-derived and re-checked mechanically.

**Partially refreshed 2026-08-14** (§4 and its sub-sections, §7's header, §8's ordering): the
arc-vs-arc rows and the far-field defect were stale, and §4's guard line numbers predated that work
by ~1400 lines. All guards and line numbers re-derived from the source.

Re-measured that day, every bucket: `maxQuaParTest` **28 / 0** (was 25 / 1), fast **203 / 0** (was
200 / 1), normal **6 / 0**, slow **111 / 4** — the slow bucket identical to its recorded baseline,
with all four failures being the documented open items of §7 and §8. No regression anywhere. The
full 26-suite census below stands at its own date. Treat every count in this file as carrying its
own date.

Last regenerated: 2026-08-02, against a full census of **300 passed / 1 failed / 0 incomplete**
across 26 suites (`CCA2_TEST_TIMEOUT=5400 bash .claude/suite.sh`, run against a snapshot of the
tree so that concurrent editing could not contaminate it). The single failure is
`biconjugateTest/biconjugateOverATwoFaceSubdivisionIsTheEnvelope`, left failing deliberately: it
pins the open defect in §7 and its comment carries the traced mechanism.

`testRegion/testCreation`, recorded here since 2026-07-28 as a longstanding toolbox-compatibility
failure, now passes — `testRegion` is 23 / 0.

Comparison runs behind that census, all on snapshots of the same revisions:

| tree | pass / fail / incomplete |
|---|---|
| `46fac7c` (mid-session) | 298 / 2 / 1 |
| `46fac7c` + the substitution batching, nothing else | 298 / 2 / 1 — **byte-identical per suite**, which is the evidence that the batching changed no answers |
| final | **300 / 1 / 0** |

The two suites that differ between `46fac7c` and the final tree are `QuaParTest` (10 → 11, the new
`evalLocatesAPointExactlyAtItsOwnVertex`) and `testMaxMultiRegion` (23/1/1 → 24/0/0, a regression
`46fac7c` introduced and the commit after it repaired — see that commit's message).

---

## 0.0 Downstream consumers — check before renaming anything

**The SCIP feasibility spike (`AI/spike/SCIP`) calls into this toolbox** through the MATLAB
Engine API for Python. Its glue is `SCIP/src/cca2ConvexEnvelope.m`, and it uses exactly one entry
point: **`convEnvCPLQ`** (the convex ENVELOPE), consuming the resulting `RatPol`'s
`V/E/F/f(:,5:10)/den` as plain arrays. It does **not** call `conj`, `biconj`, `partialConj`, or
any `QuaPar` method — per-node cut math is pure Python after one bridge call.

This matters twice over:

1. **`QuaPol.m`'s 2026-07-27 rename note is wrong.** It says no compatibility shim was left for
   `QuaPoly` "because CCA2 has no tagged release, so nothing external can depend on it". SCIP's
   glue constructs `QuaPoly(Vin, Ein, fin, Fin)`, so the rename broke that bridge silently.
   `QuaPoly.m` is now an alias, as `PLQVC.m` already was. **Do not remove it, and check this
   section before renaming any other public name.**
2. **The blocker list below is largely irrelevant to SCIP.** Every open item in section 8 lives
   in the `conj`/`biconj` pipeline. The convex-envelope path SCIP actually uses is the mature one
   — bounded faces, no `QuaPar.eval`, no cPLQ Step 3.

---

## 0.0.1 What a SCIP + QPLIB run would need from CCA2 — measured 2026-08-18

_Re-derived from `checkBoxEnvelopeForSCIP`'s own output on the current tree (log:
`.claude/boxenvelope.log`). The previous version of this section was dated 2026-08-02 and every
row of it is now stale: two ERROR rows are gone, the times went from 38–77 s to 0 s, and the
"returns no mesh" headline no longer holds for any box case._

QPLIB's viable family for this work (Sahinidis `QPLIB_1913/1922/1931/1940`) has **unit-BOX**
variable domains, so per term SCIP needs the convex envelope of a quadratic **over an axis-aligned
box**, and over a sub-box after each branching. That is `f**`, NOT `convEnvCPLQ`'s output —
`convEnvCPLQ` is Step 1, a per-TRIANGLE intermediate, and the SCIP spike established on 2026-07-11
that cuts built from it over a box domain are **invalid** (they cut off feasible points). That
finding stands.

| term over an axis-aligned box | result | worst error | returns | time |
|---|---|---|---|---|
| `x·y` on `[0,1]²` | **OK** | 5.6e-17 | `QuaPol`, meshed | 0 s |
| `3xy + 7x − 2y + 5` on `[0,1]²` | **OK** | 0 | `QuaPol`, meshed | 0 s |
| `x² − y²` on `[0,1]²` | **OK** | 1.6e-4 | `QuaPol`, meshed | 0 s |
| `(x²+y²)/2` on `[0,1]²` | **OK** | 1.6e-4 | `QuaPol`, meshed | 0 s |
| `x·y` on `[0.25,0.75]×[0,1]` (post-branching sub-box) | **OK** | 1.1e-16 | `QuaPol`, meshed | 0 s |
| `x·y` on `[−2,3]×[−1,4]` (general bounds) | **OK** | 0 | `QuaPol`, meshed | 0 s |

The two 1.6e-4 rows are the CHECK's own discretisation, not the envelope's: those two envelopes
are curved (`x² − y` and `(x²+y²)/2`), and the reference is the lower convex hull of the graph
SAMPLED on a grid, which sits slightly above a curved surface between samples. The four rows whose
envelope is piecewise affine come out at machine precision, which is what says the comparison is
sound. The verdict threshold is 2e-3.

**Why it is now 0 s: the short-circuits, not a faster pipeline.** `biconj` IS the closed convex
envelope operator, so it answers these shapes from closed forms before any conjugation happens —
convex `f` is its own answer, separable-on-a-box drops to one 1-D envelope per axis, bilinear on a
box is McCormick / Al-Khayyal-Falk, bilinear on a diamond rotates to the separable case
(`biconjCPLQ.m`; the ordering and the reasons are in `ALGORITHM.md`). Since 2026-08-18 a Step 0
normalisation (`mergeSameQuadFaces`) also merges adjacent faces carrying the same quadratic, so a
caller who hands the same box in as two triangles reaches the same short-circuits.

So the gaps, in the order they would bite, all four re-checked:

1. **No entry point — unchanged, and it is now a one-liner.** The bridge (`SCIP/src/cca2ConvexEnvelope.m`)
   calls `convEnvCPLQ`; what it wants is `q.biconj('cplq')`, which returns the box envelope in
   0.01 s (`SCIP_READINESS.md` B2 names it as the entry point). Wiring, not mathematics — and it
   is on the spike side of the fence, so it has not been re-read here.
2. **~~`biconj` returns a `QuaParCPLQ`, which has NO MESH.~~ FIXED for every box case.** All six
   rows above come back as a **meshed `QuaPol`** with `V/E/F/f` populated, which is exactly what
   the bridge reads. The mesh-less `QuaParCPLQ` is still what the GENERAL (non-box, multi-piece)
   route returns, so the value+subgradient escape hatch that `RETURN_TYPE.md` records is still
   the plan for that path — just not for SCIP's.
3. **~~A quadratic with DIAGONAL terms over a box fails.~~ FIXED.** `x² − y²` and `(x²+y²)/2` both
   return in 0 s: the first is SEPARABLE as written (its envelope is `x² − y`), the second is
   CONVEX and is its own envelope. The `MATLAB:badsubscript` came from taking the long route into
   the second conjugation, which neither now reaches. Pinned by
   `biconjCPLQTest.separableOverABoxTakesTheOneDimensionalRoute` and
   `convexOverABoxIsItsOwnBiconjugate`.
4. **~~40–60 s per term is the practical blocker.~~ NOT a blocker on the box path.** At 0 s per
   term, `QPLIB_1940`'s 288 objective off-diagonal terms are seconds and its ~27,586 constraint
   terms are minutes, offline. Per-node recomputation is no longer obviously out of reach either.
   What still costs minutes is the NON-box case (a polygon domain, whose envelope couples several
   pieces and reaches Step 3's symbolic maximum: 43 min for one quadrilateral term).

**The caveat that decides whether this is worth doing has NOT changed, and is now sharper.** For a
BILINEAR term over a BOX the convex envelope *is* McCormick, in closed form (Al-Khayyal–Falk) —
which is exactly what CCA2 returns, and `biconjugateTest/bilinearOverABoxGivesTheMcCormickEnvelope`
pins it. CCA2 is therefore a 0-second reimplementation of a formula SCIP already applies: correct,
a good validation, and **no stronger as a cut**. CCA2 can only beat McCormick where the domain is
not a box (a triangle or polygon cut out by constraints) or the piece is not bilinear — and that
is precisely the case that still costs minutes. Being fast on the box family removes the
performance objection to a QPLIB run; it does not make one scientifically interesting.

---

## 0.2 Where the cPLQ code lives — `cPLQ/` is NEVER executed

Verified 2026-08-02, because "are we running vendored cPLQ?" is a fair question and the answer is
not obvious from the file listing:

* **`cPLQ/` is a local reference clone only.** It is **gitignored** (`.gitignore:7`), so 0 files
  under it are tracked; it is **never added to the MATLAB path** (no `addpath` anywhere mentions
  it); and `which('region')` / `which('plq_1p')` resolve to `CCA2\region.m` / `CCA2\plq_1p.m`.
  **No computation uses it.** It exists to diff against.
* **The working copies are CCA2's own, in the repository root**, and they are merged and
  substantially rewritten rather than vendored-as-is:

  | file | pristine cPLQ | CCA2 root | changed lines | CCA2 commits |
  |---|---|---|---|---|
  | `region.m` | 3789 | 4741 | 8530 | 14 |
  | `plq_1p.m` | 691 | 1071 | 453 | 4 |
  | `functionNDomain.m` | 1487 | 1760 | 388 | 8 |
  | `symbolicFunction.m` | 831 | 906 | 85 | — |
  | `plq.m` | 215 | 233 | 84 | — |
  | `plq_1piece.m` | 2603 | 2621 | 28 | — |
  | `domain.m`, `conjugateExpr.m`, `yIntercept.m` | — | — | 0 | — |

* **The five `test*.m` suites are the verification role**: `testSymbolicFunction`, `testcPLQ` and
  `testfunctionNDomain` are byte-identical to the originals, `testMaxMultiRegion` and `testRegion`
  differ by 11 and 39 lines. They run in the ordinary suite and are what pins the merged code
  against the behaviour it was merged from.

So when a defect below is attributed to `plq_1p` or `functionNDomain`, that is **a defect in
CCA2's own merged code**, and CCA2 owns the fix — not a limitation of an outside dependency. The
one consequence worth acting on is §7.1: Step 1 exists TWICE in the root (`convEnvCPLQ.m` and
`plq_1p.convexEnvelope1`) with different coverage, which is the merge not being finished.

---

## 0.1 Reproducibility rule for quoted measurements

**Every number quoted in this file, in `DESIGN.md`, or in a source header must be re-derivable.**
If it came from a randomized sweep, the sweep must be COMMITTED and must set an explicit seed
(`rng(<seed>)`), and the quote must name the script and the seed. A measurement whose generator
was a throwaway script is an unreproducible claim, however carefully it was made.

This rule exists because the toolbox already has several such claims and they have cost real time:

| Quoted | Where | Status |
|---|---|---|
| ~~`QuaPar.eval` wrong at ~1.4% of polyhedral result vertices (5/356)~~ → **replaced 2026-08-02** by **225 of 1205 vertices (18.7%)** under the pre-fix exact test, **0** under the current tolerant one, 0 ring mismatches | `sweepQuaParEvalAtVertices.m`, seed **20260802**, 200 cases; pinned deterministically by `QuaParTest/evalLocatesAPointExactlyAtItsOwnVertex` | **REPRODUCIBLE.** The retired figure's population (maxQuaPar *result* vertices from unrecorded random operand pairs) cannot be reconstructed; the new sweep measures the SAME mechanism on a population that can be — the vertices of randomly generated polyhedral subdivisions, evaluated directly, with no maxQuaPar involved. Comparable in kind, not in value. Reports both the exact and the tolerant test in ONE run, so the fix is a column difference rather than a comparison between source versions |
| ~~`QuaPar.eval` wrong at ~0.8% of CURVED result vertices (9/1105)~~ → **0 of 369 as of 2026-08-02** | `sweepMaxQuaParCurvedSplit.m`, seed **20260802**, 200 cases | **REPRODUCIBLE, and it no longer reproduces.** Every one of the 369 vertices of the 30 assembled curved results evaluates exactly. That is evidence — not proof, the population is regenerated and smaller — that the same `QuaPar.eval` tolerance which closed the polyhedral half closed this one too |
| ~~115 sampled splits, 85 assembled, 340/340 edge midpoints, 5100/5100 interior points~~ → **131 sampled, 30 assembled, 112/112 straight-edge midpoints, 1800/1800 interior points, all exact** → **re-measured 2026-08-13 after the arc-vs-arc work: 142 sampled, 59 assembled, 0/1031 result vertices, 0/571 midpoints and 0/3540 interior points disagree** | `sweepMaxQuaParCurvedSplit.m`, seed **20260802**, 200 cases | **REPRODUCIBLE.** The assembly rate is far lower than the retired 85/115 because this sweep does not pre-filter: of the 131 splits, **80 never reach `maxQuaPar` at all** — Step 2 refuses the triangle with `conjPieceCPLQ:notImplemented` — and 21 more hit `maxQuaPar:notImplemented`. The retired figure evidently counted only splits that got past Step 2. Breaking the guards out by identifier is the point: "85 of 115 assembled" concealed that the dominant obstacle is upstream of `maxQuaPar` |
| 58 → 76 of 395 quadrilaterals after the arc-split work | session handoff, section 4 | **still not reproducible** — generator not kept. It is a BEFORE/AFTER pair across a code change, so re-deriving it needs the sweep run against both revisions; `sweepMaxQuaParCurvedSplit` gives the "after" half on a new population but cannot recover the original 395 |

Note the committed test suite itself IS deterministic — there is no `rand`/`randn`/`randi`
anywhere in it. The problem is entirely that the *validation* sweeps behind these figures were
never made part of the repository. When a future sweep is run, commit it under a name like
`sweepXxx.m`, seed it, and cite it here.

---

## 0. Scope — read this first

**CCA2's goal is `QuaPol` conjugate and biconjugate.** It is *not* to cover every possible
`RatPol` or `QuaPar`.

Everything downstream of a `QuaPol` is a **special case**, and that special-ness is load-bearing:

- the convex envelope of any triangle coming from a `QuaPol` is a very special case of `RatPol`;
- the conjugate of those triangles is a very special case of `QuaPar` — e.g. a parabolic edge only
  ever occurs **surrounded by two parallel rays**.

Consequently **hyperbolic edges never need to be stored**: they never arise from a `QuaPol`
conjugate/biconjugate, nor from any intermediate computation ([COAP]/[JOGO], Karmarkar & Lucet
2026 — see `DESIGN.md`'s reference list).

This is why the table below separates **"unreachable by design"** from **"gap"**. A guard that
fires only on an input the wired pipeline cannot produce is an *assertion protecting an invariant*,
not a missing feature — there is nothing there to implement. Counting such guards as gaps
overstates what is missing, and has caused repeated mis-reporting of this project's status.

Legend:

| Mark | Meaning |
|---|---|
| **OK** | Implemented and covered by tests |
| **GAP** | Genuinely not implemented; a caller can reach it with a legitimate `QuaPol` |
| **N/R** | Not reachable from a `QuaPol` conjugate/biconjugate — out of scope by design, not a gap |
| **INV** | Internal invariant; firing indicates a bug, not an unsupported input |

---

## 1. Conjugate — `conj(f, engine)`

### 1.1 By engine

| Engine | Status | Guard |
|---|---|---|
| `'cplq'` (default) | **OK** — the only implemented engine | — |
| `'pqp'` | **GAP** | `PLQ:conj:engine` — `QuaPol.m:589`, `QuaPar.m:633` |
| `'graph'` | **GAP** | same guard |
| unknown name | **INV** (input validation) | `PLQ:conj:engine` — `QuaPol.m:592`, `QuaPar.m:636` |

Two of the three engines named in the design are **not implemented at all**. This is the single
largest functional gap.

### 1.2 `'cplq'` engine, by input shape (`conjCPLQ.m`)

| Input | Status | Returns | Guard |
|---|---|---|---|
| Full-domain **strictly convex** quadratic (`nv==0, nf==1`) | **OK** | `QuaPol` | — |
| Full-domain quadratic with a **negative** eigenvalue (indefinite or concave) | **GAP**, and it is a REPRESENTATION gap, not a missing construction: conv f = -inf so f* = +inf everywhere, i.e. dom f* is EMPTY, which no mesh encodes (`nf=0` means dim<2, not empty) | — | `PLQ:conjCPLQ:conjugateHasEmptyDomain` |
| Full-domain **affine** (Q = 0) | **GAP**, representation: f* is the indicator of the single point L with value -kappa. The mesh has a needle domain (nv=1, ne=0) but `QuaPar.eval` returns +inf on it | — | `PLQ:conjCPLQ:conjugateIsAPoint` |
| Full-domain **PSD of rank 1** | **GAP**, representation: f* is finite only on the LINE `<s-L,n> = 0`, n spanning null(Q), where `f*(s) = 1/2 (s-L)' pinv(Q) (s-L) - kappa`. QuaPar's dim<2 domain is a SEGMENT or ray between two vertices, not a line | — | `PLQ:conjCPLQ:conjugateIsALine` |
| Single **bounded triangle**, Step 1 envelope has **one** face | **OK** (numeric; fast) | `QuaPar` | — |
| …envelope split with a **rational** face, 2-face split | **OK** (falls back to cPLQ's Step 2/3; slow, ~100 s) | `QuaParCPLQ` | — |
| …envelope split with a **rational** face, 4-face split | **OK** (2026-07-29; symbolic, ~27 min) | `QuaParCPLQ` | — |
| General **bounded** domain (multi-face and/or non-triangular) | **OK** — numeric when every face is convex or the fold succeeds, symbolic otherwise | `QuaPol`/`QuaPar`, else `QuaParCPLQ` | — |
| Any face carrying a **positive definite** quadratic, bounded or not | **OK** (2026-08-24, numeric, 0.01–0.1 s) — `conjConvexPolygon` takes the face WHOLE: no triangulation, any number of facets, and the result is POLYHEDRAL (`QuaPol`), so nothing curved reaches Step 3 from it | `QuaPol` | — |
| **Unbounded** multi-face domain, faces carrying **convex quadratics** | **OK** (2026-08-24, NUMERIC) — the `isDomBounded` gate on the numeric route is gone; a recession direction changes one branch of `conjConvexPolygon` and nothing else | `QuaPar` (meshed) | — |
| **Unbounded** multi-face domain, faces whose convex envelope is **affine** | **OK** (2026-07-31; still symbolic) — the conjugate of an affine piece over an unbounded polygon is a SUPPORT function, so it is `+inf` off a cone, and `maxQuaPar` refuses an operand that is not finite everywhere. `max(0,x,y)` is the canonical case; `TODO.md` G2 | `QuaParCPLQ` | `maxQuaPar:notFullDomain` |
| **Unbounded** face carrying a **non-convex** quadratic, envelope FINITE | **OK** (symbolic; measured 2026-08-25, ~14 s) — RECLASSIFIED. The previous edition read GAP, which was the NUMERIC route's decline read as the whole story. `x*y` on the first quadrant is bounded below on its own recession cone, so conv f is finite and Case C answers exactly: f* is the indicator of the third quadrant. What is missing is a NUMERIC route, i.e. a symbolic FALLBACK, not a gap. Pinned by `conjCPLQTest/unboundedNonConvexFaceWithAFiniteEnvelopeIsANSWERED` | `QuaParCPLQ` | numeric route only: `PLQ:conjCPLQ:notImplemented` |
| **Unbounded** face carrying a **non-convex** quadratic, envelope **-inf** | **GAP**, and a REPRESENTATION one: f decreases without bound along a direction of its own recession cone, so conv f = -inf and f* = +inf everywhere. dom f* is EMPTY, which no mesh encodes. The refusal is correct and names the cause. Pinned by `conjCPLQTest/unboundedNonConvexFaceUNBOUNDEDBelowIsRefusedByName` | — | `convEnvUnbounded:minusInfinity` |
| Single **bounded triangle** whose Step 1 envelope SPLITS, and the fold disagrees with the faces | **OK** (2026-08-25) — refused, not returned. `conjMaxOfSubTriangles` cross-checks the assembled Step 3 against the per-face conjugates it was built from; `f* = max_k (q_k + I_T_k)*` is an identity, so the faces are their own oracle. Measured on `TODO.md` G4's fixture: the fold OVER-estimates by 4x at (-10,0) — which the one-sided edge lower bound cannot see — and the refusal costs 2.5 s against 25+ min for the symbolic fallback. Under `CCA2_CONJ_VERIFY` like the edge bound | — | `PLQ:conjCPLQ:foldDroppedACell` |
| **Unbounded** multi-face domain, a face with a **curved convex** envelope | **GAP** | — | `plq_1p:conjugateFunction:unboundedNonAffine` |
| Step 3 with a **non-triangular** envelope piece | **GAP** | — | `PLQ:conjCPLQ:notImplemented` — `conjCPLQ.m:161` |
| Cubic (`PLC`) input | **N/R** — cubic is for `isConvex` only | — | `assertOperable`; `quaPolToPlq:cubic` |

Those three rows replace a single previously-claimed **OK**. Since Step 1 became recursive
(`convEnvCPLQ`'s `solveTriangleBF`/`splitTwoConvexEdges`, 2026-07-17/18), the envelope of an
indefinite triangle can contain genuinely rational faces, and `conjMaxOfSubTriangles` hands one to
Step 2, which has no rational branch (§1.3).

**Step 2 now falls back to cPLQ's own symbolic Step 2/3 for exactly those envelopes**
(`conjCPLQ.m`'s `conjEnvelopeViaCPLQ`, via the new `ratPolToPlq.m` adapter), which closes the
2-convex-edge case: on `f=xy` over `(2,1),(0,0),(1,0)` the result matches `sup_x ⟨s,x⟩ − f(x)` to
≤ 8.9e-16 at all 10 sampled dual points. Pinned by
`conjCPLQTest.indefiniteTriangleTwoConvexEdgesSplitViaCPLQStep2`.

Note **which half** of cPLQ is reused — its Step 2/3, on *CCA2's* Step 1 output. Running cPLQ end
to end instead does not work, and it is worth recording why, because "cPLQ already does this" is
the natural assumption:

- for a **2-convex-edge** triangle cPLQ's own envelope is the single Appendix A.4 formula, which
  CCA2 established is not always tight (`convEnvCPLQTest.bilinearTwoConvexEdgesSplitIsTight`). Its
  conjugate inherits that: on the triangle above it leaves the paper's own flagged dual point
  `s=(-0.008727,-0.999962)` covered by **no region at all** (evaluates `NaN`).
- for a **3-convex-edge** triangle cPLQ has **no branch whatsoever** — neither `convexEnvelope1`
  nor `conjugateFunction` in `plq_1p.m` has an `nCE==3` case — so it silently returns an empty
  envelope and an empty conjugate.

CCA2's Step 1 is ahead of cPLQ's on both counts (the tightness split, and [COAP] Appendix A.5's
3-convex-edge split, are CCA2's own), so the working combination is CCA2 Step 1 + cPLQ Step 2/3.

There is **no 3-convex-edge case to implement**, in CCA2 or in cPLQ. [COAP] Appendix A.5's split
reduces such a triangle to 2-convex-edge sub-triangles and Step 1 already applies it, so every
face reaching Step 2 has come through it. The edge count describes the *input*; it is not what
fails.

The **4-face** split took longest, and its last failure was in **Step 3** — cPLQ's cross-piece
maximum (`plq.maximumConjugate` → `functionNDomain.maximumP` → `region.maximum`). It is now
**closed** (see the two updates below); this is the history of how, kept because each stage
records a wrong diagnosis worth not repeating.

**Update (2026-07-29) — Step 3 runs, but its ASSEMBLY was still wrong.** Four blockers were
fixed and the 4-face pipeline completed end to end (~15 min). Steps 1 and 2 are **exact** on the
reference case `f = xy` over `conv{(0,0),(3,3),(1,2)}`: the pointwise max of the four per-piece
conjugates matches `sup_{x∈T} ⟨s,x⟩ − xy` at **all 289 points** of a 17×17 dual grid, and Step 1's
envelope reproduces `f*` to grid resolution. The four fixes were

| Site | Defect |
|---|---|
| `region.maxArray` (`region.probeAlong`/`probePerp`) | tie-break probe directions came from constraint *slopes*; a symbolically-zero-but-not-syntactically-zero slope made `d == 0` false, so `-1/d` survived as an unevaluated expression and detonated later in `subsF`'s `simplifyFraction` as `symbolic:kernel:DivisionByZero` |
| `functionNDomain.mergeL` | second accumulation loop passed the **first** region's vertices to `removeTangent`, and none at all (`Unrecognized function or variable 'nP'`) when that first region was empty |
| `region.linear3pt` | translated an edge-local vertex index into a region index **in place**, so the next iteration indexed a 3-element array out of range |
| `plq_1p.conjugateFunction`, `nCE==2` | the `grad` half-plane per edge was the hard-coded pattern `+,−,+` — a statement about one edge *ordering*. For any other ordering two edges got each other's half-plane, so `f*` came out as the **min** of the two edge candidates. Now derived from geometry: the rank-1 Hessian's kernel `k = [b,−2a]` gives `grad = ⟨s−L,k⟩`, so edge *j* owns `{grad ≥ 0}` exactly when its outward normal has `⟨k,n_j⟩ > 0`; an edge parallel to `k` is a level edge and contributes no region |

**Update (2026-07-29, later session) — the ASSEMBLY is now fixed too; this row is OK.** The final
assembled partition went from **57 of 289 wrong to 0 of 289**, exact at every one of the fold's
seven rounds, and the whole fold got *faster* (1645 s vs 1768 s) because a correct partition
carries fewer regions than a damaged one.

Both causes were the same question — "is this linear form ≤ 0 over that polyhedron" — which is an
LP: exact, and it answers unboundedness and infeasibility as first-class results, which matters
because these regions are routinely unbounded. New helpers `region.maxLinear`/`impliedBy`/
`linearForm`/`redundantSubset`/`deleteIfRedundant`/`unionIsExact`; every uncertain answer is a
refusal, since over-describing a region is harmless and under-describing it is the bug.

1. **`region.merge` over-claimed.** With `A = A' ∩ {g≤0}`, `B = B' ∩ {g≥0}`, it returns
   `M = A' ∩ B'`. `M` never *loses* a point (any `x ∈ A' ∩ B'` has `g≤0`, so is in `A`, or `g≥0`,
   so is in `B`) but it can *gain* points belonging to neither. `M = A ∪ B` exactly when
   `A ⊆ B'` and `B ⊆ A'` — equivalently, when `A ∪ B` is convex — and `unionIsExact` decides that
   by LP before any facet is deleted. Two shared facets are refused outright, since a point with
   `g₁≤0, g₂≥0` is in neither operand yet survives `M`.
2. **`region.simplifyUnboundedRegion` dropped non-redundant constraints.** "Delete any constraint
   not passing through a finite vertex" is not a redundancy test; `redundantSubset` uses the real
   one, `max{gᵢ : gⱼ ≤ 0, j≠i} ≤ 0`. Both deletion sites now route through `deleteIfRedundant`, so
   a heuristic can only ever *propose*.

Guarding (1) alone with a constraint-set-equality test was provably sound and made things **worse**
(36 → 125 wrong of 289), because refusing merges left more regions for (2) to damage — which is why
they were fixed together, and why the *exact* condition above is the right one rather than that
stronger proxy: it refuses only what it must.

Fixing these exposed one further latent bug, the usual pattern here: **`region.slopes2`** takes a
curved constraint's tangent at a region vertex lying on it, and a curved constraint need not have
one — previously unreachable because such constraints were being deleted. Skipping the constraint
"works" but costs `maxArray` its probe directions, so it returns undecided more often, `maxEqDom`
falls through to `splitmax3`, and every undecided region splits, compounding round over round —
`maximumP 3` went from 153 s to over 90 minutes. Taking the tangent at the region's finite-vertex
centroid instead runs it in 192 s. **Do not "simplify" that fallback back to skipping.**

`conjCPLQ`'s `assertStep3MatchesPieces` gate stays: it cross-checks the assembled maximum against
the per-piece max — the same `f*` by a different route — and a wrong partition fails silently by
nature, returning plausible numbers rather than erroring. Unit coverage of the two certificates is
in the new `regionTest.m`.

**Correction (2026-07-28).** An earlier edition of this section blamed the numeric/symbolic
boundary and called for exact arithmetic in Step 1. That was **wrong**, and the measurement that
refutes it is worth keeping: at the vertex where a face's denominator vanishes, the numerator
vanishes too, with ∇N parallel to ∇den — residual 0 to 1.3e-16 — so the `0/0` is cleanly
**removable** and the envelope coefficients are accurate to ~1e-16. Nothing needed to be made
exact. What was missing was *resolving* the `0/0` by a limit instead of substituting into it.
`symbolicFunction.subsF` already flagged it (returning `NaN`) and `region.funcVertices` already
used a `NaN`-then-limit idiom, but `symbolicFunction.limit` takes **iterated univariate** limits,
which return `NaN` again for a bivariate `0/0`. The new `symbolicFunction.limitDirectional`
restricts to a line through the point — one univariate limit, which a rational function always
answers — and requires several directions to agree, which is itself the check that the limit
exists. Pinned by `conjCPLQTest.indefiniteTriangleThreeConvexEdgesUsesStep3`. `biconj` is
unaffected throughout — see §3.

### 1.3 Per-piece conjugate (`conjPieceCPLQ.m`, Step 2)

| Piece | Status | Guard |
|---|---|---|
| Affine | **OK** | — |
| Strictly convex (PD) quadratic | **OK** | — |
| Convex rank-1 PSD quadratic (COAP A.4 envelope output) | **OK** | — |
| Indefinite quadratic, **0 or 1** convex edge (any frame/shift) | **OK** | — |
| Indefinite quadratic, **2** convex edges | **N/R** | `conjPieceCPLQ:notImplemented` — `conjPieceCPLQ.m:586` |
| Concave (negative semidefinite) | **N/R** — call `convEnvCPLQ` first | `conjPieceCPLQ.m:128` |
| Genuinely **rational** piece (nonzero `RatPol` denominator) | **GAP** (was listed N/R — see below) | `conjPieceCPLQ.m:107` |
| Non-triangle / unbounded piece | **N/R** — Step 1 only emits bounded triangles | `conjPieceCPLQ.m:103` |
| Cubic numerator | **N/R** | `conjPieceCPLQ.m:112` |

Notes on these rows:

- **2 convex edges** (**N/R**): Step 1 (`convEnvCPLQ`) always convexifies such a piece into a
  rank-1 PSD quadratic before Step 2 sees it. Reaching this branch requires hand-feeding a raw
  indefinite piece, bypassing Step 1. See
  `conjPieceCPLQTest.bilinearTwoConvexEdgesIsUnreachableFromTheWiredPipeline`.
- **Rational piece** (**GAP**, reclassified 2026-07-28): the previous edition claimed
  `conjSingleTriangle` conjugates the *original* quadratic piece directly and so never
  materializes a rational envelope. That reasoning covers only the branch where Step 2 accepts
  the raw piece. When it does not (concave, or ≥2 convex edges), `conjSingleTriangle` falls back
  to Step 1 — and since Step 1 became recursive, that envelope can carry rational faces, which
  `conjMaxOfSubTriangles` then feeds straight into Step 2. A legitimate single-triangle `QuaPol`
  reaches this guard, so it is a gap in `conjPieceCPLQ`, not an unreachable assertion.
  **It is no longer a gap in the pipeline**: `conjCPLQ` catches it and routes the whole envelope
  through cPLQ's Step 2/3 instead (§1.2). Implementing a native numeric rational branch here would
  buy speed — ~100 s drops to milliseconds — not new coverage.

---

## 2. Partial conjugate — `partialConj(f, idx, engine)`

| Engine | Status | Guard |
|---|---|---|
| `'cplq'` | **GAP** | `PLQ:conjCPLQ:partialNotImplemented` — `conjCPLQ.m:53` |
| `'pqp'` | **GAP** | `PLQ:partialConj:engine` — `QuaPol.m:607`, `QuaPar.m:651` |
| `RatPol` input | **GAP** | `RatPol:partialConj:notImplemented` — `RatPol.m:631` |

**`partialConj` is not implemented for any engine or any type.** The dispatch exists; every path
errors.

---

## 3. Biconjugate and convex envelope

`biconj` dispatches by input shape in `biconjCPLQ.m`, exactly as `conjCPLQ` does. Cases A and C
are the literal `conj∘conj`; **Case B is not**, and that is what closed this section's former
blocker. For a compact domain `T` and `q` continuous on it, `f = q + I_T` gives

    f** = cl conv f = conv f = conv(q + I_T)

which is precisely what **Step 1** (`convEnvCPLQ`) computes. So the biconjugate of a bounded
triangle is available in closed form with no second conjugation, and therefore without needing the
unbounded multi-face conjugate of §1.2 — which stays open, and is still the right long-term fix
for `conj` itself.

| Input | `conj` | `biconj` | Route |
|---|---|---|---|
| **A** full-domain strictly convex quadratic | **OK** → `QuaPol` | **OK** → `QuaPol` | `conj∘conj` (conjugate is again a full-domain quadratic, so Case A twice) |
| **B** single bounded triangle | **OK** → `QuaPar` (except §1.2's rational-sub-piece row) | **OK** → `RatPol` | Step 1's convex envelope, `convEnvCPLQ` |
| **C** general bounded multi-face domain | **OK** → `QuaParCPLQ` | **OK** → `QuaParCPLQ` | `conj∘conj`; `QuaParCPLQ.conj` routes through cPLQ's own symbolic machinery |
| **unbounded** domain | **GAP** | **GAP** | falls through to `conj∘conj`, which errors at `conjCPLQ.m:103` |

Case B's route is **strictly stronger** than `conj∘conj` would have been, not merely equal: it
succeeds on the very triangles whose *conjugate* Step 2 cannot compute (§1.2's rational-sub-piece
row) — Step 1 produces the envelope regardless, so `f**` is available even where `f*` is not.

Correctness is checked against a ground truth that uses **none** of the conjugate pipeline:
`f*(s)` by exact maximization of `⟨s,x⟩ − q(x)` over the triangle, then `f**(x) = sup_s ⟨s,x⟩ −
f*(s)` numerically. Agreement ≤ 8.9e-16 at interior points of 7 triangles covering every Step 1
branch (affine, convex, concave, indefinite with 0/1, 2 and 3 convex edges — i.e. 1-, 2- and
4-face envelopes). Regression tests: `biconjCPLQTest` (10 tests) and
`conjCPLQTest.biconjCoverageByInputCase`.

| Operation | Status | Notes |
|---|---|---|
| `QuaPol.biconj` | **OK** for every bounded domain | see the table above; unbounded domains still error |
| `QuaPar.biconj` | **OK** for a polyhedral bounded triangle | a triangle with a genuinely **parabolic** side is not a `convEnvCPLQ` input, so it falls through to `conj∘conj` and still errors |
| `RatPol.biconj` | **GAP** | `RatPol:biconj:notImplemented` — `RatPol.m:636` |
| `RatPol.conj` | **GAP** | `RatPol:conj:notImplemented` — `RatPol.m:625` |
| `convEnvCPLQ`, nonconvex quadratic over a **bounded triangle** | **OK** | 0/1/2/3-convex-edge splits all implemented |
| `convEnvCPLQ`, nonconvex quadratic over a non-triangle/unbounded domain | **GAP** | `convEnvCPLQ:notImplemented` — `convEnvCPLQ.m:65` |
| `convEnvCPLQ`, multi-piece with an **unbounded** face | **GAP** | `convEnvCPLQ:notImplemented` — `convEnvCPLQ.m:457` |
| `convEnvCPLQ` of a nonconvex **full-domain** quadratic | **N/R** — envelope is not finite | `convEnvCPLQ.m:60` |

---

## 4. Pointwise max of conjugates — `maxQuaPar` (Step 3)

Line numbers below were re-derived from the source on 2026-08-14; the ones this table carried
before that date predate the arc-vs-arc work and were off by ~1400 lines.

| Case | Status | Guard |
|---|---|---|
| Both inputs purely polyhedral | **OK** | — |
| Exactly **one** input with parabolic edges | **OK** (added 2026-07-27) | — |
| **Both** inputs curved | **OK** (2026-08-13) — arc-vs-arc clipping, the two-arc splits and the far-field over-extension are all done; see §4.1 | — |
| Splitting a cell that **already carries an arc** | **OK** (2026-07-30) | — |
| Splitting curve genuinely **crosses** a cell's arc | **OK** (2026-08-13) — the crossing is a third boundary hit, on the arc's own edge, and the existing two-hit split then divides the arc with everything else | — |
| Splitting an **unbounded** cell that carries an arc | **OK** (2026-08-13) — the unbounded branch restores the inherited arc exactly as the bounded one does | — |
| …and the resulting unbounded half carries **both** the inherited arc and the splitting curve | **OK** (2026-08-14) — a chord cannot close an unbounded piece, but a RAY from the vertex between the two arcs can, and the split is used only when each half passes all three exact invariants; otherwise it still refuses | `maxQuaPar:notImplemented` — `maxQuaPar.m:2864` (now unreached on every fixture) |
| An **unbounded** piece whose vertices all lie ON `{f1=f2}` | **OK** (2026-08-14) — `splitCell`'s unbounded "rest" piece can have exactly the two crossing points as its vertices, and then neither they nor their centroid can decide the winner; it was coming out of floating-point noise. Read in the piece's **recession cone** instead, sharing `assignSideFromCone`'s probe | — |
| An unbounded piece that genuinely **straddles** `{f1=f2}` | **INV**, backstop — the assigned winner must still hold at infinity along each of the piece's own rays, tested exactly via the leading coefficient. Fires on no current fixture; kept because a genuine straddle is otherwise silent | `maxQuaPar:notImplemented` — `maxQuaPar.m:3472` |
| Split curve meets an unbounded cell once and escapes as a **parabola** | **GAP** (defensive) — an unbounded curved edge | `maxQuaPar.m:2650`, `maxQuaPar.m:2904` |
| Clip line cutting one arc **twice** (arc bulging across) | **OK** (2026-08-24) — `clipPolyHalfPlane` cuts the cell along the line through the arc's own stationary point PARALLEL TO THE PARABOLA'S AXIS, which meets the conic exactly once, so each half has one sub-arc and the clip line crosses each at most once. Was a GAP claimed to need "curveAfter becoming a set"; it does not. See §4.4 | — |
| Curved cut crossing a cell's own arc twice | **GAP** (defensive) | `maxQuaPar.m:1361` |
| Curved cut that **separates** an unbounded cell | **N/R** — each component would need the cutting conic running to infinity | `maxQuaPar.m:1494` |
| Curved corner cut leaving an unbounded survivor with two arcs / one ray | **GAP** (defensive) | `maxQuaPar.m:1668`, `maxQuaPar.m:1678` |
| Two-arc lens cut `A→M→B` that leaves the cell | **GAP** — the two arcs join corners on opposite parabola branches, so the polyline exits; needs a different subdivision | `maxQuaPar.m:1751` |
| Face vertex inside the **open interior of an arc** (arc must be split) | **OK** (2026-08-24) — the same axis-parallel cut, applied at that vertex. NO LONGER "0 occurrences observed": it fires on every indefinite quadratic over a polygon, which is how it was found | — |
| Curved **ray** edge | **N/R** — `QuaPar` has no unbounded curved edge | `maxQuaPar.m:839` |
| Face with more than one curved edge | **N/R** — one `Ec` slot per piece | `maxQuaPar.m:903` |
| Input not finite everywhere | **N/R** — Step 3 combines full-domain conjugates | `maxQuaPar:notFullDomain` — `maxQuaPar.m:848` |
| Difference of the two candidate quadratics is an irreducible ellipse/hyperbola | **N/R** | `maxQuaPar:notDegenerate` — `maxQuaPar.m:2510` |
| Clipping a cell whose arc bulges, reached BELOW `clipPolyHalfPlane` | **INV** (2026-08-24) — the caller subdivides first, so `clipPolyHalfPlaneCurved` is only ever entered on a cell where the at-most-one-crossing invariant holds. Kept because it is the invariant every branch of that routine assumes | `maxQuaPar:internal` |
| Split curve meeting an unbounded cell only at the SINGULAR POINT of a line pair | **OK** (2026-08-24) -- the gradient vanishes there, so the branch direction comes from the NULL DIRECTIONS of the quadratic part instead of from the tangent. Was a SILENT WRONG ANSWER: the caller's tangency branch then read the winner at the cell's centroid, which for a cone IS that point. See section 4.5 | -- |
| Assembly after an arc split: one side CURVED, the other STRAIGHT | **GAP**, and it is the live one — `matchHalfEdges` pairs curved with curved, and the subdivision is currently consistent per face PAIR rather than globally. This is what still sends an indefinite quadratic over a POLYGON to the symbolic fallback. `DECISIONS.md` 2026-08-24 has the design | `maxQuaPar:internal` — `assemblePieces` |
| `maxQuaPar:internal` (other sites) | **INV** | assembly-topology invariants |
| `maxQuaPar:pieceInvariant` (`maxQuaPar.m:335`) | **INV**, opt-in | the three exact piece invariants; off unless `MAXQP_ASSERT` is set — see §4.2 |

### 4.5 The split direction at a line pair's SINGULAR POINT (2026-08-24)

The smallest wrong answer this file has produced, and the shape is not exotic.

    h1 = max(s1,0)^2/2      two cells, split by the s2-axis
    h2 = max(s2,0)^2/2      two cells, split by the s1-axis

`max(h1,h2)` must split the FIRST QUADRANT along `s2 = s1`, giving 5 cells. It gave **4**, with the
whole quadrant carrying `h1`, so `max(h1,h2)(2,3)` came back `2.0` against a truth of `4.5`.

`{h1 = h2}` is `s1^2 = s2^2` -- a degenerate conic, and a PAIR OF LINES crossing at the origin,
which is the quadrant cone's apex and its only vertex. `decideWinner` correctly refused (the two
rays give opposite asymptotic signs). `splitCell` found ONE boundary crossing, at the apex, and
`splitUnboundedAtOneCrossing` took the branch direction as the perpendicular to the GRADIENT --
which vanishes at a line pair's crossing point. It returned "no branch here", and the tangency
branch below it resolved the winner at the cell's CENTROID, which for a cone is that same apex,
lying exactly on the curve.

Where the gradient vanishes, `diffRow(X + t d) = (t^2/2) d'Qd`, so the branches through `X` are the
NULL DIRECTIONS of `Q` -- available in closed form from its eigen-decomposition
(`d = sqrt(-l2) u1 +- sqrt(l1) u2`), with no root-finding and no case analysis. Both are put
through the existing recession test.

**A cone whose apex is that singular point is the dual of any WEDGE face**, so every 4-cone
conjugate hit this. It is why unbounded input had to be gated away from the numeric route at all;
the 4-cone fan now assembles to 8 cells with error 0 against the definition sup.

Pinned by `maxQuaParTest/twoHalfPlaneQuadraticsSplitTheirSharedQUADRANT` (asserting the cell COUNT
as well as the values, because 4-vs-5 is the defect) and `aFourConeFanFoldsExactlyAgainstItsOwnPieces`.

### 4.4 The axis-parallel split — why the two arc guards above became OK (2026-08-24)

Both guards refused the same thing for the same stated reason: the result would be "bounded by two
separate arcs, neither of which is representable as one QuaPar face", and lifting it "means
generalizing a piece to carry several arcs (curveAfter becoming a set)". **It does not.** One arc
per face is an invariant maintained by SUBDIVIDING, exactly as `splitCell` already does, and there
is a canonical cut:

> Cut the cell along the line through the split point **parallel to the parabola's axis**. Such a
> line meets the conic EXACTLY ONCE -- which is precisely what `parabolaArcFrame`'s own header says
> makes `u` a global monotone parameter -- so it divides the arc into two sub-arcs, one per half,
> and creates no second arc. Both halves carry the same face function, so the FUNCTION is unchanged
> and only the subdivision is finer.

For the bulge case the cut point is forced rather than chosen: along the arc `nrm*x'-c` is the
quadratic `A2 u^2 + A1 u + A0` (`lineCoeffs`), so it has ONE stationary point `u* = -A1/(2 A2)` --
the quantity `arcBulgesAcross` already computes to DETECT the case -- and each sub-arc lies on one
side of it, hence is monotone. The recursion is one level deep by construction: clipping a half by
the split line itself restricts to `u - u*`, which is affine.

`maxQuaParTest` 29/0 and `addQuaParTest` 4/0 with the split in. What it does NOT yet close is the
assembly row above: the split is applied per face PAIR, so two pieces sharing an arc can disagree
about where it is divided.

### 4.1 Arc-vs-arc, measured 2026-08-13

Both operands curved used to be refused outright, then produced answers that were correct near the
arcs and wrong far from them. What that turned out to be, in the end, is one sentence: **a curved
edge is a bounded ARC and its conic is not**, so "every bounding conic, sign-oriented, ≤ 0" — the
point-location rule — is exact on a parabola's convex side and admits points arbitrarily far away on
its concave side. `QuaPar.chordCuts` derives the missing constraint (the arc's chord) per face.

Four other defects were fixed along the way, each of which had produced a silently wrong answer:
a single boundary crossing on an unbounded cell read as a tangency when the curve escapes through
the recession cone; a conic cut skipped when its only crossings are at the cell's corners (the
generic case, since conjugates of triangles sharing an edge have arcs through the same two dual
points); half-edges and boundary walks identifying an edge by its endpoints alone, which four arcs
between the same two points make ambiguous; and a cut polyline leaving one half REFLEX, which no
half-plane test can represent.

Measured now, on the two committed sweeps:

| sweep | before | now |
|---|---|---|
| `sweepMaxQuaParCurvedSplit(20260802, 200)` | 131 sampled, 30 assembled | **142 sampled, 59 assembled, 0 of 1031 vertices / 571 midpoints / 3540 interior points wrong** |
| 397 seeded quadrilateral splits, 200 directions at radii 1–500 | 7 of 64 arc-vs-arc results wrong in the far field | **0 of 64** |

**`maxQuaParTest` is 28 pass / 0 fail as of 2026-08-14** (it was 25 / 1). The red that closed was
`arcVsArcRefusesAnUnboundedTwoArcSplit`, and closing it took four defects, of which the first two
are the reason two earlier attempts at it failed:

1. **`pieceRecessionRays` got the parabola's axis from an exact discriminant.** `arcNullDirs` solves
   `d·Q·d' = 0` exactly and returns **nothing** when `b²−4ac` comes out negative — which is what a
   floating-point parabola's `Q` does about half the time, being only semidefinite up to rounding
   (measured `−2.78e-17`). The derived chord was then never emitted, the piece's constraint region
   stayed a slab open at both ends, and `reccConeViolation` reported it receding in both senses.
   `parabolaArcFrame` has always taken the smallest-magnitude eigenvector instead.
2. **`curveAfter ≠ 0` does not mean "this edge is curved".** `boundedPiece` tags every piece with
   the closing edge's index, including the straight-splitting-curve case where `curveEc` is all
   zeros — `pieceIsCurved`'s header says so, and says why the tag must stay. Five call sites read
   the tag as "is curved" anyway: `polyConstraints` emitted **no half-plane at all** for an
   ordinary straight edge; `pieceStraightEdges` skipped it, so every boundary minimisation built on
   that list was blind to it; and `containmentViolation`/`boundaryMinOf` called `parabolaArcFrame`
   on an all-zero conic, which raises `degenerateAxis` — the crash that made `MAXQP_ASSERT`
   unusable on three of the four arc-vs-arc fixtures.
3. **A whole unbounded piece could have its winner decided by floating-point noise.**
   `splitCell`'s unbounded "rest" piece can come out with exactly the two crossing points as its
   vertices — both, by construction, ON `{f1=f2}` — so `assignSide` had nothing to read the winner
   from, and its centroid fallback is on that line too. Read in the piece's **recession cone**
   instead, sharing the probe `assignSideFromCone` has used since it was written for the same
   problem in `splitUnboundedAtOneCrossing`.
   Worth recording how this was nearly mis-diagnosed: the symptom (an unbounded piece carrying the
   wrong operand, beaten by `+Inf` along its own ray) reads exactly like `{f1=f2}` being a *pair*
   of parallel lines with the second escaping through the recession cone — a real subdivision gap.
   Measuring the conic refuted it: `diffRow` there is `[0 … 0 −1.4979 −3.6486 5.4652]`, its entire
   quadratic part zero, so `{f1=f2}` is a **single straight line** and nothing straddles it.
4. The ray cut itself, now gated on all three exact invariants per half rather than on heuristics.

The six heuristics that did **not** separate the good case from the bad, and which should not be
re-derived, are in `DECISIONS.md`.

Accuracy of the curved path, measured on a randomized sweep of convex quadrilaterals split by a
diagonal (115 splits, 85 assembled): exact to ~1e-14 at **all 340** straight-edge midpoints and
**all 5100** interior sample points, with **zero** `maxQuaPar:internal` crashes and **zero**
arrangement-validity violations. See `maxQuaPar.m`'s header VALIDATION block.

### 4.2 Verification tools for the curved path — all opt-in, all exact

None of these is on in a production run; each exists because the far-field defect was a *silently*
wrong answer, and the lesson recorded in `FARFIELD_FIX_PLAN.md` is that sampling is a development
tripwire and never an acceptance test.

| tool | switch | what it decides | how |
|---|---|---|---|
| `assertPiecesWellFormed` (`maxQuaPar.m:286`) | global `MAXQP_ASSERT` (`1` = containment + winner domination, `2` = also the recession cone) | three exact invariants on every piece *before* assembly: it lies inside both source faces; its carried operand really dominates on it; its encoded region recedes exactly where it declares rays | closed-form minimisation and `pieceRecessionRays`; nothing sampled |
| `QuaPar.eval` validate mode (`QuaPar.m:171`) | global `QUAPAR_VALIDATE` | two faces admitting one point with **different** values — a shared boundary is fine, disagreement on it is an over-extended face | raises `QuaPar:eval:ambiguous` instead of letting last-admitter-wins resolve it |
| `verifyMaxIsExactSymbolically` | called explicitly | `g == max(g1,g2)` over whole **regions**: identity, domination and attainment on every `R_k ∩ F_i ∩ G_j` | minimises each quadratic difference in closed form (quadratics along segments and rays, quartics in a parabola's own frame) |
| `verifyFacesCoverThePlane` | called explicitly | the half `verifyMaxIsExactSymbolically` could not reach: that the faces leave **no hole** | see §4.3 |

`MAXQP_ASSERT=2` costs seconds per call; `1` is cheap enough to leave on in a test suite.

### 4.3 The covering proof — the last part of Phase 4 that was still sampling

`verifyMaxIsExactSymbolically` proves `g = max(g1,g2)` **on every face**. It cannot see a region
belonging to **no** face at all, because every one of its checks quantifies over a face; until
2026-08-14 the only evidence for coverage was `partitionReport`, which samples. The hole fixed on
2026-08-13 was one point wide at that sampling density and was found by accident, which is the
argument for replacing it.

`verifyFacesCoverThePlane(g)` decides coverage from four checks on the constraint data, each of
them a statement about a whole curve rather than about probe points:

| | check | decided by |
|---|---|---|
| **A** | every edge separates **two** faces — a full-domain result has no domain boundary | `F(j,:)` has no zero |
| **B** | every edge lies inside **both** its faces' constraint regions | maximising each constraint along the edge: a quadratic in the parameter on a segment or ray, a **quartic** in the parabola's own frame on an arc |
| **C** | no face's constraint region has boundary anywhere **other than its own edges** — `{c = 0} ∩ R_k` sits inside the edges lying on that curve | splitting the whole curve at the roots of every other constraint; between two consecutive roots each constraint has constant sign, so one evaluation settles a whole interval exactly |
| **D** | no face is squeezed onto a curve by two constraints on it with opposite orientations | comparing the oriented conics whose restriction to the edge vanishes identically |

Together these force `∂A = ∅` for `A = ⋃ R_k`, and since `A` is closed, nonempty and `ℝ²` is
connected, `A = ℝ²`. The argument is written out in the file's header.

Check **C** is the one that matters most: it is exactly the far-field over-extension of §4.1,
stated as a property rather than as a symptom. Before `QuaPar.chordCuts`, a face on a parabola's
concave side had `{arc conic = 0} ∩ R_k` running to infinity along both branches while its only
edge on that curve was a bounded arc — which C reports by naming the curve and the parameter range,
instead of waiting for a probe to land in the over-extended region.

**The prover can be made to fail, and that is pinned.** A check that accepts every correct result
is equally consistent with one that examines nothing, and that is not hypothetical here — an
independent review found three routes by which this one could have returned "no hole" without
looking at a constraint (an unparametrisable edge, an unparametrisable curve, and a residual-based
"this constraint is the curve I am walking" test whose shared scale one ill-conditioned sibling
could inflate until every other constraint was discarded). All three are closed, and
`maxQuaParTest/coverProofRejectsBrokenArrangements` breaks a certified result three ways — a
same-sign edge, an edge with one face, and a face whose edges no longer reach as far as its
constraints — and requires a finding each time.

The same floating-point caveat as §4.2 applies: the structure is exhaustive, the comparisons carry
relative tolerances.

---

## 5. Arithmetic and derived operators

| Operation | `QuaPol` | `QuaPar` | `RatPol` |
|---|---|---|---|
| `add` / `sub` | **OK** | **OK** (partial, see below) | **GAP** — not implemented |
| `scalarMul`, `negate` | **OK** | **OK** | **OK** |
| `addQuadratic` | **OK** | **OK** | **OK** |
| `eval` | **OK** | **OK** (see §7) | **OK** |
| `isConvex` | untested (`%TO BE TESTED` in source) | untested | untested |
| `plotDomain` | **OK** | linear edges only — `QuaPar:plotDomain:curved` (`QuaPar.m:542`) | **OK** |

`addQuaPar` carries **11 distinct `addQuaPar:notImplemented` guards** (curved rays, >1 curved edge
per face, degenerate conics, disconnected clip results, new unbounded curved rays, …). These are
**GAP**s for arbitrary `QuaPar` input, but most are **N/R** for `QuaPar`s that arise from a
`QuaPol` conjugate. They have not been individually classified — doing so requires the same
"can a QuaPol conjugate reach this?" analysis applied per guard, and is listed as follow-up work.

Derived operators — all thin compositions, all **OK** within the limits of the `conj` they call:

| Operator | Built on | Validity precondition |
|---|---|---|
| `infConv(f,g)` | `conj((conj f)+(conj g))` | `f,g` **convex** |
| `moreau(f,mu)` | a **single** `conj` (expand-the-square, [HIRIART-URRUTY-07]) | any `f`, convex or not |
| `lasryLions(f,lam,mu)` | `−moreau(−moreau(f,λ),μ)` | as `moreau` |
| `proxAverage(f,g,lam,mu)` | two `conj` around a weighted `add` | `f,g` **convex** |

---

## 6. Class hierarchy

| Item | Status |
|---|---|
| `QuaPol`, `QuaPar`, `RatPol`, `QuaParCPLQ` | **OK** |
| Common `RatPar` parent (`DESIGN.md` II.3) | **OK** — built 2026-07-27 |
| Axis markers `Rat`/`Qua` × `Par`/`Pol` | **OK** — property-less abstract traits, so either axis can be queried alone |
| `RatPar` itself as a *result* type (rational-on-parabolic) | **N/R** — abstract; no operator produces one |

Every type now inherits from `RatPar`, so `conj`/`biconj`/`convEnv` are statically typed: they
return a `RatPar`, and `kind()` reports the concrete one (`'QuaPol'` / `'RatPol'` / `'QuaPar'` /
`'QuaParCPLQ'`). Use `isMeshed()` to tell a `V/E/Ec/F/P` mesh from the still-symbolic
`QuaParCPLQ`. See `RETURN_TYPE.md` and `RatPar.m`; the contract is pinned entry by entry in
`RatParTest`.

---

## 7. Known defects (not "unimplemented" — these are wrong answers)

**The arc-vs-arc far-field defect that headed this section from 2026-08-04 is CLOSED (2026-08-13).**
It was not a family of edge cases but one sentence — *a curved edge is a bounded ARC and its conic
is not* — and `QuaPar.chordCuts` supplies the missing constraint per face. §4.1 has the mechanism
and the re-measured sweeps; §4.2 lists the tools built to keep it closed, and §4.3 the covering
proof that replaces the sampling that let it hide.

**Three more silent wrong answers were found and fixed on 2026-08-14**, all from one overloaded
flag — `curveAfter ≠ 0` read as "this edge is curved", when `boundedPiece` also sets it for a
STRAIGHT splitting curve (§4.1). `polyConstraints` emitted no half-plane at all for such an edge, so
a piece could be admitted two units outside itself; `pieceStraightEdges` skipped it, so the piece
invariants were blind along it; and two more sites called `parabolaArcFrame` on an all-zero conic,
which is the `degenerateAxis` crash that had made `MAXQP_ASSERT` unusable on three of the four
arc-vs-arc fixtures. The only `maxQuaPar` case still open is a **refusal**, not a wrong answer — an
unbounded piece straddling `{f1=f2}` when that curve is a pair of lines (§4, §8.6).

| Defect | Impact | Where |
|---|---|---|
| ~~**A region with no finite VERTEX was declared EMPTY**~~ — **FIXED 2026-08-16** | `simplifyUnboundedRegion` counted finite vertices and returned `region.empty()` when there were none. A half-plane has none; nor does a slab, nor the whole plane. **A half-plane is exactly what a TANGENT vertex produces**, so this was not a corner case: the cone at a vertex is built from the two edges meeting there, and when those are tangent — an arc and its chord touching, which is how a curvilinear piece ends — both half-planes are the SAME one. Measured on piece 9 of `f*` for `x·y` over the parallelogram `conv{(0,0),(2,0),(2.5,1),(0.5,1)}`: the cone at `(1/4,7/8)` is `{2x/3 + y ≥ 4/3}` carrying `x/4 + 7y/8 − 1/2`, and deleting it left that piece's conjugate uncovered on exactly that half-plane — 6 of 10 probe points wrong or uncovered against a brute-force sup. The verdict is now refuted by a WITNESS: a feasible point certifies non-emptiness, failing to find one proves nothing, so the old verdict stands whenever no witness turns up and this can only ever recover a region that was being deleted. Pinned by `regionTest/aHalfPlaneIsNotEmpty` and its converse | `region.simplifyUnboundedRegion` |
| ~~**A BOUNDED region with an extra constraint was read as UNBOUNDED**~~ — **FIXED 2026-08-16**, the same root cause as bug 1 without the lens collision | Both slot conventions in `conjugateOfPiecePoly` are chosen by `size(ineqs,2) == nv`, so a bounded region carrying a constraint that bounds no edge is read with the unbounded one: `endNv = nv−1`, one edge cell SHORT, and every edge read one slot along. A bounded piece's conjugate is finite everywhere, so the missing cell is a hole where the answer must be finite. Measured on the same piece 9 — 3 vertices, 3 genuine edges, plus a conic touching one vertex, so 4 constraints for 3 vertices. `functionNDomain.edgeIndexList` now supplies the list for any bounded piece the count mislabels, not only for a lens; the count-matching branches keep precedence, so a region the conventions already describe correctly is untouched. Together with the row above: **6 of 10 probe points wrong or uncovered → 2**, and the parallelogram's `biconj` went from `emptyResult` to a finite answer, exact at all four vertices and `+inf` outside. Pinned by `functionNDomainTest/aBoundedPieceWithATangentVertexConjugatesOntoTheWholePlane` | `functionNDomain.conjugateOfPiecePoly` |
| **`getInterior` on a SINGULAR quadratic** — LIVE, found 2026-08-16 | `getInterior` separates an edge cell from its neighbours by eliminating `s` between `x = ∂₁f` and `y = ∂₂f`. When `f` is a singular convex quadratic the gradient map is not invertible and that elimination returns the map's IMAGE LINE, which is not a separating curve — so two edge cells come out on the SAME region and the first one checked wins. On piece 9 of the parallelogram's `f*` (where `f` is constant along a whole edge and `∇f = 0` at two of three vertices) the chord's cell wins over the arc's and is low: `0.0176` and `−0.0138` where the truth is `0.0333` and `0.0363`. End to end that is the last 2 of 10 interior points of the parallelogram's `f**`, each about 4% LOW. **Do not attack the `isQuad` chord rewrite for this** — chording the vertices the conic actually touches makes it WORSE (2 wrong of 10 → 3) and skipping the rewrite changes nothing; both measured 2026-08-16, in `DECISIONS.md`. Reproduce with `functionNDomainTest.parallelogramPiece9` in about a minute | `functionNDomain.getInterior` |
| ~~Cross-piece maximum decided dominance from ONE probe point~~ — **FIXED 2026-08-02** | `region.maxArray` falls back to interior probes when the two candidate functions tie at every vertex of the cell, and it returned on the FIRST feasible probe — asserting global dominance from one sample. Wrong on any cell the tie line crosses. Concretely, for `f = x*y` on `[0,1]²` given as TWO triangles, the two conjugates overlap in the lens `{(s1+s2)² ≤ 4s1, (s1+s2)² ≤ 4s2}`, whose only vertices `(0,0)` and `(1,1)` both lie on `s1 = s2`; the lens came back carrying `s2` throughout, so `f*(0.66,0.18)` was `0.18` for a truth of `0.66` — **wrong at 800 of 40000 grid points, worst error 0.48**, and silent. The recorded check that had passed used 7 hand-picked dual points and missed it. Now `region.maxFromPts` requires every decisive probe to agree and reports disagreement, which makes the caller split on `f1 = f2`. Pinned by `biconjugateTest/twoFaceConjugateIsExactIncludingTheLensCell` | `region.maxArray` / `region.maxFromPts` |
| ~~`getNormalConeVertexQ` indexed a vertex outside its own guard~~ — **FIXED 2026-08-02** | Both halves of the routine pair vertex `j` with a neighbour `k` (`j∓1`) and guard `k` falling off the ends — then the `isZero` fallback that follows re-probed from vertex `k` **unguarded**, raising `MATLAB:badsubscript` on `obj.vx(0)` / `obj.vx(nv+1)`. Same shape as the already-corrected `getEdges` / `splitmax3` / `poly2orderUnbounded` cases. `region.probeVertexIndex` now clamps to `j`. (The same block also had the `py = py(1)` before `isempty(py)` inversion, corrected twice before elsewhere in the same function) | `region.getNormalConeVertexQ` |
| ~~**Conjugate of a BOUNDED piece with a parabolic edge**~~ — **FIXED 2026-08-15.** The last step was the one the four earlier attempts had all been circling: both slot conventions identify an edge by a VERTEX INDEX (`edge j is at ineqs(j)` / `at ineqs(j+1)`), and a LENS's two edges join the SAME pair of vertices — so neither convention can name them apart, and no reassignment of slots to constraints ever could. `functionNDomain.edgeIndexList` derives the edge list from the geometry instead (both endpoints on the constraint AND its own curve between them inside the region — the second half is necessary, since the lens's redundant conic passes through both vertices too), and the three edge-indexed readers take it as an argument. It turned out much smaller than the mapped-out refactor predicted: `getNormalConeEdgeQ` and `getNormalConeEdgeQ3` are the SAME routine under two slot conventions and collapse into one (`getNormalConeEdgeQE` over `region.coneNormalAt`), `getSubdiffVertexT2` and `T2Q` are identical on these inputs, and both loops live in `conjugateOfPiecePoly` itself so the list needs no field and no signature change. Entered only on the lens signature (`edgesStillCollide`); no constraint is dropped, so the two unsound drops recorded in `DECISIONS.md` stay ruled out. **Measured: both half-lenses conjugate to 3 cells — 2 vertex cones plus the arc, the chord's cell being a ray that drops out — exact against a brute-force sup at 12 points, 0 wrong, where the old code was `+inf` at `(0,0)` and `(-1,0.5)` and `0` at `(2,-1)` for a truth of `1/2`. `biconjugateTest` is 7 / 0**, including `biconjugateOverATwoFaceSubdivisionIsTheEnvelope`, which had been the repository's only red. Pinned at unit level by `functionNDomainTest/twoEdgesOnOneVertexPairAreBothKept` and `halfLensConjugateIsFiniteEverywhereAndExact` (~10 s, against 10–40 minutes through the pipeline) | `functionNDomain.edgeIndexList` / `conjugateOfPiecePoly`, `region.getNormalConeVertexQ` / `getNormalConeEdgeQE` |
| _history of the row above_ — cause CORRECTED and largely FIXED 2026-08-15 | The half-lens `{(s1+s2)² ≤ 4s1, (s1+s2)² ≤ 4s2, s2 ≤ s1}` carrying `s1` is bounded with 2 vertices and 2 edges (one arc, one segment), and its conjugate used to come back as the conjugate of the CHORD — finite only on a strip, which collapses the whole biconjugate because `f**` is a MAX and its domain is the INTERSECTION of the per-piece conjugate domains. **The cause recorded here until 2026-08-15 was wrong**: it was not the `size(d.ineqs,2) == d.nv` count. `getEdgeNosInf` numbers an edge by one of its endpoint VERTICES, and a LENS has two edges joining the SAME pair, so both get one number and the last-write-wins scatter destroys one — feed the arc first and both slots hold the chord. Three defects fixed: that numbering; `getNormalConeVertexQ` indexing its second constraint as `j+1` unwrapped, which raised `badsubscript` on any BOUNDED region and so forced its only caller to send every bounded region to the POLYHEDRAL routine (whose cones come from the chord, wrong for a curved edge); and `biconj` handing its second conjugation the curved MESH `conj` now returns, which `quaPolToPlq` refuses — it takes the symbolic form on purpose now. **Measured: the half-lens conjugates to 3 cells exact against a brute-force sup at all 10 probe points (2 identical wrong cells before), and `f**` of `x·y` over the two-face square is exact at 5 of 7 probe points, was 0 of 7.** What remains is a hole in the DOMAIN, not a wrong value: `f**` is `+inf` on `{x+y > 1}`, so one piece of `f*` still conjugates onto too small a set. Pinned (still failing) by `biconjugateTest/biconjugateOverATwoFaceSubdivisionIsTheEnvelope` | `functionNDomain.conjugateOfPiecePoly`, `region.getNormalConeVertexQ`, `biconjCPLQ` |
| `mergeL` / `removeTangent` exact-tie-point gap — **NO LIVE REPRODUCER as of 2026-08-02; do not spend time on it without one** | The only concrete case ever recorded for this entry is `f = x*y` on the unit square as two triangles, whose exact symmetric tie point `s = (0.5,0.5)` used to be excluded from `cplqAdapterTest`. It is exact now, and so is everything around it: a scan of the assembled `f*` over a 25x25 grid on `[-3,3]^2` plus both symmetry axes gives **0 uncovered and 0 wrong** (tolerance 2e-3 against the numeric sup), on `46fac7c` as well as on the current tree — so this was already closed by the `maximumP` fall-through-to-`splitmax3` fix that `functionNDomain.maximumP`'s own HISTORY comment describes, and the entry simply outlived it. What remains is a LATENT defect of the same family, now fixed: `removeTangent` chose `s0(1)` as the second point on the parabola before checking that the root set was non-empty and real (empty threw, complex made both `isAlways` tests undecidable and let the routine mark constraints off a nonsense probe point). Re-open this row only with a case that actually returns a wrong number | vendored cPLQ `functionNDomain` / `plq.biconjugateF`; inherited by `QuaParCPLQ.conj` |
| ~~`QuaPar.eval` exactly **at a vertex**~~ — **FIXED and VERIFIED 2026-08-02** | Cause was `QuaPar.eval`'s exact, no-tolerance point location: `all(vals <= 0, 2)` on a conic that is only zero in exact arithmetic, so a `+1e-17` left the point in NO face and `eval` returned its `Inf` initialization. A conic-magnitude-relative tolerance fixes it. NOW REPRODUCED AND MEASURED, which is what was missing: `sweepQuaParEvalAtVertices(20260802, 200)` finds **225 of 1205** subdivision vertices unlocatable under the exact test and **0** under the current one, with all 7230 probes on rings of radius 1e-8 around them correct — the exact signature recorded for the defect. Deterministically pinned by `QuaParTest/evalLocatesAPointExactlyAtItsOwnVertex`, whose mesh is case 2 of that sweep. Still open, separately: the CURVED half of the original claim (~0.8%) is not covered by this sweep | `QuaPar.eval` |
| ~~`testRegion/testCreation`~~ — **PASSES since 2026-07-28**, and re-measured 2026-08-15 (`testRegion` 23 / 0). Row kept struck through because this file cited it as the suite's single failure for weeks after it had stopped failing | — | — |
| ~~Step 1 ignored the face's function~~ - **FIXED 2026-08-01** | Step 1 now classifies by the SIGNS OF THE EIGENVALUES of `Q`, not by `nCE` (which tests edge slopes and so only classifies `x*y`). Convex/affine keep `q` as their own envelope; concave get the affine interpolant through the actual values of `q`; indefinite are moved by `xyFrame.m` into the frame where `q` **is** `x*y`, so cPLQ's own closed forms apply to the function they were written for. Case C values verified exact against brute force for `(x^2+y^2)/2`, `x*y`, `x^2-y^2`, `3xy+7x-2y+5` and `-(x^2+y^2)/2`, where `f*(0.3,0.4)` used to be `0.4` for a truth of `0.125` | `plq_1p.convexEnvelope1`, `xyFrame.m` |
| ~~Case C's **biconjugate**~~ — **CLOSED 2026-08-15.** The convex case has worked since 2026-08-13, and the INDEFINITE path is the two-face unit square with `f = x·y`, i.e. exactly the row above: `biconjugateOverATwoFaceSubdivisionIsTheEnvelope` passes and `checkBiconjDomainCoverage` re-measures that case as **OK, error 0**. The chain's last named link was fixed with it — `functionNDomain.getInterior` indexed `c2(2)` under a guard testing only `size(c1,2)`; the two `coeffs` results have independent sizes, so any `f` with an `s1²` term but no `s1·s2` term raised `badsubscript`. Where only one of the two equations depends on `s1` there is nothing to eliminate and the other equation IS the relation; where neither does, the old expression is identically zero, which is right for an affine `f` and is what every case reaching there relies on, so that branch is untouched. History follows | `functionNDomain.conjugateOfPiecePoly`, `getInterior` |
| _history of the row above_ | `caseC.biconj()` gives ZERO pieces on pristine `HEAD` - `f** = +inf` everywhere - for an `f` that is convex and hence its own biconjugate. The old test passed anyway: it asserted only `.kind()`, and `QuaParCPLQ(empty).kind()` is still `'QuaParCPLQ'`. It fails inside `conjugateOfPiecePoly` behind a CHAIN of latent bugs, each reachable only once the previous is fixed. Fixed so far: `region.getNormalConeVertexQ` indexed `py(1)` before its own `isempty(py)` guard (dead code); `region.splitmax3` left its output unassigned when `f1 < f2` at every vertex. Next down: `functionNDomain.getInterior` indexes `c2(2)` under a guard testing only `size(c1,2)`. Unrelated to the unbounded / general-quadratic work, which only makes the FIRST conjugation richer (11 pieces vs 9) and so carries the second far enough to reach these | `functionNDomain.conjugateOfPiecePoly` |

_Status 2026-08-15: `conjCPLQTest/biconjCoverageByInputCase` asserts `f** = f` for the convex Case C and passes, because Case B2 answers that domain in closed form and the second conjugation is handed a clean mesh. The chain described in this row is unfixed and still reachable from an INDEFINITE domain._

---

## 7.1 `biconj` domain coverage — measured, 2026-08-02

Re-derive with `checkBiconjDomainCoverage`, committed beside this file. Ground truth owes nothing
to the conjugate pipeline: for a bounded domain it is the lower convex hull of the sampled graph
(`convhulln`), for the unbounded cases it is the identity `f** = f` on a convex `f`.

| domain | `f` | result |
|---|---|---|
| triangle | convex `(x²+y²)/2` | **OK** |
| triangle | indefinite `x·y` | **OK** |
| triangle | concave `-(x²+y²)/2` | **OK** |
| triangle | affine | **OK** |
| box `[0,1]²`, ONE face | `x·y` (McCormick) | **OK**, exact, +inf outside |
| box `[0,1]²`, ONE face | indicator `0` | **OK** |
| axis-aligned box `[0,2]×[0,3]`, one face | `x·y` | **OK** |
| unbounded, full domain | `(x²+y²)/2` | **OK** |
| unbounded, 4 cones | `\|x\|+\|y\|` | **OK** |
| unbounded, 3 wedges | `max(0,x,y)` | **OK** |
| box, TWO faces sharing a diagonal | `x·y` | ~~**WRONG**~~ — **OK, re-measured 2026-08-15**: error 0, 4/4 points inside and 1/1 outside, against the `convhulln` ground truth. This was §7's open defect; see the struck-through row there |
| parallelogram, one face | `x·y` | ~~**ERROR** `QuaParCPLQ:conj:emptyResult`~~ — **COMPUTES since 2026-08-16, and is nearly right**: exact at all four vertices, `+inf` outside the domain, and 8 of 10 interior probe points correct against a brute-force double conjugate. The other two are LOW by about 4% (`0.986` for `1.031`, `1.913` for `1.950`) — an envelope that is still slightly too small, from the residual `getInterior` defect below. Two defects fixed to get here, both general: see §7 |
| general convex quadrilateral, one face | `x·y` | **ERROR** `MATLAB:badsubscript` by default; **OK with `CCA2_A45_SPLIT` set** (2026-08-16) — exact at 10 of 10 probe points against the vertex-attained sup, with no piece leaving a hole, and the fully assembled Step 3 answer exact at 8 of 8. The fix is written and tested; it is opt-in because turning it on costs `testcPLQ` 1542 s → 4728 s and makes `testRectBiconj` error. See below and §8 |

**These are not failures of the ALGORITHM, and the two remaining ones are different defects.**
Read the failing step before drawing a conclusion:

| case | which conjugation | cause |
|---|---|---|
| general quadrilateral | **FIRST** | Step 1 has no `nCE == 3` branch on the path taken |
| ~~parallelogram~~ | ~~SECOND~~ | fixed 2026-08-16 — three per-piece conjugates had holes, so the max was empty; see below |
| ~~two-face box~~ | ~~SECOND~~ | fixed 2026-08-15, §7 |

**The parallelogram — RESOLVED 2026-08-16, and the error message had been blaming the wrong
routine.** `QuaParCPLQ.conj` raised `emptyResult` saying "conjugateOfPiecePoly returned no pieces";
it did not — all **12** pieces conjugated, to **27** cells. What returned nothing was the step
after, `functionNDomain.maxOfList`, and the max was empty for an honest reason: three of the twelve
per-piece conjugates had HOLES inside the parallelogram, and since `f**` is a max its domain is the
INTERSECTION of theirs.

*How the three were identified, which is the reusable part:* `f**` of a bounded domain is finite
exactly ON that domain, so **every** per-piece conjugate must be finite there. Evaluating all 12
groups at six interior points named groups 9, 11 and 12 in one cheap run, turning "the max is
empty" into "these three pieces are wrong".

Two defects, both fixed and both general (§7):

1. **`region.simplifyUnboundedRegion` declared any region with no finite VERTEX empty.** A
   half-plane has none — and a half-plane is exactly what a TANGENT vertex produces, since the cone
   there is built from the two edges meeting at it and tangent edges give the SAME half-plane
   twice. On piece 9 the cone at `(1/4,7/8)` is `{2x/3 + y ≥ 4/3}` carrying `x/4 + 7y/8 − 1/2`.
2. **The edge list, in bug 1's other form.** Piece 9 is bounded with 3 vertices and 3 genuine edges
   plus a conic touching one vertex: 4 constraints for 3 vertices, so the COUNT called it unbounded
   and it was built one edge cell short.

All 12 groups now cover the parallelogram, and `biconj` returns a `QuaParCPLQ` instead of raising.

**What remains, and it is a distinct defect:** the last 2 of 10 interior points come back about 4%
LOW. On piece 9 the chord's edge cell and the arc's claim the SAME region and the chord's is checked
first and is wrong there. `f` is a SINGULAR convex quadratic on that piece — constant along one
whole edge, `∇f = 0` at two of its three vertices — so `functionNDomain.getInterior`, which
eliminates `s` between `x = ∂₁f` and `y = ∂₂f`, returns the gradient map's image LINE rather than a
curve separating the two cells. **Do not attack the `isQuad` chord rewrite for it**: both
alternatives were measured on 2026-08-16 and are recorded in `DECISIONS.md` — chording the vertices
the conic actually touches makes the piece WORSE (2 wrong of 10 → 3), and skipping the rewrite
changes nothing at all.

**THE GENERAL QUADRILATERAL — FIXED 2026-08-16, on the fourth attempt, and the fix is a DOMAIN
split done in EXACT SYMBOLIC arithmetic.** `splitTightTriangleSym` splits a triangle into
sub-triangles on each of which cPLQ's own closed form for THAT sub-triangle's convex-edge count IS
the convex envelope — [COAP] A.4 gives one two-convex-edge half, its own form unchanged and now
tight, plus one one-convex-edge half whose A.3 rational form the `nCE == 1` branch derives
analytically; A.5 gives two two-convex-edge halves that recurse into A.4 — and
`plq_1p.triangulate` emits them as PIECES, which keeps every one of them on a path Step 2 already
has. Sound for Step 3 because a sup over a union is the max of the sups.

**Symbolic is the whole difference from the attempt that hung.** A.4's cevian has slope
`−sqrt(mh·mw)`, so its foot is irrational; taking it from `convEnvCPLQ`'s double-precision faces
gives `2^53` denominators that grow past `1e25` downstream and `isAlways` then decides nothing.
Carried symbolically the coordinates stay compact surds. Measured: the no-split path costs 20 ms —
which is why nothing else in the suite moved — and the split paths 0.3 s (A.4) and 1.2 s (A.5).

**What it costs is Step 3's, not the split's.** Six pieces instead of two, and the cross-piece
maximum then takes **73 minutes** (folds of 93, 294, 647, 1273, 2087 s; cells 5, 14, 29, 45, 70,
86). The answer is exact at both levels. That blow-up has its own `TODO.md` item — 86 cells is
about ten times what the answer needs, and the surplus is same-function neighbours that
`region.merge` never merges.

**And the cost is paid by any `x·y` polygon whose triangles need the split**, because A.4's cevian
foot is irrational: a split sub-triangle has SURD coordinates, so every symbolic operation
downstream works in a quadratic extension rather than the rationals. Measured on `testcPLQ`, whose
domains are general polygons carrying `x·y`: **1542 s with the split off (matching its historical
1427 s), over 3100 s with it on and still unfinished when stopped**, uncontended both times — and
only two of its six domains even gain a piece (2 → 3 and 1 → 2), so this is the algebraic degree of
the coordinates, not the piece count. **The split is therefore OPT-IN, via `CCA2_A45_SPLIT`, and OFF by
default** — and off for a measured reason: with it on `testcPLQ` takes 4728 s instead of 1542 s AND
`testRectBiconj` ERRORS. It trades a documented, LOUD failure on one domain shape for a new one on
another, and until that is understood it cannot be the default. Two things to fix before it can be,
in order: Step 3's cell blow-up, and that error.

Pinned by `cplqAdapterTest/generalQuadrilateralStep1IsTheEnvelopeNotAMinorant`, which asserts only
what must hold of an envelope (it exists, it is `≤ x·y`, and where `x·y ≥ 0` it is `≥ 0`), and
`generalQuadrilateralConjugateMatchesTheSup`. The three failed attempts are in `DECISIONS.md`; what
follows is the diagnosis they produced.

**THE GENERAL QUADRILATERAL IS NOT ONLY A WIRING GAP — measured 2026-08-15 by writing the wiring
and reverting it.** The `nCE == 3` branch below was implemented (build the triangle as a one-face
`QuaPol` carrying `x·y`, call `convEnvCPLQ`, install its faces through `ratPolToPlq`) and **works at
Step 1**: 4 envelope faces come back for the offending triangle — the A.5 split, each half then
needing A.4's — all `≤ x·y`, and `conj` stops raising. The answer is then **wrong**, so it was not
committed. `triangulate` splits the test quadrilateral into piece 1 `[2.5 1.5; 2 0; 0 0]` (`nE = 2`)
and piece 2 `[2.5 1.5; 0 0; 0.5 1]` (`nE = 3`), and with the branch in:

* **piece 2 gets its 4 envelope faces and cPLQ's Step 2 returns ZERO conjugate cells for it** — the
  new envelope is computed and discarded;
* **piece 1, untouched by any of this, is wrong on its own through cPLQ's Step 2**: 6 cells,
  `f*(0,0) = 0.28647` where the truth over its own triangle is `0`, `f*(0.5,1) = 1.00464` where it
  is `1`, and NOT COVERED at `(-1,0.5)` and `(-2,-2)`.

So **the crash was masking a silent wrong answer**, and landing the wiring alone trades a loud
refusal for it. The measurement that localises it: that same triangle conjugated ON ITS OWN via
`QuaPol.conj` is exact at 7 of 7 probes, because a single bounded triangle takes the NUMERIC route
(`conjBoundedPolygon`) rather than cPLQ — **the numeric route is right on this input and the
vendored symbolic one is not**. `assertStep3MatchesPieces` correctly does not fire, since Step 3
agrees with Step 2.

**THE UNDERLYING DEFECT IS IN STEP 1, and it is a silent wrong answer in its own right.**
`plq_1p.convexEnvelope1`'s `nCE == 2` branch applies [COAP]'s single-quadratic form to the WHOLE
triangle. That form is a valid convex MINORANT, but Appendix A.4 shows it is tight only over a
sub-region — which is what `convEnvCPLQ`'s `splitTwoConvexEdges` tests for; this branch never tests.
On `conv{(2.5,1.5),(2,0),(0,0)}` carrying `x·y` it returns an envelope whose minimum over the
triangle is **−0.2835** at `(1,0)`. There `x ≥ 0` and `y ≥ 0`, so `x·y ≥ 0`, the affine minorant `0`
is admissible, and the true envelope is `≥ 0` everywhere: this one is strictly below it. A too-small
envelope gives a too-large conjugate — that is the `0.28647`. `convEnvCPLQ` on the same triangle
returns 2 faces, with minimum `0`.

**And routing Step 1 through `convEnvCPLQ` does not fix it** — tried and reverted:
`conjugateFunction`'s `nCE == 2` branch reads its envelope with `coeffs(...)` and matches monomials,
while A.4/A.5's faces are RATIONAL, so it raises `symbolic:coeffs:NotAPolynomial`. **cPLQ's Step 2
has no rational-envelope branch at all.** The split therefore belongs in the DOMAIN: have
`plq_1p.triangulate` emit the A.4/A.5 sub-triangles as PIECES — the same shape as
`conjEnvelopeViaCPLQ`, which already hands rational faces to cPLQ one piece at a time. `TODO.md`
has the plan and its cost.

**AMENDED 2026-08-20 — the obstacle was the DISPATCH, not the rational faces.**
`conjugateFunction` took its `nCE`, slopes and far vertex from the PIECE, and its per-vertex
COORDINATES from the piece's polygon, so a rational A.3 face was read by the `nCE == 2` branch's
`coeffs(...)`. It now dispatches each face on that FACE's own geometry, and the cross-face max can
split a RATIONAL pair (`region.clearedDifference`). `plq_1p.convexEnvelope1` therefore splits at
the ENVELOPE level too, by calling `splitTightTriangleSym` — the same construction the domain route
uses — for a piece that did not come through `triangulate`. Measured on `{(0,0),(1,0),(2,1)}`
carrying `x·y`: envelope `≥ 0` (was −0.0429), `f*(0,0) = 0` (was 0.0429), `f*(0,1) = 0.125` (was 0),
about 80 s; pinned by `cplqAdapterTest/twoConvexEdgeTriangleEnvelopeIsTightNotAMinorant`. The DOMAIN
split remains the route for a piece that DOES come through `triangulate`, and the cheaper one.

The wiring half of the story, which stands and is ready to re-land once Step 2 is fixed: there are
two Step 1 implementations in this repository:

* **`convEnvCPLQ.m` — CCA2's own.** It *has* the 3-convex-edge case: `splitThreeConvex` cuts the
  triangle (in the bilinear frame) through the middle vertex into two 2-convex-edge sub-triangles,
  [COAP] Appendix A.5. This is the "retriangulate first" the method calls for, and it is exactly
  what a general polyhedral set needs. `biconj`'s Case B and the SCIP bridge both use it.
* **`plq_1p.convexEnvelope1` — the one MERGED IN from cPLQ.** It branches on `nCE == 0`, `1`, `2` and
  then simply falls off the end: for `nCE == 3` it sets no envelope and never sets `lCE`, so
  `obj.envelope` stays EMPTY. `plq_1p.conjugateFunction`'s `for i = 1:max(1, size(envelope,2))`
  — a guard written for "triangles where the convex envelope is not computed" — then indexes
  `envelope(1)` and raises `MATLAB:badsubscript`.

Case C (`conjCPLQ.m`: `quaPolToPlq` → `triangulate` → `maximum`) drives Step 1 through the
**merged-in** one. So the split CCA2 already implements is simply not reachable from `conj`/`biconj`
on a multi-vertex domain. Routing Case C's Step 1 through `convEnvCPLQ`, or applying
`splitThreeConvex` to the pieces before handing them to `plq_1p`, is the fix.

`nCE` counts edges of positive finite SLOPE, which is why an axis-aligned box (all pieces
`nCE = 0`, affine envelopes) sails through while a sheared one does not — the *symptom* is
orientation-dependent even though the missing branch is not.

---

## 8. Summary — what actually blocks a general release

Ordered by how likely a downstream caller is to hit it:

1. **`partialConj` is entirely unimplemented** (§2).
2. ~~**Unbounded multi-face domains**~~ — **RESOLVED 2026-08-15.** Steps 1, 2 and 3 are all done.
   This was blocker 2 from 2026-08-01; the last two pieces closed on the same day, and both are
   below for the record.
   **(a) DONE.** `quaPolToPlq` builds ray-carrying faces from HALF-PLANES
   (`faceDomainFromHalfPlanes`), orientation off `P{k}`'s own sign convention, so the ray
   direction is no longer discarded.
   **(b) DONE.** `plq_1p.triangulate` covers an unbounded face with triangles + half-strip +
   wedge (`fanUnboundedFace`); `convexEnvelope1` gets its envelope from `convEnvUnbounded`;
   `conjugateFunction` dispatches on the ENVELOPE, not on `nCE`.
   **(c) DONE, by REPAIR rather than by the earlier workaround.** The dual regions were wrong
   because `poly2orderUnbounded` — the routine that puts an unbounded region's vertices in
   boundary cyclic order, which `getNormalConeVertex` requires — **threw** on a 2-constraint
   wedge, and `getEdges` raised `MATLAB:unassignedOutputs` instead of returning an empty list
   when no constraint is active at a point (true of the box-clip corner `(intmax,intmax)`). Both
   are fixed in place, so every vertex-based consumer benefits, and the parallel half-plane
   construction added as a stopgap has been deleted. The ray REPRESENTATION was never at fault.
   **(d) DONE.** A CURVED convex envelope over an unbounded face is conjugated by
   `conjConvexOverPiece.m` (KKT active set: vertex, edge and interior cells). The handoff's
   headline case — `(x²+y²)/2` over `{x≤0,y≥0}` — returns `min(s1,0)²/2 + max(s2,0)²/2`, exact
   at 10 probes. The same routine supplies the edge/interior cells cPLQ's `conjugateOfPiecePoly`
   omits for a convex `q` on a bounded triangle.
   **(e) RESOLVED 2026-08-15.** cPLQ's Step 3, the CROSS-PIECE maximum
   (`plq.maximumConjugate` → `functionNDomain.maximumP`), disagreed with its own per-piece
   conjugates on the 4-cone fan with convex faces. It had two halves, and **both were the same
   trap**: a routine concluding something from a direction built at a vertex where the split
   conic's GRADIENT VANISHES — the apex of a cone, which is exactly where an unbounded fan's split
   conics meet, and where no such direction exists.
   * **The DROP, fixed 2026-08-02.** The assembled maximum kept only 4 of the 16 cells, losing
     face 1's `s₂²/2` cell on `{s1≤0, s2≥0}`, so `f*(-0.5,2)` came back `1.125` for a truth of
     `2`. `region.simplifyUnboundedRegion` declared a genuine 2-D cone EMPTY from probe directions
     built out of constraint SLOPES at that vertex. Fixed by `region.witnessAwayFrom`, which
     refutes an emptiness verdict by exhibiting a feasible point — sound, because a genuinely
     empty region has none.
   * **The OVER-CLAIM, fixed 2026-08-15.** At `s = (-3,-2.4)` the assembly gave `5.130` where the
     per-piece max gives `4.500`. The cell carrying `s1²/4 + s2²/2` had lost its `−s1 ≤ 0`,
     keeping only `{s2 ≤ 0, s2²/2 − s1² ≤ 0, s1² − 2s2² ≤ 0}` — two constraints **blind to the
     sign of `s1`** — so the region was symmetric under `s1 → −s1` and claimed the mirror wedge.
     `region.removeTangent` deleted it, after building a TANGENT LINE to a quadratic at its own
     apex. Fixed by refusing to conclude anything from a vanishing gradient.
     Cleared on the way, so they are not re-checked: `redundantSubset` certifies nothing redundant
     there (correct), `simplifyUnboundedRegion` leaves the constraint alone, and Step 2 is right —
     each piece's own conjugate has the correct quadrant constraints.
   **`conjCPLQTest` is 25 / 0.** `step3DropsCellsOnSomeUnboundedAssemblies`, which pinned the gate
   FIRING, is renamed `step3UnboundedAssemblyAgreesWithItsOwnPieces` and now pins what the gate
   protects; its sibling `step3UnboundedAssemblyMatchesTheTruth` pins the closed-form values.

   The `assertStep3MatchesPieces` gate STAYS: it is an independent cross-check of the whole
   assembly against the same `f*` computed the other way, and the failure it guards against is
   silent by nature. It simply no longer fires on this input.

   **(f) DONE 2026-08-15 — the CURVED convex envelope over an unbounded face.** `convEnvUnbounded`
   computed only the AFFINE envelope and raised `convexAlongRay` as soon as `d'Qd > 0` along a ray.
   Two shapes are now derived and implemented, each with its proof in the source: a WEDGE with one
   flat and one convex ray, whose envelope is `q` with its CROSS TERM deleted; and a HALF-STRIP
   convex along the ray whose base edge is Q-ORTHOGONAL to it, where `q` separates and the envelope
   is the affine interpolant along the concave base plus the convex part along the ray. A negative
   cross term on the wedge means the envelope is `−inf`, now reported rather than answered; a
   non-orthogonal half-strip is refused loudly. `conjConvexOverPiece.m` already conjugates such an
   envelope. **`unboundedFaceTest` 18 / 0**, from 16 / 2.

3. **`'pqp'` and `'graph'` engines missing** (§1.1).
4. **`RatPol.conj`/`biconj`/`add` missing** (§3, §5).
5. ~~**One known wrong-answer defect** — `plq_1p.convexEnvelope1`'s `nCE == 2` branch returns a
   convex MINORANT rather than the ENVELOPE~~ — **RESOLVED 2026-08-16** by splitting the DOMAIN
   along [COAP] A.4/A.5 in exact symbolic arithmetic (§7.1). What is left of that entry is a
   PERFORMANCE problem, not a wrong answer, and it is now the binding one — see 5b.
5b. **Step 3's cross-piece maximum does not scale.** Measured on `x·y` over
   `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}`: Steps 1 and 2 take about 25 s for all six pieces, and
   `functionNDomain.maxOfList` then takes **73 minutes**, with the cell count running 5, 14, 29,
   45, 70, **86**. The same shape is on record for a pentagon (885 s, 41 regions). 86 cells is
   roughly ten times what the answer needs — `f*` of `x·y` over a convex quadrilateral has a cone
   per vertex and a cell per edge — and the surplus is adjacent cells carrying the SAME function
   that `region.merge` never merges. `TODO.md` has where to start.
5c. **`getInterior` on a SINGULAR quadratic** — the one remaining wrong ANSWER, and it is small:
   about 4% low at 2 of 10 interior points of the parallelogram's `f**` (§7, with a one-minute
   reproducer).
6. ~~**`maxQuaPar`: a piece whose two arcs are ADJACENT**~~ — **RESOLVED 2026-08-15.**
   `splitTwoArcPiece`'s two candidate chords join the arcs' facing endpoints, which for arcs
   sharing a vertex ARE the arcs' own edges, so both chains came out too short, the piece was
   returned unsplit with one arc flattened to its chord, and assembly reported an orphan three
   stages later. Generalised the `nv == 3` shared-vertex fallback to `nv >= 4` with the ordinary
   diagonal to a non-adjacent vertex. **The seeded sweep is now 18 exact / 0 wrong / 0 errored of
   18** (it was 17 / 0 / 1, and 16 / 0 / 2 before 2026-08-14). **`maxQuaPar` has no open case.**
   ~~an unbounded half carrying two arcs~~ and ~~an unbounded piece straddling `{f1=f2}`~~ — both
   **RESOLVED 2026-08-14**, see §4.1.
7. ~~**arc-vs-arc results are only locally correct (wrong far from the arcs)**~~ — **RESOLVED
   2026-08-13** (§4.1). This was the top blocker on this list from 2026-08-04. The cause was the
   point-location rule, not the subdivision: a curved edge is a bounded arc and its conic is not,
   so "every bounding conic, sign-oriented, ≤ 0" admits points arbitrarily far away on a
   parabola's concave side. `QuaPar.chordCuts` derives the missing chord per face. Four other
   silent-wrong-answer defects were fixed alongside it. Re-measured: 0 of 64 far-field results
   wrong where 7 of 64 had been, and 0 of 1031 vertices / 571 midpoints / 3540 interior points
   wrong on the committed sweep.
8. ~~**`maxQuaPar` cannot split a cell that already carries an arc**~~ — **RESOLVED 2026-07-30**
   (§4). It needed no conic-conic solver and no multi-arc representation, contrary to what this
   file and `maxQuaPar.m`'s own TODOs used to say. Every curved edge is a parabola
   (`QuaPar.assertParabolic`), so restricting the splitting conic to the arc via
   `parabolaArcFrame.conicCoeffs` gives one univariate quartic; and the splitting curve never
   CROSSES the arc here (measured: 19 untouched, 3 tangent, 0 crossed over 22 curved splits), so
   ONE ARC PER FACE is preserved by subdividing with a straight chord. Assembled results went
   58 → 76 of 395 sampled quadrilaterals, all exact to 2.8e-14 and arrangement-clean.
9. Performance: general bounded domains route through the symbolic pipeline (Phase 2).

**Resolved 2026-07-28:** `biconj` works for every bounded domain, including the single bounded
triangle that used to be blocker 1 — via Step 1's convex envelope rather than a second
conjugation. See §3 and `biconjCPLQ.m`. Also `conj` of a 2-convex-edge split triangle, by falling
back to cPLQ's Step 2/3 (§1.2, `ratPolToPlq.m`). That work also fixed four latent bugs in the
vendored symbolic layer, all previously unreachable because the `0/0` errored first:
`plq_1p.m`'s `nCE==2` branch left the coefficient `d` unassigned when the envelope had no linear
term; `region.vertexOfEdge` and `region.simplifyUnboundedRegion` substituted into a removable
`0/0` instead of taking its limit; and `region.plus` indexed into an empty region array.

**Resolved 2026-07-27:** `conj`'s return type no longer varies by input shape — every type now
inherits from `RatPar` and `kind()` reports the concrete one. See `RETURN_TYPE.md` and `RatPar.m`.

Items **not** on this list, deliberately: every **N/R** row above. They are unreachable from a
`QuaPol` conjugate/biconjugate and require no work.

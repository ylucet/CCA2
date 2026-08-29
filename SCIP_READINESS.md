# Readiness plan — what must be true before the SCIP/QPLIB work starts

_Written 2026-08-18. The gate has three conditions: **the bugs are ironed out**, **it is clear
what has a direct formula and what needs symbolic computation**, and **timing is understood and
acceptable**. We are not there. This file says what remains, in the order it has to happen, with
the measurement that closes each step._

The consumer is `AI/spike/SCIP` (`src/cca2ConvexEnvelope.m`, via the MATLAB Engine API for
Python). `SUPPORT_MATRIX.md` sections 0.0 and 0.0.1 are the standing description of that
interface.

---

## 0. The finding that shapes the whole plan

From `SUPPORT_MATRIX.md` 0.0.1, and it decides what is worth building:

> For a BILINEAR term over a BOX the convex envelope **is** McCormick (Al-Khayyal-Falk), in
> closed form -- which is what CCA2 returns, and the pin
> `biconjugateTest/bilinearOverABoxGivesTheMcCormickEnvelope` holds it. On QPLIB's box-domain
> bilinear terms CCA2 is a 40-second reimplementation of a formula SCIP already applies:
> correct, a good validation, **and no stronger as a cut**.

CCA2 can only beat McCormick in two places:

* the domain is **not a box** -- a triangle cut out by constraints;
* the piece is **not bilinear** -- a diagonal term such as `x^2 - y^2`.

**Both of those are exactly what fails or is slow today.** So the readiness work is not
"make the existing path production-ready"; it is "make the two cases that would justify the
approach work at all, and know what they cost". A SCIP run built only on box+bilinear terms
would be green and scientifically empty.

---

## Phase A -- correctness (the "no bugs" gate)

### A1. Re-measure 0.0.1. It is stale, and nothing should be planned on it.

Every number in that section dates from **2026-08-02**. Since then the conjugate pipeline changed
substantially: three double leaks fixed (Step 2 is exact), `merge` repaired, the A.4/A.5 split
turned on by default. The two ERROR rows may have moved in either direction, and the timings
certainly have.

* **Do:** run `checkBoxEnvelopeForSCIP` and rewrite 0.0.1 from its output.
* **Closes when:** the table is re-derived, with its date, on the current tree.
* **DONE 2026-08-18.** All six rows OK at 0 s, every one returning a MESHED `QuaPol`; both former
  ERROR rows pass. 0.0.1 is rewritten from that run (log: `.claude/boxenvelope.log`).

### A2. Fix the DIAGONAL terms over a box -- the one correctness gap on SCIP's path.

`x^2 - y^2` and `(x^2+y^2)/2` on the unit box raise `MATLAB:badsubscript` in the SECOND
conjugation (`functionNDomain.conjugateOfPiecePoly`). QPLIB objectives are sums of `x_i*x_j`
**including i = j**, so this is not optional -- and per section 0 it is one of only two cases
where CCA2 has anything to offer over McCormick.

* **Do:** reproduce at unit level (build the piece directly, as `functionNDomainTest`'s fixtures
  do -- not through `biconj`, which buries the evidence), then fix.
* **Closes when:** both shapes return an envelope matching the lower convex hull of the sampled
  graph, pinned by a test in `biconjugateTest`.
* **DONE — and there was no mathematics in it.** `x^2 - y^2` is SEPARABLE as written (envelope
  `x^2 - y`) and `(x^2+y^2)/2` is CONVEX, so both are answered by short-circuits in 0 s and
  neither reaches the second conjugation where `MATLAB:badsubscript` was raised. Pinned by
  `biconjCPLQTest.separableOverABoxTakesTheOneDimensionalRoute` and
  `convexOverABoxIsItsOwnBiconjugate`; measured again in A1's run above.

### A3. Decide the `unionIsExact` question -- an unknown, not a known bug.

At fold 3 of the quadrilateral, 52 same-function pairs reach `unionIsExact` and about 9 merge.
Whether the other ~43 refusals are **correct** is not established: two cells can share a facet,
touch along a segment, and still have a non-convex union.

* **Do:** take a handful of the 46 `sym=1 hyp=1 touch=seg` pairs from `.claude/step3adjacency.m`
  and test directly whether the union is convex.
* **Closes when:** each is classified correct or defective. If correct, the cell counts are near
  optimal and Phase C's target changes; if defective, that is the next fix.
* **Why it sits in the correctness phase:** three times this session a CORRECT result was read as
  a defect (`DECISIONS.md`). Do not optimise this gate before knowing which it is.
* **DONE 2026-08-19, and it is CORRECT, not defective.** 22 `unionIsExact` calls captured over
  three folds and scored with a DECISIVE oracle (`.claude/a3score.m`): merge is wrong exactly when
  it LOSES a point, since `M = A' ∩ B' ⊆ A ∪ B` always. **Zero defects** -- no accepted merge
  loses a point. Of the ten refusals, five are PROVEN correct by an exact witness, and the other
  five all fail for the same reason, a CURVED certificate (`certifiesNonPositive`), with no lost
  point found in 2e6 samples. So the LP verdicts are right, the cell counts are not hiding an
  over-claim, and the only sharpenable gate is the curved certificate -- worth at most 5 of 22
  calls. Details in `DECISIONS.md` 2026-08-19.

### A4. State what is deliberately NOT fixed, and confirm it is off SCIP's path.

`SUPPORT_MATRIX.md` section 8 lists open items -- `partialConj` unimplemented, the `pqp` and
`graph` engines missing, `RatPol.conj`/`biconj`/`add` missing. Section 0.0 already records that
SCIP calls none of them.

* **Do:** confirm against the current `cca2ConvexEnvelope.m` and write the exclusion into 0.0.
* **Closes when:** the SCIP-relevant surface is a short, explicit list of entry points.
* **DONE 2026-08-28.** `SUPPORT_MATRIX.md` §0.0 re-confirmed against the current tree: `biconj('cplq')`
  (the recommended entry point per gap 1 below) routes through `biconjCPLQ.m` and never reaches
  `partialConj`, `pqp`/`graph`, or `RatPol.conj`/`biconj`/`add` — same exclusion list as the current
  `convEnvCPLQ` bridge, so switching entry points (a spike-side task) does not change CCA2's surface.

---

## Phase B -- the direct-formula / symbolic map (the "clear picture" gate)

**BOTH ITEMS ANSWERED 2026-08-18.** The question was not academic: anything with a direct formula
should not be paying for the symbolic engine, and anything that needs the engine has a cost floor
no tuning removes. The map below says which is which, and B2 names the entry point.

Read Phase B as a MAP, not as a programme of work: it decides which shapes need the symbolic
engine. Replacing symbolic computation inside the one shape that does need it is Phase C (and the
separate `solve()`-removal programme in `DECISIONS.md`), not this gate.

### B1. The map -- MEASURED 2026-08-18 (`.claude/phaseBmap.m`)

One row per input SHAPE: the route `biconj` takes, whether it is closed form, where the formula
is, and what it costs. Times are the **minimum of three** repetitions and the machine was running
the slow test bucket throughout, so treat them as upper bounds (CLAUDE.md section 3); the ROUTE
and the returned CLASS are deterministic given the shape and are not contention-sensitive.

| shape | route | form | where the formula is | cost | returns |
|---|---|---|---|---|---|
| bilinear over a BOX, one face | McCormick short-circuit | closed | `biconjCPLQ.mccormickEnvelope` (Al-Khayyal-Falk) | 0.01 s | `QuaPol` |
| bilinear over a BOX, two triangles | Step 0, then McCormick | closed | `mergeSameQuadFaces` + the above | 0.01 s | `QuaPol` |
| diagonal indefinite over a BOX (`x^2-y^2`) | separable, 1-D per axis | closed | `biconjCPLQ.separableEnvelopeCoefs`/`oned` | 0.00 s | `QuaPol` |
| convex quadratic over a BOX | `co f = f` | closed | `biconjCPLQ.convexEnough` | 0.00 s | `QuaPol` |
| bilinear over a DIAMOND | rotate 45 deg, separate, rotate back | closed | `biconjCPLQ.diamondEnvelope` | 0.00 s | `QuaPol` |
| TRIANGLE, affine or convex | `co f = f` | closed | `biconjCPLQ.convexEnough` | 0.00 s | `QuaPol` |
| TRIANGLE, concave | Step 1, affine interpolant | closed | `convEnvCPLQ`, [COAP] A.2 | 0.02 s | `RatPol` |
| TRIANGLE, indefinite 1 convex edge | Step 1 | closed | `convEnvCPLQ`, [COAP] A.3 eq.16 | 0.00 s | `RatPol` |
| TRIANGLE, indefinite 2 convex edges | Step 1 | closed | `convEnvCPLQ`, [COAP] A.4 | 0.02 s | `RatPol` |
| TRIANGLE, indefinite 3 convex edges | Step 1, A.5 split into A.4 | closed | `convEnvCPLQ.splitThreeConvexEdges` | 0.01 s | `RatPol` |
| general POLYGON, indefinite | triangulate, per-piece closed form, **Step 3 symbolic max** | SYMBOLIC | `functionNDomain.maximumP` | **2579 s** (A.4/A.5 quadrilateral, cited from the profile below) | `QuaParCPLQ` |
| UNBOUNDED multi-face, convex (3 wedges) | no short-circuit fires -- full `conj(conj(f))` | SYMBOLIC | -- | **26-28 s** (three runs) | `QuaParCPLQ` |

**The expected conclusion is CONFIRMED: every single PIECE is closed form.** Nothing in the first
ten rows touches the symbolic engine, and none of them costs more than 0.02 s. All symbolic cost
is the cross-piece work -- Step 3's maximum, and the second conjugation that follows it.

**Two things the map showed that the prediction did not.**

1. **A convex MULTI-FACE input pays 26 s for the answer `f`.** `co f = f` is the largest
   short-circuit there is, but `convexEnough` can only PROVE convexity for a SINGLE piece;
   for several faces it needs the caller's `fIsConvex` flag, because the honest test requires the
   gradient jump across every shared edge (biconjCPLQ.m). Unset, the three-wedge
   `max(0,x,y)` -- convex, and its own biconjugate -- goes the whole way round. Setting the flag
   is free and is the entire fix; that is a caller-side lever, not a pipeline one.
2. **The two-triangle spelling of a box is now the one-face spelling** (Step 0, 2026-08-18), so
   the shape tests in this table are reached by the meshes callers actually build, not only by
   the tidiest one.

### B2. The consequence for SCIP -- the entry point is `biconj`, and Step 3 is never reached

**Answer: a separator computing one term's envelope over a box never enters Step 3, and needs no
per-piece entry point of its own.** The short-circuits live INSIDE `biconjCPLQ`, so the entry
point a separator should call is simply

    h = q.biconj('cplq');        % q a QuaPol: the term, on its box

and for every shape QPLIB's box-domain terms present -- bilinear over the unit box, bilinear over
a sub-box after branching, and the diagonal `x_i^2` terms -- that call returns in **0.01 s or
less**, in closed form, as a **meshed `QuaPol`** the bridge can read directly. The plan's own
premise (call the per-piece closed form rather than `biconj`) is therefore obsolete: `biconj`
IS the per-piece closed form for these shapes.

Two consequences follow, and they point in opposite directions.

* **The performance blocker on SCIP's path is gone.** The "40-60 s per term" of section 0.0.1
  (2026-08-02) predates the short-circuits. At 0.01 s per term, `QPLIB_1940`'s 288 objective
  off-diagonal terms are seconds, not 4 hours, and even its ~27,586 constraint terms are minutes.
  Phase C's target-setting must be re-scoped around that -- see C1.
* **What remains expensive is exactly what is scientifically interesting.** Step 3 is reached only
  when the domain is a genuine multi-piece polygon whose envelope COUPLES the pieces -- the
  non-box case, where CCA2 has something McCormick does not, and where one quadrilateral term
  costs 43 minutes. Section 0 said the box+bilinear run would be "green and scientifically empty";
  the map says it is now also nearly free, which changes nothing about that judgement.

* **Closes:** B1 and B2 are both answered above, with the entry point named.

---

## Phase C -- timing (the "acceptable performance" gate)

### C1. Set the target BEFORE optimising.

No threshold is stated today. The measured facts: 40-60 s per term (stale, pre-session);
`QPLIB_1940` has 288 objective off-diagonal terms, about 4 h offline at that rate, with ~27,586
more in its constraints; per-node recomputation is out of reach at any rate near this.

* **Do:** state the target as a decision -- e.g. offline pre-computation of one QPLIB instance's
  objective must fit in X minutes, and per-node recomputation is out of scope for the spike.
* **Closes when:** the number is here and `AI/spike/SCIP/PROJECT_PLAN.md` agrees.

### C2. Attribute the cost.

Split one term's time across Step 1, Step 2, Step 3 and the MATLAB-Engine round trip.

* **Closes when:** a per-stage breakdown exists for the box case and one non-box case.

### C3. Only then optimise, guided by B2.

If B2 says the per-piece closed form suffices, this is a routing change and probably ends the
performance problem outright. If Step 3 is genuinely required, the levers are the ones already
measured: cell count (86 to 60 so far) and the merge gates.

---

## Gate -- what "ready for SCIP" means

**REWRITTEN 2026-08-29** -- the version below (as of 2026-08-18) predated an entire session's
work landing three real `region.m` bug fixes, closing G1/G4/G10 in `assemblePieces`, and
correcting several load-bearing assumptions about what "ready" actually requires. All four
conditions, each with a measurement on the CURRENT tree:

1. ~~`checkBoxEnvelopeForSCIP` shows **no ERROR rows** (A1, A2).~~ **MET 2026-08-18** -- six rows,
   no ERROR, every one 0 s and meshed. Unaffected by anything since; not re-measured.
2. ~~The direct-formula / symbolic map exists and names the entry point a separator should call
   (B1, B2).~~ **MET 2026-08-18** -- the map is B1, the entry point is `q.biconj('cplq')`, and on
   SCIP's own shapes it answers in 0.01 s without reaching Step 3.
3. **Timing is MEASURED and COMPARED, not held to an arbitrary target** (revised framing,
   2026-08-28: the point is to know the cost and compare it to plain SCIP, not to clear a number
   picked in advance). Measured extensively on the reference non-box fixture (the A.4/A.5
   quadrilateral): `maximumP` on one fold went from 195.9 s (2026-08-18) to 173 s
   (2026-08-28) from work landed since, then to as low as 2226-2341 s TOTAL for the full 5-fold
   run depending on machine load (AI/CLAUDE.md sec 3 -- direction, not magnitude). **Not yet
   measured against a real QPLIB-shaped non-box constraint** -- everything above used the one
   reference fixture; that comparison is still open and is the number that actually answers
   "fast enough for QPLIB", not another pass on the same fixture.
4. **Full suite green -- NOT MET.** Known reds, both attributed and both SCIP-relevant (neither
   is a stray, unexplained failure):
   - `testcPLQ/rectBiconjugateIsAConvexUnderestimator` (verylong) -- does not finish. Confirmed
     2026-08-28 to exercise the EXACT SAME legacy Step 3 machinery
     (`functionNDomain.maxOfList`/`mergeL`/`region.maximum`) that `biconj('cplq')`'s non-box
     Case C path reuses, so this is a genuine SCIP-relevant performance ceiling, not an
     unrelated legacy-only concern. Not yet formally quarantined (§8's own rule: name it, state
     why, never leave it silently red -- see TODO below).
   - Slow bucket not re-run in full this session before 2026-08-29; a fresh run is the other
     half of closing this gate item.

Only then: wire the bridge, expose value and subgradient off whatever B2 names, and run QPLIB
(spike/SCIP's job, not this project's, per the umbrella working agreement).

---

## What changed since 2026-08-18, for anyone reading this file cold

- **Three real, independently-verified bugs fixed in `region.m`'s `maxAffineOverRegion`**
  (different-axis convex conics, same-axis-opposite-sense convex conics, and a sign-gap in the
  original mechanism) -- each previously returned `Inf` where the true value was finite. None of
  the three moved the reference fixture's own cell count (measured three times, byte-identical),
  but each is a genuine correctness fix on its own witness and may matter on QPLIB's actual
  constraint shapes, which this project has not yet tried.
- **G1/G4/G10 landed** (the parked `assemblePieces` diff, `collapseTinyEdges` +
  `matchHalfEdges`'s sagitta test) -- a deliberate, user-authorized trade-off: fixes the
  numeric-path assembly defect at the cost of one symbolic Case-C fallback input (sweep case 21)
  hitting the known Step 3 legacy gap instead of refusing fast.
- **The scaling defect's real cause found**: a high-degree hub vertex in the dual arrangement
  where many same-function cells from different pieces meet at one point without sharing a real
  edge, so pairwise folding re-tests far more candidates than the answer needs. The fix is a
  fold-STRATEGY change (resolve each hub's fan once, not pairwise) -- a genuine algorithm change,
  not attempted. This is what actually determines non-box timing for real QPLIB constraints, more
  than any single bug fix does.
- **Conic-conic exact intersection infrastructure exists** (`conicMeet.m`/`ratQ.m`, tested,
  12/0) but is wired into nothing -- built for a `QuaCon` class that exists only as design docs
  (`doc/QuaConExample.md`, `CONJ_FIELD_PROOF.md` sec 6), not as code. Corrects this file's
  earlier read of `region.getVertices`'s remaining `solve()` cost as "needs a quartic
  implementation" -- it needs `QuaCon`, or a proven-safe adapter, neither of which exists yet.

Full detail, evidence, and the mechanisms ruled unsafe along the way are in `DECISIONS.md`
2026-08-28/29 and `TODO.md`'s G1/G4/G10 and item-3 entries.

---

## Sequencing

A1 is cheap and goes first -- it may change A2's shape. B1 can run in parallel with A2, being
mostly reading and measuring rather than fixing. C is last, because B2 may remove most of the
problem C is trying to solve.

**Update 2026-08-18: B is done, and it did remove most of C's problem** -- on SCIP's own shapes
the per-term cost is 0.01 s, so C1's target must be re-scoped around the NON-box case, which is
the only one that still reaches Step 3. A1 and A2 are also close: `x^2 - y^2` over the unit box
now returns in 0.00 s through the separable route (B1's table;
`biconjCPLQTest.separableOverABoxTakesTheOneDimensionalRoute` pins it), so section 0.0.1's two
ERROR rows and its "no mesh" headline both need re-deriving.

**Update 2026-08-29: sequencing is effectively done.** A1/A2/B are all closed; what remains is
gate item 4 (quarantine the known reds, confirm the slow bucket) and the fold-strategy work that
would make gate item 3's non-box timing actually competitive, not just measured.


---

## Phase C evidence: PROFILE OF ONE FOLD (2026-08-18)

_Filed under Phase B when it was written, which made that gate read like an optimisation
programme. It is not: Phase B decides WHICH shapes need the symbolic engine (answered above), and
everything below is about making the one shape that does need it cheaper -- C2's attribution and
C3's levers._

Two folds of Step 3 on the A.4/A.5 quadrilateral, under MATLAB's own profiler so that nothing in
the repository is instrumented. **Call counts are the measurement** -- the machine is shared, so
absolute times are contended, but counts are not, and the two folds give the SCALING.

### It refuted the static prediction

The Phase B table guessed that `mergeL`'s pairwise "are these the same function" test -- an
`isAlways(simplifyFraction(...))` per pair, O(n^2) -- was the cost. **It is not.**
`simplifyFraction` runs 4735 times in fold 2 and costs **3.1 s** of a 200 s fold. The equality
test is cheap.

### What actually costs

    engine call     fold 1     fold 2      of fold 2's ~200 s
    subs              3215      11803      85.0 s   <-- the single largest
    isAlways          3664      11793      18.3 s
    simplify           914       2256      19.4 s
    solve              203        700      21.5 s
    simplifyFraction  1401       4735       3.1 s

    routine          calls (f1 -> f2)   time in fold 2
    mergeL              2  ->  2          165.2 s   (83% of maximumP's 199.8 s)
    getVertices        57  -> 131         133.1 s
    region ctor        31  ->  75          89.3 s
    merge              32  -> 130          66.2 s
    ptFeasible       1204  -> 3538          40.7 s

**The chain is: `merge` -> build a candidate `region` -> `getVertices` -> `ptFeasible` per
candidate -> `subsF` -> `subs`.** Substitution is 7.2 ms per call because each one is a round
trip to MuPAD, and there are 11803 of them in one fold.

### Scaling confirms it is quadratic in the cell count

Cells grow 12 -> 23 (1.9x) while the engine calls grow 3.2-3.7x and `merge` calls 4.1x -- i.e.
~n^2, which is what trying to merge every same-function pair predicts.

### The revised lever, in order

1. **Make `ptFeasible` numeric-first.** It substitutes symbolically and then tests a SIGN. Now
   that Step 2 is exact, evaluating in double precision and falling back to the symbolic path
   only when the value is near zero is sound -- the classic filter. 3538 calls per fold, driving
   most of the 11803 substitutions.
2. **Finish `getVertices`.** The affine x affine pair is already closed form (a determinant, done
   for exactly this reason -- it was 322 of 438 solve calls). Affine x conic and conic x conic
   still call `solve`: 700 calls, 21.5 s in fold 2, and both have textbook closed forms.
3. **Do not rebuild the region to test a merge.** 130 merge attempts per fold, few of which
   succeed, and each one constructs geometry.

Only after those is the O(n^2) pair count itself worth attacking.


### The ptFeasible numeric filter, measured (2026-08-18)

Implemented and re-profiled on the same two folds. **Counts are the measurement**; the machine was
contended in both runs, so the wall-clock column is indicative only.

    fold 2 (23 cells)        before      after     change
    subs                      11803       9445      -20%
    isAlways                  11793       8683      -26%
    ptFeasible own time        40.7 s     32.2 s     -21%
    maximumP total            199.8 s    195.9 s      -2%
    fold wall time              289 s      280 s      -3%

    fold 1 (12 cells)
    subs                       3215       2512      -22%
    isAlways                   3664       2657      -27%

**It works, and it is smaller than hoped.** The filter removes a fifth of the substitutions and a
quarter of the `isAlways` calls, exactly where predicted -- but `ptFeasible` was only 20% of the
fold, so cutting it by 21% moves the fold by ~4%, which is inside contention noise. Kept, because
the call-count reduction is real, machine-independent and sound; not claimed as a speed-up.

**What the same profile now says to do next**, fold 2, after the filter:

    getVertices     125.6 s     (still 700 solve calls -- untouched by the filter)
    region ctor      84.3 s
    merge            66.1 s
    subsF            61.2 s  (4471 calls)

`getVertices` is the target. Its affine x affine pair is already a determinant; affine x conic and
conic x conic still call `solve` 700 times per fold, and both have textbook closed forms -- the
line substituted into the conic gives a quadratic in one variable.

**Checked 2026-08-28: affine x conic is DONE** (`region.m:4527` on, `region.lineMeetsConicSym`) --
the only shape left calling `solve()` in `getVertices` (region.m:4588) is conic x conic. That is
NOT the same textbook closed form: two general conics meet in up to 4 points (Bezout), a quartic,
not a quadratic, and `DECISIONS.md` 2026-08-20 (T2b) already found one vertex of `f*` with Galois
group S4 -- "no tower of square roots is enough" -- so a from-radicals closed form here would need
more than the square-root extensions this codebase's exactness model uses elsewhere, or Ferrari's
full quartic formula (cube roots), neither attempted. `CONJ_FIELD_PROOF.md`'s naming-a-vertex-by-
its-conic-pair route sidesteps needing coordinates at all, but that is a representational change to
`region.m`, not a drop-in replacement for this one `solve()` call -- not started. Left for whoever
picks this up: it is NOT a quick win, and re-deriving that S4 finding from scratch would repeat
already-refuted work.

### Re-measured on the current tree (2026-08-28), fold 2 only (`CCA2_STEP3_FOLDS=2`)

Same fixture (the A.4/A.5 quadrilateral), same fold-2 cell count (23, unchanged) -- a clean
before/after since the input geometry did not move. `maximumP`: **173 s**, against 195.9 s the
last time this exact fold was measured (2026-08-18, after the `ptFeasible` filter). A real,
machine-independent improvement (cell count and call-count reductions from `getVertices`'
closed-form affine x affine path and the two-conic recession-ray extension both landed since),
though per `AI/CLAUDE.md` §3 not claimed as a precise percentage on a shared machine -- direction,
not magnitude. `getVertices`'s remaining `solve` calls (affine x conic, conic x conic) are still
the next lever, unchanged from the 2026-08-18 recommendation above; not attempted this session.

**Why this matters for SCIP, stated plainly (2026-08-28).** This fold is not an isolated
benchmark: `QuaParCPLQ.conj` -- reached by `biconj('cplq')` whenever a QPLIB term's domain is
non-box (Case C, `conjCPLQ.m`) -- "reuses `plq.biconjugateF`'s own recipe verbatim"
(`QuaParCPLQ.m:59`), i.e. the SAME `functionNDomain.maxOfList`/`mergeL`/`region.getVertices` fold
measured here. So `testcPLQ`'s slow-bucket reds (G11, `TODO.md`) are not an orthogonal
legacy-pipeline concern: `rectMaximumIsTheConjugateOfTheWholeDomain` (G17, fixed 2026-08-27) and
`rectBiconjugateIsAConvexUnderestimator` (still times out, `>3600 s`) are running this exact fold
machinery on a different fixture, and a fix to one moves the other. Confirmed by re-running
`rectBiconjugateIsAConvexUnderestimator` this session: it produces the identical warning
signature (`isAlways:TruthUnknown` inside `region/getEdgeNos`, `functionNDomain/mergeL`) as this
profile, not a different failure.


### getVertices: affine x conic in closed form (2026-08-18)

An affine constraint meets a conic in at most two points, both given by the quadratic formula.
Substituting the line into the conic and clearing `b^2` gives `alpha*t1^2 + beta*t1 + gamma = 0`
with the coefficients in `getVertices`' header; a vertical line swaps the roles. Symbolic
throughout -- these become the region's VERTICES, and this session showed what one ULP of vertex
error costs.

Cumulative effect on one fold, with the `ptFeasible` filter:

    fold 2 (23 cells)      orig    +ptFeas   +getVert    net
    solve                   700        700        405   -42%
    subs                  11803       9445       9843   -17%
    isAlways              11793       8683       8688   -26%
    simplify               2256       2256       2256     0%
    getVertices           133.1 s    125.6 s    123.3 s   -7%
    maximumP              199.8 s    195.9 s    181.1 s   -9%
    fold wall               289 s      280 s      265 s   -8%

**TWO THINGS WENT WRONG ON THE WAY, and both are the same mistake.** The first version guarded
the closed form with a SAMPLED REFIT -- re-evaluate at two probe points, `simplify` the residual,
`isAlways` it -- to catch a cubic being truncated by a six-point quadratic basis. That guard cost
more than the `solve` calls it saved: `getVertices` went 125.6 s -> 136.9 s even though `solve`
fell 700 -> 405. The guard is a DEGREE question and `symbolicFunction.degreeNum` already answers
it for nothing.

The second: the affine and conic paths each built their own evaluation table, 3 + 6 = 9
substitutions per region where 6 suffice. Sharing one six-point table is what took the change
from break-even to -8%.

The lesson is the one this session keeps re-learning: **a safety check added to a hot path has to
be priced like anything else on that path.** The `ptFeasible` filter survived only because its
guard is a cheap numeric test.

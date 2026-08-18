# Decisions & Dead Ends

_The record of what has been tried and ruled out, so it is not tried again._

**Append-only.** Never delete an entry. If a decision is later overturned, strike
the heading through (`## ~~...~~`) and add a new entry saying what changed.

**Nothing here is an action item.** Near-term actions live in `TODO.md`. Dead ends,
reverted attempts, rejected approaches and refuted diagnoses live here.

**Read this before attacking a hard problem** — especially one that looks like it
has been attacked before.

Newest entries at the top.

> **Seeded 2026-08-13** by copying the dead-end passages out of `TODO.md`. `TODO.md`
> was left untouched because another session was working in this project at the
> time, so **the originals are still there and this file duplicates them**. When
> CCA2 is next idle, delete the copied passages from `TODO.md` — specifically the
> `## Retired hypotheses — do not re-try` section, the `TRIED, reverted` block under
> `MAJOR FINDING 2026-08-04`, the `ATTEMPTED TWICE AND REVERTED` text in the step-4
> item, the `Tried and REVERTED` paragraph in the obsolete `decideWinner` note, the
> `--- OBSOLETE ---` block, and the `SUPERSEDED` section at the end — and leave
> `TODO.md` as action items only.

---

## 2026-08-18 (last) — the parallelogram's "remaining 1%" was MY REFERENCE, not the code

Piece 9 of `f*` for `x*y` over `conv{(0,0),(2,0),(2.5,1),(0.5,1)}` is EXACT at all ten probe
points, one cell each -- `BAD 0 of 10`. The singular-quadratic overlap at `(1/2,1/4)` was the only
real defect, and `functionNDomain.singularEdgeCut` closed it.

**The three points I reported as a lingering ~1% over-claim were a broken reference.** The
brute-force sup was a grid over the piece, and for this piece the sup is attained AT A VERTEX --
`(1/4, 7/8)`, on the conic. A grid never lands there, so it reported `2.840439` where the true sup
is `2.875`, and the correct vertex cone read as over-claiming by 1%. Adding the exact vertices to
the candidate set makes all ten agree.

**The reusable part: a grid is not a reference for a sup that a VERTEX attains.** Every conjugate
in this codebase is a sup of an affine form, so its maximiser is at an extreme point far more
often than not -- put the region's own vertices in the reference before reading any disagreement
as a defect. This is the second time in one session that a measurement artefact was reported as a
code defect (the first is the correction entry about unmeasured Step 3 numbers).

## 2026-08-18 (later) — that revert was WRONG: the code was fine, the methods were in the wrong block

**Overturns the entry below, and the mistake is worth more than the fix.** The vertex-plus-arc-
tangency maximum was reverted after `testfunctionNDomain/testMerge` and
`cplqAdapterTest/twoTriangleSquareMaxMatchesNumericSup` went red, and it was written up as a
soundness failure with a guess about redundant conic facets. It was neither. The actual error:

    The class region has no Constant property or Static method named 'quadraticParts'.

`quadraticParts` and `lineMeetsConic` were appended next to `maxAffineOverRegion`, an INSTANCE
method, and then called as `region.quadraticParts(...)`. Moving them into the `methods (Static)`
block is the whole fix. `testMerge` has no assertions, so its "failure" was that exception --
which a single look at the report would have said, and which cost a revert instead.

**The lesson, and it is the reusable part: READ THE ERROR BEFORE THEORISING ABOUT THE MATH.** A
red test in this codebase is as likely to be a MATLAB scoping rule as a defect in the geometry,
and the report says which. The `unionIsExact` gate is soundness-critical, which is exactly why the
temptation was to assume the deep explanation.

**Measured with it back in**, same three folds of the A.4/A.5 quadrilateral:

    fold 3 cells                38 -> 36
    quadFacet_exactAnotInB      63 -> 41
    merges at fold 3             7 -> 9
    TOTAL                      828 -> 804 s

Green: fast 206/0, normal 11/0, testfunctionNDomain + regionTest 17/0, cplqAdapterTest 4/0.

## 2026-08-18 — the A.4/A.5 split is now the DEFAULT

Both objections recorded on 2026-08-16 are gone, measured:

  * `testcPLQ` with the split on: **8 passed / 0 failed in 2273 s**, against 4728 s and one ERROR.
    `testRectBiconj` is one of the eight -- that exception was a casualty of the double leaks and
    nothing in the test or the split was changed to fix it.
  * assembling `f*` for the general quadrilateral: 86 cells / 73 min -> 60 cells / 43 min.

2273 s against 1542 s off is a BUCKET question by the standing rule (2026-08-17), not a blocker:
`testcPLQ` is already in the slow bucket and finishes well inside its timeout. Verified with the
flip alone: fast 206/0, normal 11/0 in 230 s.

`plq_1p.appendTriangle` now gates on `CCA2_NO_A45_SPLIT` (opt OUT); `CCA2_A45_SPLIT` is still
honoured so the two quadrilateral tests that set it keep working.

## ~~2026-08-18 — REVERTED: an exact max over a region with a CURVED facet, for unionIsExact~~ (OVERTURNED, see below)

- **Tried:** replacing `region.impliedBy`'s LP-over-the-linear-relaxation with an exact maximum
  over the region itself (`maxAffineOverRegion` + `holdsOn`). The relaxation drops a conic facet,
  so it is sound but conservative exactly when that conic is what would have cut the violating
  part away -- `quadFacet_exactAnotInB`, 98 of fold 5's refusals and the largest NAMED gate left.
- **The argument, which still looks right:** a linear form on a compact set attains its max on the
  boundary; on a straight edge that is an endpoint (a vertex), on a conic arc it is an endpoint or
  a point where `grad h` is PARALLEL to the form. That parallel condition is affine, so it meets
  the conic in at most two points, in closed form. Vertices plus those points therefore cover
  every candidate.
- **Why it was reverted:** it FAILS. `testfunctionNDomain/testMerge` and
  `cplqAdapterTest/twoTriangleSquareMaxMatchesNumericSup` both go red with it in, and both go
  green the moment it comes out (11/0 with only the A.4/A.5 default flip live). Since the whole
  point of `unionIsExact` is soundness, a gate that admits merges the old one refused and then
  changes a VALUE is refuted, not debuggable by guesswork.
- **Where the hole most likely is, for whoever picks it up:** the argument needs the region's
  `vx`/`vy` to be every corner of its boundary, and needs the region to be the compact convex set
  it is assumed to be. A REDUNDANT conic facet is enough to break the first -- it makes `lin`
  false, so the refinement engages, while the vertex list still describes only the polyhedral
  part. Test that specific case first, on a region built by hand, before touching `merge` again.
- **Do not re-run this as-is.** The next attempt should be measured on
  `testfunctionNDomain/testMerge` FIRST, which is a unit-level merge test and takes seconds.

## 2026-08-17 (last) — item 1's root cause is the MESH VERTEX TYPE, and that is a design change

The 5 same-function pairs whose shared facet `merge` cannot see are two doubles of the same exact
number, one ULP apart. Traced to the source, and it is not `conjConvexOverPiece` (which now carries
whatever it is given exactly) and not `convEnvCPLQ` (which contains no `double(` at all):

    convEnvCPLQ on x*y over conv{(0,0),(3,3),(1,2)} returns envelope vertices
        (0,0)  (1.4142,1.4142)  (1,2)  (1.5,1.5)  (3,3)  (1.5858,1.5858)

`1.4142` is `sqrt(2)` and `1.5858` is `2 - sqrt(2)`, stored as DOUBLES. The reason is one line in
`RatPar.m`:

    V (:,2){mustBeNumeric} % nv x 2 matrix storing unique vertices

The mesh vertices are CONSTRAINED to be numeric, for the whole RatPar / RatPol / QuaPol / QuaPar
lattice, whose header records that every property lives there once and deliberately. So an exact
cevian point cannot survive being stored, no matter how exactly it was computed.

**This is therefore a design change, not a fix**, and it is why it was not attempted here: relaxing
`mustBeNumeric` on `V` (and consistently on `f`, `den`, `P`, `Ec`) touches every consumer of the
lattice, and the classes' own header explains why the properties are declared in exactly one place.
It also has a cheaper alternative worth pricing first: the A.4/A.5 path already keeps its geometry
exact WITHOUT going through a RatPol mesh (`splitTightTriangleSym` -> `plq_1p.triangulate`), and it
is the path the general quadrilateral uses. If the `ratPolToPlq` route is only needed for rational
envelopes, its inexactness costs cells and time but no correctness -- price that before changing the
lattice.

**What is measured, so the size of the prize is known:** on the tri case 5 of 31 same-function pairs
at fold 1 lose a real facet to this. On the A.4/A.5 quadrilateral, which does NOT go through the
mesh, the conjugate is exact (worst denominator 56).

## 2026-08-17 (latest) — testRectBiconj PASSES with the split on: the correctness blocker is gone

`testcPLQ/testRectBiconj` ERRORED with `CCA2_A45_SPLIT` set, and that -- not the runtime -- was
the stated reason the split could not be the default. Re-run against the exactness work
(`domain.mE`/`cE`, `region.limitOfFAtVertices`, `plq_1p.quadPartsOf`, `conjConvexOverPiece`):

    RESULT passed=1 failed=0 incomplete=0

So the exception was a casualty of the double leaks -- 145-digit coefficients and comparisons
`isAlways` could not decide -- and not an independent defect. Nothing was changed in that test or
in the split to achieve it.

**What remains for the default is therefore ONLY cost**, which by the standing rule
(2026-08-17, "DECIDED") is a bucket question rather than a blocker, unless it is on the way to not
finishing. Measure `testcPLQ` with the split on against its 1542 s off / 4728 s on, and note that
the machine is shared so a single timing decides nothing.

## 2026-08-17 (last) — the blow-up, finally COUNTED: two causes, and neither is the arithmetic

`.claude/step3adjacency.m` classifies every same-function pair of fold-1 cells three ways at once:
does `merge`'s own test see a shared facet; do the two carry the same hyperplane with opposite
orientation (numeric rows, so independent of how the constraint is written); and do they actually
MEET in a segment rather than a point. Run with Step 2 exact, on the control case:

    merge sees | same hyperplane | how they meet | pairs
    no         | no              | do not touch  |  5
    no         | yes             | a point only  |  7
    NO         | yes             | A SEGMENT     |  5     <-- facet detection FAILS
    yes        | yes             | do not touch  |  1
    yes        | yes             | a point only  | 14
    yes        | yes             | A SEGMENT     |  6     <-- reaches unionIsExact; 1 merged

**Only the 11 `segment` pairs are ones a merge could ever be right about.** 21 of the 38 pairs meet
at a POINT and must not merge -- the union of two cells touching at a corner is not convex, so
those refusals are correct and were never the problem. That alone retires "578 refusals" as a
number to reason from: most of them are right.

**The 11 split evenly into two DIFFERENT defects:**

1. **5 pairs share a real facet that `ineqs(i) == -ineqs(j)` does not find** -- with exact
   arithmetic, so this is not rounding. `symbolicFunction.eq` is `if (obj1.f == obj2.f)`, a
   STRUCTURAL comparison (its own comment says "change to isAlways"), so the same constraint
   written at a different positive SCALE does not match. `region.normalize1` is supposed to
   prevent that by dividing by `abs(coeffs(f,vars))(end)` -- check whether it picks the same
   term for both operands, which for two differently-written forms it need not.
2. **6 pairs are found and then refused by `unionIsExact`** -- and the fold-1 tally's
   `lin_exactCurvedTest = 6` matches exactly, so all six are `certifiesNonPositive` declining.
   It declines by design outside its hypothesis: a rational `h`, a non-convex quadratic, or a
   relaxation with no vertex. Find out WHICH of the three before extending it.

Both are small and both are now counted rather than guessed. That is the state to start from.

## 2026-08-17 (final for the session) — Step 2 is EXACT now, and Step 3 STILL does not merge

**Read this before spending another hour on exactness.** Three double leaks are fixed and the
whole conjugate is exact -- worst denominator across the six quadrilateral pieces went
`1.2e18 / 9.7e33 / 9.0e16 / 2.6e144 / 6.1e18 / 1.4e145` to `4 / 4 / 7 / 56 / 14 / 56`. On the
control case, Step 3 is essentially UNCHANGED:

              before exactness            after exactness
    FOLD 1    20 cells, 7 fns, 1 merge    21 cells, 7 fns, 1 merge
    FOLD 2    37 cells, 11 fns, 2 merges  37 cells, 10 fns, 2 merges
    FOLD 3    57 cells, 10 fns, 4 merges  53 cells,  9 fns, 4 merges
    fold 3 refusals: noSharedFacet 475 of 511.

**So the ULP-apart doubles were REAL and were NOT the cause.** The entry below identified two
constraints one ULP apart on `4 - 2*sqrt(2)` and concluded the facet test was blinded by
arithmetic. The arithmetic is now exact and the test still does not find the facets. That
conclusion is retracted; what remains true is that the leaks were genuine defects (one of them
produced a y-intercept of `-9.06e-72` where the exact answer is `0`) and are worth having fixed on
their own.

**What is NOT yet explained, and is where to start next.** 15 of 31 same-function pairs at fold 1
carry the same hyperplane with OPPOSITE orientation, measured off `linearForm`'s numeric rows,
while `merge`'s `ineqs(i) == -ineqs(j)` sees nothing. That measurement predates the exactness work
and should be REPEATED first -- if it still holds with exact numbers, the remaining candidates are:

  * `symbolicFunction.eq` is `if (obj1.f == obj2.f)`, a STRUCTURAL test (its own comment says
    "change to isAlways"), so two forms of the same constraint differing by a positive SCALE do
    not match. `region.normalize1` divides by `abs(coeffs(f,vars))(end)`, which is the highest
    term in MATLAB's ordering -- check that it picks the same term for both operands.
  * sharing a hyperplane is necessary for adjacency but NOT sufficient: the two cells may lie on
    the same line yet not touch. The probe above does not distinguish those, and should.

Do that measurement before writing any code. The three hypotheses this session tested were each
plausible, each partly right about a real defect, and each wrong about the cell count.

## 2026-08-17 (latest, and it corrects the entry below) — the facet test cannot match two doubles of the same number

**The entry below is WRONG about the control case, and the way it is wrong is worth keeping.** It
concluded that exactness was not the lever because the "all-rational" case blew up too. That case
is not rational: `conjConvexOverPiece` converts `Q, L, c` and the piece's vertices to DOUBLE by
design (its own lines 59 and 73), whatever it is handed. There was no control.

**What the refusals actually are.** With merge's two heuristics removed, every refusal on that case
became `noSharedFacet` -- 578 of fold 3's 612. A direct check of the fold-1 cells says that is a
FAILING TEST, not geometry: of 31 same-function pairs, **15 carry the same hyperplane with opposite
orientation and `ineqs(i) == -ineqs(j)` does not see it** (12 are seen by both, 4 share nothing).

**And here is why it cannot see it.** Cells 8 and 9 meet along one facet, and carry

    s_2 - 659536895553805/562949953421312       = 5276295164430440/4503599627370496
    s_2 - 5276295164430439/4503599627370496

-- two doubles of the same exact number, `4 - 2*sqrt(2)`, **one ULP apart**. No comparison can
identify them: not the structural `==` this code uses, and not `isAlways` either, because they are
genuinely different rationals. The facet is real and the arithmetic has destroyed the evidence.

**So the chain is one chain, and exactness is the lever after all:** a double enters Step 2 ->
the same quantity acquires two different values in two different cells -> merge cannot match the
shared facet -> nothing merges -> the cell count grows without bound. `domain.mE`/`cE` was one
source of doubles and is fixed; `conjConvexOverPiece` is the other, and it is the one left.

**What the heuristic removal was worth, honestly.** Nothing yet, on the measurements: cell counts
and successful merges are unchanged (20/37/57 cells, 1/2/4 merges). It is still the right code --
two heuristics replaced by `region.certifiesNonPositive`, a sound closed-form certificate -- and it
is what will do the work once the facets can be found again, since `unionIsExact` then becomes the
gate that actually decides. But it fixed nothing on its own and should not be described as if it
had.

## 2026-08-17 (latest) — The CONTROL case: the doubles are real but they are NOT the blow-up

**This corrects the entry below.** `domain.mE`/`cE` really were double arrays and that really was a
defect -- an exact slope arrived as `0.6` and an exact zero y-intercept as `-9.06e-72`, which is a
wrong value, not merely an expensive one, and fixing it took the quadrilateral's worst conjugate
coefficient from `1e144` to `1e33`. But it changed Step 3 **not at all**: cell counts and the merge
tally came back byte-identical.

**The control that decides it.** `.claude/step3cost.m` with `CCA2_STEP3_CASE=tri` runs `x*y` over
`conv{(0,0),(3,3),(1,2)}` through `convEnvCPLQ` + `ratPolToPlq` -- four pieces, ALL RATIONAL, no
A.4/A.5 split, no surds and no doubles anywhere:

    FOLD 1: paired=17 -> cells=20 distinctF= 7   merge: okLinear=1  noSharedFacet=14 quadCutsOther=14 quadMismatch=34
    FOLD 2: paired=36 -> cells=37 distinctF=11   merge: okLinear=2  noSharedFacet=54 quadCutsOther=26 quadMismatch=58
    FOLD 3: paired=60 -> cells=57 distinctF=10   merge: okLinear=4  noSharedFacet=266 quadCutsOther=50 quadMismatch=272 lin_exactCurvedTest=19

57 cells for 10 distinct functions, and **4 successful merges out of 612 attempts**. Same blow-up,
clean numbers. So exactness is not the lever: **`region.merge`'s own gates are**, and they are what
to fix.

**Ranked by what actually fires, measured:**

1. **`quadMismatch`** -- 272 of fold 3's 612. When BOTH regions carry a quadratic constraint and
   they do not share one as a facet, merge demands that EVERY quadratic of A equal EVERY quadratic
   of B, as a cross product, and refuses otherwise. Two adjacent cells each carrying a different
   parabolic arc elsewhere have a perfectly convex union; this refuses them outright.
2. **`noSharedFacet`** -- 266. With clean numbers a large part of this is honest: 10 functions over
   57 cells means groups of ~6, most of whose pairs are genuinely not adjacent. Check before
   attacking it.
3. **`quadCutsOther`** -- 50. The other heuristic: refuse if one region's quadratic meets the other
   anywhere but at a vertex.
4. **`exactCurvedTest`** -- 19. The SOUND gate: `unionIsExact` cannot certify `A subset B'` when a
   constraint it must test is non-affine, so it refuses.

**The fix these point to, and it is one fix.** 1 and 3 are heuristics standing in for the
certificate 4 cannot supply. Give `unionIsExact` that certificate -- "is `max h <= 0` over the other
region" for a non-affine `h`, which for a CONVEX quadratic over a polyhedron is decided exactly by
its vertices plus its recession directions, the region's own curved facets being droppable because
the linear relaxation is a superset -- and then the two heuristics can go, because `unionIsExact` is
the exact criterion (`M = A' n B'` equals `A u B` iff `A subset B'` and `B subset A'`) and refusing
only ever costs compactness, never correctness.

## 2026-08-17 (later) — Step 3's blow-up MEASURED: merge stops working entirely, and doubles are why

Two measurements, on `x*y` over `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}` with `CCA2_A45_SPLIT` on.
Both are reproducible with `.claude/step3cost.m`.

**Merge succeeds ZERO times from fold 2 on.** 190 attempts, no success, while the 29 surviving
cells carry only **7 distinct functions** -- and `distinctF` is 7 at fold 1 too, so the answer
never needs more than a handful of cells. The cell sequence 5, 14, 29, 45, 70, 86 is not cells
being created, it is merges being refused. Merging is also half the cost (182 s of fold 2's
223 s). Refusals by reason: `noSharedFacet` 70, `quadMismatch` 50, `quadCutsOther` 34,
`exactCurvedTest` 30, `exactAnotInB`/`exactBnotInA` 6.

**And the reason merge cannot decide anything is that Step 2 reintroduces DOUBLES.** Worst
denominator by stage: split sub-triangle domains **12**, Step 1 envelope faces **20**, Step 2
`conjugates` **1.2e18 / 9.7e33 / 2.6e144 / 1.4e145**, `maximumConjugate` unchanged. The split is
exact and Step 1 is exact; `plq_1p.conjugate` is where it goes wrong. Piece 4's conjugate
constraints carry `7307585874000779/9007199254740992` (denominator **2^53**) alongside the exact
`30^(1/2)/12 - 15^(1/2)/6 + 3/4` -- the same kind of quantity in two forms, one exact and one its
own double -- and one coefficient is 97 digits long.

**Why this reframes the work.** `merge` finds a shared facet by `ineqs(i) == -ineqs(j)` and
compares quadratics by `~=`. Neither test can succeed when one side is exact and the other is a
double of it, which is exactly what `noSharedFacet=70` and `quadMismatch=50` look like. So the
merge gates are not established as the defect yet: the numbers they are being asked to compare
are. **Fix the Step 2 leak first, then re-measure the tally**, and only judge `unionIsExact` and
the two quadratic pre-checks against clean numbers.

This is attempt 3's pathology surviving downstream of the fix that removed it from the split --
the same lesson the quadrilateral entry records, in a place nobody had looked.

## 2026-08-17 — DECIDED: the split stays opt-in until Step 3's cost is fixed, and the principle behind it

**The user's call, and the reasoning is worth keeping because it settles a whole class of future
questions.**

**Chosen: option (a)** — leave `CCA2_A45_SPLIT` opt-in and attack Step 3's cell blow-up first,
because that is also what the parallelogram's residual 4% needs. Options (b) (diagnose
`testRectBiconj` first) and (c) (flip now) were declined.

**The principle the user attached to it, which is general and outranks the flag:**

> Ultimately all computations have to be CORRECT even if they take a long time. In that case,
> split the unit tests between fast, medium and slow. Correctness is more important than speed —
> **unless** it is so slow that the user does not wait, because getting a timeout is not helpful
> either.

So the reason the split is still opt-in is **not** that it costs 1542 s → 4728 s. A correct path
that is slow gets its test moved down a bucket (`.claude/suite.sh` already has fast / normal /
slow) and keeps running. The only live objection is the **`testcPLQ/testRectBiconj` EXCEPTION**,
which is a correctness question, plus the risk that 4728 s is on the way to "does not finish",
which is the one failure mode the principle does not tolerate either. Both of those are what
Step 3's cost work addresses.

**What this rules out for good:** proposing that a proven-correct path be reverted or gated purely
because it is slower. That is not an available answer here.

## 2026-08-16 (last) — Making the A.4/A.5 split the DEFAULT: deliberately not done

- **Tried:** shipping the split on by default, since it fixes a documented crash and a documented
  wrong answer (`plq_1p.appendTriangle`, `splitTightTriangleSym`).
- **Why it was not kept:** with it on, `testcPLQ` runs 4728 s against 1542 s off (its historical
  time is 1427 s), **and `testcPLQ/testRectBiconj` ERRORS**. That test has no assertions — it runs
  `triangulate`, `maximum`, `biconjugateF` and nothing else — so the failure is an exception, not a
  wrong value. Undiagnosed. Turning the split on therefore trades a documented, LOUD failure on one
  domain shape for a new one on another, which is not a trade to make silently.
- **Before retrying:** fix Step 3's cell blow-up (the cost) and diagnose `testRectBiconj`'s
  exception (the correctness question), in that order. Both are in `TODO.md` with their numbers.
  The switch is `CCA2_A45_SPLIT`; the two quadrilateral tests set it themselves, so the fix stays
  exercised either way.
- **Evidence:** the runs behind those two timings are the `testcPLQ`-only runs of 2026-08-16,
  uncontended in both directions; commit `ba3457d`.

## 2026-08-16 (later) — The general quadrilateral, FIXED on the fourth attempt. What the three failures were each worth.

The method the user prescribed is what worked, and it is worth stating on its own: **do it
symbolically first, and only then reach for explicit formulas.** Attempt 3 failed for exactly the
reason attempt 4 succeeded — it computed the same geometry in double precision and let `sym` turn
the result into `2^53` denominators.

**The fix.** `splitTightTriangleSym` splits a triangle into sub-triangles on each of which cPLQ's
own closed form for THAT sub-triangle's convex-edge count IS the convex envelope, and
`plq_1p.triangulate` emits them as PIECES. A.4 gives one two-convex-edge half (its own form
unchanged and now tight) plus one one-convex-edge half (A.3's rational form, which the `nCE == 1`
branch derives analytically); A.5 gives two two-convex-edge halves that recurse into A.4. Exact
symbolic arithmetic throughout, so the cevian foot — irrational in general — stays a compact surd
(`5/2 − sqrt(5)/2`, `3/2 − 3·sqrt(5)/10`).

**Two things learned while writing it that are worth keeping.**

- **`twoEdgeQuadPlain`'s ± branch search is vestigial: it always returns `s = +1`.** Substituting
  `y = mh·x + qh` into the `s`-form gives `x²` coefficient `mh·denom/denom = mh`, `x` coefficient
  `qh·denom/denom = qh` and constant `0` — for BOTH signs. So the form touches `x·y` along edge `h`
  either way and the touching test cannot separate them; `s = −1` is merely undefined when
  `mh = mw`. cPLQ hard-codes `s = +1` and is right to.
- **The "no split needed" test has a one-line geometric reading.** With `s = +1` the curvature of
  `q1` along the weak edge `AB` is `(sqrt(mh·mw)·dx + dy)² / (sqrt(mh) + sqrt(mw))²` — a perfect
  square — so it vanishes exactly when `AB` is PARALLEL to the cevian direction `−sqrt(mh·mw)`.
  That is the same degeneracy `convEnvCPLQ`'s header describes as "both candidate cevians
  degenerate", arrived at from the other side.

**What the split costs, and where the cost actually is.** The no-split path — every input cPLQ
itself ever had — is 20 ms, which is why the fast bucket did not move. The split paths are 0.3 s
(A.4) and 1.2 s (A.5). What is expensive is what comes AFTER: six pieces instead of two, and
Step 3's cross-piece maximum then takes **73 minutes**, with the cell count running 5, 14, 29, 45,
70, 86. The answer is exact at both levels — 10 of 10 on the per-piece max, 8 of 8 on the assembled
one — so this is a scaling problem in `maxOfList`, not a defect in the split. It has its own
`TODO.md` item.

**And the real reason it is expensive is the ALGEBRAIC DEGREE, not the piece count — which is worth
knowing before optimising the wrong thing.** A.4's cevian foot is irrational, so a split
sub-triangle has SURD coordinates and every symbolic operation downstream works in a quadratic
extension instead of the rationals. Measured on `testcPLQ`, whose domains are general polygons
carrying `x·y`: **1542 s with the split off (matching its historical 1427 s), over 3100 s with it on
and still unfinished when stopped**, uncontended both times — while only two of its six domains
gain a piece at all (2 → 3 and 1 → 2). Two extra pieces cannot explain a 2×+ slowdown; surds can.
`CCA2_NO_A45_SPLIT` therefore exists as an opt-out, and correctness keeps the default.

## 2026-08-16 — The parallelogram's `emptyResult`: TWO defects, both fixed, and two more measured and not taken

Traced from `QuaParCPLQ:conj:emptyResult` down to a single piece, then to a unit reproducer that
runs in about a minute. Piece 9 of `f*` for `x·y` over `conv{(0,0),(2,0),(2.5,1),(0.5,1)}` went from
**6 of 10 probe points wrong or uncovered to 2**, against a brute-force sup over the piece.

**How it was found, which is the reusable part.** `f**` of a bounded domain is finite exactly ON
that domain, and `f**` is a MAX, so EVERY per-piece conjugate must be finite there. Evaluating all
12 groups at six points inside the parallelogram showed three with holes — 9, 11 and 12 — while the
other nine covered everything. That turned "the max comes out empty" into "these three pieces are
wrong", in one cheap measurement. The accumulator and group 11 covered DISJOINT halves of the
parallelogram, so the empty intersection was the honest answer to a wrong question.

**Defect 1 — `region.simplifyUnboundedRegion` declares any region with no finite vertex empty.** A
half-plane has none; so do a slab and the whole plane. And a half-plane is exactly what a TANGENT
vertex produces: the cone at a vertex is built from the two edges meeting there, and when those are
tangent — an arc and its chord touching, which is how a curvilinear piece ends — both half-planes
are the SAME one. Here the cone at `(1/4,7/8)` is `{2x/3 + y ≥ 4/3}` carrying `x/4 + 7y/8 − 1/2`,
and it was being deleted. Now refuted by a WITNESS: a feasible point certifies non-emptiness, while
failing to find one proves nothing, so the old verdict stands whenever no witness turns up.

**Defect 2 — the edge list, again, in its other form.** The piece has 3 vertices and 3 genuine
edges plus a conic touching one vertex: 4 constraints for 3 vertices, so `size(ineqs,2) == nv` calls
a BOUNDED region unbounded, `endNv` comes out `nv−1`, and it is built one edge cell short. Same
root cause as bug 1, without the lens collision — so `conjugateOfPiecePoly` now derives the edge
list for any bounded piece the count mislabels, with the count-matching branches keeping precedence.

**Measured and NOT taken, both at the `isQuad` chord rewrite.** `biconjugateTest`'s own comment
asked for new evidence before re-trying this; here it is, and it says leave it alone.
- Chording **the vertices the conic actually touches** (`region.vertexOfEdge`) instead of
  `vx(1),vx(2)`: makes this piece WORSE, 2 wrong of 10 → 3. The arc's cell grows to swallow a
  region the vertex cone was answering correctly.
- **Skipping the rewrite entirely**: changes this piece not at all, 2 of 10 either way.

**What the last 2 of 10 are.** The chord's edge cell and the arc's claim the SAME region, and the
chord's is checked first and is wrong there (0.0176 and −0.0138 against 0.0333 and 0.0363; the
arc's cell has the right value at both). `f` is a SINGULAR convex quadratic on this piece —
constant along one whole edge, `∇f = 0` at two of the three vertices — so `functionNDomain.getInterior`,
which eliminates `s` between `x = ∂₁f` and `y = ∂₂f`, returns the gradient map's image LINE rather
than a curve separating the two cells. That is the next defect, and it is not in the chord rewrite.

## 2026-08-16 — The A.4/A.5 split as a DOMAIN split: written, measured, REVERTED. The blocker is arithmetic, not structure.

Third attempt at the general quadrilateral, and the first one that failed for a reason none of the
previous analysis predicted. **The structure was right; the ARITHMETIC is what kills it.**

**What was built.** `plq_1p.triangulate` split each indefinite triangle into the A.4/A.5
sub-triangles before emitting pieces, so every sub-piece lands on a path Step 2 already has — a
2-convex-edge sub-triangle where cPLQ's own closed form IS tight, plus a 1-convex-edge one whose
A.3 rational envelope the `nCE == 1` branch derives analytically. The sub-triangles came from
`convEnvCPLQ` itself (its FACES are the split, in x coordinates, recursion done), so nothing had to
be moved out of that file. It worked at Step 1 exactly as designed.

**Why it was reverted: it turns the crash into a HANG.** The quadrilateral's first conjugate ran
for over 45 minutes without finishing, against ~3 minutes before, and had stopped producing output
— stuck inside a symbolic call, with 3.8 MB of `isAlways:TruthUnknown` warnings behind it carrying
coefficients like

    −609298613085773108668343859 / 14507109835375549828038656 .

A hang is worse than a crash here: this repository's own suite driver exists because "a wedged
suite no longer costs the whole run."

**The cause, and it is inherent to the approach.** [COAP] A.4's cevian has slope `−sqrt(mh·mw)`, so
its foot is IRRATIONAL in general. `convEnvCPLQ` is double precision throughout, and `sym` of a
double is EXACT — a denominator near `2^53`. Snapping the new vertices to the simplest rational
within `1e-10` (which the reverted code did, and which is the same numeric→symbolic hygiene
`ratPolToPlq` already applies to coefficients) bounds the VERTEX denominators but not the
downstream ones: the conjugate is a rational function of those coordinates, so a few squarings and
additions carry `1e5` denominators to `1e25`, and MuPAD's `isAlways` then cannot decide anything.

**So the split has to be carried SYMBOLICALLY, not snapped.** The cevian foot is an exact algebraic
number — intersect the line through the far vertex of slope `−sqrt(mh·mw)` with the other convex
edge — and `sqrt` is something the symbolic layer handles natively and keeps small. That means
implementing A.4 (and A.5's smooth-fit line) in symbolic form beside the existing double version,
rather than calling `convEnvCPLQ` for the geometry. It is a bounded piece of work — the formulas
are two line intersections and a curvature test — but it is new code, not wiring.

**What NOT to do next.** Do not re-try any of: installing A.4/A.5 faces as ENVELOPE faces (Step 2
has no rational-envelope branch — 2026-08-15 entry), calling `convEnvCPLQ` for the domains and
rounding (this entry), or leaving `nCE == 2` alone in the hope the error is elsewhere (it is the
envelope itself — 2026-08-15 entry, with the −0.2835 measurement).

## 2026-08-15 — cPLQ's `nCE == 2` envelope is NOT the envelope, and the obvious fix cannot work

**The defect, derived and then measured.** `plq_1p.convexEnvelope1`'s `nCE == 2` branch applies
[COAP]'s single-quadratic form to the WHOLE triangle. That form touches `x·y` along both convex
edges and is a valid convex MINORANT — but Appendix A.4 shows it is tight only over a sub-region,
which is exactly what `convEnvCPLQ`'s `splitTwoConvexEdges` tests for and splits on. This branch
never tests.

Measured on `conv{(2.5,1.5),(2,0),(0,0)}` carrying `x·y`, the piece `triangulate` produces from the
test quadrilateral. cPLQ returns

    0.954915·y − 0.572949·x + 0.427051·x·y + 0.286475·x² + 0.159153·y²

whose minimum over the triangle is **−0.2835**, at `(1,0)`. On that triangle `x ≥ 0` and `y ≥ 0`, so
`x·y ≥ 0` and the affine minorant `0` is admissible: **the true envelope is ≥ 0 everywhere**, and
this is strictly below it. A too-small envelope gives a too-large conjugate, and that is the whole
error — `f*(0,0) = 0.28647` for a truth of `0`, `f*(0.5,1) = 1.00464` for `1`. `convEnvCPLQ` on the
same triangle returns **2 faces** — it does apply the split — with minimum `0`.

The `0.28647` is not a coincidence worth chasing: it is the form's own `x²` coefficient,
`m_h·m_w/(√m_h+√m_w)²` with the two convex slopes `3` and `0.6`.

**AND ROUTING STEP 1 THROUGH `convEnvCPLQ` DOES NOT FIX IT.** Tried, measured, reverted:
`plq_1p.conjugateFunction`'s `nCE == 2` branch reads its envelope's coefficients with
`coeffs(envelope.f.f, vars)` and matches monomials. `convEnvCPLQ`'s A.4/A.5 faces are **RATIONAL**,
so it raises `symbolic:coeffs:NotAPolynomial` outright. **cPLQ's Step 2 has no rational-envelope
branch at all**, and the dispatch keys on the PIECE's `nCE` rather than on the envelope face in
hand — which the routine's own header already complains about for a different reason.

**So the split belongs in the DOMAIN, not in the envelope.** The route that already works for
rational faces is `conjCPLQ`'s `conjEnvelopeViaCPLQ`: it hands each rational face to cPLQ as its own
PIECE via `ratPolToPlq`, and lets `plq.maximum` take the max. The same idea applies here — have
`plq_1p.triangulate` split a 2- or 3-convex-edge triangle into the A.4/A.5 sub-triangles and emit
each as a piece, recursing while `splitTwoConvexEdges` still reports `needsSplit`. Every sub-piece
is then a triangle on which cPLQ's own closed form IS tight, so Step 2 is untouched and every
envelope stays polynomial.

**What that costs, and why it was not done unattended.** `splitTwoConvexEdges`, `splitThreeConvex`
and their helpers (`classifyConvexEdges`, `solveTriangleBF`, `envelopeFromClassified`,
`bilinearFrame`, …) are all file-local to `convEnvCPLQ.m`, so exposing them means moving a
connected web of functions out of a well-tested file — and `triangulate` feeds every Case C result,
so the blast radius is the whole symbolic pipeline. That is a design change with a full
re-verification behind it, not a fix.

## 2026-08-15 — The general quadrilateral's `nCE == 3` wiring: written, measured, REVERTED

**Do not re-land this until cPLQ's Step 2 is fixed.** The `nCE == 3` branch is not the reason the
general convex quadrilateral fails; it is only the reason it fails LOUDLY.

**What was written, and it is correct as far as it goes.** In `plq_1p.convexEnvelope1`, for
`nCE == 3`, build the triangle as a one-face `QuaPol` carrying `x·y` — safe, because reaching that
line with an indefinite `q` means `plq_1p.isCanonicalXY` held, so `q` is EXACTLY `x·y` with no
linear or constant part — call `convEnvCPLQ` (which has [COAP] A.5's `splitThreeConvex`, and
recurses so each half also gets A.4's tightness check), convert with `ratPolToPlq`, and install the
faces as envelope pieces. `plq_1p.conjugate` already loops over envelope pieces and accumulates
`conjfia`, so several per input piece is the normal shape, and Step 3's max over the results is
correct because a sup over a union is the max of the sups. **Measured: 4 envelope faces for the
offending triangle, two quadratic and two rational, all ≤ `x·y`, and `conj` stops raising
`MATLAB:badsubscript`.**

**Why it was reverted.** With the branch in, `f*` of the quadrilateral is WRONG at 4 of 8 probe
points — over-claiming `0.28647` at `(0,0)` where the truth is `0`, `1.00464` at `(0.5,1)` where it
is `1`, and uncovered at `(-1,0.5)` and `(-2,-2)`. **A silent wrong answer is worse than a loud
refusal**, so the crash stays until the thing underneath it is fixed.

**Where the fault actually is, separated by measurement rather than argued.** `triangulate` splits
the test quadrilateral into piece 1 `[2.5 1.5; 2 0; 0 0]` (`nE = 2`) and piece 2
`[2.5 1.5; 0 0; 0.5 1]` (`nE = 3`). Evaluating each piece's own Step 2 conjugate inside Case C:

- **piece 2 gets its 4 envelope faces and Step 2 returns ZERO conjugate cells for it.** The new
  envelope is computed and then discarded, so the wiring buys nothing today.
- **piece 1 — untouched by any of this — carries the whole error**: 6 cells, and every wrong value
  above is its.

**The measurement that localises it, and the trap it avoids.** That same `nE = 2` triangle
conjugated ON ITS OWN via `QuaPol.conj` is exact at 7 of 7 probes — which looks like an alibi and is
not: a single bounded triangle takes the NUMERIC route (`conjBoundedPolygon`), not cPLQ. **The
numeric Step 2 is right on this input and the vendored symbolic one is not.** Checking a suspect
piece "on its own" through the public API can silently change which implementation runs; evaluate
`p.pieces(k).maxConjugate` inside the pipeline instead. `assertStep3MatchesPieces` correctly does
not fire here, because Step 3 does agree with Step 2 — the gate is doing its job, and the fault is
one stage earlier than it looks at.

## 2026-08-15 (last) — BUG 1 FIXED: the edge list, and why the four earlier attempts could not have worked

**What the refactor turned out to be, and it is smaller than the note below predicted.** Four
attempts had all tried to make the LENS fit one of the two slot conventions — free a slot, park a
claim, drop a constraint. None can work, and the reason is worth stating plainly: *the conventions
cannot express a lens at all.* Both say "edge `j` is at `ineqs(j)`" or "at `ineqs(j+1)`", i.e. they
identify an edge by a VERTEX INDEX — and a lens's two edges join the SAME pair of vertices, so no
assignment of slots to constraints can name them apart. The question is not which constraint gets
which slot; it is that slots are the wrong addressing scheme for this shape.

So `functionNDomain.edgeIndexList` derives `eIdx(j)` — the constraint bounding the edge from vertex
`j` to vertex `j+1` — from the geometry, and the three edge-indexed readers take it as an argument.
Two things made this much less invasive than the earlier note predicted:

- **Both loops are in `conjugateOfPiecePoly` itself**, so the list is built in the first and
  consumed in the second with no class field, no signature change on the pieces, nothing to carry.
- **`getNormalConeEdgeQ` and `getNormalConeEdgeQ3` are the SAME routine.** Both build the
  perpendicular to the EDGE'S OWN constraint at each endpoint, oriented by the other endpoint; they
  differ only in which slot they believe that constraint occupies. Given the list they collapse to
  one routine (`getNormalConeEdgeQE`, over `region.coneNormalAt`). Likewise `getSubdiffVertexT2` and
  `getSubdiffVertexT2Q` are identical on these inputs — `T2Q`'s extra third column is only ever
  read when `NCE(:,3)` is non-zero, and neither edge routine sets it. **Three of the "four routines
  that move together" were two routines wearing four names.**

**How the list is decided, and the trap in the obvious version.** A constraint bounds the edge
between two consecutive vertices when both vertices lie on it AND its own curve between them stays
in the region. The second half is not a refinement: on the pipeline's own lens BOTH conics pass
through both vertices, and the redundant one meets the region nowhere else. The first version of
this preferred "the slot today's convention would use" before filtering on that, which handed edge
2 the redundant conic — a constraint bounding nothing. **Filter first, prefer second.** The
preference is what keeps every piece that works today on exactly the indices it has.

**Scope, and why it is safe.** The new path is entered only when two constraints that each bound a
genuine edge STILL share an edge number after `spreadCollidingEdges` has moved everything it can
(`edgesStillCollide`) — the lens signature, reached by nothing else. No constraint is dropped: the
two unsound drops recorded below stay ruled out, and the vertex cones' own feasibility probes
(`ptFeasible`, inside `getNormalConeVertexQ`) need the full constraint set.

**Measured.** Both half-lenses conjugate to 3 cells — two vertex cones plus the arc; the chord's
cell is a ray and drops out — exact against a brute-force sup over the lens at 12 points, 0 wrong.
The three the old code got wrong: `(0,0)` and `(-1,0.5)`, `+inf` where the truth is `0`, and
`(2,-1)`, `0` where the truth is `1/2`.

**Two defects were found by re-reading the finished code, not by a test**, and both are recorded
because the second is the same shape as three other bugs in this file. (i) The two halves of
`getNormalConeVertexQ` were given `eIdx(j)` and `eIdx(j+1)` where the routine's own probe points say
they mean `eIdx(j-1)` and `eIdx(j)` — invisible at `nv = 2`, where predecessor and successor
coincide, and wrong for any larger region. **When a routine's index convention is unclear, read its
PROBE POINTS: they are the geometry it is actually using.** (ii) `edgeIndexList` built `nv-1`
entries for an unbounded region while its consumer walks all `nv` vertices — an index used outside
the guard that produced it, for the fourth time in this codebase.

## 2026-08-15 — BUGS 3 and 4 FIXED: the curved envelope over an unbounded face

Recorded because the DERIVATIONS are the reusable part, and because the shape of the argument is
different from the affine cases already in `convEnvUnbounded`'s header.

- **Wedge, one flat ray `d1` and one convex ray `d2`:** `co q` is `q` with its CROSS TERM deleted,
  `q(v) + α·g1 + β·g2 + β²·A22/2`. The derivation is the file's own method — parametrise the affine
  minorants by their gap parameters and minimise the gap at an arbitrary target — but with one
  crucial difference: **the minimiser DEPENDS on the target point** (through `β0`), which is
  exactly why this envelope is not affine and why the affine argument does not extend. `α → ∞`
  forces `A12 ≥ 0`; a negative `A12` makes `d1 + t·d2` a recession direction of negative curvature,
  so the envelope is `−inf`, and that is now reported rather than answered.
- **Half-strip convex along the ray:** it separates **only when the base edge is Q-ORTHOGONAL to
  the ray**. Then `co q = q(v1) + s·(q(v2) − q(v1)) + t·⟨∇q(v1), d⟩ + t²·(d'Q d)/2`. `w'Q d ≠ 0` is
  refused loudly, with the value in the message — the two directions interact and the minimising
  minorant moves with the target in both coordinates, so there is no separable answer.
- **Checked against the fixtures, not fitted to them:** `co(x·y + I_K) = y²` on `{0 ≤ y ≤ x}` and
  `co(−x²+y²) = −x + y²` on `{0 ≤ x ≤ 1, y ≥ 0}` are what the general formulas produce, and both
  match what `unboundedFaceTest` derives by hand. `unboundedFaceTest` 18 / 0, from 16 / 2.
- **What is still not covered:** a wedge with BOTH rays convex, and a half-strip whose base edge is
  not Q-orthogonal to its ray. Both are refused, not approximated.

## 2026-08-15 — BUG 2 FIXED: a tangent built where the gradient vanishes

**A vanishing gradient at a cone's apex is a recurring failure mode in `region.m`, and this is the
second routine to fall into it.** `SUPPORT_MATRIX.md` §8.2(e) records the first:
`simplifyUnboundedRegion` decided emptiness from probe directions built out of constraint SLOPES at
a vertex where the split conic's gradient vanishes, and `region.witnessAwayFrom` was written for
it. Same input, same trap, different routine.

- **What it was:** `removeTangent` takes a quadratic constraint active at a vertex, builds the
  TANGENT LINE to it there, and deletes any constraint equal to that tangent as redundant. At the
  APEX OF A CONE the quadratic's gradient VANISHES — there is no tangent line, every direction is
  tangent — and whatever it computes is meaningless. It then deletes a constraint that matches.
- **Measured, on the 4-cone fan:** the assembled cell carrying `s1²/4 + s2²/2` lost its `−s1 ≤ 0`,
  keeping only `{s2 ≤ 0, s2²/2 − s1² ≤ 0, s1² − 2s2² ≤ 0}` — two constraints **blind to the sign of
  `s1`** — so the region became symmetric under `s1 → −s1` and claimed the mirror wedge.
  `f*(-3,-2.4)` came back `5.130` for a truth of `4.500`. It is now `4.5`; `conjCPLQTest` 25 / 0.
- **Fix:** refuse to conclude anything from a vanishing gradient.
- **Cleared on the way, so nobody re-checks them:** `region.redundantSubset` certifies nothing
  redundant on that constraint set (correct), and `simplifyUnboundedRegion` leaves the constraint
  alone. Step 2 is right too — each primal piece's own conjugate has the correct quadrant
  constraints, and the per-piece max at that point is the truth.
- **How it was found, and this is the transferable part:** by bisecting the pipeline rather than
  reading it. Dump the cells carrying the offending quadratic immediately before and after
  `mergeL`; the constraint is present before and gone after. Then feed that exact region to each
  simplification `mergeL` applies, in turn, and see which one drops it. Three routines, one run.

## 2026-08-15 (later) — BUG 1: three defects fixed, and one attempted fix that is UNSOUND

**Fixed, and each was necessary before the next was visible.** The lever that made this tractable
was a unit-level reproducer: the half-lens `{(s1+s2)² ≤ 4·s2, s2 ≤ s1}` carrying `s1`, built by
hand as a `region` and conjugated against a brute-force sup over its own boundary. No pipeline,
seconds per run, and it went from **2 identical wrong cells to 3 cells exact at all 10 probe
points**. Build that first next time; the pipeline runs took 10–40 minutes each.

1. `getEdgeNosInf` numbers an edge by one of its endpoint VERTICES, so a LENS — two edges joining
   the same pair — gets one number for both, and the last-write-wins scatter destroys one.
2. `getNormalConeVertexQ` (the routine that builds a vertex cone from the CONSTRAINT's own tangent
   rather than the chord) indexed its second constraint as `j+1` unwrapped, so it raised
   `badsubscript` on any BOUNDED region — which is why its only caller gated it behind a
   constraint COUNT and sent every bounded region to the polyhedral routine instead. Wrapped
   cyclically; identical to `j+1` for the unbounded layout, so nothing that worked changes.
3. `biconj` handed its second conjugation the curved MESH `conj` has returned since 2026-08-13,
   and `quaPolToPlq` refuses a curved domain. It now asks `conjCPLQ` for the symbolic form.

**UNSOUND, and reverted: freeing an edge slot by dropping the constraint holding it.**

- **Tried:** the lens's two edges need slots 1 and 2, which are held by constraints with a single
  vertex on them. Dropping constraints with `nOn ≤ 1` on a bounded region frees them.
- **Why it failed:** a constraint active at exactly one vertex of a convex region can still be
  ESSENTIAL. Removing it enlarges the piece, and an enlarged piece of `f*` has a SMALLER conjugate
  domain — so `f**` loses coverage somewhere else. Measured: with the drop, `f**` is exact at
  `(0.25,0.25)` and `(0.1,0.1)` and `+inf` at `(0.9,0.6)` and `(0.6,0.6)`; without it, exactly the
  other way round. Both are 5 of 7 and only one is sound.
- **A second, independent unsoundness found on the way, and worth its own note:** the boundedness
  test written for that drop read the vertex list AFTER `removeInfV`, which deletes the `±intmax`
  box-clip vertices that are the only mark of an unbounded region — so every region looked
  bounded. Fixed on its own merits (read it before `removeInfV`, which is the codebase's own
  convention), and it did NOT rescue the drop: the harmed piece is genuinely bounded.
- **A FOURTH attempt, also reverted: PARKING the slot claim instead of dropping the constraint.**
  Since only `endNv` constraints are ever read as edges — and on a lens that is ONE — the arc and
  its chord fight for one slot, and the scatter's last-write-wins hands it to the chord. Of the
  two, only the ARC has a two-dimensional dual cell (every point of an arc has its own normal, so
  it sweeps a cone; a straight edge has one normal and contributes a ray). So: give the arc the
  slot and PARK the chord's claim above the edge slots, dropping nothing — the region as a
  constraint set is unchanged, and only the label "this is edge number k" moves.
  The reasoning is sound and the result was worse: `conjugateOfPiecePoly` returned NO pieces at
  all, so it breaks something beyond the lens. Not diagnosed further.
- **Before retrying:** do not look for a better rule for dropping, and do not re-try parking
  without first finding out which OTHER piece the parking breaks. Give `conjugateOfPiecePoly` an
  explicit EDGE LIST instead of a count with two conventions (`endNv = nv` or `nv-1`; edge `j` at
  `ineqs(j)` or `ineqs(j+1)`). It cannot be done in that routine alone — `j` indexes
  `getNormalConeEdgeQ`/`Q3`'s output and `getSubdiffVertexT2`/`T2Q`'s `subdE` simultaneously, so
  all four move together. That is why the original "derive the edge list from the geometry" note
  underestimated the job.

## 2026-08-15 — Two of the five "remaining bugs" were described WRONG. Measure before fixing.

Both descriptions had been written from a symptom and carried forward as fact. Each cost an
attempt before measurement refuted it. Corrected shapes below.

> **Both are now resolved** — bug 5 fixed the same day (see its own entry), bug 1 taken from 0 to
> 5 of 7 probe points. This entry is kept for the corrected DIAGNOSES, which is what it is for.

### BUG 1 — "conjugateOfPiecePoly returns the conjugate of the chord"

- **What the record said:** the routine decides how many edges a piece has from
  `size(d.ineqs,2) == d.nv`, a COUNT standing in for "is this region unbounded", so the half-lens
  takes the unbounded convention and its arc is never read as an edge. Fix: derive the edge list
  from geometry.
- **What is actually true, measured** (instrumented dump of the pre-scatter constraint list):
  the half-lens arrives with `nv = 2` and FIVE constraints, and `edgeNo = [3 1 2 2 2]` —
  **three** constraints claim edge slot 2, all with both vertices on them:
      con 3: (s1+s2)^2 - 4*s1     con 4: (s1+s2)^2 - 4*s2     con 5: s2 - s1
  `getEdgeNosInf` numbers an edge by one of its ENDPOINT VERTICES, and a lens has two edges
  joining the SAME pair — so they are indistinguishable to it. The scatter
  `d.ineqs(edgeNo) = d.ineqs` is last-write-wins, so one edge is destroyed and another stored
  twice. Feed the arc first and both slots end up holding the CHORD; feed the chord first and both
  hold the arc. That is the whole of the "conjugate of the chord" symptom, and the count test is a
  consequence, not the cause.
- **Also measured, and this is why fixing the numbering is not enough:** hand-build the lens with
  ONLY its two genuine edges (`{(s1+s2)^2 <= 4*s2, s2 <= s1}` carrying `s1`, `nv = 2`,
  `nineq = 2`) and `conjugateOfPiecePoly` still returns **1 cell**, not the 4 the piece needs
  (2 vertex cones + 2 edge cells). So the downstream cell generation does not handle a 2-vertex
  CURVED region either.
- **Tried, twice, both REVERTED:** (i) spreading colliding genuine edges over distinct numbers;
  (ii) the same, plus dropping constraints that bound no edge — a constraint whose curve between
  the shared vertices leaves the region (which is what con 3 is, and which LP redundancy cannot
  see because it is not linear), and constraints active at a single vertex of a bounded region.
  Both make the lens's slots correct — measured, `nineq = 2` with the arc and the chord in
  distinct slots — and both take the second conjugation from WRONG VALUES to
  `QuaParCPLQ:conj:emptyResult`, no pieces at all.
- **Before retrying:** the edge numbering is necessary and not sufficient. Do the downstream half
  FIRST — make `conjugateOfPiecePoly` produce 4 cells for a bounded 2-vertex region with one
  curved edge, checked on the hand-built lens above, which needs no pipeline. Only then re-apply
  the numbering fix. Note also that `size(d.ineqs,2) == d.nv` is used a SECOND time, to choose
  between the polyhedral and the quadratic-aware normal-cone routines
  (`getNormalConeVertex`/`getNormalConeEdgeQ3` vs `getNormalConeVertexQ`/`getNormalConeEdgeQ`),
  so changing the count changes which of those runs.
- **Separately: the pinned test no longer fails the way it says it does.** Since `conj` began
  returning a MESHED QuaPar for a bounded multi-face domain (2026-08-13),
  `biconjugateOverATwoFaceSubdivisionIsTheEnvelope` ERRORS at `quaPolToPlq:curvedEdge` — the
  second conjugation is handed a curved QuaPar and `quaPolToPlq` requires a polyhedral domain. The
  symbolic route reaches the lens defect only when forced. Both need fixing; they are different.
- **Refuted on the way:** `convEnvCPLQ` is NOT a route to `f**` here. `f** = conv f` for a compact
  domain, but `convEnvCPLQ` is Step 1, a PER-TRIANGLE envelope with no cross-piece hull — measured
  on this input, it returns the per-piece envelopes (0.25 at (0.5,0.5) where the truth is 0).

### BUG 5 — "a piece that spans TWO sub-arcs of the SAME conic"

- **What the record said** (written 2026-08-14, from the piece's vertex list): piece 4 `src [1 6]`
  has g1's arc cut twice and spans both sub-arcs, so its one curve slot holds one and the other
  becomes a chord.
- **Measured, and it is false:** evaluating piece 4's own conic at its vertices gives
      V1 (-2.960609, 0.9606088) -> 0        V5 (-2.744821, 1.372827) -> 0
      V4 (-2, 2) -> +0.149        V2 -> +0.161        V3 -> +0.327
  with a tolerance of 8.8e-07. Only V1 and V5 are on it. The two curved edges are on **different**
  conics: piece 4's own slot holds the SPLITTING curve `{g1f1 = g2f6}`, and the flattened edge
  `(-2,2) -> (-2.744821,1.372827)` is a sub-arc of **g1's arc**, which its neighbour piece 5
  `src [2 6]` carries curved on `[-0.015625 -0.03125 -0.015625 -0.25 0.25 -1]`.
- **So it is the ordinary two-DIFFERENT-arcs case**, which `splitTwoArcPiece` exists for — and the
  cell does carry the arc going in: `src [1 6]`, `nv = 4`, `arcPos0 = 2`, arc
  `(-2,2) -> (-3.125,1.125)`. But `splitTwoArcPiece` is called exactly ONCE on this input and not
  for this cell, so the loss happens before it: either the crossed-arc restoration finds no
  position for the sub-arc (`findArcPosition` returning 0) or the piece is emitted before reaching
  it.
- **AND THE CAUSE IS NOW LOCATED, by that instrumentation.** For `src [1 6]` on that shift:
      HITS: 2 hits on edges [2 3] at (-2.744821,1.372827) and (-2.960609,0.9606088), arcPos0 = 2
      RESTORE half: nv=5 curveAfter=5 **p=4** straightCut=0
  So the crossed-arc restoration works perfectly: it FINDS the inherited sub-arc at edge 4, the
  splitting curve is genuinely curved (straightCut=0), and the half is correctly handed to
  `splitTwoArcPiece(half, 4, arcEc0)`. **`splitTwoArcPiece` then returns it UNSPLIT**, which its
  own header says it may do — "if neither works the piece is returned unsplit ... the assembly's
  own arrangement check is what would catch a dropped arc". It did, three stages later, as the
  orphan.
- **Why it finds no cut: THE TWO ARCS ARE ADJACENT.** They share vertex V5 (arc at edge 4, split
  curve at edge 5, `nv = 5`). Its two candidate chords are `arcPos+1 -> c+1` = V5->V1 and
  `arcPos -> c` = V4->V5 — but those ARE edges 5 and 4 themselves, so one chain comes out with 2
  vertices and the `numel(chain) < 3` guard skips both. The `nv == 3` fallback that handles the
  shared-vertex case (cut from the shared vertex to the midpoint of the opposite straight edge)
  does not apply at `nv = 5`.
- **FIXED 2026-08-15, exactly as prescribed:** the shared-vertex case is generalised to `nv >= 4`
  with the ordinary diagonal from S to a NON-ADJACENT vertex. On this piece S = V5 and `V5 -> V2`
  gives chains {2,3,4,5} (carrying edge 4, the inherited arc) and {5,1,2} (carrying edge 5, the
  split curve): one arc each. Guarded by `insideStraightHull` like the existing candidates, and
  each half passed through `splitAtReflexVertex`. **Seeded sweep 17 exact / 0 wrong / 1 errored →
  18 / 0 / 0**; `maxQuaParTest` 29 / 0; fast bucket 203 / 0. Pinned by
  `arcVsArcSplitsTwoADJACENTArcsOnAPieceWithADiagonal`.
- **Worth keeping:** the symptom pointed at the wrong stage twice. A straight edge facing an
  identical CURVED one at 8e-16 reads as a clip dropping a conic; the actual cause was a
  subdivision that declined to cut and announced it only by leaving an arc flattened. When a piece
  reaches assembly with an edge its neighbour calls curved, ask which producer returned it UNSPLIT
  before looking at the clip.
- Do NOT write a same-conic sub-arc splitter; that shape has never been observed, and the code for
  it was written and removed on 2026-08-15.

## 2026-08-14 — RESOLVED: the two-arc ray split, and why two attempts at it failed

- **Outcome:** `arcVsArcRefusesAnUnboundedTwoArcSplit` is GREEN and `maxQuaParTest` is 28 / 0 (from
  25 / 1). The entries below about the ray split stand as history; the construction was right in
  outline all along, and **neither reverted attempt failed for the reason it was thought to**.
- **What actually blocked it, both times.** Two defects, neither in the split:
    1. `pieceRecessionRays` took the parabola's axis from `arcNullDirs`, which solves `d·Q·d' = 0`
       *exactly* and returns **nothing** when `b²−4ac` comes out negative — which is what a
       floating-point parabola's `Q` does about half the time, being only semidefinite up to
       rounding (measured `−2.78e-17` on the pinned fixture's first arc). The derived chord was
       then silently never emitted, so the half's constraint region stayed a slab open at BOTH
       ends and `reccConeViolation` refused it. This is why **check (5) of the six did not
       separate the cases**: it was answering a question about a constraint set that was not the
       piece's. `parabolaArcFrame` has always taken the smallest-magnitude eigenvector; do that.
    2. `curveAfter ≠ 0` was read as "this edge is curved" in five places, when `boundedPiece` also
       sets it for a STRAIGHT splitting curve (`curveEc` all zeros). Two of those places call
       `parabolaArcFrame` on the zero conic and raise `degenerateAxis` — **so `MAXQP_ASSERT`
       crashed on three of the four arc-vs-arc fixtures, and the invariants that would have named
       all of this never ran on the inputs that needed them.** An invariant that errors is off.
- **And the seeded shift `[1.4979 3.6486]`, blamed on the split twice:** the split's halves pass
  all three exact invariants. The wrong 0.3531 came from a *different* piece, `src [2 4]`, whose
  only finite edge `polyConstraints` had dropped (defect 2) — the refusal had merely been masking
  it. It is now refused explicitly; see the next entry.
- **Lesson worth keeping:** when a gate refuses a construction you have independent reason to
  believe is right, suspect the gate. Both reverts were correct decisions on the evidence
  available, and the evidence was wrong because the tooling was broken in a way that was silent.

## 2026-08-14 — "`{f1=f2}` must be a PAIR of parallel lines here" (diagnosis, refuted by measuring)

- **Tried:** explaining the last wrong answer on seeded shift `[1.4979 3.6486]` — piece `src [2 4]`
  carrying g1 face 2's ZERO quadratic while g2 face 4 beat it by `+Inf` along its own ray — as a
  subdivision gap: `{f1=f2}` is a degenerate conic, so it can be a pair of parallel lines, and a
  half strictly on one side of the line `splitCell` cut along could be crossed by the other, which
  would leave through the recession cone and contribute no finite crossing to find. That story fits
  every symptom, and a guard was written for it (`assertWinnerHoldsAtInfinity`).
- **Why it failed:** the conic says otherwise. For that cell
  `diffRow = [0 0 0 0 0 0 0 −1.4979 −3.6486 5.4652]` — **its entire quadratic part is zero**, so
  `{f1=f2}` is a SINGLE straight line and nothing straddles it. `delta = 0`, `det3 = 0`,
  `eig(Q) = [0 0]`.
- **What it actually was:** the piece's only two vertices ARE the two crossing points, so both lie
  ON that line; `assignSide` had nothing to read the winner from and fell back to a centroid which
  is on the line too. The winner of a whole unbounded piece came out of floating-point noise.
  `assignSide` now reads such a piece in its RECESSION CONE, sharing the probe
  `assignSideFromCone` has used since it was written for the same problem in
  `splitUnboundedAtOneCrossing`. That shift now assembles CORRECTLY: the seeded sweep is 17 exact /
  0 wrong / 1 errored of 18, from 16 / 0 / 2.
- **Before retrying:** the guard is kept as a cheap exact backstop and currently fires on nothing.
  Do not read its existence as evidence that the pair-of-lines case occurs — no input has ever been
  observed to produce it. **Classify the conic before theorising about its shape**; one line of
  `eig(Q)` would have saved the detour.

## 2026-08-14 — Why check (5) may have failed to separate the two-arc ray split's cases

- **Tried (2026-08-13):** "each half's recession cone must equal the cone its own rays span" was
  one of the six checks applied to the unbounded two-arc ray split, and it did not separate the
  pinned fixture from the seeded shift that assembles to a wrong value. It was therefore recorded
  below as exhausted, and the search was pointed at WHERE the cut starts instead.
- **Why that conclusion may be premature:** the check is implemented by `pieceRecessionRays`, and
  until 2026-08-14 that routine derived the arc's chord — the constraint that makes a concave-side
  arc-piece compact at all — by asking which side the piece's other VERTICES fall on. That is the
  same rule `DECISIONS.md` already records as unsound one level up, in `QuaPar.chordCuts`, where it
  killed two green tests: a lens has both vertices ON the chord, so they say nothing. It also had
  no gate on *when* a chord may be emitted, so a piece that genuinely straddles the chord line
  could be handed one and be reported compact when it is not. Either way the check was answering a
  question about a constraint set that was not the piece's.
- **Corrected:** both questions are now settled by the conic itself. Along the chord `X0 + t·ch`
  the conic restricts to `q(t) = A·t·(t−1)` with `A = ch·Q·ch'`, because both endpoints are on it.
  `A ≤ 0` means the chord's interior is inside the kept side, the piece straddles the line and no
  chord may be emitted (the lens, and every convex-side face). `A > 0` means the piece touches the
  line only at the arc's endpoints and lies on the side the arc's own interior points are on,
  reached along the parabola's axis from the chord midpoint. The vertex test survives as a veto.
- **Status: CONFIRMED, and it was the whole of it.** Correcting the axis derivation is exactly what
  turned both halves of the pinned fixture's cut from "cannot be proved" into "proved". See the
  RESOLVED entry at the top.
- **Evidence:** `arcVsArcRefusesAnUnboundedTwoArcSplit` green; `maxQuaParTest` 28 / 0; fast bucket
  202 / 0.

## 2026-08-14 — A newly minted OUTGOING ray was given sign +1 (a live bug, not a dead end)

- **Recorded here because the reasoning is easy to re-derive wrongly.** `polyConstraints` reads a
  ray's outward normal as `sign · rot90ccw(direction)`, and both `dirIn` and `dirOut` store the
  direction pointing from the apex OUT to infinity. A piece is walked CCW with its interior on the
  LEFT, so the incoming ray is traversed along `−dirIn` and takes `+1`, and the outgoing ray is
  traversed along `+dirOut` and takes `−1`.
- **The bug:** both branches of `splitCell` that mint an escaping ray wrote `+1` for the OUTGOING
  one. That flips the kept half-plane to the far side of the ray's line, so the piece's constraint
  region is the reflection of its true region across that line — over-extended on one side, short
  on the other, which is the shape every far-field wrong answer here has had.
- **What NOT to do:** do not "fix" the INHERITED signs the same way. A sign is a property of the
  `P{k}` entry the ray came from, not of its role, and a face whose whole boundary is two rays
  sharing one apex can legitimately carry the same sign on both — that generalisation is itself a
  recorded fix, with its own test.
- **Evidence:** derived from `polyConstraints`' own HISTORY note; spotted independently while the
  ray split was being reverted and recorded in that commit. `newRaySign` now states it once.
  **Verified 2026-08-14**: no regression on any bucket, and an independent review confirmed the
  derivation (the previous `+1` gave the two halves of a split the *same* outward normal across
  their shared ray — overlap on one side, hole on the other).

## 2026-08-13 — Making the arc's chord a REAL EDGE by splitting the neighbour (option B)

- **Tried:** Nothing was built. It was proposed — and recommended — as the way to close a face
  whose arc is concave towards it, resolving the open choice left by the entry below.
- **Why it failed:** Geometry, checked before implementing. The chord runs through the **interior
  of the neighbour** on the other side of the arc, so making it a real edge splits that neighbour
  and leaves the offending face's own edge list unchanged. It cannot supply the missing constraint
  to the face that needs it. No arrangement can: the chord is not on that face's boundary.
- **Before retrying:** Do not. The chord must be **derived per face**, which is what
  `QuaPar.chordCuts` now does — and it resolved the whole far-field defect. Two soundness rules
  came out of building it, both learned by breaking tests: the side must be read from a point just
  INSIDE the face (the arc's midpoint stepped along the inward normal), never from the face's other
  vertices — a lens has both vertices ON the chord, so they say nothing and the wrong side is
  chosen; and a chord is emitted only when every other vertex and ray direction agrees, i.e. when
  it is redundant for the face, so it can never shrink a face below its true region.
- **Evidence:** `QuaPar.chordCuts` + `SUPPORT_MATRIX.md` §4.1. The vertex-based side rule broke
  `maxQuaParCurvedMatchesGroundTruthOnRandomQuadrilaterals` and
  `maxQuaParSplitsACellThatAlreadyCarriesAnArc`; the interior-point rule turns both green.

## 2026-08-13 — Splitting an UNBOUNDED two-arc piece with a ray (implemented, reverted)

- **Tried:** An unbounded piece carrying two arcs cannot be separated by a chord (no chord closes
  an unbounded piece), so: cut with a RAY from the vertex between the two arcs, along a direction
  the piece recedes in. Each half then keeps one arc, one original ray and the new one. It works —
  `arcVsArcRefusesAnUnboundedTwoArcSplit` assembles with it.
- **Why it failed:** It makes seeded shift `[1.4979 3.6486]` assemble to a value wrong by 0.3531.
  A silent wrong answer is worse than the refusal it replaces, which is what that test's own
  comment says, so it was reverted.
- **Before retrying, fix:** Six checks were tried and NONE separates the good case from the bad:
  (1) the ray recedes every straight constraint; (2) it stays inside BOTH arcs for its whole
  infinite length — necessary, a ray did leave through an arc; (3) each half admits a point just
  inside its own arc — necessary, it caught a mis-oriented pair; (4) the new ray's SIGN (an
  OUTGOING ray needs `-1` in `polyConstraints`' `sign*rot90ccw(dir)`, not the `+1` the
  escape-to-infinity branch uses — that branch looks wrong the same way and deserves its own
  check); (5) each half's recession cone equals the cone its rays span; (6) scoping to the
  half-strip shape. Start past all six: the open question is WHERE the cut starts, not which
  direction it takes — test whether the vertex between the arcs is the right starting point at
  all, or whether the cut must start on one of the arcs.
- **Evidence:** `maxQuaParTest/arcVsArcMatchesGroundTruthOverRandomShifts` (seed 20260803, N=18)
  versus `arcVsArcRefusesAnUnboundedTwoArcSplit`; commit "Revert the unbounded two-arc ray split".

## 2026-08-13 — "Equal areas and a shared polyline prove two halves tile" (refuted reasoning)

- **Tried:** Concluding, from a measurement, that a two-arc split was NOT responsible for a
  coverage hole: the two halves shared their cut polyline bit-exactly (17 digits), were both CCW,
  and their polygon areas summed to the parent's exactly. That sent the next session's search into
  assembly instead.
- **Why it failed:** Both facts were true and the conclusion was wrong. Area says nothing about
  which side of a **bent** boundary a point falls on. The halves did not tile: the cut polyline
  left one of them REFLEX at the bend, and every point-location test here reads a face as an
  intersection of half-planes, which is exact only for a CONVEX face.
- **Before retrying:** To decide whether pieces tile, use the per-edge cross product at the
  disputed point, not area. `splitAtReflexVertex` now splits such a half along a diagonal.
- **Evidence:** `arcVsArcDoesNotCrashOnSeededQuadSplits`, fixture 1, the point
  `(0.998629534754574, -0.0523359562429444)` — outside both halves, inside the parent.

## 2026-08-13 — Subdividing a bounded arc-piece whose constraint region is non-compact

- **Tried:** Twice, on the two quadrilateral-fixture pieces `src [1 2]` and `[1 6]`.
  Cutting the piece with a line parallel to the arc's chord, just past the arc's
  deepest point, to close the non-compact side.
- **Why it failed:** Both halves came back still non-compact, and the straight half
  wrongly kept the arc tag from `clipPolyHalfPlane`. The second attempt settled the
  reason: for a piece whose arc is **concave towards it**, the constraint set is a
  wedge intersected with the **outside** of a parabola. A cut parallel to the chord
  leaves the arc-side sliver still receding along the two side edges' own direction,
  which neither the cut nor the parabola blocks.
- **Before retrying, fix:** Do not retry any subdivision. This is a **representation**
  decision, not a subdivision bug. The constraint that would close the region is the
  arc's own **chord** — redundant for the region (which lies entirely on one side of
  it) but not one of its edges. So either a face may list a redundant bounding conic,
  or the chord becomes a real edge by splitting the neighbour across it. Pick one of
  those two before writing any more clipping code.
- **Evidence:** `maxQuaParTest`, quadrilateral fixture, pieces `src [1 2]` and `[1 6]`.
  (A separate `fixArcTag` fix has since stopped the losing half carrying the arc's
  constraint, but that does not make the cut work.)

## 2026-08-04 — Guarding against non-compact arc-pieces (four guards, all reverted)

- **Tried:** Four different guards, after establishing that arc-vs-arc results are
  correct near the arcs and wrong in the far field:
    1. Piece-level compactness guard on `triHalf` output.
    2. Piece-level compactness guard on every bounded curved piece.
    3. Post-assembly check: "a bounded face admits a far point", using QuaPar's own
       EC + P signs.
    4. Post-assembly check: "two faces disagree at a far point".
- **Why it failed:** (1) and (2) used the wrong sign convention and rejected valid
  faces — (2) errored on all three fixtures. (3) and (4) **detect correctly**, but
  they also fire on `[-1 0.75]` and `[2 -0.5]`, because those results genuinely *are*
  non-compact far out. A backstop that refuses non-compact results is therefore
  correct and would turn every arc-vs-arc result red.
- **Before retrying, fix:** Do not add a guard. The fix must be **upstream**: give an
  unbounded quadratic face's arc sub-pieces their **ray** boundaries, so they stay
  unbounded and compact-as-faces, instead of closing them with straight edges. That
  is a rework of the arc-vs-arc clip/split, not a guard.
- **Evidence:** Ring sweeps of `g.eval` against the pointwise max — `[-1 0.75]` 0/60
  wrong at radius 8 but 2/60 at radius 30 (worst 6.1); `[2 -0.5]` 2/60 at radius 8,
  11/60 at radius 30 (worst 33.8). `[0.5 0.5]`: assembled FACE 15 carries g1-face-1's
  quadratic over a tiny triangle near (-2,2) yet admits (-3.98,0.61) two units away,
  overlapping the correct zero face (FACE 9).

## 2026-08-03 — `decideWinner`: sampling `diff` at 1-D stationary points

- **Tried:** Making `decideWinner` sound on unbounded cells by additionally sampling
  `diff = f1 - f2` at the 1-D stationary point of each edge and each ray.
- **Why it failed:** Sound and strict, but it does not fix the defect and was
  confirmed neutral on all three fixtures. `{diff = 0}` is a parabola (rank-1 `Q`,
  which `splitCell` already asserts) whose branches run to infinity along `Q`'s null
  direction, strictly *between* the cell's two bounding rays. The sign change is off
  the 1-skeleton entirely, so no amount of sampling along edges and rays sees it.
- **Before retrying, fix:** Do not attempt a probe-based patch. A real fix needs
  **both** a sound sign-over-the-cone test in `decideWinner` **and** a `splitCell`
  that can cut a cell along a parabola entering and leaving at infinity — the current
  one asserts exactly two finite crossings, and here there are zero.
- **Evidence:** Cell `src[4 3]` stores `winner=g2` but `g1=0 > g2=-8.648` at the
  interior point `s=(-5.4843,1.5866)`; symmetrically `src[6 5]` stores `winner=g1`
  while g2 wins at `s=(-1.1298,4.1007)`. 4/68 samples, worst 8.648. Kept out of the
  commit.

## 2026-08-03 — Returning both components from a disconnecting curved cut

- **Tried:** Replacing `clipPolyByConic`'s blanket refusal of a disconnecting curved
  cut with "return both components".
- **Why it failed:** Two components are **unrepresentable**. A separated survivor
  would need the cutting conic running to infinity, i.e. an unbounded curved edge,
  which `QuaPar` cannot hold. Separately, the premise was wrong: the cut measured
  **connected** on all three fixtures, so the single-cell construction was right.
- **Before retrying, fix:** `Do not retry` unless `QuaPar` gains unbounded curved
  edges. The blanket refusal was instead replaced with a real connectivity test,
  which now refuses only the genuinely disconnecting case — and that case does not
  arise on these fixtures.
- **Evidence:** `partitionReport` OK on all three fixtures after the connectivity
  test landed.

## 2026-08-13 — Refuted diagnoses (the analysis was wrong, not just the fix)

Each of these was written up in detail and then disproved. Re-deriving any of them
from the same symptoms is the failure mode this entry exists to prevent.

- **"Piece 5 (`src[2 2]`) is over-extended and its ray should terminate."** Refuted
  by the sign data: `evalConic(Ecut)@V = [0, -0.046, -0.015]`, so the vertex `(-2,2)`
  sits on the g2-face-1 (discard) side and the far ray is legitimately g2-face-2. The
  clip is correct. The geometric measurements in that write-up are still good; the
  conclusion is not.
- **"The missing g1f6 ∩ g2f2 mirror is the defect."** The sharper diagnosis that
  replaced the above, itself now superseded — the mirror (`src[6 2]`) is built
  correctly.
- **"`decideWinner` wrongly declares an unbounded cell decided" (for `[2 -0.5]`).**
  The real cause was `assignSide` reading the winner at `piece.V(2,:)`, which for
  `splitCell`'s unbounded "rest" piece (`curveAfter=1`) is a **crossing point** where
  `diff ≈ 0` and its sign is noise, so both halves could get the same winner. Fixed
  in `53fc9fd` by reading the vertex farthest from `{diff=0}`. The parabola-to-infinity
  write-up was the wrong lead for this fixture.
- **"Piece 16's ray is piece 5's partner."** It is not — piece 16's ray lies on
  `x+y=0.5`, a parallel but different line from piece 5's `x+y=0`. They were never
  meant to match.

## 2026-08-03 — Retired hypotheses about the orphan ray — do not re-try

- **Tried / Why each failed:**
    - *The orphan ray is one physical ray covered to different extents* — nothing
      lies on it.
    - *The cut must be restricted to the arc's own span* — it must not; the straight
      half-planes are applied first.
    - *`clipPolyByConic` emits clockwise cells* — a guard on every bounded output
      never fires.
    - *`polyL` is non-convex* — all four curved faces measured convex.
    - *Pair `(6,1)` is spurious* — it is a genuine thin cell. The "zero intersection"
      reading was a resolution artifact: a 0.0625-spaced grid cannot see a cell of
      area 0.008.
- **Before retrying, fix:** `Do not retry.` Each was measured, not argued.

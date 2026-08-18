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
* **Cost:** one run, about 10 minutes of compute.

### A2. Fix the DIAGONAL terms over a box -- the one correctness gap on SCIP's path.

`x^2 - y^2` and `(x^2+y^2)/2` on the unit box raise `MATLAB:badsubscript` in the SECOND
conjugation (`functionNDomain.conjugateOfPiecePoly`). QPLIB objectives are sums of `x_i*x_j`
**including i = j**, so this is not optional -- and per section 0 it is one of only two cases
where CCA2 has anything to offer over McCormick.

* **Do:** reproduce at unit level (build the piece directly, as `functionNDomainTest`'s fixtures
  do -- not through `biconj`, which buries the evidence), then fix.
* **Closes when:** both shapes return an envelope matching the lower convex hull of the sampled
  graph, pinned by a test in `biconjugateTest`.
* **Risk:** unknown until reproduced. This is the one genuinely open piece of mathematics here.

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

### A4. State what is deliberately NOT fixed, and confirm it is off SCIP's path.

`SUPPORT_MATRIX.md` section 8 lists open items -- `partialConj` unimplemented, the `pqp` and
`graph` engines missing, `RatPol.conj`/`biconj`/`add` missing. Section 0.0 already records that
SCIP calls none of them.

* **Do:** confirm against the current `cca2ConvexEnvelope.m` and write the exclusion into 0.0.
* **Closes when:** the SCIP-relevant surface is a short, explicit list of entry points.

---

## Phase B -- the direct-formula / symbolic map (the "clear picture" gate)

**This does not exist today in one place, and it is the deliverable most likely to change what we
build.** The question is not academic: anything with a direct formula should not be paying for the
symbolic engine, and anything that needs the engine has a cost floor no tuning removes.

### B1. Build the map.

One table, per input shape, with four columns: **route taken**, **closed form or symbolic**,
**where the formula is** (source, and the [COAP] appendix where applicable), **measured cost**.
The shapes that matter:

| shape | expected route |
|---|---|
| bilinear over a BOX | McCormick / Al-Khayyal-Falk -- closed form |
| bilinear over a TRIANGLE, 0/1/2 convex edges | cPLQ Step 1 closed forms (A.2 / A.3 / A.4) |
| bilinear over a TRIANGLE, 3 convex edges | A.5 split into A.4 -- closed form, SURD coordinates |
| bilinear over a general POLYGON | triangulate, per-piece closed form, **Step 3 symbolic max** |
| convex quadratic over a piece | `conjConvexOverPiece`, KKT active set -- closed form per cell |
| diagonal / indefinite non-bilinear quadratic | `xyFrame` change of variables, then the above |
| anything UNBOUNDED | `fanUnboundedFace` + `convEnvUnbounded` -- closed form per shape |

* **Closes when:** every row cites its routine, says closed-form or symbolic, and carries a
  measured time from the current tree.
* **The expected conclusion, to be confirmed or refuted:** every PIECE has a closed form, and
  **all** the symbolic cost is Step 3's cross-piece maximum (`functionNDomain.maximumP`) -- which
  is why the quadrilateral spends about 35 s in Steps 1+2 and about 43 minutes in Step 3.

### B2. Draw the consequence for SCIP.

If B1 confirms that a single piece is always closed-form, then a separator needing one term's
envelope over a box **never needs Step 3 at all** -- and the 40-60 s per term is being paid for
machinery the caller does not use. That turns the integration from "call `biconj`" into "call the
per-piece closed form", and it is the largest performance lever available.

* **Closes when:** the answer is written down either way, with the entry point named.

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

All of the following, each with a measurement on the current tree:

1. `checkBoxEnvelopeForSCIP` shows **no ERROR rows** (A1, A2).
2. The direct-formula / symbolic map exists and names the entry point a separator should call
   (B1, B2).
3. Per-term cost is measured, attributed, and inside a **stated** target (C1-C3).
4. Full suite green -- currently **332 / 0** -- with A2's new tests in it.

Only then: wire the bridge, expose value and subgradient off whatever B2 names, and run QPLIB.

---

## Sequencing

A1 is cheap and goes first -- it may change A2's shape. B1 can run in parallel with A2, being
mostly reading and measuring rather than fixing. C is last, because B2 may remove most of the
problem C is trying to solve.


---

## Phase B, measured: PROFILE OF ONE FOLD (2026-08-18)

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

# Session Handoff

_Last updated: 2026-08-02 (overnight session)_

## What happened this session

**A silent wrong answer in the cross-piece maximum, found by disbelieving the previous
diagnosis.** The handoff said the two-face defect lived in the second pass and that "both first
conjugates are verified correct against brute force at 7 dual points". A grid audit says the FIRST
conjugate was wrong at **800 of 40000 points, worst error 0.48**. The 7 points missed it. Root
cause: `region.maxArray` decided which of two candidate functions dominates a cell by PROBING one
interior point, and a probe is evidence, not proof. That is now replaced by a sound global sign
test, which also fixed next-step 4's 4-cone case (it no longer returns a wrong number).

**Performance: the symbolic-engine round trips roughly halved** -- 10809 -> 5591 `subs` calls on
the two-face `f*` -- by batching substitutions that were being done one constraint at a time. No
arithmetic changed, and the evidence for that is stronger than a test count: the batching-only
tree's suite results are BYTE-IDENTICAL to the baseline's, suite by suite. No wall-clock figure is
quoted, deliberately; see next-step 3.

**The two unreproducible measurement claims are now reproducible**, via two committed, seeded
sweeps -- and one of them retired a defect: `QuaPar.eval`'s at-a-vertex fix is verified.

## Where things stand

- Branch: `main`. Commits this session, oldest first:
  - `46fac7c` maxArray's single-probe dominance + `getNormalConeVertexQ`'s out-of-range index
  - `375f59d` `sweepQuaParEvalAtVertices` (seeded) -- verifies the `QuaPar.eval` fix
  - `3badfb7` `sweepMaxQuaParCurvedSplit` (seeded) -- retires maxQuaPar's curved-split numbers
  - `f861db3` batched substitutions + the sound max + the Step 3 emptiness veto (the one that
    repairs `46fac7c`, see below)
- **Suite: COMPLETE and green.** `300 passed / 1 failed / 0 incomplete over 26 suites`, the one
  failure being `biconjugateTest/biconjugateOverATwoFaceSubdivisionIsTheEnvelope`, which fails by
  design (next-step 1). Compare `46fac7c` mid-session at 298/2/1. Run against a SNAPSHOT, not the
  working tree -- see the harness note.
- Pushed: **pending** -- push before doing anything else.
- **`46fac7c` on its own is BROKEN** and the commit after it repairs the damage: its
  "all probes must agree" rule made `maxArray` answer "undecided" often enough to send
  `maxEqDom` into `splitmax3` with a RATIONAL pair, which `region.normalize1` cannot represent.
  Measured on `testMaxMultiRegion/testMax`: pre-session 1/0, `46fac7c` 0/1, final tree 1/0. If
  you ever bisect through this session, expect that one commit to fail.
- **Tag `v0.1` is PUSHED and public.** Moving it is not cheap; anything further should be `v0.1.1`.

### THE ALGORITHM (unchanged; get this right)

```
f  -> triangulate -> convex envelope of f+I_T per triangle T (splitting as needed)
                  -> conjugate of g+I_T per piece
                  -> MAX of all those conjugates                       = f*

h = f* on pieces P (P polyhedral, or with a single parabolic edge)
   -> conjugate of h+I_P per piece
   -> MAX of all those conjugates                                      = f**
```

- **Step 1's output is NOT the convex envelope of f over its domain.** It is a per-triangle
  intermediate. The envelope over a box/polygon is `f**`, i.e. two passes.
- **The max is the same code in both passes** and must not be duplicated: `functionNDomain.maxOfList`,
  called by `plq.maximumConjugate`, `plq.biconjugateF` and `QuaParCPLQ.conj`.
- `ia` from `conjugateOfPiecePoly` is only the block delimiter.
- `region.plus` is INTERSECTION. `functionNDomain.mtimes` (`*`) does NOT multiply: it intersects
  each pair of regions and stores BOTH functions there so `maximumP` has two values to compare.
- **`maxOfList` intersecting domains is CORRECT, not a bug.** `f**` is a MAX of extended-real
  functions, so its domain is the INTERSECTION of the per-piece conjugates' domains. One piece
  whose conjugate comes out on too small a set collapses the whole result -- which is exactly what
  the open defect below does.

## Next steps

0. **`biconj` does NOT work for a general polygon, and the reason is ORIENTATION.** Measured with
   `checkBiconjDomainCoverage` (committed; re-run it, it prints the table). Triangles work for
   convex/indefinite/concave/affine `f`; the unit box and an axis-aligned `[0,2]x[0,3]` box work;
   all three unbounded families tried work. But:
   - a **general convex quadrilateral** errors with `MATLAB:badsubscript`, because `triangulate`
     gives it a piece with **`nCE = 3`** and cPLQ's Step 1 has NO `nCE == 3` branch --
     `convexEnvelope` returns ZERO envelope pieces and `plq_1p.conjugateFunction`'s
     `for i = 1:max(1, size(obj.envelope,2))` (a guard written for "triangles where the envelope
     is not computed") then indexes `obj.envelope(1)`. `conjCPLQ.m`'s header already flagged the
     missing branch; what is new is that an ordinary quadrilateral reaches it.
   - a **parallelogram** gets `nCE = 1` pieces and rational envelopes, survives Step 1, and then
     dies in the SECOND conjugation with `QuaParCPLQ:conj:emptyResult`.

   `nCE` counts edges of positive finite SLOPE, so an axis-aligned box gives `nCE = 0` everywhere
   and sails through while a sheared one does not. **Rotating a working domain breaks it.** This
   is the biggest gap between what the toolbox does and what "biconjugate of a PLQ function"
   promises, and it is upstream of everything in next-step 1.

1. **The remaining half of the two-face defect: `conjugateOfPiecePoly` counts edges instead of
   looking at them.** This is the sharpest lead in the repository and the mechanism is TRACED, not
   guessed (instrumented run, notes in `biconjugateTest.m` at the failure site and
   `SUPPORT_MATRIX.md` section 7).

   `f*` for the two-face unit square is now exact and arrives as ten pieces, two of them
   half-lenses such as `{(s1+s2)^2<=4s1, (s1+s2)^2<=4s2, s2<=s1}` carrying `s1`. That piece is
   BOUNDED with 2 vertices and 2 edges (one parabolic arc, one segment), so its conjugate is finite
   on all of R^2 and needs 4 cells. It gets 2, and they are the conjugate of the straight CHORD
   between its vertices.

   Why: the routine decides how many edges a piece has from the COUNT `size(d.ineqs,2) == d.nv`,
   which is standing in for "is this region unbounded". Here `nv = 2` with FIVE stored constraints
   (the scatter duplicates `s2-s1` into slots 2 and 5 and parks the arc's conic in slot 4), so the
   count says "unbounded": `endNv = nv-1 = 1`, edge `j` is read from `ineqs(j+1)`, one edge cell is
   built from the straight edge, and **slot 4 -- the arc -- is never read as an edge at all**.

   Fix direction: derive the edge list from the geometry (`region.vertexOfEdge` already answers
   "which constraints bound an edge") instead of from that count. Alternative worth weighing first:
   pieces 1-4 of `f*` all carry `s1` and their union is the POLYHEDRAL cone `{s1>=0, s1>=s2, s2<=1}`,
   so an exact `merge` across a conic boundary would remove every curved piece here and the second
   pass would have nothing curved left to do -- that is `region.merge`/`unionIsExact`, which today
   cannot certify a union across a conic edge.

2. **Step 3's cross-piece max: the DROP is fixed, an OVER-claim remains.** The 4-cone fan
   (`conjCPLQTest/step3DropsCellsOnSomeUnboundedAssemblies`) used to return `f*(-0.5,2) = 1.125`
   for a truth of 2, because the cell carrying `s2^2/2` on `{s1<=0, s2>=0}` was lost. Traced to
   `region.simplifyUnboundedRegion` declaring the split half
   `{s1<=0, s2>=0, s1^2/2 - s2^2/4 <= 0}` EMPTY -- a genuine 2-D cone containing `(-0.5,2)`,
   `(-0.1,3)`, `(-1,4)`. It decides emptiness from probe directions built out of constraint
   SLOPES at a vertex, and the split conic's gradient VANISHES at exactly that vertex, so those
   directions are meaningless. `region.witnessAwayFrom` now refutes an emptiness verdict by
   exhibiting a feasible point away from the vertices; that is sound in one direction only (an
   empty region has no such point, so no true verdict can be overturned) and it only ever vetoes.
   8 cells now assemble and all 8 probes including `f*(-0.5,2) = 2` are exact.

   **What remains is the opposite error at a different point.** `assertStep3MatchesPieces` now
   fails at `s = (-3,-2.4)`: the assembly gives `5.130`, the per-piece max gives `4.500`, and
   `4.500` is right (the four cone suprema there are 0, 4.5, 3.69, 2.88, by hand). `5.13` is
   `s1^2/4 + s2^2/2`, i.e. FACE 4's cell -- which belongs on `{s1>=0, s2<=0}`, so some region has
   grown across `s1 = 0`. Start there: which region, and at which fold. `conj` still refuses to
   return the result (`PLQ:conjCPLQ:cplqFailed`), so nothing wrong is handed out.

   Note the same anti-pattern behind both this and next-step 1's `maxArray` bug: **a geometric
   decision taken from a handful of probe points stepped off a vertex.** `simplifyUnboundedRegion`
   still does that everywhere except the two verdicts now vetoed. If a third wrong answer shows up
   in this area, look there first.

   **A structural limit found on the way, and worth fixing on its own:** `maxArray`'s refusal to
   guess is applied only to a POLYNOMIAL pair, because a region's constraints must be polynomial
   and `splitmax3` hands `f1 - f2` straight to `region()`. On a rational pair -- which is EVERY
   pair on the biconjugate path -- refusing produced `symbolic:coeffs:NotAPolynomial` out of
   `region.normalize1`, so there the old, unsound vertex verdict still stands. Closing that means
   teaching `splitmax3` to clear denominators, which is sound only where both denominators are
   provably nonzero on the cell. Until then the biconjugate's max is still guessing.

3. **Performance -- measured in OPERATION COUNTS, not seconds; more left to take.** On the
   two-face `f*` (`profile on`, then attribute by PARENT, not self time):

   | | `subs` calls | `solve` calls |
   |---|---|---|
   | before this session | 10809 | 438 |
   | after the batching alone | 5161 | 330 |
   | after batching + the correctness work | 5591 | 378 |

   So the batching removed **~half** the round trips to the symbolic engine, and the sound max
   then added ~8% back because it splits cells the old guess did not -- that is the price of the
   right answer, and it is small. **Do not quote a wall-clock number from these runs**: three
   suite sweeps were running concurrently and the timings moved by 40% between otherwise
   identical runs. A controlled measurement still needs to be taken (CLAUDE.md section 3).

   What is left, by total time in the last profile: `getVertices` 60 s (now the largest),
   `mergeL` 55 s, `simplifyUnboundedRegion` 47 s, `region` constructor 34 s, `removeTangent` 21 s,
   `assertStep3MatchesPieces` 19 s (a self-CHECK, 13% of the run). The batching pattern that won
   here -- substitute into the WHOLE constraint vector in one call, since `subs`/`isAlways` are
   elementwise -- has not been applied to any of those yet. `region.slopes2` is now uncalled and
   its cost is gone entirely.

4. **`QuaPar.eval` at a vertex -- DONE, verified.** `sweepQuaParEvalAtVertices(20260802, 200)`:
   225 of 1205 subdivision vertices are located by NO face under the pre-fix exact test and 0 under
   the current tolerant one, with all rings of radius 1e-8 correct. Pinned deterministically by
   `QuaParTest/evalLocatesAPointExactlyAtItsOwnVertex`.

5. **Reproducibility debt -- 2 of 3 discharged.** `sweepQuaParEvalAtVertices` and
   `sweepMaxQuaParCurvedSplit` are committed and seeded; `SUPPORT_MATRIX.md` section 0.1,
   `DESIGN.md` and `maxQuaPar.m`'s header now quote them instead of the retired figures. Still
   unreproducible, and flagged: the **58 -> 76 of 395** quadrilaterals, a before/after pair across
   a code change whose population no longer exists.

   Worth acting on from the new sweep: **80 of 131 random diagonal splits never reach `maxQuaPar`
   at all** -- Step 2 refuses the triangle with `conjPieceCPLQ:notImplemented`. The old "85 of 115
   assembled" concealed that the dominant obstacle is UPSTREAM of `maxQuaPar`.

6. **`mergeL`/`removeTangent` exact-tie-point -- NO LIVE REPRODUCER; do not spend time on it
   without one.** Its one recorded case (`s=(0.5,0.5)`, `x*y` on the unit square as two triangles)
   is exact on `46fac7c` and now: 0 wrong, 0 uncovered over a 25x25 grid on `[-3,3]^2` plus both
   symmetry axes. It was closed by the `maximumP` fall-through-to-`splitmax3` fix and the entry
   outlived it. A latent same-family defect WAS found and fixed: `removeTangent` took `s0(1)` as
   the second point on the parabola before checking the root set was non-empty and real.

7. Unimplemented: `partialConj` (all engines/types), `'pqp'`/`'graph'` engines,
   `RatPol.conj`/`biconj`/`add`.

8. **`QuaParCPLQ` data-shape wart.** Uniform in TYPE (all `conj` results are `RatPar`, with
   `kind()`), but it stores `fnd` and leaves `V/E/Ec/f/F` empty, and its `add` refuses a
   `QuaPol`/`QuaPar` operand, so cross-case `infConv` errors. `RETURN_TYPE.md` records this.

## SCIP (`AI/spike/SCIP`) -- read before renaming anything

SCIP bridges in through the MATLAB Engine and calls **exactly one entry point, `convEnvCPLQ`**,
via its own glue `SCIP/src/cca2ConvexEnvelope.m`, consuming the returned `RatPol`'s
`V/E/F/f(:,5:10)/den` as plain arrays. It never calls `conj`, `biconj`, `partialConj`, or any
`QuaPar` method -- so **none of the open items above is on its path**. `QuaPoly.m` is an alias its
glue depends on; do not remove it. See `SUPPORT_MATRIX.md` section 0.0.

## Do NOT redo these -- tried, with reasons

- **Believing the previous diagnosis of the two-face defect.** It said the first pass was verified
  correct at 7 dual points and blamed the second. The first pass was wrong at 800 of 40000 grid
  points. **Audit on a grid, not on hand-picked points.**
- **Deciding a max by probing.** `maxArray` used to step 0.1 off a vertex along constraint-slope
  bisectors and take the first feasible probe's verdict. It was wrong on two independent inputs.
  A probe is evidence, never proof. The whole apparatus (`slopes2`'s only caller, `probeAlong`,
  `probePerp`, `maxFromPt`) is deleted.
- **Applying `maxArray`'s "refuse to guess" to a RATIONAL pair.** Tried, and it breaks two
  things: `splitmax3` feeds `f1 - f2` to `region()`, whose `normalize1` raises
  `symbolic:coeffs:NotAPolynomial` on a rational constraint, and it also evaluates `f1 - f2` at a
  vertex that can be a pole (`symbolic:kernel:DivisionByZero`, now guarded). Measured cost of the
  unrestricted version: `conjCPLQTest/multiFaceUnboundedConvexFacesConjugateExactly` went from
  pass to ERROR on a full sweep. The refusal is deliberately confined to polynomial pairs.
- **Asking `isAlways` about a sign WITHOUT re-declaring the variables real.** A bare `sym('s_1')`
  is COMPLEX to the engine, so `-s_1^2/2 <= 0` is simply not true and no quadratic difference can
  ever be decided. Every comparison then fell through to an unnecessary split, the 4-cone fan
  produced the vacuous constraint `-s_1^2 <= 0` as a split boundary, and `removeTangent` crashed
  downstream. `maxArray` substitutes fresh real symbols before asking.
- **Reading `convEnvCPLQ`'s output as "the convex envelope"** and concluding a box is unsupported.
  It is Step 1's per-triangle intermediate; the box's envelope is `f**`.
- **Replacing `addEq`'s intersection with a union.** Adding equations IS intersection, and that is
  correct: `f**` is finite only where EVERY piece's conjugate is finite.
- **Adding a second, private max inside `biconjugateF`.** Use `functionNDomain.maxOfList`.
- **The `isQuad` chord endpoints from the conic's own vertices** (instead of `vx(1),vx(2)`).
  Changed nothing on the two-face case; reverted. It is not even reached on that input -- the
  rewrite is guarded by `f.isQuad` and those pieces' `f` is affine.
- **Returning recession directions as `(cos t, sin t)`.** A direction of a rational half-plane is
  rational; the round trip puts `6.123e-17` where a `0` belongs. Pinned by `regionTest`.
- **Probing every envelope for affineness in `plq_1p.conjugateFunction`.** The `nCE==1` envelope
  is RATIONAL and its denominator vanishes at `(0,0)`.
- **`merge` guarded by constraint-set EQUALITY** -- provably sound and makes things WORSE
  (36 -> 125 wrong of 289).
- **`slopes2` SKIPPING a curved constraint with no vertex** -- `maximumP 3` went 153 s -> >90 min.
  (`slopes2` is now uncalled, but the trap stands for any future caller; its header says so.)
- **Replacing `conjugateOfPiecePoly`'s scatter with a sort** -- the scatter legitimately GROWS the
  array and is last-write-wins; sorting reproduces neither.
- **PARKING vertexless constraints above the real edge slots** -- the `isQuad` branch then chords
  them from unrelated vertices.

## Environment / harness

- MATLAB R2024b is **network-licensed only**. Off the UBC VPN, `matlab -batch` dies with
  `Licensing error: -96,7`. **Connect the VPN before any session that runs anything.**
- **`.claude/suite.sh`** -- one MATLAB per suite under `timeout`; default 3600 s,
  `CCA2_TEST_TIMEOUT` to override, `0` disables, `TIMEOUT` reported separately from failure.
  Also takes suite names: `bash .claude/suite.sh biconjugateTest conjCPLQTest`.
- **Run long sweeps against a SNAPSHOT, not the working tree.** `git archive HEAD | tar -x -C <dir>`
  (or `cp *.m`), then `bash <dir>/.claude/suite.sh` -- `suite.sh` takes its root from its own
  location. A sweep against the live tree is invalidated the moment you edit a file, because each
  suite starts a fresh MATLAB that reloads whatever is on disk. That happened this session and cost
  a full sweep. Snapshots also let two or three sweeps run concurrently; the machine has 16 cores.
- Long waiters must watch a **per-run** log filename; a reused name makes the exit condition
  unsatisfiable.
- Profiling: `profile on -historysize 40000000`, then attribute by PARENT
  (`ft(j).Parents`) rather than by self time -- everything expensive is inside the Symbolic Math
  Toolbox, so self time tells you nothing about which of YOUR routines to fix.

## Relevant files

- `region.m` -- `maxArray` (the sound sign test, and the header explaining why probing is gone);
  `ptFeasible`, `linearForm`, `getVertices` (batched substitutions + closed-form linear/linear
  intersection); `probeVertexIndex` (the `getNormalConeVertexQ` index guard); `removeTangent`
  (the root-set guard); `slopes2` (now uncalled).
- `functionNDomain.m` -- `maxOfList`, `maximumP`, `conjugateOfPiecePoly` (the `endNv` edge-count
  defect of next-step 1 lives here), `getInterior` (indexes `c2(2)` under a guard testing only
  `size(c1,2)` -- next in the latent chain, NOT fixed).
- `evalFunctionNDomain.m`, `symbolicFunction.m` (`subsF`) -- batched substitutions.
- `biconjugateTest.m` -- 7 tests, one failing by design. **Read the failure-site comment before
  touching the parabolic branch**; it carries the traced mechanism.
- `sweepQuaParEvalAtVertices.m`, `sweepMaxQuaParCurvedSplit.m` -- the two seeded sweeps. Neither
  matches `*Test.m`/`test*.m`, so neither joins the suite; run them by hand.
- `SUPPORT_MATRIX.md` -- section 0.0 downstream consumers, 0.1 the reproducibility rule (now with
  two discharged rows), 7 known wrong-answer defects, 8 the blocker list.

# Session Handoff

_Last updated: 2026-08-14_

## What happened this session

**`maxQuaParTest` is GREEN — 28 pass / 0 fail, from a baseline of 25 / 1.** The last arc-vs-arc
red, `arcVsArcRefusesAnUnboundedTwoArcSplit`, is closed, and so is the covering half of
`FARFIELD_FIX_PLAN.md` Phase 4.

The headline is not the fix, it is why the fix took three attempts:

> **The tooling that judged the two previous attempts was itself broken, in two ways, and
> silently.** `MAXQP_ASSERT` was *crashing* (`degenerateAxis`) on three of the four arc-vs-arc
> fixtures, so the invariants that eventually named three defects had never run on the inputs that
> needed them. An invariant that errors is an invariant that is off, and nothing was noticing.

Both earlier reverts were correct decisions on the evidence available. The evidence was wrong.

## What was fixed

1. **`pieceRecessionRays` took the parabola's axis from an exact discriminant.** `arcNullDirs`
   solves `d·Q·d' = 0` exactly and returns **nothing** when `b²−4ac` comes out negative — which is
   what a floating-point parabola's `Q` does about half the time, being only semidefinite up to
   rounding (measured `−2.78e-17`). The derived chord was then never emitted, the piece's
   constraint region stayed a slab open at both ends, and `reccConeViolation` refused it. This is
   why check (5) of the six recorded heuristics never separated the good case from the bad.
2. **`curveAfter ≠ 0` does not mean "this edge is curved"** — `boundedPiece` also sets it for a
   STRAIGHT splitting curve, where `curveEc` is all zeros, and `pieceIsCurved` exists to say so.
   Five sites read the tag as "is curved" anyway: `polyConstraints` emitted **no half-plane at all**
   for an ordinary straight edge (a piece was admitted two units outside itself, answering `0`
   where the truth was `0.35310191`); `pieceStraightEdges` skipped it, blinding every boundary
   minimisation built on that list; and `containmentViolation`/`boundaryMinOf` called
   `parabolaArcFrame` on the zero conic — the `degenerateAxis` crash above.
3. **A whole unbounded piece could have its winner decided by floating-point noise.**
   `splitCell`'s unbounded "rest" piece can come out with exactly the two crossing points as its
   vertices — both, by construction, ON `{f1=f2}` — so `assignSide` had nothing to read from, and
   its centroid fallback is on that line too. Now read in the piece's RECESSION CONE. The seeded
   sweep went **16 exact / 0 wrong / 2 errored → 17 / 0 / 1** of 18.
   This one was nearly mis-diagnosed as a subdivision gap (`{f1=f2}` a *pair* of parallel lines);
   classifying the conic refuted it in one line — its whole quadratic part is zero. `DECISIONS.md`
   records that, because the symptom is a convincing fit for the wrong story.
4. **A newly minted OUTGOING ray was given sign `+1`** where `polyConstraints`' convention needs
   `−1`, giving the two halves of a split the same outward normal across their shared ray.
   `newRaySign` states the derivation once.
5. **The two-arc ray split is restored**, gated on all three exact invariants per half
   (containment, recession cone, winner domination) rather than on heuristics; when nothing can be
   proved it returns `[]` and the caller refuses exactly as before.
6. **`verifyFacesCoverThePlane`** — the covering proof. Four checks on the constraint data that
   together force the boundary of the union of the faces to be empty. An independent review found
   three routes by which it could have passed *vacuously*; all three are closed, and
   `coverProofRejectsBrokenArrangements` breaks a certified result three ways and requires a
   finding each time.

## Where things stand

- Branch: `overnight/2026-08-13` (not merged, not pushed).
- **`maxQuaParTest` 28 / 0. Fast bucket 203 / 0** (was 200 / 1). **Normal bucket 6 / 0.** The slow
  bucket had not finished when this was written — check it before merging.
- `MAXQP_ASSERT=1` and `=2` now run clean on all four arc-vs-arc fixtures.

## Next steps

1. **Finish the slow bucket and merge.** Nothing on this branch touches the symbolic path, but the
   slow bucket is the only thing not yet re-measured.
2. **The one open `maxQuaPar` case: an orphan boundary edge in `assemblePieces`**
   (`maxQuaPar:internal`). An ERROR, not a wrong answer, and it was masked until now by the two-arc
   refusal upstream. Reproducer: seeded shift `[-2.6434 -1.8066]` of the two-curved fixture —
   piece 4, `src [1 6]`, straight edge `(-2,2)→(-2.744821,1.372827)`, no matching neighbour.
   `checkOrphanHalfEdges` prints the closest candidates; `insertGlobalPassthrough` is what handles
   the T-junction form of this.
3. **Turn `MAXQP_ASSERT=1` on in the test suite.** It is cheap, and this session is the argument:
   it was off *and* crashing, and three defects lived behind that for weeks.
4. **Then SCIP/QPLIB**, in the order that bites: wire `biconj` into `SCIP/src/cca2ConvexEnvelope.m`
   → expose value+subgradient off `QuaParCPLQ` → fix diagonal terms over a box (`x²−y²`,
   `(x²+y²)/2` on `[0,1]²` still error in the second conjugation) → performance (~40–60 s/term).

## Relevant files

- `DECISIONS.md` — dead ends and refuted reasoning. The top entry has the lesson: when a gate
  refuses a construction you have independent reason to believe is right, suspect the gate.
- `TODO.md` — live action items.
- `SUPPORT_MATRIX.md` §4.1 (what the arc-vs-arc work was), §4.2 (the verification tools), §4.3
  (the covering proof).
- `maxQuaPar.m` — `newRaySign`, `assertWinnerHoldsAtInfinity`, `splitUnboundedTwoArcPiece`,
  `pieceIsCurved` (read its header before touching `curveAfter` anywhere).
- `verifyFacesCoverThePlane.m`, `verifyMaxIsExactSymbolically.m`, `pieceRecessionRays.m`.

# Morning report — 2026-08-13 overnight run

Branch: `overnight/2026-08-13`

Task as given: items 1–4 from `TODO.md` / the session handoff —
1. `arcVsArcRefusesAnUnboundedTwoArcSplit` (the last `maxQuaPar` red),
2. the covering proof for `verifyMaxIsExactSymbolically`,
3. bounded arc-pieces whose constraint region is non-compact,
4. refresh `SUPPORT_MATRIX.md` §7/§8 and `FARFIELD_FIX_PLAN.md`.

Explicit instruction: do not pause for input, and keep attacking a bug past the
usual three-strikes rule.

## Headline

`maxQuaParTest` **26 pass / 1 fail**, against a `main` baseline of **25 / 1**. One more test
passing, the same single red — `arcVsArcRefusesAnUnboundedTwoArcSplit`, which is **still red** and
which item 1 did not fix.

For the first ~half of the run MATLAB was unusable (the UBC licence server is on an internal domain
and the VPN was down, so `matlab -batch` failed with License Manager Error -96 on every attempt).
That half produced documentation and code written blind. The VPN came back later and **everything
below has now been run**; two of the blind changes were wrong and were fixed or reverted on
evidence.

## What changed

- **Item 4 (docs) — DONE.** `SUPPORT_MATRIX.md` §4's table was stale in both directions: three rows
  marked GAP are now handled, several real guards were missing from it, and every line-number
  citation predated the arc-vs-arc work by ~1400 lines, so each pointed at unrelated code. All
  re-derived from the source. §4.1's closing sentence described the ray split that was *reverted*
  rather than the code that is there. New §4.2 (the verification tools, which existed only in the
  session handoff) and §4.3 (the covering proof). §7 now says the far-field defect is closed; §8
  promotes the one remaining `maxQuaPar` case to its own entry and marks the far-field blocker
  resolved. `FARFIELD_FIX_PLAN.md` gets an OUTCOME section scoring its five phases against what
  actually happened — its own diagnosis was wrong, and that gap is the useful part of the record.
- **Item 2 (covering proof) — DONE and GREEN.** `verifyFacesCoverThePlane.m` decides, from the
  constraint data rather than from probe points, that the faces leave no hole; four checks that
  together force the boundary of the union to be empty. `maxQuaParTest/arcVsArcResultsCoverThePlane`
  passes on all three arc-vs-arc fixtures. This is the last part of `FARFIELD_FIX_PLAN.md` Phase 4
  that was resting on `partitionReport`'s sampling.
  It took two corrections after its first run: a row-vs-column indexing slip that made the
  tolerance non-scalar, and — the real one — a tolerance read off the curve **parameter** instead
  of the **point**. On an arc the parameter enters as a quartic, so the threshold reached ~0.1 and
  the check reported four far-field over-extensions on a result that is exact. Probing the ranges
  it named settled that they were false. It now uses `QuaPar.eval`'s own admission rule, which is
  the thing being decided.
- **A live bug fixed: a newly minted OUTGOING ray was given sign `+1`.** `polyConstraints` reads a
  ray's outward normal as `sign · rot90ccw(direction)`; an outgoing ray needs `−1`, and both
  branches of `splitCell` that mint an escaping ray wrote `+1`. That flips the kept half-plane to
  the far side of the ray's line, so the piece's constraint region is the reflection of its true
  region. Measured no regression.
- **Item 3 — the TODO item was stale, and what was left of it is fixed.** Its representation
  question was answered on 2026-08-13 by deriving the chord per face (`QuaPar.chordCuts`); what
  remained was that `pieceRecessionRays`, the piece-level analogue, still read the chord's side off
  the piece's other **vertices** — the rule `DECISIONS.md` already records as unsound one level up —
  and had no gate on when a chord may be emitted at all. Both are now settled by the conic itself.
  Measured no regression.

## What is broken

- **`arcVsArcRefusesAnUnboundedTwoArcSplit` is still red** — unchanged from `main`. Item 1 was
  attempted and **reverted for the second time**, on measurement; see below.
- Nothing else. No regression anywhere I ran.

## Needs a decision

**Item 1 is the one thing that needs you, and it is now much better characterised.** The ray split
was restored with an acceptance gate built from the two exact invariants (`reccConeViolation`,
`winnerDominationViolation`) instead of a seventh heuristic, and it measured as the worst of both
worlds: it still **refuses** on the pinned fixture, and it **regresses** the seeded sweep on the
same shift `[1.4979 3.6486]` by the same `0.3531` as the first attempt. So it was reverted.

What the run established, and it changes the diagnosis:

- The wrong point is admitted by exactly **one** face, which carries a zero quadratic where the
  truth is `0.35310191`. It is **not** two faces overlapping with different values — not the shape
  the far-field defect had, and not the shape `reccConeViolation` is built to catch.
- `verifyMaxIsExactSymbolically` names it exactly, and the beater is **g2 face 4** — a *different*
  g2 face from the cell's own source. That makes it a **containment** failure, and containment is
  the one invariant of the three that was left out of the gate.
- `winnerDominationViolation` cannot see it by construction: it compares the carried row against
  the other row *of the same cell*, which is the right comparison only while the piece stays inside
  both source faces.

So the question is no longer "where does the cut start" (last session's guess, now superseded) but
**whether the parent cell already fails containment, with the old guard's error having masked it.**
Completing the gate needs `facePoly(g1,k)`/`facePoly(g2,l)` inside `splitCell`, which has only
`f1row`/`f2row` — not a one-line change, which is why it stops here rather than being guessed at.

**Blocking that measurement:** with `MAXQP_ASSERT` on, this shift dies in `parabolaArcFrame` with
`degenerateAxis` — some piece carries a `curveEc` that is not a genuine parabola in its own frame.
The invariant tooling cannot currently run on the one input that most needs it. That is worth
fixing first.

## Where I stopped

Items 2, 3 and 4 are done and green. Item 1 is reverted with its diagnosis recorded in
`DECISIONS.md` and `TODO.md`. The next concrete step is the `degenerateAxis` crash, because it
blocks the measurement item 1 now depends on.

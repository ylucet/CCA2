# Session Handoff

_Last updated: 2026-08-10_

## What happened this session

Worked the arc-vs-arc `maxQuaPar` defect on branch `overnight/2026-08-02`. Fixed **two of the three**
arc-vs-arc pins earlier in the run (T-junction passthrough + `raySideVector`; `assignSide` winner),
then got the third (`twoCurvedWhereTheSplitCurveCrossesAnArc`, shift `[0.5 0.5]`) to **assemble**
via six chained fixes. Then discovered — and this is the headline — that the arc-vs-arc result is
only **locally correct**: it is right near the arcs but **WRONG in the far field**, pervasively (even
the two "fixed" pins are wrong far out; they pass only because their samples sit near the arcs).
Drafted a proof-based plan, and implemented its sound foundation.

## Where things stand

- Branch: `overnight/2026-08-02` @ `5e704ba` — "Phase 0: sound localisation of the over-extenders
  (recession cone over 18 seeded shifts)"
- Pushed: **no** — branch is ahead of `origin/overnight/2026-08-02` (~17 commits). NOT merged to
  `main`. Do not merge or push without the user's explicit say-so.
- `maxQuaParTest`: 18 pass / 2 fail. The 2 reds are `twoCurvedWhereTheSplitCurveCrossesAnArc`
  ([0.5 0.5], assembles-but-wrong) and `arcVsArcMatchesGroundTruthOverRandomShifts` (the new seeded
  far-field sweep, deliberately red — pins the defect).
- **Methodology decided (do NOT regress to sampling):** detection/verification must be SYMBOLIC
  proofs, not numerical probing. Sampling far points can miss the bad direction forever (the
  `maxArray` "decide from probes" anti-pattern). See `FARFIELD_FIX_PLAN.md` "Methodology".
- **No unrepresentable case exists.** conj/biconj of a nonconvex QuaPol is proven representable as a
  QuaPar (`proveStageA-D`). Any "needs an unbounded curved edge" / non-compact-face signal is a
  SUBDIVISION BUG, not a limitation. Do not add an "unrepresentable escape hatch".

## Next steps

- **Immediate open question (the user's last Socratic point, and it is right):** a bounded region of
  a parabolic arc + 2 segments IS a valid compact QuaPar face — so the over-extension is a
  **construction bug** (a dropped/incorrect bounding edge, or a flipped arc orientation), NOT a
  representational impossibility. My earlier "arc-triangle is non-compact / Phase 2 is a rework"
  framing was too pessimistic and was retracted.
- **Next concrete action:** verify whether the offending piece's OWN cell (straight out of
  `clipByFace`, BEFORE `splitCell`) is already non-compact. That pins the bug to `clipByFace` vs
  `splitCell`. Use the capture hook `MAXQP_CAPTURE`/`MAXQP_PIECES` (in `maxQuaPar.m`) + the sound
  primitive `pieceRecessionRays.m`. On `[0.5 0.5]`, the culprit cell is `g1f6 ∩ g2f1`: `g1f6` is a
  cone at `(-2,2)` opening UP, but the piece admits `(-3.98,0.61)` (down-left, OUTSIDE that cone) —
  so a bounding edge that should exclude the down-left region is missing.
- Then Phase 2: restore the missing bounding edge / fix orientation in the culprit producer.
  Phase 4 acceptance is SYMBOLIC (region-partition + per-cell polynomial identity `q_g ≡ max`), not
  the seeded sweep (that is a dev tripwire only).
- Reproducer is deterministic: `runtests('maxQuaParTest/arcVsArcMatchesGroundTruthOverRandomShifts')`,
  seed `20260803`, N=18, near + radius-25 ring; ~9/18 assemble to a wrong value.

## Relevant files

- `FARFIELD_FIX_PLAN.md` — the phased, proof-based plan; Phase 0 RESULTS recorded there. Start here.
- `pieceRecessionRays.m` — the sound recession-cone primitive (closed-form, exact `sym` sign tests).
- `maxQuaPar.m` — `splitCell` (incl. the session's `triHalf` two-arc split + escape-to-infinity
  branch), `assignSide`, `clipByFace`, `clipPolyByConic`, `insertGlobalPassthrough`; capture hook
  near `assemblePieces`.
- `maxQuaParTest.m` — `arcVsArcMatchesGroundTruthOverRandomShifts` + `sweepTwoCurvedShifts` (seeded,
  far-ring), and the three `twoCurved...` pins.
- `TODO.md` — "MAJOR FINDING 2026-08-04" (far-field pervasive) and the 17/2→18/2 arc-vs-arc history.

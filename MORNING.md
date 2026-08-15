# Morning report — 2026-08-15 overnight run

Branch: **`main`** (merged, and pushed as the run went — you authorised pushing for this run)

Task as given: housekeeping (merge the previous branch, turn the piece invariants on, retire the
stale TODO entry, commit and push), then fix bugs 1–5. Later, explicitly: fix bug 5. Standing
instructions: parallel agents permitted, do not stop after three failures, do not wait for input.

## Headline

**Two bugs fixed, one localised, one re-scoped.** `maxQuaPar` now has **no open case**: the seeded
arc-vs-arc sweep is **18 exact / 0 wrong / 0 errored of 18**, from 16 / 0 / 2 two sessions ago.
`maxQuaParTest` 29 / 0; fast bucket 204 / 0.

## What changed

- **Housekeeping — DONE and pushed** (`5f22486`). Previous branch merged into `main`.
  `MAXQP_ASSERT = 1` is now on for every `maxQuaParTest` test, via a `TestMethodSetup` that
  restores the previous value on teardown. The stale `twoCurvedWhereTheSplitCurveCrossesAnArc`
  entry is retired.

- **BUG 5 — FIXED** (`db4a188`). `splitTwoArcPiece` found no cut when a piece's two arcs are
  ADJACENT: its two candidate chords join the arcs' facing endpoints, which for arcs sharing a
  vertex ARE the arcs' own edges, so both chains came out too short and the piece was returned
  unsplit with one arc flattened to its chord — surfacing three stages later as an orphan
  boundary edge. Generalised the `nv == 3` shared-vertex fallback to `nv >= 4` with the ordinary
  diagonal to a non-adjacent vertex. **Sweep 17 / 0 / 1 → 18 / 0 / 0.** Pinned by
  `arcVsArcSplitsTwoADJACENTArcsOnAPieceWithADiagonal`, which exists because the sweep test
  asserts `nWrong == 0` and would have counted this input in `nErr`.

- **BUG 1 — three defects fixed; `f**` exact at 5 of 7 probe points, was 0 of 7**
  (`b105814`, `508b336`). In the order they had to be found:
  1. `getEdgeNosInf` numbers an edge by one of its endpoint VERTICES, so a LENS — two edges
     joining the SAME pair — gets one number for both and the last-write-wins scatter destroys
     one. That, not the `size(ineqs,2) == nv` count the record blamed, is the "conjugate of the
     chord".
  2. `getNormalConeVertexQ` indexed its second constraint as `j+1` unwrapped, raising
     `badsubscript` on any BOUNDED region — which is why its only caller sent every bounded region
     to the POLYHEDRAL routine, whose cones come from the CHORD and are wrong for a curved edge.
  3. `biconj` handed its second conjugation the curved MESH `conj` has returned since 2026-08-13;
     `quaPolToPlq` refuses a curved domain, so it died before reaching any of the above. It now
     asks `conjCPLQ` for the symbolic form on purpose.
  Unit level: the half-lens conjugates to 3 cells exact against a brute-force sup at all 10 probe
  points (2 identical wrong cells before).

- **BUG 2 — LOCALISED** (`1f02101`), and two suspects cleared. Step 2 is right (each piece's own
  conjugate has the correct quadrant constraints, and the per-piece max at `(-3,-2.4)` is the
  truth). The assembled cell is `s1²/4 + s2²/2` on `{s2 ≤ 0, s2²/2 − s1² ≤ 0, s1² − 2s2² ≤ 0}`:
  the sign constraint `−s1 ≤ 0` is gone and the quadratics that replaced it are **blind to the
  sign of `s1`**, so the region is symmetric and claims the mirror wedge.
  `region.redundantSubset` is exonerated — asked directly, it certifies nothing as redundant.

- **BUGS 3–4 — re-scoped, not attempted.** They are a missing ALGORITHM, not a defect:
  `convEnvUnbounded` computes only the AFFINE envelope over an unbounded face and refuses the
  rest by design, and both fixtures have envelopes that are convex and not affine.

## What is broken

The four documented slow-bucket reds, unchanged in kind:
`biconjugateOverATwoFaceSubdivisionIsTheEnvelope` (now 5 of 7 points right instead of 0),
`step3UnboundedAssemblyMatchesTheTruth`, and two in `unboundedFaceTest`.

## Needs a decision

1. **Bug 1's remainder is a refactor, and I stopped rather than force it.** The lens's two edges
   need slots 1 and 2, held by constraints that bound no edge. Freeing them by DROPPING those was
   tried and is **unsound** — a constraint active at one vertex of a convex region can still be
   essential, and removing it enlarges the piece, which SHRINKS its conjugate domain. Measured:
   with the drop, `f**` is exact at two points and `+inf` at two others; without it, exactly the
   other way round. Both 5 of 7, one of them sound. The real fix is an explicit EDGE LIST spanning
   `conjugateOfPiecePoly`, `getNormalConeEdgeQ`/`Q3` and `getSubdiffVertexT2`/`T2Q` together,
   because the loop variable indexes all of them.
2. **Bugs 3–4 belong off the bug list**, as research items.

## Where I stopped

All work is committed and pushed on `main`. `DECISIONS.md` carries the unsound approaches so they
are not re-derived; `TODO.md` and `.claude/SESSION_HANDOFF.md` carry the next concrete step for
each open item.

The transferable lesson, now in the handoff: **build a unit-level reproducer before touching the
symbolic layer.** The hand-built half-lens runs in seconds against a brute-force sup; the same
defect through the pipeline takes 10–40 minutes per run. That is what made bug 1 tractable after
two sessions of failed attempts on it.

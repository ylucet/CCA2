# Morning report — 2026-08-15 overnight run

Branch: **`main`** (merged, and pushed as the run went — you authorised pushing for this run)

Task as given: housekeeping (merge the previous branch, turn the piece invariants on, retire the
stale TODO entry, commit and push), then fix bugs 1–5. Later, explicitly: fix bug 5. Standing
instructions: parallel agents permitted, do not stop after three failures, do not wait for input.

## Headline

**Bugs 2 and 5 FIXED, bug 1 taken from 0 to 5 of 7, bugs 3–4 re-scoped and their maths derived.**
`maxQuaPar` has **no open case** (seeded sweep **18 / 0 / 0**, from 16 / 0 / 2 two sessions ago),
and **cPLQ's Step 3 cross-piece maximum is closed** — the last blocker of `SUPPORT_MATRIX.md`
§8.2. `maxQuaParTest` 29 / 0, `conjCPLQTest` 25 / 0, fast bucket 204 / 0.

The two fixes turned out to be the same lesson twice over, one per subsystem:

> **A degenerate geometric object is not a geometric object.** Bug 5: two arcs that share a vertex
> have no separating chord, because the "chords" are the arcs' own edges. Bug 2: a quadratic has no
> tangent line at its own apex, where the gradient vanishes — and `SUPPORT_MATRIX.md` §8.2(e)
> already recorded a *different* routine falling into that exact trap on that exact input.

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

- **BUG 2 — FIXED** (`1d04c75`). `region.removeTangent` builds the TANGENT LINE to a quadratic
  constraint at a vertex and deletes any constraint matching it. At the APEX OF A CONE the
  gradient vanishes, there is no tangent, and what it computes is meaningless — and that apex is
  exactly where an unbounded fan's Step 3 split conics meet. It deleted `−s1 ≤ 0`, leaving two
  constraints **blind to the sign of `s1`**, so the region went symmetric and claimed the mirror
  wedge. `f*(-3,-2.4)` is now **4.5**, the truth; it was 5.130. **`conjCPLQTest` 25 / 0.**
  Cleared on the way, so nobody re-checks them: `redundantSubset`, `simplifyUnboundedRegion`, and
  Step 2 itself. The old `step3DropsCellsOnSomeUnboundedAssemblies` pinned the GATE firing; it no
  longer does, so it is renamed and rewritten to pin what the gate protects.

- **BUGS 3–4 — re-scoped, and their mathematics derived.** They are a missing ALGORITHM, not a
  defect: `convEnvUnbounded` computes only the AFFINE envelope and refuses the rest by design.
  Both envelopes are now in `TODO.md` rather than left to the next attempt:
  `co(x·y + I_K) = y²` on `K = {0 ≤ y ≤ x}` — with its proof — and `co(−x²+y²) = −x + y²` on the
  half-strip, plus the pattern they share and a warning not to ship a formula that merely matches
  the two fixtures.

## What is broken

Three reds, down from four:
`biconjugateOverATwoFaceSubdivisionIsTheEnvelope` (now 5 of 7 points right instead of 0) and the
two in `unboundedFaceTest` (bugs 3–4, the missing algorithm).
`step3UnboundedAssemblyMatchesTheTruth` is GREEN.

The full slow bucket was still running against the final tree when this was written — check it
before merging anything on top. Every suite re-run individually against these changes was clean:
`conjCPLQTest` 25 / 0, `testMaxMultiRegion` 24 / 0, `testcPLQ` 8 / 0, `testRegion` 23 / 0,
`biconjCPLQTest` 10 / 0.

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

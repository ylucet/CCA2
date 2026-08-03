# TODO

_Seeded 2026-08-02 at the start of the overnight run, from the task given when it
was launched. The repository had no TODO.md; the acceptance criterion is precise
(three named tests green), so the run works from this list._

## Status 2026-08-03 (Opus 4.8)

Three arc-vs-arc pins were red; **now two**, and they are THREE DIFFERENT defects, not one.
`maxQuaParTest` 16/3 -> **17/2**, no regression (random-quadrilateral sweep + all arrangement
tests stay green). Committed `96aad61` on `overnight/2026-08-02`.

- **[-1 0.75] `twoCurvedThatMustAssembleAcrossRays` -- FIXED.** Was a cross-piece T-JUNCTION:
  a decided cone's ray ran to infinity while the neighbour side was split (segment+ray) at a
  point P lying exactly on that ray (perp ~2e-15). `insertGlobalPassthrough` now re-inserts
  every piece vertex that lands on another piece's ray/segment; a companion fix to
  `raySideVector` (its old adjacent-vertex representative is COLLINEAR with a just-subdivided
  ray, so `oppositeSides` decided on sign noise -- now falls back to the other ray's direction).
  Green and exact at all 68 ground-truth samples.

- **[0.5 0.5] `twoCurvedWhereTheSplitCurveCrossesAnArc` -- STILL RED (errors).** This is the
  ORIGINAL piece-5 arc-clip defect below, on the `clipPolyByConic` path, NOT a T-junction (the
  post-passthrough scan finds nothing collinear with piece 5's ray). The analysis below stands.

- **[2 -0.5] `twoCurvedAssembleAcrossRaysSecondFixture` -- STILL RED, now a WRONG VALUE not an
  error.** The T-junction fix lets assembly COMPLETE, which uncovers a coverage/winner defect:
  `g` returns the wrong operand over a region -- 4/68 samples, e.g. `s=(-5.4843,1.5866)` gives
  `-8.648` (=g2) where truth is `0` (=g1), and `s=(-1.1298,4.1007)` gives `0` (=g1) where truth
  is `1.768` (=g2). Two regions with swapped winners. All accepted RAY matches are geometrically
  correct adjacent-`src` pairs, and the passthrough only inserts collinear vertices (never touches
  the winner row), so this is a SEPARATE coverage defect the completing assembly newly reaches --
  likely a wrong FACE built in `buildFinalEdgesAndFaces`, or a cell whose winner/region is wrong.
  Next: probe the 4 bad points against the pieces PRE-assembly (does the correct-winner piece
  cover them there?) to localise assembly-vs-decideWinner.

## Next up (the [0.5 0.5] defect)

- [ ] **Piece 5 (src `[2 2]`) emits a RAY where its boundary should terminate.**
      Localised to the cell, the edge and the reason; this is the last defect.

      * Its unmatched ray: apex `(-2.03125, 2.03125)`, direction `(-1,1)`, lying
        on the line `x+y=0`, which is g1's face-2/face-6 edge.
      * CORRECTION to an earlier note here: piece 16's ray is NOT its partner.
        That one lies on `x+y=0.5`, a parallel but different line, so the two were
        never meant to match.
      * Sampling across the apex, the local structure is three cells:
        `(2,2) | (6,1) | (6,2)`, with g2's ARC separating the last two. So piece
        5's neighbour along `x+y=0` is a `(6,1)` cell — and the `(6,1)` pieces
        (13 and 14) are BOUNDED slivers of area 0.008 and 0.004.
      * A bounded neighbour means the shared boundary is a finite SEGMENT, not a
        ray. So piece 5 is over-extended: its boundary along `x+y=0` should stop
        where g2's arc crosses that line a second time, and instead it runs to
        infinity. `matchHalfEdges` pairs rays with rays and segments with
        segments, so a ray facing a segment can never match — which is exactly
        the reported symptom.

      Next: find why `clipPolyByConic` cuts cell `(2,2)` at the first crossing
      with g2's arc but not at the second. Note the restriction of that conic to
      this ray came out with `A = 1.7e-18` — numerically degenerate, so the
      quadratic is treated as linear and yields ONE root. Check whether the true
      second crossing is being lost there, or whether it lies on a different
      boundary element of the cell.

## Done recently

- [x] `clipPolyByConic`: replaced the blanket refusal of a disconnecting curved
      cut with a real CONNECTIVITY TEST. Measured CONNECTED on all three
      fixtures, so the single-cell construction was right and "return both
      components" was not the fix — see Blocked note below for why two
      components would be unrepresentable anyway.
- [x] The clip and split stages now produce a VALID PARTITION on all three
      fixtures (`partitionReport` OK), where before they overlapped.
- [x] Fixed the no-crossing keep/drop decision, which evaluated the conic at the
      centroid of the vertices — not necessarily inside a cell bounded by an
      inward-bulging arc.
- [x] Retracted the hole/overlap evidence and made the partition diagnostic
      sound (it omitted the arc edge, so a curved sliver appeared to contain a
      point 3 units away).
- [x] Tagged pieces with their source `(k,l)` pair.

## Blocked

- Nothing is blocked on a decision. The note that "return both components" is
  unrepresentable stands as a finding, not a blocker: a separated survivor would
  need the cutting conic running to infinity, i.e. an unbounded curved edge,
  which QuaPar cannot hold. The connectivity test now refuses only that case,
  and it does not arise on these fixtures.

## Retired hypotheses — do not re-try

- The orphan ray is one physical ray covered to different extents (nothing lies
  on it).
- The cut must be restricted to the arc's own span (it must not; the straight
  half-planes are applied first).
- `clipPolyByConic` emits clockwise cells (a guard on every bounded output never
  fires).
- `polyL` is non-convex (all four curved faces measured convex).
- Pair `(6,1)` is spurious (it is a genuine thin cell; a 0.0625-spaced grid
  simply cannot see a cell of area 0.008 — the "zero intersection" reading was a
  resolution artifact, not evidence).

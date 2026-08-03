# TODO

_Seeded 2026-08-02 at the start of the overnight run, from the task given when it
was launched. The repository had no TODO.md; the acceptance criterion is precise
(three named tests green), so the run works from this list._

## Status 2026-08-03 (Opus 4.8) -- 18/1

Three arc-vs-arc pins were red; **TWO are now fixed**, they were THREE DIFFERENT defects.
`maxQuaParTest` 16/3 -> **18/1**, no regression (random-quadrilateral sweep + all arrangement
tests stay green). Commits `96aad61` (T-junction) and `53fc9fd` (assignSide) on
`overnight/2026-08-02`.

- **[-1 0.75] `twoCurvedThatMustAssembleAcrossRays` -- FIXED.** Was a cross-piece T-JUNCTION:
  a decided cone's ray ran to infinity while the neighbour side was split (segment+ray) at a
  point P lying exactly on that ray (perp ~2e-15). `insertGlobalPassthrough` now re-inserts
  every piece vertex that lands on another piece's ray/segment; a companion fix to
  `raySideVector` (its old adjacent-vertex representative is COLLINEAR with a just-subdivided
  ray, so `oppositeSides` decided on sign noise -- now falls back to the other ray's direction).
  Green and exact at all 68 ground-truth samples.

- **[2 -0.5] `twoCurvedAssembleAcrossRaysSecondFixture` -- FIXED (`53fc9fd`).** Was NOT decideWinner
  (that correctly flags the cell undecided) -- it was `assignSide` reading the winner at `piece.V(2,:)`,
  which for splitCell's UNBOUNDED "rest" piece (curveAfter=1) is a CROSSING point where diff~0 and its
  sign is noise, so both halves could get the same winner. Now reads the vertex farthest from
  `{diff=0}`. Exact at all 68 samples. (The long note below about decideWinner/parabola-to-infinity was
  the WRONG lead -- the real cell here is a strip cut by a LINE at two finite boundary points.)

- **[0.5 0.5] `twoCurvedWhereTheSplitCurveCrossesAnArc` -- STILL RED (errors). NEW, sharper
  diagnosis; the old "piece 5 over-extended" note below is WRONG.** Piece 5 (`src[2 2]` =
  g1f2 n g2f2, an **arc-vs-arc** clip via `clipPolyByConic`) emits an unmatched ray on `x+y=0`
  (g1's face-2/face-6 edge) from apex `(-2.03125,2.03125)`, dir `(-1,1)`. The clip is CORRECT, not
  over-extended: the sign data at that cell is `evalConic(Ecut)@V = [0, -0.046, -0.015]`, so the
  vertex `(-2,2)` sits on the g2-face-1 (discard) side and the far ray is genuinely g2-face-2. What
  is missing is the MIRROR piece across `x+y=0`: `src[6 2]` = g1f6 n g2f2. g1 face 6 does NOT carry
  the arc (the arc edge borders g1 faces 2 and 1 only), so `src[6 2]` is clipped by the
  **arc-vs-HALF-PLANE** path (`clipPolyHalfPlaneCurved`), a DIFFERENT code path -- and it does not
  produce the matching `x+y=0` ray (its `src[6 2]` pieces 15/16 sit at `x+y=0.5` and `1.0`, and
  piece 15 has a bounded vertex exactly at `(-2.03125,2.03125)` but no ray along `(-1,1)`). So piece
  5's ray has no neighbour because the two clip paths DISAGREE along the shared g1 mesh edge.
  Post-passthrough T-junction scan finds nothing collinear with the ray, confirming it is not a
  subdivision miss. FIX DIRECTION (not yet done -- delicate, high regression risk to the arc
  machinery): make the arc-vs-half-plane clip of g2f2 by g1f6 produce the same `x+y=0` boundary the
  arc-vs-arc clip of g2f2 by g1f2 does, i.e. reconcile `clipPolyHalfPlaneCurved` with
  `clipPolyByConic` along a shared straight g1 edge.
  DECISIVE CHECK NOW DONE (confirms "f6 clip under-producing"): g1's mesh for [0.5 0.5] has
  `E(6,:)=[2 7 0]` a RAY from V2=(-2,2) in dir (-0.707,0.707) lying on `x+y=0`, with `F(6,:)=[2 6]`
  -- so it is exactly the g1 face-2/face-6 boundary and it runs to INFINITY. Therefore `src[6 2]`
  (=g1f6 n g2f2) must have a matching `x+y=0` ray, and `clipByFace(facePoly(g1,6), facePoly(g2,2))`
  -- which swaps to clip the curved g2f2 by g1f6's straight half-planes -- is dropping it. Start in
  `clipPolyHalfPlaneCurved` / `clipByFace`'s swap path: the g2f2 boundary ray on `x+y=0` is being
  lost when g2f2 (curved) is clipped by g1f6's `x+y=0` half-plane (edge 6). g1 face 6 is a cone at
  (-2,2) between edge 6 (dir (-0.707,0.707), on x+y=0) and edge 7 (dir (0.707,0.707)).

  --- OBSOLETE (kept for history; superseded by the sharper diagnosis above) ---
  The old `decideWinner`/parabola-to-infinity write-up for [2 -0.5]. The T-junction
  fix lets assembly COMPLETE, uncovering the real defect: `decideWinner` wrongly declares an
  UNBOUNDED cell "decided". Traced with a pre-assembly coverage probe: cell `src[4 3]` stores
  `winner=g2`, but at `s=(-5.4843,1.5866)` -- an interior point of that cell -- `g1=0 > g2=-8.648`,
  so g1 wins there; symmetrically cell `src[6 5]` stores `winner=g1` while g2 wins at
  `s=(-1.1298,4.1007)`. 4/68 samples, worst 8.648.
  WHY `decideWinner` misses it: it proves domination by sampling finite vertices + the ARC midpoint
  + the ASYMPTOTIC sign along the two bounding rays (t->inf). But `diff = f1-f2` is a PARABOLIC
  quadratic (rank-1 Q -- splitCell already asserts this), so `{diff=0}` is a parabola whose two
  branches run to infinity ALONG Q's null direction, which lies strictly BETWEEN the cell's two
  bounding rays. That parabola bounds a region of the opposite winner entirely inside the cell,
  touching neither the finite boundary nor either bounding ray -- so every point `decideWinner`
  samples is on one side and it "decides" wrongly.
  Tried and REVERTED (kept out of the commit): sampling `diff` at the 1-D stationary point of each
  edge and each ray. Sound and strict, but it does NOT fix this -- the sign change is off the
  1-skeleton entirely. Confirmed neutral on all three fixtures.
  Note this cell is ALSO beyond `splitCell`: even if `decideWinner` returned undecided, `{diff=0}`
  here makes ZERO finite boundary crossings (both parabola branches escape to infinity inside the
  recession cone), and `splitCell` asserts exactly two. So a real fix needs BOTH a sound
  sign-over-the-cone test in `decideWinner` AND a `splitCell` that can cut a cell along a parabola
  that enters and leaves at infinity. Same family as next-step 2 in the session handoff (the
  Step-3 unbounded over-claim); do not attempt a probe-based patch.

## Next up (the [0.5 0.5] defect) -- SUPERSEDED, read the Status section above first

The "over-extended, should terminate" framing below is REFUTED by the sign data (see Status):
piece 5's ray is legitimately g2-face-2; the real bug is the missing g1f6 n g2f2 mirror from the
arc-vs-half-plane clip path. Kept below only for its geometric measurements.

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

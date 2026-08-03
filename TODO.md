# TODO

_Seeded 2026-08-02 at the start of the overnight run, from the task given when it
was launched. The repository had no TODO.md; the acceptance criterion is precise
(three named tests green), so the run works from this list._

## Next up

- [ ] **Two cells disagree about where their SHARED ray starts.** Fully localised;
      this is the last defect.

      * Orphan A: piece 5, src `[2 2]`, apex `(-2.03125, 2.03125)`, dir `(-1,1)`.
      * Orphan B: piece 16, src `[6 2]`, apex `(-1.7578, 2.2578)`, dir `(-1,1)`.
      * Same direction, and the ray separates cell `(2,2)` from cell `(6,2)`, so
        these two ARE each other's intended partner — they simply disagree on the
        apex by 0.355, far above `matchHalfEdges`' 1e-3 tolerance.
      * Which is right: walking the line, the face-pair transition
        `(1,1)/(5,1) -> (2,2)/(6,2)` happens at `(-2.03125, 2.03125)`. So **piece
        5 is correct and piece 16 starts 0.355 too late** — cell `(6,2)` is
        missing the stretch between the two apexes.
      * Why they can differ: the two are built by DIFFERENT paths. `(2,2)` is
        both-curved and goes through `clipPolyByConic`; `(6,2)` is swapped and
        clipped by half-planes only (`cutConic` is empty once `polyL` is the
        straight face), so its apex comes from `clipArcByHalfPlane` instead.

      Next: find why the half-plane path stops cell `(6,2)` at the sliver corner
      `(-1.7578, 2.2578)` rather than at the true crossing of g1's face-2/6 edge
      with g2's arc.

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

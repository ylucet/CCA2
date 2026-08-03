# TODO

_Seeded 2026-08-02 at the start of the overnight run, from the task given when it
was launched. The repository had no TODO.md; the acceptance criterion is precise
(three named tests green), so the run works from this list._

## Next up

- [ ] **`clipPolyByConic` misses a crossing and keeps the cell whole.** This is
      the single remaining defect and it is now located exactly, not inferred
      from the orphan ray three stages downstream:

      * The orphan ray separates cell `(2,2)` from cell `(6,2)`; both exist.
      * Walking across the ray, the true face-pair transition happens at
        `(-2.03125, 2.03125)`, so piece 5's apex is CORRECT and its neighbour's
        `(-2, 2)` is over-extended by 0.044.
      * `(-2, 2)` is `g1.V(2)`, an ORIGINAL mesh vertex; `(-2.03125, 2.03125)` is
        a vertex of NEITHER mesh — it is where a g1 edge crosses g2's ARC, i.e. a
        genuine arrangement vertex that only the conic cut can produce.
      * Traced: `clipPolyByConic` runs 4 times on this fixture. Twice it finds 2
        crossings and cuts. **Twice it finds 0 crossings on an unbounded 2-vertex
        cell and keeps it whole** — those are the cells that stay extended to
        `g1.V(2)`.

      So the fix is in crossing DETECTION, not in the surgery: find why the
      arc-vs-boundary crossing at t ≈ 0.044 along a ray is not detected. Prime
      suspects, in order: the sign-change filter in `conicRootsAlong` rejecting
      it, and the `hasArc && p == cePair` skip firing on the wrong pair.

- [ ] Only after that: `assemblePieces`. The orphan ray should disappear once
      both neighbours share the true apex — it is a symptom, not a second bug.
- [ ] Once assembled: confirm the three arc-vs-arc tests against the pointwise
      max, then the fast and normal buckets for regression.

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

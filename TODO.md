# TODO

_Seeded 2026-08-02 at the start of the overnight run, from the task given when it
was launched. The repository had no TODO.md; the acceptance criterion is precise
(three named tests green), so the run works from this list._

## Next up

- [ ] Decide how to handle a curved cut that DISCONNECTS an unbounded cell — see
      **Blocked** below. Every remaining item is downstream of this decision.
- [ ] `splitCell`: whatever the disconnection decision surfaces there.
- [ ] `assemblePieces`: orphan ray / half-edge matching.
- [ ] All three arc-vs-arc tests green.
- [ ] No regression: `maxQuaParTest` full, then the fast and normal buckets.

## Done recently

- [x] Made the partition diagnostic sound, and RETRACTED the hole/overlap
      evidence it had produced — `pieceContainsPt` omitted the arc edge, so a
      3-vertex curved sliver was reported as containing a point 3 units away.
- [x] Fixed a real bug in `clipPolyByConic`'s no-crossing decision: it evaluated
      the conic at `interiorSample`, the centroid of the vertices, which is not
      necessarily inside a cell bounded by an inward-bulging arc. It now decides
      from the boundary signs.
- [x] Tagged pieces with their source `(k,l)` face pair, which is what made an
      offending cell's provenance readable.

## Blocked

- **"Return both components" may be the wrong fix, and that is a decision for
  Yves.** If the curved cut genuinely disconnects an unbounded cell, each
  component contains one of the two rays AND is bounded in part by the cutting
  conic running out to infinity. That is an unbounded curved edge, which QuaPar
  cannot represent at all (`assertCurvedEdgesAreArcs` requires bounded arcs). So
  returning two components is not merely unimplemented — for the disconnected
  case it is not representable, and the current loud refusal is the correct
  behaviour rather than a gap.

  What is NOT yet established is whether the configuration reached on these three
  fixtures is genuinely disconnected, or is connected and simply built wrongly by
  a branch copied from the half-plane path. Deciding that is the next step, and
  it is a measurement (sample the survivor, test connectivity between the two ray
  ends), not a code change.

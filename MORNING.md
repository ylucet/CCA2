# Morning report — 2026-11-01 overnight run

Branch: overnight/2026-11-01

Task: (1) cut Example 2's overlay face count; (2) step d — cut triangles after the overlay,
biconjQ subdivides Example 2, slow test against EnvelopeMin on 4000 points.

## What changed
- Item 1 (6467b6e): Example 2's overlay sheds about a third of its faces — first 30 of 302
  tuples 1236 -> 852 faces, values identical (6000 points). Three certificates: tie rows proved
  false, faces below the max of the others in a pair's cell, algebraic rows on final pieces.
  Pairs 43/134/62 -> 42/98/59 faces. Exact polynomial divisions filtered in doubles (5476 -> 70
  on one tuple) — that was 60% of the overlay's time.
- Item 2 (b856958): triangles are cut AFTER the overlay (no cut lines needed); the output type
  holds the degree-4 triangle; biconjQ subdivides pairs + triangles. Fixed on the way: a piece
  split for one triangle was never checked against the next.
- A piece budget on the overlay: Example 2 (11567 pieces > 4000) is declined by name, so
  biconjQ still returns the EnvelopeMin by default; a new verylong test lifts the budget.
- VERIFIED: with the budget lifted, Example 2's biconjugate comes back SUBDIVIDED (AlgCrv, the
  degree-4 triangle as its own face), equal to the EnvelopeMin on 4000 points and passing the
  definition check -- erratumBiconjEx2SubdivisionTest 3/0 in 54 min (shared machine).
- AlgCrv.eval runs face by face over all points (f4ecf2f), needed for thousands of faces.
- Speed (fa541f1): two heavy Example 2 tuples 780 s -> 147 s, same faces (single runs, shared
  machine): quadratic fibres solved vectorised, line pairs and row polynomials remembered.
- Example 2's WHOLE overlay now: 3299 faces in ~21 min (one run, shared machine). Yesterday's
  projection was ~12000 faces and ~4 h.
- Tests on this branch: fast 814/0, normal 32/0, slow 97/0 (+2 known-incomplete legacy tests),
  subdivision slow suites and slow Example 2 suite 19/0, verylong Example 2 subdivision 3/0.

## What is broken
Nothing known.

## Needs a decision
None. DECIDED 2026-11-02 by Yves: keep declining Example 2 by default (option a); branch merged to main.

- (resolved) Example 2 by default: subdivide (~30 min in biconjQ now) or decline (current)? I kept DECLINE
  with a piece budget, reversible by one number. Options: (a) keep declining [my pick until the
  tangent faces are cut]; (b) raise the default budget so biconjQ subdivides it (the slow
  Example 2 suite would take ~30 min more); (c) opt-in only, skipping the ~3 min of pair
  subdivisions spent before the refusal.

## Where I stopped
Both items done; the run is finished. Left for item 1: faces TANGENT to the winner (most of what
Example 2's pair cells still keep). The exact certificate is derived in DECISIONS 2026-11-01 but
needs a big-integer square-free part of a degree-20 polynomial -- not started. Branch not merged
or pushed.

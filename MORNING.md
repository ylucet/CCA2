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
- Speed (fa541f1): two heavy Example 2 tuples 780 s -> 147 s, same faces (single runs, shared
  machine): quadratic fibres solved vectorised, line pairs and row polynomials remembered.
- Tests: fast 814/0, normal 32/0 on fa541f1; subdivision slow suites 15/0 and the slow
  Example 2 suite 4/0 on b856958.

## What is broken
Nothing known.

## Needs a decision
- Example 2 by default: subdivide (hours) or decline (now)? I kept DECLINE with a budget,
  reversible by one number. Options: (a) keep the budget [my pick]; (b) make subdividing opt-in
  only and skip even the ~3 min of pair subdivisions before the refusal; (c) lift it once the
  overlay is fast enough.

## Where I stopped
Running the verylong Example 2 subdivision test (budget lifted) from a snapshot of fa541f1, plus
a full 302-tuple overlay count. Left for item 1: faces TANGENT to the winner (most of what remains);
the exact certificate is derived in DECISIONS 2026-11-01 but needs big-integer bivariate gcds.

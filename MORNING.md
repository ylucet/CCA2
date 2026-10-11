# Morning report — 2026-10-10 overnight run

Branch: overnight/2026-10-10 (from main @ 08f1391)

## What changed
- **biconjQ now computes bounded inputs with a non-convex piece that curves up along an edge**
  (it used to refuse). New `biconjAtoms`: split each non-convex piece into its vertices and convex
  edges (erratum Thm 1), take the min over atom-pair envelopes and triple tangent planes. Erratum
  Example 1 exact in all 4 cells; matches the legacy triangle envelope to 1e-14 on 4 triangles.
  Result type is a new `EnvelopeMin` (eval only, not meshed) — see Needs a decision.
- Three bugs fixed on the way (each red first): whole-plane conjugate cells were skipped in pair
  enumeration; rational bridge faces had no membership test; corners where a cell facet meets a
  parabolic tie once were all dropped (made the envelope +inf).
- envelopeEval: a point's value depended on the other points in the call (tolerance scaled by
  the batch); fixed per point, red first.
- Exact-sum audit: 49 more integer sums that could silently get the wrong sign past 2^53 now
  refuse instead (17 found by a widened lint, 32 hidden in matrix/dot products).
- Two pre-existing envelope defects found on the poly2 kink (both made values WRONG, not just
  refused): a cell pair tying on irrational parallel lines got no corners (envelope +inf on 64%
  of its domain); and every corner face was evaluated through a branch index that picks the
  OTHER corner's value past a crossing line (0.23 too high). Both fixed, red first.
- Two corners of the poly2 kink that had no region rows now have them (a wedge edge that is a
  rational line times an irrational number). Verified against the corner's own KKT test.
  Consequence: the kink now reaches the global mesh, which has 85 curves and did not finish in
  90 min; the mesh now refuses past 40 curves by name and biconjQ returns the evaluable face list
  (verified against the definition), exactly as before.
- Erratum Example 2's biconjugate gets much further: three overflows on intermediate numbers
  fixed (each red first, from the captured input) and a witness search that only looked in a
  small box around the origin replaced by a closed form. It now stops at an honest limit: one
  face's polynomial genuinely needs 59 bits, past what double storage holds.
- **biconjQ now handles unbounded domains on the atom route** (erratum Case E). Erratum Example 3
  comes out exactly x1^2 (its quarantine is lifted), and a cone fixture refused since September
  is exactly its closed form. Both are values no finite convex combination attains.
- Daily slow gate: 101 pass / 0 fail (2 known quarantined).

## What is broken
- Erratum Example 2 biconj: one algebraic face's stored polynomial needs 59 bits — 3 tests
  quarantined on that. Needs wider storage for implicit polynomials (TODO).
- poly2 kink biconj is still not MESHED (evaluable only): the single global arrangement has 85
  curves. Fix is a per-cell-pair arrangement (TODO).
- Unbounded inputs with a CONVEX unbounded piece are still refused by name (TODO).

## Needs a decision
- **TODO item "biconjCPLQ scaling on the poly2 kinks" (legacy symbolic route, runs > 4 h) was
  skipped.** The exact route is its replacement, and tonight's fixes go there. Options: (a) drop
  the item and retire the two quarantined legacy tests in favour of the exact route's tests
  [my pick]; (b) keep it for a dedicated session; (c) keep the legacy tests quarantined forever.
- **Return type of biconjQ on the new route.** It returns `EnvelopeMin` (a min of face lists,
  answers `eval`), not an AlgCon node. Options: (a) keep it as a separate type until meshing lands
  [what I did — reversible]; (b) mesh before returning (needs envelopeMesh on several pairs);
  (c) return it only under an explicit engine flag. I'd keep (a).

## Where I stopped
Row-less corners done; Example 2 integer width scope-reduced to the stored-polynomial width; unbounded atom route done.

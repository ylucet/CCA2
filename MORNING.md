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
- Daily slow gate: 101 pass / 0 fail (2 known quarantined).

## What is broken
- Erratum Example 2 biconj: `ratQ:overflow` (integers past 2^53) — 3 tests quarantined. Needs a
  wider number type in the face constructors (TODO).
- Erratum Example 3 biconj: unbounded, atom route is bounded-only — 1 test quarantined (TODO).

## Needs a decision
- **Return type of biconjQ on the new route.** It returns `EnvelopeMin` (a min of face lists,
  answers `eval`), not an AlgCon node. Options: (a) keep it as a separate type until meshing lands
  [what I did — reversible]; (b) mesh before returning (needs envelopeMesh on several pairs);
  (c) return it only under an explicit engine flag. I'd keep (a).

## Where I stopped
Item 1 done (scope-reduced); regression buckets running before commit.

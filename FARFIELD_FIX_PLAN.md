# Plan: fix the arc-vs-arc far-field over-extension

_Drafted 2026-08-04. Companion to TODO.md "MAJOR FINDING 2026-08-04"._

## The bug in one paragraph

`maxQuaPar(g1,g2)` for two curved operands is correct NEAR the arcs and WRONG in the far field.
Root: an UNBOUNDED quadratic conjugate face (e.g. g1 face 1 carries rays) gets emitted by the
arc-vs-arc subdivision as a BOUNDED arc-piece — one parabola arc plus straight edges, no rays.
The piece sits on the parabola's OPEN side, so it is not a compact QuaPar face: `QuaPar.eval`
(locate by "every bounding conic, sign-oriented, ≤ tol") admits points arbitrarily far away, the
piece overlaps the true face out there, and eval's last-admitter-wins returns the wrong value.
The fix is UPSTREAM: a sub-piece that inherits an unbounded direction must keep the RAY edge that
bounds it; it must not be closed off with a straight chord.

## Reproducer (seeded, deterministic)

- `runtests('maxQuaParTest/arcVsArcMatchesGroundTruthOverRandomShifts')` — seed 20260803, N=18,
  samples near the arcs AND a radius-25 ring. Currently **9/18 assemble to a wrong value, worst ~58**.
  The failure message lists the exact wrong shifts.
- Single-shift check: build `[g1,g2]=buildCurvedG1G2(T1,T1+shift)`, `g=maxQuaPar(g1,g2)`, compare
  `g.eval(x)` to `max(g1.eval(x),g2.eval(x))` on rings of radius 8 and 30. Confirmed wrong on the
  three "pinned" shifts too: [-1 0.75] 2/60 @ R=30, [2 -0.5] 11/60 @ R=30, ~0 @ R=8.

## Phase 0 — Localise which subdivision path drops the rays (measure, don't guess)

For each wrong shift from the seeded sweep, on the ASSEMBLED result:
1. Find the over-extending face: probe far points; the face that admits one AND disagrees in value
   with another admitter there is the culprit (the prototype in this session's git history —
   `assertNoOverExtendedFace`, reverted — does exactly this using QuaPar's own EC + P signs; lift it
   back as a DIAGNOSTIC, not a hard error).
2. Map that face back to its PIECE (g face i ↔ piece i, since `g.f(p,:)=pieces(p).f`) and its `src`
   `(k,l)` and `curveAfter`.
3. Tag which producer emitted it: `clipByFace` (arc-vs-half-plane or the straight path),
   `clipPolyByConic` (arc-vs-arc), or `splitCell` (which branch: bounded, unbounded-escape-to-inf,
   `splitTwoArcPiece`, or the `triHalf` two-arc-triangle fallback added this session).
Output: a table {shift → culprit piece src, curveAfter, producer}. Expect a small number of
producers. On [0.5 0.5] the culprit was a `triHalf` piece (src[6 1]); the pinned shifts use OTHER
producers, so this MUST be done across several shifts, not just [0.5 0.5].

## Phase 1 — State and check the invariant

INVARIANT (the thing every producer must preserve): a piece is a valid QuaPar face iff its
constraint region {arc side} ∩ {straight half-planes} ∩ {ray half-planes} is COMPACT, or it is
legitimately unbounded via ray edges. Equivalently: a piece bounded partly by an arc must lie on
the arc's CLOSED (convex) side, OR carry the ray(s) that bound the open side.

Add `pieceIsValidFace(piece)` (reliable version = probe far points using the SAME orientation the
QuaPar constructor will give it; the raw-piece sign convention got this wrong this session — build
the tiny QuaPar for the single piece and reuse `assertNoOverExtendedFace`'s logic, or normalise
signs exactly as `buildFinalEdgesAndFaces` does). Assert it on every piece right before
`assemblePieces`. This turns silent-wrong into a LOUD, localised failure and is the regression gate
for the rest of the plan.

## Phase 2 — Fix ray propagation, producer by producer

Hypothesis to confirm in Phase 0, then fix wherever it holds: when a producer intersects/splits a
region that is UNBOUNDED in some direction (because one operand's face is an unbounded quadratic
face, i.e. `facePoly` returned rays), the resulting sub-piece that still contains that direction
must inherit the corresponding RAY, not be closed with a straight edge.

Likely fixes, in rough order of suspicion:
- **`triHalf` / the two-arc-triangle split** (added this session): it cuts from the shared vertex to
  the opposite-edge MIDPOINT, producing bounded triangles. If the parent two-arc region is
  unbounded, the correct cut is a RAY to infinity along the recession direction, not a chord to a
  midpoint. Rework `triHalf` to emit an unbounded half (with a ray) when the parent has rays.
- **`splitCell` bounded/unbounded branches**: verify `cellRest`/`cellMid` keep `dirIn`/`dirOut`
  whenever the corresponding side is unbounded; the crossed-arc restoration (`splitTwoArcPiece`)
  must not flatten an unbounded side into a chord.
- **`clipPolyByConic`**: the arc-vs-arc clip builds `mkPoly(...)`; confirm it carries the rays for an
  unbounded survivor and never emits an arc-piece on the open side without them.
- **`clipByFace` arc-vs-half-plane path** (`clipPolyHalfPlaneCurved`): same check.

For each, the concrete question is: "does this construction ever produce a bounded arc-piece from a
parent that had rays in the arc's open direction?" If yes, carry the ray.

## Phase 3 — The genuinely-unrepresentable escape hatch

Some intersections may need an UNBOUNDED curved edge (a parabola running to infinity as a face
boundary), which QuaPar cannot represent (`assertCurvedEdgesAreArcs`). Where Phase 2 cannot produce
a compact-or-ray-bounded piece, ERROR LOUDLY (`maxQuaPar:notImplemented`) rather than emit a wrong
face. `conj`'s verification then turns it into a clean failure. Distinguish this from Phase 2 bugs
by whether the recession direction is bounded by a straight edge (fixable) or only by the conic
itself (unrepresentable).

## Phase 4 — Lock it in with tests

- Keep the seeded far-field sweep (already committed); target it to **0 wrong** (only exact or loud
  error).
- Add a radius-25/30 far ring to the THREE existing arc-vs-arc pins
  (`twoCurved...`) so they can never again pass while being wrong far out.
- Keep `pieceIsValidFace`/`assertNoOverExtendedFace` as a permanent assertion (gated so it is cheap)
  so any future producer that drops a ray fails loudly at the source.

## Decision point (raise before Phase 2)

Incremental (patch ray-propagation in each producer, guarded by the Phase-1 assertion) vs. a rework
of the arc-vs-arc assembly. Recommend INCREMENTAL: the invariant + assertion make each producer's
fix independently verifiable, and Phase 0 should show only a few producers are at fault. Escalate to
a rework only if Phase 0 shows the ray-loss is spread across most producers.

## Do NOT repeat (tried this session, reverted — see git log around 0b5eee8)

- Raw-piece compactness guards (`arcFaceIsBounded`): got the arc/half-plane sign convention wrong and
  rejected valid faces. Any compactness check must use the ASSEMBLED QuaPar's own EC + P signs.
- A hard post-assembly "no non-compact face" error: correct, but errors on ~every arc-vs-arc result
  because the defect is pervasive — it is a diagnostic/gate, not a shippable guard, until Phase 2
  lands.

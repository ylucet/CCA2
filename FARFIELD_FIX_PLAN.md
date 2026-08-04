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

## Methodology: PROVE, do not sample

Numerical sampling (probing far points, rings, the seeded sweep) is EVIDENCE, never proof — one can
sample forever and miss the one direction a face over-extends. It is the same "decide from probes"
anti-pattern the handoff removed from `maxArray`. So detection, the invariant, and acceptance are all
SYMBOLIC; sampling is only a cheap tripwire during development, never the acceptance test.

- **Boundedness / over-extension is decidable EXACTLY via the recession cone.** A face is a region
  `{ sign_i * evalConic(c_i, .) <= 0 }`. Direction `d` is recessive iff for every bounding edge the
  LEADING term of `sign_i * evalConic(c_i, x0 + t d)` as `t->inf` is `<= 0`: for a line
  `sign*(a dx + b dy) <= 0`; for a parabola `sign*Q(d) <= 0` on the quadratic form, tie-broken to the
  linear term where `Q(d)=0`. A face with NO ray edge is valid iff its recession cone is `{0}`; a
  face WITH rays is valid iff its recession cone equals exactly the cone those rays span. The
  over-extension is then PROVED ("this arc-triangle has recession direction along the parabola axis
  yet carries no ray"), not sampled. Build on `region`'s recession-direction machinery (rational
  half-planes) and `isAlways`, not on a far-point probe.
- **Correctness `g = max(g1,g2)` is proved REGION BY REGION.** Overlay the three subdivisions into
  common cells (finite). On each cell `C` the overlay's switching curve `{g1=g2}` fixes which operand
  dominates, so `g = max(g1,g2)` on `C` reduces to a POLYNOMIAL IDENTITY `q_g ≡ q_{g_i}` on `C`
  (exact: `isAlways` / coefficient equality), plus an exact proof the cells PARTITION the plane (no
  gap, no overlap) via region algebra. Finitely many proved identities cover ALL points at once —
  unlike any ring of samples.

## Reproducer (seeded numerical tripwire ONLY — not the proof)

- `runtests('maxQuaParTest/arcVsArcMatchesGroundTruthOverRandomShifts')` — seed 20260803, N=18,
  samples near the arcs AND a radius-25 ring. Currently **9/18 assemble to a wrong value, worst ~58**.
  The failure message lists the exact wrong shifts.
- Single-shift check: build `[g1,g2]=buildCurvedG1G2(T1,T1+shift)`, `g=maxQuaPar(g1,g2)`, compare
  `g.eval(x)` to `max(g1.eval(x),g2.eval(x))` on rings of radius 8 and 30. Confirmed wrong on the
  three "pinned" shifts too: [-1 0.75] 2/60 @ R=30, [2 -0.5] 11/60 @ R=30, ~0 @ R=8.

## Phase 0 — Localise which subdivision path drops the rays (PROVE, don't probe)

For each wrong shift (the seeded sweep only tells us WHICH shifts are wrong quickly; the localisation
itself is symbolic), on the pieces feeding `assemblePieces`:
1. Compute each piece's recession cone `R` in closed form (see Methodology). The culprit is a piece
   with NO ray edge but `R != {0}` (over-extended), or a piece WITH rays whose `R` is strictly
   larger than the cone its rays span (dropped a ray). This is a PROOF the piece is malformed, not a
   far-point probe.
2. Independently, that malformed piece's recession cone `R` should equal the recession cone of its
   TRUE region = `recc(facePoly(g1,k)) ∩ recc(facePoly(g2,l))` intersected with the split half it
   belongs to. Comparing the piece's `R` to this EXPECTED cone (both symbolic) names exactly which
   ray was dropped.
3. Tag the producer that built it: `clipByFace` (arc-vs-half-plane or straight path),
   `clipPolyByConic` (arc-vs-arc), or `splitCell` (bounded / escape-to-infinity / `splitTwoArcPiece`
   / the `triHalf` two-arc-triangle fallback added this session).
Output: {shift -> culprit src, curveAfter, producer, dropped ray direction}. Expect a small number
of producers. On [0.5 0.5] the culprit was a `triHalf` piece (src[6 1]); the pinned shifts use OTHER
producers, so do this across several wrong shifts, not just [0.5 0.5].

## Phase 1 — State and PROVE the invariant (recession cone, not probing)

INVARIANT every producer must preserve: a piece's recession cone equals the recession cone of its
TRUE region (the intersection of the two operand faces' recession cones, cut by any split). In
particular a piece with no ray edge must have recession cone `{0}` (compact), and a piece with rays
must have recession cone equal to the span of exactly those rays.

Add `assertPieceReccConeCorrect(piece, reccExpected)`: compute the piece's recession cone from its
own edges (leading-term test above) and prove — with `isAlways`, not sampling — that it equals
`reccExpected`. Assert it on every piece right before `assemblePieces`. This turns a silent wrong
answer into a LOUD, localised, SOUND failure (it cannot miss a bad direction the way a far-point
probe can) and is the regression gate for Phase 2. Do NOT resurrect the reverted far-point
`assertNoOverExtendedFace` as anything but a throwaway smoke-test — it is unsound.

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

## Phase 3 — There is NO unrepresentable case; a "needs unbounded curved edge" signal is a BUG

CORRECTION (2026-08-04): do NOT add an "unrepresentable escape hatch" here — that is the general-case
trap. It is ESTABLISHED (`proveStageA–D`: 1140 dual boundaries classified, zero exceptions) that the
conjugate/biconjugate of a nonconvex QuaPol is representable as a QuaPar: every dual boundary is a
line, a line pair, or a PARABOLA, and QuaPar holds unbounded faces via RAYS and parabolic edges as
BOUNDED arcs. So the true `f*`/`f**` never contains an unbounded parabolic edge. Any place the code
appears to need one (an arc running to infinity as a face boundary; a piece that is neither compact
nor ray-bounded) is my SUBDIVISION computing the wrong geometry, not a limitation of QuaPar.

Consequence for the fix: the Phase-1 assertion and the `splitCell` "escape to infinity" guard are
BUG DETECTORS, not permanent escape hatches. Their firing means "the ray/arc bounding is wrong,
Phase 2 is incomplete" — the goal is that they NEVER fire on any conj/biconj input. Do not convert a
firing into a supported-input error; treat it as a failing assertion to be driven to zero.

The only genuinely-unsupported inputs remain the ones already out of scope and proven not to arise
here (hyperbolic edges, ellipses) — those keep their existing loud `notImplemented` errors, but they
are NOT what the far-field defect is about.

## Phase 4 — Acceptance is a PROOF, not a passing sample

The acceptance criterion is a SYMBOLIC theorem, checked exactly, that covers all points at once:
- **Partition:** the faces of `g` cover the plane with no gap and no overlap — proved with region
  algebra (pairwise face intersections have empty interior; union is the plane). The existing
  `maxQuaParResultsAreValidArrangements` is the structural half; complete it to a set-algebra proof.
- **Per-cell identity:** overlay `g`, `g1`, `g2`; on each common cell `C` (finite set), prove the
  polynomial identity `q_g ≡ max(q_{g1}, q_{g2})` — which, since the switching curve fixes the
  winner on `C`, is `q_g ≡ q_{g1}` (or `q_{g2}`), verified by `isAlways` / coefficient equality on
  `C`. This proves `g = max(g1,g2)` EVERYWHERE, far field included.
- Add such a `verifyMaxIsExactSymbolically(g1,g2,g)` helper and assert it on the three arc-vs-arc
  pins (plus a couple of the previously-wrong seeded shifts pinned by their exact shift value).
- Keep the recession-cone assertion (`assertPieceReccConeCorrect`) as a permanent, cheap gate so any
  future producer that drops a ray fails loudly and SOUNDLY at the source.
- The seeded numerical far-ring sweep stays only as a fast development tripwire; it is NOT the
  acceptance test and must never be cited as proof of correctness.

## Decision point (raise before Phase 2)

Incremental (patch ray-propagation in each producer, guarded by the Phase-1 assertion) vs. a rework
of the arc-vs-arc assembly. Recommend INCREMENTAL: the invariant + assertion make each producer's
fix independently verifiable, and Phase 0 should show only a few producers are at fault. Escalate to
a rework only if Phase 0 shows the ray-loss is spread across most producers.

## Do NOT repeat (tried this session, reverted — see git log around 0b5eee8)

- Deciding compactness / correctness by NUMERICAL SAMPLING (far-point probes, rings, the seeded
  sweep). Unsound — it can miss the one bad direction forever. Use the recession-cone / per-cell
  identity PROOFS above. Sampling is a development tripwire only.
- Raw-piece compactness guards (`arcFaceIsBounded`): got the arc/half-plane sign convention wrong AND
  are sampling-based. Replace with the symbolic recession-cone test, not a "fixed" probe.
- A hard post-assembly far-point check (`assertNoOverExtendedFace`): correct-ish but sampling-based
  and errors on ~every arc-vs-arc result because the defect is pervasive. Superseded by the symbolic
  recession-cone assertion.

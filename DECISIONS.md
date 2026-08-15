# Decisions & Dead Ends

_The record of what has been tried and ruled out, so it is not tried again._

**Append-only.** Never delete an entry. If a decision is later overturned, strike
the heading through (`## ~~...~~`) and add a new entry saying what changed.

**Nothing here is an action item.** Near-term actions live in `TODO.md`. Dead ends,
reverted attempts, rejected approaches and refuted diagnoses live here.

**Read this before attacking a hard problem** — especially one that looks like it
has been attacked before.

Newest entries at the top.

> **Seeded 2026-08-13** by copying the dead-end passages out of `TODO.md`. `TODO.md`
> was left untouched because another session was working in this project at the
> time, so **the originals are still there and this file duplicates them**. When
> CCA2 is next idle, delete the copied passages from `TODO.md` — specifically the
> `## Retired hypotheses — do not re-try` section, the `TRIED, reverted` block under
> `MAJOR FINDING 2026-08-04`, the `ATTEMPTED TWICE AND REVERTED` text in the step-4
> item, the `Tried and REVERTED` paragraph in the obsolete `decideWinner` note, the
> `--- OBSOLETE ---` block, and the `SUPERSEDED` section at the end — and leave
> `TODO.md` as action items only.

---

## 2026-08-15 — Two of the five "remaining bugs" were described WRONG. Measure before fixing.

Both descriptions had been written from a symptom and carried forward as fact. Each cost an
attempt before measurement refuted it. Corrected shapes below; neither is fixed.

### BUG 1 — "conjugateOfPiecePoly returns the conjugate of the chord"

- **What the record said:** the routine decides how many edges a piece has from
  `size(d.ineqs,2) == d.nv`, a COUNT standing in for "is this region unbounded", so the half-lens
  takes the unbounded convention and its arc is never read as an edge. Fix: derive the edge list
  from geometry.
- **What is actually true, measured** (instrumented dump of the pre-scatter constraint list):
  the half-lens arrives with `nv = 2` and FIVE constraints, and `edgeNo = [3 1 2 2 2]` —
  **three** constraints claim edge slot 2, all with both vertices on them:
      con 3: (s1+s2)^2 - 4*s1     con 4: (s1+s2)^2 - 4*s2     con 5: s2 - s1
  `getEdgeNosInf` numbers an edge by one of its ENDPOINT VERTICES, and a lens has two edges
  joining the SAME pair — so they are indistinguishable to it. The scatter
  `d.ineqs(edgeNo) = d.ineqs` is last-write-wins, so one edge is destroyed and another stored
  twice. Feed the arc first and both slots end up holding the CHORD; feed the chord first and both
  hold the arc. That is the whole of the "conjugate of the chord" symptom, and the count test is a
  consequence, not the cause.
- **Also measured, and this is why fixing the numbering is not enough:** hand-build the lens with
  ONLY its two genuine edges (`{(s1+s2)^2 <= 4*s2, s2 <= s1}` carrying `s1`, `nv = 2`,
  `nineq = 2`) and `conjugateOfPiecePoly` still returns **1 cell**, not the 4 the piece needs
  (2 vertex cones + 2 edge cells). So the downstream cell generation does not handle a 2-vertex
  CURVED region either.
- **Tried, twice, both REVERTED:** (i) spreading colliding genuine edges over distinct numbers;
  (ii) the same, plus dropping constraints that bound no edge — a constraint whose curve between
  the shared vertices leaves the region (which is what con 3 is, and which LP redundancy cannot
  see because it is not linear), and constraints active at a single vertex of a bounded region.
  Both make the lens's slots correct — measured, `nineq = 2` with the arc and the chord in
  distinct slots — and both take the second conjugation from WRONG VALUES to
  `QuaParCPLQ:conj:emptyResult`, no pieces at all.
- **Before retrying:** the edge numbering is necessary and not sufficient. Do the downstream half
  FIRST — make `conjugateOfPiecePoly` produce 4 cells for a bounded 2-vertex region with one
  curved edge, checked on the hand-built lens above, which needs no pipeline. Only then re-apply
  the numbering fix. Note also that `size(d.ineqs,2) == d.nv` is used a SECOND time, to choose
  between the polyhedral and the quadratic-aware normal-cone routines
  (`getNormalConeVertex`/`getNormalConeEdgeQ3` vs `getNormalConeVertexQ`/`getNormalConeEdgeQ`),
  so changing the count changes which of those runs.
- **Separately: the pinned test no longer fails the way it says it does.** Since `conj` began
  returning a MESHED QuaPar for a bounded multi-face domain (2026-08-13),
  `biconjugateOverATwoFaceSubdivisionIsTheEnvelope` ERRORS at `quaPolToPlq:curvedEdge` — the
  second conjugation is handed a curved QuaPar and `quaPolToPlq` requires a polyhedral domain. The
  symbolic route reaches the lens defect only when forced. Both need fixing; they are different.
- **Refuted on the way:** `convEnvCPLQ` is NOT a route to `f**` here. `f** = conv f` for a compact
  domain, but `convEnvCPLQ` is Step 1, a PER-TRIANGLE envelope with no cross-piece hull — measured
  on this input, it returns the per-piece envelopes (0.25 at (0.5,0.5) where the truth is 0).

### BUG 5 — "a piece that spans TWO sub-arcs of the SAME conic"

- **What the record said** (written 2026-08-14, from the piece's vertex list): piece 4 `src [1 6]`
  has g1's arc cut twice and spans both sub-arcs, so its one curve slot holds one and the other
  becomes a chord.
- **Measured, and it is false:** evaluating piece 4's own conic at its vertices gives
      V1 (-2.960609, 0.9606088) -> 0        V5 (-2.744821, 1.372827) -> 0
      V4 (-2, 2) -> +0.149        V2 -> +0.161        V3 -> +0.327
  with a tolerance of 8.8e-07. Only V1 and V5 are on it. The two curved edges are on **different**
  conics: piece 4's own slot holds the SPLITTING curve `{g1f1 = g2f6}`, and the flattened edge
  `(-2,2) -> (-2.744821,1.372827)` is a sub-arc of **g1's arc**, which its neighbour piece 5
  `src [2 6]` carries curved on `[-0.015625 -0.03125 -0.015625 -0.25 0.25 -1]`.
- **So it is the ordinary two-DIFFERENT-arcs case**, which `splitTwoArcPiece` exists for — and the
  cell does carry the arc going in: `src [1 6]`, `nv = 4`, `arcPos0 = 2`, arc
  `(-2,2) -> (-3.125,1.125)`. But `splitTwoArcPiece` is called exactly ONCE on this input and not
  for this cell, so the loss happens before it: either the crossed-arc restoration finds no
  position for the sub-arc (`findArcPosition` returning 0) or the piece is emitted before reaching
  it.
- **AND THE CAUSE IS NOW LOCATED, by that instrumentation.** For `src [1 6]` on that shift:
      HITS: 2 hits on edges [2 3] at (-2.744821,1.372827) and (-2.960609,0.9606088), arcPos0 = 2
      RESTORE half: nv=5 curveAfter=5 **p=4** straightCut=0
  So the crossed-arc restoration works perfectly: it FINDS the inherited sub-arc at edge 4, the
  splitting curve is genuinely curved (straightCut=0), and the half is correctly handed to
  `splitTwoArcPiece(half, 4, arcEc0)`. **`splitTwoArcPiece` then returns it UNSPLIT**, which its
  own header says it may do — "if neither works the piece is returned unsplit ... the assembly's
  own arrangement check is what would catch a dropped arc". It did, three stages later, as the
  orphan.
- **Why it finds no cut: THE TWO ARCS ARE ADJACENT.** They share vertex V5 (arc at edge 4, split
  curve at edge 5, `nv = 5`). Its two candidate chords are `arcPos+1 -> c+1` = V5->V1 and
  `arcPos -> c` = V4->V5 — but those ARE edges 5 and 4 themselves, so one chain comes out with 2
  vertices and the `numel(chain) < 3` guard skips both. The `nv == 3` fallback that handles the
  shared-vertex case (cut from the shared vertex to the midpoint of the opposite straight edge)
  does not apply at `nv = 5`.
- **The fix, not yet written:** generalise that shared-vertex case to `nv >= 4` — cut from the
  shared vertex S to a NON-ADJACENT vertex. On this piece S = V5 and the diagonal V5->V2 gives
  chains {2,3,4,5} (carrying edge 4, the inherited arc) and {5,1,2} (carrying edge 5, the split
  curve): one arc each, which is the invariant. Guard it with `insideStraightHull` exactly as the
  existing candidates are, and fall through unsplit when no diagonal is interior.
- Do NOT write a same-conic sub-arc splitter; that shape has never been observed, and the code for
  it was written and removed on 2026-08-15.

## 2026-08-14 — RESOLVED: the two-arc ray split, and why two attempts at it failed

- **Outcome:** `arcVsArcRefusesAnUnboundedTwoArcSplit` is GREEN and `maxQuaParTest` is 28 / 0 (from
  25 / 1). The entries below about the ray split stand as history; the construction was right in
  outline all along, and **neither reverted attempt failed for the reason it was thought to**.
- **What actually blocked it, both times.** Two defects, neither in the split:
    1. `pieceRecessionRays` took the parabola's axis from `arcNullDirs`, which solves `d·Q·d' = 0`
       *exactly* and returns **nothing** when `b²−4ac` comes out negative — which is what a
       floating-point parabola's `Q` does about half the time, being only semidefinite up to
       rounding (measured `−2.78e-17` on the pinned fixture's first arc). The derived chord was
       then silently never emitted, so the half's constraint region stayed a slab open at BOTH
       ends and `reccConeViolation` refused it. This is why **check (5) of the six did not
       separate the cases**: it was answering a question about a constraint set that was not the
       piece's. `parabolaArcFrame` has always taken the smallest-magnitude eigenvector; do that.
    2. `curveAfter ≠ 0` was read as "this edge is curved" in five places, when `boundedPiece` also
       sets it for a STRAIGHT splitting curve (`curveEc` all zeros). Two of those places call
       `parabolaArcFrame` on the zero conic and raise `degenerateAxis` — **so `MAXQP_ASSERT`
       crashed on three of the four arc-vs-arc fixtures, and the invariants that would have named
       all of this never ran on the inputs that needed them.** An invariant that errors is off.
- **And the seeded shift `[1.4979 3.6486]`, blamed on the split twice:** the split's halves pass
  all three exact invariants. The wrong 0.3531 came from a *different* piece, `src [2 4]`, whose
  only finite edge `polyConstraints` had dropped (defect 2) — the refusal had merely been masking
  it. It is now refused explicitly; see the next entry.
- **Lesson worth keeping:** when a gate refuses a construction you have independent reason to
  believe is right, suspect the gate. Both reverts were correct decisions on the evidence
  available, and the evidence was wrong because the tooling was broken in a way that was silent.

## 2026-08-14 — "`{f1=f2}` must be a PAIR of parallel lines here" (diagnosis, refuted by measuring)

- **Tried:** explaining the last wrong answer on seeded shift `[1.4979 3.6486]` — piece `src [2 4]`
  carrying g1 face 2's ZERO quadratic while g2 face 4 beat it by `+Inf` along its own ray — as a
  subdivision gap: `{f1=f2}` is a degenerate conic, so it can be a pair of parallel lines, and a
  half strictly on one side of the line `splitCell` cut along could be crossed by the other, which
  would leave through the recession cone and contribute no finite crossing to find. That story fits
  every symptom, and a guard was written for it (`assertWinnerHoldsAtInfinity`).
- **Why it failed:** the conic says otherwise. For that cell
  `diffRow = [0 0 0 0 0 0 0 −1.4979 −3.6486 5.4652]` — **its entire quadratic part is zero**, so
  `{f1=f2}` is a SINGLE straight line and nothing straddles it. `delta = 0`, `det3 = 0`,
  `eig(Q) = [0 0]`.
- **What it actually was:** the piece's only two vertices ARE the two crossing points, so both lie
  ON that line; `assignSide` had nothing to read the winner from and fell back to a centroid which
  is on the line too. The winner of a whole unbounded piece came out of floating-point noise.
  `assignSide` now reads such a piece in its RECESSION CONE, sharing the probe
  `assignSideFromCone` has used since it was written for the same problem in
  `splitUnboundedAtOneCrossing`. That shift now assembles CORRECTLY: the seeded sweep is 17 exact /
  0 wrong / 1 errored of 18, from 16 / 0 / 2.
- **Before retrying:** the guard is kept as a cheap exact backstop and currently fires on nothing.
  Do not read its existence as evidence that the pair-of-lines case occurs — no input has ever been
  observed to produce it. **Classify the conic before theorising about its shape**; one line of
  `eig(Q)` would have saved the detour.

## 2026-08-14 — Why check (5) may have failed to separate the two-arc ray split's cases

- **Tried (2026-08-13):** "each half's recession cone must equal the cone its own rays span" was
  one of the six checks applied to the unbounded two-arc ray split, and it did not separate the
  pinned fixture from the seeded shift that assembles to a wrong value. It was therefore recorded
  below as exhausted, and the search was pointed at WHERE the cut starts instead.
- **Why that conclusion may be premature:** the check is implemented by `pieceRecessionRays`, and
  until 2026-08-14 that routine derived the arc's chord — the constraint that makes a concave-side
  arc-piece compact at all — by asking which side the piece's other VERTICES fall on. That is the
  same rule `DECISIONS.md` already records as unsound one level up, in `QuaPar.chordCuts`, where it
  killed two green tests: a lens has both vertices ON the chord, so they say nothing. It also had
  no gate on *when* a chord may be emitted, so a piece that genuinely straddles the chord line
  could be handed one and be reported compact when it is not. Either way the check was answering a
  question about a constraint set that was not the piece's.
- **Corrected:** both questions are now settled by the conic itself. Along the chord `X0 + t·ch`
  the conic restricts to `q(t) = A·t·(t−1)` with `A = ch·Q·ch'`, because both endpoints are on it.
  `A ≤ 0` means the chord's interior is inside the kept side, the piece straddles the line and no
  chord may be emitted (the lens, and every convex-side face). `A > 0` means the piece touches the
  line only at the arc's endpoints and lies on the side the arc's own interior points are on,
  reached along the parabola's axis from the chord midpoint. The vertex test survives as a veto.
- **Status: CONFIRMED, and it was the whole of it.** Correcting the axis derivation is exactly what
  turned both halves of the pinned fixture's cut from "cannot be proved" into "proved". See the
  RESOLVED entry at the top.
- **Evidence:** `arcVsArcRefusesAnUnboundedTwoArcSplit` green; `maxQuaParTest` 28 / 0; fast bucket
  202 / 0.

## 2026-08-14 — A newly minted OUTGOING ray was given sign +1 (a live bug, not a dead end)

- **Recorded here because the reasoning is easy to re-derive wrongly.** `polyConstraints` reads a
  ray's outward normal as `sign · rot90ccw(direction)`, and both `dirIn` and `dirOut` store the
  direction pointing from the apex OUT to infinity. A piece is walked CCW with its interior on the
  LEFT, so the incoming ray is traversed along `−dirIn` and takes `+1`, and the outgoing ray is
  traversed along `+dirOut` and takes `−1`.
- **The bug:** both branches of `splitCell` that mint an escaping ray wrote `+1` for the OUTGOING
  one. That flips the kept half-plane to the far side of the ray's line, so the piece's constraint
  region is the reflection of its true region across that line — over-extended on one side, short
  on the other, which is the shape every far-field wrong answer here has had.
- **What NOT to do:** do not "fix" the INHERITED signs the same way. A sign is a property of the
  `P{k}` entry the ray came from, not of its role, and a face whose whole boundary is two rays
  sharing one apex can legitimately carry the same sign on both — that generalisation is itself a
  recorded fix, with its own test.
- **Evidence:** derived from `polyConstraints`' own HISTORY note; spotted independently while the
  ray split was being reverted and recorded in that commit. `newRaySign` now states it once.
  **Verified 2026-08-14**: no regression on any bucket, and an independent review confirmed the
  derivation (the previous `+1` gave the two halves of a split the *same* outward normal across
  their shared ray — overlap on one side, hole on the other).

## 2026-08-13 — Making the arc's chord a REAL EDGE by splitting the neighbour (option B)

- **Tried:** Nothing was built. It was proposed — and recommended — as the way to close a face
  whose arc is concave towards it, resolving the open choice left by the entry below.
- **Why it failed:** Geometry, checked before implementing. The chord runs through the **interior
  of the neighbour** on the other side of the arc, so making it a real edge splits that neighbour
  and leaves the offending face's own edge list unchanged. It cannot supply the missing constraint
  to the face that needs it. No arrangement can: the chord is not on that face's boundary.
- **Before retrying:** Do not. The chord must be **derived per face**, which is what
  `QuaPar.chordCuts` now does — and it resolved the whole far-field defect. Two soundness rules
  came out of building it, both learned by breaking tests: the side must be read from a point just
  INSIDE the face (the arc's midpoint stepped along the inward normal), never from the face's other
  vertices — a lens has both vertices ON the chord, so they say nothing and the wrong side is
  chosen; and a chord is emitted only when every other vertex and ray direction agrees, i.e. when
  it is redundant for the face, so it can never shrink a face below its true region.
- **Evidence:** `QuaPar.chordCuts` + `SUPPORT_MATRIX.md` §4.1. The vertex-based side rule broke
  `maxQuaParCurvedMatchesGroundTruthOnRandomQuadrilaterals` and
  `maxQuaParSplitsACellThatAlreadyCarriesAnArc`; the interior-point rule turns both green.

## 2026-08-13 — Splitting an UNBOUNDED two-arc piece with a ray (implemented, reverted)

- **Tried:** An unbounded piece carrying two arcs cannot be separated by a chord (no chord closes
  an unbounded piece), so: cut with a RAY from the vertex between the two arcs, along a direction
  the piece recedes in. Each half then keeps one arc, one original ray and the new one. It works —
  `arcVsArcRefusesAnUnboundedTwoArcSplit` assembles with it.
- **Why it failed:** It makes seeded shift `[1.4979 3.6486]` assemble to a value wrong by 0.3531.
  A silent wrong answer is worse than the refusal it replaces, which is what that test's own
  comment says, so it was reverted.
- **Before retrying, fix:** Six checks were tried and NONE separates the good case from the bad:
  (1) the ray recedes every straight constraint; (2) it stays inside BOTH arcs for its whole
  infinite length — necessary, a ray did leave through an arc; (3) each half admits a point just
  inside its own arc — necessary, it caught a mis-oriented pair; (4) the new ray's SIGN (an
  OUTGOING ray needs `-1` in `polyConstraints`' `sign*rot90ccw(dir)`, not the `+1` the
  escape-to-infinity branch uses — that branch looks wrong the same way and deserves its own
  check); (5) each half's recession cone equals the cone its rays span; (6) scoping to the
  half-strip shape. Start past all six: the open question is WHERE the cut starts, not which
  direction it takes — test whether the vertex between the arcs is the right starting point at
  all, or whether the cut must start on one of the arcs.
- **Evidence:** `maxQuaParTest/arcVsArcMatchesGroundTruthOverRandomShifts` (seed 20260803, N=18)
  versus `arcVsArcRefusesAnUnboundedTwoArcSplit`; commit "Revert the unbounded two-arc ray split".

## 2026-08-13 — "Equal areas and a shared polyline prove two halves tile" (refuted reasoning)

- **Tried:** Concluding, from a measurement, that a two-arc split was NOT responsible for a
  coverage hole: the two halves shared their cut polyline bit-exactly (17 digits), were both CCW,
  and their polygon areas summed to the parent's exactly. That sent the next session's search into
  assembly instead.
- **Why it failed:** Both facts were true and the conclusion was wrong. Area says nothing about
  which side of a **bent** boundary a point falls on. The halves did not tile: the cut polyline
  left one of them REFLEX at the bend, and every point-location test here reads a face as an
  intersection of half-planes, which is exact only for a CONVEX face.
- **Before retrying:** To decide whether pieces tile, use the per-edge cross product at the
  disputed point, not area. `splitAtReflexVertex` now splits such a half along a diagonal.
- **Evidence:** `arcVsArcDoesNotCrashOnSeededQuadSplits`, fixture 1, the point
  `(0.998629534754574, -0.0523359562429444)` — outside both halves, inside the parent.

## 2026-08-13 — Subdividing a bounded arc-piece whose constraint region is non-compact

- **Tried:** Twice, on the two quadrilateral-fixture pieces `src [1 2]` and `[1 6]`.
  Cutting the piece with a line parallel to the arc's chord, just past the arc's
  deepest point, to close the non-compact side.
- **Why it failed:** Both halves came back still non-compact, and the straight half
  wrongly kept the arc tag from `clipPolyHalfPlane`. The second attempt settled the
  reason: for a piece whose arc is **concave towards it**, the constraint set is a
  wedge intersected with the **outside** of a parabola. A cut parallel to the chord
  leaves the arc-side sliver still receding along the two side edges' own direction,
  which neither the cut nor the parabola blocks.
- **Before retrying, fix:** Do not retry any subdivision. This is a **representation**
  decision, not a subdivision bug. The constraint that would close the region is the
  arc's own **chord** — redundant for the region (which lies entirely on one side of
  it) but not one of its edges. So either a face may list a redundant bounding conic,
  or the chord becomes a real edge by splitting the neighbour across it. Pick one of
  those two before writing any more clipping code.
- **Evidence:** `maxQuaParTest`, quadrilateral fixture, pieces `src [1 2]` and `[1 6]`.
  (A separate `fixArcTag` fix has since stopped the losing half carrying the arc's
  constraint, but that does not make the cut work.)

## 2026-08-04 — Guarding against non-compact arc-pieces (four guards, all reverted)

- **Tried:** Four different guards, after establishing that arc-vs-arc results are
  correct near the arcs and wrong in the far field:
    1. Piece-level compactness guard on `triHalf` output.
    2. Piece-level compactness guard on every bounded curved piece.
    3. Post-assembly check: "a bounded face admits a far point", using QuaPar's own
       EC + P signs.
    4. Post-assembly check: "two faces disagree at a far point".
- **Why it failed:** (1) and (2) used the wrong sign convention and rejected valid
  faces — (2) errored on all three fixtures. (3) and (4) **detect correctly**, but
  they also fire on `[-1 0.75]` and `[2 -0.5]`, because those results genuinely *are*
  non-compact far out. A backstop that refuses non-compact results is therefore
  correct and would turn every arc-vs-arc result red.
- **Before retrying, fix:** Do not add a guard. The fix must be **upstream**: give an
  unbounded quadratic face's arc sub-pieces their **ray** boundaries, so they stay
  unbounded and compact-as-faces, instead of closing them with straight edges. That
  is a rework of the arc-vs-arc clip/split, not a guard.
- **Evidence:** Ring sweeps of `g.eval` against the pointwise max — `[-1 0.75]` 0/60
  wrong at radius 8 but 2/60 at radius 30 (worst 6.1); `[2 -0.5]` 2/60 at radius 8,
  11/60 at radius 30 (worst 33.8). `[0.5 0.5]`: assembled FACE 15 carries g1-face-1's
  quadratic over a tiny triangle near (-2,2) yet admits (-3.98,0.61) two units away,
  overlapping the correct zero face (FACE 9).

## 2026-08-03 — `decideWinner`: sampling `diff` at 1-D stationary points

- **Tried:** Making `decideWinner` sound on unbounded cells by additionally sampling
  `diff = f1 - f2` at the 1-D stationary point of each edge and each ray.
- **Why it failed:** Sound and strict, but it does not fix the defect and was
  confirmed neutral on all three fixtures. `{diff = 0}` is a parabola (rank-1 `Q`,
  which `splitCell` already asserts) whose branches run to infinity along `Q`'s null
  direction, strictly *between* the cell's two bounding rays. The sign change is off
  the 1-skeleton entirely, so no amount of sampling along edges and rays sees it.
- **Before retrying, fix:** Do not attempt a probe-based patch. A real fix needs
  **both** a sound sign-over-the-cone test in `decideWinner` **and** a `splitCell`
  that can cut a cell along a parabola entering and leaving at infinity — the current
  one asserts exactly two finite crossings, and here there are zero.
- **Evidence:** Cell `src[4 3]` stores `winner=g2` but `g1=0 > g2=-8.648` at the
  interior point `s=(-5.4843,1.5866)`; symmetrically `src[6 5]` stores `winner=g1`
  while g2 wins at `s=(-1.1298,4.1007)`. 4/68 samples, worst 8.648. Kept out of the
  commit.

## 2026-08-03 — Returning both components from a disconnecting curved cut

- **Tried:** Replacing `clipPolyByConic`'s blanket refusal of a disconnecting curved
  cut with "return both components".
- **Why it failed:** Two components are **unrepresentable**. A separated survivor
  would need the cutting conic running to infinity, i.e. an unbounded curved edge,
  which `QuaPar` cannot hold. Separately, the premise was wrong: the cut measured
  **connected** on all three fixtures, so the single-cell construction was right.
- **Before retrying, fix:** `Do not retry` unless `QuaPar` gains unbounded curved
  edges. The blanket refusal was instead replaced with a real connectivity test,
  which now refuses only the genuinely disconnecting case — and that case does not
  arise on these fixtures.
- **Evidence:** `partitionReport` OK on all three fixtures after the connectivity
  test landed.

## 2026-08-13 — Refuted diagnoses (the analysis was wrong, not just the fix)

Each of these was written up in detail and then disproved. Re-deriving any of them
from the same symptoms is the failure mode this entry exists to prevent.

- **"Piece 5 (`src[2 2]`) is over-extended and its ray should terminate."** Refuted
  by the sign data: `evalConic(Ecut)@V = [0, -0.046, -0.015]`, so the vertex `(-2,2)`
  sits on the g2-face-1 (discard) side and the far ray is legitimately g2-face-2. The
  clip is correct. The geometric measurements in that write-up are still good; the
  conclusion is not.
- **"The missing g1f6 ∩ g2f2 mirror is the defect."** The sharper diagnosis that
  replaced the above, itself now superseded — the mirror (`src[6 2]`) is built
  correctly.
- **"`decideWinner` wrongly declares an unbounded cell decided" (for `[2 -0.5]`).**
  The real cause was `assignSide` reading the winner at `piece.V(2,:)`, which for
  `splitCell`'s unbounded "rest" piece (`curveAfter=1`) is a **crossing point** where
  `diff ≈ 0` and its sign is noise, so both halves could get the same winner. Fixed
  in `53fc9fd` by reading the vertex farthest from `{diff=0}`. The parabola-to-infinity
  write-up was the wrong lead for this fixture.
- **"Piece 16's ray is piece 5's partner."** It is not — piece 16's ray lies on
  `x+y=0.5`, a parallel but different line from piece 5's `x+y=0`. They were never
  meant to match.

## 2026-08-03 — Retired hypotheses about the orphan ray — do not re-try

- **Tried / Why each failed:**
    - *The orphan ray is one physical ray covered to different extents* — nothing
      lies on it.
    - *The cut must be restricted to the arc's own span* — it must not; the straight
      half-planes are applied first.
    - *`clipPolyByConic` emits clockwise cells* — a guard on every bounded output
      never fires.
    - *`polyL` is non-convex* — all four curved faces measured convex.
    - *Pair `(6,1)` is spurious* — it is a genuine thin cell. The "zero intersection"
      reading was a resolution artifact: a 0.0625-spaced grid cannot see a cell of
      area 0.008.
- **Before retrying, fix:** `Do not retry.` Each was measured, not argued.

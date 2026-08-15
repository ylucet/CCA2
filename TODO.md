# TODO

## 2026-08-13 -- the far-field defect, worked steps 1-6 of the plan

**Where it stands (updated 2026-08-13, end of the second pass).** `maxQuaParTest` **25 pass / 1 fail**; fast bucket 200 pass / 1 fail. The one red is `arcVsArcRefusesAnUnboundedTwoArcSplit` (first item below). `sweepMaxQuaParCurvedSplit(20260802,200)` went 30 -> 59 assembled of 142 sampled, with 0 of 1031 result vertices, 571 midpoints and 3540 interior points wrong. What the whole defect turned out to be is one sentence, and it is recorded in `SUPPORT_MATRIX.md` 4.1: **a curved edge is a bounded ARC and its conic is not**, so the point-location rule admits far-away points on a parabola's concave side; `QuaPar.chordCuts` derives the missing constraint. Original note follows.

**Where it stood mid-session.** `maxQuaParTest` 21 pass / 5 fail. The seeded far-field sweep
(`arcVsArcMatchesGroundTruthOverRandomShifts`) and the unit-square pin are GREEN, and on a
397-quadrilateral seeded sweep the arc-vs-arc results wrong in the far field went **7 of 64 to 0
of 64** (200 directions at radii 1, 5, 50, 500).

**Three exact piece invariants** (`assertPiecesWellFormed`, global `MAXQP_ASSERT`, off by default;
1 = containment + winner domination, 2 = also the symbolic recession cone): containment in the
source cell (vertices, ray directions, straight edges against the source arc, and the piece's own
arc against every source constraint), the carried operand actually dominating on the piece, and the
encoded region receding exactly where the piece declares rays. Nothing sampled.

**Two producers fixed.**
1. `splitCell` read "one finite crossing" as a tangency and returned the cell intact. On an
   UNBOUNDED cell the curve can enter there and leave through the recession cone, cutting it in
   two; `splitUnboundedAtOneCrossing` decides the two cases from whether the branch direction
   recedes the cell, and each half keeps one original ray plus the escaping branch.
2. `clipPolyByConic` skipped a cut whose only crossings are AT the cell's corners -- the generic
   arc-vs-arc case, since conjugates of triangles sharing a primal edge have arcs through the same
   two dual points. Now decided from one representative point per boundary element; a two-arc
   survivor is cut along `A -> M -> B` with M between the arcs and each half given the midpoint as
   a third vertex (two-vertex lunes sharing a chord defeat assembly: the chord never becomes an
   edge).

**Two new tools.** `QuaPar.eval` validate mode (`QUAPAR_VALIDATE`) errors when two faces admit a
point with DIFFERENT values; `verifyMaxIsExactSymbolically` proves `g = max(g1,g2)` over whole
regions by closed-form minimisation over each face-pair intersection. The three arc-vs-arc fixtures
verify with zero findings.

### Fixed since (2026-08-13, second pass)

- The **unit square now VERIFIES**: zero findings from `verifyMaxIsExactSymbolically`, zero
  ambiguous points from eval's validate mode over 900 ring points. Two causes, both "an edge
  identified by its endpoints alone", which an arc-vs-arc arrangement makes ambiguous because FOUR
  arcs run between the same two dual points: `matchHalfEdges` now also requires two curved
  half-edges to lie on the same conic, and `QuaPar.orderEdges` now looks the next boundary edge up
  among the FACE's own edges instead of the whole edge list.
- `fixArcTag`: a clipped half that does not keep the arc no longer carries the arc's CONSTRAINT.
- `splitTwoArcPiece` re-locates both curve indices from geometry, so a stale index from a rebuilt
  vertex list can no longer index off the end (the two seeded crash fixtures no longer crash).
- `dropDegeneratePieces`: collapsed pieces (2 vertices, no arc, bounded) no longer reach assembly.

### 2026-08-15 -- FOUR of the five bugs fixed; ONE red left in the whole repository

`maxQuaParTest` 29 / 0, `conjCPLQTest` 25 / 0, `unboundedFaceTest` 18 / 0, fast bucket 204 / 0,
slow bucket **114 / 1**. Bugs 2, 3, 4 and 5 are fixed; bug 1 is at 5 of 7 probe points, was 0 of 7,
and is the only remaining red. `maxQuaPar` has NO open case -- the seeded arc-vs-arc sweep is
**18 exact / 0 wrong / 0 errored of 18**, from 16 / 0 / 2 on 2026-08-13. What it took is in `SUPPORT_MATRIX.md` 4.1 and `DECISIONS.md`; the short version is that
**neither of the two earlier attempts at the last red failed for the reason it was thought to** --
the tooling that judged them was itself broken, in two ways, and silently.

### Still open, in the order they should be taken

- [x] **`MAXQP_ASSERT` is now ON in `maxQuaParTest`** (2026-08-15), at level 1, via a
      `TestMethodSetup` that restores the previous value on teardown. 28 / 0 with it on. Level 2
      stays opt-in -- it costs seconds per call, and the tools that want it call
      `pieceRecessionRays` directly.

### Bugs, in the order they should be taken (2026-08-15)

- [ ] **BUG 1 -- TWO defects, and the description this item used to carry was wrong about the
      cause.** Measured 2026-08-15; full write-up in `DECISIONS.md`.
      **(1a) The pinned test no longer fails the way it says.** Since `conj` began returning a
      MESHED QuaPar for a bounded multi-face domain, `biconj`'s second conjugation is handed a
      CURVED QuaPar and `quaPolToPlq` refuses it -- the test ERRORS at `quaPolToPlq:curvedEdge`
      before any of (1b) is reached. Either the second conjugation learns to take a curved QuaPar,
      or `biconj` routes a bounded multi-face domain through the symbolic form on purpose.
      **(1b) The lens defect, reachable by forcing the symbolic route.** The cause is NOT the
      `size(d.ineqs,2) == d.nv` count: it is that `getEdgeNosInf` numbers an edge by one of its
      ENDPOINT VERTICES, and a lens has two edges joining the SAME pair, so they get the same
      number and the last-write-wins scatter destroys one. Fixing the numbering is NECESSARY and
      NOT SUFFICIENT -- measured: hand-build the lens with only its two genuine edges and
      `conjugateOfPiecePoly` still returns 1 cell where the piece needs 4. **Do the downstream half
      first** (4 cells for a bounded 2-vertex region with one curved edge, checked on that
      hand-built region, which needs no pipeline); the numbering fix is written up in
      `DECISIONS.md` and can be re-applied after.
      **PROGRESS 2026-08-15 (uncommitted at the time of writing): 5 of 7 probe points now right,
      was 0 of 7.** Three defects fixed, in the order they had to be:
        1. `spreadCollidingEdges` -- give a lens's two edges distinct numbers instead of letting
           the scatter destroy one. Scoped to fire only on that signature.
        2. `getNormalConeVertexQ` indexed its second constraint as `j+1` UNWRAPPED, so it raised
           `badsubscript` on any BOUNDED region -- which is why the only caller sent every bounded
           region to the POLYHEDRAL `getNormalConeVertex`, whose cones come from the CHORD and are
           wrong for a curved edge. Wrapped cyclically (identical to `j+1` for the unbounded
           layout, so nothing that worked changes), and the dispatch now asks whether a constraint
           is quadratic rather than comparing a constraint COUNT to the vertex count.
        3. `biconj` on a bounded multi-face domain takes its FIRST conjugate in symbolic form
           (`conjCPLQ(..., 'symbolic')`), because the second conjugation cannot take the curved
           MESH `conj` now returns and died at `quaPolToPlq:curvedEdge`.
      Unit level: the half-lens conjugates to 3 cells matching a brute-force sup at all 10 probe
      points (2 identical wrong cells before).
      **REMAINING, and it is a REFACTOR, not a fix.** `f**` is exact at 5 of 7 probe points; the
      other two come back `+Inf` -- a hole in the DOMAIN, not a wrong value, since `f**`'s domain
      is the INTERSECTION of the per-piece conjugate domains and one piece of `f*` still
      conjugates onto too small a set.
      What that piece needs is its two lens edges in slots 1 and 2, and those slots are held by
      constraints that bound no edge. **Freeing them by DROPPING those constraints was tried and
      is unsound** -- a constraint active at exactly one vertex of a convex region can still be
      essential, and removing it enlarges the piece, which SHRINKS its conjugate domain. Measured:
      with the drop, `f**` is exact at `(0.25,0.25)` and `(0.1,0.1)` and `+inf` at `(0.9,0.6)` and
      `(0.6,0.6)`; without it, exactly the other way round. Both 5 of 7, one of them sound.
      The real fix is to give `conjugateOfPiecePoly` an explicit EDGE LIST instead of a count with
      two conventions (`endNv = nv` or `nv-1`; edge `j` at `ineqs(j)` or `ineqs(j+1)`). It cannot
      be done in that routine alone: the loop variable `j` indexes `getNormalConeEdgeQ`/`Q3`'s
      output and `getSubdiffVertexT2`/`T2Q`'s `subdE` at the same time, so all four move together.
      Start from the hand-built lens in the probe -- it needs no pipeline and runs in seconds.

- [x] **BUG 2 -- FIXED 2026-08-15.** `region.removeTangent` built a TANGENT LINE to a quadratic
      constraint at a vertex where that quadratic's GRADIENT VANISHES -- the apex of a cone, which
      is exactly where the Step 3 split conics of an unbounded fan meet. There is no tangent line
      there, every direction is tangent, and whatever it computed was meaningless; it then deleted
      a constraint matching that "tangent". On the 4-cone fan it deleted `-s1 <= 0` from the cell
      carrying `s1^2/4 + s2^2/2`, leaving only `{s2 <= 0, s2^2/2 - s1^2 <= 0, s1^2 - 2*s2^2 <= 0}`
      -- two constraints BLIND TO THE SIGN of `s1` -- so the region became symmetric under
      `s1 -> -s1` and claimed the mirror wedge. `f*(-3,-2.4)` is now **4.5**, the truth; it was
      5.130. `conjCPLQTest` 25 / 0.
      This is the SIBLING of the 2026-08-02 fix: `simplifyUnboundedRegion` fell into the same trap
      on the same input, and `region.witnessAwayFrom` was written for it. One input, one trap, two
      routines -- **a vanishing gradient at a cone's apex is a recurring failure mode here.**
      Cleared on the way, so nobody re-checks them: `redundantSubset` (certifies nothing there,
      correctly) and `simplifyUnboundedRegion` (leaves the constraint alone).
      `step3DropsCellsOnSomeUnboundedAssemblies` asserted the GATE firing; it no longer does, so
      that test is renamed `step3UnboundedAssemblyAgreesWithItsOwnPieces` and now pins what the
      gate protects.

- [x] **BUGS 3 and 4 -- FIXED 2026-08-15.** `convEnvUnbounded` computed only the AFFINE envelope
      and raised `convexAlongRay` as soon as `d'Qd > 0` along a ray. Two shapes are now derived and
      implemented, each with its proof in the source:
        * **WEDGE, one flat ray and one convex ray.** `co q` is `q` with its CROSS TERM deleted:
          `q(v) + alpha*g1 + beta*g2 + beta^2*A22/2`. A negative `A12` means `d1 + t*d2` recedes
          with negative curvature, so the envelope is `-inf` -- now reported rather than answered.
        * **HALF-STRIP convex along the ray, base edge Q-ORTHOGONAL to it.** `q` separates, so
          `co q = q(v1) + s*(q(v2)-q(v1)) + t*<grad q(v1),d> + t^2*(d'Q d)/2` -- the affine
          interpolant along the concave base plus the convex part along the ray.
      `w'Q d ~= 0` is deliberately not handled and is refused loudly.
      **`unboundedFaceTest` 18 / 0**, from 16 / 2. Fast bucket 204 / 0.

- [x] **BUG 5 -- FIXED 2026-08-15.** `splitTwoArcPiece` found no cut when the two arcs are
      ADJACENT: its two candidate chords join the arcs' facing endpoints, which for adjacent arcs
      ARE the arcs' own edges, so both chains come out with two vertices, the `numel(chain) < 3`
      guard skips them, and the piece was returned unsplit with one arc flattened to its chord.
      The `nv == 3` shared-vertex fallback did not apply at `nv = 5`.
      Generalised to `nv >= 4` with the ordinary DIAGONAL from the shared vertex to a non-adjacent
      one -- the two arcs leave that vertex in opposite directions around the boundary, so any such
      diagonal puts one arc in each half by construction. Same `insideStraightHull` guard as the
      existing candidates, and each half goes through `splitAtReflexVertex`.
      **Measured: the seeded sweep goes 17 exact / 0 wrong / 1 errored -> 18 / 0 / 0.**
      `maxQuaParTest` 29 / 0, fast bucket 203 / 0. Pinned by
      `arcVsArcSplitsTwoADJACENTArcsOnAPieceWithADiagonal` -- by VALUE, and with its own test
      because `arcVsArcMatchesGroundTruthOverRandomShifts` asserts `nWrong == 0` and would have
      counted this input in `nErr`.
      (The "two sub-arcs of the same conic" description this item used to carry was refuted by
      measurement; `DECISIONS.md` records it.)


- [x] **FIXED** `arcVsArcDoesNotCrashOnSeededQuadSplits` -- the last piece was a REFLEX vertex left by the
      bent two-arc cut: half-plane point location cannot represent a non-convex face, so the notch at the
      bend belonged to neither half. `splitAtReflexVertex` splits such a half along a diagonal. NOTE the
      retracted reasoning below -- equal areas and a bit-identical shared polyline do NOT imply two
      halves tile, because area says nothing about which side of a BENT boundary a point falls on.
      History follows. -- **both fixtures now ASSEMBLE** (three defects
      fixed: a duplicated consecutive vertex shifting the curve labels, `numel` where
      `size(...,1)` was meant in two index guards, and a STRAIGHT splitting curve being routed
      through the two-arc split, which flattened the inherited arc). Fixture 2 is exact on all
      1080 ring points. Fixture 1 has **one** point left, and it is a COVERAGE HOLE, not a wrong
      value: at `s = (0.998629534754574, -0.0523359562429444)` no face admits, and the three
      nearest faces all carry the right quadratic and miss by 0.0019, 0.0052 and 0.0072 in
      normalised conic units -- far above any tolerance. The point belongs to cell (2,1)
      (`g1` face 2 n `g2` face 1); both pieces from that cell exist (pieces 3 and 4) and neither
      covers it. **MEASURED, and it exonerates the two-arc split:** the two pieces from cell (2,1)
      share their cut polyline BIT-EXACTLY (both carry `M = 1.1254915141491897,
      0.074667480226358884` and the same two endpoints, to all 17 digits), their orientations are
      both CCW, and their polygon areas sum to the parent's exactly (0.03699 each way). So the
      pieces tile. The slit therefore appears AFTER the split -- in assembly: the vertex merge, the
      half-edge pairing, or the face edge-lists `orderEdges` builds. Next: take the assembled faces
      that miss the point (they carry the right quadratic and miss by 0.0019, 0.0052 and 0.0072 in
      normalised conic units) and compare each one's `P` list against the piece it came from.
      OLD NOTE, now fixed, kept for the trail: the orphan error reported a STRAIGHT boundary edge
      of piece 4 facing an IDENTICAL CURVED edge of its neighbour, at distance 0.
        * fixture 1: piece 4 `src[2 1]` (1.297862,0.278742)->(0.915534,-0.078641) straight, versus
          piece 5 `src[2 2]` curved.
        * fixture 2: piece 4 `src[1 1]` (1.163109,-0.285096)->(-1.244161,-1.161034) straight,
          versus piece 6 `src[1 2]` curved.
      So a piece that must carry TWO arcs -- its own operand's, plus the other operand's arc along
      the shared face boundary -- has one of them represented by its chord. matchHalfEdges then
      correctly refuses to pair a straight half-edge with a curved one, and the arrangement fails.
      Piece 4 has nv=4 with an arc already on edge 1, so the flattening happens where the SECOND
      arc is introduced: `clipPolyByConic` (including its corner-cut branch) or
      `clipPolyHalfPlaneCurved`, not in `splitTwoArcPiece` -- three fixes there (index re-location,
      degenerate-piece filter, polyline cut when no chord is interior) each held and none changed
      this symptom. Next: instrument where piece 4 acquires its arc, and check whether the cut
      conic was applied as an EDGE or silently dropped to a chord.
- [ ] `splitTwoArcLens` refuses when the cut `A -> M -> B` leaves the cell (the seeded far-field
      fixture). The two arcs there join corners on OPPOSITE branches of their parabolas, so the arc
      between them swings far out and the polyline exits the cell; a different subdivision is needed.
- [x] **FIXED** `twoCurvedWhereTheSplitCurveCrossesAnArc` -- the test passes, and `MAXQP_ASSERT=2`
      is clean on that fixture. The four non-compact-arc-piece findings this entry recorded were
      closed by `QuaPar.chordCuts` (2026-08-13) and the corrected chord derivation in
      `pieceRecessionRays` (2026-08-14); the entry outlived them. Original text follows.
      -- 2 of 68 sample points wrong, and it is the same defect as the item below, not a separate
      one. At both bad points `QuaPar.eval` reported `region 0`, i.e. SEVERAL faces admitted them
      (`(-3.9811,0.6115)` gave 0.468 and `(-5.0954,0.1351)` gave 1.229, truth 0 at both); the
      verifier named face 13, carrying g1 face 2, beaten by g2 faces 5 and 6 by +inf along one of
      its own rays; and `MAXQP_ASSERT=2` listed four pieces with non-compact constraint regions
      (`src[2 1]` and three `src[6 1]`), all BOUNDED arc-pieces.
- [x] **Step 4 of the plan: bounded arc-pieces whose CONSTRAINT region is non-compact** --
      **the REPRESENTATION question this item poses was answered on 2026-08-13, and the item as
      written is stale.** Neither (A) nor (B) below is what happened: option B was checked before
      implementing and refuted (the chord runs through the NEIGHBOUR's interior, so making it a
      real edge splits the neighbour and leaves the offending face's own edge list unchanged -- see
      `DECISIONS.md`). The answer is that the chord is **derived per face**, which is what
      `QuaPar.chordCuts` does, and it resolved the whole far-field defect.
      What was left of this item on 2026-08-14 was that `pieceRecessionRays` -- the piece-level
      analogue, which decides the same question before assembly -- was still using the weaker rule:
      it read the chord's side off the piece's other VERTICES and had no gate on when a chord may
      be emitted at all. Corrected to read both off the conic. **Unrun**; the residual
      `MAXQP_ASSERT=2` findings on `src [1 2]`, `[1 6]` and on
      `twoCurvedWhereTheSplitCurveCrossesAnArc` are the measurement that closes this.
      Original text follows, for the record: for a piece whose arc is CONCAVE towards it, the
      constraint set is a wedge intersected with the OUTSIDE of a parabola; a cut parallel to the
      chord leaves the arc-side sliver still receding along the two side edges' own direction,
      which neither the cut nor the parabola blocks. The two options considered were (A) let a face
      list a REDUNDANT bounding conic in its `P`, and (B) make the chord a REAL edge by splitting
      the neighbour across the arc.
- [x] The verifier does not prove the faces COVER the plane; `partitionReport` only samples.
      **`verifyFacesCoverThePlane` (2026-08-14) does, in four checks on the constraint data:** every
      edge separates two faces; every edge lies inside both of them; no face's constraint region
      has boundary anywhere but on its own edges; no face is squeezed onto a curve. Together they
      force the boundary of the union of the faces to be empty, so the union is the plane. The
      argument is in the file's header, the summary in `SUPPORT_MATRIX.md` 4.3. **Unrun.**


_Seeded 2026-08-02 at the start of the overnight run, from the task given when it
was launched. The repository had no TODO.md; the acceptance criterion is precise
(three named tests green), so the run works from this list._

## 2026-08-13: four arc-vs-arc failures pinned from ORDINARY polygon splits (not the shift fixture)

The arc-vs-arc defect does not need the translated-triangle fixture: `f = x*y` over a quadrilateral
handed in as two triangles either side of a diagonal reaches it, with a closed-form reference that
owes nothing to the pipeline (`supBilinearOverPoly` over each triangle; `sup` over a union is the
max of the sups, so overlapping triangles are fine). All four are now RED tests in `maxQuaParTest`,
vertices written out literally:

- `unitSquareSplitByItsDiagonalIsExactNearTheArc` — the SMALLEST failing case: `x*y` on `[0,1]^2`
  split by the main diagonal. Truth is closed form and POLYHEDRAL (`max(0,s1,s2,s1+s2-1)`, since the
  objective is bilinear on a box), so the two operands' arcs must cancel in the max; they do not.
  Wrong at 17 of 1080 ring points, ALL at radius <= 1 — worst 0.437 at `s=(0.45399,0.891007)`, where
  it returns the `s1` face (bounded by the arc conic `(s1+s2)^2/4 - s1 = 0`) in `s2` territory.
  So on this fixture the over-extension is NEAR the arcs, not in the far field.
- `arcVsArcIsExactFarFromTheArcsOnASeededQuadSplit` — the FAR-FIELD form, worst 5.28e4 at
  `s = 100*(cos(-57deg), sin(-57deg))`: got 53126, truth 309. Error growing like `|s|^2` is the
  non-compact bounded arc-piece signature of the MAJOR FINDING below.
- `arcVsArcDoesNotCrashOnSeededQuadSplits` — two fixtures, both `MATLAB:badsubscript` at the same
  site, `splitTwoArcPiece` line 2216 (via `clipPolyByConic` -> `clipByFace`). A raw badsubscript is
  never a designed refusal, so this is a bug independent of how the far-field work turns out.
- `arcVsArcRefusesAnUnboundedTwoArcSplit` — `maxQuaPar:notImplemented` from `splitCell`: "an
  UNBOUNDED half carries both the inherited arc and the splitting curve". Per FARFIELD_FIX_PLAN.md
  Phase 3 that guard is a bug detector, not a supported-input error.

`maxQuaParTest` is now 18 pass / 6 fail (the 2 pre-existing reds plus these 4), 40 s, still fast.

## MAJOR FINDING 2026-08-04: arc-vs-arc results are only LOCALLY correct (wrong far from the arcs)

The random-shift sweep's "silent WRONG" is NOT a handful of edge cases -- it is PERVASIVE, and it
afflicts even the two pins this session marked FIXED. Measured directly (g.eval vs pointwise-max on
rings around the origin):
  - `[-1 0.75]`: 0/60 wrong on a radius-8 ring, **2/60 wrong on radius 30** (worst 6.1).
  - `[2 -0.5]`:  2/60 wrong on radius 8,     **11/60 wrong on radius 30** (worst 33.8).
Their suite tests pass only because `curvedSamplePoints` samples NEAR the arcs. So `maxQuaPar` on two
curved operands is correct locally and WRONG in the far field, generally.

ROOT: a quadratic conjugate face is genuinely UNBOUNDED (e.g. g1 face 1 carries rays). The arc-vs-arc
subdivision, though, emits the sub-pieces of such a face as BOUNDED arc-pieces (one parabola arc plus
straight edges, no rays). A bounded piece left on the parabola's OPEN side is not a compact QuaPar
face: `QuaPar.eval` (locate by "every bounding conic, sign-oriented, <= tol") then admits points
arbitrarily far away into it, so it OVERLAPS the true face out there and, by eval's last-admitter-wins
rule, returns the wrong value. Verified on `[0.5 0.5]`: assembled FACE 15 carries g1-face-1's quadratic
over a tiny triangle near (-2,2) but admits (-3.98,0.61) two units away, overlapping the correct zero
face (FACE 9).

TRIED, reverted (all four break the suite or don't help):
  - piece-level compactness guard on triHalf output -- WRONG sign convention, rejected valid faces.
  - piece-level compactness guard on every bounded curved piece -- same, errored on all 3 fixtures.
  - post-assembly "bounded face admits a far point" check (QuaPar's own EC + P signs) -- CORRECT
    detection, but errors on `[-1 0.75]`/`[2 -0.5]` too, because they ARE non-compact far out.
  - post-assembly "two faces disagree at a far point" check -- also errors on all three, since the
    disagreement is real (they are wrong far out). Confirms the issue is systematic, not per-fixture.
A safety backstop that refuses non-compact results is thus correct but would turn every arc-vs-arc
result RED -- so the real fix must be UPSTREAM: give an unbounded quadratic face's arc sub-pieces their
RAY boundaries (so they stay unbounded and compact-as-faces) instead of closing them with straight
edges. That is a rework of the arc-vs-arc clip/split, not a guard.

NB: this reframes the whole session. The six arc bugs fixed are real and make the results correct
NEAR the arcs; but "18/2, two pins fixed" means "fixed where the tests sample", not "fixed".

## Genericity of the arc-vs-arc fixes (measured 2026-08-03)

`arcVsArcMatchesGroundTruthOverRandomShifts` sweeps seeded random shifts (not the 3 hand-picked
ones). Over 60 shifts: **65% exact, 15% assemble-to-WRONG, 20% error.** So the fixes generalise
well (not ad-hoc), but the 15% wrong is the TOP next task: it is the SAME pre-existing far-field
over-extension as [0.5 0.5]'s residual 2/68 (a decided unbounded polyhedral cell reaching past its
g1 face near a mesh vertex), which the arc-assembly fixes newly EXPOSE by letting those cases
assemble instead of erroring. conj's verification catches it in production. Fixing that far-field
coverage bug should turn most of the 15% wrong (and [0.5 0.5]'s 2/68) green.

## Status 2026-08-03 (Opus 4.8) -- 18/1, [0.5 0.5] now assembles

Three arc-vs-arc pins were red; **TWO fully fixed**, the third now ASSEMBLES (was erroring) and is
red only by a 2/68 value error from a SEPARATE pre-existing far-field bug. `maxQuaParTest` 16/3 ->
**18/1** throughout, no regression (random-quadrilateral sweep + all arrangement tests green).
Commits on `overnight/2026-08-02`: `96aad61` (T-junction), `53fc9fd` (assignSide), `be1a31f`
(arcEdge off-by-one + escape-to-infinity split + dedup normalisation), `3e1a6b2` (triangle
two-arc split). The [0.5 0.5] fixture chained ~9 distinct arc-handling bugs; six are fixed.

- **[0.5 0.5] `twoCurvedWhereTheSplitCurveCrossesAnArc` -- now ASSEMBLES, red by VALUE (2/68).**
  Six bugs fixed to get here (see the four commits above); the remaining error is UNRELATED to the
  arc machinery. The 2 wrong samples are `(-3.98109,0.61148)` and `(-5.09537,0.13508)` -- both far
  out in the lower-left, both essentially ON g1's mesh vertex `V5=(-3.98249,0.610504)` and its
  face-3/face-4 boundary ray (dir `(-0.857,-0.514)`). Piece `src[4 3]` (a DECIDED unbounded
  polyhedral strip, untouched by any session change) over-extends a hair past g1's edge 4 there and
  evaluates g1-face-4's quadratic (~0.47/1.23) where g1's real value via the neighbouring face is 0.
  A PRE-EXISTING far-field coverage/precision issue at a g1 vertex, newly reachable now that assembly
  completes; conj's own verification would catch such a result in production. NEXT: why
  facePoly(g1,4)/clipByFace lets src[4 3] cross g1's edge 4 near V5. The long "arc-clip mirror"
  analysis below is SUPERSEDED -- the mirror (src[6 2]) is now built correctly.

- **[-1 0.75] `twoCurvedThatMustAssembleAcrossRays` -- FIXED.** Was a cross-piece T-JUNCTION:
  a decided cone's ray ran to infinity while the neighbour side was split (segment+ray) at a
  point P lying exactly on that ray (perp ~2e-15). `insertGlobalPassthrough` now re-inserts
  every piece vertex that lands on another piece's ray/segment; a companion fix to
  `raySideVector` (its old adjacent-vertex representative is COLLINEAR with a just-subdivided
  ray, so `oppositeSides` decided on sign noise -- now falls back to the other ray's direction).
  Green and exact at all 68 ground-truth samples.

- **[2 -0.5] `twoCurvedAssembleAcrossRaysSecondFixture` -- FIXED (`53fc9fd`).** Was NOT decideWinner
  (that correctly flags the cell undecided) -- it was `assignSide` reading the winner at `piece.V(2,:)`,
  which for splitCell's UNBOUNDED "rest" piece (curveAfter=1) is a CROSSING point where diff~0 and its
  sign is noise, so both halves could get the same winner. Now reads the vertex farthest from
  `{diff=0}`. Exact at all 68 samples. (The long note below about decideWinner/parabola-to-infinity was
  the WRONG lead -- the real cell here is a strip cut by a LINE at two finite boundary points.)

- **[0.5 0.5] `twoCurvedWhereTheSplitCurveCrossesAnArc` -- STILL RED (errors). NEW, sharper
  diagnosis; the old "piece 5 over-extended" note below is WRONG.** Piece 5 (`src[2 2]` =
  g1f2 n g2f2, an **arc-vs-arc** clip via `clipPolyByConic`) emits an unmatched ray on `x+y=0`
  (g1's face-2/face-6 edge) from apex `(-2.03125,2.03125)`, dir `(-1,1)`. The clip is CORRECT, not
  over-extended: the sign data at that cell is `evalConic(Ecut)@V = [0, -0.046, -0.015]`, so the
  vertex `(-2,2)` sits on the g2-face-1 (discard) side and the far ray is genuinely g2-face-2. What
  is missing is the MIRROR piece across `x+y=0`: `src[6 2]` = g1f6 n g2f2. g1 face 6 does NOT carry
  the arc (the arc edge borders g1 faces 2 and 1 only), so `src[6 2]` is clipped by the
  **arc-vs-HALF-PLANE** path (`clipPolyHalfPlaneCurved`), a DIFFERENT code path -- and it does not
  produce the matching `x+y=0` ray (its `src[6 2]` pieces 15/16 sit at `x+y=0.5` and `1.0`, and
  piece 15 has a bounded vertex exactly at `(-2.03125,2.03125)` but no ray along `(-1,1)`). So piece
  5's ray has no neighbour because the two clip paths DISAGREE along the shared g1 mesh edge.
  Post-passthrough T-junction scan finds nothing collinear with the ray, confirming it is not a
  subdivision miss. FIX DIRECTION (not yet done -- delicate, high regression risk to the arc
  machinery): make the arc-vs-half-plane clip of g2f2 by g1f6 produce the same `x+y=0` boundary the
  arc-vs-arc clip of g2f2 by g1f2 does, i.e. reconcile `clipPolyHalfPlaneCurved` with
  `clipPolyByConic` along a shared straight g1 edge.
  DECISIVE CHECK NOW DONE (confirms "f6 clip under-producing"): g1's mesh for [0.5 0.5] has
  `E(6,:)=[2 7 0]` a RAY from V2=(-2,2) in dir (-0.707,0.707) lying on `x+y=0`, with `F(6,:)=[2 6]`
  -- so it is exactly the g1 face-2/face-6 boundary and it runs to INFINITY. Therefore `src[6 2]`
  (=g1f6 n g2f2) must have a matching `x+y=0` ray, and `clipByFace(facePoly(g1,6), facePoly(g2,2))`
  -- which swaps to clip the curved g2f2 by g1f6's straight half-planes -- is dropping it. Start in
  `clipPolyHalfPlaneCurved` / `clipByFace`'s swap path: the g2f2 boundary ray on `x+y=0` is being
  lost when g2f2 (curved) is clipped by g1f6's `x+y=0` half-plane (edge 6). g1 face 6 is a cone at
  (-2,2) between edge 6 (dir (-0.707,0.707), on x+y=0) and edge 7 (dir (0.707,0.707)).

  --- OBSOLETE (kept for history; superseded by the sharper diagnosis above) ---
  The old `decideWinner`/parabola-to-infinity write-up for [2 -0.5]. The T-junction
  fix lets assembly COMPLETE, uncovering the real defect: `decideWinner` wrongly declares an
  UNBOUNDED cell "decided". Traced with a pre-assembly coverage probe: cell `src[4 3]` stores
  `winner=g2`, but at `s=(-5.4843,1.5866)` -- an interior point of that cell -- `g1=0 > g2=-8.648`,
  so g1 wins there; symmetrically cell `src[6 5]` stores `winner=g1` while g2 wins at
  `s=(-1.1298,4.1007)`. 4/68 samples, worst 8.648.
  WHY `decideWinner` misses it: it proves domination by sampling finite vertices + the ARC midpoint
  + the ASYMPTOTIC sign along the two bounding rays (t->inf). But `diff = f1-f2` is a PARABOLIC
  quadratic (rank-1 Q -- splitCell already asserts this), so `{diff=0}` is a parabola whose two
  branches run to infinity ALONG Q's null direction, which lies strictly BETWEEN the cell's two
  bounding rays. That parabola bounds a region of the opposite winner entirely inside the cell,
  touching neither the finite boundary nor either bounding ray -- so every point `decideWinner`
  samples is on one side and it "decides" wrongly.
  Tried and REVERTED (kept out of the commit): sampling `diff` at the 1-D stationary point of each
  edge and each ray. Sound and strict, but it does NOT fix this -- the sign change is off the
  1-skeleton entirely. Confirmed neutral on all three fixtures.
  Note this cell is ALSO beyond `splitCell`: even if `decideWinner` returned undecided, `{diff=0}`
  here makes ZERO finite boundary crossings (both parabola branches escape to infinity inside the
  recession cone), and `splitCell` asserts exactly two. So a real fix needs BOTH a sound
  sign-over-the-cone test in `decideWinner` AND a `splitCell` that can cut a cell along a parabola
  that enters and leaves at infinity. Same family as next-step 2 in the session handoff (the
  Step-3 unbounded over-claim); do not attempt a probe-based patch.

## Next up (the [0.5 0.5] defect) -- SUPERSEDED, read the Status section above first

The "over-extended, should terminate" framing below is REFUTED by the sign data (see Status):
piece 5's ray is legitimately g2-face-2; the real bug is the missing g1f6 n g2f2 mirror from the
arc-vs-half-plane clip path. Kept below only for its geometric measurements.

- [ ] **Piece 5 (src `[2 2]`) emits a RAY where its boundary should terminate.**
      Localised to the cell, the edge and the reason; this is the last defect.

      * Its unmatched ray: apex `(-2.03125, 2.03125)`, direction `(-1,1)`, lying
        on the line `x+y=0`, which is g1's face-2/face-6 edge.
      * CORRECTION to an earlier note here: piece 16's ray is NOT its partner.
        That one lies on `x+y=0.5`, a parallel but different line, so the two were
        never meant to match.
      * Sampling across the apex, the local structure is three cells:
        `(2,2) | (6,1) | (6,2)`, with g2's ARC separating the last two. So piece
        5's neighbour along `x+y=0` is a `(6,1)` cell — and the `(6,1)` pieces
        (13 and 14) are BOUNDED slivers of area 0.008 and 0.004.
      * A bounded neighbour means the shared boundary is a finite SEGMENT, not a
        ray. So piece 5 is over-extended: its boundary along `x+y=0` should stop
        where g2's arc crosses that line a second time, and instead it runs to
        infinity. `matchHalfEdges` pairs rays with rays and segments with
        segments, so a ray facing a segment can never match — which is exactly
        the reported symptom.

      Next: find why `clipPolyByConic` cuts cell `(2,2)` at the first crossing
      with g2's arc but not at the second. Note the restriction of that conic to
      this ray came out with `A = 1.7e-18` — numerically degenerate, so the
      quadratic is treated as linear and yields ONE root. Check whether the true
      second crossing is being lost there, or whether it lies on a different
      boundary element of the cell.

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

# TODO

## 2026-08-24 — SYM-FREE `conj`: what is left, measured. READ THIS FIRST.

The premise changed on 2026-08-24, and it changes the whole plan below: **vertices are stored as
INTERSECTIONS OF CONICS with floating-point coordinates, not exactly.** By
`CONJ_FIELD_PROOF.md` Theorem 1 the face functions and the edge conics are always rational, so the
only thing that ever needed an extension field was the vertex layer -- and naming a vertex by the
pair of conics it solves removes the need for one entirely.

**Two items below are CANCELLED by that premise, not deferred**: the degree-<=4 real algebraic
kernel, and the interim "factor the vertex quartic and refuse" gap. `exactQ` is likewise not on the
path to anything; it is referenced by nothing but its own test.

### What was measured, and it is the thing to know

The numeric conjugate path was **already 100% sym-free** before this work: `conjCPLQ`,
`conjPieceCPLQ`, `convEnvCPLQ`, `maxQuaPar`, `clipArcByConic` and `mergeSameQuadFaces` contain
**zero** `sym`/`subs`/`simplify`/`solve`/`isAlways`/`coeffs` between them, counting non-comment
lines. Every symbolic call on the conj route is behind ONE dispatch -- the fallback to Case C. So
the work is not "rewrite Step 3"; it is **shrink the set of inputs that fall back**, and
`checkConjSymFree` measures that set with the reason for each.

    2026-08-24 baseline, 16 fixtures:   SYMBOLIC 2 of 16, both maxQuaPar:notImplemented
    after conjConvexPolygon landed:     the unbounded CONVEX family joined the numeric route

### The gaps that remain, in the order they should be closed

- [ ] **G1 -- a missing LENS in the overlay.** `maxQuaPar` can now SPLIT a cell whose arc is cut
      twice, or whose arc carries a neighbour's vertex in its interior (`bulgeSplit` / the
      passthrough split, 2026-08-24), and `maxQuaParTest` is green with it. The two fixtures that
      needed it then die two stages later in `assemblePieces`, and the piece dump says why:

          piece 1 src[1 1] carries g2's parabola from (0,0) to (1,1);
          its would-be neighbours are STRAIGHT edges through (0.5,0.5),
          which is not on that parabola.

      MEASURED with each operand's own point location: every point of the lens between g1's straight
      diagonal and g2's arc lies in **g1 face 4 and g2 face 2**, and the fold produces **no piece
      with src [4 2]** -- while the control point one step away, in [1 2], is there. So
      `clipByFace` returned nothing for a pair whose intersection is not empty, and the orphaned
      arc in `assemblePieces` is the symptom rather than the defect. `polyConstraints` already
      skips a curved edge's chord, so the candidates are the operand SWAP at the top of
      `clipByFace`, `clipPolyByConic`, and the three reduction passes. **Instrument that one pair.**
      `DECISIONS.md` 2026-08-24 (later) has the numbers; `SUPPORT_MATRIX.md` 4.1 records the same
      family of defect from 2026-08-13.
      Closing this removes the LAST measured fallback of the bounded family.

- [ ] **G2 -- an AFFINE face over an UNBOUNDED polygon.** `max(0,x,y)` is the canonical example and
      it still falls back. Two things are missing and they are separable:
      1. the construction: `(L'x + c + I_P)*` is the support function `sigma_P(s-L) - c`, which for
         `P = conv(W) + cone(d1,d2)` is `max_i <s-L, w_i> - c` restricted to the polar cone
         `{t : <t,d1> <= 0, <t,d2> <= 0}` -- a QuaPol with affine faces and an INDICATOR domain;
      2. `maxQuaPar.assertFullDomain` refuses any operand that is not finite everywhere, and every
         piece of this family is `+inf` off a cone.
      **Worth pricing a third route first.** For an input whose pieces are ALL affine, `f*` is the
      max of finitely many affine functions restricted to one polyhedron -- the upper envelope of
      planes, i.e. a 3-D hull computation -- so a direct `conjAffinePLQ` would cover every
      piecewise-LINEAR input in one construction and never enter `maxQuaPar` at all. That is a new
      operator rather than a gap fix, which is why it is not started.

- [x] **G2b -- DONE 2026-08-24. `maxQuaPar` dropped a cell on some unbounded folds.** Found 2026-08-24 and it is the
      one silent wrong answer in the session: a 4-cone fan assembled to 4 cells and returned 2.0 at
      `s=(-2,-3)` where the definition sup is 4.5. Reproducer:
      `conjCPLQTest/step3UnboundedAssemblyMatchesTheTruth`'s fixture in its `F = [3 2;2 1;1 4;4 3]`
      orientation; `PARTITION OVERLAP` fires on it. A cross-check now catches this and falls back
      **Fixed the same night.** The cause was one line: at the SINGULAR POINT of a degenerate
      conic (a pair of lines crossing) the gradient vanishes, and
      `splitUnboundedAtOneCrossing` took its branch direction as the perpendicular to that
      gradient. It bailed, and the caller's tangency branch then read the winner at the cell's
      centroid -- which for a cone IS that same point, on the curve. The branches are now taken
      from the NULL DIRECTIONS of the quadratic part. `DECISIONS.md` 2026-08-24 (last) and
      `SUPPORT_MATRIX.md` section 4.5 have it; the two-line reproducer is
      `maxQuaParTest/twoHalfPlaneQuadraticsSplitTheirSharedQUADRANT`. The cross-check that caught
      it stays.

- [ ] **G4 -- `conj` of `xy` on some triangles is wrong, and it is the FOLD.** RE-MEASURED
      2026-08-25; the description below it replaced was wrong on both counts, so read this one.
      Step 1 splits sweep case 21's triangle into **four** faces (the nCE==3 cevian split with each
      half re-split), two of them slivers of area 8.7e-05 against 2.7e-02. **Every face's own Step 2
      conjugate is EXACT** at the bad point -- 1.032507658472 to twelve digits, four times over --
      and folding faces 2 and 3 keeps it; the FOURTH fold returns 1.005089907622. Pairwise folds of
      any two faces are all exact: only folding a sliver into the ALREADY-ACCUMULATED mesh loses it,
      which is why every two-operand `maxQuaParTest` passes.
      **It is also far worse than 2.742e-02**, which was only the worst on a probe grid of radius
      <= 6: the fold cross-check finds `f*(-10,0) = 47.10181578` against a true 10.86895777, an
      OVER-estimate by a factor of four. The G6 edge bound is one-sided and cannot see that.
      **Contained, not fixed.** `conjMaxOfSubTriangles` now cross-checks the fold against its own
      pieces (the identity `f* = max_k (q_k + I_T_k)*` makes them their own oracle) and refuses by
      name in 2.5 s as `PLQ:conjCPLQ:foldDroppedACell`, under `CCA2_CONJ_VERIFY` like the edge
      bound. It REFUSES rather than falling back because Case C did not finish in 25 minutes on it.
      **What is left is the `maxQuaPar` defect.** With `MAXQP_ASSERT = 2` this input raises two
      invariant violations, both real: a piece whose declared rays are the NEGATIVE of the direction
      its constraint region recedes along, and a piece carrying one operand's quadratic where the
      other is larger by `Inf` along a ray. The second is what produces the wrong value.
      `DECISIONS.md` 2026-08-25 (G4).

- [x] **G5 -- DONE 2026-08-25. `MATLAB:badsubscript`, and it was never 5-gon-specific.** The old
      entry said "SOME indefinite 5-gons" from sweep case 17; a 40-case re-run found case 29, a
      **4-gon**, crashing identically, so a fix aimed at 5-gons would have missed half of it.
      `plq_1p.conjugate`'s FRAME-CHANGE branch copied the z-frame object's ENVELOPE (2 faces) while
      replacing its `conjfia` -- the per-face block boundaries into `conjugates` -- with the single
      block `[1 nConj+1]`; `maximumConjugate` sizes its loop from the envelope, so face 2 asked for
      `conjfia(3)` of a 2-element array. Measured: `envelope faces = 2, numel(conjfia) = 2,
      nConjugates = 11`. The branch now carries `objT.conjfia`, and `maximumConjugate` takes its
      block count from `numel(conjfia)-1` -- the two arrays are not interchangeable, and the
      SEPARABLE branch legitimately returns one block over a multi-face envelope.
      **What the crash was hiding:** with the index repaired the cross-face symbolic max this input
      always needed actually runs, and `maximumP` on its two envelope faces did **not finish in 25
      minutes** -- measured both in the z-frame and after the read-back, with no difference. So case
      29 goes from a fast wrong crash to a correct computation nobody can wait for; the route that
      would make it fast is the NUMERIC one, which declines it today with `maxQuaPar:notImplemented`
      (clipPolyByConic separating an unbounded cell). Pinned by
      `conjCPLQTest/frameChangedPieceKeepsItsEnvelopeBLOCKS`, which asserts the INVARIANT rather
      than the value, because a value assertion here is a test nobody can run.
      `DECISIONS.md` 2026-08-25 (G5).

- [x] **G6 -- DONE 2026-08-25. The EDGE lower bound is a DEFAULT refusal.** Along each
      edge of `dom f` the objective `<s,x> - q(x)` is a quadratic in the segment parameter, so its
      maximum is closed form, and `f*` must be at least that. On the 24-case random sweep it fires
      on **exactly** case 21 (G4, at -2.742e-02) and on nothing else -- every other case sits at
      ~1e-15. The VERTEX-only version of the same idea catches NOTHING (0 of 24, including case
      21, because that sup is not attained at a vertex), so build the edge one.
      It is cheaper than the fold cross-check already in `conjPolygonalDomain` and, unlike it,
      covers the SINGLE-PIECE route. Built as `conjEdgeLowerBound.m`, raised by `conjCPLQ` as
      `PLQ:conjCPLQ:belowEdgeBound`, and ON by default (`CCA2_CONJ_VERIFY = 0` turns it off).
      Gated: fast 303/0, slow 88/0, verylong 26/7/1 -- the last IDENTICAL to a pristine `b9243d3`.
      It RAISES rather than falling back because the symbolic route returns the same wrong value on
      the known-bad case. `DECISIONS.md` 2026-08-25 (final).

- [x] **G7 -- DONE 2026-08-25. `plqStage`'s cache raced under `--verylong -j N`.** `save` now
      writes a unique temporary name in the SAME directory and `movefile`s it onto the real one, so
      a reader sees either the old file or the complete new one; and `load` is wrapped so an
      unreadable or half-written cache is treated exactly as a missing one -- it recomputes.
      NOT verified under contention: that needs the contended `--verylong` run this session was
      asked not to make. The change is one-directional (it can only turn a throw into a recompute),
      and the uncontended path is exercised by every staged test.

- [ ] **G3 -- a non-convex face over an UNBOUNDED polygon.** Declines by name today
      (`the fan-triangulation route needs a BOUNDED face`). Needs Step 1 or Step 2 for an unbounded
      indefinite piece; `convEnvUnbounded`/`fanUnboundedFace` are the symbolic side's answer and
      are the only remaining symbolic helpers reachable from a conj-shaped input.

- [ ] **G4 -- `QuaCon` storage (H-form).** Now optional rather than blocking, and that is the
      premise's doing. `Conic`/`RatCon` landed 2026-08-24 so the LATTICE can hold an elliptical
      edge; `ratQ` and `conicMeet` are the exact coefficient layer and the vertex-naming primitive
      it would use. Build it when a `conj` result actually needs a non-parabolic edge -- i.e. when
      G1 lands and a three-piece input with two non-adjacent pieces reaches Step 3.

- [ ] **G8 -- `testPCE2`'s CONVEX ENVELOPE is wrong, one stage before the conjugate.** Found
      2026-08-25 by splitting `testMaxMultiRegion`'s PCE family by stage: on the triangle
      {(0,0),(1,0),(2,1)} carrying `x*y`, `pce2EnvelopeUnderestimates` fails in **26 s** --
      Step 1's envelope does not underestimate `f` / touch it at the vertices. The known red
      everyone quoted for this fixture was the CONJUGATE one (`f*(0,0) = 0.0429` against a sup of
      0); `pce2ConjugateMatchesItsOwnSup` reproduces it in **36 s** and reports it as
      `f*(1,1) = 1.125 BELOW the sampled sup 1.24999381`. Fix the envelope first: a conjugate
      computed from a wrong envelope is not evidence of anything.

- [ ] **G9 -- SIX of the seven `verylong` reds are still unnamed.** No verylong log was ever kept
      and only `testMaxMultiRegion/testPCE2` is named anywhere in the repository, so "seven
      pre-existing failures" has never been a list. Naming them costs one `--verylong -j 1` run.
      Until that run happens no one can say whether any of the six is a `conj` defect.
      The PCE family is now split by stage (below), which is what makes each of them cheap to
      attack once named; `testBiconjugate` through `testBiconjugate8` are NOT split yet and are the
      remaining monolithic block -- each still runs `.maximum` then `.biconjugateF` inline in one
      method, uncached.

- [ ] **G10 -- the `maxQuaPar` accumulated fold, which is what G4 actually is.** Two invariant
      violations, both raised by `MAXQP_ASSERT = 2` on the G4 fixture and neither fixed: a piece
      whose declared rays are the NEGATIVE of the direction its constraint region recedes along,
      and a piece carrying one operand's quadratic where the other is larger by `Inf` along a ray.
      Pairwise folds are exact; only folding into an accumulated mesh fails. `SUPPORT_MATRIX.md`
      section 4.6's "assembly after an arc split: one side CURVED, the other STRAIGHT" is the
      design note for the same area.

### Tools built on 2026-08-24, so they are not rebuilt

- `checkConjSymFree.m` -- the fallback RATE and its reasons, per fixture. Run it before and after
  any change to the dispatch.
- `conjCPLQ(q, [], 'numeric')` -- refuses to fall back, so a test can pin the ROUTE. `conjSymFreeTest`
  does exactly that, including for the two gaps above (a gap test going GREEN is good news: promote
  it and delete the entry here).
- `conjConvexPolygon.m` -- a convex quadratic over ANY convex polygon, bounded or not, closed form,
  no triangulation, returns a QuaPol.
- `checkConjAgainstDefinition.m` -- randomized `conj`-vs-definition sweep over random polygons and
  all four sign classes. 22 of 24 exact; the two that are not are G4 and G5, both pre-existing.
  Slow by design (the reference is a scan plus a pattern search per probe), so it is a CHECK, not a
  bucket member. Run it after any change to the conjugate, and against a snapshot of the old tree
  when something fails -- that is what established G4/G5 were not introduced on 2026-08-24.
- `conjPolygonalDomain`'s fold cross-check -- the assembled Step 3 result is verified against
  `max_k (q_k + I_P_k)*`, the identity it was built from, and DECLINES on a disagreement so the
  caller falls back. One-sided by construction: it can miss a defect, it cannot invent one.
- `ratQ.m` / `conicMeet.m` -- exact rational coefficient arithmetic, and conic-conic intersection
  through an exact integer quartic.

---

## 2026-08-20 — THE SYM-FREE PORT, RE-PLANNED (superseded above, kept for its measurements)

The plan changed on 2026-08-20 because two things were measured, both in `DECISIONS.md` under that
date: one quadratic extension is not enough (a single A.5 triangle needs `sqrt(15)` and `sqrt(30)`
in ONE coordinate), and no tower of square roots is enough either (an `f*` vertex can have degree 4
with Galois group S4). Order matters below; each item says what decides it.

Terminology used throughout: **Row 7** is the recommendation of `CONJ_FIELD_PROOF.md` section 8.2 —
store the mesh in **H-form** (faces as sign conditions on rational conics, edges as rational conics,
a vertex named by the PAIR OF CONICS it solves rather than by coordinates) and run every predicate
through an interval **filter** first, dropping to the exact kernel only when the interval straddles
zero. That file is UNTRACKED at the time of writing — another session left it in the tree — so this
section states what it needs rather than pointing at it.

- [ ] **T3 — `symbolicFunction`'s payload becomes a coefficient vector, over RATIONALS.**
      The face and conic layers need no extension field at all (that document's Theorem 1: face
      functions and edge conics of `f*` are always rational), so `Rat`/int64 covers them and
      `exactQ` is not load-bearing there. This is where the engine calls actually are — live
      counts, non-comment: `subs` 72 in `region.m`, 68 in `plq_1piece.m`, 21 in `plq_1p.m`;
      `coeffs` 18 + 10; `simplify`/`simplifyFraction` 11 + 18 + 10; `hessian`/`gradient`/`dfdx`
      about 45. Every one becomes arithmetic on coefficients.

- [x] **CANCELLED 2026-08-24 by the H-form premise — a degree-<= 4 real algebraic kernel**: a rational quartic plus an isolating
      interval, sign of a rational polynomial at it, and comparison by resultant or Sturm sequence.
      `exactQ` is now multiquadratic and still NOT enough: of twelve continuous three-piece
      configurations the vertex quartic is irreducible over Q in ten of them, and the S4 case is
      proved reachable. Bounded work, because the degree is capped at 4 for both `conj` and the
      envelope.

- [x] **CANCELLED 2026-08-24 with the kernel above — DETECTED refusal.** Factor the vertex quartic and refuse
      by name when it does not split into rational or quadratic factors. That turns a reachable
      wrong answer into one nameable `SUPPORT_MATRIX.md` GAP, which is the discipline this project
      already has for unreachable branches.

- [ ] **Row 7 itself — H-form storage plus filtered predicates.** Biggest rewrite: V-form is baked
      into `RatPar`'s `V/E/F/P`, `eval`, `createP`, `orderEdges`, plotting, and every test that
      names a vertex. Plotting and user output still need coordinates, but then as an OUTPUT
      convenience, not the stored truth. The filter must be certified (a real error bound, not a
      tuned tolerance) with a terminating exact fallback, or it is the refuted double-plus-tolerance
      design wearing a disguise.

- [ ] **When Row 7 lands, INDEX the mesh for evaluation — do not linear-scan.**
      Measured (`.claude/evalbench.c`, gcc -O2 -march=native, min of three runs, value AND gradient
      every call; the caveats are in `DECISIONS.md`):

          baseline: 20-node expression tree, forward-mode gradient      44 ns
          linear scan          9 / 81 / 1024 cells             38 / 130 / 1670 ns
          uniform bucket       9 ... 1024 cells                24 ...    27 ns
          slab, two binary searches                            15 ...    26 ns
          cached cell + bucket, COHERENT queries                        11.5 ns
          no mesh, per-piece closed-form sup   3 / 6 / 24 pieces  105 / 215 / 822 ns

      Take the **uniform bucket plus a last-cell cache**: flat in the cell count, O(N) to build,
      and SCIP's queries move in small steps so the previous cell usually still holds. Skip the
      slab method — its 15-26 ns is measured on a grid of axis-aligned boxes, which flatters it,
      and on a real conic arrangement each of its ~log2(N) steps costs a conic evaluation. Keep
      unbounded cells in a separate short list, tested only when the point is outside the global
      bounding box. The index is pure preprocessing, microseconds to build.

- [ ] **T6 — delete `plq_1piece`** (57 live engine calls, 13 of them `solve`). RE-RUN THE FIXTURE
      SWAP FIRST: of the three regressions recorded on 2026-08-19, one may already be gone, since
      `testPCE2`'s domain is the triangle A.4 now gets right. The other two are specific —
      `testFractional` needs `conjugateFunction` to dispatch on the envelope's KIND (a rational
      face to the A.3 branch) and not only on its domain, and `testConvex` feeds a four-vertex
      polygon to a routine that wants triangles.

- [ ] **T4 / T5 — the remaining global sign questions and `solve()` sites.** `signEverywhere`,
      `impliedBy`/`holdsOn`, `certifiesNonPositive` are the only `isAlways` uses that ask about a
      FUNCTION rather than a number; degree <= 2 has a closed-form PSD/minimum test and half of it
      is already written. Live `solve`: 6 in `region.m` (most measured dead), 3 in
      `functionNDomain.m`, 5 in `symbolicFunction.m`.

**Facts about the SCIP side worth keeping, so they are not re-derived:** `f*` is convex, so the
gradient at any point is already a global linear underestimator, the max over a box is attained at
one of four corners, and the min is the global min when it is inside — the interval evaluation that
usually dominates propagation is nearly free. Evaluation is pure double precision: the exact kernel
is a BUILD-time cost and the finished mesh compiles to flat arrays, so the degree-4 vertex question
never reaches the solver. The gradient jumps across cell boundaries; the value does not, and either
side is a valid subgradient, so separation stays correct but an Ipopt-based heuristic will feel it.
Building the mesh is minutes to an hour, so it must be once-per-constraint preprocessing, never
reached from a node.

## 2026-08-23 — RETURN TYPES: `QuaCon` for `conj`, `AlgCon` for `biconj`

Why the current types cannot hold the answers is `DECISIONS.md` 2026-08-23 (envelope face type) and
2026-08-21 (`f*`'s elliptical edge). Both axes of the `RatPar` lattice grow by one level ABOVE the
present top, so nothing existing changes behaviour:

      subdivision:   Pol  <  Par  <  Con        `Con` drops b^2 - 4ac = 0
      function:      Qua  <  Rat  <  Alg        `Alg` = root of a rational quartic in z

- [x] **DONE 2026-08-24 (as `Conic`; `CON` is a Windows device name)** — a new trait plus a new
      data-holding parent `RatCon`. Cheapest
      item here, and the elliptical edge already forces it. `QuaPar` becomes a real specialization
      instead of a type that cannot hold the values `conj` produces.

- [ ] **`QuaCon` = the return of `conj`.** Rational faces `f` and edge conics `Ec` as int64
      `[a b c d e f | den]`, primitive and sign-normalized so equality is bitwise (the
      `4 - 2 sqrt 2`-as-two-doubles failure is what canonical integers prevent). Vertices become
      NAMES: `Vname(nv,3) = [edgeA edgeB rootIdx]`, the point where two edge conics meet, `rootIdx`
      canonical among the real intersections. `Vx` (double) and `Vbox` (rational isolating box) are
      CACHES — rebuildable, deleting them changes runtime only. This is Row 7; it is the expensive
      item, because `V`, `eval`, `createP`, `orderEdges`, plotting and every vertex-naming test read
      coordinates today. Do it on `QuaCon` alone; the four legacy types keep coordinate `V`.

- [ ] **`AlgCon` = the return of `biconj`, stored as a DECORATION of the `QuaCon`, not standalone.**
      Forced, not stylistic: an affine cell is `<p,x> - f*(p)` with `p` a dual vertex of degree <= 4,
      so its coefficients are irrational unless the cell NAMES the dual vertex. Payload:
      `src` (the `QuaCon`), `fkind` in {QUAD, RAT, RULED, AFFINE}, `fref` (dual cell pointer plus
      which adjacent face plays `i`), `fdeg` in {1,2,4}, `froot`; and `Ekind` in {CONIC, RULING}
      because the extreme rulings join two degree-<=4 points and their supporting line is rational
      only over `Q(p)`. `f`/`den` stay as caches for the QUAD and RAT faces.

- [ ] **`RatPol.m`'s header is now wrong** — "quadratic over linear on a polyhedral subdivision,
      proven not to need a square root" holds for ONE piece. Fix the comment when `AlgCon` lands;
      `RatPol` becomes `fkind = RAT, fdeg = 1`.

Evaluation stays pure double: locate the cell (uniform bucket + last-cell cache, the measured
choice), then dispatch on `fdeg` — coefficients, one square root, or one quartic. `grad f** = s` is
free and the Hessian is the rank-one closed form. The degree-<=4 exact kernel is a BUILD-time service
behind the interval filter and never reaches `eval` or SCIP.

---

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

### Next up (2026-08-18) — where Step 3 actually stands

**MEASURED end to end on the A.4/A.5 quadrilateral, after this session's work:**

    cells per fold   5, 14, 29, 45, 70, 86   ->   5, 12, 23, 38, 51, 60
    total            73 min                  ->   43 min (2579 s)

The machine was running three MATLABs throughout, so that timing is pessimistic, and a single
timing settles nothing anyway (see `CLAUDE.md` §3).

**What got it there, all measured:** three double leaks fixed so Step 2 is exact
(`domain.mE`/`cE`, `region.limitOfFAtVertices`, `plq_1p.quadPartsOf` + `conjConvexOverPiece`); a
sound certificate for a curved constraint (`region.certifiesNonPositive`) replacing merge's two
quadratic heuristics; `quadprog` deciding the CONCAVE conics, which is what the conics here
actually are; and `functionNDomain.singularEdgeCut` closing the singular-quadratic overlap.

**AND `testcPLQ/testRectBiconj` NOW PASSES with `CCA2_A45_SPLIT` on** -- `passed=1 failed=0
incomplete=0`, nothing changed in the test or the split. That exception was the stated correctness
blocker for making the split the default, and it was a casualty of the double leaks.

- [x] **1. THE A.4/A.5 SPLIT IS THE DEFAULT since 2026-08-18.** `testcPLQ` 8 passed / 0 failed in
      2188 s against 4728 s and one ERROR; full suite 332 / 0 with it on. `CCA2_NO_A45_SPLIT` opts
      out. The slow bucket went ~92 -> ~113 minutes, which is the bucket cost the standing rule
      says to accept.

- [x] **2. `region.impliedBy` over a region with a CURVED facet -- DONE 2026-08-18.**
      `region.holdsOn` + `maxAffineOverRegion` take the max over the region ITSELF (vertices, plus
      arc tangencies found in closed form) instead of over its linear relaxation. Measured on three
      folds: `quadFacet_exactAnotInB` 63 -> 41, fold-3 cells 38 -> 36, merges 7 -> 9.
      **An earlier revert of this was MY ERROR** -- two helpers in the instance block called as
      static -- see `DECISIONS.md`.

- [x] **3. `noSharedFacet` -- MEASURED 2026-08-18 and mostly HONEST.** At fold 3, of 137
      same-function pairs, 74 do not touch or touch at a single POINT, where merging would be
      wrong. Only 11 are genuine misses: 4 share an affine hyperplane the symbolic test does not
      match, 7 meet along a CONIC that neither facet search identifies.
      **The open question moved to `unionIsExact`**: 52 pairs reach it and about 9 merge, so ~43
      are refused by the sound gate -- and whether those are right is NOT established. Two cells
      can share a facet, touch along a segment, and still have a non-convex union. Check a handful
      of them directly before optimising, or risk reading a correct refusal as a defect for the
      third time this session.

- [x] **4. The parallelogram's piece 9 -- DONE 2026-08-18, `BAD 0 of 10`.** The singular-quadratic
      overlap was the only defect; the "remaining 1%" was a grid reference that missed the vertex
      where the sup is attained. See `DECISIONS.md`.

- [x] **5. `RatPar`'s `V (:,2){mustBeNumeric}` -- DECIDED 2026-08-18: leave it.** It costs cells
      and time, not correctness; the default path no longer goes through the mesh; and the change
      is lattice-wide against a design the classes state deliberately. `DECISIONS.md` records what
      would reverse it.

### Measurements that stand (2026-08-16)

- [ ] **STEP 3's CROSS-PIECE MAXIMUM DOES NOT SCALE, and it is now the binding cost.** Measured
      2026-08-16 on `x*y` over `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}`, which the A.4/A.5 split turns
      into 6 pieces: Steps 1 and 2 take about 25 s for all six, and `functionNDomain.maxOfList`
      then takes **73 minutes** -- folds of 93, 294, 647, 1273 and 2087 s, with the cell count
      running 5, 14, 29, 45, 70, **86**. `SUPPORT_MATRIX.md` records the same shape for a pentagon
      (885 s, 41 regions).
      **86 cells is roughly ten times what the answer needs.** `f*` of `x*y` over a convex
      quadrilateral has a cone per vertex and a cell per edge -- of the order of a dozen. The
      surplus is adjacent cells carrying the SAME function that are never merged, which is
      `region.merge` / `unionIsExact` refusing a union it cannot certify across a conic edge.
      **Where to start:** merge same-function neighbours after each fold, not only at the end. Many
      of these cells are POLYHEDRAL (the vertex cones), where `unionIsExact` already decides
      exactly, so a large fraction should collapse without needing the conic case at all. Measure
      the cell count per fold before and after -- that is the number that predicts the time.
      **This is the FIRST of two things standing between the A.4/A.5 split and being the default**:
      `testcPLQ` goes from 1542 s to 4728 s with the split on. The second is that `testRectBiconj`
      then ERRORS -- a separate, undiagnosed failure, and the reason the split is opt-in.

### Then (2026-08-15, after bug 1 closed the last red)

- [x] **THE GENERAL CONVEX QUADRILATERAL — FIXED 2026-08-16, on the fourth attempt, by doing the
      split SYMBOLICALLY.** `splitTightTriangleSym` splits a triangle into sub-triangles on each of
      which cPLQ's own closed form for THAT sub-triangle's convex-edge count IS the convex
      envelope, and `plq_1p.triangulate` emits them as PIECES.
      **Measured:** `f*` of `x*y` over `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}` is exact at 10 of 10
      probe points against the vertex-attained sup, with no piece leaving a hole; the fully
      assembled Step 3 answer is exact too, at 8 of 8. Pinned by
      `cplqAdapterTest/generalQuadrilateralStep1IsTheEnvelopeNotAMinorant` (the envelope must
      exist, be a minorant, and be `>= 0` where `x*y >= 0`) and
      `generalQuadrilateralConjugateMatchesTheSup`.
      **What made attempt 4 work where attempt 3 hung:** exact symbolic arithmetic. A.4's cevian
      foot is irrational, and taking it from `convEnvCPLQ`'s doubles gives `2^53` denominators that
      grow to `1e25` downstream; carried symbolically the coordinates stay compact surds
      (`5/2 - sqrt(5)/2`, `3/2 - 3*sqrt(5)/10`).
      **A COST, not a defect, and it is Step 3's:** assembling the cross-piece maximum for this
      input takes **73 minutes** (folds of 93, 294, 647, 1273, 2087 s, cells running 5, 14, 29, 45,
      70, 86). The per-piece conjugates take about 25 s in total. See the next item.
      **The cost is paid by any x*y polygon whose triangles need the split**, because A.4's cevian
      foot is IRRATIONAL: a split sub-triangle has SURD coordinates, and every symbolic operation
      downstream then works in a quadratic extension instead of the rationals. Measured on
      `testcPLQ`, whose domains are general polygons carrying `x*y`: **1542 s with the split off
      (matching its historical 1427 s), over 3100 s with it on and still unfinished when stopped**,
      uncontended both times. Only two of its six domains even gain a piece (2 -> 3 and 1 -> 2), so
      this is the algebraic degree of the coordinates, not the piece count.
      **The split is OPT-IN, via `CCA2_A45_SPLIT`, and OFF by default** -- and off for a measured
      reason: with it on, `testcPLQ` takes 4728 s instead of 1542 s AND `testRectBiconj` ERRORS.
      So it trades a documented, LOUD failure on one domain shape for a new one on another, and
      until that is understood it cannot be the default. The two tests turn it on themselves.
      History of the three failed attempts follows.
      **Attempt 3 (2026-08-16), and why it failed.** The domain split was built exactly as attempt
      2's write-up prescribed, taking the sub-triangles from `convEnvCPLQ`'s own faces, and it
      WORKED at Step 1. It then turned the crash into a **HANG** -- the first conjugate ran 45+
      minutes with no output, stuck in a symbolic call behind 3.8 MB of `isAlways:TruthUnknown`
      warnings carrying denominators around `1e25`. The cause is inherent to taking the geometry
      from `convEnvCPLQ`: that routine is double precision, `sym` of a double is EXACT (denominator
      near `2^53`), and snapping the new vertices to the simplest rational within `1e-10` bounds
      the VERTEX denominators but not the downstream ones -- the conjugate is a rational function of
      those coordinates, so a few squarings carry `1e5` to `1e25`. A hang is worse than a crash.
      **Attempts 1 and 2 (2026-08-15), for the record.** Wiring the missing `nCE == 3` branch alone
      leaves the answer wrong, because the `nCE == 2` branch returns a MINORANT and not the
      ENVELOPE (measured: its envelope reaches `-0.2835` on a triangle where the truth is `>= 0`,
      so `f*(0,0) = 0.28647` for a truth of `0`). Routing the envelope through `convEnvCPLQ`
      instead raises `symbolic:coeffs:NotAPolynomial`, because A.4/A.5's faces are RATIONAL and
      **cPLQ's Step 2 has no rational-envelope branch at all**.
      Original framing follows.
      What was on record: `p.conj('cplq')` raises `MATLAB:badsubscript` because
      `plq_1p.convexEnvelope1` branches on `nCE == 0, 1, 2` and falls off the end, so
      `obj.envelope` stays EMPTY and `conjugate`'s `for i = 1:max(1,size(envelope,2))` indexes
      `envelope(1)`. All true, and CCA2 already has the missing algorithm ([COAP] A.5,
      `convEnvCPLQ`'s `splitThreeConvex`).
      **The wiring was written and it WORKS at Step 1.** Build the triangle as a one-face `QuaPol`
      carrying `x*y` (safe: reaching that line with an indefinite `q` means `isCanonicalXY` held,
      so `q` is EXACTLY `x*y`), call `convEnvCPLQ`, convert with `ratPolToPlq` and install the
      faces as envelope pieces -- `plq_1p.conjugate` already loops over them. Measured on the
      offending triangle: **4 envelope faces** come back (the A.5 split, each half then needing
      A.4's), two quadratic and two rational, all `<= x*y`. `conj` no longer raises.
      **AND THE ANSWER IS THEN WRONG, so it was reverted.** `triangulate` splits the test
      quadrilateral into piece 1 = `[2.5 1.5; 2 0; 0 0]` (`nE = 2`) and piece 2 =
      `[2.5 1.5; 0 0; 0.5 1]` (`nE = 3`). With the branch in:
        * **piece 2 gets 4 envelope faces and cPLQ's Step 2 returns ZERO conjugate cells for it**
          -- the new envelope is computed and then discarded, so the wiring buys nothing today;
        * **piece 1, untouched by any of this, is WRONG on its own through cPLQ's Step 2**: 6
          cells, `f*(0,0) = 0.28647` where the truth over its own triangle is `0`,
          `f*(0.5,1) = 1.00464` where it is `1`, and NOT COVERED at `(-1,0.5)` and `(-2,-2)`.
      So the crash was masking a **silent wrong answer**, and landing the wiring alone trades a
      loud refusal for it. That is why it is not committed.
      **The one measurement that tells you where to look.** That same triangle conjugated ON ITS
      OWN via `QuaPol.conj` is exact at 7 of 7 probe points -- because a single bounded triangle
      takes the NUMERIC route (`conjBoundedPolygon`), not cPLQ. **The numeric Step 2 is right on
      this input and the vendored symbolic one is not.** `assertStep3MatchesPieces` correctly does
      NOT fire: Step 3 agrees with Step 2: the fault is inside Step 2.
      **THE UNDERLYING DEFECT IS FOUND, and it is in STEP 1, not Step 2.** `plq_1p`'s `nCE == 2`
      branch applies [COAP]'s single-quadratic form to the WHOLE triangle. That form is a valid
      convex MINORANT but A.4 shows it is tight only over a sub-region -- and this branch never
      tests. Measured on `conv{(2.5,1.5),(2,0),(0,0)}` carrying `x*y`, it returns an envelope whose
      minimum over the triangle is **-0.2835** at `(1,0)`; on that triangle `x >= 0` and `y >= 0`,
      so `x*y >= 0`, the affine minorant `0` is admissible, and the TRUE envelope is `>= 0`
      everywhere. A too-small envelope gives a too-large conjugate: that is the `0.28647`.
      `convEnvCPLQ` on the same triangle returns 2 faces -- it does apply the split -- with
      minimum `0`.
      **And routing Step 1 through `convEnvCPLQ` does NOT fix it -- tried, measured, reverted.**
      `conjugateFunction`'s `nCE == 2` branch reads its envelope with `coeffs(...)` and matches
      monomials; `convEnvCPLQ`'s A.4/A.5 faces are RATIONAL, so it raises
      `symbolic:coeffs:NotAPolynomial`. **cPLQ's Step 2 has no rational-envelope branch at all.**
      **So the split belongs in the DOMAIN.** The route that already works for rational faces is
      `conjCPLQ`'s `conjEnvelopeViaCPLQ`, which hands each rational face to cPLQ as its own PIECE
      through `ratPolToPlq`. Do the same here: have `plq_1p.triangulate` split a 2- or
      3-convex-edge triangle into the A.4/A.5 sub-triangles and emit each as a piece, recursing
      while `splitTwoConvexEdges` still reports `needsSplit`. Every sub-piece is then a triangle on
      which cPLQ's own closed form IS tight, Step 2 is untouched, and every envelope stays
      polynomial.
      **The cost, which is why it was not done unattended:** `splitTwoConvexEdges`,
      `splitThreeConvex` and their helpers are file-local to `convEnvCPLQ.m`, so exposing them
      means moving a connected web of functions out of a well-tested file -- and `triangulate`
      feeds every Case C result. Design change plus a full re-verification, not a fix.

- [x] **The PARALLELOGRAM's `emptyResult` — TWO defects found and FIXED 2026-08-16**, taking its
      worst piece from **6 of 10** probe points wrong or uncovered to **2**, against a brute-force
      sup. Both are of this codebase's recurring kind, and both are general fixes rather than
      special cases.
      1. **`region.simplifyUnboundedRegion` declared any region with no finite VERTEX empty** --
         a half-plane, a slab, the whole plane. A half-plane is exactly what a TANGENT vertex
         produces: the cone there is built from the two edges meeting at it, and when those are
         tangent (an arc and its chord touching, how a curvilinear piece ends) both half-planes
         are the SAME one. Now refuted by a WITNESS. `regionTest/aHalfPlaneIsNotEmpty`.
      2. **The edge list, in bug 1's other form.** 3 vertices, 3 genuine edges, plus a conic
         touching one vertex: 4 constraints for 3 vertices, so the count called a BOUNDED region
         unbounded and it was built one edge cell short. `conjugateOfPiecePoly` now derives the
         list for any bounded piece the count mislabels, not just for a lens.
      `functionNDomainTest/aBoundedPieceWithATangentVertexConjugatesOntoTheWholePlane` pins it.
      **How it was found, worth reusing:** `f**` of a bounded domain is finite exactly on that
      domain and is a MAX, so EVERY per-piece conjugate must be finite there. Evaluating all 12
      groups at six interior points named the three bad ones in one cheap run.

- [ ] **The parallelogram's LAST 2 of 10 — `getInterior` on a SINGULAR quadratic.** The chord's
      edge cell and the arc's claim the SAME region and the chord's is checked first and is wrong
      there (0.0176 and -0.0138 against 0.0333 and 0.0363; the arc's cell is right at both).
      `f` is a singular convex quadratic on that piece -- constant along one whole edge, `∇f = 0`
      at two of the three vertices -- so `functionNDomain.getInterior`, which eliminates `s`
      between `x = ∂1f` and `y = ∂2f`, returns the gradient map's image LINE instead of a curve
      separating the two cells.
      **Do NOT attack the `isQuad` chord rewrite for this.** Both alternatives were measured
      2026-08-16 and are recorded in `DECISIONS.md`: chording the vertices the conic actually
      touches makes the piece WORSE (2 wrong of 10 -> 3), and skipping the rewrite changes nothing.
      Start from `functionNDomainTest.parallelogramPiece9`, which reproduces it in about a minute.

### Bugs, in the order they should be taken (2026-08-15)

- [x] **BUG 1 -- FIXED 2026-08-15**, by the edge-list refactor mapped below. Everything from
      "Original framing" down is HISTORY, kept because the four failed attempts are recorded
      against it; the live write-up is the newest `DECISIONS.md` entry.
      **What it was.** Both slot conventions identify an edge by a VERTEX INDEX, and a LENS's two
      edges join the SAME pair of vertices -- so neither convention can name them apart, and no
      reassignment of slots to constraints ever could. That is why all four earlier attempts
      failed, and it is a stronger statement than "the count is the wrong test".
      **What it took.** `functionNDomain.edgeIndexList` derives `eIdx(j)` -- the constraint
      bounding the edge from vertex `j` to vertex `j+1` -- from the geometry (both endpoints on
      the constraint AND its own curve between them inside the region; the second half is
      necessary, since the lens's redundant conic passes through both vertices too).
      `region.getNormalConeVertexQ` takes the list as an optional argument;
      `region.getNormalConeEdgeQE` replaces `getNormalConeEdgeQ` and `Q3`, which turned out to be
      the same routine under two slot conventions; `getSubdiffVertexT2` and `T2Q` are identical on
      these inputs. **The refactor was smaller than this item predicted: three of the "four
      routines that move together" were two routines wearing four names, and both loops are in
      `conjugateOfPiecePoly` itself, so the list needs no field and no signature change.**
      **Scope.** Entered only when two constraints that each bound a genuine edge still share an
      edge number after `spreadCollidingEdges` (`edgesStillCollide`) -- the lens signature. No
      constraint is dropped, so the two unsound drops stay ruled out.
      **Measured.** Both half-lenses conjugate to 3 cells (2 vertex cones + the arc; the chord's
      cell is a ray and drops out), exact against a brute-force sup at 12 points, 0 wrong -- where
      the old code was `+inf` at `(0,0)` and `(-1,0.5)` and `0` at `(2,-1)` for a truth of `1/2`.
      Pinned by `functionNDomainTest/twoEdgesOnOneVertexPairAreBothKept` and
      `halfLensConjugateIsFiniteEverywhereAndExact`, ~10 s where reaching the same piece through
      `biconjugateTest` takes 10-40 minutes.
      Original framing follows.
      **TWO defects, and the description this item used to carry was wrong about the
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

      **THE INDEXING CONTRACT, mapped 2026-08-15 -- this is what the refactor has to preserve.**
      Two conventions run through four routines at once, chosen by the COUNT
      `size(d.ineqs,2) == d.nv`:
        * BOUNDED layout: `nv` slots, edge `j` at `ineqs(j)`, `endNv = nv`, cones from
          `getNormalConeVertexQ` / `getNormalConeEdgeQ3` and `getSubdiffVertexT2`.
        * UNBOUNDED layout: `nv+1` slots with slot 1 reserved for the ray, edge `j` at
          `ineqs(j+1)`, `endNv = nv-1`, cones from `getNormalConeEdgeQ` and `getSubdiffVertexT2Q`.
        * Confirmed in the source: `getNormalConeEdgeQ` loops `j = 1:nv-1` and reads
          `slopeIneq(j+1, vertex j)`, so its output ROW `j` is the edge at slot `j+1`.
          `getNormalConeVertexQ` uses the same "vertex `j` lies on constraints `j` and `j+1`" rule
          (now wrapped cyclically, so it works in the bounded layout too).
      The loop variable `j` indexes `NCE`'s rows, `subdE`'s rows and `d.ineqs` SIMULTANEOUSLY, so
      the count cannot simply be replaced -- all four move together.

      **THE LENS NEEDS THE BOUNDED LAYOUT AND CANNOT GET IT.** It has `nv = 2` and 2 genuine edges,
      so the bounded layout is right; but it arrives with 5 constraints, so the count picks the
      unbounded one, `endNv = 1`, and only one edge cell is built.
      Two ways to give it the bounded layout, and the trap in each:
        1. DROP the three non-edge constraints. **Unsound** -- a constraint active at exactly one
           vertex of a convex region can still be essential; measured, see `DECISIONS.md`. It would
           be sound if guarded by a REDUNDANCY PROOF, but `redundantSubset` is an LP and cannot
           certify redundancy in the presence of a conic. A closed-form "maximise this linear form
           over the region's own boundary" test would -- the machinery exists in
           `verifyMaxIsExactSymbolically`, for QuaPar faces rather than `region`s.
        2. PARK them above slot `nv` and dispatch on boundedness instead of the count.
           `wasBounded` is now computed soundly (before `removeInfV`, which deletes the evidence).
           The trap: the `f.isQuad` chord rewrite iterates ALL slots, so a parked QUADRATIC
           constraint gets a chord derived from vertices it does not join -- the harm rule (2) in
           `conjugateOfPiecePoly` already records. Scope any parking to `~obj(i).f.isQuad`.
      Start from the hand-built lens probe -- it needs no pipeline and runs in seconds, and with
      the three fixes already in it produces the correct 3 cells when handed exactly its 2 edges.

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

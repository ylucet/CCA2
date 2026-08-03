# Session Handoff

_Last updated: 2026-08-03 07:15 PDT_

## What happened this session

**The triangle-case enumeration was proved and then used to drive the code.** `proveStageA–D.m`
discharge symbolically that the case list is exhaustive, what envelope each case yields, and — the
payoff — that **every dual boundary is a line, a line pair or a parabola, never an ellipse or
hyperbola** (1140 boundaries classified, zero exceptions). The parabolic edge originates in exactly
**one** case: indefinite with **one convex edge**.

**Two closed-form wins on `main`.** Step 2 needs no rational-piece conjugate at all (conjugating the
ORIGINAL quadratic over the sub-triangle is exact), and bounded polygons now take the closed-form
path (Case B2), which also fixed Case C's biconjugate (`f** = +inf` → `f** = f`).

**Arc-vs-arc is built but not finished.** `clipArcByConic` is done and validated standalone;
clip and split now produce a **valid partition**. Three tests still fail, all in `assemblePieces`,
and the cause is localised to one cell and one edge (below).

## Where things stand

- Branch: `overnight/2026-08-02` @ `7d145d2` — "Localise the last defect to one cell and one edge"
- 8 commits ahead of `main`. **`main` is fully pushed** (18 commits went to `origin/main` this
  session).
- Pushed: **no** for this branch — it has never been pushed, so it has no tracking branch yet.
  (`origin` exists; `git push -u origin overnight/2026-08-02` would create it.)
- `maxQuaParTest` **16 pass / 3 fail** — the three arc-vs-arc pins. Elsewhere: two deliberate reds
  in `unboundedFaceTest`, one by-design `biconjugateTest` failure. Nothing else red.
- Suite buckets: `.claude/suite.sh --fast` (189 tests, **100 s**), `--normal` (266 s), `--slow`
  (the symbolic suites, hours). Slow == symbolic, exactly.

## Next steps

1. **THE ONE REMAINING DEFECT — read `TODO.md` first, it has the full evidence.**
   Piece 5 (src `[2 2]`) emits a **ray** where its boundary should **terminate**. Sampling across
   its apex `(-2.03125, 2.03125)` shows three cells locally — `(2,2) | (6,1) | (6,2)` — with g2's
   arc separating the last two, so piece 5's neighbour is a **bounded** sliver (area 0.008). A
   bounded neighbour means a finite SEGMENT, and `matchHalfEdges` pairs rays with rays and segments
   with segments — so a ray facing a segment can never match. That is the reported symptom.
   **Lead:** the restriction of g2's arc conic to that ray has `A = 1.7e-18`, numerically
   degenerate, so `conicRootsAlong` takes its linear branch and returns ONE root. Check whether the
   genuine second crossing is lost there, or sits on a different boundary element.
2. **Do not re-try the retired hypotheses** listed at the bottom of `TODO.md` — six things already
   refuted by measurement, including three claims of mine that had to be retracted. The sharpest:
   a "zero intersection" reading from a 0.0625-spaced grid was a resolution artifact for cells of
   area 0.008, and acting on it would have deleted legitimate cells.
3. Decide whether to merge `overnight/2026-08-02` into `main` and push it.
4. Older open items: the Step-3 over-claim on unbounded assemblies, and the curved-envelope
   derivation for a nonconvex q over an unbounded piece. Both have red tests pinning them.

## Relevant files

- `TODO.md` — the live defect analysis and the retired-hypotheses list. Start here.
- `MORNING.md` — the overnight run's account.
- `maxQuaPar.m` — `clipPolyByConic`, `conicRootsAlong` (just fixed), `clipByFace`, `splitCell`,
  `assemblePieces`/`matchHalfEdges`, `partitionReport` (sound only for polyhedral pieces).
- `clipArcByConic.m` + `clipArcByConicTest.m` — the arc-vs-arc primitive, 7 tests, one per
  geometric configuration.
- `proveStageA.m` … `proveStageD.m`, `sweepCaseEnumeration.m` — the enumeration and its proof.
- `conjCPLQ.m` — `conjFaceOrOriginal`, `conjBoundedPolygon` (Case B2).

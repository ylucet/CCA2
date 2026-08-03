# Session Handoff

_Last updated: 2026-08-02 23:40 PDT_

## What happened this session

**The triangle-case enumeration was made exhaustive, proved, and then used to drive the code.**
Stages A–D (`proveStageA/B/C/D.m`) discharge symbolically that the case list is complete, what
envelope each case yields, and — the payoff — that **every dual boundary is a line, a line pair or
a parabola, never an ellipse or hyperbola** (1140 boundaries classified, zero exceptions). The
parabolic edge originates in exactly **one** case: indefinite with **one convex edge**.

**Two closed-form wins.** Step 2 no longer needs a rational-piece conjugate at all — conjugating
the ORIGINAL quadratic over the sub-triangle is exact, which made the 2-convex-edge case fully
numeric. And bounded polygons now take the closed-form path (Case B2), which also fixed Case C's
biconjugate (`f** = +inf` → `f** = f`).

**Arc-vs-arc is half-built.** `clipArcByConic` is done and validated standalone (7/7, one test per
geometric configuration) and wired into `clipByFace`; clip and split now produce a **valid
partition**. The remaining failure is confined to `assemblePieces`.

## Where things stand

- Branch: `overnight/2026-08-02` @ `57d2206` — "Overnight run: clip/split now partition correctly"
- `main` @ `88d482b`, **18 commits ahead of `origin/main`, none pushed.** The overnight branch is
  4 further commits on top of `main`.
- Pushed: **no** — nothing was pushed this session (the overnight run's guardrail is local-only).
- **Suite split into buckets:** `.claude/suite.sh --fast` (16 suites, 189 tests, **100 s**),
  `--normal` (3 suites, 266 s), `--slow` (7 symbolic suites, hours). Slow == symbolic, exactly.
- `maxQuaParTest` 16 pass / **3 fail** — the three arc-vs-arc pins. Elsewhere: two deliberate reds
  in `unboundedFaceTest` (curved envelope over an unbounded piece) and one by-design
  `biconjugateTest` failure. Nothing else is red.

## Next steps

1. **`assemblePieces`: a boundary ray with no matching neighbour.** All three arc-vs-arc fixtures
   now fail here and nowhere else. The arrangement contains genuine thin cells (area ~0.008 beside
   cells of ~10) whose corners are the orphan apexes; those apexes differ by ~0.044, far above
   `matchHalfEdges`' 1e-3 tolerance, so they are genuinely different rays — **not** the
   near-degenerate cluster case `checkOrphanHalfEdges` already handles for segments.
2. **Read `TODO.md`'s "Retired hypotheses" before touching this.** Five things are already refuted
   by measurement, including two of my own retracted claims. The sharpest lesson: a "zero
   intersection" reading from a 0.0625-spaced grid was a resolution artifact for cells of area
   0.008, and acting on it would have deleted legitimate cells.
3. **Decide whether to merge `overnight/2026-08-02` into `main`**, then whether to push the 18+4
   commits.
4. Still open from before: the Step-3 over-claim on unbounded assemblies, and the curved-envelope
   derivation for a nonconvex q over an unbounded piece (both have red tests pinning them).

## Relevant files

- `MORNING.md`, `TODO.md` — the overnight run's full account and the next concrete step.
- `proveStageA.m` … `proveStageD.m`, `sweepCaseEnumeration.m` — the enumeration and its proof.
- `clipArcByConic.m` + `clipArcByConicTest.m` — the arc-vs-arc primitive (parabolas only, so the
  cut is always one representable `Ec` arc).
- `maxQuaPar.m` — `clipPolyByConic`, `clipByFace`, `splitCell`, `assemblePieces`,
  `partitionReport` (sound only for polyhedral pieces — it says so).
- `conjCPLQ.m` — `conjFaceOrOriginal` (the rational sidestep), `conjBoundedPolygon` (Case B2).
- `.claude/suite.sh` — the fast/normal/slow buckets.

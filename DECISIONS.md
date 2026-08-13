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

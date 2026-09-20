# Morning report — 2026-09-19 overnight run

Branch: overnight/2026-09-19 (7 commits, all local; nothing pushed)

## What changed

- **Compaction: 321 faces → 160, 41 edges instead of 69** (item 2), on the convex fixture, still
  reproducing f exactly — 0.0e0 at 300 points, no domain mismatch. `compactCells` drops a row that
  the cell's other rows already imply, decided by `cellInterior` so CURVED and BRANCH rows are
  candidates too, and alternates with `mergeAdjacentCells` to a fixpoint. The older linear-only
  `dropRedundantRows` is now unreferenced except in a comment — delete it after review.

- **`cellInterior` decides a TWO-PAIR cell instead of keeping it** (item 3). The case was recorded as
  never having occurred; the probe that asked built one immediately — four concurrent lines from the
  pairs `x² = 2y²` and `x² = 3y²` leave 8 of 16 sign patterns undecidable in either quadratic field
  alone, and all eight are EMPTY cells that were being kept. `ratQ.qs2*` is the compositum
  Q(t₁, t₂): exact zero test, then interval bisection on exact rational midpoints. Two pairs of the
  SAME field are refused by name rather than answered, because the zero test is invalid there.

- **Bridge faces now enter the arrangement for COVERAGE** (item 1, the half that needs no decision).
  Their outer bounds join the curve table as curves, and membership is asked ONCE per cell at an
  interior point — sound because the predicate is constant on a cell (50000 points, 40 cells, none
  mixed). The blanket refusal is gone: a curved-tie envelope is now refused by
  `envelopeMesh:bridgeCompetes` only where a bridge's VALUE competes with a quadratic face.

- **The bridge item was split by measurement, not by construction** (item 1). Two constructions over
  sign conditions were already refuted; asked of the arrangement instead, the REGION half is settled
  and the VALUE half turns out to need a curve the conic level cannot hold. DECISIONS.md 2026-09-19.

## What is broken

Nothing. Both daily gates ran once, at the start of the run and mid-run:

- **slow**: 98 passed, 0 failed, 5 suites.
- **verylong**: 37 passed, 0 failed, 8 incomplete — all of them the documented `assumeFail`
  quarantines (the G11 biconjugate fixture and the legacy `plq_1piece` pins), not new.
- **fast**: 688 passed, 0 failed over 48 suites at the last run.

Two reds appeared and were fixed inside the run, both mine and both instructive: `compactCells`
stripped THIN faces (five reds — a point domain has no interior, so every row's violation is
vacuously empty), and each new guard had to be added to `guardExerciseTest`'s inventory.

## Needs a decision

**Does the SUBDIVISION axis grow?** This is the only thing left on the list, and it is yours.

A bridge face carries the envelope's value — measured strictly above the best quadratic face at 22
of 60 points, by up to 1.2e-1 — so where it competes the mesh must draw `FP(x, y, q(x,y)) = 0`, a
rational curve of degree well above two. Measured on the c = -2 fixture: 2 of 40 cells, and in both
the competitors are faces of DIFFERENT cell pairs, so it is a genuine max between pairs rather than
a partition inside one. A conic subdivision cannot hold that edge.

1. **Grow the subdivision axis** to carry a higher-degree curve — the `AlgAlg` level set aside on
   2026-09-13 because only the function axis seemed to need it. Biggest change; gives a real mesh, a
   boundary, and the paper's figures. **This is what I would pick**: every other route leaves
   `biconj` returning something that cannot draw its own domain.
2. **Return the face list** for a curved tie — what `envelopeEval` already answers from — and keep
   the mesh for straight and tangential ties. Cheap and correct, but no boundary.
3. **Mesh only where no bridge competes**, refusing otherwise. That is what is now implemented, and
   on this fixture it still refuses, so it buys little on its own.

Smaller, and not blocking: `dropRedundantRows` is now dead code, and `envelopeMesh:bridgeCompetes`
is on the guard inventory as `todo` — pinning it needs a normal-bucket test, since building a bridge
face costs the whole assembly.

## Where I stopped

All three items are off the list: 2 and 3 are done and committed with tests; 1 is split, its region
half implemented, and its value half is the decision above. Nothing else under **Next up** is
actionable without you, which is the stop condition.

The run is on `overnight/2026-09-19`, seven commits, nothing pushed. `git log main..overnight/2026-09-19`
is the whole night; `git switch main` abandons it.

# Session Handoff

_2026-09-19_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main`, well ahead of `origin/main`, tree clean
- Tests: fast bucket over 46 suites · `envelopeAssembleTest` 8/0 (normal) ·
  slow 98/0 · verylong 37/0 (09-06/07, not re-run since)
- Known reds: none

## Where the envelope stands
The corner faces' REGIONS are exact rows now (`cornerWedgeRows`, 6400/6400 against
the KKT predicate), and the arrangement can take them: `cellInterior` decides a
cell with branch rows in three stages — exact rational tests, a margin LP in
doubles, then an exact Farkas certificate for the ties. Measured on the c = -2
fixture: 486 live cells where ignoring the branch rows gives 31104, and NOTHING
left undecided.

## Next
1. **The corner faces' VALUES.** A corner's value is affine with quadratic-
   irrational coefficients and the mesh emits a QuaCon of quadratic faces, so
   `envelopeMesh` still refuses an envelope containing them — now purely an
   OUTPUT TYPE question, the geometry being done.
2. **The BRIDGE face's region** stays a predicate (its cut is refuted twice,
   DECISIONS 2026-09-17).
3. Compaction: 321 faces where a handful is right.

## Files
- `cellInterior.m` — the emptiness test, and `cellInteriorTest.m`'s six-sector count
- `cornerWedgeRows.m` — a corner's wedge as conjugate-pair rows
- `ratQ.qsRed/qsAdd/qsMul/qsSign` — the Q(t) arithmetic, in one place
- `.claude/branchRowPruning.m` — the measurement behind all of this

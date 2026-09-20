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
The CORNER faces are in the mesh — carried for their REGION, not their value:
measured, a corner's value never exceeds the best quadratic face owning the same
point (largest excess exactly 0.0e0), so what a corner decides is COVERAGE, and a
pair with no covering face is +infinity. A BRIDGE face is the opposite and is now
the only blocker.

The corner faces' REGIONS are exact rows now (`cornerWedgeRows`, 6400/6400 against
the KKT predicate), and the arrangement can take them: `cellInterior` decides a
cell with branch rows in three stages — exact rational tests, a margin LP in
doubles, then an exact Farkas certificate for the ties. Measured on the c = -2
fixture: 486 live cells where ignoring the branch rows gives 31104, and NOTHING
left undecided.

## Next
1. **The BRIDGE face's region** — now the ONLY thing between a curved tie and a
   meshed `biconj`. Its cut is refuted twice (DECISIONS 2026-09-17); either solve
   it or carry a bridge as a predicate face in a type that allows one. Its VALUE
   already has a home in `AlgCon`'s `FP`/`FB`, and it does need one: measured, a
   bridge is strictly above the best quadratic face at 22 of 60 points.
2. Compaction: 321 faces where a handful is right.
3. `cellInterior` skips a certificate whose branch rows name two different pairs
   and keeps the cell; it has not occurred.

## Files
- `cellInterior.m` — the emptiness test, and `cellInteriorTest.m`'s six-sector count
- `cornerWedgeRows.m` — a corner's wedge as conjugate-pair rows
- `ratQ.qsRed/qsAdd/qsMul/qsSign` — the Q(t) arithmetic, in one place
- `.claude/branchRowPruning.m` — the measurement behind all of this

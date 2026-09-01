# Session Handoff

_2026-08-31_

## Blocked
- 80% fast+normal coverage — needs YOUR call on two things: `plq_1piece.m`
  (migrate T6/G14 fixtures to `plq_1p`, or exclude from the metric — see
  DECISIONS 2026-08-31) and whether print/plot lines count (21% of the gap).

## State
- Branch `main` @ `e3a0671` — "fix: limitDirectional crashed on a direction
  tangent to the denominator's zero set"
- Pushed: pending
- Tests (2026-08-31): fast 363/0 (30 suites, ~175 s) · normal 13/0 (369 s) ·
  slow 98/0 (measured before today's deletions; re-run before a tag)
- Known reds: none. G11 still quarantined via `assumeFail`.
- Coverage fast+normal, production: 59.3% → **68.5%** (74.9% excl
  `plq_1piece.m`); +615 lines to 80%.

## Next
1. `plq_1p` conjugate branches — 268 uncovered lines, the largest single
   remaining block and the hardest.
2. `fanUnboundedFace.cutNonPointed` (19) + `boundaryOrder` (16);
   `quaPolToPlq.faceDomainFromHalfPlanes` (14). All have clear contracts.
3. Re-run `--verylong` before any tag: it uses `plqCheck`, extended today.

## Files
- `plqCheck.m` — region fixtures + definition checks, memoized
- `meshPredicateTest` `frameAndFanTest` `regionUnitTest` `regionMinusTest`
  `symbolicFunctionUnitTest` `QuaParCPLQUnitTest` `biconjCPLQUnitTest` — new
- `.claude/delete-dead-functions.py` — reproduces the dead-code sweep

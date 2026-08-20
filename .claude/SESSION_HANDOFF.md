# Session Handoff

_2026-08-20_

## Blocked

- **A.4 cevian split in `plq_1p`** — needs `maximumP` to max two RATIONAL
  conjugates; `region.normalize1` raises `NotAPolynomial` on one. Not external:
  the fix is stated at `biconjugateTest.m:250` (clear denominators where both
  are provably nonzero on the cell). Three attempts; reverted code kept at
  `.claude/a4split_attempt.m.txt`. See DECISIONS 2026-08-19 (night, later).
- **T6, delete `plq_1piece`** — a migration, not a deletion; behind A.4.
- **Phase C1 per-term cost target** — EXTERNAL, needs Yves. Box terms are now
  0.01 s, so the old 40-60 s figure is stale; C2/C3 wait on the number.

## State

- Branch `main` @ `2802f40` — "WIP: slopeAtVertex must return a REAL slope"
- Pushed: yes — 18 commits, `801ee1f..bbb4680`, 2026-08-20
- Tests (2026-08-20): fast 217/0 (93 s) · normal 11/0 (215 s) · slow 88/0
  (601 s) · verylong one red
- Known reds: `testMaxMultiRegion/testPCE2` — `plq_1p`'s 2-convex-edge branch
  returns a convex MINORANT, not the envelope. On `{(0,0),(1,0),(2,1)}`,
  `f = x*y`, it dips to -0.0429 where `f >= 0`. Real defect, pre-dates the
  session.
- **`2802f40` is UNVERIFIED** — `regionTest` was interrupted before reporting.

## Next

1. Verify `2802f40`: `bash .claude/suite.sh regionTest` (~45 s). Third guard on
   `slopeAtVertex`; each earlier one moved the crash two lines down the same
   caller (pole → `double()` → `atan2`). A fourth move ⇒ umbrella §5 rung 3,
   state and test the contract.
2. Teach the cross-face max to handle a rational pair — unblocks A.4 → T6 → T3.
3. T3-T5 sym-free port, once one per-piece class remains. `exactQ` (T2) is green.

## Files

- `DECISIONS.md` — read the A.4 and T6 headings first
- `.claude/a4split_attempt.m.txt` — reverted split; working cevian geometry
- `region.m` — `slopeAtVertex`, `probeOnConstraint`, `rootsIn`, `isFeasible`
- `plqCheck.m` / `plqStage.m` — definition verifiers, stage cache

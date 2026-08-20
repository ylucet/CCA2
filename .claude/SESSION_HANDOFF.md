# Session Handoff

_Last updated: 2026-08-20_

## What happened this session

**The symbolic-removal programme finished for `region.m` and then hit a real wall.** Live `solve()`
went **16 → 3** (two of the three a deliberate fallback, measured never hit; one kept for conic
pairs). `getNormalConeVertexQ`'s specification was established and made a test, A3 was answered
(merge never over-claims), Step 0 and Phase B landed, and `SUPPORT_MATRIX.md` §0.0.1 was re-derived.

**The test suites were the bigger change.** `testcPLQ` and `testMaxMultiRegion` held 32 tests and
**zero assertions** — they ran the pipeline, printed, and returned. They now verify against
definitions (`plqCheck`), split by stage with caching (`plqStage`), and live in a new `--verylong`
bucket. `suite.sh` gained `-j` (2 h → 73 min). That work immediately paid: it exposed three real
defects that had been invisible.

**Then the sym-free port (T1–T8) started and stopped at a genuine blocker** — see below.

## READ THIS FIRST — what is actually blocked, and why

1. **`testMaxMultiRegion/testPCE2` fails on `main`, and did before this session.** It is the
   executable statement of a real defect: `plq_1p`'s 2-convex-edge branch returns [COAP] A.4's
   single quadratic over the WHOLE triangle, which is a valid convex MINORANT but not the
   envelope. Measured on `{(0,0),(1,0),(2,1)}`, `f = x*y`: dips to −0.0429 where `f ≥ 0`.
2. **The fix is written and does not work yet — and the reason is NOT in `plq_1p`.** Three
   attempts (`DECISIONS.md`, and the code is kept at `.claude/a4split_attempt.m.txt`). The cevian
   geometry is correct, both envelope faces build, both conjugate, both reach the max
   (`conjfia = [1 6 12]`). What loses the answer is the max ACROSS them: face 1's conjugate is
   quadratic, face 2's is RATIONAL, and `biconjugateTest.m:246` already records that `splitmax3`
   hands `f1 − f2` to `region()`, whose `normalize1` raises `NotAPolynomial` on a rational.
   **Adding the split before that is fixed trades an overshoot for an UNDERSHOOT, which is worse.**
3. **T6 (delete `plq_1piece`, −75 symbolic calls) is a MIGRATION, not a deletion.** Swapping the
   fixtures regresses three tests. Blocked behind item 1.

## Where things stand

- Branch: `main` @ `2802f40` — "WIP: slopeAtVertex must return a REAL slope, not merely a convertible one"
- Pushed: **yes** — 18 commits, `801ee1f..bbb4680`, on 2026-08-20
- **fast 217 / 0** (93 s), **normal 11 / 0** (215 s), **slow 88 / 0** (601 s, five suites)
- `--verylong` (`testcPLQ`, `testMaxMultiRegion`): **one known red**, `testPCE2`, cause above
- **`2802f40` IS UNVERIFIED** — `regionTest` was interrupted before reporting. Run it first.

## Next steps

1. **Verify `2802f40`**: `regionTest` (~45 s), then
   `testcPLQ/rectBiconjugateIsAConvexUnderestimator` (~50 min). This is the third guard on
   `slopeAtVertex`; each previous one moved the crash two lines down the same caller (pole →
   `double()` → `atan2`). **If it moves a fourth time, stop patching symptoms** — state and test
   `slopeAtVertex`'s contract directly.
2. **The unblocking item, now well-defined:** teach the cross-face max to handle a RATIONAL pair —
   clear denominators where both are provably nonzero on the cell, the fix `biconjugateTest.m:250`
   itself proposes. That unblocks A.4 → which unblocks T6 → which unblocks T3.
3. **T3–T5, the sym-free port**, once one per-piece class remains. `exactQ` (T2) is built and
   green; `symbolicFunction`'s payload becomes a coefficient vector in it, then `region`'s
   predicates, then the 96 `isAlways` calls.
4. Phase C1 needs **one decision from the user**: the per-term cost target. Box terms are now
   0.01 s, so the old 40–60 s figure is stale and C2/C3 wait on the number.

## Relevant files

- `DECISIONS.md` — the six newest entries are this session's; read the A.4 and T6 ones before
  touching either. One entry carries a **correction** to an overstated claim made earlier the
  same day.
- `.claude/a4split_attempt.m.txt` — the reverted A.4 split: working cevian geometry
  (`a4Split`, `twoEdgeQuad`) plus a real fix for `conjugateFunction` dispatching on the PIECE's
  domain instead of the FACE's. Both reusable once the rational max lands.
- `plqCheck.m` / `plqStage.m` — the definition-based verifiers and the stage cache
- `exactQ.m` / `exactQTest.m` — T2, exact `a + b√d` arithmetic, 19 tests
- `.claude/suite.sh` — four buckets and `-j`; the throttle is `xargs -P` for a measured reason
- `region.m` — `slopeAtVertex` (item 1 above), `probeOnConstraint`, `rootsIn`, `isFeasible`

# Session Handoff

_2026-09-05_

## Blocked
Nothing.

## State
- Branch `main` @ `1a4185c` — "docs: record the CRT refutation in DECISIONS.md"
- Pushed: yes
- Tests (2026-09-05): fast 488/0 (35 suites, ~190 s) · normal 13/0 (324 s) ·
  slow 98/0 (963 s, run once overnight)
- Known reds: none. Coverage, exact pipeline, fast bucket: 88.8%.

## Next
1. **UNANSWERED QUESTION, answer it first.** Which algorithm for the exact
   piecewise convexity test, and does it follow Singh & Lucet (SIOPT 2021) or
   Singh's 2019 thesis Algorithm 3, which `QuaPol.isEdgeConvex` cites? NEITHER
   IS IN `reference/` — I have not read them. Do not restate the derivation
   sketched in TODO.md as if it were theirs.
2. `biconjQ` exact convexity test — unblocks all 8 fixtures the legacy envelope
   answers and this one refuses. See TODO.md 2026-09-05 (biconj tasks), which
   also names two LATENT smells in the existing `isConvex`/`isEdgeConvex`.
3. Then: unbounded non-convex envelope; then edge-convex + `AlgAlg`.

## Files
- `TODO.md` 2026-09-05 — biconj tasks, overflow options, coverage
- `conjQ.m` `biconjQ.m` — the exact operators
- `ratQ.m` — arithmetic + degree-4 sign kernel (Sturm, signAt)
- `.claude/*-probe.m` `*-sweep.m` `*-diff.m` — reproduce every number quoted

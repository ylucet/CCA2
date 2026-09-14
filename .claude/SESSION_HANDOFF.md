# Session Handoff

_2026-09-14_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — delete / track /
  note? Your call. MORNING.md.
- `maxQuaPar:notDegenerate` reachability — an in-contract pair needs the symbolic
  pipeline (hours). DECISIONS 2026-09-06.

## State
- Branch `main` @ `e110130` — "chore: update session handoff"
- Pushed: yes — `origin/main` @ `e110130`
- Tests: fast 630/0 (09-13) · normal 13/0 (09-11) · slow 98/0 · verylong 37/0
  (09-06/07, not re-run since)
- Known reds: none

## Next
1. Promote the BRIDGE faces: `.claude/bridgeImplicitProbe.m` and
   `bridgeRegionBuild/Sides.m` are already exact — fill `AlgCon.FP`/`.FB` from them.
2. Fold the per-piece faces into one `AlgCon`; route `biconjQ`'s
   `convexPieceCoupling` through it so `biconj(QuaPol)` stops refusing.
3. The paper's two classification tables — which input gives which conic in
   `conj`, which return type in `biconj` — one witness fixture per row.

## Files
- `envelopeCellPairs.m`, `envelopeCellRegions.m` — assembly steps 1-2, exact
- `AlgCon.m`, `algBranch.m` — implicit-face type and its branch evaluator
- `.claude/bridgeImplicitProbe.m` — the elimination to promote next

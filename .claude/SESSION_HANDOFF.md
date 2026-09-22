# Session Handoff

_2026-09-22_

## Blocked
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — your call.
- `maxQuaPar:notDegenerate` reachability — needs an in-contract symbolic input
  (hours). DECISIONS 2026-09-06.

## State
- Branch `overnight/2026-09-20` @ `d32bd18` — "feat: algCrvIncidence -- an AlgCrv
  gets its edge and face adjacency"
- Pushed: no upstream; also UNMERGED to `main` (10 commits ahead)
- Tests: fast 707/0 (09-22) · normal 8/0 (09-20) · slow 98/0 and verylong 37/0
  with 8 documented `assumeFail` quarantines (09-19, not re-run since)
- Known reds: none

## Next
1. Order vertices ALONG a curve. 31 of the c = -2 object's 56 edges carry more
   than two, so `E` describes only part of them (`algCrvIncidence`,
   `info.oversubdivided`). Needs comparison of algebraic numbers.
2. Move `fN, fD, EcQ, Vname, FC` from QuaCon up to AlgCon. The five-line version
   is REFUTED — 692 -> 363 passing, 326 incomplete. DECISIONS 2026-09-20.
3. Merge the branch into `main`, then run `bash .claude/suite.sh --slow`
   (~20 min; stop the podman machine first, it commits 4 GB).

## Files
- `Crv.m`, `AlgCrv.m` — the fourth subdivision level and the (Alg, Crv) type
- `envelopeMeshAlg.m`, `assembleAlgCrv.m` — the curved-tie mesh `biconjQ` returns
- `algCurveSide.m` — the new edge's side, by Sturm root count
- `algCrvVertices.m`, `algCrvIncidence.m` — vertices, then adjacency

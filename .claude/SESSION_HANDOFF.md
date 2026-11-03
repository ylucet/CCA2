# Session Handoff

_2026-11-02 (overnight run)_

## Blocked
Nothing.

## State
- Branch `overnight/2026-11-01` merged into `main` 2026-11-02 (fast-forward, not pushed) -- see `MORNING.md`
- Tests: fast 814/0, normal 32/0, slow 97/0 + subdivision slow suites (2026-11-02 on this branch)
- Known reds: none

## Next
1. Example 2's default DECIDED (2026-11-02): biconjQ keeps declining its overlay by the piece
   budget; `overlayBudget = Inf` subdivides it (~30 min, verylong test).
2. Tangent faces (most of what Example 2's pair cells keep): exact certificate derived in DECISIONS
   2026-11-01, needs a big-integer square-free part of a degree-20 polynomial.

## Files
- `envelopeSubdivisionAlg.m` -- dropBelow (below the max), minOverlay (tieEmpty, budget)
- `algFaceBelow.m`, `algRowFails.m` -- several upper faces; division filter; row memo
- `cutTriangles.m`, `AlgCrv.m` (table 3, FT, vectorised eval), `assembleAlgCrv.m`, `biconjQ.m`
- `erratumBiconjEx2SubdivisionTest.m` (verylong)

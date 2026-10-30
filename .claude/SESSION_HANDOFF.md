# Session Handoff

_Last updated: 2026-10-20T17:56-07:00_

## What happened this session
The biconjugate list (TODO items 1-7) was worked through: active atom-pair prune, subdividing of one-pair
and several-pair (min of maxes) EnvelopeMins, the poly2 kink subdivided (~3 min, was >90 min), `biconj`
now defaults to 'exact' (Yves's call), suites rebalanced so the normal bucket fits its budget, and a
certified Bernstein test drops empty faces caused by algebraic rows. Every change and its
measurements is in DECISIONS.md 2026-10-16..20.

## Where things stand
- Branch: `main` @ `e09b1e2` -- "feat(subdivision): certified emptiness for algebraic rows; normal bucket within budget"
- Pushed: no -- GitHub credentials missing/invalid on this machine (`gh auth login` needed); ~10 commits ahead of origin
- Tests (2026-10-20): fast 766/0; normal 30/0 in 535 s (budget 600); slow subdivision suites green
  (minSubdivisionTest 3/0, kinkSubdivisionTest 3/0, algCrvSubdivisionTest 3/0, erratumBiconjEx2Test 4/0,
  envelopeSubdivisionSlowTest 1/0). Full `--slow` last run whole on 2026-10-17 (106/0).

## Next steps
- Re-auth (`gh auth login`) and push main.
- Yves: Example 2's subdivision needs a type for triangle faces with degree-4 irrational corners; it stays
  an EnvelopeMin (eval only) until then.
- TODO item 3 remainder: ~1600 empty faces on the 3-convex-edge triangle whose cell the curve
  R = FP(x,y,q) = 0 crosses -- needs that curve's arrangement to prove them empty (values unaffected).
- Optional: planar Minkowski-Weyl in Lean (lean/README "Not formalized").

## Relevant files
- `envelopeSubdivisionAlg.m` -- the subdivider: header, minMaxFaces (several pairs), liveCell, pieceLive, algEmpty
- `biconjQ.m` (subdividePair, hullRows), `biconjAtoms.m` (activeAtomPairs), `envelopeEval.m` (crack retry)
- New helpers: `polySignOn.m`, `cellSamplePoint.m`, `envelopeRationalCorners.m`
- New slow suites: `kinkSubdivisionTest.m`, `minSubdivisionTest.m`, `erratumBiconjEx2Test.m`, `envelopeSubdivisionSlowTest.m`
- `DECISIONS.md` (2026-10-16..20), `TODO.md` (biconjugate list at the top)

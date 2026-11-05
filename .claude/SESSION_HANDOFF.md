# Session Handoff

_2026-11-04_

## Blocked
Nothing.

## State
- `overnight/2026-11-02` merged to main 2026-11-03 (fast-forward, not pushed) -- see `MORNING.md`
- Tests on main (e45e39e): fast 826/0, normal 35/0, slow 122/0 (+2 known-incomplete), verylong Ex2Subdivision 3/0
- 2026-11-03: triangle 140 -> 84 faces, zero-area 12 -> 2 (mirror pair at (1/(2 sqrt5), 1 - 1/(2 sqrt5)));
  algBranch double-root window exact; face row vs a rational face; bridgeMember fixed beside the
  edge so erratumBiconjTest checks every point again (DECISIONS 2026-11-03 (second)-(fourth))
- 2026-11-04: surd face row (faces over Q(sqrt 55) compared directly): Example 2 pairs 29/91/47 -> 20/79/35;
  fast 826/0, normal 35/0, slow 125/0 (+2 known), Ex2Subdivision 3/0
- 2026-11-04 (later): rows against a quartic face's TOP root (pairs 20/79/35 -> 20/62/35) and surd rows
  divided by the region's edge factors (triangle 84 -> 60 -> 56, no holes); overlay default budget
  4000 -> 3000 (Example 2's bound fell to 3736; keeps the default fast). DECISIONS 2026-11-04.
  Tests on main (53b0bb3): fast 826/0, normal 35/0, slow 127/0 (+2 known), Ex2Subdivision 3/0
- Known reds: none

## Next
1. (done) Branch merged to main 2026-11-03.
2. Example 2 faces left: pairs 20/62/35 (winners 15/32/15). Pair 2's open calls: corner touches
   (gaps 1e-6..1e-12), places with no top root, z_k below several roots (DECISIONS 2026-11-04).
3. Triangle at 56 faces; ~20 own no point of 400000 samples (small or zero-area) -- not analysed.
4. Guards still `todo`: bridgeFaceAlg bigZContent/branchOverflow/witnessOverflow,
   envelopeAssemble:irrationalStraightFacet (no fixtures).

## Files
- `algRowFails.m` -- topRowPolys / topFaceRow (top-root rows, no corner blow-up), divideOutSurd
- `algRowFails.m` -- surdFaceRow / surdPoly / oddPartZ (exact BigInteger blow-ups); `algRowFailsSurdTest.m` (slow, new)
- `algRowFails.m` -- face rows (.W), Yun odd part, symbolic square-free past 2^53
- `algFaceBelow.m` -- face-row fallback, `own` rows; `envelopeSubdivisionAlg.m` dropBelow passes them
- `guardExerciseSymbolicTest.m` (slow, new); `addQuaPar.m`, `conjConvexOverPiece.m` guard fixes

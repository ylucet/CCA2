# Session Handoff

_2026-08-28_

## Blocked
- Nothing external. Two-conic extension for item 3 and the `assemblePieces`
  defect (G1/G4/G10) are both open but merely hard, not blocked — see Next.

## State
- Branch `main` @ `2e80373` — "chore: final overnight report update -- stale-entry
  sweep and closing note"
- Pushed: yes — `origin/main`
- Tests (2026-08-28): fast 309/0 · normal 12/0 · slow 94/1 · verylong 29/0/7/1-timeout
- Known reds: `theEMPTYDomainConjugateRoundTripsToMinusInfinity` (slow, pre-existing
  per B3) · `testcPLQ` (verylong, G11's known unresolved timeout)

## Next
1. Two-conic extension of item 3's `maxAffineOverRegion` fix (60% of remaining
   `exactAnotInB` refusals need it) — a harder proof, not a continuation; see
   `DECISIONS.md` 2026-08-27 (item 3 follow-up) for the sketch.
2. `assemblePieces` curved/straight matching (G1/G4/G10) — any retry of the
   parked diff needs the acceptance test named in `TODO.md`: failure mode must
   not regress from fast/named refusal to slow/unrelated crash.
3. `splitTwoArcLens` and the piece-5-ray item (`TODO.md`) — same family as
   above, no quick reproducer; swept for staleness 2026-08-27, still open.

## Files
- `region.m` — tonight's 2 real fixes: `certifiesNonPositive` (G17), `maxAffineOverRegion` (item 3)
- `.claude/suite.sh` — new `--coverage` flag (Cobertura report)
- `MORNING.md` — full overnight-run narrative, only read once

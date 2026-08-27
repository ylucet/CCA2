# Session Handoff

_2026-08-27_

## Blocked
- Phase C1 per-term cost target — EXTERNAL, needs Yves. Box terms are 0.01 s;
  the 40–60 s figure is stale.

## State
- Branch `main` @ `2e18333` — "docs: G17 LOCATED -- the hole appears in fold 4's
  maximumP, after the pairing"
- Pushed: yes — `origin/main`
- Tests: fast 309/0 (08-27) · slow 93/0 (08-26, STALE: predates B3's LINE+EMPTY)
  · verylong 36/8/1 (08-25, STALE: predates the legacy pins)
- Known reds: two, both `plq_1p` — G17 (`rectMaximumIsTheConjugateOfTheWholeDomain`,
  a hole at (-0.5,2)) and `rectBiconjugateIsAConvexUnderestimator` (the scaling
  defect). The seven legacy `plq_1piece` reds are pinned → incomplete, not failed.

## Next
1. G17, one level deeper: count coverage across `maximumP`'s own split/merge
   stages at FOLD 4 — the fold is already identified. ~15 min per run.
2. Run `--slow` and `--verylong`; neither has run since B3's LINE/EMPTY landed.
3. Scaling (`quadFacet_exactAnotInB`, 374 vs 37): hand-check ONE pair on PRect3
   (44 s) before coding, then extend `maxAffineOverRegion` to unbounded regions
   via the recession cone.

## Files
- `TODO.md` — G1/G4/G10/G16/G17, each with its measurement
- `DECISIONS.md` — 08-26/27 entries; FOUR refutations, read before proposing
- `.claude/assembly_attempt_2026-08-25.diff` — parked assembly prerequisites

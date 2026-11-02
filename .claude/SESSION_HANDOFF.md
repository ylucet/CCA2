# Session Handoff

_2026-11-02_

## Blocked
- Push: `gh auth login` needed (no GitHub credentials on the old machine).

## State
- Branch `main` @ `c0bd4f6` — "perf(subdivision): overlay witnesses in doubles; corner charts only where a row vanishes"
- Pushed: no (credentials)
- Tests: fast 813/0, normal 32/0 (2026-11-02); slow last whole 2026-10-17, the subdivision slow suites 2026-10-31
- Known reds: none

## Next
1. Example 2 overlay: ~4 h, ~12000 faces (`.claude/ex2Overlay.m`). Cut the face COUNT first:
   13 of 21 pieces on tuple 1 have no witness (tie rows unproved). DECISIONS 2026-11-02 (both).
2. Step d (TODO item 1): cut the triangles AFTER the overlay (no cut lines in the arrangement);
   biconjQ subdivides Example 2; slow test vs EnvelopeMin on 4000 points.
3. Daily slow run (`bash .claude/suite.sh --slow`) — not run since the local arrangement.

## Files
- `envelopeSubdivisionAlg.m` — dropBelow, minOverlay, branchAt
- `algFaceBelow.m` — one face certified below another
- `algRowFails.m` — big-integer rows; no D sign
- `.claude/ex2Overlay.m` — regenerates Example 2's run
- `TODO.md` item 1 c'/d; `DECISIONS.md` 2026-11-02

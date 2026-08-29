# Session Handoff

_2026-08-29_

## Blocked
Nothing external.

## State
- Branch `main` @ `dde2a93` — "docs: record the getVertices memoization
  win in TODO.md's G4 entry"
- Pushed: yes — `origin/main`. Tag `v0.2` cut and pushed (before the two
  memoizations below — consider `v0.2.1` if that matters).
- Tests (2026-08-29): fast 312/0 · normal 12/0 · slow 96/0 (543s, was
  798s) · regionTest 27/0 · testcPLQ (verylong) 8/0/1 (G11 quarantined).
- Known reds: none silent. G11 quarantined via `assumeFail`.

## Next
1. Fold-strategy redesign (TODO G4): root cause found (high-degree hub
   vertex). TWO safe memoizations landed (`unionIsExact`, `getVertices`) —
   reference fixture TOTAL time 2944s -> 2008s, cell count unchanged (58).
   The real fix (resolve each hub's fan once, not pairwise) is still open.
2. `splitTwoArcLens`: reproducer FOUND. Target cell is an INTERMEDIATE
   `clipByFace` product — trace forward from operands, not backward.
3. Conic-conic closed form: `conicMeet.m`/`ratQ.m` exist, tested, unused.
4. QPLIB-shaped timing — spike/SCIP's side now (notified, scip-09, no
   reply yet).

## Files
- `region.m` — `unionIsExact` + `getVertices` both memoized this session;
  5 real `Inf`/needle bugs fixed in `maxAffineOverRegion` earlier
- `maxQuaPar.m` — G1/G4/G10 landed
- `testcPLQ.m` — G11 quarantined via `assumeFail`
- `SCIP_READINESS.md` — gate rewritten, items 1/2/4 met
- `DECISIONS.md`, `TODO.md` — extensive, all evidenced

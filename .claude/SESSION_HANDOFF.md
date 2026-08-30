# Session Handoff

_2026-08-29_

## Blocked
Nothing external.

## State
- Branch `main` @ `89648d4` — "docs: region.plus has zero redundancy --
  memoization thread stops on a clean negative"
- Pushed: yes — `origin/main`. Tag `v0.2` predates all three
  memoizations below — consider `v0.2.1`.
- Tests (2026-08-29): fast 312/0 · normal 12/0 · slow 96/0 (493s, was
  798s) · regionTest 27/0 · testcPLQ (verylong) 8/0/1 (G11 quarantined).
- Known reds: none silent. G11 quarantined via `assumeFail`.

## Next
1. Fold-strategy redesign (TODO G4): hub-vertex root cause found. THREE
   safe memoizations landed — reference fixture 2944s -> 1830s (~38%),
   cell count unchanged (58). `region.plus`/`mtimes` CHECKED, zero
   redundancy (0/246 dup keys) — not a memoization candidate. Real fix
   (resolve each hub's fan once) still open.
2. `splitTwoArcLens`: reproducer FOUND. Target cell is an INTERMEDIATE
   `clipByFace` product — trace forward, not backward.
3. Conic-conic closed form: `conicMeet.m`/`ratQ.m` exist, unused.
4. QPLIB-shaped timing — spike/SCIP's side (notified, scip-09, no reply).

## Files
- `region.m` — THREE memoizations this session (`unionIsExact`,
  `getVertices`, `simplifyUnboundedRegion`); 5 `Inf`/needle bugs fixed
  earlier in `maxAffineOverRegion`
- `maxQuaPar.m` — G1/G4/G10 landed
- `testcPLQ.m` — G11 quarantined via `assumeFail`
- `SCIP_READINESS.md` — gate rewritten, items 1/2/4 met
- `DECISIONS.md`, `TODO.md` — extensive, all evidenced

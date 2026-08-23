# Session Handoff

_2026-08-23_

## Blocked

- **Phase C1 per-term cost target** — EXTERNAL, needs Yves. Box terms are 0.01 s,
  so the old 40–60 s figure is stale.

## State

- Branch `main` @ `db7d13b` — "The convex envelope is not rational"
- Pushed: pending — 12 commits `da2c5b9..db7d13b`, 10 of them `proof/` from
  another session
- Tests (2026-08-20): fast 249/0 · normal 12/0 · slow 88/0 · verylong NOT run.
  Nothing run since — no VPN, and no `.m` file changed this session.
- Known reds: `testMaxMultiRegion/testPCE2` — `plq_1piece`'s envelope route

## Next

1. **`Con` trait** — relax `RatPar`'s `set.Ec`; cheapest item, and `f*`'s
   elliptical edge already forces it.
2. **T3** — `symbolicFunction`'s payload to a rational coefficient vector;
   order and counts in `TODO.md` top section.
3. `bash .claude/suite.sh --verylong` before the next tag (~73 min).

## Files

- `DECISIONS.md` 2026-08-23 — envelope face degree 1/2/4; `f**` is not C1
- `TODO.md` 2026-08-23 — `QuaCon`/`AlgCon` return types, in order
- `proof/` — another session's Lean proof of `conj_isQuaCon`, 0 sorry
- `RatPol.m` header — now wrong, it describes the one-piece case only

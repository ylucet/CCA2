# Session Handoff

_2026-09-06_

## Blocked
- The IMPLICIT-FACE TYPE (`p(x,y,z)=0` + branch selector). External: your call on
  adding a third value to `RETURN_TYPE.md`'s function axis. DECISIONS 2026-09-06.

## State
- Branch `main` @ `ca208da` — "measure: ASSEMBLY's cost is the cell-pair count"
- Pushed: yes
- Tests: fast 563/0 (09-06) · normal 13/0 (09-05) · slow 98/0 (09-06, **stale**:
  predates the `ratQ` big-layer rewrite)
- Known reds: none

## Next
1. `bash .claude/suite.sh --slow -j 4` (~13 min). Stale gate, run it first.
2. REGION LOGIC — which primal region a bridge face owns. Untouched, and the
   only part of assembly with real content. Then one full envelope on a triangle
   pair (~24 eliminations x 1-19 s each).
3. If that is too slow: second filter (does `g1-g2` change sign on the cell
   intersection), then cache cells by type. First filter cut only 1/3 to 1/2.

## Files
- `ratQ.m` — resultant, BigInteger fallback, exact division
- `.claude/bridgeCellProbe.m` — all four cell-pair kinds, verified
- `.claude/cellPairFilterProbe.m` — cell-pair counts and the filter
- `DECISIONS.md` — six refutations recorded this session

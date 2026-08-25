# Session Handoff

_2026-08-25_

## Blocked
- Phase C1 per-term cost target — EXTERNAL, needs Yves. Box terms are 0.01 s;
  the 40–60 s figure is stale.

## State
- Branch `main` @ `cd66cc2` — "docs: REFUTED -- branching cannot isolate
  two runs that share a working tree"
- Pushed: yes — `origin/main`. Every other branch deleted; all were merged.
- Tests (2026-08-25): fast 303/0 · slow 88/0 · verylong 26 pass / 7 fail /
  1 timeout, that last figure IDENTICAL to a pristine `b9243d3`.
- Known reds: the seven verylong ones, all pre-existing (`testPCE2` among them).
- Hazard: `--verylong -j N` races on `plqStage`'s cache and can fake a red.
  Re-run a suspect at `-j 1`. `TODO.md` G7 has the fix.

## Next
1. G1 — `clipByFace` returns nothing for a face pair whose intersection is not
   empty (measured: g1 face 4 × g2 face 2). Last bounded fallback.
2. G2 — affine face over an unbounded polygon (`max(0,x,y)`); `TODO.md` prices
   an all-affine route that never enters `maxQuaPar`.
3. G4 — `conj` of `xy` on some triangles computes a MINORANT. It now raises
   instead of returning it, so the refusal is the safety net, not the fix.

## Files
- `TODO.md` — opens with the measured gap list G1–G7.
- `MORNING.md` — the overnight run's report; `proof/MORNING.md` is the other's.
- `DECISIONS.md` — nine entries dated 2026-08-24/25; read the headings.
- `checkConjSymFree.m` — the symbolic-fallback rate, with reasons.

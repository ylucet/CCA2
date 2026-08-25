# Session Handoff

_2026-08-24_

## Blocked

- **Phase C1 per-term cost target** — EXTERNAL, needs Yves. Box terms are 0.01 s,
  so the old 40–60 s figure is stale.

## State

- Branch `overnight/2026-08-24`, an unmerged fast-forward off `main`. **Read
  `MORNING.md` first** — it is the report for that run. A parallel `proof/`
  session commits to `main`.
- Pushed: no. The overnight run never pushes.
- Tests (2026-08-24, on the branch): fast **296 / 0**, slow **88 / 0** (identical
  to its pre-change baseline), verylong NOT run.
- Known reds: none.

## What changed, in one line

`conj`'s numeric path was ALREADY sym-free; the work was shrinking the set of
inputs that fall back to the symbolic Case C. Measured with `checkConjSymFree`:
3 of 17 fixtures still fall back, and the unbounded CONVEX family moved from
"no numeric route at all" to 0.16 s.

## Next

1. **G2b** — `maxQuaPar` DROPS a cell on some unbounded folds. The one silent
   wrong answer of the session; a cross-check catches it and falls back, so it
   cannot reach a caller, but it is unfixed. Reproducer in `TODO.md`.
2. **G1** — the arc split works; the pieces then fail to pair up in
   `assemblePieces`. Attack the MATCHING, not the refusal.
3. **G2** — an affine face over an unbounded polygon (`max(0,x,y)`).
4. `bash .claude/suite.sh --verylong` before the next tag.

## Files

- `MORNING.md` — the overnight report; the gap list is `TODO.md` G1–G3.
- `DECISIONS.md` 2026-08-24 — four entries: the precision budget, the `syms`
  assumption leak, the arc-split geometry, the dropped cell.
- `SUPPORT_MATRIX.md` §1.2 and §4 refreshed; §4.4 is new.
- `proof/` — another session's Lean proof; do not stage its files.

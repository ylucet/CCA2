# Session Handoff

_2026-08-24_

## Blocked

- **Phase C1 per-term cost target** — EXTERNAL, needs Yves. Box terms are 0.01 s,
  so the old 40–60 s figure is stale.

## State

- Branch `main`. **Read `MORNING.md` first** — it is the report for the overnight
  run, and its first section explains why the branch that run created was a
  mistake and has been folded back in. A parallel `proof/` session also commits
  to `main`; do not stage its files.
- Pushed: no. The overnight run never pushes.
- Tests (2026-08-25): fast **303 / 0**, slow **88 / 0**, verylong
  **26 pass / 7 fail / 1 timeout** — and that verylong figure is IDENTICAL to a
  pristine `b9243d3`, so its reds are all pre-existing (`testPCE2` among them).
- Known reds: the seven pre-existing verylong ones. Nothing this run caused.
- One tooling hazard: `--verylong -j N` RACES on `plqStage`'s shared cache and
  can produce a spurious red. Re-run a suspicious test at `-j 1` before believing
  it. `TODO.md` G7 has the three-line fix.

## What changed, in one line

`conj`'s numeric path was ALREADY sym-free; the work was shrinking the set of
inputs that fall back to the symbolic Case C. Measured with `checkConjSymFree`:
3 of 17 fixtures still fall back, and the unbounded CONVEX family moved from
"no numeric route at all" to 0.16 s.

## Next

1. **G1 — build the missing LENS.** Two operands can have boundaries between the
   SAME two dual points, one straight and one curved; the arrangement is missing
   the cell between them, and the orphaned arc in `assemblePieces` is looking for
   a partner that was never built. Do NOT chase the matching. This is the last
   fallback of the bounded family.
2. **G2** — an affine face over an unbounded polygon (`max(0,x,y)`). `TODO.md`
   prices a direct route for all-affine input that never enters `maxQuaPar`.
3. G2b is **DONE**: the silent wrong answer on unbounded folds is fixed (the
   split direction at a line pair's singular point).

## Files

- `MORNING.md` — the overnight report; the gap list is `TODO.md` G1–G3.
- `DECISIONS.md` 2026-08-24 — four entries: the precision budget, the `syms`
  assumption leak, the arc-split geometry, the dropped cell.
- `SUPPORT_MATRIX.md` §1.2 and §4 refreshed; §4.4 is new.
- `proof/` — another session's Lean proof; do not stage its files.

# Morning report — 2026-09-06 overnight run

Branch: `overnight/2026-09-06`

Started from `main` @ `a57a3e2`, tree clean.

**On the branch.** Stored feedback says "commit on main, don't branch" for this repo. This run
branches anyway because `/overnight`'s safety model rests on it. Say the word and I will merge
straight across to `main` — the history is linear.

**Daily gate.** `--slow -j 4` already ran earlier in this session: **98 passed, 0 failed**, which
clears the stale gate the handoff flagged. `--verylong -j 4` was started in the background at the
top of this run, once.

The four items came from Yves directly rather than from `TODO.md`'s **Next up**, and were written
into `TODO.md` as a **Next up** block at the start of the run.

## What changed

- **Item 1 DONE — two `N/R` guards are reachable, one survives.** `SUPPORT_MATRIX.md` classified
  guards as unreachable on a scope claim that was refuted yesterday. Re-tested each guard's
  CONDITION on real `conjQ` output (`.claude/guardConditionProbe.m`), over eight two-piece kink
  fixtures with both single-piece fixtures as a control: **curved ray edges 52** (control 0),
  **irreducible ellipse/hyperbola differences 44** (control 0), **faces with more than one curved
  edge 0**. The first two rows are reclassified; the third survives the re-check.
  `d3e0cba`

## What is broken

## Needs a decision

## Where I stopped

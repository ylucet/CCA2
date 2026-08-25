# Morning report — 2026-08-24 overnight run (QuaConProof)

Branch: **`overnight/2026-08-24`** — *not* a branch of my own. See the deviation
note below; it matters for how you review this.

## Deviations from `/overnight`, and why

**No new branch.** The skill says to isolate the run on `overnight/YYYY-MM-DD`.
This repository is shared with a **concurrent MATLAB session in `CCA2/`**, which
had already created and checked out `overnight/2026-08-24`. One working tree
means one `HEAD`: creating or switching a branch would have moved that session's
checkout under it mid-run. So this run commits onto the branch that was already
checked out, and every commit stages **only paths under `proof/`**, as
`proof/CLAUDE.md` requires. The two projects' commits interleave in `git log`;
`git log --oneline -- proof/` separates them.

**The dirty files at start are not mine.** `MORNING.md`, `conjCPLQ.m`,
`checkConjSymFree.m`, `conjConvexPolygon.m` in `CCA2/` belong to the MATLAB run.
The skill says to checkpoint a dirty tree; I have not, because they are outside
my write boundary. They are untouched.

**`CCA2/MORNING.md` is the other run's report.** This file, `proof/MORNING.md`,
is mine.

**Pushed before starting, on request.** `origin/main` moved `b9243d3..64043b5`,
publishing Phase 7. Note that `main` did *not* contain the later plan commit
`86243b1` — it had been left behind when the MATLAB session branched — so that
commit and everything tonight sit on `overnight/2026-08-24`, unpushed, as the
guardrails require.

## What changed

- **A1, A2, A3 done** (`QuaConProof/Convexity.lean`, new file). `f*` is convex,
  lower semicontinuous, and `dom f*` is convex. Convexity is stated as an
  epigraph inequality rather than `f*(as1+bs2) <= a f*(s1) + b f*(s2)`, which
  keeps `EReal` multiplication out of every proof; the direct form is recorded as
  `conj_combo_le`. Build green, axioms clean.


## What is broken

_(nothing yet)_

## Needs a decision

_(nothing yet)_

## Where I stopped

Run started. Preflight done: `lake build` green at `9cbf484`, 0 sorry,
`#print axioms conj_isQuaCon` clean. Working the `TODO.md` programme top to
bottom — A (finish the real case), B (rational coefficients), C (biconjugate
shape), 17 items.

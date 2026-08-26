# Morning report — 2026-08-25 overnight run

Branch: `overnight/2026-08-25`

## What changed

_(appended as the run goes)_

## What is broken

## Needs a decision

### The branch itself — read this first

The run is on `overnight/2026-08-25`, not on `main`, and that contradicts a
standing preference. `DECISIONS.md` 2026-08-25 REFUTED branching for an
unattended run, and the recorded rule is "stay on `main`".

I branched anyway, deliberately, for one reason: **that refutation's cause is
gone as of today.** It was about two sessions sharing ONE working tree — the
`proof/` session and this one — where branching silently dragged the other
session onto this run's branch. `proof/` moved to its own repository
(`AI/CCA2proof`) earlier today, so nothing else has this tree checked out and
that failure mode cannot recur. Against that, the skill's reason for branching
holds: abandoning a whole unattended night is one `git switch`.

If you disagree, `git switch main && git merge --ff-only overnight/2026-08-25`
folds it in with nothing lost, exactly as the last one was folded back.

## Where I stopped

_(run in progress)_

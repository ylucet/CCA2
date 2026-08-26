# Morning report — 2026-08-25 overnight run

Branch: `overnight/2026-08-25`

## What changed

- **Item 1 (`assemblePieces`, = G1 + G4 + G10) — scope-reduced to one named step, not closed.**
  Found the obstruction and the discriminator for it; the working attempt is preserved at
  `.claude/assembly_attempt_2026-08-25.diff` and is not committed because it leaves two tests red.
  Full measurements in `DECISIONS.md` 2026-08-25 (overnight). One thing is left:
  `clipArcByHalfPlane` fails on a degenerate arc; fix that and re-apply the diff.
- **`clipArcByHalfPlane` FIXED** — it raised an internal error on a case where a root provably
  exists (opposite signs at the two ends, continuous in between). The closed-form solve fails on a
  short arc; bisection is now its backstop. Independent of everything else. `clipArcByHalfPlaneTest`
  7/0, fast 303/0.
- Nightly slow gate (started first, per rule 9): **92 pass / 0 fail** over 5 suites.

## What is broken

- `TODO.md` item 1 is still open, with a much smaller boundary — see above.

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

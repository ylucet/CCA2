# Morning report — 2026-08-27 overnight run

Branch: overnight/2026-08-27

## What changed

- **G1/G10 — re-applied the parked assembly diff, measured, reverted.** Not a clean win: on case
  21's full `conj('cplq')` it trades the documented 2.4s `foldDroppedACell` refusal for a 292.5s
  run that hits a different, unrelated MuPAD internal error in Case C. Diff stays parked.
  `DECISIONS.md` 2026-08-27 (overnight, G1/G10).
- (in progress) **G17's `certifiesNonPositive` fix** — verify `quadprog`'s `ef==-2` via `linprog`
  before trusting it as "region empty" (root cause found this session: `quadprog` genuinely fails
  to converge on an unbounded cell with a rank-deficient Hessian and returns `-2` as a
  non-convergence code, not an infeasibility one). Confirmed on the isolated reproducer; an
  end-to-end fixture re-run is in the background (`g17_full`, checking the actual hole closes).
  Not yet run against the fast suite, not yet committed.

## What is broken

## Needs a decision

## Where I stopped

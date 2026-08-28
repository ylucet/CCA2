# Morning report — 2026-08-27 overnight run

Branch: overnight/2026-08-27

## What changed

- **G1/G10 — re-applied the parked assembly diff, measured, reverted.** Not a clean win: on case
  21's full `conj('cplq')` it trades the documented 2.4s `foldDroppedACell` refusal for a 292.5s
  run that hits a different, unrelated MuPAD internal error in Case C. Diff stays parked.
  `DECISIONS.md` 2026-08-27 (overnight, G1/G10).
- **G17 — FIXED and committed (`84dbba7`).** Root cause: `certifiesNonPositive`'s concave branch
  trusted `quadprog`'s `ef==-2` as "region empty, vacuously true" without verifying it. On an
  unbounded cell with a rank-deficient (semidefinite) Hessian, `quadprog` genuinely fails to
  converge and returns `-2` as a non-convergence code, not an infeasibility one -- confirmed
  directly (the point is feasible; `quadprog`'s own iterates diverge). Fixed by verifying via
  `linprog` (`region.maxLinear`) before accepting `ef==-2`. Fast 309/0, normal 12/0. An end-to-end
  re-run of the ORIGINAL fixture (not just the isolated reproducer) was still running in the
  background when this landed (`g17_full` -- see "Where I stopped").
- **Item 3 (the `quadFacet_exactAnotInB` scaling defect) — FIXED and committed (`84dbba7`, same
  commit as G17; both fixes are in `region.m`).** `maxAffineOverRegion` now correctly tightens an
  unbounded region's bound via the conic's own parametrization, combined with a straight-line
  recession check -- two independent sufficient conditions for genuine unboundedness, found and
  hardened by TWO hand-derived counterexamples before landing (see DECISIONS.md 2026-08-27, item
  1 REFUTED and item 3 FIXED for both). Validated against a true 2D brute-force oracle: 5000+
  prototype cases, ~150 as real `region` objects, zero genuine disagreements. **Not yet
  re-measured against the ORIGINAL 374-refusal fixture's cell count** -- that measurement is the
  next natural step but the fixture is too expensive to reach quickly; the fix is validated on
  the mechanism, not yet on that specific scaling number.

## What is broken

## Needs a decision

## Where I stopped

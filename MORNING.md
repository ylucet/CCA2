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
- **Item 3 (the `quadFacet_exactAnotInB` scaling defect) — the FIX is correct and committed
  (`84dbba7`), but does NOT close the scaling defect it targets.** `maxAffineOverRegion` now
  correctly tightens an unbounded region's bound via the conic's own parametrization combined
  with a straight-line recession check -- two independent sufficient conditions for genuine
  unboundedness, found and hardened by TWO hand-derived counterexamples before landing, validated
  against a true 2D brute-force oracle (5000+ prototype cases, ~150 as real `region` objects,
  zero genuine disagreements). **Then re-measured overnight against the ORIGINAL fixture
  (`.claude/step3cost.m`) and it is essentially a wash: 58 cells vs the 56-cell baseline, no
  material improvement**, even though the `exactAnotInB` refusal RATE genuinely dropped (6.9% vs
  11.6%). Likely cause: the fix only handles a region with EXACTLY ONE curved facet by design,
  and this fixture's surplus cells plausibly carry two or more. Recorded honestly rather than
  overclaimed -- `DECISIONS.md` 2026-08-27 (overnight, item 3 vs the ORIGINAL fixture), commit
  `ca4c025`. The scaling defect stays OPEN in `TODO.md` with the next concrete step (instrument
  `unionIsExact` for `numel(qidx)>1` before attempting the two-conic case).

## What is broken

- **The `quadFacet_exactAnotInB` scaling defect (Step 3's cross-piece maximum on the A.4/A.5
  quadrilateral) is still open.** Item 3's fix closed a real, independently-validated soundness
  gap but did not move this fixture's cell count. Boundary: the fix's own scope is exactly one
  curved facet per region; extending to two-or-more is a materially harder proof (combined
  recession behaviour of two conics) and was not attempted.

## Needs a decision

## Where I stopped

# Morning report — 2026-08-27 overnight run

Branch: overnight/2026-08-27

## What changed

- **G17 — FIXED and committed (`84dbba7`).** Root cause: `certifiesNonPositive`'s concave branch
  trusted `quadprog`'s `ef==-2` as "region empty, vacuously true" without verifying it. On an
  unbounded cell with a rank-deficient (semidefinite) Hessian, `quadprog` genuinely fails to
  converge and returns `-2` as a non-convergence code, not an infeasibility one -- confirmed
  directly (the point is feasible; `quadprog`'s own iterates diverge). Fixed by verifying via
  `linprog` (`region.maxLinear`) before accepting `ef==-2`. Fast 309/0, normal 12/0, and the exact
  isolated reproducer (cells 12/17/21) now correctly refuses the bad merge. **The end-to-end
  re-run of the ORIGINAL fixture was ABANDONED after 12 continuous CPU-bound hours** (killed, not
  hung -- confirmed via CPU-time tracking) -- not a correctness concern, but the honest fix's real
  verification work is far more expensive on this fixture than the old wrong-but-fast merge was.
  Left `UNVERIFIED` at the full-pipeline level; the isolated reproducer is what stands as
  verification. `DECISIONS.md` 2026-08-27 (G17, ROOT CAUSE FOUND) and (overnight, G17 end-to-end).

- **Item 3 (the `quadFacet_exactAnotInB` scaling defect) — the FIX is correct and committed
  (`84dbba7`), but does NOT close the scaling defect it targets.** `maxAffineOverRegion` now
  correctly tightens an unbounded region's bound via the conic's own parametrization combined
  with a straight-line recession check -- two independent sufficient conditions for genuine
  unboundedness, found and hardened by TWO hand-derived counterexamples before landing, validated
  against a true 2D brute-force oracle (5000+ prototype cases, ~150 as real `region` objects,
  zero genuine disagreements). **Re-measured overnight against the ORIGINAL fixture
  (`.claude/step3cost.m`) and it is essentially a wash: 58 cells vs the 56-cell baseline**, even
  though the `exactAnotInB` refusal RATE genuinely dropped (6.9% vs 11.6%). Recorded honestly
  rather than overclaimed -- `DECISIONS.md` 2026-08-27 (overnight, item 3 vs the ORIGINAL
  fixture), commit `ca4c025`.
  **Follow-up, CONFIRMED the reason:** instrumented `unionIsExact` directly (temp probe, reverted
  immediately after, nothing shipped) -- of 141 bare `exactAnotInB` refusals on the target
  fixture, 60% (84) carry TWO OR MORE curved facets, outside `tightenUnboundedFacet`'s deliberate
  one-conic scope. The 39% single-conic share is genuinely helped (`okQuadFacet` successes rose).
  `DECISIONS.md` 2026-08-27 (overnight, item 3 follow-up), commit `464a04c`. The two-conic case
  is scoped as a separate, harder future item in `TODO.md` (no single global arc parameter once
  two conics are both active constraints).

- **G1/G10 — re-applied the parked assembly diff, measured, reverted.** Not a clean win: on case
  21's full `conj('cplq')` it trades the documented 2.4s `foldDroppedACell` refusal for a 292.5s
  run that hits a different, unrelated MuPAD internal error in Case C. Diff stays parked.
  `DECISIONS.md` 2026-08-27 (overnight, G1/G10).

- **Stale TODO.md entry found and checked off.** "The parallelogram's LAST 2 of 10" (dated
  2026-08-16) described a bug that `functionNDomain.singularEdgeCut` (commit `699326d`, the very
  next day) already fixed; nobody had checked the box. Confirmed still passing tonight
  (`aBoundedPieceWithATangentVertexConjugatesOntoTheWholePlane`, 1/0/0). Commit `378a45b`. Swept
  the two remaining open items after this (`splitTwoArcLens`, the piece-5-ray item under an
  explicitly `SUPERSEDED` 2026-08-03 status section) -- both are the same difficult
  `assemblePieces` family already explored tonight via G1/G10, without a quick reproducer the way
  the parallelogram item had; not attempted fresh tonight rather than starting a multi-hour
  investigation this late in the run.

- **Code coverage tooling — DONE.** `.claude/suite.sh --fast --coverage` / `--normal --coverage`
  via `matlab.unittest.plugins.CodeCoveragePlugin`, Cobertura-XML report to
  `.claude/coverage/cobertura.xml` (gitignored). Verified: `--fast --coverage` matches plain
  `--fast`'s 309/0/0 with negligible overhead, and produced a real baseline -- **47.4% line
  coverage from the fast bucket alone** (9825/20719 lines). Commit `9ad2acc`.

## What is broken

- **The `quadFacet_exactAnotInB` scaling defect (Step 3's cross-piece maximum on the A.4/A.5
  quadrilateral) is still open.** Item 3's fix closed a real, independently-validated soundness
  gap but did not move this fixture's cell count -- confirmed why (60% of remaining refusals need
  the two-conic case, out of scope). Boundary for a future attempt: the region's true recession
  cone becomes the INTERSECTION of both conics' recessive-direction sets, and there is no single
  global parameter analogous to `parabolaArcFrame`'s `u` once two arcs are both active
  constraints -- a materially harder proof than the one just built.
- **G17's full-pipeline fixture is UNVERIFIED**, not because of any known defect, but because the
  honest fix's real verification cost made a 12-hour run inconclusive on wall clock alone. The
  fix itself (the `certifiesNonPositive` routine) is verified by the isolated reproducer and both
  regression buckets; only this one specific downstream fixture's END-TO-END completion time is
  open.
- **G1/G4/G10** (the `assemblePieces` curved/straight matching defect) remains open. The parked
  diff (`.claude/assembly_attempt_2026-08-25.diff`) does not close it cleanly on the actual
  `conj('cplq')` call, only on the narrower assembly-only test it was originally measured
  against -- see the acceptance-test requirement recorded in `TODO.md` for the next attempt.

## Needs a decision

None outstanding -- every choice made tonight had a reversible option and was taken (see
`DECISIONS.md` for the record of each: G1/G10's revert, G17's kill-after-12h, item 3's honest
non-improvement report).

## Where I stopped

Worked through, in order: G1/G10 (measured, reverted), G17 (root cause found and fixed), item 3
(fixed, then honestly re-measured against its own target and found wanting, then the shortfall's
cause confirmed), code coverage (built and verified), a stale TODO.md entry (found and checked
off). 15 commits deep on `overnight/2026-08-27`, tree clean, fast 309/0 and normal 12/0 both
green as of the last real code change (commit `84dbba7`; every commit after that is docs/tooling
only). Stopped here after sweeping the two remaining open items and finding neither had a quick
reproducer the way the stale one did -- both are the same difficult family already worked
tonight without a clean win, and starting a fresh multi-hour investigation into either did not
look like a good use of what was left of the run.

**Next, in the order I would take it:**
1. The two-conic extension to item 3's fix, if the scaling defect is worth more time -- it is a
   new, harder proof, not a continuation of tonight's work, and `TODO.md`/`DECISIONS.md` both
   sketch what it would need.
2. G1/G4/G10's `assemblePieces` defect, with the acceptance test tonight's measurement named
   (failure mode must not regress from fast/named to slow/unrelated) as the bar for any future
   attempt at the parked diff or a fresh one.
3. Re-running `--slow`/`--verylong` against tonight's `region.m` changes -- both ran once earlier
   in this session (before tonight's fixes), so a fresh run would catch anything the fast/normal
   buckets cannot see. Not started tonight; would cost the better part of an hour given `--slow`'s
   own budget, and the night's remaining time went to the measurements above instead.

No blockers reached `TODO.md`'s `## Blocked` section; nothing here needed a decision only Yves
can make.

# Session Handoff

_Last updated: 2026-07-19T00:00:00Z_

## What happened this session

Fixed the tie-point coverage gap the prior session's handoff flagged as the one remaining
biconjugate bug — but the actual root cause was different from what was documented: not
`functionNDomain.mergeL`/`region.removeTangent`, but `region.maxArray`'s dominance-probing
heuristic, which gives up (`lSing`) when a region's only vertices lie exactly on the `f1==f2` tie
line (e.g. two independent triangles' conjugates overlapping in a lens bounded by two parabolas
meeting exactly at the tie point). `functionNDomain.maximumP` was silently discarding that whole
overlap region instead of falling back to `splitmax3`, producing a real coverage gap at the tie
point, not a tolerance artifact. Fixed by falling through to `splitmax3` on `lSing`, mirroring what
`maxEqDom` already does. Fixing that exposed 3 further latent empty-region crashes (in
`mergeL`'s merge-accumulation/copy-through paths and `region.merge` itself) that were previously
unreachable because the tie-point bug quietly dropped exactly the cases that trigger them; all
guarded with `isempty` checks. `cplqAdapterTest.m`/`conjCPLQTest.m` no longer exclude
`s=(0.5,0.5)`. Ran the full `testMaxMultiRegion` suite (previously only ~12/24 cases had ever been
exercised) against a captured pre-fix baseline to separate genuine regressions from pre-existing
failures — all resolved regressions confirmed fixed, no new ones introduced.

## Where things stand

- Branch: `main` @ `f715df2` — "Fix real coverage gap at exact conjugate tie points -- cPLQ
  multi-face maximum now handles degenerate overlaps"
- Pushed: no upstream configured
- `testMaxMultiRegion` full suite: 19/24 pass. The 5 failures (`testPCE2`, `testFractional`,
  `testMaxThesis`, `testMaxThesis2`, `testOpenconvex`) are confirmed pre-existing — they fail
  identically on the pre-fix baseline, unrelated to this session's changes. **These are now the
  next thing to fix** (see Next steps).
- `cplqAdapterTest.m`/`conjCPLQTest.m`: 20/20 pass, including both tie-point assertions.
- `testRegion`/`testfunctionNDomain`/`testSymbolicFunction`: pass except
  `testRegion/testCreation`, confirmed failing identically pre-fix (pre-existing, the
  known `isequal(sym,double)` toolbox-compat gap — unrelated).
- `cPLQ/` (the original reference clone) remains intentionally untracked, per explicit user
  instruction — do not `git add` it.

## Next steps

- **Fix the 5 remaining `testMaxMultiRegion` failures**: `testPCE2`, `testFractional`,
  `testMaxThesis`, `testMaxThesis2`, `testOpenconvex`. Error messages captured this session
  (see scratch logs referenced below, though those are session-temp and may be gone — rerun
  `runtests('testMaxMultiRegion')` to reproduce). `testOpenconvex`'s error was captured directly:
  `functionNDomain/maxEqFun` (line 608, `eqFun.*eqFun`) → `region.minus` (line 352,
  `v1(i,:,:) = obj1.getEndpoints(i)`) → `sym/privsubsasgn` — a `1x2x2` vs `4x2` size mismatch
  assigning into a sym array, likely a different bug in the same "unguarded shape assumption"
  family as this session's fixes. The other 4 weren't diagnosed as deeply — start with a
  fresh `runtests` capture per case, same before/after-fix-diff method this session used to
  distinguish real regressions from pre-existing failures.
- **Give Case C a proper `QuaPar`-like return type** (convert `evalFunctionNDomain`'s underlying
  `functionNDomain` result into something with its own `.conj()`), so `biconj`/`infConv`/`moreau`/
  `proxAverage` can compose with a genuinely multi-face nonconvex conjugate, not just evaluate it.
- **Phase 2 (later, not started)**: once Phase 1 is fully solid, replace `cPLQ`'s symbolic
  computation with closed-form numeric formulas incrementally, one case/step at a time, validating
  each against the working Phase-1 symbolic result before moving to the next — this is what
  actually fixes the slowness (`isAlways` "truth unknown" retries; a single biconjugate/maximum
  run over a nontrivial multi-piece domain can cost 5-20+ minutes).
- Unrelated, longstanding, lower-priority items (unaffected by this session, still open): exact
  `[LOCATELLI]` citation in `DESIGN.md`; 2/741 residual `maxQuaPar:internal` crashes;
  `QuaPar.orderEdges`/`createP`'s near-degenerate-triangle error; `partialConj` for
  `'cplq'`/`'pqp'`; `add` for `RatPol`/`RatPar`; conjugate engines `'pqp'`/`'graph'`.

## Relevant files

- `functionNDomain.m` — this session's core fix (`maximumP`'s `lSing` handling) plus the 3
  empty-region guards in `mergeL`/`maximumP`; search "HISTORY (tie-point fix..." comments.
- `region.m` — defensive `isempty` guard added to `merge`; search its own "HISTORY" comment.
- `cplqAdapterTest.m`/`conjCPLQTest.m` — tie-point tests, no longer excluding `s=(0.5,0.5)`.
- `testMaxMultiRegion.m` — the untested-case backlog is now down to the 5 listed above; each
  test method's fixture setup is in `setUpTestData` at the top of the file.
- `conjCPLQ.m` — Case C (general bounded multi-face/non-triangular domain via the `cPLQ`
  pipeline); see its own header for the `functionNDomain`-vs-`QuaPoly`/`QuaPar` return-type caveat
  (still open, see Next steps).
- `DESIGN.md` §0 point 3-4 and "Next planned"/II.5.1 — the full two-phase plan and current status.

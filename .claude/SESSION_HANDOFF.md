# Session Handoff

_Last updated: 2026-07-18T17:00:00Z_

## What happened this session

Redirected the nonconvex-PLQ-conjugate work from deriving closed-form numeric formulas to
**integrating the reference `cPLQ` symbolic package into CCA2** (Phase 1 of a two-phase plan;
Phase 2, later, replaces the symbolic computation with closed-form numerics once Phase 1 is solid
— see `DESIGN.md` §"Next planned"/II.5.1). Copied the 9-file runtime dependency closure of `cPLQ`
plus its test suites into the repo root, fixed 2 toolbox-version-compat bugs to get a passing
baseline, built a `QuaPoly`<->`plq` adapter (`quaPolyToPlq.m`/`evalFunctionNDomain.m`), wired it
into `conjCPLQ.m`'s previously-unimplemented general (multi-face) case, and fixed 2 real bugs in
the vendored `region.m` that were blocking biconjugate (`poly2orderUnbounded` index overflow,
`getNormalConeVertexQ` dimension mismatch). Biconjugate now passes on every test tried.

## Where things stand

- Branch: `main` @ `06f18a0` — "Fix region.getNormalConeVertexQ dimension-mismatch bug -- biconjugate now works"
- Pushed: no upstream configured
- A genuinely multi-triangle nonconvex `QuaPoly` can now be conjugated end to end via
  `q.conj('cplq')` (previously errored unconditionally for `nf>1`). The result (`g`) for that case
  is a cPLQ `functionNDomain` array, NOT a `QuaPoly`/`QuaPar` — evaluate with
  `evalFunctionNDomain(g, s)`, not `g.eval(s)`; composition (`biconj`/`infConv`/`moreau`/...) does
  not support this case yet for that reason (documented in `conjCPLQ.m`'s own header).
- The only known open bug from this integration: `functionNDomain.mergeL`/`region.removeTangent`
  drop an exact symmetric tie point (`s=(0.5,0.5)` in the test example) when merging two
  independent triangles' conjugates during `.maximum` — narrow (exact ties only), documented and
  excluded in `cplqAdapterTest.m`/`conjCPLQTest.m`, not fixed.
- `cPLQ/` (the original reference clone this session copied files FROM) is intentionally left
  untracked, per explicit user instruction — do not `git add` it.

## Next steps

- **Fix the `mergeL`/`removeTangent` tie-point gap** — the one remaining known bug from this
  integration; Step 3 (`.maximum`), not biconjugate.
- **Give Case C a proper `QuaPar`-like return type** (convert `evalFunctionNDomain`'s underlying
  `functionNDomain` result into something with its own `.conj()`), so `biconj`/`infConv`/`moreau`/
  `proxAverage` can compose with a genuinely multi-face nonconvex conjugate, not just evaluate it.
- **Run the ~12 still-untested `testMaxMultiRegion` cases** (`testMaxR3`/`testMaxT`/`testMaxP`/
  `testMax3`/`testPSqroot`/`testFractional`/`testMaxThesis`/`testMaxThesis2`/`testOpenconvex`,
  plus 6 of the 8 `testBiconjugate*` not yet exercised) to see what else is lurking before calling
  Phase 1 done.
- **Phase 2 (later, not started)**: once Phase 1 is fully solid, replace `cPLQ`'s symbolic
  computation with closed-form numeric formulas incrementally, one case/step at a time, validating
  each against the working Phase-1 symbolic result before moving to the next — this is what
  actually fixes the slowness (`isAlways` "truth unknown" retries; a single biconjugate run costs
  5-9 minutes).
- Unrelated, longstanding, lower-priority items (unaffected by this session, still open): exact
  `[LOCATELLI]` citation in `DESIGN.md`; 2/741 residual `maxQuaPar:internal` crashes;
  `QuaPar.orderEdges`/`createP`'s near-degenerate-triangle error; `partialConj` for
  `'cplq'`/`'pqp'`; `add` for `RatPol`/`RatPar`; conjugate engines `'pqp'`/`'graph'`.

## Relevant files

- `DESIGN.md` §0 point 3-4 and "Next planned"/II.5.1 — the full two-phase plan and current status,
  including the corrected `cPLQ` runtime-dependency list (search "IN PROGRESS, core adapter DONE").
- `conjCPLQ.m` — Case C (general bounded multi-face/non-triangular domain via the `cPLQ`
  pipeline); see its own header for the `functionNDomain`-vs-`QuaPoly`/`QuaPar` return-type caveat.
- `quaPolyToPlq.m` / `evalFunctionNDomain.m` — the adapter (input conversion / numeric eval).
- `region.m` — both bug fixes this session (`poly2orderUnbounded`, `getNormalConeVertexQ`); search
  "HISTORY" comments there for the reasoning.
- `plq.m`, `plq_1p.m`, `plq_1piece.m`, `functionNDomain.m`, `domain.m`, `symbolicFunction.m`,
  `conjugateExpr.m`, `yIntercept.m` — the rest of the copied `cPLQ` runtime dependency closure,
  otherwise unmodified.
- `cplqAdapterTest.m`, `conjCPLQTest.m`, `testMaxMultiRegion.m`, `testcPLQ.m`,
  `testRegion.m`/`testSymbolicFunction.m`/`testfunctionNDomain.m` — tests; the last three include
  the toolbox-compat fix (`isequal(sym,double)` no longer works in the current Symbolic Math
  Toolbox) in `testRegion.m`'s `testCreation`.

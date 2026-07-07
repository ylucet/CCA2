# Session Handoff

_Last updated: 2026-07-07 (local session end)_

## What happened this session
Continued the `cplq` conjugate engine (Step 2, `conjPieceCPLQ.m`) for the CCA2 MATLAB toolbox.
Added the zero-convex-edge bilinear conjugate case: `conjBilinearXYzeroCE` computes the conjugate
of `f = x*y` over a triangle with no convex edge (sup attained at a vertex → 3-cone piecewise-linear
QuaPar). Refactored the shared vertex-max geometry out of `conjLinearTriangle` into a reusable
`conjVertexMax(V, ellv)` helper. Replaced the old stub test with a real numeric+analytic
verification test, and added a new not-implemented test pinning the still-open two-convex-edge case.
Verified full suite 99/99 on Frances (MATLAB R2024a), committed (`9e35f09`) and pushed. Also added
`CLAUDE.md` + this handoff file (`624f461`, pushed). Started nCE==2 (two convex edges, COAP
Appendix B.3) but only got as far as locating the reference implementation
(`/home/ylucet/FLTW/codeOld/cPLQ/conjugateExpr.m`, `plq_1p.m` has `conjugateFunction` around line
~540+ per earlier notes) before being interrupted — no code written yet for nCE==2.

## Where things stand
- Branch: `cplq-engine` @ `624f461` — "Add CLAUDE.md wiring and session handoff notes"
- Pushed: yes (working tree clean, local == origin/cplq-engine)

## Next steps
- **In progress**: nCE==2 (two convex edges, COAP Appendix B.3) in `conjPieceCPLQ.m`. Next action
  is to read `/home/ylucet/FLTW/codeOld/cPLQ/conjugateExpr.m` and the `conjugateFunction` method in
  `plq_1p.m` (~line 540+) for the closed-form structure, then port it the same way
  `conjBilinearXYoneCE` ported the nCE==1 case (COAP B.2) — derive closed-form formulas (not
  symbolic solve), build the QuaPar subdivision, verify by numeric sup-sampling on a fine grid
  before trusting the analytic form, then wire into the `conjPieceCPLQ` dispatch (currently in the
  `else` branch that errors with `notImplemented` for `size(ce,1) > 1`). Test triangle already
  chosen and pinned as a not-implemented case in `conjPieceCPLQTest.m`
  (`bilinearTwoConvexEdgesNotImplemented`, triangle `(0,0),(2,1),(1,2)`) — reuse it, but convert it
  from a `verifyError` check into a real numeric-sup verification once implemented.
- Then nCE==3 (three-convex-edge split).
- Then general indefinite bilinear quadratics (scale/rotate via `bilinearFrame` + affine
  substitution back to the general case, not just pure `x*y`).
- Then conjugate of rational pieces (parabolic edges — curved-edge `orderEdges` support already
  landed in commit `ccbaba6`, so this is unblocked).
- Finally, Step 3: max of per-piece conjugates → `conj(QuaPoly)` for nonconvex inputs, wiring the
  full `conjCPLQ` pipeline (Step1 `convEnvCPLQ` → RatPol pieces; Step2 `conjPieceCPLQ` each →
  QuaPar; Step3 max).
- Test command (run on Frances, MATLAB path is fiddly re: Maple/Symbolic toolbox conflict):
  ```
  matlab -batch "restoredefaultpath; rehash toolboxcache; cd('/home/ylucet/CCA2/CCA2'); \
    res=runtests({'RatPolTest','convEnvCPLQTest','QuaParTest','conjCPLQTest','conjPieceCPLQTest','PLQVCTest'}); \
    fprintf('\n==== SUMMARY ====\nTOTAL %d PASSED %d FAILED %d\n',numel(res),sum([res.Passed]),sum([res.Failed])); \
    exit(sum([res.Failed])>0)"
  ```
- Do not commit or push without asking first — the author has said not to commit unless
  explicitly asked. Prior commit+push actions were one-off explicit requests, not a standing
  authorization; ask again each time.

## Relevant files
- `conjPieceCPLQ.m` — Step 2 of the pipeline (this session's main file); dispatch for bilinear
  pieces is inside the main `conjPieceCPLQ` function, `else` branch near the bottom errors on
  `size(ce,1) > 1` (nCE==2/3) — that's where the nCE==2 branch needs to be added.
- `conjPieceCPLQTest.m` — its tests; `bilinearTwoConvexEdgesNotImplemented` is the pinned nCE==2
  test triangle to convert once implemented.
- `convEnvCPLQ.m` / `convEnvCPLQTest.m` — Step 1 (per-piece convex envelope), already complete.
- `QuaPar.m`, `RatPol.m`, `QuaPoly.m` — the class hierarchy the pipeline builds on.
- `/home/ylucet/FLTW/codeOld/cPLQ/plq_1p.m` (`conjugateFunction`, ~line 540+) and
  `conjugateExpr.m` — unintegrated research code with the nCE==2 reference formulas.
- `/home/ylucet/CCA2/DESIGN.md` — design proposal (outside this repo).
- `/home/ylucet/CCA2/s10589-026-00781-5.pdf` (COAP, convex envelope) and
  `/home/ylucet/CCA2/s10898-025-01503-7.pdf` (JOGO, conjugate) — authoritative spec, esp. Appendix
  A (convex envelope, pp. 25-34) and Appendix B (conjugate, B.3 for nCE==2) for the remaining cases.

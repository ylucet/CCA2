# Session Handoff

_Last updated: 2026-07-07T00:00:00 (local session end)_

## What happened this session
Continued the `cplq` conjugate engine (Step 2, `conjPieceCPLQ.m`) for the CCA2 MATLAB toolbox.
Added the zero-convex-edge bilinear conjugate case: `conjBilinearXYzeroCE` computes the conjugate
of `f = x*y` over a triangle with no convex edge (sup attained at a vertex → 3-cone piecewise-linear
QuaPar). Refactored the shared vertex-max geometry out of `conjLinearTriangle` into a reusable
`conjVertexMax(V, ellv)` helper. Replaced the old stub test with a real numeric+analytic
verification test, and added a new not-implemented test pinning the still-open two-convex-edge case.
Verified full suite 99/99 on Frances (MATLAB R2024a), committed (`9e35f09`), and pushed to
`origin/cplq-engine`.

## Where things stand
- Branch: `cplq-engine` @ `9e35f09` — "conjPieceCPLQ: add zero-convex-edge bilinear xy conjugate"
- Pushed: yes

## Next steps
- Implement nCE==2 (two convex edges, COAP Appendix B.3) in `conjPieceCPLQ.m`.
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
  explicitly asked. This session's commit+push was a one-off explicit request, not a standing
  authorization.

## Relevant files
- `conjPieceCPLQ.m` — Step 2 of the pipeline (this session's main file).
- `conjPieceCPLQTest.m` — its tests.
- `convEnvCPLQ.m` / `convEnvCPLQTest.m` — Step 1 (per-piece convex envelope), already complete.
- `QuaPar.m`, `RatPol.m`, `QuaPoly.m` — the class hierarchy the pipeline builds on.
- `/home/ylucet/CCA2/DESIGN.md` — design proposal (outside this repo).
- `/home/ylucet/CCA2/s10589-026-00781-5.pdf` (COAP, convex envelope) and
  `/home/ylucet/CCA2/s10898-025-01503-7.pdf` (JOGO, conjugate) — authoritative spec, esp. Appendix
  A (convex envelope, pp. 25-34) and Appendix B (conjugate) for the remaining cases.

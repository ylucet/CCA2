# Session Handoff

_Last updated: 2026-07-13T00:00:00-07:00_

## What happened this session

Reviewed the `cplq-engine` branch (35 commits ahead of a 2024-06-24 `main`) and, on request,
fast-forward merged it into `main` and deleted `cplq-engine`. Then, on request, audited which
operators/engines from `DESIGN.md` actually exist as MATLAB code (answer: only the `'cplq'`
conjugate engine, plus `QuaPoly`/`QuaPar`/`RatPol`/`maxQuaPar`; `'pqp'`/`'graph'`/`RatPar`/most
derived operators are still just proposal). Updated `DESIGN.md` with a prominent "Implementation
status" section reflecting this, and corrected several authorship/attribution errors (parametric-QP
is Jakee Khan's work, not Deepak Kumar's; added Bryan Gardiner and Tasnuva Haque references;
clarified `infConv`'s convexity requirement and that `moreau` is a single conjugate, not built on
`infConv`). Implemented `scalarMul`/`negate` (all three classes) and `add` for `QuaPoly`
(`addQuaPoly.m`, adapting `maxQuaPar.m`'s polygon-clipping machinery), with a new test file
verified against hand-computed ground truth. While building `add`, found and fixed a real bug:
`QuaPoly.m`/`RatPol.m` had unfused copies of the `orderEdges` pivot-vertex bug fixed in `QuaPar.m`
earlier this session, AND that earlier fix was itself incomplete (it corrupted the closing-angle
convexity check for bounded faces) — fixed identically in all three files, caught by `QuaPoly`'s
own pre-existing `testNonconvex1-4`. Full suite: 123/123 PASS. All work committed and pushed.

## Where things stand

- Branch: `main` (the `cplq-engine` branch no longer exists — merged and deleted this session).
- Pushed: yes, `origin/main` @ `425625c` — "Implement scalarMul, negate, and add (QuaPoly);
  document pqp/graph as unimplemented".

## Next steps

1. **`add` is QuaPoly-only.** Extending it to `QuaPar` (needs curved/`Ec`-edge clipping in
   `clipByFace`/`clipPolyHalfPlane`) or `RatPol` (needs a common-denominator sum, not just adding
   numerators) is open — see `DESIGN.md`'s Implementation status section.
2. **`'pqp'`/`'graph'` conjugate engines are unimplemented** (`conj(f,'pqp'|'graph')` errors
   explicitly, doesn't silently give a wrong answer) — porting Jakee Khan's parametric-QP code or
   Tasnuva Haque's entity-graph algorithm is future work, not started.
3. Still-open items from before this session (untouched): the standalone `RatPol.conj` gap; the
   `delta>0, Delta=0` degeneracy open question (see `/home/ylucet/CCA2/3-edge.tex`'s Conclusion,
   outside this repo).
4. Worth a quick audit: are there other unfused copies of buggy shared code across
   `QuaPoly.m`/`QuaPar.m`/`RatPol.m` beyond `orderEdges`? All three files are largely
   parallel/copy-pasted, so a bug found in one is worth checking in the other two (as happened
   twice this session already).
5. Test command (Frances, prefix every batch call — Maple toolbox conflicts with Symbolic):
   ```
   matlab -batch "restoredefaultpath; rehash toolboxcache; cd('/home/ylucet/CCA2/CCA2'); \
     res=runtests({'RatPolTest','convEnvCPLQTest','QuaParTest','conjCPLQTest','conjPieceCPLQTest','PLQVCTest','maxQuaParTest','addQuaPolyTest'}); \
     fprintf('\n==== SUMMARY ====\nTOTAL %d PASSED %d FAILED %d\n',numel(res),sum([res.Passed]),sum([res.Failed])); \
     exit(sum([res.Failed])>0)"
   ```

## Relevant files

- `DESIGN.md` — new "Implementation status" section (right after the TL;DR) is the authoritative
  "what's real vs proposed" pointer; also has corrected authorship attributions throughout
  (search for [JAKEE-13], [HAQUE-17]/[HAQUE-18], [GARDINER-09/11/13/14], [KUMAR-20],
  [HIRIART-URRUTY-07]) and the `infConv`/`moreau` convexity clarifications in II.6.
- `reference/` — added this session (via git, not by this session's own commits): PDFs for all
  the above references (Jakee Khan's and Tasnuva Haque's MSc theses, Gardiner's four papers, etc.).
- `addQuaPoly.m` / `addQuaPolyTest.m` — new this session; the `QuaPoly.add` implementation.
- `QuaPoly.m`, `QuaPar.m`, `RatPol.m` — each got `scalarMul`/`negate` methods, and each had its
  `orderEdges` fixed (see "What happened" above for the two-part bug).
- `maxQuaPar.m` — untouched this session, but `addQuaPoly.m`'s geometry helpers were adapted from
  it (facePoly/clipByFace/clipPolyHalfPlane/polyConstraints/insertPassthroughVertices/
  assemblePieces-equivalent) — the two files' geometry cores are meant to stay in sync if either
  is debugged further.

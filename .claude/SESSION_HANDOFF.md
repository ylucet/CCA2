# Session Handoff

_Last updated: 2026-07-14T23:59:00Z_

## What happened this session

Diagnosed (did not fix) the prior session's #1-priority "finite but wrong `g.eval`" bug. Traced it
to `convEnvCPLQ.m`'s 2-convex-edge quadratic (`envelopeFromClassified` case 2 / `twoEdgeQuadPlain` /
`buildTwoEdge`, implementing [COAP] Appendix A.4) not always being the true tightest convex
envelope, even though it is a valid minorant that correctly touches `f=u1u2` along both classified
convex edges. This is a much deeper finding than expected: it reproduces on **the paper's own
Appendix A.4.3 worked example** (`V=(2,1),(0,0),(1,0)`), confirmed via three independent methods
(hand-derived exact boundary analysis, the real `convEnvCPLQ`+`conjPieceCPLQ` pipeline, and an
independent biconjugate reconstruction `conv(f)(x0)=sup_s[s.x0 - supBilinearOverPoly(s,V)]`). Root
cause and the true envelope's boundary conditions along the triangle's third ("weak") edge are now
characterized, but a correct general fix needs an actual extension of the published derivation
(likely a genuine sub-partition of the triangle, not a simple pointwise combination) — this is
research-level work, not a code bug fix, and was intentionally left open this session per the
user's direction. Full derivation trail is in `DESIGN.md`.

## Where things stand

- Branch: `main` @ `b69dfec` — "docs: diagnose deep Appendix A.4 envelope-tightness gap in
  convEnvCPLQ".
- Pushed: pending (awaiting user confirmation).
- No code changed this session — only `DESIGN.md` (documentation of the finding). Full test suite
  therefore still 145/145 PASS (unchanged from end of prior session).

## Next steps

- **Top priority, and now known to be genuinely hard**: derive the correct general convex envelope
  for a triangle with exactly 2 convex edges (COAP Appendix A.4), covering the case where the
  paper's own "checking the bounds we get the domain as the entire triangle" check (stated but not
  proven general, and never implemented as an actual runtime check) fails. Established boundary
  conditions to build from:
  - Along the 2 classified convex edges: the envelope must equal `f` exactly (already correct in
    the current `twoEdgeQuadPlain` formula).
  - Along the 3rd ("weak", non-convex) edge: the true envelope equals the AFFINE CHORD between its
    2 endpoint `f`-values (verified to 6 decimals via biconjugate reconstruction, at multiple
    points, on 2 different triangles).
  - A naive `max(q1, affine-interpolation-of-all-3-vertex-values)` is NOT a valid general fix: the
    affine piece exceeds `f` (and `q1`) in a region nearer the common vertex (by up to 1.46 on the
    stress-test triangle), and a systematic re-check over random dual points `s` found the
    combination's conjugate wrong in the OPPOSITE direction too (undershooting truth by >1). A
    correct fix likely needs a genuine second piece/sub-partition (analogous to how the
    3-convex-edge case already splits into two sub-triangles), not a simple max.
  - Repro triangles or points for verification once a candidate fix exists (compare against
    `maxQuaParTest.supBilinearOverPoly`, the exact/trustworthy ground truth):
    - Paper's own example: `V=(2,1),(0,0),(1,0)`, bad `s=(-0.008727,-0.999962)` (exact truth `0`, current
      code gives `0.03864091`); the "weak-edge dip" point is `x0=(0.474343,0)`, `q1(x0)=-0.042780`
      vs true `conv(f)(x0)≈0`.
    - Stress-test triangle: `T=(3.8398,5.0413),(8.8152,7.2338),(8.7969,5.6447)`, bad
      `s=(6.457616,8.384129)` (exact truth `54.477034`, current code gives `54.532567`); dip point
      `x0=(8.360846,5.984837)` (on the weak/seam edge, `t≈0.3367` of the way from V2 to V3),
      `q1(x0)=49.636211` vs true `conv(f)(x0)=49.746078` (matches the affine chord there to 6
      decimals).
  - This session's scratch MATLAB probes (not committed — see
    `C:\Users\ylucet\AppData\Local\Temp\claude\...\scratchpad\`, session-specific temp dir, likely
    gone by next session) contain the biconjugate-reconstruction technique
    (`conv(f)(x0)=sup_s[s.x0-supBilinearOverPoly(s,V)]`, solved via a coarse-then-refined 2D grid
    search since `fminsearch`/`fminunc` were numerically unreliable here) — reuse this technique to
    probe candidate fixes against ground truth at arbitrary points before trusting any closed form.
  - Once a fix is found: add a regression test in `convEnvCPLQTest.m` (or a new dedicated test)
    using the paper's own example at the bad `s` above, since the existing
    `bilinearTwoConvexEdgesQuadratic` test only checks the primal formula's VALUE at a few interior
    points against the published closed form — it does NOT check conjugate/tightness, which is
    exactly why this bug went undetected until now.
- Related, still unconfirmed whether same root cause: `QuaPar.orderEdges`/`createP` rejects the
  face topology for some near-degenerate/thin triangles with a DIFFERENT error ("Face k has ..." /
  "expected 2 but got N").
- Lower priority / untouched: `partialConj` for `'cplq'`/`'pqp'`; `add` for `RatPol` and the
  `RatPar` parent class; conjugate engines `'pqp'`/`'graph'` (unimplemented, error explicitly); the
  standalone `RatPol.conj` gap.

## Relevant files

- `DESIGN.md` — `maxQuaPar.m`'s Implementation-status bullet has the full new-finding paragraph
  ("Open research question found while diagnosing the small-magnitude wrong-answer cases"),
  including the exact repro, the 3-method verification, and why the naive `max(q1,affine)` fix
  doesn't work.
- `convEnvCPLQ.m` — `envelopeFromClassified` case 2 / `twoEdgeQuadPlain` / `buildTwoEdge`: the
  functions implementing [COAP] Appendix A.4, now known to need an additional piece/case-split.
  NOT modified this session (diagnosis only).
- `convEnvCPLQTest.m` — `bilinearTwoConvexEdgesQuadratic` (line ~78): the existing test using the
  paper's own example, which only checks the primal formula's value (not tightness/conjugate) — the
  gap this reveals is exactly why the bug was invisible until this session.
- `maxQuaParTest.m` — `supBilinearOverPoly` is the exact, trustworthy ground-truth function used
  throughout this diagnosis; `buildG1G2ForTriangle` is the real end-to-end pipeline entry point used
  to reproduce the bug via `maxQuaPar`.
- `QuaPar.m` — NOT modified this session; still a candidate location for a DIFFERENT, unconfirmed
  issue (the `orderEdges`/`createP` face-topology rejection above).

# Session Handoff

_Last updated: 2026-07-14T23:10:00Z_

## What happened this session

Fixed the prior session's #1-priority residual `maxQuaPar:internal` crash ("no matching
neighbour", ~1/800 valid random triangles). Diagnosed it first: the orphaned half-edge is always
part of a genuinely ambiguous 3-way (or more) vertex cluster, where several pieces meeting near
one point each compute a slightly different position (order ~1e-4 — too coarse for
`matchHalfEdges`' cross-arithmetic-noise tolerance, too fine to be a distinct feature) for what is
mathematically ONE vertex, so the best-first greedy matcher can pair up only 2 of the 3-plus
mutually-close half-edges and always orphans one. The prior handoff guessed this needed full
vertex PROVENANCE (tagging each vertex with which original `g1`/`g2` face boundary produced it);
that turned out to be unnecessary. Diagnosis (via targeted instrumentation, since run and removed)
showed the orphaned edge's own two endpoints ALWAYS resolve to the very same global vertex once
the OTHER, already-confirmed matches on its own piece's boundary are accounted for — i.e. the
orphaned edge is provably zero-length in the resolved geometry, so dropping it (emitting no edge
for it) is exactly correct, not a guess. Fixed via a new `checkOrphanHalfEdges` function, called
after global vertex identity is built: it drops an orphan only when its two endpoints already
coincide globally; otherwise it still raises the original error (a genuinely unresolved gap, or an
orphaned ray — rays keep the strict old behaviour, no evidence yet they can be legitimately
degenerate this way). See `checkOrphanHalfEdges`'s own header in `maxQuaPar.m`, the new HISTORY
entry in `assemblePieces`, and `maxQuaParTest.checkOrphanHalfEdgesDropsProvablyDegenerateOrphanEdges`
(regression test using 5 of the 6 randomly-found repro triangles, checked end-to-end against
ground truth).

**Important — a second, more serious issue was found while verifying this fix, and is NOT fixed:**
one of the 6 repro triangles (a 5-piece cluster, more complex than the simple 3-way case) no longer
crashes but instead SILENTLY returns the wrong answer (`g.eval` returns `Inf`, i.e. a genuine
domain-coverage gap) over a substantial region of the plane. Confirmed this is a real bug (not a
borderline-precision artifact): checked via `polyConstraints` that the original, pre-assembly piece
geometrically DOES cover the missing region with a comfortable margin (not just barely inside), so
the piece-generation loop in `maxQuaPar.m` (`clipByFace`/`splitCell`) is fine — the bug is somewhere
in final assembly or in `QuaPar.eval`'s per-edge membership test. Likely candidate (not confirmed):
a sign/orientation issue in `QuaPar.eval`'s `edgeConics`/`evalConic`-based membership test
(`vals(:,t) = evalConic(...)*sign(Pe(t)); flags = all(vals<=0,2)`, with NO tolerance margin)
specific to an unbounded face whose two rays point in nearly the SAME direction (a degenerate
"parallel-strip" shape rather than a normal wedge) — this repro's face 2 has exactly that shape.
This bug is NOT caused by this session's `checkOrphanHalfEdges` fix — it pre-dates it and was
simply masked by the crash for this specific triangle (execution never reached final assembly
before). A broader correctness-focused stress test (not just crash-counting) found this
silent-wrong-answer pattern in **~2/300 valid random triangles** — the same order of magnitude as
the crash rate this session's fix addressed, and worse in kind (silent, not loud). **This should be
the top priority for the next session** — see "Repro for the still-open issue" below.

## Where things stand

- Branch: `main` @ `bcc2851` — "Fix residual maxQuaPar:internal crash by dropping
  provably-degenerate orphan edges".
- Pushed: pending.
- Full test suite: 144/144 PASS (143 pre-existing + 1 new regression test,
  `checkOrphanHalfEdgesDropsProvablyDegenerateOrphanEdges`). Working tree clean (this session's
  work is committed).
- Targeted stress test (3000 random triangles, same generator as before): the specific
  "no matching neighbour" crash is now **0/455 valid samples** (down from ~4-6/800 before this
  session's fix). 5 unrelated `orderEdges`-topology-rejection errors also appeared in that run —
  pre-existing, already known (see "Related" below), unaffected by this session.

## Next steps

- **Highest priority (new, found this session)**: the silent domain-coverage-gap bug described
  above. Repro triangle: `T=(7.8665,4.6784),(2.6908,1.9477),(0.3892,0.7130)` via
  `maxQuaParTest.buildG1G2ForTriangle(T)` then `maxQuaPar(g1,g2)` — e.g.
  `g.eval([3.380265 -0.644943])` returns `Inf` but the true value is `2.69462`
  (`maxQuaParTest.supBilinearOverPoly`). Diagnostic trail (not left in the code — re-derive via
  temporary `fprintf`/`getenv` instrumentation the same way this session did, in
  `matchHalfEdges`/`checkOrphanHalfEdges`/`buildFinalEdgesAndFaces`, or via a standalone script
  calling `maxQuaParTest.buildG1G2ForTriangle` + inspecting `g.P{2}`/`g.E`/`g.F` directly):
  1. `g1.eval(s)`/`g2.eval(s)` both correctly return finite values at the failing point (each
     INPUT arrangement has no gap) — the bug is specific to the ASSEMBLED max.
  2. The original piece (from `clipByFace(facePoly(g1,1), facePoly(g2,4))`, an unbounded
     3-real-vertex cell) DOES geometrically contain the failing region
     (`polyConstraints` margin ~1e-4 to 1e-5, comfortably inside, not borderline).
  3. The final `g.P{2}` reconstructs what LOOKS combinatorially correct on manual inspection
     (ray-in at the right apex, two segments, ray-out at the right apex, matching the piece's own
     true CCW order once the CW-to-CCW reversal is accounted for) — so the bug is suspected to be
     in `QuaPar.eval`'s membership test itself (or possibly `orderEdges`' sign/`P` construction,
     not yet ruled out), not in `maxQuaPar.m`'s assembly.
  4. This repro's face 2 is unusual: its two rays (`dirIn`/`dirOut`) point in NEARLY THE SAME
     direction (a thin, skewed "infinite strip" rather than a normal diverging wedge) — a
     plausible trigger for a sign bug not exercised by the existing test suite's ray-containing
     cases (which all have more normally-diverging wedges).
  5. Suggested next step: write a minimal, direct `QuaPar` (V/E/F/P by hand, 2 faces, one an
     unbounded double-same-direction-ray strip) to isolate whether the bug is in `eval` or in how
     `orderEdges`/`createP` build `P`/sign conventions for this specific shape, independent of the
     whole `maxQuaPar` pipeline.
- Once the above is fixed, re-run the correctness-focused stress test in this session's scratch
  work (not committed — regenerate similarly: sample random triangles, call `maxQuaPar`, spot-check
  `g.eval` at several random points against `maxQuaParTest.supBilinearOverPoly`) to confirm the
  ~2/300 silent-wrong-answer rate drops to 0.
- Related (unchanged this session): `QuaPar.orderEdges`/`createP` rejects the face topology for
  some near-degenerate/thin triangles with a DIFFERENT error ("Face k has ..." / "expected 2 but
  got N") — still not confirmed whether it shares a root cause with either issue above. Appeared
  5/455 times in this session's stress test, consistent with the prior session's ~4% figure.
- Lower priority / untouched this session: `partialConj` for `'cplq'`/`'pqp'`; `add` for `RatPol`
  and the `RatPar` parent class; conjugate engines `'pqp'`/`'graph'` (unimplemented, error
  explicitly); the standalone `RatPol.conj` gap.

## Relevant files

- `maxQuaPar.m` — new `checkOrphanHalfEdges` function (this session's fix), called from
  `assemblePieces` right after `buildGlobalVertices`; `matchHalfEdges` no longer errors on an
  unmatched half-edge itself (moved to `checkOrphanHalfEdges`). Full derivation in
  `assemblePieces`' HISTORY comment (new entry) and `checkOrphanHalfEdges`'s own header.
- `maxQuaParTest.m` — `checkOrphanHalfEdgesDropsProvablyDegenerateOrphanEdges`, new regression test
  covering 5 of the 6 repro triangles found this session (all near-degenerate/thin, e.g.
  `T=(8.5697,2.6142),(5.0151,1.8051),(1.3296,0.9185)`).
- `DESIGN.md` — `maxQuaPar.m`'s Implementation-status bullet has a new "Reliability fix
  (2026-07-14, later session)" paragraph for this session's fix, plus a paragraph documenting the
  new still-open silent-wrong-answer issue.
- `QuaPar.m` — NOT modified this session, but `eval` (and possibly `orderEdges`/`createP`) is the
  prime suspect for the new still-open issue above; not yet touched.
- `convEnvCPLQ.m` — `splitThreeConvex` (prior-prior session's correctness fix; UNCHANGED this
  session).

# Session Handoff

_Last updated: 2026-07-17T23:40:00Z_

## What happened this session

Attempted to fix the 2-convex-edge tightness bug diagnosed (but not fixed) last session, using
the counterexample triangle `T=(-9.95506,3.70366),(-9.345,-5.34552),(1.29049,5.31738)`. Result: a
real, verified partial fix (Part 1), plus a much more precise diagnosis of a second, deeper,
still-open gap (Part 2) that the fix does not touch. Full derivation trail is in `DESIGN.md`'s
`convEnvCPLQ.m` entry (search "Diagnosis (2026-07-17, later session)").

**Part 1 (fixed, verified, kept)**: `splitTwoConvexEdges`'s "other" sub-triangle used a bespoke
quadratic (`buildEdgeAffinePiece`) that was never actually the tight envelope, just a
valid-on-the-boundary guess (proved via a real residual, ~0.36 max, fitting a general quadratic to
independent ground truth). Fix: that sub-triangle's own edges always reclassify (checked on ~17
random triangles) as an ordinary ONE-convex-edge triangle, so it should reuse the existing,
already-validated Appendix A.3 rational formula (`envelopeFromClassified` case 1) directly, no
bespoke construction -- the same pattern the 3-convex-edge split already uses, now unified in
`convEnvCPLQ.m`. Verified against two independent ground-truth methods (numerically-maximized
biconjugate `sup_s[s.x0-supBilinearOverPoly(s,V)]`, and the 3D convex hull of the graph
`(x,y,xy)`) to grid resolution on the counterexample and 16/17 fresh random triangles. Consequence:
the "other" piece is now genuinely rational, which `conjPieceCPLQ` cannot conjugate yet (a
pre-existing, separately-tracked gap) -- `maxQuaParTest.m`'s `extractTriFace` used to silently drop
the denominator there, building the wrong `QuaPoly`; now errors clearly instead
(`extractTriFace:rationalFaceNotSupported`), and the one affected regression test was updated to
expect that clear error. Full MATLAB suite: 79/79, no regressions.

**Part 2 (found this session, NOT fixed, more fundamental)**: fixing Part 1 does not make the
counterexample's wrong answer go away, because `q1` (the formula kept UNCHANGED for the
sub-triangle containing shared vertex `P`) is itself not tight throughout that sub-triangle --
confirmed via ground truth, undershooting by up to ~2.55 along the internal seam's interior. A
1500-triangle `rng(2026)` resample found 508 genuine-split 2-convex-edge triangles, 231 (~45%)
still showing a real gap, IDENTICALLY before and after Part 1's fix -- Part 1 is real but isolated,
Part 2 drives the aggregate wrong-answer rate. **Important correction made this session**: initial
numerical exploration (tracing where `q1` stops matching truth) fit a general conic far better
than any straight line, and was reported as possible evidence the true domain partition might need
a curved boundary. The user correctly pushed back: **Locatelli proved the convex envelope of a
quadratic over a polyhedron always admits a polyhedral subdivision** -- so that curve-fit reflects
an incomplete search (the wrong split/pairing), not a counterexample to the theorem. Recorded as an
established fact in `DESIGN.md` (new `[LOCATELLI]` reference entry, exact citation not yet filled
in) specifically so the next session does not re-derive this and does not pursue a curved-`RatPol`
extension. The correct polyhedral split has not yet been found; several specific attempts were
ruled out this session (see `DESIGN.md` for details: naive max of two candidates, a 3-cevian split
via a third point on the weak edge, an unconstrained general-rational fit) -- kept as negative
results so they aren't retried.

## Where things stand

- Branch: `main` @ `acef5b1` -- "Fix convEnvCPLQ 2-convex-edge split: other sub-triangle is a
  genuine 1CE piece". Contains this session's `convEnvCPLQ.m` fix, `maxQuaParTest.m` update, and
  `DESIGN.md`'s full diagnosis (including the Locatelli correction).
- Pushed: pending (see below).
- Full MATLAB suite: 79/79 passing (`convEnvCPLQTest`, `conjPieceCPLQTest`, `conjCPLQTest`,
  `maxQuaParTest`, `RatPolTest`, `QuaParTest`).
- Untracked, not part of this repo: `cPLQ/` (clone of the original reference implementation, see
  prior handoff); not touched this session.
- Scratch diagnostic scripts (MATLAB + Python) live under the OS-temp scratchpad
  (session-specific, not persisted) -- regenerable from `DESIGN.md`'s description if needed again,
  not guaranteed to survive to a future session. Note: use
  `C:\Users\ylucet\AppData\Local\miniforge3\python.exe` for `numpy`/`scipy`/`sympy` -- the bare
  `python3` on PATH is a different interpreter without them.

## Next steps

- **Top priority**: find the correct POLYHEDRAL (straight-edge) split for Part 2. Per Locatelli
  one is guaranteed to exist; `q1` + a single Appendix A.3 rational piece (any placement/pairing
  tried this session) is simply the wrong combinatorial split. Likely directions: more than 2
  pieces; a different choice of which edge/vertex anchors each piece; a genuinely different split
  geometry than a single cevian from a weak-edge endpoint or from `P`. Do NOT pursue a curved
  `RatPol` domain extension -- that was a mistaken interim conclusion, corrected this session.
- Decide whether/how to unblock `conjPieceCPLQ` for a genuinely rational (quadratic/linear) piece
  -- Part 1's fix makes this reachable via the ordinary single-triangle 2-convex-edge path (not
  just the previously-known multi-face path), so it's now a more common need. Not attempted this
  session.
- Fill in the exact `[LOCATELLI]` citation in `DESIGN.md`'s reference list.
- Still open, unchanged from before: the 2/741 (0.3%) residual `maxQuaPar:internal` crashes (see
  prior handoff for repro triangles) -- not investigated this session.
- Lower priority / untouched: `QuaPar.orderEdges`/`createP`'s "Face k has ... expected 2 but got
  N" error on near-degenerate/thin triangles; `partialConj` for `'cplq'`/`'pqp'`; `add` for
  `RatPol`/`RatPar`; conjugate engines `'pqp'`/`'graph'`; the standalone `RatPol.conj` gap.

## Relevant files

- `convEnvCPLQ.m` -- `splitTwoConvexEdges` (Part 1 fix: the "other" sub-triangle now uses
  `envelopeFromClassified` on its own reclassified edges, not `buildEdgeAffinePiece`'s result;
  `buildEdgeAffinePiece` kept, now only as a seam-finding aid). `nCE==2` branch unified with the
  `nCE==3` branch's per-sub-triangle `envelopeFromClassified` loop. Part 2 is NOT fixed here.
- `maxQuaParTest.m` -- `extractTriFace` now errors clearly on a rational face instead of silently
  building a wrong `QuaPoly`; `matchHalfEdgesRejectsSameSideRayPairingAndDropsSubsumedPieces`
  updated to expect that error for `V=[2 1;0 0;1 0]`.
- `DESIGN.md` -- full derivation trail for both Part 1 and Part 2, including the Locatelli
  correction; search "Diagnosis (2026-07-17, later session)" under `convEnvCPLQ.m`'s entry, and
  the new `[LOCATELLI]` reference entry near the top.
- `doc/bug.svg`, `doc/bug_pipeline.svg` -- prior session's figures for the same counterexample
  triangle; still accurate background, not updated this session.
- `cPLQ/` -- untracked clone of the original reference implementation; not touched this session.

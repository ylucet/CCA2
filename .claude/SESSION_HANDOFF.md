# Session Handoff

_Last updated: 2026-07-14T23:55:00Z_

## What happened this session

Fixed the prior session's #1-priority silent-wrong-answer bug (`g.eval` returning `Inf` over a
real, comfortably-covered region of the plane). Root cause: `insertPassthroughVertices` in
`maxQuaPar.m` shared its "already a vertex" check with `onOpenSegment`/`onOpenRay`'s own
edge-matching tolerance (1e-7) — tighter than the ~3e-5 cross-arithmetic-noise floor between an
original face vertex and the already-present cell vertex it geometrically coincides with. A point
that geometrically WAS an already-present vertex narrowly failed that check and got inserted as a
brand-new, near-duplicate vertex instead, creating a near-zero-length sliver edge whose line
equation was dominated by floating-point noise in the tiny direction vector — `QuaPar.eval`'s exact
(no-tolerance) membership test then wrongly excluded a real region behind that noise-direction
line. Fixed by decoupling the "already a vertex" pre-check into its own, wider tolerance
(`tolSnap=1e-4`), leaving the edge-matching tolerance untouched. Widening the SHARED tolerance was
tried first and reverted: it fixed the target bug but broke two other regression tests at different
noise scales (a genuinely distinct ~7.9e-4-away vertex wrongly absorbed into the wrong edge in one;
a genuine ~1e-4 feature of a near-degenerate triangle wrongly merged in the other) — no single
shared value separated all three cases, but `tolSnap` only ever recognizes a point as coincident
with a vertex the cell construction already produced, so it can't manufacture that kind of wrong
topology. See `insertPassthroughVertices`'s own header in `maxQuaPar.m`,
`maxQuaParTest.insertPassthroughVerticesDropsNearDuplicateCrossingPoint`, and the new DESIGN.md
paragraph for the full derivation.

## Where things stand

- Branch: `main` @ `5133c17` — "Fix silent maxQuaPar wrong-answer bug from unmerged near-duplicate
  passthrough vertices".
- Pushed: pending (awaiting user confirmation).
- Full test suite: 145/145 PASS (144 pre-existing + 1 new regression test,
  `insertPassthroughVerticesDropsNearDuplicateCrossingPoint`). Working tree clean.
- 3000-triangle correctness stress test (sample random triangles, call `maxQuaPar`, spot-check
  `g.eval` at several points against `maxQuaParTest.supBilinearOverPoly`): 0/491 valid samples gave
  the silent-`Inf` domain-gap pattern (down from ~2/300 documented before this session's fix).

## Next steps

- **New issue found this session while stress-testing the fix (not yet diagnosed)**: a handful of
  valid triangles produce a FINITE but WRONG `g.eval` value (not `Inf` — a different failure mode
  from the bug just fixed), with small relative error (~0.1%-2%). Confirmed these are PRE-EXISTING
  (reproduce identically on the code from before this session's fix, so unrelated to today's
  change) — likely a genuinely separate, still-open bug. Repro triangles/points found this session
  (from `maxQuaParTest.buildG1G2ForTriangle(T)` then `maxQuaPar(g1,g2)`, evaluated at `s`, compared
  to `maxQuaParTest.supBilinearOverPoly(s,T)`):
  - `T=(3.8398,5.0413),(8.8152,7.2338),(8.7969,5.6447)`, `s=(6.457616,8.384129)`: actual=54.532567,
    expected=54.477034.
  - `T=(2.4884,8.4297),(1.3254,1.3537),(8.2373,8.4915)`, `s=(7.687732,4.914175)`: actual=40.267523,
    expected=39.578708 (largest diff seen, ~0.69).
  - 5 more found in the same stress run (not narrowed down further) — see this session's scratch
    `stress.m`/`check_one.m` (not committed) for the generator (`rng` seeds 999, 555555) to
    reproduce more.
  Not yet isolated to a specific mechanism; a plausible place to start is a point that lands very
  close to a REAL (non-degenerate) face boundary and gets assigned to a face whose coefficients are
  slightly the wrong one, or an issue in `orderEdges`/`createP`'s `P` construction for some face
  shape — not confirmed. **Should be the next priority.**
- Related (unchanged this session): `QuaPar.orderEdges`/`createP` rejects the face topology for
  some near-degenerate/thin triangles with a DIFFERENT error ("Face k has ..." / "expected 2 but
  got N") — still not confirmed whether it shares a root cause with the issue above. One instance
  appeared in this session's 3000-triangle stress test.
- Lower priority / untouched this session: `partialConj` for `'cplq'`/`'pqp'`; `add` for `RatPol`
  and the `RatPar` parent class; conjugate engines `'pqp'`/`'graph'` (unimplemented, error
  explicitly); the standalone `RatPol.conj` gap.

## Relevant files

- `maxQuaPar.m` — `insertPassthroughVertices` now uses a separate `tolSnap=1e-4` for its "already a
  vertex" pre-check, decoupled from `onOpenSegment`/`onOpenRay`'s own `tol=1e-7` matching tolerance.
  Full derivation in the function's own header comment.
- `maxQuaParTest.m` — new regression test
  `insertPassthroughVerticesDropsNearDuplicateCrossingPoint`, using the exact repro triangle
  `T=(7.8665,4.6784),(2.6908,1.9477),(0.3892,0.7130)` from the prior handoff.
- `DESIGN.md` — `maxQuaPar.m`'s Implementation-status bullet has a new "Correctness fix
  (2026-07-14, later session)" paragraph for this session's fix, including why the shared-tolerance
  approach was tried and reverted, and noting the new still-open small-error issue above.
- `QuaPar.m` — NOT modified this session; still a candidate location for the new still-open
  small-error issue (unconfirmed).
- `convEnvCPLQ.m` — `splitThreeConvex` (prior-prior session's correctness fix; UNCHANGED this
  session).

# Session Handoff

_Last updated: 2026-07-18T03:00:00Z_

## What happened this session

Resolved the decision left open by the prior session's WIP commit (`e6e936a`): whether to keep
the recursive `nCE==3` generalization of the Part 2c tightness fix. **Resolved: KEEP it — it was
never actually an open judgment call.** `tangentCevian`'s tightness criterion is proven to depend
only on a 2-convex-edge triangle's own `mh,mw`, independent of how it arose, so a
`splitThreeConvex` sub-triangle can't be exempt from it; reverting would mean knowingly
reinstating a formula already proven untight. Corrected DESIGN.md/`conjCPLQ.m`/`conjPieceCPLQ.m`
to remove the false "undecided" framing and the now-false "3CE already works end to end" claims
(keeping the fix means the 3CE single-triangle case now also hits the pre-existing
`conjPieceCPLQ` rational-piece gap). Also corrected that gap's own framing: the reference `cPLQ/`
package (untracked, opened and read this session) already has a complete symbolic recipe for it
(`plq_1piece.m`'s "type 1" branch) — the remaining work is deriving the closed-form numeric
formula matching every other `conjPieceCPLQ` case, not open research.

Fixed all 8 tests broken by keeping the fix (the original 6, plus 2 more found by re-running the
full suite: `conjPieceCPLQTest.pickRepFindsThinEdgeStripFace`,
`convEnvCPLQTest.threeConvexEdgesSplit`) — mostly by freezing `g1,g2` `QuaPar` fixtures from the
pre-Part-2c pipeline or updating hardcoded piece-count assertions. While diagnosing the last
failure, found and fixed a real bug in `RatPol.eval`: a genuinely rational piece has a removable
(0/0) singularity at its own anchor vertex, and the eval loop let that NaN silently clobber an
already-correct finite value from another matching face at a shared vertex. Full suite: 147/147
passing. Committed as `083de95`.

## Where things stand

- Branch: `main` @ `083de95` — "Keep Part 2c nCE==3 fix (decision resolved), fix all 8 broken
  tests, fix RatPol.eval NaN bug".
- Pushed: pending (asking user this session).

## Next steps

- **Highest value**: derive `conjPieceCPLQ`'s closed-form rational-piece conjugate formula from
  the `cPLQ/plq_1piece.m` recipe (polynomial-divide num/den, parametrize a dual point by scalar
  `t`, eliminate `t` to get the parabola equation, build vertex/edge subdifferential regions —
  see `conjPieceCPLQ.m`'s own TODO comment for the full recipe). This is the single remaining gap
  confining `infConv`/`moreau`/`lasryLions`/`proxAverage` to full-domain-quadratic inputs.
- Lower priority, unchanged from before: exact `[LOCATELLI]` citation in DESIGN.md; 2/741 (0.3%)
  residual `maxQuaPar:internal` crashes; `QuaPar.orderEdges`/`createP`'s near-degenerate-triangle
  error; `partialConj` for `'cplq'`/`'pqp'`; `add` for `RatPol`/`RatPar`; conjugate engines
  `'pqp'`/`'graph'`; standalone `RatPol.conj` gap; a proper random-sample stress test of how often
  a 3-convex-edge triangle actually needs the further split (currently a 9/9 non-random sample —
  a frequency question, separate from whether the fix is correct, which is settled).

## Relevant files

- `convEnvCPLQ.m` — `tangentCevian`/`solveTriangleBF`/`assemblePiecesBF` KEPT, decision resolved,
  unchanged this session (already committed in `e6e936a`).
- `RatPol.m` — `eval`'s main polytope-membership loop fixed: a later face's NaN (removable
  singularity) no longer overwrites an already-found finite value; genuine disagreement still
  correctly produces NaN.
- `conjPieceCPLQ.m` — TODO comment rewritten (no functional change) to correctly point at the
  `cPLQ/plq_1piece.m` recipe instead of describing an open derivation.
- `DESIGN.md` — Part 2c entry rewritten; `maxQuaPar.m`/`conjCPLQ`/"Next planned" bullets updated.
- `conjCPLQTest.m`, `maxQuaParTest.m`, `conjPieceCPLQTest.m`, `convEnvCPLQTest.m` — the 8 fixed
  tests.
- `cPLQ/` — untracked reference clone; opened and read this session, not modified.

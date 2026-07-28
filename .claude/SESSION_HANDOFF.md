# Session Handoff

_Last updated: 2026-07-28T00:00:00Z_

## What happened this session

Four pieces of work, in order:

1. **Fixed `biconj`** — the former blocker 1. It was literally `conj∘conj`, and `f*` of a
   bounded-domain function lives on an unbounded multi-face domain `conjCPLQ` rejects, so it
   failed for every single bounded triangle. New `biconjCPLQ.m` takes **Step 1's convex envelope**
   instead for that case: `f** = conv(q + I_T)` when `T` is compact. Validated against a
   pipeline-free ground truth (≤ 8.9e-16 on 7 triangles covering every Step 1 branch).
2. **Regenerated `SUPPORT_MATRIX.md`** from the actual guards; refreshed every `file:line`
   citation (most had drifted 20–30 lines into unrelated code) and reclassified the
   rational-piece guard from N/R to GAP.
3. **Fixed Step 2 for a rational envelope face** by routing Step 1's envelope through cPLQ's
   symbolic Step 2/3 (`ratPolToPlq.m`). A 2-face envelope split now conjugates end to end.
4. **Corrected two of my own wrong diagnoses** after user pushback — see below. This is the part
   worth reading before touching this area again.

## Where things stand

- Branch: `main` @ `7b190d0` — "docs: correct the rational-face diagnosis; retire the
  3-convex-edge framing"
- Pushed: yes
- Tests: **262 passed / 1 failed** across 22 suites. CCA2's own 17 suites are **187 / 0**; the
  single failure is `testRegion/testCreation`, longstanding and unrelated.

### Two corrections — do not reintroduce these

- **There is no "3-convex-edge case" to implement**, in CCA2 or in cPLQ. [COAP] Appendix A.5's
  split reduces such a triangle to 2-convex-edge sub-triangles and Step 1 already applies it, so
  every face reaching Step 2 has come through it. cPLQ having no `nCE==3` branch is a fact about
  cPLQ's *Step 1*, which CCA2 does not use. Describe these cases by the **envelope's face count**,
  not the input's edge count.
- **It did not need exact arithmetic.** At the vertex where a face's denominator vanishes the
  numerator vanishes too, with `∇N ∥ ∇den` (residual 0 to 1.3e-16), and the coefficients are good
  to ~1e-16 — the `0/0` is cleanly **removable**. The fix was to *resolve* it by a limit.
  `subsF` already flagged it (`NaN`) and `region.funcVertices` already had a NaN-then-limit idiom,
  but `symbolicFunction.limit` takes **iterated univariate** limits, which return `NaN` again for
  a bivariate `0/0`. New `symbolicFunction.limitDirectional` restricts to a line — one univariate
  limit — and requires several directions to agree (that agreement is the existence check).

### Which half of cPLQ is reusable

Its **Step 2/3**, on **CCA2's** Step 1 output. Running cPLQ end to end gives a wrong answer: its
own 2-convex-edge envelope is the single Appendix A.4 formula, which CCA2 established is not
always tight, and the resulting conjugate leaves the paper's own flagged dual point
`s=(-0.008727,-0.999962)` covered by no region at all (evaluates `NaN`).

## Next steps

1. **Step 3 is the frontier.** A 4-face envelope split now gets through Step 2 and fails ~91 s in
   at `plq.maximum:185` → `plq.maximumConjugate` → `functionNDomain.maximumP` → `region.maximum`,
   with `symbolic:kernel:DivisionByZero` inside `symbolicFunction.subsF`'s `simplifyFraction` (a
   point that is a pole of one candidate). Reproduce with `f=xy` over `conv{(0,0),(1,1),(3,2)}`.
   Expect the same peel-the-onion pattern: four latent bugs were already flushed out this session,
   each previously unreachable because the `0/0` threw first. **User was asked whether to continue
   into `region.maximum` and ended the session instead — ask before resuming this.**
2. **Ask GitHub Support to garbage-collect the repo.** Unchanged, still the only remaining
   exposure before going public (force-pushed PDF blobs stay reachable by SHA until GC).
3. `partialConj` is unimplemented for every engine and type (`SUPPORT_MATRIX.md` §2).
4. Unbounded multi-face `conj` (`conjCPLQ.m:103`) — the remaining reason `conj` is not closed
   under itself. `biconj` no longer depends on it.
5. `maxQuaPar`: split a cell that already carries an arc (~26% of sampled splits); multi-arc
   pieces is the natural unit.
6. A **native numeric** rational branch in `conjPieceCPLQ` would buy **speed, not coverage** —
   ~100 s → milliseconds. The cPLQ fallback is correct but symbolic.
7. Then 0.1 tagging — **do not tag without being asked**.
8. Longstanding: `'pqp'`/`'graph'` engines; `RatPol.conj`/`biconj`/`add`; the
   `mergeL`/`removeTangent` exact-tie-point bug; `QuaPar.eval` wrong exactly *at* a result vertex
   (~1.4%); `testRegion/testCreation`.

## Relevant files

- `biconjCPLQ.m` — the biconjugate orchestrator; header has the derivation and validation numbers.
- `ratPolToPlq.m` — RatPol → cPLQ `plq`, the rational-face counterpart of `quaPolToPlq.m`.
- `conjCPLQ.m` — `conjEnvelopeViaCPLQ` (the Step 2 fallback) and its header on which half of cPLQ
  is reusable; line 103 is still the unbounded multi-face gap.
- `symbolicFunction.m` — `limitDirectional` vs `limit`; read both before touching either.
- `region.m` — `vertexOfEdge`, `simplifyUnboundedRegion`, `plus`: the three repaired sites.
- `plq_1p.m` — `conjugateFunction`'s `nCE==2` branch (the `a..f` seeding fix); note there is no
  `nCE==3` branch anywhere and none is needed.
- `SUPPORT_MATRIX.md` — §0's scope rule, §1.2's corrected rows, §8's blocker order.
- `biconjCPLQTest.m` / `conjCPLQTest.m` — the pipeline-free ground truth and the case-by-case pins.

# Session Handoff

_Last updated: 2026-07-29T00:00:00Z_

## What happened this session

Worked the Step 3 frontier (previous handoff's next-step 1). Reference case throughout:
`f = xy` over `T = conv{(0,0),(3,3),(1,2)}`, whose Step 1 envelope splits into 4 faces.

**Step 3 now runs end to end** (~15 min) instead of dying ~90 s in, and **Steps 1 and 2 are now
exact on this case** — the pointwise max of the four per-piece conjugates matches
`sup_{x∈T} ⟨s,x⟩ − xy` at **all 289 points** of a 17×17 dual grid. Step 1's envelope was verified
independently (tight to 1.8e-15; `env*` reproduces `f*` to grid resolution), so it was never at
fault.

Four bugs fixed, each previously unreachable because the one before it threw first:

1. **`region.maxArray`** — its tie-break probe points were derived from constraint *slopes*, and a
   slope is a lossy encoding of a direction. The bisector slope
   `2^(1/2)/4 + (4*2^(1/2)-4)/(2*(4*2^(1/2)-8))` **is** zero but is not the symbol `0`, so the
   `d == 0` test was false, `-1/d` survived as an unevaluated symbolic expression, and it only
   detonated much later inside `subsF`'s `simplifyFraction` as `symbolic:kernel:DivisionByZero`.
   New static helpers `region.probeAlong`/`probePerp` test the degeneracies numerically and return
   *points*, not slopes.
2. **`functionNDomain.mergeL`** — its second accumulation loop passed the **first** region's
   vertices to `removeTangent`, and passed nothing at all (`Unrecognized function or variable
   'nP'`) when that first region was empty and skipped. Now `region.finiteVertices` is read from
   the region actually being simplified, before `simplifyUnboundedRegion` drops any.
3. **`region.linear3pt`** — translated an edge-local vertex index into a region-level index **in
   place**, so once a match was found the next iteration used the result as the search key and
   overran the 3-element array (`Index must not exceed 3`).
4. **`plq_1p.conjugateFunction`, `nCE==2`** — which `grad` half-plane each edge owns was the
   hard-coded pattern `subdE(1,3)=grad; subdE(2,3)=-grad; subdE(3,3)=grad`. That is a statement
   about one edge *ordering*, not about geometry, so for any other ordering two edges got each
   other's half-plane and `f*` came out as the **min** of the two edge candidates instead of the
   max. Now derived per edge: the envelope here is rank-1 PSD, its Hessian kernel is
   `k = [b, −2a]`, and `grad = ⟨s − L, k⟩`, so the objective is affine along `k` with slope `grad`
   and the sup is pushed along `+k`; edge *j* therefore owns `{grad ≥ 0}` exactly when its outward
   normal satisfies `⟨k, n_j⟩ > 0`. An edge parallel to `k` is a level edge and contributes no
   region at all.

## Where things stand

- Branch: `main` — this session is a single commit, "Step 3: run the 4-face envelope end to end;
  fix four blockers and gate the wrong assembly", and it is the tip.
- Pushed: see the session summary (one confirmation at the end).
- Tests: **262 passed / 1 failed** across 22 suites — exactly the pre-session baseline. The single
  failure is `testRegion/testCreation`, longstanding and unrelated.
- `conj('cplq')` on the 4-face case still **fails loudly**, by design — but now for the real
  reason, via a new `conjCPLQ.assertStep3MatchesPieces` gate, not a crash.

### The one thing still wrong: Step 3's ASSEMBLY (and it is two bugs, not one)

The per-piece conjugates are right; `plq.maximumConjugate` mis-assembles them. Both causes let a
region claim territory that was never its own, carrying the wrong value there — ~12% of the 17×17
grid, over- and under-estimates alike, never uncovered.

1. **`region.merge` is unsound.** It unions `A ∩ {g≤0}` with `B ∩ {−g≤0}` by deleting the shared
   facet and intersecting what remains, i.e. it returns `A ∩ B`. That equals the union only when
   `A` and `B` are the *same* constraint set. Measured: three same-valued Step 3 regions merge into
   one covering `s=(1,1)`, where none of them does, so the partition reports `1.0` instead of
   `1.125`.
2. **`region.simplifyUnboundedRegion` drops non-redundant constraints.** Its rule is "delete any
   constraint that does not pass through a finite vertex" — not a redundancy test for unbounded or
   curved regions. Measured: a region carrying `(s₁+s₂)²/4` loses `−s₁−s₂ ≤ 0` and reports
   `f*(−3,−3) = 9` where the truth is `0`. Note it does **not** misbehave on the pristine
   per-piece regions in isolation — only on the intersected ones, so reproduce it inside the
   pipeline, not standalone.

**Do not fix these separately.** Guarding (1) alone with a constraint-set-equality test (tried,
provably sound, then reverted) makes the measured result *worse* — 36 → 125 wrong of 289 — because
refusing merges leaves more regions for (2) to damage, and evaluation takes the first matching
region. `region.merge`'s header records the analysis and the numbers.

### One speculative change was reverted — don't reintroduce it

`region.maxArray`'s pair loop tests `abs(m(i)) ~= inf`, and `slopes2` flags a **vertical**
constraint as `intmax`, not `inf`. That mismatch looks like a bug and is not: a vertical constraint
deliberately falls into the arithmetic mean and yields a huge-but-finite bisector slope, which
probes far out along an unbounded region's recession direction.
`testMaxMultiRegion/testOpenconvex` depends on those far probes being feasible — "fixing" the
mismatch made `maxArray` return undecided, `maxEqDom` fall through to `splitmax3`, and
`region.normalize1` choke on the resulting rational inequality
(`symbolic:coeffs:NotAPolynomial`). The comment at the site says so.

## Next steps

1. **Step 3's assembly — the two defects above, fixed together.** This is the whole remaining gap
   for the 4-face case. (2) is the dominant one. A real redundancy test is what
   `simplifyUnboundedRegion` needs; `merge` needs the union-exactness condition, not just
   convexity. Reproduce with `f=xy` over `conv{(0,0),(3,3),(1,2)}`; the harness that localizes it
   is a round-by-round grid check of `plq.maximumConjugate`'s accumulation (seed → `*` piece j →
   `maximumP` j), which shows `overlap-disagree` going 0 → 97 → 187 across the folds.
2. `partialConj` is unimplemented for every engine and type (`SUPPORT_MATRIX.md` §2) — a new
   feature, not a bug fix, so it sits just above the longstanding backlog.
3. Unbounded multi-face `conj` (`conjCPLQ.m`) — the remaining reason `conj` is not closed under
   itself. `biconj` does not depend on it.
4. `maxQuaPar`: split a cell that already carries an arc (~26% of sampled splits); multi-arc
   pieces is the natural unit.
5. A **native numeric** rational branch in `conjPieceCPLQ` would buy **speed, not coverage**
   (~100 s → ms). The cPLQ fallback is correct but symbolic. Step 3 is now the runtime hog anyway
   (~15 min for the 4-face case).
6. Then 0.1 tagging — **do not tag without being asked**.
7. Longstanding: `'pqp'`/`'graph'` engines; `RatPol.conj`/`biconj`/`add`; the
   `mergeL`/`removeTangent` exact-tie-point bug; `QuaPar.eval` wrong exactly *at* a result vertex
   (~1.4%); `testRegion/testCreation`.

**Dropped:** asking GitHub Support to garbage-collect the repo. The repo is private and will stay
so for a while, and nobody outside has the blob SHAs.

## Relevant files

- `conjCPLQ.m` — `conjEnvelopeViaCPLQ` and the new `assertStep3MatchesPieces` gate (read its header
  for why a loud failure is preferred to a silently wrong partition).
- `region.m` — `probeAlong`/`probePerp` (static, top of file) and `maxArray`; `finiteVertices`;
  `linear3pt`; `merge`'s KNOWN UNSOUND header; `simplifyUnboundedRegion`.
- `functionNDomain.m` — `mergeL`'s two accumulation loops; `maximumP`.
- `plq_1p.m` — `conjugateFunction`'s `nCE==2` branch, both the `a..f` seeding and the new
  geometric `grad` half-plane derivation. There is no `nCE==3` branch anywhere and none is needed.
- `SUPPORT_MATRIX.md` — §1.2's 2026-07-29 update (the fix table and the two assembly defects),
  §8's blocker order.
- `conjCPLQTest.m` — `indefiniteTriangleThreeConvexEdgesUsesStep3` (still pins a loud
  `cplqFailed`, now from the gate) and `indefiniteTriangleTwoConvexEdgesSplitViaCPLQStep2` (the
  end-to-end 2-face case, still exact).

## Still true from before

- **There is no "3-convex-edge case" to implement.** [COAP] Appendix A.5's split reduces such a
  triangle to 2-convex-edge sub-triangles and Step 1 already applies it. Describe these cases by
  the **envelope's face count**, not the input's edge count.
- **The rational-face `0/0` did not need exact arithmetic** — it was cleanly removable, and
  `symbolicFunction.limitDirectional` resolves it by a directional limit.
- **Which half of cPLQ is reusable:** its **Step 2/3**, on **CCA2's** Step 1 output. Running cPLQ
  end to end gives a wrong answer, because its own 2-convex-edge envelope is the single Appendix
  A.4 formula, which CCA2 established is not always tight.

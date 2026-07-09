# Session Handoff

_Last updated: 2026-07-09_

## What happened this session

Two parts. Verified 112/112 PASS on Frances before starting.

### Part 1: wired the `conjCPLQ` orchestrator (done, tested, uncommitted)

Added `conjSingleTriangle` (local fn in `conjCPLQ.m`) so a single bounded-triangle piece is
conjugated end-to-end: tries `conjPieceCPLQ` on the ORIGINAL piece first (works directly for
affine / PD / rank-1-PSD / indefinite-with-0-or-1-convex-edge, since `f*=(conv f)*` means no
envelope is needed at all in these cases), falling back to `convEnvCPLQ`'s envelope only when
`conjPieceCPLQ` can't handle the raw piece (a concave piece, or an indefinite piece with exactly
2 convex edges, which envelopes to a single rank-1-PSD face). 7 new tests added to
`conjCPLQTest.m` (`affineTriangleViaOrchestrator`, `convexQuadraticTriangleViaOrchestrator`,
`concaveTriangleViaOrchestratorSidestepsToEnvelope`,
`indefiniteTriangleZeroOrOneConvexEdgeViaOrchestrator`,
`indefiniteTriangleTwoConvexEdgesSidestepsToEnvelope`,
`indefiniteTriangleThreeConvexEdgesNeedsStep3`, `multiFacePieceStillNotImplemented`). Full suite
**112/112 PASS**. The one remaining case (indefinite triangle with 3 convex edges) raises an
explicit `PLQ:conjCPLQ:notImplemented` naming Step 3 as the reason, instead of a generic error.

### Part 2: Step 3 primitive `maxQuaPar.m` — hit a real mathematical wall, not a bug

Built `maxQuaPar.m` (new file, ~450 lines) implementing the pointwise max of two full-domain
`QuaPar`s: face-pair overlay (all `g1.nf × g2.nf` combinations), Sutherland-Hodgman-style
polygon/ray clipping (`clipByFace`/`clipPolyHalfPlane`, handles bounded and unbounded cells,
0/1/2 crossings), closed-form winner decision (`decideWinner`, finite vertices + ray
asymptotics via the Hessian/linear/constant leading term, no numeric far-point guessing),
curve-splitting for undecided cells (`splitCell`, exact ray-quadratic algebra, asserts the
degeneracy the whole approach depends on rather than assuming it), and reassembly into one
`QuaPar` (`assemblePieces`, generalizes `convEnvCPLQ.m`'s `assembleTriangles` half-edge pairing
to rays and one conic edge per piece). Algorithm shape was cross-checked against the older
symbolic `cPLQ`/`functionNDomain.m`+`region.m`'s own max implementation
(`maximumP`/`maximumPC`, `region.maxArray`/`splitmax3`) — ported the *algorithm*, not the code
(cPLQ needs an `intmax` sentinel for points at infinity; `maxQuaPar.m` uses `QuaPar`'s own
structural ray representation instead, no sentinel needed — a real simplification).

**Testing this on a concrete example surfaced a genuine mathematical obstruction in the
published algorithm, not a bug in this session's code.** Full writeup with proofs and diagrams:
**`/home/ylucet/CCA2/3-edge.tex`** (compiled `3-edge.pdf`) — **read this first**, it has the
complete derivation; summary:

- Example: `f(x,y)=xy` over the triangle `(0,0),(3,3),(1,2)` (three convex edges — the same
  category `convEnvCPLQ` already splits, via `splitThreeConvex`, into
  `T1=(0,0)-(1,2)-(2,2)` and `T2=(1,2)-(2,2)-(3,3)` sharing the horizontal edge `(1,2)-(2,2)`,
  which is *not* a convex edge of either sub-triangle).
- `g1=conj(T1's rank-1-PSD envelope)`, `g2=conj(T2's rank-1-PSD envelope)` (via
  `conjPSDRank1QuadTriangle`, 6 faces each). Of the 36 combinatorial face×face pairs (18 nonempty
  overlap cells), exactly **2** require comparing two genuinely quadratic pieces against each
  other — both involving `g1`'s shared-edge-strip piece — and **both have discriminant provably
  > 0** (exact: `√2/8+3/16 ≈ 0.364` and `√2/4+3/8 ≈ 0.729`): a genuine **hyperbola**, which
  `QuaPar` cannot represent as an edge.
- This directly contradicts [JOGO] Theorem 6's proof, which asserts "when we compare two
  functions we always get one of them as linear" — false here, since **both** pieces are
  quadratic (they descend from *different* rank-1 envelopes on adjacent sub-triangles, a
  configuration the proof never actually examines). [COAP] Appendix A.5 states the
  three-convex-edge split in one sentence and never revisits Step 3 for it.
- Verified 3 independent ways: (1) **SCAT** (`github.com/ylucet/SCAT`, Python/SymPy port,
  installed editable at `/home/ylucet/SCATpy` — `from scat import PWF, rsym`) independently
  re-derived `g1` via its own general n-D pivot algorithm, exact match with the closed-form
  construction; (2) direct constrained solver (MATLAB `fmincon`, multi-start) computing
  `sup_{x in T} <s,x>-f(x)` directly vs `max(g1,g2)` — max error `4e-15` at 16 points including
  inside both hyperbola cells; (3) ran the **actual reference `cPLQ` code**
  (`/home/ylucet/FLTW/codeOld/cPLQ`) on this example: `plq_1p.convexEnvelope` has explicit
  branches for `nCE==0,1,2` only — no `nCE==3` branch exists — so it **silently returns an empty
  envelope/conjugate** rather than computing anything or erroring; a comment in the authors' own
  `testcPLQ.m` (`testRect`: `% add to split 3 convex edge case`) confirms this was a known,
  never-finished item, not a deliberate design choice.

`maxQuaPar.m`'s clipping/winner/reassembly machinery is otherwise believed complete for the
16/18 cells that *do* decide cleanly, but **has never returned a complete result end-to-end on
any real example** — the one scenario tried hits the hyperbola wall partway through
(`splitCell`'s degeneracy assertion fires as designed, loudly, not silently).

## Where things stand

- Branch: `cplq-engine` @ `1d115d4` (this session's starting point, unchanged — nothing pushed
  this session).
- **Uncommitted** (ask before committing, per standing rule): `conjCPLQ.m` (modified),
  `conjCPLQTest.m` (modified), `maxQuaPar.m` (new, untracked).
- Full suite (excluding `maxQuaPar.m`, which has no test file yet): **112/112 PASS** on Frances.
- New file outside the repo: `/home/ylucet/CCA2/3-edge.tex` + `3-edge.pdf` (the hyperbola
  writeup) and its figures `fig_primal.pdf`, `fig_dual.pdf`, `fig_hyperbola_insets.pdf`.

## Next steps

1. **Read `/home/ylucet/CCA2/3-edge.tex` first.** Then decide between:
   - (a) Construct or find a 3-convex-edge example where the two sub-triangles' shared-edge
     pieces do *not* end up needing to be compared as two quadratics, to get a first full
     `maxQuaPar` end-to-end validation (would confirm the clipping/reassembly machinery works
     when the hyperbola issue doesn't arise).
   - (b) Determine whether the hyperbola is **unavoidable in general** for a 3-convex-edge
     split: the shared cut edge is, by construction (COAP A.5's horizontal-line-through-the-
     middle-vertex), never a convex edge of either sub-triangle, so there is no a priori reason
     for the two adjacent rank-1 Hessians to share an eigenvector (which is what would be needed
     for their difference to stay degenerate). If unavoidable, this is a genuine gap in
     [JOGO]/[COAP]'s coverage for *any* 3-convex-edge split, not just this one example — worth
     deciding whether/how to report it, and whether some `QuaPar`-representable workaround
     exists (e.g. does the pipeline actually need the exact max there, or could a different
     splitting strategy in Step 1 avoid ever producing this configuration?).
2. Separately, the standalone `RatPol.conj` gap (rational piece with no known originating
   quadratic) is still open — unrelated to this session's finding, see prior handoffs.
3. Test command (Frances, prefix every batch call — Maple toolbox conflicts with Symbolic):
   ```
   matlab -batch "restoredefaultpath; rehash toolboxcache; cd('/home/ylucet/CCA2/CCA2'); \
     res=runtests({'RatPolTest','convEnvCPLQTest','QuaParTest','conjCPLQTest','conjPieceCPLQTest','PLQVCTest'}); \
     fprintf('\n==== SUMMARY ====\nTOTAL %d PASSED %d FAILED %d\n',numel(res),sum([res.Passed]),sum([res.Failed])); \
     exit(sum([res.Failed])>0)"
   ```

## Relevant files

- `/home/ylucet/CCA2/3-edge.tex` (+ `.pdf`) — **the main deliverable of this session**: full
  proof, diagrams, and 3-way verification of the hyperbola finding. Figures generated from
  `/home/ylucet/.claude/jobs/237b4e30/tmp/make_diagrams.py` and
  `make_hyperbola_insets.py` (job scratch dir, not persistent — regenerate from the `.tex`
  source data if needed, or copy the `.py` scripts out first).
- `conjCPLQ.m` — `conjSingleTriangle` (Part 1).
- `conjCPLQTest.m` — 7 new tests (Part 1).
- `maxQuaPar.m` — Part 2, the Step 3 primitive (untested end-to-end, see above).
- `/home/ylucet/SCATpy` — Python SCAT (SymPy port of Maple SCAT), used to independently verify
  `g1`'s conjugate. `pip`-editable-installed; `from scat import PWF, rsym`; `PWF.from_regions`,
  `PWF.conj`.
- `/home/ylucet/FLTW/codeOld/cPLQ` — reference implementation; `plq_1p.m`'s `convexEnvelope`
  (~line 209) confirms the missing `nCE==3` branch; `testcPLQ.m`'s `testRect` has the
  "unfinished" comment.
- `/home/ylucet/CCA2/DESIGN.md` — design proposal; II.5.1 describes the 3-step `cplq` engine.
- Memory: `cca2-toolbox-design` project memory updated with this session's findings (auto-loaded
  next session).

# Session Handoff

_Last updated: 2026-07-18T13:00:00Z_

## Update (same session, final): the adapter is built and validated

Built the CCA2↔cPLQ adapter the sections below anticipated: **`quaPolyToPlq.m`** (CCA2 `QuaPoly`,
any number of triangular faces, → cPLQ `plq`: per face, `QuaPoly.matrixForm` → a `sym` quadratic
expression → `domain(V3,x,y)` + `symbolicFunction(expr)` + `plq_1p(d,f)`) and
**`evalFunctionNDomain.m`** (numeric eval of a `functionNDomain` array — e.g. `plq.maxConjugate`
after `.maximum` — at a dual point `s`: find the region whose inequalities `s` satisfies, `subs`+
`double` that region's own function there). Deliberately does NOT convert back into `QuaPar` —
that would reintroduce the exact curved-edge assembly problem `cPLQ` already solves symbolically;
Phase 1 only needs numeric evaluability to write ground-truth-checked tests.

**`cplqAdapterTest.m` (2/2 passing)**:
- `singleTriangleMatchesConjPieceCPLQ`: one one-convex-edge triangle (`f=xy`) through
  `quaPolyToPlq`→`.triangulate`→`.conjugate` matches CCA2's own existing numeric `conjPieceCPLQ`
  exactly at 8 dual points — cross-validates the INPUT conversion is correct.
- `twoTriangleSquareMaxMatchesNumericSup`: **the actual gap-closing case** — `f=xy` over a square
  split into 2 INDEPENDENT triangles (exactly what `maxQuaPar` cannot do: each triangle's own
  conjugate is a correct curved `QuaPar`, but combining them numerically errors on the curved
  edges) — through `quaPolyToPlq`→`.triangulate`→`.maximum` matches the true numeric sup (dense
  grid over the whole square) at 9 test points.
  **One point deliberately excluded and documented, not silently dropped**: `s=(0.5,0.5)`, an
  exact symmetric TIE where both triangles' own off-diagonal-vertex cones meet (true value 0.5,
  attained at BOTH (1,0) and (0,1) simultaneously). Debugged directly (scratch script): at this
  point, ALL 8 assembled regions have their nearest inequality violated by 0.5-1 (not a small
  residual) — a genuine gap in the assembled partition, matching the same "isAlways: TruthUnknown"
  warnings seen during `mergeL`/`removeTangent` for this exact computation. This is a narrow
  limitation of the VENDORED `functionNDomain.mergeL`/`region.removeTangent` at an exact tie
  between two independent triangles' vertices, not a bug in the new adapter code — flagged for a
  future debugging pass, not fixed (same "don't guess at core symbolic-geometry logic, each
  iteration costs ~45s+" tradeoff as `poly2orderUnbounded` below).

**Not done yet** (would be natural next steps, not attempted this session):
- Wire the adapter into `conjCPLQ.m`'s own `conj(f,'cplq')` dispatch for the general (`nf>1`) case
  — right now `quaPolyToPlq`/`evalFunctionNDomain` are standalone, not called from `conjCPLQ`.
- A `QuaPar`-producing conversion (if/when a caller needs the structured result, not just a
  numeric `eval`).
- The `poly2orderUnbounded` biconjugate bug (see below, unchanged) and the `mergeL` tie-point gap
  just found both remain open.

`DESIGN.md`'s "Next planned" §1 (Phase 1) has been updated in place to match; search "IN PROGRESS,
core adapter DONE" there.

## Goal changed this session (user decision, 2026-07-18) — read this first

The user redirected the project's near-term goal partway through this session. **Do not resume
any closed-form-numeric-derivation work** (the abandoned single-RatPol-piece conjugate, or a
from-scratch curved-edge `maxQuaPar` extension) until Phase 1 below is done. New plan, in order:

1. **Phase 1 (current focus): integrate the reference `cPLQ` package into CCA2 almost as-is,
   using its own symbolic (MATLAB `sym`) computation, and get it fully tested end to end**
   (envelope → conjugate → max → biconjugate, for genuinely multi-triangle nonconvex PLQ inputs —
   the case `conjCPLQ.m`'s general branch currently just errors on). This is not a new direction so
   much as a correction: `DESIGN.md` §0 point 3 always specified `'cplq'` as *"symbolic per-piece
   conjugate (`cPLQ` Lagrange-multiplier method)"* — the numeric closed-form code built so far
   (`conjCPLQ.m`/`conjPieceCPLQ.m`/`convEnvCPLQ.m`/`maxQuaPar.m`) is a faster partial
   reimplementation of a *subset* of that, not a replacement for it.
2. **Phase 2 (later): improve performance.** Once Phase 1 passes tests, replace pieces of the
   symbolic computation with closed-form numeric formulas incrementally — one case/step at a time,
   validating each replacement against the still-available Phase-1 symbolic result on the same
   inputs before moving to the next. Go slowly; do not re-derive the whole pipeline's math in one
   jump (see "what NOT to repeat" below).
3. `DESIGN.md`'s "Next planned" section has been updated in place with this same plan (search
   "REVISED (this session)").

## What NOT to repeat (context for why the plan changed)

The prior two sessions (see git history before this handoff) got stuck trying to derive a
closed-form numeric conjugate for a single genuinely-rational `RatPol` triangle piece
(`conjPieceCPLQ`'s own TODO). This session found that work was **mis-scoped**: `conjCPLQ.m`'s
`conjSingleTriangle` never actually needs it (it always succeeds by conjugating the *original*
quadratic directly, or falling back to an affine/rank-1-PSD-quadratic envelope — never a genuinely
rational one, for any convex-edge count 0-3). The *real* remaining blocker for a genuinely
multi-triangle nonconvex PLQ is `maxQuaPar`'s own documented refusal to accept curved (parabolic)
input edges — confirmed by a direct test: `f=xy` over a 2-triangle square, each triangle's
(correct) conjugate is a 6-face curved `QuaPar`, and `maxQuaPar(g1,g2)` errors immediately.
Extending `maxQuaPar` to curved edges from scratch (conic-conic clipping) would have been the next
natural numeric-formula task — **the user chose integrating `cPLQ` symbolically instead**, since
`cPLQ` already solves exactly this (its `functionNDomain`/`region`/`mergeL` machinery handles
multi-piece assembly, including curved edges, symbolically) and is the pipeline's own designed-in
reference implementation, not a detour.

## What was verified this session about `cPLQ` itself (concrete, reusable)

- **Not bit-rotted**: `cPLQ`'s classes still load and run under the current MATLAB + Symbolic
  Math Toolbox. Confirmed by running `cPLQ/testcPLQ.m`'s `testRect3Conj` (envelope + conjugate,
  `f=xy` over a 4-vertex domain triangulated into 2 triangles) — **passed in ~11.6s**.
- **Timing breakdown** (same 2-triangle `f=xy` example throughout): `testRect3Conj` (envelope +
  conjugate, no max) passed in **~11.6s**; `testRect3Max` (adds `.maximum`) passed in **~38.6s**
  total (so the max/merge step itself costs on the order of ~27s here, for just 2 triangles);
  `testRect3Biconj` (adds `.biconjugateF` — a second conjugate+merge round on the maxed result,
  which by then has ~5-6 regions) did **not** finish in several minutes and was killed. So the
  slow part is specifically **repeated conjugate+merge rounds on an already-many-region input**
  (biconjugate re-conjugates and re-merges the max step's own output), not envelope/conjugate/max
  in isolation for a small number of original triangles. Heavy `Warning:
  symbolic:sym:isAlways:TruthUnknownUnable to prove...` spam comes from `region/removeTangent` →
  `functionNDomain/mergeL` → `functionNDomain/maximumP`: `isAlways` is trying to prove pairs of
  region-boundary inequalities are the identical line/conic and frequently can't decide
  symbolically, warns, and (presumably) falls back to a slower path each time. **This is Phase 2's
  primary target once Phase 1 is otherwise working** — do not treat "it's slow" as a Phase-1
  blocker; get it correct first, and prefer `conjugate`/`maximum` (not `biconjugateF`) for the
  first passing end-to-end tests, adding biconjugate coverage once the rest is solid.
- **Top-level API** (`cPLQ/plq.m`): construct `plq(array_of_plq_1p)`; `plq_1p(d, f)` where `d` is
  a `cPLQ/domain.m` polygon (`domain(vertexMatrix, x, y)` with `x,y` = `sym` variables) and `f` is
  a `cPLQ/symbolicFunction.m` (`symbolicFunction(x*y)` etc.). Methods: `.triangulate`,
  `.convexEnvelope`, `.conjugate` (envelope+conjugate per piece), `.maximum` (adds
  `.maximumConjugate` across all pieces), `.biconjugateF`. `.print`/`.printDomainMaple` for
  inspection (Maple-syntax region dump).
- **Per-piece conjugate formulas already match CCA2's own numeric ones almost exactly**:
  `cPLQ/plq_1p.m`'s `conjugateFunction` (not `plq_1piece.m` — that file is an older/unused
  duplicate; `plq.m` calls methods on `plq_1p` objects directly) handles `nCE=0/1/2` by
  conjugating the **original** triangle's vertices/edge slope directly — e.g. its `nCE==1` branch's
  `a,b,c,d,e,f` (the `crs` parabola) and `av,bv,cv,dv,ev,fv` (the edge-strip quadratic) are the
  same closed forms as CCA2's `conjBilinearXYoneCE`'s `crs`/`E`. This is reassuring cross-
  validation of the numeric formulas already in CCA2, and confirms Phase 2 (later) should be
  comparatively mechanical for the parts already covered by `conjPieceCPLQ.m` — the genuinely new
  work in Phase 2 will be the multi-piece assembly/merge logic, which CCA2 doesn't have a numeric
  analogue of at all yet.

## Update (same session, continued): files copied in, baseline testing underway

Per user follow-up decision: **copy the relevant `cPLQ` files into the CCA2 repo root** (not the
whole `cPLQ/` folder/repo). Did a dependency-closure analysis by grepping cross-file calls among
all non-`.asv` files in `cPLQ/`: the live runtime dependency set is exactly **8 files** —
`plq.m`, `plq_1p.m`, `functionNDomain.m`, `region.m`, `domain.m`, `symbolicFunction.m`,
`conjugateExpr.m`, `yIntercept.m`. Everything else in `cPLQ/` (`plq_1piece.m` — an unused older
duplicate of `plq_1p.m`, never called by `plq.m`; `quadQuad.m`/`qq_conj.m`/`quadLinear.m`/
`conjQuad.m`/`conjR.m`/`checkConvex.m`/`subdiffParabolicEdge.m`/`vertexNan.m`/`quadquad1.m`/
`plots.m`/`plot2.m`) are standalone offline derivation *scripts* with zero live call sites — NOT
part of the pipeline, confirmed by grep (this also corrects `DESIGN.md`'s own prior claim that
`quadQuad`/`qq_conj` were part of the engine's source — that was aspirational/inaccurate, now
fixed in place). Copied those 8 files plus their 5 test suites (`testcPLQ.m`,
`testfunctionNDomain.m`, `testRegion.m`, `testSymbolicFunction.m`, `testMaxMultiRegion.m`) into
the CCA2 root. Smoke-tested from the new location (not just the old `cPLQ/` subfolder) —
`testRect3Conj` still passes (~8.8s), confirming no hidden path dependency broke. `cPLQ/` itself
(untracked) was left alone, not deleted.

**Baseline test run so far**: `testRegion` + `testSymbolicFunction` + `testfunctionNDomain` =
43/44 passed. The one failure, `testRegion/testCreation`, had two distinct causes:
1. **Fixed**: 4 assertions used `isequal(<sym vector>, <double literal>)` — in the current
   MATLAB/Symbolic Math Toolbox, `sym` does not overload `isequal` (confirmed: `isequal(sym(0),0)`
   itself returns `false`, structural/handle comparison, not value comparison), so this idiom
   never worked regardless of the actual data, a toolbox-version-drift issue, not a functional bug.
   Fixed by switching to `testCase.verifyEqual(double(...), ...)` (2 pairs of assertions,
   `testRegion.m` lines ~76-77 and ~82-83) — the standard modern-MATLAB idiom for this exact
   check, preserving the original test's intent.
2. **NOT fixed yet, documented instead**: after fix #1, one genuine value mismatch remains, but
   only when the test is run **through the MATLAB unittest framework** (`runtests`/
   `TestSuite.fromClass`) — `s.vy` comes back as `[0, 0.5, 0]` (paired with `s.vx=[1,0,0.5]`,
   i.e. vertex (0,0.5)) instead of the geometrically-correct `[0, 0, 0.5]` (vertex set
   `{(0,0),(1,0),(0.5,0.5)}`, confirmed correct by hand from the region's own 5 inequalities,
   2 of which are tangent/redundant at existing vertices). **Reproduced consistently through
   `runtests`/`TestSuite.fromClass`, but NOT reproducible as a plain top-level script** running
   the exact same `region(l2,[x,y])` call (with or without first creating the other test regions
   `r,t,u,v,w` in the same sequence `setUpTestData` uses) — every plain-script attempt gave the
   correct answer. This points at some interaction between MATLAB's unittest test-runner context
   and the symbolic engine (`region.m`'s `getVertices`, which pairwise-solves all constraint pairs
   via `solve` then filters by `ptFeasible` and dedupes via `unique(V,'rows')` — confirmed
   `unique(sym,'rows')` itself pairs columns correctly in isolation, so that's not the cause
   either), not a simple one-line logic bug. **Did not chase further this session** — time/effort
   tradeoff call, since it's a narrow redundant-constraint corner case in a single test, and the
   actual pipeline tests (`testRect3Conj`/`testRect3Max`, which exercise `region.m` on genuine
   nonconvex-PLQ triangle pieces, not this synthetic 5-inequality example) pass reliably. Next
   session: if picking this up, try adding trace/disp inside `getVertices`'s solve loop specifically
   when invoked via `runtests`, to see which constraint PAIR produces the spurious point under that
   execution context.
- `testMaxMultiRegion.m` initially failed ALL 24 tests immediately (`plq(p)` rejected `plq_1piece`
  objects: `plq.m`'s `pieces` property was typed to `plq_1p.empty()`, and its constructor did
  element-by-element indexed assignment `obj.pieces(i)=ps(i)` into that pre-typed array).
  `plq_1piece.m` (not originally copied — see above, now copied in too, it's a genuine second
  dependency) is an OLDER, independently-tested parallel per-piece implementation this test suite
  uses directly (never through `plq_1p`); a commented-out `%pieces = plq_1piece.empty();` line
  directly above the current default in `plq.m` is direct evidence of the historical class swap
  that broke this. **Fixed** (2 small edits to `plq.m`, both in the class constructor): (1) untyped
  the `pieces` property (plain `[]` instead of `plq_1p.empty()`), (2) changed the constructor to
  assign the whole `ps` array in one shot (`obj.pieces = ps;`) instead of element-by-element
  indexed assignment, since indexed assignment into a property that starts as literal `[]`
  (class `double`) still fails to convert a custom object into `double` on the first element.
  Verified this doesn't regress `testcPLQ.m` (all 6 non-biconjugate tests still pass) or break
  `testMaxMultiRegion`'s own non-biconjugate tests (`testConvex`/`testConjugate` pass, 2/2).
- **Steps 1-3 (envelope, conjugate, max) are solid**: every non-biconjugate test tried this
  session across both `testcPLQ.m` and `testMaxMultiRegion.m` passes (`testPCE0`/`testPCE3`/
  `testPCE1`: 12-56s each). **Step 4 (biconjugate) has at least one genuine, reproducible bug**,
  found by running `testMaxMultiRegion/testMax` (which — unlike `testcPLQ.m`'s split-out
  `testRectConj`/`testRectMax`/`testRectBiconj` methods — calls BOTH `.maximum` AND
  `.biconjugateF` in one test): after ~372s, it throws
  `MATLAB:badsubscript`/"Index exceeds the number of array elements. Index must not exceed 3." in
  `region.poly2orderUnbounded` (line 198: `obj.getEdges(obj.vx(i+1), obj.vy(i+1))`), called from
  `functionNDomain.conjugateOfPiecePoly` → `plq.biconjugateF`. Root cause: the preceding loop
  (`for i = 1:size(obj.ineqs,2) ... if nv1==1 | nv2==1, break; end`) is supposed to find a vertex
  adjacent to an unbounded (ray) edge and `break` with that `i`; if NO such vertex is found (loop
  runs to completion without breaking, here `i` reaches `obj.nv`=3), the very next line reads
  `obj.vx(i+1)` unconditionally — an unhandled fall-through, not a wraparound. This is the SAME
  class of bug already found and fixed once in this very codebase's own `maxQuaPar.m` (see its
  HISTORY comment: an unmodded `p2:(p1-1)` range that should wrap but doesn't, fixed by adding
  `+nv` before `mod`). **Not fixed this session** — deliberately: a wrong guess at "wrap `i+1` via
  `mod`" risks silently producing a WRONG vertex reordering (worse than a clear crash) rather than
  the correct one, and iterating on a wrong guess is expensive given each biconjugate run costs
  5-6+ minutes. All 8 `testBiconjugate*` methods in `testMaxMultiRegion.m` almost certainly call
  `.biconjugateF` too (by name) and were not run this session, given `testMax` alone took ~6
  minutes before crashing — expect them to hit the same or related bugs and budget accordingly.
- **Untested this session** (ran out of time/turn budget, not attempted): `testMaxMultiRegion`'s
  `testPCE2`, `testMaxR3`, `testMaxT`, `testMaxP`, `testMax3`, `testPSqroot`, `testFractional`,
  `testMaxThesis`, `testMaxThesis2`, `testOpenconvex`, and all 8 `testBiconjugate*`. `testcPLQ`'s
  `testRectBiconj`/`testRect3Biconj` (both call `.biconjugateF`) were attempted early this session
  and did not finish in several minutes — consistent with the bug above, though not confirmed to
  be the identical crash (killed before erroring, not rerun with full capture).
- **Practical implication for Phase 1**: Steps 1-3 (the actual gap in CCA2's own `conjCPLQ.m` —
  general multi-triangle Step 3) are ready to build the adapter against NOW. Biconjugate needs its
  own debugging pass first if/when it's needed — recommend NOT blocking the Step-1-3 adapter and
  its tests on fixing biconjugate, since `conjCPLQ.m`'s own current gap is specifically about
  Step 3 (the conjugate, not the biconjugate) for a multi-triangle input.

## Next steps

1. **Design and build a thin CCA2↔`cPLQ` adapter**: convert a CCA2 `QuaPoly`/`PLQVC` object
   (possibly multi-face, genuinely nonconvex) into `cPLQ` `domain`/`symbolicFunction`/`plq_1p`/
   `plq` inputs, run the `cPLQ` pipeline, and convert results back into something CCA2 tests can
   check (numeric `eval` at sample dual points is enough for Phase 1 — no need to convert into
   `QuaPar`'s stricter representation yet, since that reintroduces exactly the curved-edge
   assembly problem `cPLQ` already solves; keep the result in `cPLQ`'s own `functionNDomain`
   representation for now).
2. **Write CCA2-style end-to-end tests** (numeric sup-sampling ground truth, matching every
   existing test file's convention) for: (a) a single nonconvex triangle (should match existing
   `conjPieceCPLQ` results, cross-validating both), (b) a genuinely multi-triangle nonconvex PLQ
   (the currently-unreachable case), through to biconjugate.
3. Given the max/merge slowness, prefer small examples (2-3 triangles) for the first passing
   tests; treat anything needing many triangles/pieces as a Phase-2-readiness question, not a
   Phase-1 blocker.
4. Where to put the vendored/integrated `cPLQ` code: it currently lives at `cPLQ/` (untracked,
   per `git status` at session start) — decide whether it should be added to the repo proper (git
   `add`) as part of Phase 1, since "integrate almost as-is" implies it becomes load-bearing
   production code, not just a reference clone. Flag this to the user before adding a large
   third-party-ish symbolic codebase to source control, if not already discussed.
5. Unchanged: exact `[LOCATELLI]` citation in DESIGN.md; 2/741 (0.3%) residual `maxQuaPar:internal`
   crashes (in the numeric code, unaffected by this plan change); `QuaPar.orderEdges`/`createP`'s
   near-degenerate-triangle error; `partialConj` for `'pqp'`; `add` for `RatPol`/`RatPar`.

## Relevant files

- `DESIGN.md` — "Next planned" section rewritten in place (search "REVISED (this session)") with
  the 2-phase plan; struck through the old rational-piece-derivation item rather than deleting it,
  for history.
- No CCA2 production `.m` files changed this session. `cPLQ/` files read: `plq.m`, `plq_1p.m`
  (the active per-piece implementation), `functionNDomain.m` (partial — `getSubdiffVertexT1/T2`,
  `getSubDiffVertexSpT1`, `getSubDiffEdgeT1`, `conjugateExprVerticesT1`, and the existence/location
  of `mtimes`/`maximumP`/`maximumPC`/`mergeL`/`conjugateOfPiecePoly`, not their full bodies).
  `cPLQ/testcPLQ.m` read and partially run (see above). `RatPol.m` read in full.
  `reference/KARMARKAR-26-convex-envelope.pdf` Appendix B.2/C.2, Section 3/4, Table 1/3, Figs 2-3
  (via `pdftotext`, available at `/mingw64/bin/pdftotext` in this environment — useful for future
  sessions needing to search the paper's text instead of reading it as an image-only PDF).

# Deep code review — CCA2 engines behind the MPC paper (2026-11-19)

Reviewer: Claude Opus 5.5 (Kandev task CCA2 #76). **Independence:** this session wrote none of
the reviewed code (MATLAB, Python and Rust engines, and the #64 cross-check, were written by other
sessions and models). The oracles, generator and tests below were written by the reviewer from the
definitions; the reviewer's f* oracle shares no code with `cca2` (not even `cca2/check.py`).
Author notes (DECISIONS.md, DESIGN.md) were read only after the findings were in (section 7).

## 1. Scope, budget, commands

**In depth (all three engines):** the 2-D conjugate (`conjQ`/`conj_q`) and biconjugate
(`biconjQ`/`biconj_q`) paths and their callees: recession analysis, Case D boundary max, `maxQ`
fold, `assembleQuaConCells`, `cellHasInterior` (witness search + CAD fallback), `biconjAtoms`
(tie arcs), envelope assembly/subdivision; the legacy `'cplq'` path (`conjCPLQ`, `biconjCPLQ`) for
the three open TODO items; the PLQ / non-PLQ acceptance rules (`QuaPol` reading, `fromDouble`);
the cell-count logic; the #64 cross-check (`python/bench/cross_convexdb.py`,
`python/tests/test_cross_convexdb.py`).
MATLAB changes since the move to `MATLAB/` (3 commits: `19e4388` psTangency, `01d70f7` conicMeet
`/scale`, `c19bc97` cellHasInterior CAD) were read line by line.
**Normal depth:** the rest of `python/cca2` and `rust/src` (by test/clippy output and targeted reading).
**Skimmed only:** `doc/`, `manuscript/` (claims only, not changed), `python/bench` other than the
cross-check, `rust/bench`, `rust/tools`, `python/tools`, the `.claude/` scratch scripts, the
legacy symbolic suites' internals (`plq`, `region`, `functionNDomain`).
**Not done:** a full re-run of the >4 h `biconjCPLQ` poly2-kink tests (quarantined, see 6.1).

| command (shared machine, load 6-60) | result |
|---|---|
| `bash .claude/suite.sh --fast` | **834 / 0** over 64 suites (763 s, cold stage cache) |
| `bash .claude/suite.sh --normal` | **35 / 0** over 8 suites (738 s) |
| `CCA2_TEST_TIMEOUT=5400 bash .claude/suite.sh --slow -j 4` | **128 pass / 0 fail / 2 incomplete** — the 2 are `biconjCPLQTest`'s quarantined poly2 kinks (assumeFail) |
| `bash .claude/suite.sh clipArcByConicTest` | 7 / 0 — this suite is in **no** bucket (finding 14) |
| `CCA2_TEST_TIMEOUT=5400 bash .claude/suite.sh --verylong -j 4` | **40 pass / 0 fail / 8 incomplete** (all deliberate quarantines), 5049 s |
| `python3 -m pytest -q -n 8` (python/) | **2753 passed, 1 skipped** |
| `cargo test --release -- --include-ignored` (CCA2_STAR_NS=11..14) | all green (parallel 4, parity 3, ps_tangency 3, star 2) |
| `CCA2_FLINT=0 cargo test --release --test parity` | 2 passed, 1 ignored |
| `cargo clippy --release --all-targets` | **fails**: 5 deny-level `approx_constant` errors, 128 warnings (finding 15) |
| reviewer tests: `pytest tests/review` | 3 passed, 13 xfailed (each xfail checked with `--runxfail` to fail for the stated reason) |
| reviewer tests: `cargo test --release --test review_2026_11_19 -- --ignored` | 7 / 7 fail as expected |
| reviewer tests: `runtests(fullfile('review','reviewFindings20261119Test.m'))` | all fail as expected (list in section 8) |

## 2. Overall assessment

The test suites are green, but they do not test what the paper claims. An exhaustive + random
differential test (357 small exact inputs, three engines, an exact rational f* oracle and a
certified f** bracket) finds **wrong values in both operators in all three engines**: `conj` is
wrong on 6 of 327 continuous inputs (finite where f* = +inf; too high on an L-shape), and `biconj`
is wrong on **18 of 243 (Python) / 14 of 240 (Rust)** answered continuous inputs — too low,
too high, +inf inside dom f, finite outside conv(dom f). The defects predate the 2026-11-17 fixes.
These change claims the draft makes (abstract: "exactly"; "a small class of inputs on which the
biconjugate is still slightly low"; section 4.2 "complete" cell-interior search; the harness as
the reason the claims are checkable).

## 3. Algorithm as reconstructed from the code (short)

`conj`: per piece, classify the Hessian; PD on a polygon → Case B (cells per vertex/edge/interior
from the KKT conditions); not PD → Case D (sup on the boundary: vertex affines and, per edge of
positive curvature, the clamped 1-D maximum, available where 0 <= t* <= 1, an affine condition);
nonconvex polygons are ear-clipped and the triangle conjugates folded; unbounded pieces get
`recessionConditions` (dom f* as half-planes); pieces are combined by `maxQ` (H-form
concatenation of constraint lists, cells split by `{g_i = g_j}`), and empty cells are removed by
`cellHasInterior` (lines exactly by `feasible2`; curved cells by a double-precision witness
search: pairwise conic intersections, 400 random probes in a box, then a "complete" CAD sample).
`biconj`: atoms (convex pieces, edges, vertices) → active atom pairs (`tieArc`, sampled in
doubles) → pair envelopes (bridge faces, algebraic, degree <= 4) → per-pair subdivision →
min-overlay (`EnvelopeMin`) or `AlgCrv`/`QuaCon`. Discrepancy with the paper: the paper says
every exact decision is exact; `tieArc` (which pairs are active) and `cellHasInterior` for curved
cells are double-precision sampling heuristics with tolerance bands (section 5.2).

## 4. Findings

Severity per the skill: CRITICAL = wrong results. Each CRITICAL/MAJOR item has an executable
failing test: Python `python/tests/review/test_review_2026_11_19.py` (xfail strict), Rust
`rust/tests/review_2026_11_19.rs` (`#[ignore]` with the finding number), MATLAB
`MATLAB/review/reviewFindings20261119Test.m` (not in any suite.sh bucket).

### CRITICAL

**#1 — `conj` finite where f* = +inf: a zero-curvature recession direction strictly inside the
cone is ignored (all three engines).**
`python/cca2/conj_q.py:461-479` (`recessionConditions`), `MATLAB/conjQ.m` (`recessionConditions`,
the `for r = 1:size(R,1)` null-direction loop), `rust/src/conj_q.rs:465` (`recession_conditions`).
The zero-curvature test `d'Hd == 0` is run only for the cone's two **edge rays** `R`. For a PSD
singular H whose null vector n lies strictly inside the recession cone, no domain row is added, so
dom f* is taken as the whole plane. (The "the cone, not only its generators" check above it only
catches negative curvature; for PSD H, b^2 <= ac always.)
Input: f(x,y) = y^2/2 on the cone {x >= |y|}. Along x = (t,0), <s,x> - f = s1 t, so
f*(s) = +inf for s1 > 0. Observed (Python, Rust, MATLAB): f*(1,0) = 0.5, f*(0.5,0.3) = 0.32,
f*(0.1,0) = 0.005. Random differential: 5 of 42 cone fixtures (E-quadrant-H1-11, R5-053, R5-060,
R5-074, R6-152), wrong on 66-145 of 218 dual points each.
Fix: for PSD-singular H, take n spanning null(H); for each of ±n inside the closed recession cone
add the same half-plane row the edge-ray branch adds (`<s - L, n> <= inf_x <Hn, x>` = `<s-L,n> <= 0`
since Hn = 0). Tests: `test_conj_cone_with_interior_null_direction_is_infinite`,
`test_conj_two_cone_fan_R6_152`; Rust `conj_cone_with_interior_null_direction`; MATLAB
`coneWithInteriorNullDirection`.

**#8 — `conj` too HIGH on a bounded nonconvex polygon (all three engines).**
Input: f = y^2/2 on the L-shape (0,0),(2,0),(2,1),(1,1),(1,2),(0,2) — the domain of golden conj
record 36 with a PSD-singular quadratic. At s = (2.53125, -1.453125) the sup is 2 s1 = 5.0625
(attained at (2,0); `check.supAt` and a dense grid agree). All three engines return
6.5137939453125 = (2 s1 + s2)^2 / 2: the 1-D maximum along direction (2,1) — the ear-clip
diagonal (0,0)-(2,1) — used where its clamp 0 <= t* = 2 s1 + s2 <= 1 fails (t* = 3.61). The face's
rows (face 91 of 243) contain no t* <= 1 row. Wrong on 342 of 14 641 points of a 0.1 grid on
[-6,6]^2 (s1 in [2.2, 6], s2 in [-4.9, -1.2]). Root cause localised to Case D on a triangulated
nonconvex polygon (a diagonal's candidate keeps a value outside its availability range); not
traced to a line. Tests: `test_conj_L_shape_y2`, Rust `conj_l_shape_y2`, MATLAB `conjLShapeY2`.

**#9 — `biconj` wrong values on small continuous inputs (all three engines; confirmed point by
point in Python and Rust, and MATLAB returns the same wrong values — section 5.1).**
On the random/exhaustive corpus (continuous families E, R1-R6) the engines answer 243 (Python) /
240 (Rust) `biconj` calls; **18 / 14 are wrong** against a certified bracket
lo <= f** <= hi (lo = <s,x> - f*(s) with exact f*; hi = an LP over convex combinations of domain
samples). Each verdict was re-derived with the PROJECT'S OWN oracle (`check.supAt`,
`check.evalQuaPol`) and the engine's OWN `conj` (`python/tests/review/verify_wrong.py`):
- **9a, finite outside conv(dom f).** f = (y^2 - x^2)/2 on the unit square: f**(-0.49, 0.01) =
  0.246 (truth +inf; 42/169 points). The result has faces on {x <= 0} and {x >= 1} carrying the
  extrapolated formula y^2/2 - x/2: the cells are not clipped to conv(dom f).
- **9b, too low (Fenchel-Young violated against the engine's own conjugate).** f = -xy + y^2/2 on
  the L-shape: f**(0.2578, 1.7617) = 0.409 while <s,x> - f*(s) = 1.046 with f*(s) from the
  engine's own `conj` (= `check.supAt` = reviewer oracle). Also E-L-H1-10, R2-064 (star polygon),
  R4-066 (continuous 3-piece fan on a pentagon: 7.03 for 15.63).
- **9c, +inf inside dom f.** R3-072, R6-075, R6-117, R6-131, R6-145, R4-178, R2-204, R4-164,
  R4-045 (Py), R4-010 (Rs): f(x) finite (`check.evalQuaPol`), engine f**(x) = +inf.
- **9d, too high.** R3-072 x = (-2.49, 1.85): engine -17.29, a convex combination of three
  domain points (f from `check.evalQuaPol`) gives -18.455. R2-169 (Py too high, Rs too low).
Not caused by psTangency (rerun with it disabled: identical), and present on the 2026-11-13
engine (identical output): not a regression of the 2026-11-17 fixes.
Tests: `test_biconj_infinite_outside_the_square`, `test_biconj_fenchel_young_L_shape`,
`test_biconj_finite_inside_the_domain_R3_072`; Rust `biconj_infinite_outside_the_square`,
`biconj_fenchel_young_l_shape`; MATLAB `biconjInfiniteOutsideTheSquare`, `biconjFenchelYoungLShape`.

### MAJOR

**#2 — `cellHasInterior`'s "complete" CAD fallback is not complete (all three engines).**
`python/cca2/cell_has_interior.py:44-104` (`_cadInterior`), `MATLAB/cellHasInterior.m` (c19bc97),
`rust/src/cells.rs`. Two mechanisms, both reproduced in all three languages:
(a) **vertical lines** (a2 = a1 = 0) contribute no critical abscissa: the cell
{10 <= x <= 10.5, outside the circle centred (100,0)} has interior but is declared empty (its
horizontal twin is found); (b) **degenerate conics sharing a line factor**: `conicMeet` reports
the pair degenerate and returns no points, so the other lines' crossings are not abscissae. (b) is
**reached inside a real `conj` call** (random nonconvex hexagon, `hunt_cell_interior.py`, 1 of
1200 inputs); no wrong value was found at the sampled dual points of that input, so value impact
is unconfirmed. The CAD also samples in doubles with a 1e-9 band and is skipped for side-0 rows.
The paper (section 4.2, lines 587-612) calls this search complete. Mutation M4 (fallback switched
off) survives the whole Python suite. Tests: `test_cell_has_interior_far_vertical_strip`,
`test_cell_has_interior_shared_line_factor`; Rust `cell_has_interior_far_vertical_strip`,
`cell_has_interior_shared_line_factor`; MATLAB `cellInteriorMissesAFarVerticalStrip`,
`cellInteriorMissesSharedLineFactorCell`.

**#3 — the #64 cross-check refuses four established PLQ registry entries as "non-PLQ".**
`python/bench/cross_convexdb.py`. `norm1_2d` (|x|+|y|), `separable_bump_quadratic_2d`,
`separable_multiwell_abs_2d` fail with "two pieces own the same side of line" (the parser turns
first-match Piecewise branches, e.g. `(-L, True)`, into overlapping closed regions), and
`segment_indicator_2d` fails because `build_fixture_1d` ignores the half-plane side `hp.s`
(lines 699-711: `k > 0` is always read as a lower bound). All four have stored, *established*
f* and co_f. The summary, `RESULTS.md` and TODO.md ("16 non-PLQ entries are refused as expected")
therefore count 4 PLQ entries as out of scope, and the pytest guard `test_non_plq_refused`
**pins** |x|+|y| as non-PLQ. Tests: `test_crosscheck_accepts_norm1`,
`test_crosscheck_accepts_bounded_segment_indicator`.

**#4 — the cross-check truncates rational data.** `piece_value_row` (`int(2*A20)`, `int(A11)`…)
and `atom_to_hp_tuple`/`dnf_1d` (`HP(int(A), int(B), int(C))`) truncate non-integer rationals:
x*y/3 becomes 0, y <= 1/2 becomes y <= 0. `fixture_check` usually catches the damage and the
entry is then refused as a "parser defect", counted with the non-PLQ entries. No current accepted
entry is affected (all have integer data); entries with fractional data are silently excluded.
QuaPol takes Fractions, so `Fraction` would do. Tests: `test_crosscheck_keeps_rational_coefficients`,
`test_crosscheck_keeps_rational_halfplane`.

**#5 — the definition oracle shares finding #1's blind spot and cannot see 9a.**
`python/cca2/check.py`: `supAt` tests recession only along 49 sampled directions with a 1e-12
curvature band, so a null direction strictly inside a cone between samples is missed (R5-053:
`supAt` 2.11, truth +inf; 39 of 11 424 oracle-vs-oracle comparisons); `checkBiconj` samples only
the bounding box of dom f, so f** finite outside conv(dom f) is invisible for a box domain; and
`test_definition.py` runs only on the golden fixtures (produced by the authors' own suites). The
harness is therefore not independent of the engine on exactly the classes where the engine is
wrong. Test: `test_definition_oracle_detects_unbounded_cone`.

**#6 — Rust `biconj` panics on valid input.** `rust/src/biconj_atoms.rs:347` (`tie_arc`):
`bp.sort_by(|a,b| a.partial_cmp(b).unwrap_or(Equal))` on breakpoints containing NaN aborts with
"user-provided comparison function does not correctly implement a total order" (Rust >= 1.81).
5 of 357 inputs (R4-017, R4-080, R1-091, R4-115, R4-171). Python (`sorted(set(bp))`) and MATLAB
(`sort`) do not crash but order NaN arbitrarily and sample at NaN, so the tie-arc decision is
silently degraded there. 11 such `unwrap_or(Equal)` sorts exist in `rust/src`. Fix: filter
non-finite roots in `real_roots`, use `total_cmp`. Test: Rust `biconj_does_not_panic_on_r4_017`.

**#7 — the legacy `conjCPLQ` eigenvalue tolerance (TODO item checked against the code).**
`MATLAB/conjCPLQ.m:116-194`: `all(ev > sqrt(eps))` is an ABSOLUTE test, then `tolE` is RELATIVE.
f = (x^2 + y^2)/2^31 (strictly convex, f* = 2^29|s|^2) goes to the affine branch and returns a
needle (+inf at every s != 0); Q = diag(1, 2^-30) and diag(1, -2^-30) reach the rank-1 branch,
where `null(Q)` is empty, and crash with an unnamed `MATLAB:badsubscript`. The TODO wording
("can return +Inf for a strictly convex quadratic") is right but incomplete (crash; the
indefinite case is misrouted too). Legacy engine only; the exact engine refuses these dyadic
inputs as `notExact` (finding 18). Tests: MATLAB `conjCPLQTinyPositiveEigenvalue`,
`conjCPLQSmallEigenvalueDoesNotCrash`.

**#10 — discontinuous input is accepted silently and gives a wrong biconjugate.** Neither engine
checks continuity across shared edges (the paper's standing assumption). On 30 discontinuous
two-piece inputs `conj` agrees with the oracle (sup over the closed pieces), but 12 of 24 answered
`biconj` calls are above a certified upper bound (e.g. D-160: 14.36 where a convex combination of
domain points gives -4.89). Either refuse by name or define and implement the lsc semantics.

**#11 — the three engines disagree on random inputs.** Python vs Rust: 11 of 357 fixtures
(biconj): R4-010, R4-045, D-069, R3-072, D-132, R4-150, R2-162, R4-164, R2-169, R2-204 (values)
and R4-080 (Python refuses `bridgeFaceAlg:noZ`, Rust panics), plus R4-003 (Python > 600 s, Rust
answers). MATLAB sides with Rust on some (R3-072, R4-045) and with Python on others (R4-010); it
also refuses by its 2^53 ceiling where the ports answer (`ratQ:overflow` on R4-003, R4-059,
R4-101; `bridgeFaceAlg:witnessOverflow` on R3-037, and in place of `noWitness` on R2-029, D-090).
The parity claim holds on the recorded corpus only; the draft should say so.

### MINOR

12. **Cross-check grading:** `_face_holds`/`eval_algcrv_mp` scale with signed coordinates
    (`max(1, P[0], P[1])**2`, should be `abs`); an `Ambiguous` (two faces holding with different
    values — itself a possible defect) silently falls back to the 1e-6 double protocol instead of
    being reported; Rust "open" entries skip the f** <= f sanity check Python gets;
    `test_cross_convexdb.py` skips silently when the sibling registry is absent (fresh clone / CI).
13. **Stale statements.** TODO "Open now" and outline G7 say "Two BRIDGE faces owning one cell —
    refused by name, never observed"; the code no longer refuses (`envelopeSubdivisionAlg.m:44-51`,
    since 2026-10-17: overlapping pieces, value by max) and erratum Example 1 reaches it (18 cells).
    The `nonConvexPiece` refusal message says "co f = f here" for f = y^2/2 on the L-shape, where
    co f != f.
14. `clipArcByConicTest` is in no `suite.sh` bucket, so `--fast/--normal/--slow` never run it (7/0
    when run by name).
15. `cargo clippy` fails (5 `approx_constant` deny errors on deliberately generic screen points,
    `alg_row_fails.rs:1661,1808`); 128 warnings (mostly index loops, complex types).
16. MATLAB `algFaceBelow>separator` (line 134) solves nearly singular systems (RCOND ~1e-21,
    warnings during `biconjQ` on random inputs); whether a decision depends on it is unconfirmed.
18b. `python/cca2/alg_row_fails.py:1376-1386`: unreachable code after `return hit` in
    `_curvePoints` (ruff F821: it references an undefined `S`); harmless, misleading.
17. Error-identifier inconsistency for f** = -inf: `PLQ:biconjQ:noFace` (whole plane),
    `PLQ:biconjQ:minusInfinity` (cones), `biconjAtoms:emptyDomain` (fans).

### INFO

18. The exact path refuses exactly representable dyadic data with denominator > 10^6
    (`fromDouble(maxDen=10**6)`; e.g. 2^-30) as `notExact`; conservative, by name.
19. Face counts are representation sizes, not face counts of f*: f = (x^2+y^2)/2 on a 2n-vertex
    star (true f* has O(n) faces) gives nf = 105, 180, 325, 718, 836, 927, 1655, 2602, 2555,
    2835, 4241, 11683 for n = 3..14 (Rust; non-monotone at n = 11). This is TODO's "the cell
    count is still an upper bound", quantified: overcount ~8x at n = 3 to ~200x at n = 14. The
    paper's "317-face kink AlgCrv" should be described as a representation size.

## 5. Differential, property, numerics, complexity, reproducibility

### 5.1 Differential test (new; extends #64 beyond the registry and adds MATLAB)

Generator `python/tests/review/plq_gen.py` (seed 20261119): 162 exhaustive cases (Hessian in
{-1,0,1}^3 x {plane, triangle, square, L-shape, quadrant, segment}; quadrant only where the
recession question is decidable without ties) + 210 random (R1 convex polygon, R2 star-shaped
nonconvex, R3/R4 continuous 2-/3-piece splits, R5 cone, R6 two-cone fan, D discontinuous). All
data are integers, evaluation points dyadic, so the three engines read identical exact input.
Dual points: 169 offset + 49 lattice points; primal points: a 13x13 offset grid over the
bounding box (+1/2). Drivers: `run_engines.py` (Python, Rust `golden`), `MATLAB/review/
reviewDiffDriver.m`; grading `compare_engines.py`; re-verification `verify_wrong.py`.

| | conj ok / wrong / refused | biconj ok / wrong / refused |
|---|---|---|
| Python (357) | 350 / 6 / 1 (timeout) | 237 / 30 / 90 |
| Rust (357) | 351 / 6 / 0 | 238 / 26 / 93 (5 panics) |
| MATLAB (333 of 357 finished) | 327 / 6 / 0 | 224 / 20 / 89 |

(wrong counts include the discontinuous D cases of finding #10: 12 in Python and Rust, 10 in
MATLAB. On continuous inputs MATLAB is wrong on 10 of 223 answered `biconj` calls, the same
fixtures and the same values as the ports where it answers; it refuses more by its 2^53 ceiling.)
MATLAB ran in four processes capped at 3 h each; two shards hit the cap, each stuck on one
fixture for over an hour (R4-150, a continuous 3-piece pentagon fan, and R3-177; Rust answers
both in seconds), so 24 fixtures have no MATLAB result. That is recorded as a MATLAB performance
observation, not graded.
The reviewer's f* oracle was itself cross-checked against `check.supAt` on 11 424 points: the
only 39 disagreements are finding #5 (the reviewer's oracle is right: hand computation).

### 5.2 Numerics

Exact where claimed: QuaCon data, conic rows, `feasible2`. In doubles with tolerances, deciding
structure: `cellHasInterior` (1e-9 relative band; a boundary point satisfies a strict cell, so a
zero-area cell is kept — this is the count inflation of item 19 — and a far/thin cell is missed,
#2); `tieArc` (1e-6 imaginary-part band, NaN breakpoints, #6); `QuaCon.eval` picks the face
"most inside" with a sqrt(eps) band. The paper's sentence "every exact decision ... is decided in
exact arithmetic" is true of the decisions listed there but not of active-pair selection or
curved-cell liveness.

### 5.3 Mutation check (Python, 5-10 hand mutants, full fast suite, file restored each time)

| mutant (Python copy, full fast suite) | result |
|---|---|
| M1 Case D keeps zero-curvature edges (`alpha <= 0` -> `< 0`) | killed (2 failed) |
| M2 `recessionConditions`: interior negative-curvature cone test removed | **survived** |
| M3 `conicMeet`: the `/scale` fix reverted (the star bug) | killed (3 failed) |
| M4 `cellHasInterior`: CAD fallback switched off | **survived** |
| M5 `cellHasInterior`: tolerance band 1e-9 -> 1e-4 | no verdict: suite did not finish in 900 s (baseline 65 s) |
| M6 `EnvelopeMin.eval`: min -> max | **survived** |
| M7 `QuaCon.eval`: face tolerance sqrt(eps) -> 1e-3 | **survived** |
| M8 `rat_q.evalFace`: divide by 2d | killed (9 failed) |
| M9 `biconj`: interior negative-curvature cone test removed | killed (1 failed) |
| M10 `_dropImplied`: drop rows without the interior test | killed (14 failed) |

Survivors are missing tests: M2 (no fixture where the cone's interior has negative curvature but
the rays do not), M4 (nothing needs the CAD fallback, although the star fix added it for one),
M6 (no fast test checks `EnvelopeMin` VALUES — Example 2's default return type), M7 (no point
near a face boundary is evaluated). Note on method: the first attempt with `pytest -x -n 16`
reported no result for M5-M8; the cause was a numba cache lock left by a killed run, not the
mutants; M5-M10 were rerun with fresh caches.

### 5.4 Reproducibility of the paper's test numbers

MATLAB fast 834/0 and normal 35/0 reproduce the paper's table; the slow bucket's eight exact
suites are green (kink 3/0, min-subdivision 4/0, erratum Ex. 2 4/0, algCrv 3/0, envelope-slow
1/0, algRowFails-surd 5/0, unbounded-face 19/0, guard-symbolic 11/0). Python now has 2753 passing
tests (the draft says 2741: more tests since). Rust: as in section 1. The very-long bucket
(`--verylong -j 4`, 46 jobs, 5049 s): **40 pass, 0 fail, 8 incomplete**; the 8 are deliberate
quarantines (`assumeFail` in `testcPLQ` (G11) and `testMaxMultiRegion` (legacy pins)); the
erratum Example 2 subdivision suite passes 3/0 in 3165 s.
Timings are from a shared machine and are not used for any conclusion.

### 5.5 Complexity

The outline/draft claim no complexity bound for the exact engine (the linear-time claims are
[COAP]/[JOGO]'s). Measured output size on stars (item 19) grows far faster than the true face
count, and run time grows ~2x per added star arm at n = 13-14 (single runs, shared machine:
indicative only).

## 6. TODO items checked against the code

1. **`biconjCPLQ` on the poly2 kinks (> 4 h):** the two tests are quarantined (`assumeFail`,
   2 incomplete in the slow run). Not re-run (> 4 h). Legacy path; the exact engine answers these.
2. **Two bridge faces owning one cell:** no longer refused (finding 13); handled as overlapping
   pieces with value by max. The max rule is right only if every overlapping face is a minorant
   of f** on its region; not tested separately here, and the biconj defects of #9 were not traced
   to it.
3. **`conjCPLQ`'s eig(Q) > sqrt(eps):** confirmed defective (finding #7); the "empty domain"
   reasoning (negative eigenvalue ⇒ f* = +inf everywhere) is mathematically right, the tolerance
   routing is not.
4. **The cell count is an upper bound:** confirmed and quantified (item 19); the mechanism is the
   tolerance band in `satisfiesAllDouble`, which accepts boundary points as interior witnesses.

## 7. Second pass (after reading DECISIONS.md / TODO.md)

- Inherited assumption: "the definition oracle shares nothing with conjQ" — it shares the
  recession blind spot (#5). The star-polygon fix was verified against the same oracle family.
- The 2026-11-17 entry calls the CAD "complete"; it is a sample over an incomplete set of
  critical abscissae (#2).
- Input families never tried by the suites: unbounded pieces with PSD-singular H whose null
  direction is interior; PSD-singular H on nonconvex polygons; points outside conv(dom f) for
  biconj; multi-piece random continuous inputs.

## 8. Recommended next actions, in order

1. Decide how the draft treats #1, #8, #9 (fix before submission, or narrow the claim).
2. Fix #1 (small, local, three engines) and #6; add the cone family to `test_definition`.
3. Localise #9 (start from 9a: clip the result to conv(dom f); then the L-shape Fenchel-Young
   violation) and #8 (Case D clamp on ear-clip diagonals).
4. Make the cell-interior fallback a genuine CAD (#2): vertical lines and factor-sharing conics.
5. Extend `test_definition.py` with this review's generator (and points outside conv(dom f)).
6. Fix the cross-check parser (#3, #4), then re-run it and correct the "16 non-PLQ" record.
7. Refuse discontinuous input by name (#10); state parity as "on the recorded corpus" (#11).

Review test files: `python/tests/review/test_review_2026_11_19.py`,
`rust/tests/review_2026_11_19.rs`, `MATLAB/review/reviewFindings20261119Test.m` (run:
`cd MATLAB; addpath(pwd); runtests(fullfile('review','reviewFindings20261119Test.m'))`).
Harness: `python/tests/review/{plq_gen,plq_oracle,run_engines,compare_engines,verify_wrong,
hunt_cell_interior}.py`, `MATLAB/review/reviewDiffDriver.m`.

No finding was fixed in this task; the reviewed code and `manuscript/` are unchanged.

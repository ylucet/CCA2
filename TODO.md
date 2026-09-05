# TODO

## 2026-09-03 — Target venue decided: Mathematical Programming Computation

Author decision, 2026-09-03. CCA2's paper targets **Mathematical
Programming Computation** (backup: Computational Optimization and
Applications; SVVA if the paper ends up led by the mathematics rather
than the algorithms).

**Why MPC:** of the four manuscripts in the group — `FLTW`, `SCATpy`,
`convexdb`, `CCA2` — CCA2 is the only one that really is a *solver*, and
MPC exists for practical computation in mathematical optimization. Its
mandatory code-and-data submission with independent verification of
computational results suits exact closed-form conjugacy for bivariate PLQ
better than any other venue on the list. All four had been aimed at ACM
TOMS; they were split across communities because the same five to ten
people would referee all four whatever the masthead. FLTW keeps TOMS,
SCATpy moved to the Journal of Symbolic Computation. Full reasoning:
`../convexdb/DECISIONS.md` 2026-09-03 and
`../convexdb/manuscript/PLAN.md` "Cross-project allocation".

**Compliance and cost:** free. Springer's 12-month green route
(self-archive the accepted manuscript) clears the tri-agency bar at no
APC; gold is ~$3,290 and unnecessary.

**Two things to handle deliberately when the paper is written:**

- [ ] **The MATLAB + Symbolic Math Toolbox dependency is real friction at
  MPC**, whose results verification assumes a referee can run the code —
  and this one needs a licence on the UBC VPN. **Declare it in the cover
  letter** and offer MPC's "detailed observation / guest access" route
  rather than letting a referee discover it. COAP has no code mandate if
  it turns out fatal.
- [ ] **Cite `convexdb` by its Zenodo DOI**, not a repo URL, wherever this
  project uses that registry as ground truth. The DOI is deliberately
  being minted before any of the four papers is submitted, so the citation
  does not depend on the convexdb paper's own fate.


## 2026-09-04 (audit, closed) — non-convex faces are SPLIT, so the audit's silent defect is gone

The audit below found two defects; the crash was fixed there and the non-convex face was turned
into a named refusal. It is now ANSWERED instead.

**conj turns a union into a MAX**, so for `P = union of P_i` we have `f* = max_i (q + I_{P_i})*`:
triangulating a non-convex face and folding the pieces with `maxQ` is EXACT, not an approximation.
The triangulation is exact too -- **ear clipping**, whose every decision is the sign of an integer
cross product (does this corner turn left; does another corner lie inside this candidate triangle).
A floating-point triangulator could return a different mesh for the same polygon depending on
rounding, and each triangle's conjugate would then be exact about the wrong thing.

Two details that had to be right:
- **Collinear vertices are dropped first.** A corner that turns neither way lies on the segment
  between its neighbours and changes nothing about the region; kept, it becomes a zero-area ear.
- **The boundary cycle is walked from the EDGES and oriented by its own signed area**, not taken
  from `P`, so no clockwise/counter-clockwise convention is consulted -- reading one backwards
  would reverse the orientation the ear test depends on.

Verified two-sided on the L-shaped face against a sup sampled over the TRUE L on a 1201^2 grid:
**0 below** it and at most **5.8e-07 above**, which is the grid's own resolution error (h = 1.7e-3,
so h^2/2 ~ 1.4e-6). The fold is convex, as it must be. The answer is FAT -- 133 faces for four
triangles folded pairwise -- which is the known cell-count item, not a correctness one.

**Still refused, and it is a separate case:** a non-convex piece that is also UNBOUNDED. Ear
clipping needs a closed boundary, so splitting one means cutting along its recession directions
first.

**One more oracle bug, the third of this kind today, now fixed in the checked-in probe.**
`stress-probe.m` sampled convex combinations of the vertices -- i.e. the convex HULL -- which for a
non-convex face includes points OUTSIDE the domain, so its "sampled sup" was over too large a set
and a correct `f*` looked too low. That artefact produced the original 44-of-201 and a residual
1-of-201. It now filters samples through `inpolygon`, which is floating point and fine in an oracle:
it only has to produce points that really are in the domain.

## 2026-09-04 (audit) — "can it handle ANY QuaPol?" NO. The table was an ENUMERATION, not a proof.

The coverage table varies DOMAIN SHAPE x HESSIAN CLASS and holds everything else fixed. Asked
whether `conjQ` handles any input at all, the honest move was to vary the other axes --
`.claude/stress-probe.m` does: face convexity, vertex count, piece count, collinear vertices,
coefficient magnitude, and vertex denominators. It found **two real defects the table could not
see**, one of them silent.

**1. A NON-CONVEX face was answered WRONGLY, and silently.** A piece is described here as an
intersection of half-planes, which cannot express a REFLEX corner: for an L-shaped face the
half-plane intersection is strictly SMALLER than the L, so the sup was taken over the wrong set and
came out BELOW the true one at **54 of 303** dual points. Wrong in the direction a one-sided oracle
treats as impossible, so the existing sweeps would never have caught it -- they all use convex
polygons. NOW REFUSED by name (`PLQ:conjQ:nonConvexPiece`), detected exactly while building the
half-planes, with the remedy in the message: split the face into convex pieces, since `conj` turns
a union into a MAX and `maxQ` recombines them. **The mesh already knew** -- `QuaPol.orderEdges`
reports `isConvex = 0` for such a face.

**2. A FOUR-PIECE input CRASHED**, `MATLAB:badsubscript`, an unnamed crash rather than a named
refusal. The corner-naming loop in `assembleQuaConCells` used the PRE-merge cell count while
reading the POST-merge list. It could not fire until `mergeAdjacentCells` actually started removing
cells -- which it never did until Case D was restructured, and the merge is recorded above as having
fired ZERO times before that. Fixed; the four-quadrant fan now conjugates at 0 wrong of 401 against
a closed form.

**What the probe found CLEAN**, so these are now known rather than assumed: a 12-gon, collinear
vertices on an edge, coefficients up to 1e7 (no `ratQ` overflow), and vertices with denominators
7, 11 and 13 (no denominator blow-up).

**So the accurate claim is:** `conjQ` handles any QuaPol whose pieces are CONVEX polyhedra -- which
is what this toolbox's subdivisions are, since a face is read as an intersection of half-planes
throughout (`ALGORITHM.md` says so of `eval` too) -- and refuses a non-convex piece by name rather
than answering it. Splitting non-convex faces automatically is the obvious next item and is not
built.

**Method note, and it is the general lesson.** A coverage table is an enumeration along the axes
its author thought of. Both defects lived on axes the table held FIXED, and one was silent. Widening
the axes is a different activity from filling in the table, and it is the one that found the bugs.

## 2026-09-04 (closing) — the `conjQ` coverage table is COMPLETE. Every entry is OK.

    domain \ Hessian    PD    PSD-sing  indefinite  ND      NSD-sing  affine
    full plane          OK    OK        OK          OK      OK        OK
    bounded triangle    OK    OK        OK          OK      OK        OK
    bounded square      OK    OK        OK          OK      OK        OK
    unbounded wedge     OK    OK        OK          OK      OK        OK
    half-strip          OK    OK        OK          OK      OK        OK
    needle (dim 0) OK   segment (dim 1) OK   multi-face bounded OK   multi-face unbounded OK

`.claude/coverage-probe.m` reproduces it. The only inputs that still raise are the two the toolbox
declines BY DESIGN and neither is a gap: a CUBIC numerator (`assertOperable`, out of scope since
`CLAUDE.md` says the design target is quadratic) and an INEXACT input (`PLQ:QuaPol:notExact`,
refused deliberately because exact arithmetic on rounded coefficients yields exactly the wrong
number).

**(3) The thin dual domain is now BUILT**, and it cost one value in an existing column exactly as
predicted -- no new class, and `AlgAlg` remains unrelated (that type is for irrational
COEFFICIENTS, not thin domains). `FC`'s side takes `0` meaning ON the curve:
- **q AFFINE on the plane** -> `dom f*` is the single POINT `s = L`, value `-kappa`: two equality
  rows, one constant face.
- **q PSD-SINGULAR on the plane** -> `dom f*` is a LINE. With `H = lambda*w*w'` the objective is a
  concave parabola along `w` and LINEAR along `w`-perp, so the sup is `+inf` unless `u = fD*s - Ln`
  is parallel to `w`; on that line `f* = (u.w)^2/(2 lambda) - kappa`. Stored as one equality row
  (`u . null(H) = 0`) and a quadratic face, using `pinv(H) = m m'/(m'Hm)` for any `m` spanning the
  range -- all rational.

**Three touch points, as forecast**: `QuaCon`'s constructor validation (sides are now `-1/0/+1`),
`eval` (a side of 0 tests `|q| <= tol` instead of `sign*q >= -tol`), and the producer.

**One bug found while building it, and it is a good illustration of why the filters are separate
from the geometry.** `dropRedundantRows` reasons about HALF-PLANES: it asks whether the other rows
make a row impossible to violate. Handing it an EQUALITY turns the row into `0*c`, i.e. all-zero,
which `ratQ.feasible2` correctly reads as having no interior -- so every equality was deleted as
"redundant", the thin face silently became the whole plane, and the conjugate of an affine function
evaluated to `-kappa` EVERYWHERE. Thin cells are now skipped by that pass and by the
nonempty-interior filter in `assembleQuaConCells`, both for the same reason: a thin cell has no
interior BY CONSTRUCTION, so a filter that deletes interior-less cells deletes exactly the faces
this case exists to build.

**The caveat stands and is worth repeating**: evaluating a line- or point-supported function in
floating point answers `+inf` almost surely, because a random double is never exactly on the line.
The mesh is exactly right and every predicate that built it is exact; it is simply not an object to
probe interactively. `QuaCon.eval` says so at the point where it tests an equality.

## 2026-09-04 (final) — items (1) and (2) DONE. Only the THIN dual domain is left.

    domain \ Hessian    PD    PSD-sing  indefinite  ND      NSD-sing  affine
    full plane          OK    LINE      OK nf=0     OK nf=0 OK nf=0   POINT
    bounded triangle    OK    OK        OK          OK      OK        OK
    bounded square      OK    OK        OK          OK      OK        OK
    unbounded wedge     OK    OK        OK nf=0     OK nf=0 OK nf=0   OK
    half-strip          OK    OK        OK nf=0     OK nf=0 OK        OK
    needle (dim 0)  OK nf=1     segment (dim 1)  OK nf=2
    multi-face bounded OK    multi-face unbounded OK

**(1) The empty dual domain now RETURNS the zero-face QuaCon** instead of raising. `nf = 0` is how
"+infinity everywhere" is spelled, and it needed no change to the type -- only the willingness to
return it. `assembleQuaConCells` and `conjPieceQ` both do. `caseAFullDomain` was raising one
identifier for two different situations and now separates them: `~isPSD2(H)` means some direction
has negative curvature, so `f*` is +infinity everywhere -- three of the five non-PD full-plane
cells -- while a PSD-singular H leaves a genuinely thin domain.

**(2) A domain of dimension < 2 is now READ.** The conjugate of such an input is FULL-dimensional,
so this needed no storage work either: a needle at `p` with value `c` gives `<s,p> - c`, affine on
all of `R^2`; a segment gives the clamped one-dimensional maximum, which is exactly `caseD`'s
candidate set. The only new code is `degenerateShape`, because such a mesh has `nf = 0` and an empty
`F`, so there is no face to read the geometry from and every edge belongs to the single piece.
`caseD` is used unconditionally there: the problem is ONE-dimensional, so what decides its shape is
the curvature ALONG the edge, `d'Hd`, not the definiteness of H in the plane. A chain carrying more
than one function is refused by name (`ambiguousChain`) since `F` cannot say which edge takes which.

**LEFT: the thin dual domain, two cells, and it now says which kind it is** --
`PLQ:conjQ:domainIsAPoint` (q affine on the plane: `dom f* = {L}`) and `PLQ:conjQ:domainIsALine`
(q PSD-singular: the sup is finite only for s in the range of H). Storing either needs the H-form's
side column to carry `0` for "on the curve" -- one value in an existing field, three touch points
(`QuaCon`'s constructor validation, `eval`, the producer), and NOT a new class. Worth doing when
something consumes it; the caveat stands that evaluating a line-supported function in floating point
answers +infinity almost surely.

## 2026-09-04 (later still) — what `dim<2` actually needs: almost nothing, and NOT a new type

Asked directly, so measured. The label covered THREE different situations and they need three
different (small) answers.

**1. `dom f*` is EMPTY -- needs NOTHING. Already representable, verified.** A `QuaCon` with ZERO
faces constructs today and evaluates to `+inf` at every point, which IS the function `f* = +inf`:

    z = QuaCon(zeros(0,3), zeros(0,3), zeros(0,6), zeros(0,10), zeros(0,1), zeros(0,2), {});
    z.eval([0 0; 1 2; -3 4])   ->   Inf Inf Inf

So every `emptyDomain` refusal is a routine declining to return an answer it can already store.
This is THREE of the five full-plane cells too -- ND, indefinite and NSD-singular on `R^2` all give
`f* = +inf` everywhere, not a lower-dimensional domain. The only change needed is to stop raising:
`assembleQuaConCells` errors on an empty cell list, and `caseAFullDomain` raises
`notStrictlyConvex` without distinguishing which case it is in.

**2. The INPUT has dimension < 2 (needle, segment, ray) -- needs NOTHING either, and the label was
misleading.** The conjugate of a low-dimensional input is FULL-dimensional: a needle at `p` with
value `c` conjugates to `<s,p> - c`, affine and finite on all of `R^2`; a segment carrying an affine
`f` conjugates to a max of two affine functions, again finite everywhere. Both are ordinary
`QuaCon`s. What is missing is the INPUT side -- `pieceShape` looks for a face and these meshes have
`nf = 0` -- so this is reading work, not representation work.

**3. `dom f*` is a POINT or a LINE -- the only case that needs an extension, and it is one value in
an existing column.** Just two cells: a full-plane AFFINE `q` gives `dom f* = {L}`, a single point;
a full-plane PSD-SINGULAR `q` gives a line. A `QuaCon` cell is a list of sign conditions `>= 0`, and
a lower-dimensional set is the same list with some conditions as EQUALITIES -- so the side column in
`FC`, which today holds `+1` or `-1`, gains the value `0` meaning "on the curve". Three touch
points: `QuaCon`'s constructor validation (which currently insists on `+-1`), `eval` (which must
test `|q| <= tol` rather than `sign*q >= -tol`, since side 0 currently reads as vacuous), and
whatever produces such a cell. **No new class, and `AlgAlg` is unrelated to this** -- that type is
for faces whose COEFFICIENTS are irrational, not for domains of low dimension.

**The honest caveat on 3**, and it is why this is worth doing deliberately rather than by reflex:
evaluating such a function in floating point is nearly meaningless, because a random double is never
exactly on a line, so `eval` would answer `+inf` almost surely. The exact predicates would still be
right, and the object would still be the correct answer to store, print and conjugate back. That is
a real use, but it is not the interactive one.

**Order, if built:** (1) is a few lines and removes a whole column of refusals; (2) is
straightforward input handling; (3) is the only one that touches the type, and it should wait until
something actually consumes it.

## 2026-09-04 (later) — `conjQ` IS COMPLETE for every 2-D domain. The refusals left are the ANSWER's shape.

Supersedes the "NO" measured earlier the same day, which is kept below for its table. Re-run
`.claude/coverage-probe.m` rather than reasoning about what it would say.

    domain \ Hessian    PD    PSD-sing  indefinite  ND        NSD-sing  affine
    full plane          OK    dim<2     dim<2       dim<2     dim<2     dim<2
    bounded triangle    OK    OK        OK          OK        OK        OK
    bounded square      OK    OK        OK          OK        OK        OK
    unbounded wedge     OK    OK        emptyDom    emptyDom  emptyDom  OK
    half-strip          OK    OK        emptyDom    emptyDom  OK        OK

    multi-face bounded   OK        multi-face unbounded (oneNorm)  OK

**Every OK cell is verified**, not merely non-throwing: 9060 dual points over all six Hessian
classes at 0 wrong (`.claude/allH-sweep.m`), the unbounded family against a one-sided sampled sup
plus a growth test (`.claude/unbounded-sweep.m`), and the new cells against separable closed forms.

**Every `emptyDom` cell is CORRECT**, not a gap: the piece recedes along a direction of negative
curvature, so `f*` really is `+infinity` everywhere and `dom f*` is empty. The only thing missing is
somewhere to store a function with no domain -- a `QuaCon` carries at least one face. The same
representational gap `conjCPLQ` records as `conjugateHasEmptyDomain`.

**Every `dim<2` cell is likewise the ANSWER's shape**: a full-plane quadratic that is not strictly
convex has `dom f*` equal to a point (affine `q`), a line (PSD-singular) or nothing. A needle or
segment DOMAIN is the same story one level down.

**So what closed the gap**, and it was two things, both of which turned a refusal into a
computation:
1. **The piece's geometry became an EDGE LIST plus HALF-PLANES** (`pieceShape`), replacing the
   vertex-cycle walk that an unbounded face cannot supply. Nothing downstream needed the ordering:
   normal cones come from the edges incident to a vertex, cells from pairwise comparisons. A ray is
   then just an edge clamped at ONE end -- `t* >= 0` instead of `0 <= t* <= 1` -- which is one
   affine condition instead of two, so `caseB` and `caseD` generalised without new machinery.
2. **The null-recession condition is a linear program, not an obstacle.** Along a recession
   direction of zero curvature the slope is `<s-L,d> - <Hd,x>`, so finiteness is
   `<s-L,d> <= inf_P <Hd,x>` -- and that infimum over a polyhedron is minus infinity when some
   recession direction decreases it (then `dom f*` is empty) and otherwise the MINIMUM OVER THE
   VERTICES. Exact integer arithmetic throughout. This is what makes `conj(xy)` on the first
   quadrant come out as the indicator of the third quadrant instead of a refusal.

**One regression was introduced and caught by the existing tests**, worth recording because the
cause is a classic: the outward normal of an edge was found by matching the half-plane list on
PERPENDICULARITY, which is ambiguous whenever two edges are PARALLEL -- on the unit square the top
edge picked up the bottom edge's normal, and its cell came out as `s2 <= 1` where it must be
`s2 >= 1`. 74 of 307 dual points wrong, 53 in no cell. The normal is now built from the edge's own
direction and oriented against the piece's vertices and recession directions, which refers to no
other edge and so cannot confuse two parallel ones.

## 2026-09-04 — IS `conjQ` COMPLETE? NO. Measured, with the table.

`.claude/coverage-probe.m` enumerates the input space along the axes the dispatch branches on --
DOMAIN shape x HESSIAN class -- and reports OK or the refusal identifier for each cell. Re-run it
rather than reasoning about what it would say.

    domain \ Hessian    PD        PSD-sing  indefinite  ND        NSD-sing  affine
    full plane          OK        notSC     notSC       notSC     notSC     notSC
    bounded triangle    OK        OK        OK          OK        OK        OK
    bounded square      OK        OK        OK          OK        OK        OK
    unbounded wedge     unbounded unbounded unbounded   emptyDom  emptyDom  OK
    half-strip          unbounded unbounded unbounded   emptyDom  OK        OK

    multi-face bounded   OK        multi-face unbounded (oneNorm)  OK
    needle / segment     noFace    cubic  unsupportedType          inexact  notExact

**COMPLETE: a bounded polygonal piece, every Hessian class.** All six classes, verified over 9060
dual points at 0 wrong (`.claude/allH-sweep.m`). Multi-piece bounded inputs follow, since Step 3 is
a fold.

**ONE REAL ALGORITHMIC GAP, and it is not representational.** An UNBOUNDED piece whose `q` is NOT
concave or affine. Measured: `q = (x^2+y^2)/2` on the first quadrant has a FINITE conjugate at every
`s` (the maximiser is the projection of `s` onto the quadrant -- `s=(1,1)` gives 1, `s=(-3,4)` gives
8), and `conjQ` refuses it with `PLQ:conjQ:unbounded`. So a legitimate answer exists and is not
computed. This shape -- a convex quadratic on a cone -- is common, so this is the next item.
What it needs: on an unbounded piece the finiteness test is a QUADRATIC form on the recession CONE
rather than the linear one the concave case enjoys (`d'Hd >= 0` for all `d` in the cone, plus the
linear condition on the null directions). For a cone generated by `d1,d2` that is an exact integer
test -- `a,c >= 0` and (`b >= 0` or `b^2 <= ac`) with `a=d1'Hd1`, `b=d1'Hd2`, `c=d2'Hd2` -- so it is
tractable; it is simply not built.

**THE REST ARE REPRESENTATIONAL, not algorithmic** -- the answer is known and there is nowhere to
put it, because a `QuaCon` face is two-dimensional:
- FULL PLANE, `H` not positive definite: `dom f*` is a point (affine `q`), a line (PSD-singular) or
  EMPTY (ND, NSD-singular, indefinite). Raised as `notStrictlyConvex`, which names the test rather
  than the reason -- worth renaming when the dim-<2 representation lands.
- `emptyDomain` on an unbounded piece: correct, `f*` really is `+inf` everywhere.
- `noFace` on a needle or a segment: a domain of dimension < 2.

**Not gaps at all:** a cubic numerator (`assertOperable`, out of scope by design) and an inexact
input (`PLQ:QuaPol:notExact`, refused deliberately).

## 2026-09-03 — EXACT `conj` landed for Case A; one legacy defect found on the way

- [ ] **`conjCPLQ`'s `eig(Q) > sqrt(eps)` can claim an EMPTY domain for a strictly convex
      quadratic.** Reproducer, and it is two lines: `k=16384; N=k^2+1; QuaPol([1,k,N,0,0,0]).conj()`
      returns a QuaPar evaluating to `+Inf` everywhere, where `f*` is finite on all of R^2 (`det H`
      is exactly 1). NOT a refusal -- a silent wrong answer. `conjQ` gets it right, against the
      definition, to ten digits. Not fixed: every alternative tolerance is wrong in a different
      direction and the correct fix IS the exact test, i.e. `conjQ` itself; changing `conjCPLQ`'s
      branch would alter dispatch across three buckets to fix an input nothing currently produces.
      Full trace and the table: `DECISIONS.md` 2026-09-03 (second entry).

- [x] **DONE 2026-09-03: `conjQ` Case B -- a strictly convex quadratic on a bounded convex
      polygon**, exactly, as the KKT active-set subdivision (one affine cell per polygon vertex,
      one quadratic cell per edge, one interior cell). Every cell boundary is a straight line, so
      the answer is a `QuaCon` whose edge conics are all lines -- the shape the SCIP bridge wants.
      Ports `conjConvexOverPiece`'s decomposition and removes its ten `sym` calls. Verified against
      an independent candidate-enumeration maximisation: 40 random polygons, 8040 dual points,
      0 wrong, 0 uncovered, worst relative error 4.6e-15.

- [x] **DONE 2026-09-03: `conjQ` Case C -- a CONCAVE or AFFINE quadratic on a bounded polygon**,
      exactly. `<s,x> - q(x)` is then convex in x, so the max sits at an extreme point and
      `f*` is the max of the affine functions `<s,v> - q(v)`: one cell per vertex that wins, and
      no curvature anywhere in the answer. Verified against that max, written directly: 40 random
      polygons, 4221 dual points, 0 wrong, 0 uncovered, worst 1.8e-15. This is the conjugate half
      of `ALGORITHM.md`'s concave-envelope case, and the cells ARE the lower hull's normal fan --
      so `ratQ.hullQ` is not needed for it after all; the hull is implicit in the cell
      construction plus the feasibility filter.

- [x] **DONE 2026-09-03: `ratQ.feasible2`** -- exact 2-D polyhedral feasibility (and nonempty
      INTERIOR) by one Fourier-Motzkin elimination. Needed to drop cells that describe the empty
      set: a collinear vertex under a strictly concave q is DOMINATED, and its cell is empty
      rather than small. Cross-checked against `linprog` on 1939 random systems, 0 disagreements.

- [x] **DONE 2026-09-03: `maxQ` -- Step 3, the exact pointwise max of two `QuaCon`s**, and with it
      `conjQ` on MULTI-FACE inputs (`f* = max_k (f|P_k)*`, so Step 3 is a fold rather than a
      special algorithm). The H-form makes the overlay a CONCATENATION of two constraint lists --
      no geometric intersection, no vertex merging, no arrangement -- against `maxQuaPar.m`'s 4654
      lines of double geometry. This is also the first thing in the repository that PRODUCES a
      genuine elliptical edge: the split boundary is `{g_i = g_j}`, whose quadratic part is
      `(H_i - H_j)/2`. Verified against the pointwise max of its operands, and by asserting that
      every face of the fold is EXACTLY (as a rational) a face of one operand.

      Two gaps carried forward, both stated in `assembleQuaConCells` and neither able to produce a
      wrong VALUE: a cell empty only because of a CURVED constraint is not detected (so `nf` is an
      upper bound -- 38 cells where the true count is smaller), and corners involving a curved edge
      are not named. Both need Phase 2c's exact degree-4 sign kernel.

- [x] **DONE 2026-09-04: `conjQ` Case D -- EVERY quadratic that is not positive definite**, on a
      bounded polygon: indefinite, negative definite, PSD-singular, NSD-singular, affine. This
      closes the bounded-piece case analysis, and it did NOT need [COAP] A.2-A.5 or the `x*y`
      frame change. One fact replaces all of it: when H is not positive DEFINITE the sup of
      `<s,x> - q(x)` over P is attained on the BOUNDARY (an interior maximiser would force the
      objective's Hessian `-H` to be NSD, i.e. H PSD; and a PSD-SINGULAR H leaves a direction of
      zero curvature along which one can walk to the boundary without decreasing). So the answer is
      the max of the vertex affines and, per edge of positive curvature, that edge's clamped 1-D
      maximum -- folded with `maxQ`. The concave case falls out with no qualifying edge, so it is
      ONE branch, not two, and `caseCConcaveOnPolygon` became `maxOverVerticesQuaCon`.
      Verified over all six Hessian classes: 9060 dual points, 0 wrong, 0 uncovered, worst 4.6e-15
      (`.claude/allH-sweep.m`).

- [ ] **THE CELL COUNT IS STILL AN UPPER BOUND, but 4x less so.** Tracked on one PSD-singular
      pentagon, which is the worst class: **501 -> 274 -> 121 faces**, 70 of them ever occupied over
      200k samples, carrying 10 distinct functions. Values were never affected (0 wrong at every
      stage). What each step bought, because the order is the lesson:
      1. Three sound FILTERS -- infeasible linear part (`ratQ.feasible2`), contradictory sides on
         one canonical curve, and a constant-sign conic (`ratQ.conicSign`) -- took 501 to 274.
      2. Exact REDUNDANCY elimination (`dropRedundantRows`) shrank the constraint lists from 14
         rows to 6.3 and removed **not one cell**, and `mergeAdjacentCells` then fired **zero**
         times. Both are still in and both are right; they simply could not reach the problem,
         because what separated those cells was CURVED and no exact test available today can see it.
      3. RESTRUCTURING Case D took 274 to 121 and the runtime from 10 s to 4 s. That is the actual
         fix, and the two attempts above are what identified it: the pairwise `maxQ` fold was
         generating the conic splits in the first place. Every availability condition (`0 <= t* <= 1`)
         is AFFINE, so the arrangement is a LINE arrangement that `feasible2` decides exactly and
         completely; build that first, pruning as it grows, and let a conic appear only when asking
         which of a cell's few surviving candidates is largest.
      **What is left is genuinely curved** and needs Phase 2c's kernel: 121 reported against 70
      occupied, and 70 occupied carrying 10 functions. `CLAUDE.md` section 5's ladder, rungs 1-3,
      is what this sequence was.

- [x] **DONE 2026-09-04: `conjQ` Case E -- a CONCAVE or AFFINE quadratic on an UNBOUNDED piece.**
      The first case where `dom f*` is a PROPER SUBSET, which is the mathematics rather than a gap.
      With H negative semidefinite the objective already diverges along a recession direction `d`
      unless `Hd = 0`, and then `<grad q, d> = <L,d>` independently of where on the ray it is taken,
      so the finiteness condition is `<s - L, d> <= 0` -- LINEAR in `s` AND in `d`, which is why
      testing the extreme rays settles the whole cone. On that cone the objective is convex, so the
      max is at an extreme point and the value is the same vertex max as the bounded case.
      Covers the whole elementary unbounded family, `q` affine being the `H = 0` instance:
      indicators, support functions, norms, `max(0,x,y)`. Cross-checked against
      `QuaPol.oneNormConjugate` -- a fixture that predates all of this work and states independently
      that `|x|_1` conjugates to the indicator of the unit infinity-ball.

- [x] **DONE 2026-09-04: `biconjQ` -- the exact BICONJUGATE, first two cases.** DIRECT, not
      `conj(conj(f))`: `ALGORITHM.md` records that the double conjugation cost 436 s on an input
      whose answer was "return f", and it would additionally need the conjugate of a `QuaCon`,
      whose return type is still open (Phase 6).
      * **f CONVEX -> f itself.** Decided exactly by `ratQ.isPSD2` on every piece.
      * **a CONCAVE or AFFINE piece on a bounded polygon -> the lower hull of the lifted
        vertices.** A concave function on a polytope has its whole envelope fixed by the VERTEX
        VALUES, so this is a hull of m points in R^3 with rational coordinates: a candidate facet
        is the plane through three of them (an integer cross product) and "is it a lower facet" is
        the sign of an integer for every other point. No hull library, no orientation predicate,
        no tolerance.
      Verified against the definition in three parts -- `co f <= f` with equality at the EXTREME
      points, convexity, and the largest-convex-minorant property -- over 38 random polygons, 0
      failures on all three (`.claude/envelope-sweep.m`).

- [ ] **`AlgAlg` is still not built, and its trigger is now precise: the INDEFINITE piece.** Both
      implemented envelope cases have RATIONAL faces (f itself, or affine functions interpolating
      rational vertex values), so nothing yet produces a face that must NAME a dual vertex. The
      indefinite envelope is where the first one appears -- an affine cell of a general `f**` is
      `<p,x> - f*(p)` with `p` of degree up to 4. Build `AlgAlg` when that case is built, not
      before, or the type is speculative.

- [ ] **What is left of the DOMAIN gap**: an unbounded piece whose `q` is NOT concave/affine (the
      finiteness test becomes a quadratic form on the recession CONE rather than a linear one --
      `PLQ:conjQ:unbounded`), a piece whose sup diverges everywhere so `dom f*` is EMPTY (correct
      answer, nowhere to put it, since a `QuaCon` carries at least one face --
      `PLQ:conjQ:emptyDomain`, the same representational gap `conjCPLQ` records as
      `conjugateHasEmptyDomain`), and a domain of dimension < 2. All three refuse by name.

- [ ] **A CONVENTION WORTH KNOWING, learned the hard way 2026-09-04.** Three of the five fixtures
      written for the unbounded sweep were MALFORMED, and each looked like a defect in `conjQ`
      first. `F(j,:)` is `[left, right]` of edge `j` with 0 meaning +inf, so the two rays bounding a
      wedge need OPPOSITE columns -- `[1 0; 1 0]` puts the face on the left of both, describes no
      convex wedge, and `QuaPol.eval` then returns `+Inf` inside its own domain (measured:
      `q(1,1) = Inf` on the first quadrant). A fourth signal was the oracle evaluating a multi-piece
      `f` with face 1's coefficients. Check a hand-built unbounded fixture with `f.eval` at an
      interior point BEFORE using it to judge anything.

      the per-piece closed forms first -- `convEnvCPLQ`, `conjPieceCPLQ`, `conjConvexOverPiece`,
      `conjConvexPolygon`, `conjAffinePLQ` are ALREADY 100% sym-free, so porting them is replacing
      double arithmetic with `ratQ` calls rather than rewriting an algorithm -- then Step 3
      (`maxQuaPar`), which is the large item, then the filtered predicates behind it.

- [ ] **`ratQ.hullQ` is deferred to B1** (the concave-quadratic-on-a-polygon envelope), where it
      can be validated against a real envelope rather than in isolation. It is the one Layer 0
      routine from the plan not built.

- [x] **DONE 2026-09-03: `QuaCon` (G16) exists.** Exact faces `fN`/`fD`, primitive integer edge
      conics `EcQ`, vertices as NAMES `[edgeA edgeB rootIdx]` realised on demand through
      `conicMeet`, faces in H-form as signed conic constraints. `conicMeet`, built 2026-08-24 and
      wired to nothing since, is now wired to this. `QuaConTest` 15/0.

## 2026-08-24 — SYM-FREE `conj`: what is left, measured. READ THIS FIRST.

The premise changed on 2026-08-24, and it changes the whole plan below: **vertices are stored as
INTERSECTIONS OF CONICS with floating-point coordinates, not exactly.** By
`CONJ_FIELD_PROOF.md` Theorem 1 the face functions and the edge conics are always rational, so the
only thing that ever needed an extension field was the vertex layer -- and naming a vertex by the
pair of conics it solves removes the need for one entirely.

**Two items below are CANCELLED by that premise, not deferred**: the degree-<=4 real algebraic
kernel, and the interim "factor the vertex quartic and refuse" gap. `exactQ` is likewise not on the
path to anything; it is referenced by nothing but its own test.

### What was measured, and it is the thing to know

The numeric conjugate path was **already 100% sym-free** before this work: `conjCPLQ`,
`conjPieceCPLQ`, `convEnvCPLQ`, `maxQuaPar`, `clipArcByConic` and `mergeSameQuadFaces` contain
**zero** `sym`/`subs`/`simplify`/`solve`/`isAlways`/`coeffs` between them, counting non-comment
lines. Every symbolic call on the conj route is behind ONE dispatch -- the fallback to Case C. So
the work is not "rewrite Step 3"; it is **shrink the set of inputs that fall back**, and
`checkConjSymFree` measures that set with the reason for each.

    2026-08-24 baseline, 16 fixtures:   SYMBOLIC 2 of 16, both maxQuaPar:notImplemented
    after conjConvexPolygon landed:     the unbounded CONVEX family joined the numeric route

### The gaps that remain, in the order they should be closed

- [x] **G1/G4/G10 -- LANDED 2026-08-28.** The parked `assemblePieces` diff (`collapseTinyEdges`
      + `matchHalfEdges`' sagitta test) is committed. Deliberate trade-off, user-authorized: it
      fixes the assembly-level defect below, at the cost of case 21 (a symbolic Case-C fallback
      input) changing from a fast named refusal to hitting the already-tracked, pre-existing
      Step 3 legacy gap (`SUPPORT_MATRIX.md` 1.2) instead of refusing early -- confirmed a
      pre-existing issue, not a new one introduced by this diff. Full suite verified clean
      (fast/normal/maxQuaParTest/conjCPLQTest/conjEdgeLowerBoundTest). See the commit and
      `DECISIONS.md` 2026-08-28 for the full trade-off record. Original entry follows, for the
      trail:
- [ ] **G1 -- a missing LENS in the overlay.** `maxQuaPar` can now SPLIT a cell whose arc is cut
      twice, or whose arc carries a neighbour's vertex in its interior (`bulgeSplit` / the
      passthrough split, 2026-08-24), and `maxQuaParTest` is green with it. The two fixtures that
      needed it then die two stages later in `assemblePieces`, and the piece dump says why:

          piece 1 src[1 1] carries g2's parabola from (0,0) to (1,1);
          its would-be neighbours are STRAIGHT edges through (0.5,0.5),
          which is not on that parabola.

      MEASURED with each operand's own point location: every point of the lens between g1's straight
      diagonal and g2's arc lies in **g1 face 4 and g2 face 2**, and the fold produces **no piece
      with src [4 2]** -- while the control point one step away, in [1 2], is there. So
      `clipByFace` returned nothing for a pair whose intersection is not empty, and the orphaned
      arc in `assemblePieces` is the symptom rather than the defect. `polyConstraints` already
      skips a curved edge's chord, so the candidates are the operand SWAP at the top of
      `clipByFace`, `clipPolyByConic`, and the three reduction passes. **Instrument that one pair.**
      `DECISIONS.md` 2026-08-24 (later) has the numbers; `SUPPORT_MATRIX.md` 4.1 records the same
      family of defect from 2026-08-13.
      Closing this removes the LAST measured fallback of the bounded family.

- [x] **G2 -- DONE 2026-08-25 (overnight). `conjAffinePLQ`, the all-affine route.** `max(0,x,y)`
      now conjugates NUMERICALLY, in well under a second, and never enters `maxQuaPar`.
      The construction is three lines of mathematics and is written out in the file: on face i,
      `f*(s) = max_i [sigma_{P_i}(s - a_i) - b_i]`, so
        * dom f* is the intersection over every face and every RECESSION direction of the
          half-planes `<d,s> <= <d,a_i>` -- which is why the answer has a BOUNDED domain although
          the input does not (`max(0,x,y)`'s three cones give s1>=0, s2>=0, s1+s2<=1: the simplex);
        * on it f* is the max of the affine functions `<s,w> - <a_i,w> - b_i`, one per (face,
          vertex), since a linear functional bounded on a polyhedron attains its sup at a vertex.
      The subdivision is then the ordinary one for a max of affine functions, by half-plane
      clipping. Wired into `conjPolygonalDomain` ahead of the fan route and guarded: it declines by
      name when dom f* is unbounded (the case Step 3 already does) and falls through, so it is a
      strict addition. Pinned by `conjAffinePLQTest` (3 tests, both fixtures against closed forms
      derived independently -- the simplex, and the L1 ball with a genuine 4-cell subdivision from
      `max(0,|x|-1,|y|-1)`) and end to end by
      `conjSymFreeTest/anAFFINEUnboundedFaceIsNowNUMERICAndCorrect`, which replaces the test that
      used to pin the gap. fast 306/0.

- [x] **G2b -- DONE 2026-08-24. `maxQuaPar` dropped a cell on some unbounded folds.** Found 2026-08-24 and it is the
      one silent wrong answer in the session: a 4-cone fan assembled to 4 cells and returned 2.0 at
      `s=(-2,-3)` where the definition sup is 4.5. Reproducer:
      `conjCPLQTest/step3UnboundedAssemblyMatchesTheTruth`'s fixture in its `F = [3 2;2 1;1 4;4 3]`
      orientation; `PARTITION OVERLAP` fires on it. A cross-check now catches this and falls back
      **Fixed the same night.** The cause was one line: at the SINGULAR POINT of a degenerate
      conic (a pair of lines crossing) the gradient vanishes, and
      `splitUnboundedAtOneCrossing` took its branch direction as the perpendicular to that
      gradient. It bailed, and the caller's tangency branch then read the winner at the cell's
      centroid -- which for a cone IS that same point, on the curve. The branches are now taken
      from the NULL DIRECTIONS of the quadratic part. `DECISIONS.md` 2026-08-24 (last) and
      `SUPPORT_MATRIX.md` section 4.5 have it; the two-line reproducer is
      `maxQuaParTest/twoHalfPlaneQuadraticsSplitTheirSharedQUADRANT`. The cross-check that caught
      it stays.

- [ ] **G4 -- `conj` of `xy` on some triangles is wrong, and it is the FOLD.** RE-MEASURED
      2026-08-25; the description below it replaced was wrong on both counts, so read this one.
      Step 1 splits sweep case 21's triangle into **four** faces (the nCE==3 cevian split with each
      half re-split), two of them slivers of area 8.7e-05 against 2.7e-02. **Every face's own Step 2
      conjugate is EXACT** at the bad point -- 1.032507658472 to twelve digits, four times over --
      and folding faces 2 and 3 keeps it; the FOURTH fold returns 1.005089907622. Pairwise folds of
      any two faces are all exact: only folding a sliver into the ALREADY-ACCUMULATED mesh loses it,
      which is why every two-operand `maxQuaParTest` passes.
      **It is also far worse than 2.742e-02**, which was only the worst on a probe grid of radius
      <= 6: the fold cross-check finds `f*(-10,0) = 47.10181578` against a true 10.86895777, an
      OVER-estimate by a factor of four. The G6 edge bound is one-sided and cannot see that.
      **Contained, not fixed.** `conjMaxOfSubTriangles` now cross-checks the fold against its own
      pieces (the identity `f* = max_k (q_k + I_T_k)*` makes them their own oracle) and refuses by
      name in 2.5 s as `PLQ:conjCPLQ:foldDroppedACell`, under `CCA2_CONJ_VERIFY` like the edge
      bound. It REFUSES rather than falling back because Case C did not finish in 25 minutes on it.
      **What is left is the `maxQuaPar` defect.** With `MAXQP_ASSERT = 2` this input raises two
      invariant violations, both real: a piece whose declared rays are the NEGATIVE of the direction
      its constraint region recedes along, and a piece carrying one operand's quadratic where the
      other is larger by `Inf` along a ray. The second is what produces the wrong value.
      `DECISIONS.md` 2026-08-25 (G4).

- [x] **G5 -- DONE 2026-08-25. `MATLAB:badsubscript`, and it was never 5-gon-specific.** The old
      entry said "SOME indefinite 5-gons" from sweep case 17; a 40-case re-run found case 29, a
      **4-gon**, crashing identically, so a fix aimed at 5-gons would have missed half of it.
      `plq_1p.conjugate`'s FRAME-CHANGE branch copied the z-frame object's ENVELOPE (2 faces) while
      replacing its `conjfia` -- the per-face block boundaries into `conjugates` -- with the single
      block `[1 nConj+1]`; `maximumConjugate` sizes its loop from the envelope, so face 2 asked for
      `conjfia(3)` of a 2-element array. Measured: `envelope faces = 2, numel(conjfia) = 2,
      nConjugates = 11`. The branch now carries `objT.conjfia`, and `maximumConjugate` takes its
      block count from `numel(conjfia)-1` -- the two arrays are not interchangeable, and the
      SEPARABLE branch legitimately returns one block over a multi-face envelope.
      **What the crash was hiding:** with the index repaired the cross-face symbolic max this input
      always needed actually runs, and `maximumP` on its two envelope faces did **not finish in 25
      minutes** -- measured both in the z-frame and after the read-back, with no difference. So case
      29 goes from a fast wrong crash to a correct computation nobody can wait for; the route that
      would make it fast is the NUMERIC one, which declines it today with `maxQuaPar:notImplemented`
      (clipPolyByConic separating an unbounded cell). Pinned by
      `conjCPLQTest/frameChangedPieceKeepsItsEnvelopeBLOCKS`, which asserts the INVARIANT rather
      than the value, because a value assertion here is a test nobody can run.
      `DECISIONS.md` 2026-08-25 (G5).

- [x] **G6 -- DONE 2026-08-25. The EDGE lower bound is a DEFAULT refusal.** Along each
      edge of `dom f` the objective `<s,x> - q(x)` is a quadratic in the segment parameter, so its
      maximum is closed form, and `f*` must be at least that. On the 24-case random sweep it fires
      on **exactly** case 21 (G4, at -2.742e-02) and on nothing else -- every other case sits at
      ~1e-15. The VERTEX-only version of the same idea catches NOTHING (0 of 24, including case
      21, because that sup is not attained at a vertex), so build the edge one.
      It is cheaper than the fold cross-check already in `conjPolygonalDomain` and, unlike it,
      covers the SINGLE-PIECE route. Built as `conjEdgeLowerBound.m`, raised by `conjCPLQ` as
      `PLQ:conjCPLQ:belowEdgeBound`, and ON by default (`CCA2_CONJ_VERIFY = 0` turns it off).
      Gated: fast 303/0, slow 88/0, verylong 26/7/1 -- the last IDENTICAL to a pristine `b9243d3`.
      It RAISES rather than falling back because the symbolic route returns the same wrong value on
      the known-bad case. `DECISIONS.md` 2026-08-25 (final).

- [x] **G7 -- DONE 2026-08-25. `plqStage`'s cache raced under `--verylong -j N`.** `save` now
      writes a unique temporary name in the SAME directory and `movefile`s it onto the real one, so
      a reader sees either the old file or the complete new one; and `load` is wrapped so an
      unreadable or half-written cache is treated exactly as a missing one -- it recomputes.
      NOT verified under contention: that needs the contended `--verylong` run this session was
      asked not to make. The change is one-directional (it can only turn a throw into a recompute),
      and the uncontended path is exercised by every staged test.

- [x] **G3 -- RECLASSIFIED 2026-08-25, and it was two things, neither of them what this said.**
      The old entry read "Declines by name today (`the fan-triangulation route needs a BOUNDED
      face`)" as though the case were unimplemented. That is the NUMERIC route's decline. Measured:

      1. **Envelope FINITE -> already answered, exactly.** `x*y` on the first quadrant is bounded
         below on its own recession cone, so conv f is finite and Case C returns the right function
         in ~14 s: f* is the INDICATOR of the third quadrant, 0 there and +inf elsewhere, checked
         against the closed form at ten dual points. So this half is a **symbolic fallback**, i.e.
         a cost, not a gap -- it belongs on the fallback list with G1 and G2, and what it needs is
         a numeric route.
      2. **Envelope -inf -> correctly refused.** `x^2 - y^2` on the same quadrant runs to -inf
         along the y-axis, which is IN the recession cone, so conv f = -inf and f* = +inf
         everywhere: dom f* is EMPTY. `convEnvUnbounded:minusInfinity` says exactly that. It is the
         RIGHT answer, blocked only by having nowhere to put it -- the same representation blocker
         as Case A's `conjugateHasEmptyDomain`, not a missing construction.

      Both are pinned by tests in `conjCPLQTest`. What remains under this heading is the numeric
      route for (1), which is a performance item.

- [ ] **G16 -- `QuaCon` storage (H-form). NOT STARTED, DELIBERATELY -- reviewed 2026-08-27.**
      (Renumbered from a SECOND G4 on 2026-08-26; two items carried that label and it caused
      confusion twice. The conjugate-of-`xy` G4 keeps the name.)
      Its own trigger is "build it when a `conj` result actually needs a non-parabolic edge -- i.e.
      when G1 lands and a three-piece input with two non-adjacent pieces reaches Step 3". **G1 is
      open**, so nothing can produce that edge today, and building the storage now would be
      speculative work with no way to test it against a real result.
      **Checked, rather than assumed, that nothing is silently broken meanwhile:**
      `QuaPar.assertParabolicEdges` runs in the constructor on every nonzero conic row, and its
      message already names this exact case -- "A conjugate CAN have an elliptical edge -- between
      two NON-adjacent pieces of f, see CONJ_FIELD_PROOF.md 7.4b -- and such a function is a
      QuaCon, not a QuaPar". So an elliptical edge RAISES `QuaPar:notParabola` by name; it cannot
      be stored as a wrong parabola. `Conic.m` and `RatCon.m` exist (the lattice can hold it);
      `QuaCon.m` does not, which is the item.
      **Do not start this before G1.** The first real elliptical edge is what tells you the storage
      is right, and until then there is nothing to validate against.
      **RE-CHECKED, 2026-08-30, now that G1 IS landed: still not triggerable, but for a DIFFERENT
      reason.** Ran `doc/QuaConExample.md`'s minimal 3-piece witness through `conj('cplq')`
      directly -- it dies on the exact same `PLQ:conjCPLQ:cplqFailed` at `s=(-121,-121)` this doc
      already recorded in 2026-08-20/21, unchanged by G1. That is `SUPPORT_MATRIX.md` 1.2's
      general Step-3 assembly gap, unrelated to `matchHalfEdges`/sagitta -- the run never reaches
      far enough to test `maxQuaPar` against a real elliptical edge. **The real trigger is 1.2
      closing, not G1** (G1 was necessary but the text overstated it as sufficient). Do not re-run
      this check again until 1.2 moves. `DECISIONS.md` 2026-08-30 (later) has the run.
      **REAL BUG FOUND AND FIXED tracing 1.2, 2026-08-30 (later still): `conjConvexOverPiece`'s
      vertex cone could collapse to a single point.** Traced the `cplqFailed` disagreement past
      assembly entirely -- piece 1 ALONE (no folding) returned `NaN` at `s=(-1,-1)`, one unit
      from the origin, where an independent brute-force oracle proves `f*=0`. Root cause:
      `edgeDirsAt`'s feasibility-probe step was the fixed constant `1e-6`, uncorrelated with
      `tolA` (scale-dependent); on `T1=(0,0),(60,10),(15,10)`, `Q=[3 -1;-1 5]` the resulting
      cross-facet violation (4.9e-7) came out smaller than `tolA` (1e-6), so both signs of one
      facet's tangent direction passed, collapsing the vertex-(0,0) cell from a 2D wedge to a
      point. FIXED: `step = max(1e-6, 1e4*tolA)`. Verified: 0/2000 mismatches on a random sweep
      against the oracle (was failing before); `unboundedFaceTest` 19/0 (new regression test
      `anisotropicQOnALargeTriangleCoversTheVertexCone` added); fast 312/0 and normal 12/0
      unaffected. `conjConvexOverPiece`'s first bug since it landed 2026-08-24 -- every existing
      test used small, roughly-unit-scale geometry that never hit the `step~tolA` knife-edge.
      **Does not by itself close 1.2 or unblock G16** -- the 3-piece witness needs re-running end
      to end to see how far it gets now; not done this session. `DECISIONS.md` 2026-08-30 (later
      still) has the full trace.
      **RE-RAN END TO END WITH THE FIX, 2026-08-30 (final): real progress, not closed.** No longer
      hits `cplqFailed` at `s=(-121,-121)` at all -- gets further, then dies on a NEW, different
      `MATLAB:badsubscript` (an array indexed one past its length). Not traced -- needs a stack
      trace first. 1.2/G16 still open, moved to a later pipeline stage. `DECISIONS.md` 2026-08-30
      (final) has the run.
      **SECOND BUG FOUND AND FIXED, same day: `symbolicFunction.tangent` crashed on a degenerate
      conic missing one ambient variable** (`obj.vars` from `symvar` can be length 1; `tangent`
      assumed 2). Fixed by an explicit `vars` argument, backward compatible. `DECISIONS.md`
      2026-08-30 (final, second bug).
      **CLOSED for this witness, 2026-08-30 (final, closure): the 3-piece elliptical-edge input
      now COMPLETES AND VERIFIES end to end.** `f*(s*)=2.9278190688` matches
      `doc/QuaConExample.md`'s own recorded value exactly; 0/200 mismatches against an
      independent oracle over a wide sweep. Pinned by
      `conjCPLQTest/threePieceEllipticalEdgeWitnessNowCompletes` (1/0, 467s, slow bucket).
      **G16's own precondition is genuinely met for the first time** -- a real elliptical edge
      reaches Step 3 and is handled correctly by `QuaParCPLQ`'s existing representation, so
      `QuaCon` remains an EFFICIENCY project (H-form), not a correctness gap. `SUPPORT_MATRIX.md`
      1.2 is a FAMILY of inputs, not fully closed by one witness, but this member of it -- the
      one this session traced -- is done. `DECISIONS.md` 2026-08-30 (final, closure) has the
      full chain: three real bugs found and fixed tracing one reported gap.

- [x] **G8 -- SUBSUMED 2026-08-25 by the legacy pins.** It recorded that `testPCE2`'s convex
      envelope is wrong one stage before the conjugate. True, and it is `plq_1piece`'s envelope:
      that fixture is one of the seven legacy reds now pinned by name (G14), and `conj` is exact on
      the same triangle. Nothing here is a `conj` defect. It comes back if the fixtures migrate to
      `plq_1p`.

- [x] **G9 -- FALSE as written, struck 2026-08-26.** It says six of the seven `verylong` reds are
      unnamed. G11 named all seven on 2026-08-25 (`--verylong -j 4`, 3699 s), with per-test costs
      and the class each one exercises, and G14 then pinned them. Read G11 instead.

- [ ] **G10 -- the `maxQuaPar` accumulated fold, which is what G4 actually is.** Half closed
      2026-08-25. `pieceRecessionRays` was deciding a piece's recession cone from floating-point
      noise -- its two sign tests compare against a mathematical ZERO (`A(d) = d*Q*d'` along a
      parabola's null direction, `cross(E,d)` along a ray edge) and were exact on doubles lifted
      with `sym()`. Measured +4.4e-20 vs -2.5e-17 for the SAME direction stored twice, and a
      half-strip reported receding along the negative of its own declared ray. Both now use a
      scaled tolerance. That gate is not cosmetic: `reccConeViolation` decides
      `halfIsProvedWellFormed`, so a false violation refuses a legitimate two-arc ray cut.
      **What is left is the winner on a straddling unbounded cell**, and the fix is NOT the
      obvious one -- routing `splitCell`'s single-tangency branch through `assignSide` fixes G4
      and breaks `conjSymFreeTest`'s A.5 triangle, which is correct to 3.6e-15 today. REFUTED,
      recorded in `DECISIONS.md` 2026-08-25 (G10); do not retry it.
      **Go upstream instead:** that branch is reached only when `splitCell` found ONE boundary
      crossing, and for an unbounded cell that is not a tangency -- the curve leaves through the
      RECESSION CONE, so the cell has one finite crossing and two parts.
      `splitUnboundedAtOneCrossing` declined, and HALF OF THAT IS NOW FIXED: its
      "must not run along one of the cell's own rays" test conflated PARALLEL with ALONG, so
      for a HALF-STRIP (dirIn == dirOut) -- where the escaping branch is parallel to both rays
      by construction -- every candidate was rejected. It now asks whether the crossing sits ON
      that ray edge, which is when a cut really would cut off nothing. The A.5 triangle splits
      and stays exact to 3.6e-15; `splitDeclined` records the reason under `MAXQP_ASSERT`.
      **SCOPE-REDUCED overnight 2026-08-25 to one named step.** The obstruction is
      `matchHalfEdges` rejecting a curved half-edge against a straight one that is the SAME
      boundary (both orphans had ZERO candidates -- rejected at generation, not consumed). The
      discriminator is the SAGITTA, measured to separate the cases by five orders of magnitude
      (true pair 1.172e-06, arc-vs-its-own-chord 1.127e-01), used RELATIVE to the endpoints' own
      residuals rather than against a fixed tolerance. Sub-tolerance edges must be collapsed first
      -- a piece whose own two vertices are 4.5e-06 apart against tolPos = 1e-3 cannot be matched.
      **The working attempt is `.claude/assembly_attempt_2026-08-25.diff`** and it gets case 21
      past assembly with case 8 still exact to 3.6e-15. `clipArcByHalfPlane`'s `internal` error on
      a degenerate arc (one endpoint outside so a crossing must exist, but the discriminant on a
      5.76e-05 arc came out negative) is FIXED already -- the bisection backstop landed
      2026-08-25 and is in the tree, so that is no longer what blocks re-applying the diff.
      **RE-APPLIED AND MEASURED, overnight 2026-08-27: not a clean win.** The FULL `conj('cplq')`
      call on case 21 (not just the assembly stage) goes from the documented 2.4 s
      `foldDroppedACell` refusal to 292.5 s and a DIFFERENT, unrelated failure -- an internal
      MuPAD error inside Case C's symbolic fallback. Reverted; the diff stays parked, unchanged.
      `DECISIONS.md` 2026-08-27 (overnight, G1/G10). **Whoever lands this next needs an
      acceptance test that checks the FAILURE MODE doesn't regress** (fast named refusal ->
      slow unrelated crash), not just "does assembly complete."
      **DONE 2026-08-28: the acceptance test is landed**,
      `conjCPLQTest.sweepCase21FailsFastAndNamedNotSlowAndUnrelated` -- pins
      `PLQ:conjCPLQ:foldDroppedACell` under 30 s on sweep case 21's exact fixture (reconstructed
      from seed 20260824). Any retry of the parked diff runs this test FIRST; a red here before
      touching anything else means the diff regressed the failure mode again, no re-measurement
      needed to find that out. The diff itself is still parked and unattempted this session.

      **PROVENANCE IS RULED OUT (2026-08-27), and so is tolerance tuning.** Measured on this
      fixture, 30 pieces: edge lengths min 4.484e-06, **median 2.987e-03**, max 9.982e-01, with
      **9 of 42 edges SHORTER than `matchHalfEdges`' tolPos = 1e-3** and 4 pieces whose whole
      diameter is under it. The orphaned piece 1 has its two vertices 3.0e-03 apart -- the piece is
      barely above the resolution at which assembly claims to identify points.
      Provenance answers "which candidate neighbour is right"; the difficulty here is prior to that
      -- the pieces' own coordinates, computed along different arithmetic paths, disagree by more
      than the features they must distinguish. Tagging an edge does not make two disagreeing copies
      of it agree.
      **The two directions that do address it:** (1) stop PRODUCING sub-tolerance pieces -- collapse
      them before assembly, which is the parked diff's first half and gets two stages further; or
      (2) make the coordinates agree, i.e. exact rational vertices for the polyhedral part
      (`ratQ`/`conicMeet`), which is large. **Start any future attempt from the edge-length
      distribution: if a change does not move it, it will not fix the assembly.**
      `DECISIONS.md` 2026-08-27 (item 2).

      **The hole is LOCATED, and G10 and G1 are the SAME defect (2026-08-25).** The uncovered
      point lies in g1 face 21 and g2 face 6 and the fold produced no piece with src [21 6] --
      G1's signature verbatim. But it is NOT `clipByFace`: tracked with `MAXQP_PROBE`, the cell is
      produced, survives every straight clip, the curved cut, and all four reduction passes, and
      is still covering the point when `assemblePieces` is called. **The loss is inside
      `assemblePieces`** -- `SUPPORT_MATRIX.md` 4.6's "live one", where `matchHalfEdges` pairs
      curved with curved and the subdivision is consistent per face PAIR rather than globally.
      That is the documented redesign, not a patch to `matchHalfEdges`.

- [x] **G17 -- DONE overnight 2026-08-27. `rectMaximumIsTheConjugateOfTheWholeDomain` holes at (-0.5,2). CAUSE UNKNOWN --
      `removeTangent` was accused twice and is EXONERATED.** Read the refutation before proposing
      anything: `DECISIONS.md` 2026-08-27 (G17, third pass). Deletion cannot cause a hole
      (it enlarges regions), and neither of the two regions `removeTangent` DISCARDS contains the
      point -- measured against all their constraints, misses of +9.10 and +3.98. The only evidence
      ever behind the identification was a stack trace of a symbolic WARNING, and `removeTangent` is
      where cPLQ's `isAlways` calls are noisiest, so it appears in traces from any run.
      **What is known:** 32 cells, (-0.5,2) uncovered, nearest cell misses by 1.668523e-02 -- four
      orders above the `maxQuaPar` assembly's scale, so a genuinely absent region, not a tolerance
      artefact.
      **BISECTED 2026-08-27: the loss is in FOLD 4's `maximumP`, after the pairing.**

          start:  3 cells, covered by 1
          fold 2: before=1  after mtimes=1  after maximumP=1  ( 8 cells)
          fold 3: before=1  after mtimes=1  after maximumP=1  (19 cells)
          fold 4: before=1  after mtimes=1  after maximumP=0  (25 cells)   <-- LOST
          fold 5: before=0  after mtimes=0  after maximumP=0  (32 cells)

      Covered entering fold 4 AND after `mtimes`; gone after `maximumP`. Every operand covers the
      point on its own (the five groups give 37.5, 37.5, 2.5, 2.5, 2.5; 37.5 is the truth), so no
      per-piece conjugate is at fault -- it is purely a fold defect.
      **ONE LEVEL DEEPER, 2026-08-27: it is `mergeL`'s FIRST pass, not the split loop.**

          entering fold 4:          19 cells, covered=1
          fold 4 after mtimes:      35 cells, covered=1
          fold 4 after SPLIT loop:  39 cells, covered=1
          fold 4 after mergeL #1:   27 cells, covered=0   <-- LOST
          fold 4 after mergeL #2:   25 cells, covered=0

      The split creates 4 cells and loses nothing (covered stays 1 through `mtimes` and the split
      loop), so `maximumP`'s own split/simplify path is ruled out.
      **Probe note that cost two attempts:** after `mtimes` each region carries TWO functions, so
      `evalFunctionNDomain` cannot read a paired object (it errors in `subs`) -- count coverage from
      the CONSTRAINTS at that stage.

      **ROOT CAUSE FOUND, 2026-08-27.** Cached the 39-cell pre-merge state
      (`fold4_objSplit.mat`) and replicated `mergeL`'s two accumulation loops by hand. The
      covering cell (12) merges with cell 17 (still covers), then merges with cell 21 --
      `region.merge` reports SUCCESS and the result no longer covers `s`. **Not a discard: a
      WRONG "the union is exact" certificate.** `unionIsExact` needs `r21`'s one curved
      constraint (`h = 40s2-4s1s2-s1^2-4s2^2-40 <= 0`) to hold on all of `r1217` (cells 12+17,
      by then purely polyhedral); `certifiesNonPositive` says yes, and it is wrong -- `h(s)=27.75
      > 0` and `s` is directly confirmed feasible for `r1217`. `h` is CONCAVE (Hessian eigenvalues
      `{-10,0}`, rank 1), so `certifiesNonPositive` maximises it via `quadprog` and treats exit
      code `ef==-2` as "the region is empty, vacuously true." Run by hand with `Display,'iter'`:
      `quadprog` genuinely fails to converge (`Fval` diverges past `-6.4e7`) on this unbounded
      cell with a rank-deficient Hessian, and `ef=-2` here means NUMERICAL FAILURE, not
      infeasibility -- the two share an exit code and the code cannot tell them apart. This is
      the actual defect, three levels deeper than the first bisection (`maximumP` -> `mergeL`'s
      first pass -> the 12+17+21 merge -> `certifiesNonPositive`'s `ef==-2` branch), and it is
      NOT `removeTangent` (exonerated twice) and NOT a discard anywhere in the split loop
      (exonerated this session). `DECISIONS.md` 2026-08-27 (G17, ROOT CAUSE FOUND).
      **Not fixed here -- deserves a clean-headed attempt, not a rushed one** (item 1's
      recession-cone fix was JUST refuted this same session by an analogous "trusted a
      plausible-looking answer without a second probe" mistake). **What a fix needs:** never
      trust `ef==-2` alone as "empty" -- verify emptiness independently first (e.g.
      `region.maxLinear` on the same `A,b`, which uses `linprog` and is not subject to this
      failure mode), or detect "unbounded objective, degenerate Hessian" directly (Q's null space
      intersected with the region's recession cone, checking `L'd` there) and refuse rather than
      guess. Reproduce this EXACT case (`fold4_objSplit.mat` cells 12, 17, 21 -- scratchpad,
      session-local, not committed; regenerate via `g17f_save.m`'s method if needed) as the
      regression test before touching `certifiesNonPositive`.

      (Superseded text follows, kept for the reasoning.) ~~Identified 2026-08-27.~~ The nearest cell to the
      uncovered point (-0.5,2) misses by **1.668523e-02** -- a real gap, four orders above the
      1e-6..1e-3 scale at which the `maxQuaPar` assembly fails, so it is NOT that defect and no
      tolerance or snapping scheme touches it. The failing stack is
      `region.removeTangent <- functionNDomain.mergeL <- maximumP <- maxOfList <-
      plq.maximumConjugate`, which is the gap `DESIGN.md` records in five places as "the remaining
      known bug".
      **MECHANISM CORRECTED the same day -- it is NOT deletion.** Deleting a constraint makes a
      region BIGGER, and no enlargement can create a hole, so "deletes a constraint the region
      needs" was wrong. With every branch instrumented: the three `deletedLinear` events are all far
      from the hole, and what matters is `removeTangent`'s `lin & ~tin` branch, which does
      `obj = region.empty` -- **it discards the WHOLE region**. Two regions were discarded, and the
      second one's quadratic is satisfied at the hole
      (`4(-0.5)(2) - 40(2) + 0.25 + 16 + 40 = -27.75 <= 0`).
      **Next step, cheap:** capture the full discarded region at that branch and test (-0.5,2)
      against ALL its constraints -- containment is currently established only for the one named
      constraint, so this is a strong candidate rather than a proof. Then decide whether
      `lin & ~tin` should discard at all, or refuse to conclude when the probe sits at a tie point.
      **Work it as itself, not as part of the assembly.** `DECISIONS.md` 2026-08-27 (G17, corrected).

      **ROOT CAUSE FOUND AND FIXED, overnight 2026-08-27.** Traced three levels deeper than the
      corrected version above: the loss is a real merge (cells 12+17, then +21) that
      `region.merge`/`unionIsExact` wrongly certify as exact. `certifiesNonPositive`'s concave
      branch trusted `quadprog`'s `ef==-2` as "region empty, vacuously true" without verifying
      it -- on this unbounded cell with a rank-deficient Hessian, `quadprog` genuinely fails to
      converge and returns `ef==-2` as a non-convergence code, not an infeasibility one (confirmed
      directly: the point is feasible, `quadprog`'s own iterates diverge). Fixed by verifying via
      `linprog` (`region.maxLinear`) before accepting `ef==-2`. Fast 309/0, normal 12/0, no
      regressions; the isolated reproducer (cells 12/17/21) now correctly refuses the bad merge.
      **The end-to-end re-run of the original fixture was ABANDONED after 12 hours of continuous
      CPU-bound execution** (killed, not hung -- confirmed via CPU time tracking the whole way).
      Not a correctness concern (every other check is green); the honest fix's real verification
      work is apparently far more expensive on this fixture than the old wrong-but-fast merge was.
      **Left `UNVERIFIED` at the full-pipeline level, correctly reported as such.** The isolated
      three-cell reproducer stands as the verification this fix has. `DECISIONS.md` 2026-08-27
      (G17, ROOT CAUSE FOUND) and 2026-08-27 (overnight, G17 end-to-end).

- [ ] **G11 -- the seven `verylong` reds, NAMED at last (2026-08-25, `--verylong -j 4`, 3699 s).**
      45 jobs, 35 pass, 9 fail, 1 timeout -- and the 9 are the 7 pre-existing ones, with `testPCE2`
      now reporting as three because the split separates its stages. No regression.

          testcPLQ/rectConjugateMatchesTheSup                     GREEN 2026-08-25 (was FAIL)
          testcPLQ/rectMaximumIsTheConjugateOfTheWholeDomain      FAIL   1562 s   plq_1p
          testcPLQ/rectBiconjugateIsAConvexUnderestimator         TIMEOUT >3600 s plq_1p
          testMaxMultiRegion/testMaxBiconjugate                   FAIL   1360 s   plq_1piece
          testMaxMultiRegion/testMax3                             FAIL    557 s   plq_1piece
          testMaxMultiRegion/testPSqroot                          FAIL     75 s   plq_1piece  G13
          testMaxMultiRegion/testOpenconvex                       FAIL    113 s   plq_1piece  G14
          testMaxMultiRegion/pce2{Envelope,Conjugate,Step3}       FAIL  57/83/76 s plq_1piece

      SECOND GATE 2026-08-25, after the plqCheck oracle fix: **36 pass / 8 fail / 1 timeout** in
      3695 s, against 35/9/1 before. One red fixed, NO new reds -- so the sharper lower bound
      exposed no hidden minorant anywhere in the two suites.

      **SCIP-RELEVANCE, checked 2026-08-28 (SCIP_READINESS.md Phase C).** These two are not an
      orthogonal legacy-pipeline concern for SCIP: `QuaParCPLQ.conj` -- reached by `biconj('cplq')`
      whenever a QPLIB term's domain is non-box (Case C) -- reuses this exact
      `functionNDomain.maxOfList`/`mergeL` fold verbatim (`QuaParCPLQ.m:59`). Re-ran
      `rectBiconjugateIsAConvexUnderestimator` this session (killed after ~15 min, not to
      completion -- see AI/CLAUDE.md §9 on verylong budget): same warning signature
      (`isAlways:TruthUnknown` in `region/getEdgeNos`, `functionNDomain/mergeL`) as
      `.claude/step3cost.m`'s profiled fold, not a different failure. So optimizing THAT fold
      (SCIP_READINESS.md Phase C's `getVertices` closed-form lever, not yet done) is the same work
      as fixing this test's timeout, not a separate task.
      **Read G14 before taking any of these.** Every red measured so far belongs to `plq_1piece`,
      and `conj` is exact on the same inputs -- so they are not `conj` bugs. The two that still
      matter for `conj` are the two `plq_1p` ones: `rectMaximumIsTheConjugateOfTheWholeDomain`
      (1562 s, needs the stage split the PCE family got) and
      `rectBiconjugateIsAConvexUnderestimator`, which does not finish at all and where the first
      job is to find out where the time goes.
      **G7 is verified by this run:** `-j 4` shards per test, consecutive stages of a fixture ran
      concurrently, and no cache error or spurious red appeared.

- [x] **B5 -- STALE, corrected 2026-08-25.** `SUPPORT_MATRIX.md` listed "unbounded multi-face
      domain, a face with a CURVED CONVEX envelope" as a GAP against
      `plq_1p:conjugateFunction:unboundedNonAffine`. That case is DONE and has a green test:
      `conjConvexOverPiece`'s KKT active-set decomposition covers the bounded triangle and the
      unbounded wedge/half-strip with the same code -- a ray contributes its direction exactly as a
      bounded edge does -- and `unboundedFaceTest/convexFaceOverAWedgeGoesThroughTheWholePipeline`
      records that it *used to* raise that identifier.
      What the guard actually tests is `envUnbounded && ~envIsAffine`, reached only AFTER the
      convex branch has returned -- i.e. an unbounded face whose envelope is RATIONAL (the nCE==1 /
      A.4 form). **No input has been found that reaches it**, and there is a reason to expect none:
      over an unbounded face the recession directions force an indefinite quadratic's envelope to
      be affine (`x*y` on the first quadrant, and on the half-strip `0<=x<=1, y>=0`, both give
      conv f = 0) or `-inf` (`x*y` on `y>=0`), and both are handled. Mark it N/R when someone
      proves it; until then it is a guard over an unproven-unreachable branch, not a feature.

- [x] **B6 -- NOT A GAP, corrected 2026-08-25.** Listed as "Step 3 with a non-triangular envelope
      piece". The guard is in `conjPolygonalDomain`'s fan loop and says something else: it fires
      when `conjSingleTriangle` hands back a `QuaParCPLQ`, i.e. that ONE TRIANGLE's Step 2 fell back
      to symbolic, and Step 3's numeric max cannot mix a mesh with cPLQ's symbolic form. Refusing to
      mix them is right. And it is not a gap: `conjCPLQ` lists `PLQ:conjCPLQ:notImplemented` in its
      fallback catch, so under `auto` the domain routes to Case C and gets an answer. A symbolic
      FALLBACK, and one that shrinks automatically as the Step 2 fallbacks (G1, G4) close --
      which is why it was scheduled last.

- [x] **G12 -- DONE 2026-08-25 (overnight). `route='symbolic'` now has a destination for a single
      triangle.** Case B ignored the route entirely, so a triangle came back as the same numeric
      mesh whatever the caller asked for -- and `biconjCPLQ` asks for the symbolic form exactly
      when the numeric first conjugate is a CURVED mesh, because the second conjugation cannot take
      one. The escape did nothing. The destination is cPLQ's own PER-TRIANGLE form,
      `conjEnvelopeViaCPLQ`, which `conjSingleTriangle` already falls back to -- NOT Case C, which
      does not cover a single triangle and raises `cplqFailed` after 102 s (that route was tried
      and is refuted in `DECISIONS.md`).
      Measured: `x*y` over {(0,0),(1,0),(0,1)} returns a QuaParCPLQ in 14.4 s, exact to 0.000e+00
      against the closed form `max(0,s1,s2)`. On the {(-1,1),(-3,-3),(-4,-3)} triangle cPLQ itself
      raises `cplqFailed` -- a pre-existing limit of that pipeline on that fixture, now named
      rather than silently answered with the numeric mesh. Pinned by
      `conjSymFreeTest/routeSYMBOLICHasSomewhereToGoForASingleTriangle`.

- [x] **G13 -- SUBSUMED 2026-08-25 by the legacy pins.** Its measurement stands and is worth
      keeping (`conj` returns 1.75 where the legacy route returns -5, the objective at the vertex
      (-4,-3)); the ACTION is G14's, which pinned `testPSqroot` with exactly that evidence in the
      message. Nothing separate left to do.

- [x] **G14 -- DONE 2026-08-25: the seven legacy reds are PINNED, not counted. Measured
      2026-08-25, and this is the finding that should decide what happens to them.** Each red was
      taken to its fixture and the SAME input put through `QuaPol.conj`:

          red                                  class used     conj on the same input
          ---------------------------------    ------------   -----------------------------------
          pce2* (was testPCE2)                 plq_1piece     EXACT at all 9 dual points
          testPSqroot                          plq_1piece     EXACT (1.75; legacy gives -5)
          testOpenconvex                       plq_1piece     plq_1p REFUSES, correctly
          rectConjugateMatchesTheSup           plq_1p         was a FALSE red -- oracle, now green

      * `testPSqroot`: legacy takes the VERTEX value at (-4,-3) where the sup is strictly inside
        the edge to (-1,1). `conj` gets 1.75.
      * `testOpenconvex`: legacy's envelope contains **2147483647** -- `intmax('int32')`, the ray
        DIRECTION MARKER read as an ordinary coordinate -- and exceeds f by 2.147e+09. `plq_1p`
        raises `convEnvUnbounded:minusInfinity`, which is right: `x*y` on that half-strip runs to
        -inf along x = -1, so conv f = -inf. This is the exact defect `quaPolToPlq`'s header says
        was fixed *in plq_1p*; `plq_1piece` never got it.

      **So `conj` is not implicated in any red measured so far.** The unmeasured ones are
      `testMax3` and `testMaxBiconjugate` (both `plq_1piece`, so expect the same) and testcPLQ's
      `rectMaximumIsTheConjugateOfTheWholeDomain` / `rectBiconjugateIsAConvexUnderestimator`, which
      DO use `plq_1p` and are therefore the two worth measuring next.
      **DECIDED: pinned as legacy.** All seven now call `testMaxMultiRegion.legacyPin`, which
      `assumeFail`s with the measurement in the message, so they report INCOMPLETE rather than
      FAILED and cost 1-3 s instead of 75-1360 s apiece. The bodies are untouched: whoever migrates
      these fixtures to `plq_1p` (T6) deletes one line per test to get the checks back.
      `testMax3` and `testMaxBiconjugate` are pinned on the CLASS rather than a measurement -- their
      pins say so, and confirming them is part of the migration.

- [x] **B3 -- DONE 2026-08-27. All three cases are ANSWERED; none was blocked on the return type.**
        * **POINT** (Q = 0) -- 2026-08-25. A needle; `QuaPar.eval` had no branch for it.
        * **LINE** (PSD rank 1) -- 2026-08-26. Two opposite rays; `eval` ignored the ray marker and
          `belongToEdge`'s ray test was coordinate-wise.
        * **EMPTY** (a negative eigenvalue) -- 2026-08-27. f* = +inf everywhere, which a
          full-domain mesh with constant `+inf` represents exactly. Verified safe downstream
          (eval/isMeshed/kind/isConvex/addition) before adopting.
      **One consequence, pinned:** conjugating the EMPTY result gives -inf everywhere -- correct
      ((+inf)* = -inf) and outside the PLQ class. `biconjCPLQ` never reaches it today.
      **The transferable lesson:** three for three, "blocked on the return type" was wrong and the
      representation already existed. Build the object by hand and evaluate it before recording
      that phrase again. `DECISIONS.md` 2026-08-27 (B3).

- [x] **G15 -- ANSWERED 2026-08-26, and the hypothesis was WRONG.** The
      `rectBiconjugateIsAConvexUnderestimator` timeout is NOT item 2's hole. Measured:
      `conjugateOfPiecePoly` -- the step that consumes the holed mesh, and the only one a hole
      could plausibly derail -- takes **191 s** of the 3600 s budget and TERMINATES, turning 32
      cells into **111 cells in 32 blocks**. The remaining cost is the pointwise max across those
      111, i.e. `functionNDomain.maxOfList`.
      That is the defect already recorded below as "STEP 3's CROSS-PIECE MAXIMUM DOES NOT SCALE",
      whose own measurements run 5, 14, 29, 45, 70, 86 cells against folds of 93 ... 2087 s. At 111
      cells a 3600 s timeout is the expected cost, not a hang. **G15 folds into that entry**, and
      its fix -- merge same-function neighbours after each fold, polyhedral cells first -- is the
      one to try. `DECISIONS.md` 2026-08-26 (G15).

- [x] **Code coverage, industry-standard measurement. DONE overnight 2026-08-27.** MATLAB, not
      Python, so `pytest-cov` does not apply; wired `matlab.unittest.plugins.CodeCoveragePlugin`
      (through `TestRunner`, not the plain `runtests` convenience function) into
      `.claude/suite.sh --fast --coverage` / `--normal --coverage`, emitting a Cobertura-XML
      report to `.claude/coverage/cobertura.xml` (gitignored, like `plqStage`'s cache). Scoped to
      the fast bucket by default (cheap, re-runnable); `--normal --coverage` also works.
      `--slow`/`--verylong --coverage` is rejected outright rather than silently doing something
      partial -- sharding coverage per-suite would need merging multiple Cobertura reports, not
      attempted. Verified: `--fast --coverage` matches plain `--fast`'s 309/0/0 with negligible
      overhead (147s vs ~150s) and produced a real report -- **47.4% line coverage from the fast
      bucket alone** (9825/20719 lines). Two MATLAB/bash interop snags recorded in the commit
      message for anyone touching this again: the report path must come from MATLAB's own `pwd`
      (bash's `$CCA2DIR` is a Git-Bash path native MATLAB silently fails to write to), and
      `TestSuite.fromClass` needs array concatenation for multiple classes, not a comma list.

### Tools built on 2026-08-24, so they are not rebuilt

- `checkConjSymFree.m` -- the fallback RATE and its reasons, per fixture. Run it before and after
  any change to the dispatch.
- `conjCPLQ(q, [], 'numeric')` -- refuses to fall back, so a test can pin the ROUTE. `conjSymFreeTest`
  does exactly that, including for the two gaps above (a gap test going GREEN is good news: promote
  it and delete the entry here).
- `conjConvexPolygon.m` -- a convex quadratic over ANY convex polygon, bounded or not, closed form,
  no triangulation, returns a QuaPol.
- `checkConjAgainstDefinition.m` -- randomized `conj`-vs-definition sweep over random polygons and
  all four sign classes. 22 of 24 exact; the two that are not are G4 and G5, both pre-existing.
  Slow by design (the reference is a scan plus a pattern search per probe), so it is a CHECK, not a
  bucket member. Run it after any change to the conjugate, and against a snapshot of the old tree
  when something fails -- that is what established G4/G5 were not introduced on 2026-08-24.
- `conjPolygonalDomain`'s fold cross-check -- the assembled Step 3 result is verified against
  `max_k (q_k + I_P_k)*`, the identity it was built from, and DECLINES on a disagreement so the
  caller falls back. One-sided by construction: it can miss a defect, it cannot invent one.
- `ratQ.m` / `conicMeet.m` -- exact rational coefficient arithmetic, and conic-conic intersection
  through an exact integer quartic.

---

## 2026-08-20 — THE SYM-FREE PORT, RE-PLANNED (superseded above, kept for its measurements)

The plan changed on 2026-08-20 because two things were measured, both in `DECISIONS.md` under that
date: one quadratic extension is not enough (a single A.5 triangle needs `sqrt(15)` and `sqrt(30)`
in ONE coordinate), and no tower of square roots is enough either (an `f*` vertex can have degree 4
with Galois group S4). Order matters below; each item says what decides it.

Terminology used throughout: **Row 7** is the recommendation of `CONJ_FIELD_PROOF.md` section 8.2 —
store the mesh in **H-form** (faces as sign conditions on rational conics, edges as rational conics,
a vertex named by the PAIR OF CONICS it solves rather than by coordinates) and run every predicate
through an interval **filter** first, dropping to the exact kernel only when the interval straddles
zero. That file is UNTRACKED at the time of writing — another session left it in the tree — so this
section states what it needs rather than pointing at it.

- [ ] **T3 — `symbolicFunction`'s payload becomes a coefficient vector, over RATIONALS.**
      The face and conic layers need no extension field at all (that document's Theorem 1: face
      functions and edge conics of `f*` are always rational), so `Rat`/int64 covers them and
      `exactQ` is not load-bearing there. This is where the engine calls actually are — live
      counts, non-comment: `subs` 72 in `region.m`, 68 in `plq_1piece.m`, 21 in `plq_1p.m`;
      `coeffs` 18 + 10; `simplify`/`simplifyFraction` 11 + 18 + 10; `hessian`/`gradient`/`dfdx`
      about 45. Every one becomes arithmetic on coefficients.

- [x] **CANCELLED 2026-08-24 by the H-form premise — a degree-<= 4 real algebraic kernel**: a rational quartic plus an isolating
      interval, sign of a rational polynomial at it, and comparison by resultant or Sturm sequence.
      `exactQ` is now multiquadratic and still NOT enough: of twelve continuous three-piece
      configurations the vertex quartic is irreducible over Q in ten of them, and the S4 case is
      proved reachable. Bounded work, because the degree is capped at 4 for both `conj` and the
      envelope.

- [x] **CANCELLED 2026-08-24 with the kernel above — DETECTED refusal.** Factor the vertex quartic and refuse
      by name when it does not split into rational or quadratic factors. That turns a reachable
      wrong answer into one nameable `SUPPORT_MATRIX.md` GAP, which is the discipline this project
      already has for unreachable branches.

- [ ] **Row 7 itself — H-form storage plus filtered predicates.** Biggest rewrite: V-form is baked
      into `RatPar`'s `V/E/F/P`, `eval`, `createP`, `orderEdges`, plotting, and every test that
      names a vertex. Plotting and user output still need coordinates, but then as an OUTPUT
      convenience, not the stored truth. The filter must be certified (a real error bound, not a
      tuned tolerance) with a terminating exact fallback, or it is the refuted double-plus-tolerance
      design wearing a disguise.

- [ ] **When Row 7 lands, INDEX the mesh for evaluation — do not linear-scan.**
      Measured (`.claude/evalbench.c`, gcc -O2 -march=native, min of three runs, value AND gradient
      every call; the caveats are in `DECISIONS.md`):

          baseline: 20-node expression tree, forward-mode gradient      44 ns
          linear scan          9 / 81 / 1024 cells             38 / 130 / 1670 ns
          uniform bucket       9 ... 1024 cells                24 ...    27 ns
          slab, two binary searches                            15 ...    26 ns
          cached cell + bucket, COHERENT queries                        11.5 ns
          no mesh, per-piece closed-form sup   3 / 6 / 24 pieces  105 / 215 / 822 ns

      Take the **uniform bucket plus a last-cell cache**: flat in the cell count, O(N) to build,
      and SCIP's queries move in small steps so the previous cell usually still holds. Skip the
      slab method — its 15-26 ns is measured on a grid of axis-aligned boxes, which flatters it,
      and on a real conic arrangement each of its ~log2(N) steps costs a conic evaluation. Keep
      unbounded cells in a separate short list, tested only when the point is outside the global
      bounding box. The index is pure preprocessing, microseconds to build.

- [ ] **T6 — delete `plq_1piece`** (57 live engine calls, 13 of them `solve`). RE-RUN THE FIXTURE
      SWAP FIRST: of the three regressions recorded on 2026-08-19, one may already be gone, since
      `testPCE2`'s domain is the triangle A.4 now gets right. The other two are specific —
      `testFractional` needs `conjugateFunction` to dispatch on the envelope's KIND (a rational
      face to the A.3 branch) and not only on its domain, and `testConvex` feeds a four-vertex
      polygon to a routine that wants triangles.

- [ ] **T4 / T5 — the remaining global sign questions and `solve()` sites.** `signEverywhere`,
      `impliedBy`/`holdsOn`, `certifiesNonPositive` are the only `isAlways` uses that ask about a
      FUNCTION rather than a number; degree <= 2 has a closed-form PSD/minimum test and half of it
      is already written. Live `solve`: 6 in `region.m` (most measured dead), 3 in
      `functionNDomain.m`, 5 in `symbolicFunction.m`.

**Facts about the SCIP side worth keeping, so they are not re-derived:** `f*` is convex, so the
gradient at any point is already a global linear underestimator, the max over a box is attained at
one of four corners, and the min is the global min when it is inside — the interval evaluation that
usually dominates propagation is nearly free. Evaluation is pure double precision: the exact kernel
is a BUILD-time cost and the finished mesh compiles to flat arrays, so the degree-4 vertex question
never reaches the solver. The gradient jumps across cell boundaries; the value does not, and either
side is a valid subgradient, so separation stays correct but an Ipopt-based heuristic will feel it.
Building the mesh is minutes to an hour, so it must be once-per-constraint preprocessing, never
reached from a node.

## 2026-08-23 — RETURN TYPES: `QuaCon` for `conj`, `AlgCon` for `biconj`

Why the current types cannot hold the answers is `DECISIONS.md` 2026-08-23 (envelope face type) and
2026-08-21 (`f*`'s elliptical edge). Both axes of the `RatPar` lattice grow by one level ABOVE the
present top, so nothing existing changes behaviour:

      subdivision:   Pol  <  Par  <  Con        `Con` drops b^2 - 4ac = 0
      function:      Qua  <  Rat  <  Alg        `Alg` = root of a rational quartic in z

- [x] **DONE 2026-08-24 (as `Conic`; `CON` is a Windows device name)** — a new trait plus a new
      data-holding parent `RatCon`. Cheapest
      item here, and the elliptical edge already forces it. `QuaPar` becomes a real specialization
      instead of a type that cannot hold the values `conj` produces.

- [ ] **`QuaCon` = the return of `conj`.** Rational faces `f` and edge conics `Ec` as int64
      `[a b c d e f | den]`, primitive and sign-normalized so equality is bitwise (the
      `4 - 2 sqrt 2`-as-two-doubles failure is what canonical integers prevent). Vertices become
      NAMES: `Vname(nv,3) = [edgeA edgeB rootIdx]`, the point where two edge conics meet, `rootIdx`
      canonical among the real intersections. `Vx` (double) and `Vbox` (rational isolating box) are
      CACHES — rebuildable, deleting them changes runtime only. This is Row 7; it is the expensive
      item, because `V`, `eval`, `createP`, `orderEdges`, plotting and every vertex-naming test read
      coordinates today. Do it on `QuaCon` alone; the four legacy types keep coordinate `V`.

- [ ] **`AlgCon` = the return of `biconj`, stored as a DECORATION of the `QuaCon`, not standalone.**
      Forced, not stylistic: an affine cell is `<p,x> - f*(p)` with `p` a dual vertex of degree <= 4,
      so its coefficients are irrational unless the cell NAMES the dual vertex. Payload:
      `src` (the `QuaCon`), `fkind` in {QUAD, RAT, RULED, AFFINE}, `fref` (dual cell pointer plus
      which adjacent face plays `i`), `fdeg` in {1,2,4}, `froot`; and `Ekind` in {CONIC, RULING}
      because the extreme rulings join two degree-<=4 points and their supporting line is rational
      only over `Q(p)`. `f`/`den` stay as caches for the QUAD and RAT faces.

- [ ] **`RatPol.m`'s header is now wrong** — "quadratic over linear on a polyhedral subdivision,
      proven not to need a square root" holds for ONE piece. Fix the comment when `AlgCon` lands;
      `RatPol` becomes `fkind = RAT, fdeg = 1`.

Evaluation stays pure double: locate the cell (uniform bucket + last-cell cache, the measured
choice), then dispatch on `fdeg` — coefficients, one square root, or one quartic. `grad f** = s` is
free and the Hessian is the rank-one closed form. The degree-<=4 exact kernel is a BUILD-time service
behind the interval filter and never reaches `eval` or SCIP.

---

## 2026-08-13 -- the far-field defect, worked steps 1-6 of the plan

**Where it stands (updated 2026-08-13, end of the second pass).** `maxQuaParTest` **25 pass / 1 fail**; fast bucket 200 pass / 1 fail. The one red is `arcVsArcRefusesAnUnboundedTwoArcSplit` (first item below). `sweepMaxQuaParCurvedSplit(20260802,200)` went 30 -> 59 assembled of 142 sampled, with 0 of 1031 result vertices, 571 midpoints and 3540 interior points wrong. What the whole defect turned out to be is one sentence, and it is recorded in `SUPPORT_MATRIX.md` 4.1: **a curved edge is a bounded ARC and its conic is not**, so the point-location rule admits far-away points on a parabola's concave side; `QuaPar.chordCuts` derives the missing constraint. Original note follows.

**Where it stood mid-session.** `maxQuaParTest` 21 pass / 5 fail. The seeded far-field sweep
(`arcVsArcMatchesGroundTruthOverRandomShifts`) and the unit-square pin are GREEN, and on a
397-quadrilateral seeded sweep the arc-vs-arc results wrong in the far field went **7 of 64 to 0
of 64** (200 directions at radii 1, 5, 50, 500).

**Three exact piece invariants** (`assertPiecesWellFormed`, global `MAXQP_ASSERT`, off by default;
1 = containment + winner domination, 2 = also the symbolic recession cone): containment in the
source cell (vertices, ray directions, straight edges against the source arc, and the piece's own
arc against every source constraint), the carried operand actually dominating on the piece, and the
encoded region receding exactly where the piece declares rays. Nothing sampled.

**Two producers fixed.**
1. `splitCell` read "one finite crossing" as a tangency and returned the cell intact. On an
   UNBOUNDED cell the curve can enter there and leave through the recession cone, cutting it in
   two; `splitUnboundedAtOneCrossing` decides the two cases from whether the branch direction
   recedes the cell, and each half keeps one original ray plus the escaping branch.
2. `clipPolyByConic` skipped a cut whose only crossings are AT the cell's corners -- the generic
   arc-vs-arc case, since conjugates of triangles sharing a primal edge have arcs through the same
   two dual points. Now decided from one representative point per boundary element; a two-arc
   survivor is cut along `A -> M -> B` with M between the arcs and each half given the midpoint as
   a third vertex (two-vertex lunes sharing a chord defeat assembly: the chord never becomes an
   edge).

**Two new tools.** `QuaPar.eval` validate mode (`QUAPAR_VALIDATE`) errors when two faces admit a
point with DIFFERENT values; `verifyMaxIsExactSymbolically` proves `g = max(g1,g2)` over whole
regions by closed-form minimisation over each face-pair intersection. The three arc-vs-arc fixtures
verify with zero findings.

### Fixed since (2026-08-13, second pass)

- The **unit square now VERIFIES**: zero findings from `verifyMaxIsExactSymbolically`, zero
  ambiguous points from eval's validate mode over 900 ring points. Two causes, both "an edge
  identified by its endpoints alone", which an arc-vs-arc arrangement makes ambiguous because FOUR
  arcs run between the same two dual points: `matchHalfEdges` now also requires two curved
  half-edges to lie on the same conic, and `QuaPar.orderEdges` now looks the next boundary edge up
  among the FACE's own edges instead of the whole edge list.
- `fixArcTag`: a clipped half that does not keep the arc no longer carries the arc's CONSTRAINT.
- `splitTwoArcPiece` re-locates both curve indices from geometry, so a stale index from a rebuilt
  vertex list can no longer index off the end (the two seeded crash fixtures no longer crash).
- `dropDegeneratePieces`: collapsed pieces (2 vertices, no arc, bounded) no longer reach assembly.

### 2026-08-15 -- FOUR of the five bugs fixed; ONE red left in the whole repository

`maxQuaParTest` 29 / 0, `conjCPLQTest` 25 / 0, `unboundedFaceTest` 18 / 0, fast bucket 204 / 0,
slow bucket **114 / 1**. Bugs 2, 3, 4 and 5 are fixed; bug 1 is at 5 of 7 probe points, was 0 of 7,
and is the only remaining red. `maxQuaPar` has NO open case -- the seeded arc-vs-arc sweep is
**18 exact / 0 wrong / 0 errored of 18**, from 16 / 0 / 2 on 2026-08-13. What it took is in `SUPPORT_MATRIX.md` 4.1 and `DECISIONS.md`; the short version is that
**neither of the two earlier attempts at the last red failed for the reason it was thought to** --
the tooling that judged them was itself broken, in two ways, and silently.

### Still open, in the order they should be taken

- [x] **`MAXQP_ASSERT` is now ON in `maxQuaParTest`** (2026-08-15), at level 1, via a
      `TestMethodSetup` that restores the previous value on teardown. 28 / 0 with it on. Level 2
      stays opt-in -- it costs seconds per call, and the tools that want it call
      `pieceRecessionRays` directly.

### Next up (2026-08-18) — where Step 3 actually stands

**MEASURED end to end on the A.4/A.5 quadrilateral, after this session's work:**

    cells per fold   5, 14, 29, 45, 70, 86   ->   5, 12, 23, 38, 51, 60
    total            73 min                  ->   43 min (2579 s)

The machine was running three MATLABs throughout, so that timing is pessimistic, and a single
timing settles nothing anyway (see `CLAUDE.md` §3).

**What got it there, all measured:** three double leaks fixed so Step 2 is exact
(`domain.mE`/`cE`, `region.limitOfFAtVertices`, `plq_1p.quadPartsOf` + `conjConvexOverPiece`); a
sound certificate for a curved constraint (`region.certifiesNonPositive`) replacing merge's two
quadratic heuristics; `quadprog` deciding the CONCAVE conics, which is what the conics here
actually are; and `functionNDomain.singularEdgeCut` closing the singular-quadratic overlap.

**AND `testcPLQ/testRectBiconj` NOW PASSES with `CCA2_A45_SPLIT` on** -- `passed=1 failed=0
incomplete=0`, nothing changed in the test or the split. That exception was the stated correctness
blocker for making the split the default, and it was a casualty of the double leaks.

- [x] **1. THE A.4/A.5 SPLIT IS THE DEFAULT since 2026-08-18.** `testcPLQ` 8 passed / 0 failed in
      2188 s against 4728 s and one ERROR; full suite 332 / 0 with it on. `CCA2_NO_A45_SPLIT` opts
      out. The slow bucket went ~92 -> ~113 minutes, which is the bucket cost the standing rule
      says to accept.

- [x] **2. `region.impliedBy` over a region with a CURVED facet -- DONE 2026-08-18.**
      `region.holdsOn` + `maxAffineOverRegion` take the max over the region ITSELF (vertices, plus
      arc tangencies found in closed form) instead of over its linear relaxation. Measured on three
      folds: `quadFacet_exactAnotInB` 63 -> 41, fold-3 cells 38 -> 36, merges 7 -> 9.
      **An earlier revert of this was MY ERROR** -- two helpers in the instance block called as
      static -- see `DECISIONS.md`.

- [x] **3. `noSharedFacet` -- MEASURED 2026-08-18 and mostly HONEST.** At fold 3, of 137
      same-function pairs, 74 do not touch or touch at a single POINT, where merging would be
      wrong. Only 11 are genuine misses: 4 share an affine hyperplane the symbolic test does not
      match, 7 meet along a CONIC that neither facet search identifies.
      **The open question moved to `unionIsExact`**: 52 pairs reach it and about 9 merge, so ~43
      are refused by the sound gate -- and whether those are right is NOT established. Two cells
      can share a facet, touch along a segment, and still have a non-convex union. Check a handful
      of them directly before optimising, or risk reading a correct refusal as a defect for the
      third time this session.

- [x] **4. The parallelogram's piece 9 -- DONE 2026-08-18, `BAD 0 of 10`.** The singular-quadratic
      overlap was the only defect; the "remaining 1%" was a grid reference that missed the vertex
      where the sup is attained. See `DECISIONS.md`.

- [x] **5. `RatPar`'s `V (:,2){mustBeNumeric}` -- DECIDED 2026-08-18: leave it.** It costs cells
      and time, not correctness; the default path no longer goes through the mesh; and the change
      is lattice-wide against a design the classes state deliberately. `DECISIONS.md` records what
      would reverse it.

### Measurements that stand (2026-08-16)

- [ ] **STEP 3's CROSS-PIECE MAXIMUM DOES NOT SCALE, and it is now the binding cost.** Measured
      2026-08-16 on `x*y` over `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}`, which the A.4/A.5 split turns
      into 6 pieces: Steps 1 and 2 take about 25 s for all six, and `functionNDomain.maxOfList`
      then takes **73 minutes** -- folds of 93, 294, 647, 1273 and 2087 s, with the cell count
      running 5, 14, 29, 45, 70, **86**. `SUPPORT_MATRIX.md` records the same shape for a pentagon
      (885 s, 41 regions).
      **86 cells is roughly ten times what the answer needs.** `f*` of `x*y` over a convex
      quadrilateral has a cone per vertex and a cell per edge -- of the order of a dozen. The
      surplus is adjacent cells carrying the SAME function that are never merged, which is
      `region.merge` / `unionIsExact` refusing a union it cannot certify across a conic edge.
      **NUMBERS RE-MEASURED 2026-08-26 -- the ones below are stale.** Same fixture today:
      14, 29, 42, 60, 63 -> **56 cells in 49 min**, against the 5/14/29/45/70/86 and 73 min recorded
      here from 2026-08-16. Intervening work has taken a third off both.
      **And the cause is now located.** The SPLIT is not the growth (`afterSplit` exceeds `in` by
      1,3,3,8,0 -- fifteen cells created in the whole run) and the merge works (fold 4: 68 -> 44).
      The growth is `mtimes`, the pairwise overlay, which is inherent. What should contain it is the
      merge, and its refusals break down as: `noSharedFacet` 2743 (pairs that are NOT ADJACENT --
      expected at 63 cells), and **`quadFacet_exactAnotInB` 374 against 37 successes**. That is the
      target: adjacent pairs sharing a CURVED facet whose union `unionIsExact` will not certify.
      **ANSWERED 2026-08-27: CONSERVATIVE, and located.** The exact closed-form bound that would
      decide these pairs (`maxAffineOverRegion`'s vertex-plus-tangency candidate set, added
      2026-08-17 for this very conservatism) DOES NOT RUN on most of them. Instrumented on PRect3:
      `tighten_polyhedral 6` (no conic facet -- fine) and **`tighten_unboundedVertex 2`** -- the
      guard "every vertex must be finite for the compactness argument" returns as soon as a vertex
      is at `intmax`, the ray marker. A conjugate is full of CONES, so this fires on a large share
      of exactly the cells with conic facets, and they fall back to the loose LP whose conservatism
      IS `exactAnotInB`.
      **The guard is sound, its argument is too narrow.** A linear form on an unbounded region need
      not attain its max -- true -- but it is decidable whenever the form is bounded on the region's
      RECESSION CONE, which is when `unionIsExact`'s question has a finite answer at all. Add the
      recession directions to the candidate set instead of rejecting the region;
      `pieceRecessionRays` already computes them.
      **HAND-CHECKED 2026-08-27: CONFIRMED conservative, with a minimal witness.** PRect3 turned
      out to have NO `quadFacet_exactAnotInB` pairs at all (only 7 cells; the 374-count fixture is
      the bigger quadrilateral above, and it is too expensive to reach one from cold -- over two
      hours without hitting a single case, killed). Built the smallest possible witness directly
      from `region`'s constructor instead: `A = {s1>=s2^2, s2>=0}` (unbounded), `B = {s1<=s2^2,
      s1>=0, s2>=0}`, shared on the parabola. `unionIsExact` returns `(false, 'exactAnotInB')`; by
      hand and by a 200,000-point sample (0 mismatches), `A union B` is exactly the convex
      quarter-plane `{s1>=0, s2>=0}`. **So the refusal is wrong.** The mechanism is sharper than
      "gives up": `maxAffineOverRegion([-1 0])` on A returns `st=1` (unbounded above, DECIDED, not
      undecided) because the unbounded-region LP fallback drops the conic facet entirely, leaving
      only `s2>=0` -- which alone bounds nothing on `s1`, so the relaxation reports unbounded even
      though the true region's one recession direction cannot decrease `s1`.
      `DECISIONS.md` 2026-08-27 (item 3).
      **THE RECESSION-CONE FIX IS REFUTED -- caught by a second probe on the SAME witness before
      it ever reached the test suite.** Built it (mirroring `pieceRecessionRays`'s straight-line
      extreme rays), and it made the witness above exact -- and then wrongly certified
      `maxAffineOverRegion(A,[0 1])` (max of `s2` over A) as `0, decided`, when the true answer is
      `+Inf` (`(100,10) in A`, and `s2` grows without bound the same way for any `t`). The
      straight-line recession-cone theorem is for POLYHEDRA; A's boundary is curved, and its
      recession cone (the single ray `(1,0)`) does not see growth ALONG the arc `(t^2,t)`, whose
      direction converges in ANGLE to `(1,0)` without ever being witnessed by a fixed straight
      ray. REVERTED before any test ran; nothing shipped. `DECISIONS.md` 2026-08-27 (item 1,
      REFUTED) has the full counterexample.
      **The correct ingredient, not yet built:** parametrize the conic's own unbounded branch
      directly (a parabola is one coordinate quadratic in the other) and check whether `cRow`
      composed along that parametrization is bounded above as the parameter -> +-infinity -- a
      different and more specific question than the straight-line recession cone, which stays
      valid and necessary only for the LINEAR facets. Hand-verify on the SAME witness first
      (`cRow=(0,1)` composed is `t`, unbounded -- correctly predicts the counterexample;
      `cRow=(-1,0)` composed is `-t^2`, bounded by 0 -- correctly predicts the case that was
      right) before writing it into `region.m`.
      `DECISIONS.md` 2026-08-27 (item 1), 2026-08-27 (item 1, REFUTED), 2026-08-26 (item 1, third
      pass).

      **FIXED, overnight 2026-08-27.** Built the arc-parametrization ingredient, then found a
      SECOND gap on paper before writing anything: the arc-only check misses growth along a
      STRAIGHT edge unrelated to the arc (a region with both an arc and a straight edge reaching
      infinity, where the true maximiser sits on the straight edge). The landed fix combines TWO
      independent sufficient conditions for unboundedness -- a straight recession ray (admitted by
      the conic's own recession condition too, not just the linear facets) and arc growth via
      `parabolaArcFrame` -- and only tightens when NEITHER fires. Validated against a true 2D
      brute-force oracle: 5000+ prototype cases plus ~150 built as real `region` objects, zero
      genuine disagreements (the few flagged were oracle-threshold artifacts or a pre-existing,
      unrelated `region.linearForm` fragility on degenerate random test input). Fast 309/0, normal
      12/0. `DECISIONS.md` 2026-08-27 (item 3, FIXED) has the derivation and the counterexample.
      **RE-MEASURED overnight against the ORIGINAL fixture: NO material win.**
      `.claude/step3cost.m` on `x*y` over `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}`: final cell count
      **58, against the baseline's 56** -- unchanged, if anything slightly worse; distinctF stayed
      at 8 both times, so the ~7x surplus this defect is named for is still there. The
      `quadFacet_exactAnotInB` refusal RATE did drop (139/2010 = 6.9% here vs 374/3236 = 11.6%
      before), but not enough to move the final count, and the two runs are not a clean diff
      (intermediate paired counts differ fold by fold, so Step 1/2's own use of
      `maxAffineOverRegion` likely produced slightly different upstream cells too).
      **The fix is correct and stays** (independently validated by a brute-force oracle and a
      hand-derived counterexample, fast/normal green) -- it closes a real soundness gap, it is
      just not what closes THIS scaling defect. `DECISIONS.md` 2026-08-27 (overnight, item 3 vs
      the ORIGINAL fixture).
      **CONFIRMED overnight 2026-08-27: 60% of the remaining refusals have 2+ curved facets, out
      of the fix's scope.** Instrumented `unionIsExact` directly (temp probe, reverted
      immediately after -- nothing shipped): of 141 bare `exactAnotInB` refusals on this fixture,
      `nq1=55` (39%, the fix's target, and it IS helping there), `nq2=84` (60%, outside scope),
      `nq3=2` (1%). `DECISIONS.md` 2026-08-27 (overnight, item 3 follow-up).
      **If picked up again:** the two-conic case is a MATERIALLY HARDER proof, not a small
      generalisation -- the region's true recession cone is the INTERSECTION of both conics'
      recessive-direction sets (not decidable from either alone), and there is no single global
      parameter analogous to `parabolaArcFrame`'s `u` once two arcs are both active constraints.
      Scope it as its own item rather than extending this one; DECISIONS.md's entry has a fuller
      sketch of what would be needed.
      **PARTIALLY CLOSED 2026-08-28, in two steps, both real (not the scaling defect itself).**
      (1) Mechanism 1 (straight recession ray) generalized to two conics -- sound, but MEASURED
      to move NOTHING on this fixture (139 `quadFacet_exactAnotInB` refusals before and after,
      byte-identical), confirming rather than refuting its own stated scope: it only confirms
      genuine UNBOUNDEDNESS, and this fixture's two-conic refusals are all bounded regions.
      (2) **A real, separate CLOSED-FORM result found and landed:** when BOTH curved facets are
      CONVEX (`trace(Qi)>0`) with DIFFERENT axes, the region is provably BOUNDED -- a rank<=1 PSD
      conic can only be receded along its own null direction, so two different axes share no
      common receding direction, full stop, no arc parametrization needed. This is NOT the
      general two-conic proof (same-axis and any-concave-facet configurations still abstain,
      correctly), but it closed a genuine, previously-undetected BUG: `maxAffineOverRegion`
      returned `Inf` -- unconditionally wrong, not conservative -- on `region([y^2-x,x^2-y],[x
      y])`, whose true max of `y-x` is `0.25` at a smooth tangency point, confirmed by brute
      force. `regionTest.twoDifferentAxisConvexConicsAreProvablyBounded` pins it on 6 directions.
      (3) Generalized further same session: fixed a SIGN GAP in mechanism 1 itself (only one
      sign of each conic's null direction was ever tested, latent since 2026-08-27), which let
      mechanism 3 generalize from "different axes only" to "any two convex facets" -- a
      same-axis pair receding in OPPOSITE senses also bounds the region. Second genuine,
      previously-undetected `Inf` bug closed (`region([x^2-y-1,x^2+y-1],[x y])`, true max(x)=1,
      max(y)=1). `regionTest.sameAxisOppositeSenseConvexConicsAreProvablyBounded` pins it.
      **MEASURED against this fixture three times (mechanism 1, mechanism 3, the sign fix) --
      byte-identical numbers every time (139 `quadFacet_exactAnotInB`, cells=58).**
      **THE REASON WHY, found 2026-08-28 (final) by instrumenting `holdsOn` directly rather than
      guessing a fourth mechanism: `unionIsExact`'s `exactAnotInB` refusals on this fixture are
      almost entirely (14/14 sampled on fold 2) GENUINELY DECIDED VIOLATIONS
      (`maxAffineOverRegion` returns a correct finite value that legitimately exceeds the
      bound), not undecided proof gaps -- and mostly `nq=1`, not `nq=2` (12/14). See
      `DECISIONS.md` 2026-08-28 (final).** This means 2026-08-25's original framing of this
      scaling defect as a merge-CERTIFICATE gap was itself the wrong diagnosis: the three fixes
      above are real, correct, and worth keeping (each closes a genuine `Inf` bug on its own
      witness), but NONE of them, nor any further two-conic proof, can close THIS fixture's
      surplus -- there is nothing wrong with the certificate to fix. **The scaling defect's real
      cause is upstream of `unionIsExact`**: the cells Step 1/2/`maxQuaPar` hand it are
      genuinely too fragmented to pairwise-union convexly. Any future attempt should instrument
      WHY the arrangement is this fine before touching `region.m` again.
      **FIRST WITNESS, 2026-08-28: a genuine SLIVER** (area 0.038, diameter 1.06 -- ~15x thinner
      than round) refused against four distinct candidates on its curved edge.
      **N-ARY HYPOTHESIS TESTED AND REFUTED, same session:** built a 2,000,000-sample brute-force
      oracle (Python, independent of MATLAB) and checked whether the sliver lies inside the
      UNION of all four candidates -- essentially ZERO coverage (1 of 799,743 in-sliver samples).
      **THE REAL MECHANISM, found by comparing bounding boxes: all four candidates share only
      ONE VERTEX with the sliver, not an edge.** They are different ARCS of the SAME underlying
      parabola (`region.merge`'s candidate test checks CURVE-EQUATION equality only, via
      `obj.ineqs(mq1(i)) == -obj2.ineqs(mq2(j))`), meeting the sliver at one point and running off
      in unrelated directions. `unionIsExact` is correctly refusing candidates that were never
      real neighbours -- not a proof gap, one-conic or N-ary, at all.
      **CORRECTION (checked before implementing, same session): the "cheap pre-filter" is NOT
      safe as stated -- do not implement it.** `region.merge`'s own HISTORY names the exact
      precedent: `quadCutsOther` (removed 2026-08-17) was the SAME CLASS of plausible-looking
      geometric pre-check ("refused when one region's conic met the other anywhere but at a
      vertex") and it wrongly refused valid merges, since two cells CAN union convexly while
      sharing a curve only partially or elsewhere. `unionIsExact`'s soundness is purely
      algebraic, not dependent on segment adjacency, so an edge-sharing requirement could reject
      a merge `unionIsExact` would have correctly approved. Six refused samples is evidence, not
      a proof of safety. `DECISIONS.md` 2026-08-28 (correction) has the full argument. Any future
      attempt at a candidate-generation speedup here needs a soundness proof first, not another
      plausible heuristic.
      **RESOLVED, same session: the sliver has NO valid merge partner in this fold, and that is
      not a bug.** Captured every `exactAnotInB` refusal in fold 2 (14 total): all trace to just
      TWO cells, the sliver and one apparent edge-adjacent neighbour ("A2"). Checked directly
      whether they share a function value (`mergeL`'s own equality test, instrumented) -- **NO,
      confirmed 4/4 times.** So A2 is not a missed merge either; it is an ordinary different
      piece that happens to share an edge. That leaves the sliver's only same-function
      candidates as B1-B4 (the same-curve-different-arc false positives), none of which share a
      real boundary segment with it. **The sliver genuinely has no valid partner in fold 2's
      piece list** -- not because one was missed, but because none exists yet at this stage.
      **REFRAMES THE WHOLE ITEM.** `.claude/step3cost.m`'s own premise -- "cells above distinctF
      ought to have merged" -- assumes every function's true argmax region is a single connected
      convex piece. The sliver is a concrete counterexample to that assumption AT THIS FOLD: it
      may need a piece not yet produced, or may simply be one of several genuinely-separate
      components. **GENERALIZED, 2026-08-29: it is not an outlier.** A compact signature probe
      over the FULL 5-fold run's 141 `exactAnotInB` refusals found only 12 DISTINCT cells behind
      all of them, the top 5 accounting for 89%, and the original sliver persisting unmerged
      from fold 2 through fold 5. **The 58-vs-8 surplus is dominated by a handful (~5-6) of
      specific, persistently-stuck cells, not a broad, diffuse problem.** `unionIsExact`,
      candidate generation, and every boundedness mechanism have all been independently checked
      and cleared this session -- what remains is tracing these FEW cells' actual histories
      (which fold created each, what its true intended neighbour would have to be), a
      tractable, bounded investigation rather than a general redesign.
      `DECISIONS.md` 2026-08-29 (item 1, generalized) has the full signature breakdown.
      **ROOT CAUSE FOUND, 2026-08-29 (later): a genuine HIGH-DEGREE HUB VERTEX, present before
      any folding.** The sliver's own vertex is already shared by 4 cells in piece 1's own
      Step-2 output and 4 in piece 2's, INDEPENDENTLY, before fold 1 even runs -- plausibly one
      of the "cone per vertex" apexes `.claude/step3cost.m`'s own header names. **This is a
      structural consequence of pairwise-folding many pieces that share a hub vertex, not a
      missing merge**: each new piece's fan re-tests against every existing same-curve cell
      touching that point, and most were never going to succeed (the same-curve-different-arc
      false positives measured all session). **Reframes the fix, one more time and probably
      finally: not "find cell X's partner" (there may be none, if the correct decomposition near
      a fan genuinely needs several wedge pieces) but "does the pairwise-fold STRATEGY generate
      more candidates near a hub than the answer needs" -- e.g. resolving each hub vertex's fan
      once, grouped, rather than re-testing it against every incoming piece. A genuine
      architectural change to the fold order/strategy, not a bug fix. Not attempted -- this is
      the natural stopping point for item 1's diagnosis. `DECISIONS.md` 2026-08-29 (item 1, root
      cause) has the full trace.
      **PARTIAL MITIGATION LANDED, 2026-08-29: `unionIsExact` memoized.** A safe, verified,
      MODEST win -- not the fold-strategy fix above. `region.unionIsExact` is a pure function of
      its inputs, so caching its result cannot change any answer; measured 34% of calls were
      exact repeats (the SAME candidate re-tried in every later fold). Full 5-fold measurement:
      cell counts and every merge-tally count IDENTICAL to every prior run (correctness
      unaffected); TOTAL time 2186s vs the 2226-2546s range before. Real but partial because
      `unionIsExact` is only part of a fold's cost (`getVertices`, region construction, and
      pairing are the rest, untouched by this).
      **SECOND MITIGATION LANDED, same day: `getVertices` memoized too, bigger win.** Same
      pure-function argument (region's only properties are ineqs/nv/vx/vy/vars, and nv/vx/vy
      are exactly what this computes); measured 29% duplicate calls first. TOTAL time 2186s ->
      2008s; slow bucket 612s -> 543s.
      **THIRD MITIGATION LANDED, same day: `simplifyUnboundedRegion` too, highest duplicate rate
      of the three (37%).** Same argument (value class, no side effects, checked directly).
      TOTAL time 2008s -> **1830s**; slow bucket 543s -> **493s**. Cumulative from the 2944s
      original baseline: ~38% faster. Still IDENTICAL cell counts/tallies at every one of the
      three steps. **The fold-strategy redesign is still the actual fix** for the cell-count
      surplus itself (58, unchanged by all three) -- these reduce wasted recomputation, not the
      surplus. `mtimes` (untouched by any of the three) is now comparable to or larger than
      `maximumP` in later folds.
      **CHECKED AND RULED OUT, same day: `region.plus` (mtimes' own dominant cost) has ZERO
      redundancy** -- a 3-fold probe found 246 unique keys of 246 calls, confirming its inner
      loop tests every pair exactly once and cross-fold operands always differ. Not a
      memoization candidate; a geometric pre-filter would be a DIFFERENT, riskier kind of
      optimization (must never reject a genuinely non-empty pair), not attempted. This session's
      memoization thread stops here, on a clean negative result.
      `DECISIONS.md` 2026-08-29 ("landed" x3, "checked, ruled out") has all four measurements.
      **CHEAPEST POSSIBLE FIX CHECKED AND RULED OUT, 2026-08-30: extra `mergeL` passes on the
      FINAL fold's output find NOTHING.** 58 -> 58 -> 58 across two more calls beyond the two
      `maximumP` already makes. Rules out a merge-ORDER artifact; confirms the hub's same-function
      cells genuinely cannot pairwise-merge in any order because most pairs share only the single
      hub point, not a real edge. **The real fix needs an N-ARY simultaneous fan merge (its own
      new soundness proof -- `unionIsExact`'s pairwise argument does not generalize automatically)
      or avoiding the fragmentation upstream in Steps 1/2.** Both are genuine redesigns; neither
      attempted. `DECISIONS.md` 2026-08-30 has the measurement.

      **THE "WHERE TO START" BELOW IS STALE -- read this first (2026-08-26).** Merging after each
      fold is ALREADY what happens: `maxOfList` calls `maximumP(true)` per fold and `maximumP` calls
      `mergeL` twice. Measured on PRect3 with the existing `CCA2_TRACE_MAXP`:
      `[maxP] in=11 afterSplit=11 merge1=8 merge2=7` -- 11 cells in, 7 out.
      **The refusals are what to attack, and `region.mergeTally` names them:**
      `noSharedFacet 7, okLinear 3, emptyOperand 1` -- so 7 of 11 attempts die in `sharesFacet`,
      which is a LOCAL NECESSARY CONDITION that returns false when it cannot PROBE the edge, not a
      decision that the cells are unmergeable. Its own header says "a false here costs compactness,
      never correctness". Those pairs are never put to `unionIsExact` at all.
      **Do this instead:** let a failed probe fall through to `unionIsExact` rather than
      short-circuiting to "no" -- strictly safe, since it only asks a question currently skipped.
      **Iterate on PRect3** (44 s cross-piece max, 68 s biconjugateF), NOT PRect (1016 s, >3600 s).
      `DECISIONS.md` 2026-08-26 (item 1).

      ~~**Where to start:** merge same-function neighbours after each fold, not only at the end. Many
      of these cells are POLYHEDRAL (the vertex cones), where `unionIsExact` already decides
      exactly, so a large fraction should collapse without needing the conic case at all. Measure~~
      the cell count per fold before and after -- that is the number that predicts the time.
      **This is the FIRST of two things standing between the A.4/A.5 split and being the default**:
      `testcPLQ` goes from 1542 s to 4728 s with the split on. The second is that `testRectBiconj`
      then ERRORS -- a separate, undiagnosed failure, and the reason the split is opt-in.

### Then (2026-08-15, after bug 1 closed the last red)

- [x] **THE GENERAL CONVEX QUADRILATERAL — FIXED 2026-08-16, on the fourth attempt, by doing the
      split SYMBOLICALLY.** `splitTightTriangleSym` splits a triangle into sub-triangles on each of
      which cPLQ's own closed form for THAT sub-triangle's convex-edge count IS the convex
      envelope, and `plq_1p.triangulate` emits them as PIECES.
      **Measured:** `f*` of `x*y` over `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}` is exact at 10 of 10
      probe points against the vertex-attained sup, with no piece leaving a hole; the fully
      assembled Step 3 answer is exact too, at 8 of 8. Pinned by
      `cplqAdapterTest/generalQuadrilateralStep1IsTheEnvelopeNotAMinorant` (the envelope must
      exist, be a minorant, and be `>= 0` where `x*y >= 0`) and
      `generalQuadrilateralConjugateMatchesTheSup`.
      **What made attempt 4 work where attempt 3 hung:** exact symbolic arithmetic. A.4's cevian
      foot is irrational, and taking it from `convEnvCPLQ`'s doubles gives `2^53` denominators that
      grow to `1e25` downstream; carried symbolically the coordinates stay compact surds
      (`5/2 - sqrt(5)/2`, `3/2 - 3*sqrt(5)/10`).
      **A COST, not a defect, and it is Step 3's:** assembling the cross-piece maximum for this
      input takes **73 minutes** (folds of 93, 294, 647, 1273, 2087 s, cells running 5, 14, 29, 45,
      70, 86). The per-piece conjugates take about 25 s in total. See the next item.
      **The cost is paid by any x*y polygon whose triangles need the split**, because A.4's cevian
      foot is IRRATIONAL: a split sub-triangle has SURD coordinates, and every symbolic operation
      downstream then works in a quadratic extension instead of the rationals. Measured on
      `testcPLQ`, whose domains are general polygons carrying `x*y`: **1542 s with the split off
      (matching its historical 1427 s), over 3100 s with it on and still unfinished when stopped**,
      uncontended both times. Only two of its six domains even gain a piece (2 -> 3 and 1 -> 2), so
      this is the algebraic degree of the coordinates, not the piece count.
      **The split is OPT-IN, via `CCA2_A45_SPLIT`, and OFF by default** -- and off for a measured
      reason: with it on, `testcPLQ` takes 4728 s instead of 1542 s AND `testRectBiconj` ERRORS.
      So it trades a documented, LOUD failure on one domain shape for a new one on another, and
      until that is understood it cannot be the default. The two tests turn it on themselves.
      History of the three failed attempts follows.
      **Attempt 3 (2026-08-16), and why it failed.** The domain split was built exactly as attempt
      2's write-up prescribed, taking the sub-triangles from `convEnvCPLQ`'s own faces, and it
      WORKED at Step 1. It then turned the crash into a **HANG** -- the first conjugate ran 45+
      minutes with no output, stuck in a symbolic call behind 3.8 MB of `isAlways:TruthUnknown`
      warnings carrying denominators around `1e25`. The cause is inherent to taking the geometry
      from `convEnvCPLQ`: that routine is double precision, `sym` of a double is EXACT (denominator
      near `2^53`), and snapping the new vertices to the simplest rational within `1e-10` bounds
      the VERTEX denominators but not the downstream ones -- the conjugate is a rational function of
      those coordinates, so a few squarings carry `1e5` to `1e25`. A hang is worse than a crash.
      **Attempts 1 and 2 (2026-08-15), for the record.** Wiring the missing `nCE == 3` branch alone
      leaves the answer wrong, because the `nCE == 2` branch returns a MINORANT and not the
      ENVELOPE (measured: its envelope reaches `-0.2835` on a triangle where the truth is `>= 0`,
      so `f*(0,0) = 0.28647` for a truth of `0`). Routing the envelope through `convEnvCPLQ`
      instead raises `symbolic:coeffs:NotAPolynomial`, because A.4/A.5's faces are RATIONAL and
      **cPLQ's Step 2 has no rational-envelope branch at all**.
      Original framing follows.
      What was on record: `p.conj('cplq')` raises `MATLAB:badsubscript` because
      `plq_1p.convexEnvelope1` branches on `nCE == 0, 1, 2` and falls off the end, so
      `obj.envelope` stays EMPTY and `conjugate`'s `for i = 1:max(1,size(envelope,2))` indexes
      `envelope(1)`. All true, and CCA2 already has the missing algorithm ([COAP] A.5,
      `convEnvCPLQ`'s `splitThreeConvex`).
      **The wiring was written and it WORKS at Step 1.** Build the triangle as a one-face `QuaPol`
      carrying `x*y` (safe: reaching that line with an indefinite `q` means `isCanonicalXY` held,
      so `q` is EXACTLY `x*y`), call `convEnvCPLQ`, convert with `ratPolToPlq` and install the
      faces as envelope pieces -- `plq_1p.conjugate` already loops over them. Measured on the
      offending triangle: **4 envelope faces** come back (the A.5 split, each half then needing
      A.4's), two quadratic and two rational, all `<= x*y`. `conj` no longer raises.
      **AND THE ANSWER IS THEN WRONG, so it was reverted.** `triangulate` splits the test
      quadrilateral into piece 1 = `[2.5 1.5; 2 0; 0 0]` (`nE = 2`) and piece 2 =
      `[2.5 1.5; 0 0; 0.5 1]` (`nE = 3`). With the branch in:
        * **piece 2 gets 4 envelope faces and cPLQ's Step 2 returns ZERO conjugate cells for it**
          -- the new envelope is computed and then discarded, so the wiring buys nothing today;
        * **piece 1, untouched by any of this, is WRONG on its own through cPLQ's Step 2**: 6
          cells, `f*(0,0) = 0.28647` where the truth over its own triangle is `0`,
          `f*(0.5,1) = 1.00464` where it is `1`, and NOT COVERED at `(-1,0.5)` and `(-2,-2)`.
      So the crash was masking a **silent wrong answer**, and landing the wiring alone trades a
      loud refusal for it. That is why it is not committed.
      **The one measurement that tells you where to look.** That same triangle conjugated ON ITS
      OWN via `QuaPol.conj` is exact at 7 of 7 probe points -- because a single bounded triangle
      takes the NUMERIC route (`conjBoundedPolygon`), not cPLQ. **The numeric Step 2 is right on
      this input and the vendored symbolic one is not.** `assertStep3MatchesPieces` correctly does
      NOT fire: Step 3 agrees with Step 2: the fault is inside Step 2.
      **THE UNDERLYING DEFECT IS FOUND, and it is in STEP 1, not Step 2.** `plq_1p`'s `nCE == 2`
      branch applies [COAP]'s single-quadratic form to the WHOLE triangle. That form is a valid
      convex MINORANT but A.4 shows it is tight only over a sub-region -- and this branch never
      tests. Measured on `conv{(2.5,1.5),(2,0),(0,0)}` carrying `x*y`, it returns an envelope whose
      minimum over the triangle is **-0.2835** at `(1,0)`; on that triangle `x >= 0` and `y >= 0`,
      so `x*y >= 0`, the affine minorant `0` is admissible, and the TRUE envelope is `>= 0`
      everywhere. A too-small envelope gives a too-large conjugate: that is the `0.28647`.
      `convEnvCPLQ` on the same triangle returns 2 faces -- it does apply the split -- with
      minimum `0`.
      **And routing Step 1 through `convEnvCPLQ` does NOT fix it -- tried, measured, reverted.**
      `conjugateFunction`'s `nCE == 2` branch reads its envelope with `coeffs(...)` and matches
      monomials; `convEnvCPLQ`'s A.4/A.5 faces are RATIONAL, so it raises
      `symbolic:coeffs:NotAPolynomial`. **cPLQ's Step 2 has no rational-envelope branch at all.**
      **So the split belongs in the DOMAIN.** The route that already works for rational faces is
      `conjCPLQ`'s `conjEnvelopeViaCPLQ`, which hands each rational face to cPLQ as its own PIECE
      through `ratPolToPlq`. Do the same here: have `plq_1p.triangulate` split a 2- or
      3-convex-edge triangle into the A.4/A.5 sub-triangles and emit each as a piece, recursing
      while `splitTwoConvexEdges` still reports `needsSplit`. Every sub-piece is then a triangle on
      which cPLQ's own closed form IS tight, Step 2 is untouched, and every envelope stays
      polynomial.
      **The cost, which is why it was not done unattended:** `splitTwoConvexEdges`,
      `splitThreeConvex` and their helpers are file-local to `convEnvCPLQ.m`, so exposing them
      means moving a connected web of functions out of a well-tested file -- and `triangulate`
      feeds every Case C result. Design change plus a full re-verification, not a fix.

- [x] **The PARALLELOGRAM's `emptyResult` — TWO defects found and FIXED 2026-08-16**, taking its
      worst piece from **6 of 10** probe points wrong or uncovered to **2**, against a brute-force
      sup. Both are of this codebase's recurring kind, and both are general fixes rather than
      special cases.
      1. **`region.simplifyUnboundedRegion` declared any region with no finite VERTEX empty** --
         a half-plane, a slab, the whole plane. A half-plane is exactly what a TANGENT vertex
         produces: the cone there is built from the two edges meeting at it, and when those are
         tangent (an arc and its chord touching, how a curvilinear piece ends) both half-planes
         are the SAME one. Now refuted by a WITNESS. `regionTest/aHalfPlaneIsNotEmpty`.
      2. **The edge list, in bug 1's other form.** 3 vertices, 3 genuine edges, plus a conic
         touching one vertex: 4 constraints for 3 vertices, so the count called a BOUNDED region
         unbounded and it was built one edge cell short. `conjugateOfPiecePoly` now derives the
         list for any bounded piece the count mislabels, not just for a lens.
      `functionNDomainTest/aBoundedPieceWithATangentVertexConjugatesOntoTheWholePlane` pins it.
      **How it was found, worth reusing:** `f**` of a bounded domain is finite exactly on that
      domain and is a MAX, so EVERY per-piece conjugate must be finite there. Evaluating all 12
      groups at six interior points named the three bad ones in one cheap run.

- [x] **STALE -- already FIXED 2026-08-17, found overnight 2026-08-27 while sweeping the list.**
      The parallelogram's LAST 2 of 10 — `getInterior` on a SINGULAR quadratic. This entry is
      dated 2026-08-16; the fix (`functionNDomain.singularEdgeCut`, commit `699326d`, "Two
      fixes: the concave-conic certificate, and the SINGULAR-quadratic edge cut") landed the very
      next day and was never checked off. Confirmed still passing tonight:
      `functionNDomainTest/aBoundedPieceWithATangentVertexConjugatesOntoTheWholePlane` (the
      correct name for what this entry called `parallelogramPiece9` -- that is the fixture
      builder, not the test) is 1/0/0. The fix computes the interior cut from the KKT multiplier
      directly (`mu = <n, s - grad f(x*)>`, affine in `s`, no inversion needed) instead of the
      gradient-map elimination that fails when the Hessian is singular -- see `getInterior`'s own
      header for the derivation. No action needed; kept only as a note that this entry was stale.

### Bugs, in the order they should be taken (2026-08-15)

- [x] **BUG 1 -- FIXED 2026-08-15**, by the edge-list refactor mapped below. Everything from
      "Original framing" down is HISTORY, kept because the four failed attempts are recorded
      against it; the live write-up is the newest `DECISIONS.md` entry.
      **What it was.** Both slot conventions identify an edge by a VERTEX INDEX, and a LENS's two
      edges join the SAME pair of vertices -- so neither convention can name them apart, and no
      reassignment of slots to constraints ever could. That is why all four earlier attempts
      failed, and it is a stronger statement than "the count is the wrong test".
      **What it took.** `functionNDomain.edgeIndexList` derives `eIdx(j)` -- the constraint
      bounding the edge from vertex `j` to vertex `j+1` -- from the geometry (both endpoints on
      the constraint AND its own curve between them inside the region; the second half is
      necessary, since the lens's redundant conic passes through both vertices too).
      `region.getNormalConeVertexQ` takes the list as an optional argument;
      `region.getNormalConeEdgeQE` replaces `getNormalConeEdgeQ` and `Q3`, which turned out to be
      the same routine under two slot conventions; `getSubdiffVertexT2` and `T2Q` are identical on
      these inputs. **The refactor was smaller than this item predicted: three of the "four
      routines that move together" were two routines wearing four names, and both loops are in
      `conjugateOfPiecePoly` itself, so the list needs no field and no signature change.**
      **Scope.** Entered only when two constraints that each bound a genuine edge still share an
      edge number after `spreadCollidingEdges` (`edgesStillCollide`) -- the lens signature. No
      constraint is dropped, so the two unsound drops stay ruled out.
      **Measured.** Both half-lenses conjugate to 3 cells (2 vertex cones + the arc; the chord's
      cell is a ray and drops out), exact against a brute-force sup at 12 points, 0 wrong -- where
      the old code was `+inf` at `(0,0)` and `(-1,0.5)` and `0` at `(2,-1)` for a truth of `1/2`.
      Pinned by `functionNDomainTest/twoEdgesOnOneVertexPairAreBothKept` and
      `halfLensConjugateIsFiniteEverywhereAndExact`, ~10 s where reaching the same piece through
      `biconjugateTest` takes 10-40 minutes.
      Original framing follows.
      **TWO defects, and the description this item used to carry was wrong about the
      cause.** Measured 2026-08-15; full write-up in `DECISIONS.md`.
      **(1a) The pinned test no longer fails the way it says.** Since `conj` began returning a
      MESHED QuaPar for a bounded multi-face domain, `biconj`'s second conjugation is handed a
      CURVED QuaPar and `quaPolToPlq` refuses it -- the test ERRORS at `quaPolToPlq:curvedEdge`
      before any of (1b) is reached. Either the second conjugation learns to take a curved QuaPar,
      or `biconj` routes a bounded multi-face domain through the symbolic form on purpose.
      **(1b) The lens defect, reachable by forcing the symbolic route.** The cause is NOT the
      `size(d.ineqs,2) == d.nv` count: it is that `getEdgeNosInf` numbers an edge by one of its
      ENDPOINT VERTICES, and a lens has two edges joining the SAME pair, so they get the same
      number and the last-write-wins scatter destroys one. Fixing the numbering is NECESSARY and
      NOT SUFFICIENT -- measured: hand-build the lens with only its two genuine edges and
      `conjugateOfPiecePoly` still returns 1 cell where the piece needs 4. **Do the downstream half
      first** (4 cells for a bounded 2-vertex region with one curved edge, checked on that
      hand-built region, which needs no pipeline); the numbering fix is written up in
      `DECISIONS.md` and can be re-applied after.
      **PROGRESS 2026-08-15 (uncommitted at the time of writing): 5 of 7 probe points now right,
      was 0 of 7.** Three defects fixed, in the order they had to be:
        1. `spreadCollidingEdges` -- give a lens's two edges distinct numbers instead of letting
           the scatter destroy one. Scoped to fire only on that signature.
        2. `getNormalConeVertexQ` indexed its second constraint as `j+1` UNWRAPPED, so it raised
           `badsubscript` on any BOUNDED region -- which is why the only caller sent every bounded
           region to the POLYHEDRAL `getNormalConeVertex`, whose cones come from the CHORD and are
           wrong for a curved edge. Wrapped cyclically (identical to `j+1` for the unbounded
           layout, so nothing that worked changes), and the dispatch now asks whether a constraint
           is quadratic rather than comparing a constraint COUNT to the vertex count.
        3. `biconj` on a bounded multi-face domain takes its FIRST conjugate in symbolic form
           (`conjCPLQ(..., 'symbolic')`), because the second conjugation cannot take the curved
           MESH `conj` now returns and died at `quaPolToPlq:curvedEdge`.
      Unit level: the half-lens conjugates to 3 cells matching a brute-force sup at all 10 probe
      points (2 identical wrong cells before).
      **REMAINING, and it is a REFACTOR, not a fix.** `f**` is exact at 5 of 7 probe points; the
      other two come back `+Inf` -- a hole in the DOMAIN, not a wrong value, since `f**`'s domain
      is the INTERSECTION of the per-piece conjugate domains and one piece of `f*` still
      conjugates onto too small a set.

      **THE INDEXING CONTRACT, mapped 2026-08-15 -- this is what the refactor has to preserve.**
      Two conventions run through four routines at once, chosen by the COUNT
      `size(d.ineqs,2) == d.nv`:
        * BOUNDED layout: `nv` slots, edge `j` at `ineqs(j)`, `endNv = nv`, cones from
          `getNormalConeVertexQ` / `getNormalConeEdgeQ3` and `getSubdiffVertexT2`.
        * UNBOUNDED layout: `nv+1` slots with slot 1 reserved for the ray, edge `j` at
          `ineqs(j+1)`, `endNv = nv-1`, cones from `getNormalConeEdgeQ` and `getSubdiffVertexT2Q`.
        * Confirmed in the source: `getNormalConeEdgeQ` loops `j = 1:nv-1` and reads
          `slopeIneq(j+1, vertex j)`, so its output ROW `j` is the edge at slot `j+1`.
          `getNormalConeVertexQ` uses the same "vertex `j` lies on constraints `j` and `j+1`" rule
          (now wrapped cyclically, so it works in the bounded layout too).
      The loop variable `j` indexes `NCE`'s rows, `subdE`'s rows and `d.ineqs` SIMULTANEOUSLY, so
      the count cannot simply be replaced -- all four move together.

      **THE LENS NEEDS THE BOUNDED LAYOUT AND CANNOT GET IT.** It has `nv = 2` and 2 genuine edges,
      so the bounded layout is right; but it arrives with 5 constraints, so the count picks the
      unbounded one, `endNv = 1`, and only one edge cell is built.
      Two ways to give it the bounded layout, and the trap in each:
        1. DROP the three non-edge constraints. **Unsound** -- a constraint active at exactly one
           vertex of a convex region can still be essential; measured, see `DECISIONS.md`. It would
           be sound if guarded by a REDUNDANCY PROOF, but `redundantSubset` is an LP and cannot
           certify redundancy in the presence of a conic. A closed-form "maximise this linear form
           over the region's own boundary" test would -- the machinery exists in
           `verifyMaxIsExactSymbolically`, for QuaPar faces rather than `region`s.
        2. PARK them above slot `nv` and dispatch on boundedness instead of the count.
           `wasBounded` is now computed soundly (before `removeInfV`, which deletes the evidence).
           The trap: the `f.isQuad` chord rewrite iterates ALL slots, so a parked QUADRATIC
           constraint gets a chord derived from vertices it does not join -- the harm rule (2) in
           `conjugateOfPiecePoly` already records. Scope any parking to `~obj(i).f.isQuad`.
      Start from the hand-built lens probe -- it needs no pipeline and runs in seconds, and with
      the three fixes already in it produces the correct 3 cells when handed exactly its 2 edges.

- [x] **BUG 2 -- FIXED 2026-08-15.** `region.removeTangent` built a TANGENT LINE to a quadratic
      constraint at a vertex where that quadratic's GRADIENT VANISHES -- the apex of a cone, which
      is exactly where the Step 3 split conics of an unbounded fan meet. There is no tangent line
      there, every direction is tangent, and whatever it computed was meaningless; it then deleted
      a constraint matching that "tangent". On the 4-cone fan it deleted `-s1 <= 0` from the cell
      carrying `s1^2/4 + s2^2/2`, leaving only `{s2 <= 0, s2^2/2 - s1^2 <= 0, s1^2 - 2*s2^2 <= 0}`
      -- two constraints BLIND TO THE SIGN of `s1` -- so the region became symmetric under
      `s1 -> -s1` and claimed the mirror wedge. `f*(-3,-2.4)` is now **4.5**, the truth; it was
      5.130. `conjCPLQTest` 25 / 0.
      This is the SIBLING of the 2026-08-02 fix: `simplifyUnboundedRegion` fell into the same trap
      on the same input, and `region.witnessAwayFrom` was written for it. One input, one trap, two
      routines -- **a vanishing gradient at a cone's apex is a recurring failure mode here.**
      Cleared on the way, so nobody re-checks them: `redundantSubset` (certifies nothing there,
      correctly) and `simplifyUnboundedRegion` (leaves the constraint alone).
      `step3DropsCellsOnSomeUnboundedAssemblies` asserted the GATE firing; it no longer does, so
      that test is renamed `step3UnboundedAssemblyAgreesWithItsOwnPieces` and now pins what the
      gate protects.

- [x] **BUGS 3 and 4 -- FIXED 2026-08-15.** `convEnvUnbounded` computed only the AFFINE envelope
      and raised `convexAlongRay` as soon as `d'Qd > 0` along a ray. Two shapes are now derived and
      implemented, each with its proof in the source:
        * **WEDGE, one flat ray and one convex ray.** `co q` is `q` with its CROSS TERM deleted:
          `q(v) + alpha*g1 + beta*g2 + beta^2*A22/2`. A negative `A12` means `d1 + t*d2` recedes
          with negative curvature, so the envelope is `-inf` -- now reported rather than answered.
        * **HALF-STRIP convex along the ray, base edge Q-ORTHOGONAL to it.** `q` separates, so
          `co q = q(v1) + s*(q(v2)-q(v1)) + t*<grad q(v1),d> + t^2*(d'Q d)/2` -- the affine
          interpolant along the concave base plus the convex part along the ray.
      `w'Q d ~= 0` is deliberately not handled and is refused loudly.
      **`unboundedFaceTest` 18 / 0**, from 16 / 2. Fast bucket 204 / 0.

- [x] **BUG 5 -- FIXED 2026-08-15.** `splitTwoArcPiece` found no cut when the two arcs are
      ADJACENT: its two candidate chords join the arcs' facing endpoints, which for adjacent arcs
      ARE the arcs' own edges, so both chains come out with two vertices, the `numel(chain) < 3`
      guard skips them, and the piece was returned unsplit with one arc flattened to its chord.
      The `nv == 3` shared-vertex fallback did not apply at `nv = 5`.
      Generalised to `nv >= 4` with the ordinary DIAGONAL from the shared vertex to a non-adjacent
      one -- the two arcs leave that vertex in opposite directions around the boundary, so any such
      diagonal puts one arc in each half by construction. Same `insideStraightHull` guard as the
      existing candidates, and each half goes through `splitAtReflexVertex`.
      **Measured: the seeded sweep goes 17 exact / 0 wrong / 1 errored -> 18 / 0 / 0.**
      `maxQuaParTest` 29 / 0, fast bucket 203 / 0. Pinned by
      `arcVsArcSplitsTwoADJACENTArcsOnAPieceWithADiagonal` -- by VALUE, and with its own test
      because `arcVsArcMatchesGroundTruthOverRandomShifts` asserts `nWrong == 0` and would have
      counted this input in `nErr`.
      (The "two sub-arcs of the same conic" description this item used to carry was refuted by
      measurement; `DECISIONS.md` records it.)


- [x] **FIXED** `arcVsArcDoesNotCrashOnSeededQuadSplits` -- the last piece was a REFLEX vertex left by the
      bent two-arc cut: half-plane point location cannot represent a non-convex face, so the notch at the
      bend belonged to neither half. `splitAtReflexVertex` splits such a half along a diagonal. NOTE the
      retracted reasoning below -- equal areas and a bit-identical shared polyline do NOT imply two
      halves tile, because area says nothing about which side of a BENT boundary a point falls on.
      History follows. -- **both fixtures now ASSEMBLE** (three defects
      fixed: a duplicated consecutive vertex shifting the curve labels, `numel` where
      `size(...,1)` was meant in two index guards, and a STRAIGHT splitting curve being routed
      through the two-arc split, which flattened the inherited arc). Fixture 2 is exact on all
      1080 ring points. Fixture 1 has **one** point left, and it is a COVERAGE HOLE, not a wrong
      value: at `s = (0.998629534754574, -0.0523359562429444)` no face admits, and the three
      nearest faces all carry the right quadratic and miss by 0.0019, 0.0052 and 0.0072 in
      normalised conic units -- far above any tolerance. The point belongs to cell (2,1)
      (`g1` face 2 n `g2` face 1); both pieces from that cell exist (pieces 3 and 4) and neither
      covers it. **MEASURED, and it exonerates the two-arc split:** the two pieces from cell (2,1)
      share their cut polyline BIT-EXACTLY (both carry `M = 1.1254915141491897,
      0.074667480226358884` and the same two endpoints, to all 17 digits), their orientations are
      both CCW, and their polygon areas sum to the parent's exactly (0.03699 each way). So the
      pieces tile. The slit therefore appears AFTER the split -- in assembly: the vertex merge, the
      half-edge pairing, or the face edge-lists `orderEdges` builds. Next: take the assembled faces
      that miss the point (they carry the right quadratic and miss by 0.0019, 0.0052 and 0.0072 in
      normalised conic units) and compare each one's `P` list against the piece it came from.
      OLD NOTE, now fixed, kept for the trail: the orphan error reported a STRAIGHT boundary edge
      of piece 4 facing an IDENTICAL CURVED edge of its neighbour, at distance 0.
        * fixture 1: piece 4 `src[2 1]` (1.297862,0.278742)->(0.915534,-0.078641) straight, versus
          piece 5 `src[2 2]` curved.
        * fixture 2: piece 4 `src[1 1]` (1.163109,-0.285096)->(-1.244161,-1.161034) straight,
          versus piece 6 `src[1 2]` curved.
      So a piece that must carry TWO arcs -- its own operand's, plus the other operand's arc along
      the shared face boundary -- has one of them represented by its chord. matchHalfEdges then
      correctly refuses to pair a straight half-edge with a curved one, and the arrangement fails.
      Piece 4 has nv=4 with an arc already on edge 1, so the flattening happens where the SECOND
      arc is introduced: `clipPolyByConic` (including its corner-cut branch) or
      `clipPolyHalfPlaneCurved`, not in `splitTwoArcPiece` -- three fixes there (index re-location,
      degenerate-piece filter, polyline cut when no chord is interior) each held and none changed
      this symptom. Next: instrument where piece 4 acquires its arc, and check whether the cut
      conic was applied as an EDGE or silently dropped to a chord.
- [ ] `splitTwoArcLens` refuses when the cut `A -> M -> B` leaves the cell (the seeded far-field
      fixture). The two arcs there join corners on OPPOSITE branches of their parabolas, so the arc
      between them swings far out and the polyline exits the cell; a different subdivision is needed.
      **Still no reproducer (checked again 2026-08-28, two passes).** First pass wrongly guessed
      the sign convention; corrected by reading `maxQuaPar.m:917`: `curveEc` is SIGN-NORMALIZED so
      `evalConic(curveEc,.) > 0` on the face's OWN interior (not `<=0` as first assumed). Rebuilt
      three hand-picked two-parabola lenses with that convention right -- (a) `y>=x^2` meeting
      `y<=-(x-4)^2+20`, (b) `x>=y^2` meeting `y>=x^2` (a 90-degree axis mismatch), (c) a
      curvature-mismatched pair (`y=x^2/2` vs `x=y^2/50`, corners near the origin and far out at
      `(5.85,17.1)`) -- and NONE trigger `splitTwoArcLens`'s guard: `M` and both edge midpoints
      stay strictly positive on both conics every time, well inside the lens. So the failure is
      not "any two opposite-opening parabolas", nor just an axis or curvature mismatch by itself --
      whatever makes the real fixture's `M` swing outside is a more specific geometric
      configuration than these three tried. The real `curveEc` pair is produced by
      `clipPolyByConic`/`poly.curveEc` from an actual `maxQuaPar` fold, which is what a faithful
      reproducer needs to trace instead of guessing conic equations by hand -- not attempted
      further this session (three configurations is a reasonable bounded search; a fourth guess
      without new information would be thrashing, AI/CLAUDE.md sec 5).
      **REPRODUCER FOUND, 2026-08-28 (later): a SYSTEMATIC random search (not another hand
      guess) hit it on the FIRST trial.** Two upward-opening parabolas of very different
      curvature (`y = 1.174 x^2 - 6.348 x + 9.975` and `y = 0.518 x^2 + 2.139 x - 0.444`),
      intersecting at `A=(1.374,3.470)` (near the smaller-curvature parabola's own vertex
      region) and `B=(11.571,93.667)` (far out, where the LARGER-curvature term dominates and
      the two curves' values diverge sharply). `splitTwoArcLens`'s own guard fires:
      `evalConic` at `M`, `A->M`'s midpoint, and `M->B`'s midpoint is `< -sc` for the arc with
      the SMALLER curvature -- confirmed against the exact guard condition copied verbatim from
      `maxQuaPar.m:2033-2040`. Four more hits found in the same 20000-trial search (seed 42,
      script not committed -- trivial to rebuild: two random `y=a(x-h)^2+k` or `x=a(y-h)^2+k`
      parabolas per trial, intersect symbolically, test every point-pair as (A,B)).
      **What makes it fire, now that a real example exists to generalise from:** a LARGE
      curvature ratio between the two parabolas (here ~2.3x, and the search's other 4 hits are
      similar), NOT axis mismatch or opposite-opening (my three earlier failed attempts all
      used comparable curvatures). The smaller-curvature arc's own parametrization stretches
      much further in `u` between the same two endpoints than the larger-curvature one's does,
      so their arc-midpoints `P1`,`P2` diverge sharply and `M = (P1+P2)/2` lands far from
      either true arc.
      **Not yet wired into a maxQuaParTest regression** -- `splitTwoArcLens` is a local function
      in `maxQuaPar.m`, not callable from outside it, so a real test needs to reach it through
      the public `maxQuaPar(g1,g2)` entry point with two operands engineered to produce this
      exact two-vertex lens, which was not attempted (nontrivial reverse-construction, separate
      from finding the failing geometry itself).
      **ONE BOUNDED ATTEMPT, 2026-08-29: confirms the reverse-construction is real work, not
      just tedium.** Built `g1`,`g2` as simple bounded triangles sharing the curved edge A-B
      (`ec1`,`ec2` from the found reproducer) via `QuaPar`'s own constructor -- rejected outright,
      `maxQuaPar:notFullDomain`, "F has a 0 entry". Every `maxQuaPar` input must be FULL-DOMAIN
      (finite everywhere, matching a real conjugate's `+inf` outside its proper domain via
      genuine unbounded rays, not a plain bounded polygon) -- confirmed against
      `maxQuaParTest.m`'s own unbounded examples (`dirIn`/`dirOut` apex-plus-direction markers).
      Getting this right needs the SAME ray-construction machinery those examples use, correctly
      oriented so clipping produces the exact two-vertex lens -- real, careful, separate work,
      not attempted further (the risk flagged before: an incorrectly-built full-domain input
      could silently pass or fail for the WRONG reason, which is worse than no test).
      **SECOND ATTEMPT, same session: studied the unbounded-QuaPar convention further and found
      a DEEPER obstacle, not just a representation detail.** `splitTwoArcLens`'s own two-vertex
      lens (corners A, B) is itself BOUNDED -- it is an INTERMEDIATE cell `clipByFace` produces
      while clipping one genuinely unbounded, full-domain operand against another, not something
      a caller hands in directly. Building `g1`,`g2` correctly (each covering the whole plane,
      via real rays) is necessary but not SUFFICIENT: there is no direct control over where
      `clipByFace`'s internal clipping happens to leave a two-vertex remainder with exactly the
      found `ec1`,`ec2` -- that emerges from the interaction of both operands' full boundaries,
      not from A, B, ec1, ec2 alone. Constructing operands that provably produce THIS exact
      intermediate cell (rather than some other shape that happens to also error) needs tracing
      `clipByFace` forward from candidate operands, not backward from the target cell -- a
      genuinely different, larger task than "build two curved triangles". Stopping here rather
      than guess further; the standalone reproducer (verified against the real guard condition)
      remains the deliverable for this item.
      **THIRD ATTEMPT, 2026-08-30, a genuinely different method: SEARCH real operands instead of
      hand-building them. 0/500 hits.** Reused `maxQuaParTest`'s own trusted machinery
      (`buildCurvedG1G2`/`curvedFixtureTriangles`, real `conjPieceCPLQ` output) and swept 500
      random affine perturbations (shift+rotation+scale) of the fixture, checking every
      `maxQuaPar` error for `splitTwoArcLens` in its stack -- none hit. **Likely structural**:
      this family always conjugates `f=xy`, whose Hessian is fixed rank-1, so both arcs'
      intrinsic curvature share one source and may be unable to reach the large curvature RATIO
      (~2.3x) the real reproducer needs, however the triangle is transformed. `DECISIONS.md`
      2026-08-30 (splitTwoArcLens, third attempt) has the reasoning. Three attempts, three
      different methods, one answer: needs the ray/apex full-domain construction machinery built
      from scratch, not a search over an existing single-function fixture family. Stopping here.
- [x] **FIXED** `twoCurvedWhereTheSplitCurveCrossesAnArc` -- the test passes, and `MAXQP_ASSERT=2`
      is clean on that fixture. The four non-compact-arc-piece findings this entry recorded were
      closed by `QuaPar.chordCuts` (2026-08-13) and the corrected chord derivation in
      `pieceRecessionRays` (2026-08-14); the entry outlived them. Original text follows.
      -- 2 of 68 sample points wrong, and it is the same defect as the item below, not a separate
      one. At both bad points `QuaPar.eval` reported `region 0`, i.e. SEVERAL faces admitted them
      (`(-3.9811,0.6115)` gave 0.468 and `(-5.0954,0.1351)` gave 1.229, truth 0 at both); the
      verifier named face 13, carrying g1 face 2, beaten by g2 faces 5 and 6 by +inf along one of
      its own rays; and `MAXQP_ASSERT=2` listed four pieces with non-compact constraint regions
      (`src[2 1]` and three `src[6 1]`), all BOUNDED arc-pieces.
- [x] **Step 4 of the plan: bounded arc-pieces whose CONSTRAINT region is non-compact** --
      **the REPRESENTATION question this item poses was answered on 2026-08-13, and the item as
      written is stale.** Neither (A) nor (B) below is what happened: option B was checked before
      implementing and refuted (the chord runs through the NEIGHBOUR's interior, so making it a
      real edge splits the neighbour and leaves the offending face's own edge list unchanged -- see
      `DECISIONS.md`). The answer is that the chord is **derived per face**, which is what
      `QuaPar.chordCuts` does, and it resolved the whole far-field defect.
      What was left of this item on 2026-08-14 was that `pieceRecessionRays` -- the piece-level
      analogue, which decides the same question before assembly -- was still using the weaker rule:
      it read the chord's side off the piece's other VERTICES and had no gate on when a chord may
      be emitted at all. Corrected to read both off the conic. **Unrun**; the residual
      `MAXQP_ASSERT=2` findings on `src [1 2]`, `[1 6]` and on
      `twoCurvedWhereTheSplitCurveCrossesAnArc` are the measurement that closes this.
      Original text follows, for the record: for a piece whose arc is CONCAVE towards it, the
      constraint set is a wedge intersected with the OUTSIDE of a parabola; a cut parallel to the
      chord leaves the arc-side sliver still receding along the two side edges' own direction,
      which neither the cut nor the parabola blocks. The two options considered were (A) let a face
      list a REDUNDANT bounding conic in its `P`, and (B) make the chord a REAL edge by splitting
      the neighbour across the arc.
- [x] The verifier does not prove the faces COVER the plane; `partitionReport` only samples.
      **`verifyFacesCoverThePlane` (2026-08-14) does, in four checks on the constraint data:** every
      edge separates two faces; every edge lies inside both of them; no face's constraint region
      has boundary anywhere but on its own edges; no face is squeezed onto a curve. Together they
      force the boundary of the union of the faces to be empty, so the union is the plane. The
      argument is in the file's header, the summary in `SUPPORT_MATRIX.md` 4.3. **Unrun.**


_Seeded 2026-08-02 at the start of the overnight run, from the task given when it
was launched. The repository had no TODO.md; the acceptance criterion is precise
(three named tests green), so the run works from this list._

## 2026-08-13: four arc-vs-arc failures pinned from ORDINARY polygon splits (not the shift fixture)

The arc-vs-arc defect does not need the translated-triangle fixture: `f = x*y` over a quadrilateral
handed in as two triangles either side of a diagonal reaches it, with a closed-form reference that
owes nothing to the pipeline (`supBilinearOverPoly` over each triangle; `sup` over a union is the
max of the sups, so overlapping triangles are fine). All four are now RED tests in `maxQuaParTest`,
vertices written out literally:

- `unitSquareSplitByItsDiagonalIsExactNearTheArc` — the SMALLEST failing case: `x*y` on `[0,1]^2`
  split by the main diagonal. Truth is closed form and POLYHEDRAL (`max(0,s1,s2,s1+s2-1)`, since the
  objective is bilinear on a box), so the two operands' arcs must cancel in the max; they do not.
  Wrong at 17 of 1080 ring points, ALL at radius <= 1 — worst 0.437 at `s=(0.45399,0.891007)`, where
  it returns the `s1` face (bounded by the arc conic `(s1+s2)^2/4 - s1 = 0`) in `s2` territory.
  So on this fixture the over-extension is NEAR the arcs, not in the far field.
- `arcVsArcIsExactFarFromTheArcsOnASeededQuadSplit` — the FAR-FIELD form, worst 5.28e4 at
  `s = 100*(cos(-57deg), sin(-57deg))`: got 53126, truth 309. Error growing like `|s|^2` is the
  non-compact bounded arc-piece signature of the MAJOR FINDING below.
- `arcVsArcDoesNotCrashOnSeededQuadSplits` — two fixtures, both `MATLAB:badsubscript` at the same
  site, `splitTwoArcPiece` line 2216 (via `clipPolyByConic` -> `clipByFace`). A raw badsubscript is
  never a designed refusal, so this is a bug independent of how the far-field work turns out.
- `arcVsArcRefusesAnUnboundedTwoArcSplit` — `maxQuaPar:notImplemented` from `splitCell`: "an
  UNBOUNDED half carries both the inherited arc and the splitting curve". Per FARFIELD_FIX_PLAN.md
  Phase 3 that guard is a bug detector, not a supported-input error.

`maxQuaParTest` is now 18 pass / 6 fail (the 2 pre-existing reds plus these 4), 40 s, still fast.

## MAJOR FINDING 2026-08-04: arc-vs-arc results are only LOCALLY correct (wrong far from the arcs)

The random-shift sweep's "silent WRONG" is NOT a handful of edge cases -- it is PERVASIVE, and it
afflicts even the two pins this session marked FIXED. Measured directly (g.eval vs pointwise-max on
rings around the origin):
  - `[-1 0.75]`: 0/60 wrong on a radius-8 ring, **2/60 wrong on radius 30** (worst 6.1).
  - `[2 -0.5]`:  2/60 wrong on radius 8,     **11/60 wrong on radius 30** (worst 33.8).
Their suite tests pass only because `curvedSamplePoints` samples NEAR the arcs. So `maxQuaPar` on two
curved operands is correct locally and WRONG in the far field, generally.

ROOT: a quadratic conjugate face is genuinely UNBOUNDED (e.g. g1 face 1 carries rays). The arc-vs-arc
subdivision, though, emits the sub-pieces of such a face as BOUNDED arc-pieces (one parabola arc plus
straight edges, no rays). A bounded piece left on the parabola's OPEN side is not a compact QuaPar
face: `QuaPar.eval` (locate by "every bounding conic, sign-oriented, <= tol") then admits points
arbitrarily far away into it, so it OVERLAPS the true face out there and, by eval's last-admitter-wins
rule, returns the wrong value. Verified on `[0.5 0.5]`: assembled FACE 15 carries g1-face-1's quadratic
over a tiny triangle near (-2,2) but admits (-3.98,0.61) two units away, overlapping the correct zero
face (FACE 9).

TRIED, reverted (all four break the suite or don't help):
  - piece-level compactness guard on triHalf output -- WRONG sign convention, rejected valid faces.
  - piece-level compactness guard on every bounded curved piece -- same, errored on all 3 fixtures.
  - post-assembly "bounded face admits a far point" check (QuaPar's own EC + P signs) -- CORRECT
    detection, but errors on `[-1 0.75]`/`[2 -0.5]` too, because they ARE non-compact far out.
  - post-assembly "two faces disagree at a far point" check -- also errors on all three, since the
    disagreement is real (they are wrong far out). Confirms the issue is systematic, not per-fixture.
A safety backstop that refuses non-compact results is thus correct but would turn every arc-vs-arc
result RED -- so the real fix must be UPSTREAM: give an unbounded quadratic face's arc sub-pieces their
RAY boundaries (so they stay unbounded and compact-as-faces) instead of closing them with straight
edges. That is a rework of the arc-vs-arc clip/split, not a guard.

NB: this reframes the whole session. The six arc bugs fixed are real and make the results correct
NEAR the arcs; but "18/2, two pins fixed" means "fixed where the tests sample", not "fixed".

## Genericity of the arc-vs-arc fixes (measured 2026-08-03)

`arcVsArcMatchesGroundTruthOverRandomShifts` sweeps seeded random shifts (not the 3 hand-picked
ones). Over 60 shifts: **65% exact, 15% assemble-to-WRONG, 20% error.** So the fixes generalise
well (not ad-hoc), but the 15% wrong is the TOP next task: it is the SAME pre-existing far-field
over-extension as [0.5 0.5]'s residual 2/68 (a decided unbounded polyhedral cell reaching past its
g1 face near a mesh vertex), which the arc-assembly fixes newly EXPOSE by letting those cases
assemble instead of erroring. conj's verification catches it in production. Fixing that far-field
coverage bug should turn most of the 15% wrong (and [0.5 0.5]'s 2/68) green.

## Status 2026-08-03 (Opus 4.8) -- 18/1, [0.5 0.5] now assembles

Three arc-vs-arc pins were red; **TWO fully fixed**, the third now ASSEMBLES (was erroring) and is
red only by a 2/68 value error from a SEPARATE pre-existing far-field bug. `maxQuaParTest` 16/3 ->
**18/1** throughout, no regression (random-quadrilateral sweep + all arrangement tests green).
Commits on `overnight/2026-08-02`: `96aad61` (T-junction), `53fc9fd` (assignSide), `be1a31f`
(arcEdge off-by-one + escape-to-infinity split + dedup normalisation), `3e1a6b2` (triangle
two-arc split). The [0.5 0.5] fixture chained ~9 distinct arc-handling bugs; six are fixed.

- **[0.5 0.5] `twoCurvedWhereTheSplitCurveCrossesAnArc` -- now ASSEMBLES, red by VALUE (2/68).**
  Six bugs fixed to get here (see the four commits above); the remaining error is UNRELATED to the
  arc machinery. The 2 wrong samples are `(-3.98109,0.61148)` and `(-5.09537,0.13508)` -- both far
  out in the lower-left, both essentially ON g1's mesh vertex `V5=(-3.98249,0.610504)` and its
  face-3/face-4 boundary ray (dir `(-0.857,-0.514)`). Piece `src[4 3]` (a DECIDED unbounded
  polyhedral strip, untouched by any session change) over-extends a hair past g1's edge 4 there and
  evaluates g1-face-4's quadratic (~0.47/1.23) where g1's real value via the neighbouring face is 0.
  A PRE-EXISTING far-field coverage/precision issue at a g1 vertex, newly reachable now that assembly
  completes; conj's own verification would catch such a result in production. NEXT: why
  facePoly(g1,4)/clipByFace lets src[4 3] cross g1's edge 4 near V5. The long "arc-clip mirror"
  analysis below is SUPERSEDED -- the mirror (src[6 2]) is now built correctly.

- **[-1 0.75] `twoCurvedThatMustAssembleAcrossRays` -- FIXED.** Was a cross-piece T-JUNCTION:
  a decided cone's ray ran to infinity while the neighbour side was split (segment+ray) at a
  point P lying exactly on that ray (perp ~2e-15). `insertGlobalPassthrough` now re-inserts
  every piece vertex that lands on another piece's ray/segment; a companion fix to
  `raySideVector` (its old adjacent-vertex representative is COLLINEAR with a just-subdivided
  ray, so `oppositeSides` decided on sign noise -- now falls back to the other ray's direction).
  Green and exact at all 68 ground-truth samples.

- **[2 -0.5] `twoCurvedAssembleAcrossRaysSecondFixture` -- FIXED (`53fc9fd`).** Was NOT decideWinner
  (that correctly flags the cell undecided) -- it was `assignSide` reading the winner at `piece.V(2,:)`,
  which for splitCell's UNBOUNDED "rest" piece (curveAfter=1) is a CROSSING point where diff~0 and its
  sign is noise, so both halves could get the same winner. Now reads the vertex farthest from
  `{diff=0}`. Exact at all 68 samples. (The long note below about decideWinner/parabola-to-infinity was
  the WRONG lead -- the real cell here is a strip cut by a LINE at two finite boundary points.)

- **[0.5 0.5] `twoCurvedWhereTheSplitCurveCrossesAnArc` -- SUPERSEDED, now FIXED (see line ~1401;
  `QuaPar.chordCuts` 2026-08-13 + `pieceRecessionRays` 2026-08-14). Kept below for the diagnosis
  that led there.** Piece 5 (`src[2 2]` =
  g1f2 n g2f2, an **arc-vs-arc** clip via `clipPolyByConic`) emits an unmatched ray on `x+y=0`
  (g1's face-2/face-6 edge) from apex `(-2.03125,2.03125)`, dir `(-1,1)`. The clip is CORRECT, not
  over-extended: the sign data at that cell is `evalConic(Ecut)@V = [0, -0.046, -0.015]`, so the
  vertex `(-2,2)` sits on the g2-face-1 (discard) side and the far ray is genuinely g2-face-2. What
  is missing is the MIRROR piece across `x+y=0`: `src[6 2]` = g1f6 n g2f2. g1 face 6 does NOT carry
  the arc (the arc edge borders g1 faces 2 and 1 only), so `src[6 2]` is clipped by the
  **arc-vs-HALF-PLANE** path (`clipPolyHalfPlaneCurved`), a DIFFERENT code path -- and it does not
  produce the matching `x+y=0` ray (its `src[6 2]` pieces 15/16 sit at `x+y=0.5` and `1.0`, and
  piece 15 has a bounded vertex exactly at `(-2.03125,2.03125)` but no ray along `(-1,1)`). So piece
  5's ray has no neighbour because the two clip paths DISAGREE along the shared g1 mesh edge.
  Post-passthrough T-junction scan finds nothing collinear with the ray, confirming it is not a
  subdivision miss. FIX DIRECTION (not yet done -- delicate, high regression risk to the arc
  machinery): make the arc-vs-half-plane clip of g2f2 by g1f6 produce the same `x+y=0` boundary the
  arc-vs-arc clip of g2f2 by g1f2 does, i.e. reconcile `clipPolyHalfPlaneCurved` with
  `clipPolyByConic` along a shared straight g1 edge.
  DECISIVE CHECK NOW DONE (confirms "f6 clip under-producing"): g1's mesh for [0.5 0.5] has
  `E(6,:)=[2 7 0]` a RAY from V2=(-2,2) in dir (-0.707,0.707) lying on `x+y=0`, with `F(6,:)=[2 6]`
  -- so it is exactly the g1 face-2/face-6 boundary and it runs to INFINITY. Therefore `src[6 2]`
  (=g1f6 n g2f2) must have a matching `x+y=0` ray, and `clipByFace(facePoly(g1,6), facePoly(g2,2))`
  -- which swaps to clip the curved g2f2 by g1f6's straight half-planes -- is dropping it. Start in
  `clipPolyHalfPlaneCurved` / `clipByFace`'s swap path: the g2f2 boundary ray on `x+y=0` is being
  lost when g2f2 (curved) is clipped by g1f6's `x+y=0` half-plane (edge 6). g1 face 6 is a cone at
  (-2,2) between edge 6 (dir (-0.707,0.707), on x+y=0) and edge 7 (dir (0.707,0.707)).

  --- OBSOLETE (kept for history; superseded by the sharper diagnosis above) ---
  The old `decideWinner`/parabola-to-infinity write-up for [2 -0.5]. The T-junction
  fix lets assembly COMPLETE, uncovering the real defect: `decideWinner` wrongly declares an
  UNBOUNDED cell "decided". Traced with a pre-assembly coverage probe: cell `src[4 3]` stores
  `winner=g2`, but at `s=(-5.4843,1.5866)` -- an interior point of that cell -- `g1=0 > g2=-8.648`,
  so g1 wins there; symmetrically cell `src[6 5]` stores `winner=g1` while g2 wins at
  `s=(-1.1298,4.1007)`. 4/68 samples, worst 8.648.
  WHY `decideWinner` misses it: it proves domination by sampling finite vertices + the ARC midpoint
  + the ASYMPTOTIC sign along the two bounding rays (t->inf). But `diff = f1-f2` is a PARABOLIC
  quadratic (rank-1 Q -- splitCell already asserts this), so `{diff=0}` is a parabola whose two
  branches run to infinity ALONG Q's null direction, which lies strictly BETWEEN the cell's two
  bounding rays. That parabola bounds a region of the opposite winner entirely inside the cell,
  touching neither the finite boundary nor either bounding ray -- so every point `decideWinner`
  samples is on one side and it "decides" wrongly.
  Tried and REVERTED (kept out of the commit): sampling `diff` at the 1-D stationary point of each
  edge and each ray. Sound and strict, but it does NOT fix this -- the sign change is off the
  1-skeleton entirely. Confirmed neutral on all three fixtures.
  Note this cell is ALSO beyond `splitCell`: even if `decideWinner` returned undecided, `{diff=0}`
  here makes ZERO finite boundary crossings (both parabola branches escape to infinity inside the
  recession cone), and `splitCell` asserts exactly two. So a real fix needs BOTH a sound
  sign-over-the-cone test in `decideWinner` AND a `splitCell` that can cut a cell along a parabola
  that enters and leaves at infinity. Same family as next-step 2 in the session handoff (the
  Step-3 unbounded over-claim); do not attempt a probe-based patch.

## Next up (the [0.5 0.5] defect) -- SUPERSEDED, read the Status section above first

The "over-extended, should terminate" framing below is REFUTED by the sign data (see Status):
piece 5's ray is legitimately g2-face-2; the real bug is the missing g1f6 n g2f2 mirror from the
arc-vs-half-plane clip path. Kept below only for its geometric measurements.

- [x] **SUPERSEDED -- already closed, see line ~1401.** `twoCurvedWhereTheSplitCurveCrossesAnArc`
      is the same fixture this whole section discusses; it is green (`MAXQP_ASSERT=2` clean),
      closed by `QuaPar.chordCuts` (2026-08-13) and `pieceRecessionRays` (2026-08-14). The "Piece
      5 emits a RAY" framing below was the pre-fix diagnosis; kept for its geometric measurements,
      not as an open item.
      **Piece 5 (src `[2 2]`) emits a RAY where its boundary should terminate.**
      Localised to the cell, the edge and the reason; this was the last defect before the fix.

      * Its unmatched ray: apex `(-2.03125, 2.03125)`, direction `(-1,1)`, lying
        on the line `x+y=0`, which is g1's face-2/face-6 edge.
      * CORRECTION to an earlier note here: piece 16's ray is NOT its partner.
        That one lies on `x+y=0.5`, a parallel but different line, so the two were
        never meant to match.
      * Sampling across the apex, the local structure is three cells:
        `(2,2) | (6,1) | (6,2)`, with g2's ARC separating the last two. So piece
        5's neighbour along `x+y=0` is a `(6,1)` cell — and the `(6,1)` pieces
        (13 and 14) are BOUNDED slivers of area 0.008 and 0.004.
      * A bounded neighbour means the shared boundary is a finite SEGMENT, not a
        ray. So piece 5 is over-extended: its boundary along `x+y=0` should stop
        where g2's arc crosses that line a second time, and instead it runs to
        infinity. `matchHalfEdges` pairs rays with rays and segments with
        segments, so a ray facing a segment can never match — which is exactly
        the reported symptom.

      Next: find why `clipPolyByConic` cuts cell `(2,2)` at the first crossing
      with g2's arc but not at the second. Note the restriction of that conic to
      this ray came out with `A = 1.7e-18` — numerically degenerate, so the
      quadratic is treated as linear and yields ONE root. Check whether the true
      second crossing is being lost there, or whether it lies on a different
      boundary element of the cell.

## Done recently

- [x] `clipPolyByConic`: replaced the blanket refusal of a disconnecting curved
      cut with a real CONNECTIVITY TEST. Measured CONNECTED on all three
      fixtures, so the single-cell construction was right and "return both
      components" was not the fix — see Blocked note below for why two
      components would be unrepresentable anyway.
- [x] The clip and split stages now produce a VALID PARTITION on all three
      fixtures (`partitionReport` OK), where before they overlapped.
- [x] Fixed the no-crossing keep/drop decision, which evaluated the conic at the
      centroid of the vertices — not necessarily inside a cell bounded by an
      inward-bulging arc.
- [x] Retracted the hole/overlap evidence and made the partition diagnostic
      sound (it omitted the arc edge, so a curved sliver appeared to contain a
      point 3 units away).
- [x] Tagged pieces with their source `(k,l)` pair.

## Blocked

- Nothing is blocked on a decision. The note that "return both components" is
  unrepresentable stands as a finding, not a blocker: a separated survivor would
  need the cutting conic running to infinity, i.e. an unbounded curved edge,
  which QuaPar cannot hold. The connectivity test now refuses only that case,
  and it does not arise on these fixtures.

## Retired hypotheses — do not re-try

- The orphan ray is one physical ray covered to different extents (nothing lies
  on it).
- The cut must be restricted to the arc's own span (it must not; the straight
  half-planes are applied first).
- `clipPolyByConic` emits clockwise cells (a guard on every bounded output never
  fires).
- `polyL` is non-convex (all four curved faces measured convex).
- Pair `(6,1)` is spurious (it is a genuine thin cell; a 0.0625-spaced grid
  simply cannot see a cell of area 0.008 — the "zero intersection" reading was a
  resolution artifact, not evidence).

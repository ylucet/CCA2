# Decisions & Dead Ends

_The record of what has been tried and ruled out, so it is not tried again._

**Append-only.** Never delete an entry. If a decision is later overturned, strike
the heading through (`## ~~...~~`) and add a new entry saying what changed.

**Nothing here is an action item.** Near-term actions live in `TODO.md`. Dead ends,
reverted attempts, rejected approaches and refuted diagnoses live here.

**Read this before attacking a hard problem** — especially one that looks like it
has been attacked before.

Newest entries at the top.

> **Seeded 2026-08-13** by copying the dead-end passages out of `TODO.md`. `TODO.md`
> was left untouched because another session was working in this project at the
> time, so **the originals are still there and this file duplicates them**. When
> CCA2 is next idle, delete the copied passages from `TODO.md` — specifically the
> `## Retired hypotheses — do not re-try` section, the `TRIED, reverted` block under
> `MAJOR FINDING 2026-08-04`, the `ATTEMPTED TWICE AND REVERTED` text in the step-4
> item, the `Tried and REVERTED` paragraph in the obsolete `decideWinner` note, the
> `--- OBSOLETE ---` block, and the `SUPERSEDED` section at the end — and leave
> `TODO.md` as action items only.

---

## 2026-08-23 — the BICONJUGATE's face type is NOT rational: degree 1, 2 or 4; and `f**` is NOT C1

Three questions closed on the envelope side. `CONJ_FIELD_PROOF.md` §5.1 had computed ONE degree-4
ruled cell and left the general face type open; `RatPol.m`'s header still claims quadratic-over-linear
on a polyhedral subdivision, citing [COAP] Step 4 and [JOGO] Prop. 1.

**The ruled-cell face is ALGEBRAIC of degree <= 4, not rational.** On a 2-cell of `f**` coming from a
1-cell of `f*`, `f**(x) = <s,x> - f*(s)` with `s` solving TWO CONICS: the arc `C(s) = g_i - g_j = 0`
(fixed, rational) and the incidence

      Delta(s;x) = det[ x - grad g_i(s) ,  grad g_j(s) - grad g_i(s) ] = 0

whose coefficients are AFFINE in `x` (and which is symmetric in `i,j`: the two versions differ by
`det[u,u] = 0`). So `deg_{Q(x)}(z)` divides 4; and when the Galois group of the four intersections is
`S4` or `A4`, the only subgroups containing a point stabiliser are the stabiliser and the whole group,
so **degree 2 is impossible generically — it is 4 or 1**. Implicit form `Psi(x1,x2,z)`: degree 4 in
`z`, total degree <= 8 (Macaulay count on three ternary quadrics). Do NOT store `Psi`; it is far larger
than the parametrisation, can pick up extraneous factors, and discards the root selection.

`RatPol` is the case `H_i = 0` AND `H_j` rank 1 — a vertex branch against an edge branch, i.e. ONE
piece. There `Delta` degenerates to a LINE lying along the parabola's axis and meets the arc ONCE,
which is why the single-piece answer is rational and not a square root.

**Degree 2 IS reachable, and from two ADJACENT pieces.** Verified instance, integer data:

      f = x^2 - y^2  on T1 = [-1,1]x[0,3]        f = x^2 + y^2  on T2 = [-4,4]x[-3,0]
      continuous: both equal x^2 on the shared edge [-1,1]x{0}

`Q1 = diag(2,-2)` is indefinite, so `T1` has no interior branch. `f*` has a genuine edge where `T1`'s
VERTEX branch at `v = (1,3)` (`f(v) = -8`, rank 0) meets `T2`'s INTERIOR branch `|s|^2/4` (rank 2), on
the circle `(s1-2)^2 + (s2-6)^2 = 72` — central, non-degenerate. On the ruled cell

      f**(x) = 2x1 + 6x2 - 28 + 6 sqrt(2) sqrt( (x1-1)^2 + (x2-3)^2 )
      Psi    = (z - 2x1 - 6x2 + 28)^2 - 72[ (x1-1)^2 + (x2-3)^2 ]        degree 2 in z, irreducible

Checked against `sup_s <x,s> - f*(s)` at 9 points of the cell: max error `3.6e-15`; on the arc both
pieces attain `f*` to the last digit against a brute-force grid; `grad f** = s` exactly; the rank-one
Hessian formula below matches the analytic Hessian to `1.3e-15` at 12 points.

**Theorem 3 does NOT forbid this.** Its identity carries `det Q_i det Q_j` and is about the two
INTERIOR branches of an adjacent pair. A vertex-branch/interior-branch pair is outside its scope, so
the adjacency degeneracy never reaches it. Consequence: degree-2 faces are cheap — two adjacent pieces
— while degree-4 faces still need non-adjacent pieces.

**`f**` is NOT C1, and the gradient does NOT extend across the kink.** `|x1| + x2^2` is a legal convex
`QuaPol`: `f* = s2^2/4` on the strip `|s1| <= 1`, ONE 2-cell whose Hessian `diag(0,1/2)` is RANK 1, so
its image under `grad` is the LINE `x1 = 0` and `d f**(0,x2) = [-1,1] x {2x2}`. The kinks of `f**` are
exactly the images of the RANK-DEFICIENT dual faces — rank 1 gives a curve, rank 0 a point. That set is
lower-dimensional, which does NOT imply the gradient is continuous across it; an earlier claim in this
session that it did was wrong. The same example's two half-planes come from BOUNDARY edges of
`dom f*`, whose rulings are RAYS rather than segments; the value there is quadratic.

**The parametrisation, which is what should be stored.** With `g_k = 1/2 s'H_k s + a_k's + c_k`,
`H(lam) = (1-lam)H_i + lam H_j` and likewise `a, c`:

      x(s,lam) = (1-lam) grad g_i(s) + lam grad g_j(s) = grad g_lam(s)
      z(s,lam) = <s,x> - f*(s) = 1/2 s'H(lam)s - c(lam)          (g_i = g_j on the arc)
      grad f**(x) = s          Hess f**(x) = (Ju)(Ju)' / ( u' adj H(lam) u ),  u = grad g_j - grad g_i

Rank one with kernel `u`, PSD automatically (`H(lam)` is a convex combination of PSD matrices and the
2x2 adjugate of a PSD matrix is PSD), blowing up only where rulings converge — the corners of `f**`.
Plotting needs only the FORWARD map: sweep `s` along the arc and `lam` across `[0,1]`; no root-finding.

**Consequence for the return type.** `RatPol` stops being `biconj`'s return type and becomes the
degenerate single-piece case. Structure proposal: `TODO.md`, the `QuaCon` / `AlgCon` section.

Derived and checked this session with an independent numeric oracle; nothing was run in MATLAB (no
VPN). Both instances above are fully specified, so they re-check in MATLAB directly.

---

## 2026-08-20 (evaluation) — MEASURED: a mesh conjugate is competitive as a SCIP nonlinear constraint, and the only thing that ruins it is linear scanning

The Row 7 decision (H-form plus filtered predicates) was about to be taken with no number for the
thing the target application actually pays: SCIP calls a nonlinear constraint's value-and-gradient
millions of times per solve. Written in C rather than estimated, because MATLAB's scalar-loop
timings do not transfer to a solver's inner loop: `.claude/evalbench.c`, `gcc -O2 -march=native`,
minimum of three runs on the shared machine, value AND gradient on every call.

    baseline: a 20-node expression tree walked with a forward-mode gradient      44 ns
      -- the honest comparison: what SCIP already pays for an ordinary small nonlinear constraint

    A  linear scan over cells        9 / 81 / 1024 cells        38 / 130 / 1670 ns
    B  uniform bucket, O(1)          9 ... 1024 cells           24 ...    27 ns
    E  slab, two binary searches     9 ... 1024 cells           15 ...    26 ns
    F  cached cell + bucket, COHERENT queries (random walk, step 1e-3)         11.5 ns
    C  no mesh, max of per-piece closed-form sups   3 / 6 / 24 pieces   105 / 215 / 822 ns

**The answer is yes, with one condition.** An indexed mesh evaluates FASTER than SCIP's own walk
over a small expression, and the cost is flat in the cell count — so the measured cell counts (11
for the A.4 fixture, 41 for a pentagon, 60-86 for the quadrilateral) do not matter. A LINEAR SCAN
does matter: 130 ns at 81 cells is three times the baseline and it grows without bound.

**Read E with its caveat, which is why the recommendation is B+F and not E.** The synthetic mesh is
a grid of axis-aligned boxes, so the slab method's two binary searches are compares against stored
doubles and it never evaluates a cell predicate at all. On a real conic arrangement the y-search
inside a slab compares against ARCS, so each of its ~log2(N) steps costs a conic evaluation. The
bucket's flat ~25 ns is the transferable number, it is O(N) to build against sorting plus per-slab
arc lists, and F shows the last-cell cache is worth more than either: SCIP's queries move in small
steps, so the previous cell usually still holds and confirming it is four predicates.

**C is worth keeping even though it loses.** `f*(s) = max_k sup_{x in T_k} <s,x> - q_k(x)` with each
sup in closed form (interior critical point when the objective is concave and it is feasible, else
the three edges as 1-D quadratics), and the argmax IS the gradient by Danskin. It needs no
arrangement, no point location and no exact arithmetic, so it is both a differential oracle for the
mesh and a usable fallback before Row 7 exists.

**What the measurement settles about the architecture:** nothing on the evaluation path is exact or
symbolic. The exact kernel — degree-4 vertices included — is a BUILD-time cost, and the finished
mesh compiles once into flat double arrays. So the algebraic-number question and the SCIP
performance question are independent, and neither constrains the other.

Also recorded because they are consequences of convexity that a future session would otherwise
re-derive: the gradient at any point is a global linear underestimator, the max over a box is
attained at a corner (four evaluations), and the box minimum is the global minimum when it is
inside — so bound propagation, usually the expensive part, is nearly free. The gradient JUMPS
across cell boundaries while the value does not; either side is a valid subgradient, so separation
stays correct, but a smooth-NLP heuristic will feel it.

## 2026-08-20 (after T2b) — VERIFIED INDEPENDENTLY: a vertex of `f*` can have degree 4 with Galois group S4, so NO tower of square roots is enough

A parallel session left an untracked `CONJ_FIELD_PROOF.md` in the working tree claiming this. It
is not taken on trust and it is not the source of what follows: `.claude/s4VertexCheck.m` takes
only the INPUT from it and recomputes everything.

**The input.** Three triangles fanned from `(4,0)` inside `conv{(24,10),(-16,10),(4,-20)}`, each
carrying a POSITIVE DEFINITE rational quadratic — so every piece is convex and every piece's
conjugate is its interior-critical branch near the point of interest.

**Recomputed here, from `q*(s) = (s-beta)'Q^{-1}(s-beta)/2 - gamma` alone:**

    g1 = 2*s1^2 + s1*s2 + s2^2/2 + 4*s1 - 2*s2
    g2 = 2*s1^2 +         s2^2/2 + 2*s1 + 2*s2
    g3 = s1^2/2 - s1*s2 + s2^2          -   s2 + 3

matching that document's three exactly. Eliminating `s2` from `g1-g2` and `g1-g3` gives

    (3/2)*s1^4 - 12*s1^3 + 5*s1^2 + 80*s1 - 48      = (1/2) * (3t^4 - 24t^3 + 10t^2 + 160t - 96)

— the same quartic up to a constant. `factor` returns it unchanged, so it is IRREDUCIBLE over Q;
its resolvent cubic `9t^3 - 30t^2 - 2688t - 11008` factors no further; and the discriminant
`32853221376` has square root `181254.576...`, not an integer. Irreducible quartic + irreducible
resolvent + non-square discriminant is the standard certificate for **Galois group S4, order 24**.
24 is not a power of 2, so the point lies in NO iterated quadratic extension of Q.

**And it is a genuine vertex, not an artefact.** At the root `s1 = 0.608050881512364`, with
`s2 = 0.358525944978488`, all three branches agree to the last digit
(`g_k(p) = 2.736875828608988`), and each piece's argmax `Q_k^{-1}(p - beta_k)` is strictly INSIDE
its own triangle — barycentric `[0.739 0.121 0.140]`, `[0.764 0.129 0.107]`, `[0.674 0.188 0.138]`.
So three full-dimensional faces of `f*` meet there.

**What it costs and what it does not.**

* It does NOT overturn T2b (the multiquadratic type, committed the same day). The A.4/A.5 cevian
  feet are degree-2 and PRIMAL — `sqrt(30)/12 - sqrt(15)/6 + 5/4` is a vertex of a sub-triangle of
  the DOMAIN, and that is measured and unavoidable. The multiquadratic field is necessary.
* It does mean the multiquadratic field is NOT SUFFICIENT for the vertices of `f*`, and that no
  enlargement by square roots ever will be. The mechanism is plain once seen: a vertex where three
  faces meet is a base point of a pencil of two rational CONICS, and two conics meet in degree 4.
* The structural half is the useful part, and it is the other document's, not this entry's:
  face functions and edge conics of `f*` are claimed to be RATIONAL always, with only the vertex
  LIST irrational. That claim is not verified here.

**The design choice this forces** — recorded as open, since it is not one measurement away:
(a) exact real-algebraic arithmetic to degree 4; (b) store regions in H-form and never write a
vertex coordinate down, deciding predicates by resultants/Sturm instead; (c) keep the
multiquadratic type and refuse the degree-4 cases BY NAME, in `SUPPORT_MATRIX.md`'s discipline.

## 2026-08-20 (last) — T1 is REFUTED AS STATED: one quadratic extension is not enough. A SINGLE triangle needs two.

T1 decided the sym-free number type is `Q(sqrt(d))` -- one squarefree `d`, with mixing two of them
an ERROR rather than a promotion to a tower (`exactQ`, 2026-08-19). That decision was taken from
A.5's cevian foot `5/2 - sqrt(5)/2`, which is a single extension. It does not survive contact with
a second triangle, let alone a polygon.

**Measured (`.claude/t1RadicandProbe.m`), the radicands appearing in the SUB-TRIANGLE COORDINATES
`splitTightTriangleSym` produces:**

    conv{(0,0),(2,0),(2.5,1.5)}     2 sub-triangles     sqrt(5)
    conv{(2.5,1.5),(0,0),(0.5,1)}   4 sub-triangles     sqrt(15) AND sqrt(30)
    conv{(0,0),(1,0),(2,1)}         2 sub-triangles     sqrt(2)
    conv{(0,0),(1,1),(3,2)}         4 sub-triangles     sqrt(3)  AND sqrt(6)

The first two are the two triangles of ONE quadrilateral -- the A.4/A.5 fixture the whole Step 3
cost story is measured on. So:

* **One A.5 triangle already needs two extensions.** Its split recurses into two A.4 halves whose
  cevian slopes are `-sqrt(mh*mw)` for different edge pairs; nothing makes those two products
  differ by a square. `exactQ` would RAISE inside a single piece's Step 1.
* **Step 3 makes it worse by construction.** The cross-piece max subtracts a cell carrying
  `sqrt(5)` from a cell carrying `sqrt(15)`. That subtraction is the operation `exactQ` is designed
  to refuse.

**Why the refusal was the right design and still is.** Silently building a tower is how an exact
type turns back into a symbolic engine. What the measurement changes is the FIELD, not the rule:
the extensions that arise are square roots of SQUAREFREE INTEGERS, and products of them stay in the
same family (`sqrt(15)*sqrt(30) = 15*sqrt(2)`). So the type to build is MULTIQUADRATIC --
`Q(sqrt(p1),...,sqrt(pk))` over the primes actually seen, an element being a rational combination
of `sqrt(m)` for squarefree `m | p1...pk` -- not a general algebraic number field and not a tower
of arbitrary degree.

**And it keeps the two properties `exactQ` was built for.** Those `sqrt(m)` are linearly
INDEPENDENT over Q, so an element is zero exactly when every coefficient is zero -- zero-testing
stays trivial and exact, which is what `region`'s `isAlways(dt == 0)`-style questions need. Sign
then follows from refining a rational interval until it excludes zero, which terminates precisely
because zero-testing is exact. Neither needs floating point.

**What this does NOT overturn:** doubles are still refused (one ULP made a shared facet invisible
to `merge`), rational snapping is still refused (1e5 became 1e25), and `exactQ`'s other two rules
-- raise on int64 overflow, refuse `fromDouble` of something unrepresentable -- carry over
unchanged. `exactQ` is a correct implementation of a field that is too small; it is the base case
of the one to build, not a dead end.

## 2026-08-20 (later still) — the A.4 split is LANDED, at the ENVELOPE level, by reusing the split the domain route already had

The item filed as Blocked below is closed. `plq_1p.convexEnvelope1` now splits a piece whose
triangle needs it and installs the sub-triangles' envelopes as several FACES, so a piece built
DIRECTLY -- not through `triangulate` -- gets the envelope rather than a minorant.

**Measured, `{(0,0),(1,0),(2,1)}` with `f = x*y`, pinned by
`cplqAdapterTest/twoConvexEdgeTriangleEnvelopeIsTightNotAMinorant`:**

    envelope       >= 0 on the whole triangle   (was dipping to -0.0429)
    f*(0,0)        0                            (was 0.0429, the OVERSHOOT)
    f*(0,1)        0.125                        (was 0, the UNDERSHOOT of attempts 1-3)
    cost           about 80 s

**Four things had to be true at once, which is why three earlier attempts each failed differently:**

1. The cross-face max must be able to split a RATIONAL pair -- 2026-08-20, above.
2. `conjugateFunction` must dispatch each face on that FACE's geometry. The reverted attempt did
   half of this (`nCE`, `mE`, `cE`, `V` from a rebuilt face domain) and left the per-vertex
   COORDINATES reading `obj.d.polygon.vx(j)` -- the PARENT's triangle -- in three places, including
   the affine interpolant `a,b,c` the vertex cells subtract. With the rational max in place that
   attempt gets as far as `f*(0,0) = 0` and `f*(0,1) = 0.125` and is then WRONG at `s = (-2,-2)`:
   -1.9289 where the sup over the domain is 0, attained at the origin. Reading the face's own
   vertices is the rest of the fix.
3. `region.signEverywhere` and `region.getVertices` -- the two defects above.
4. **The geometry must not be a second implementation.** The attempt carried its own `a4Split`,
   `twoEdgeQuad` and `triFromVertices`. `splitTightTriangleSym` is the same construction, already
   measured, already refusing what it cannot certify, and already recursive through A.5 -- so the
   branch calls THAT and sends each sub-triangle back through `convexEnvelope1`, which classifies
   its own convex edges instead of inheriting the parent's. Also 13x faster on this fixture (81 s
   against 1041 s), and it fixes `nCE == 3` on the same path, where the dispatch used to fall off
   the end and leave an EMPTY envelope.

**One claim in `splitTightTriangleSym`'s header is now struck** and amended in place: it said an
envelope-level split "cannot work" because A.4/A.5's faces are rational. The obstacle was the
DISPATCH, not the faces. The DOMAIN split stays the route for a piece that comes through
`triangulate` -- it is cheaper, and it is what makes the sub-triangle's `nCE` its own.

**Still open, and NOT closed by this:** `testMaxMultiRegion/testPCE2` builds its fixture with
`plq_1piece`, whose envelope comes from the eta/`solveC` route rather than from this branch. That
red is the plq_1piece path, so it is T6's (the migration), not this item's.

## 2026-08-20 (later) — two `region` defects the A.4 split exposed, and one refuted diagnosis

Landing the A.4 split (below) stopped raising `NotAPolynomial` and started raising
`symbolic:kernel:DivisionByZero` out of `region.getVertices`, from inside `maximumP`. Two
independent defects, each now a fast unit test, and both are on the ordinary path — the split only
made them REACHABLE, it did not create them.

### 1. `isAlways` does not complete the square, so a tangency was read as undecided

Two conjugate cells that meet TANGENTIALLY have a perfect-square difference. Measured, on the
cells this split produces:

    isAlways((s1*s2)/2 - s2/2 - s1/2 + s1^2/4 + s2^2/4 + 1/4 >= 0)    UNKNOWN
    isAlways(((s1+s2) - 1)^2/4 >= 0)                                  TRUE

Same function; `simplifyFraction` leaves the first form, and `maxArray` asked about that. The
undecided answer makes the caller SPLIT a cell that must not be split, and on a cell whose vertices
all lie on the tangency line every vertex ties, so the vertex comparison then returned **f2, the
smaller one** — a wrong value, not a loose one.

`region.signEverywhere` now carries both things this test kept getting wrong — REAL variables (the
existing reason, recorded in `maxArray`'s header) and a retry on the SIMPLIFIED form — as one
routine used for the raw difference and for the cleared polynomial. The `simplify` runs only where
the alternative is splitting a cell, which costs more.

### 2. `simplify` was deciding whether a vertex EXISTS

The split boundary that same refusal emitted is `-(s1+s2-1)^2 <= 0`, which is VACUOUS -- and a line
meets that degenerate conic in a DOUBLE root, so `lineMeetsConicSym`'s radical is an unevaluated
form of `sqrt(0)`:

    t1 = 3 - 2*sqrt(2) - 2*((sqrt(2)-3/2)^2 - (sqrt(2)-2)^2 - sqrt(2) + 7/4)^(1/2)
       = 0.1715728752538099            <- `double` reads it; `simplify` RAISES on it

`getVertices` called `simplify` on the way to STORING the coordinate, so the whole `region`
constructor failed. Keeping the unsimplified expression when `simplify` raises is exact — same
number, not in normal form — and dropping the vertex instead would lose a real corner of the
region. `simplify` normalises a coordinate here; it does not decide whether the vertex exists.

### The diagnosis that was WRONG, and cost a run to refute

The first reading was "the candidate is INFINITE and `isreal(double(Inf))` is true, so a parallel
pair the `isAlways(dt == 0)` test could not prove parallel slips through" — the same shape as this
session's `slopeAtVertex` pole. Instrumented (`fprintf` on every non-finite candidate): **nothing
printed**. The candidate is finite, real, and correct; only `simplify` fails on it. Recorded
because that diagnosis is the plausible one and will suggest itself again — the guard it implies
(`isfinite` alongside `isreal`) is harmless but would have fixed nothing here.

## 2026-08-20 — the cross-face max HANDLES A RATIONAL PAIR now: clear the denominator where the cell certifies its sign

The blocker recorded below (2026-08-19, night, later) is closed. `region.maxArray` refused to
guess only for a POLYNOMIAL pair, because the refusal makes the caller split on `f1 = f2` and a
region's constraints must be polynomial -- `region.normalize1` raises
`symbolic:coeffs:NotAPolynomial` on a rational one. Every second-pass conjugate is rational, so on
the biconjugate path the refusal never applied: the vertex comparison decided, and on an
all-vertices-tied cell that means "f2 wins" for no better reason than operand order.

**The arithmetic, and it is one line.** With `f1 - f2 = N/D`, the sets `{f1 >= f2}` and
`{sign(D)*N >= 0}` are the SAME set wherever `D` has a constant strict sign. So on a cell that can
certify `D`'s sign, `sign(D)*N` is a POLYNOMIAL standing in for the rational difference, with the
same sign at every point of the cell, vertices included -- which is exactly what both the
orientation test in `splitmax3` and the constraints it returns need.

**Three pieces, all in `region.m`:**

* `signOnRegion(expr)` -- +1 / -1 / 0, where 0 is a REFUSAL in the same standing as
  `certifiesNonPositive`'s. A PRODUCT is decided factor by factor, which is the shape that
  actually arrives (`D = q1*q2` for two rational conjugates); each affine factor is decided by the
  LP over the region, STRICTLY both ways. Anything non-affine left after factoring is refused --
  `certifiesNonPositive` proves `h <= 0`, not `h < 0`, so it cannot supply the strict statement.
* `clearedDifference(f1,f2)` -- the polynomial, or a refusal with a reason.
* `maxArray` gains a SECOND sign attempt on the cleared polynomial before it refuses, and extends
  its refusal to any pair `clearedDifference` can represent; `splitmax3` splits on the cleared
  polynomial. A rational pair whose denominator the cell cannot certify falls through exactly as
  before.

**Why the second attempt is worth having on its own:** `isAlways` settles far fewer questions about
a rational difference than about a polynomial one -- it has to reason about the denominator's zero
set as well -- so clearing first also DECIDES cases that used to split, which costs cells.

**Refusing costs a cell; clearing by a sign-changing denominator costs a wrong answer.** Hence the
strictness, and hence `x + y - 1` on the unit simplex (zero on a whole facet) is a refusal, not a
sign. Four tests in `regionTest`, fast bucket, red before the change.

## ~~2026-08-19 (night, later) — the A.4 split is BLOCKED behind a KNOWN defect: maximumP cannot max two RATIONAL conjugates~~ (CLOSED 2026-08-20, see above)

Three genuine attempts at landing the A.4 cevian split in `plq_1p`. The failure moved each time
and the third attempt located it in something that was already written down elsewhere. Recording
under Blocked rather than attacking it a fourth time.

**Attempt 1 -- guessed anchor.** Passed the A.3 formula an anchor vertex and edge derived from the
PARENT triangle. Overshoot at `s = (0,0)` fixed (0.0429 -> 0); undershoot appeared at `s = (0,1)`
(0 against a true 0.125).

**Attempt 2 -- self-classifying sub-triangle.** Built the second sub-triangle as a piece and sent
it back through `convexEnvelope1`, so `domain` decides its own convex edges and far vertex instead
of the parent guessing. Same undershoot. Ruled the anchor out.

**Attempt 3 -- dispatch on the FACE's domain.** Found a real second defect on the way:
`conjugateFunction` reads the face's REGION from `obj.envelope(i).d` but took its convex-edge
count, slopes and far vertex from the PIECE's domain `obj.d`. For one envelope face those are the
same object, which is why it survived; for a genuine two-face envelope face 2 is dispatched on the
wrong geometry entirely. Fixed by rebuilding the face's `domain` from its own vertices, and only
when it differs. **Same undershoot.**

**What attempt 3's instrumentation shows, and it is conclusive:**

    parent nE = 2, nV = 0
    envelope faces = 2
      face 1  quadratic with sqrt(2) coefficients, over {P, A, R}
      face 2  y^2/(y - x/2 + 1/2)   -- RATIONAL, over {A, R, B}
    conjugate cells = 11,  conjfia = [1 6 12]      <- BOTH blocks populated, 5 cells and 6 cells

So the split is right, both faces are built, both are conjugated, and both blocks reach
`maximumConjugate`. What loses the answer is the MAX ACROSS them -- and that is a defect already
recorded in `biconjugateTest.m` (lines 246-251):

> `splitmax3` hands `f1 - f2` straight to `region()`, whose `normalize1` raises
> `symbolic:coeffs:NotAPolynomial` on a rational one. Every second-pass conjugate is rational, so
> the old vertex verdict stands: on an all-vertices-tied cell it picks f2 for no better reason
> than operand order.

Face 1's conjugate is quadratic and face 2's is rational, so their max is exactly the pair
`region.maxArray` cannot decide. It picks one, and the cells carrying the 0.125 at `s = (0,1)` are
the ones dropped.

**Therefore the A.4 fix is not a `plq_1p` change at all -- it is blocked on teaching the
cross-face max to handle a RATIONAL pair** (clear denominators where both are provably nonzero on
the cell, which is the fix `biconjugateTest` itself proposes). Until then, adding the split
replaces an OVERSHOOT with an UNDERSHOOT, and an undershoot is worse: a minorant is at least a
valid convex underestimator, while a conjugate below the sup is not a conjugate of anything.

**Reverted, and the fuller attempt is kept at `.claude/a4split_attempt.m.txt`** -- it contains the
working cevian geometry (`a4Split`, `twoEdgeQuad`, `triFromVertices`) and the face-domain dispatch
fix, both of which are correct and reusable the moment the rational max works.

**Left in place: one known-failing test**, `testMaxMultiRegion/testPCE2`, which fails on `main` and
did before today. It is now the executable statement of this gap.

## 2026-08-19 (night) — plq_1p's A.4 branch computes a MINORANT, not the envelope. Diagnosed exactly; the fix is started and not finished.

Chasing T6's failures found a defect that has nothing to do with T6: **`testPCE2` fails under BOTH
per-piece classes**, so it is not a migration regression but a wrong answer that no test could see
until this session's assertions existed.

**The defect.** `plq_1p.convexEnvelope1`'s `nCE == 2` branch applies [COAP] Appendix A.4's single
quadratic over the WHOLE triangle. `convEnvCPLQ.splitTwoConvexEdges` has carried the fix for that
since the 2026 sessions and its header states the reason: the single quadratic touches `u1*u2`
along both convex edges and is a valid convex MINORANT, but it is tight only over a sub-region
containing the two edges' common vertex. This branch never got the fix.

**Measured, on the triangle {(0,0),(1,0),(2,1)} with f = x*y:**

    unsplit A.4 envelope dips to -0.0429 inside the domain
    but f >= 0 there, so the true envelope's minimum is 0
    hence f*(0,0) = 0.0429 where the definition gives 0

The paper's own A.4.3 example has the same defect (`q1(0.474343,0) = -0.042780` for a true 0),
which is what `convEnvCPLQ`'s header already records. So this is the KNOWN A.4 gap, still present
on the symbolic per-piece path.

**The fix, and how far it got.** The cevian is forced, not chosen: the line through one convex
edge's far vertex with slope `-sqrt(mh*mw)`, met with the other convex edge -- the unique direction
along which the single quadratic and the A.3 formula agree in value AND gradient. Implemented
symbolically (so the foot stays an exact surd rather than a double), with `twoEdgeQuad`,
`oneEdgeRational` and `a4Split` factored out and the sub-triangles built through `domain` so they
get the same half-plane orientation rule as every other region.

**It half-works, which is why it is not committed.** The overshoot is GONE -- `f*(0,0)` is now 0 --
but a new UNDERSHOOT appears: `f*(0,1) = 0` where the sup over the triangle is 0.125, attained at
(0.5, 0.25) on the edge `y = x/2`. So the two faces do not yet cover what they should; the likely
causes are the choice of which far vertex anchors the A.3 face, or the vertex order of the second
sub-triangle. **Reverted to the committed baseline** rather than leaving a half-correct envelope
in place: the attempt is kept verbatim at `.claude/a4split_attempt.m.txt`.

**Status: one known-failing test, cause identified.** `testMaxMultiRegion/testPCE2` fails on `main`
and did so before today; what is new is that it now SAYS so instead of printing and passing.

## 2026-08-19 (night) — T6 REFUTED: plq_1p and plq_1piece are NOT interchangeable, so deleting the old class is not free

T6 was listed as the cheap win on the way to a sym-free CCA2: `plq_1piece.m` carries 75
symbolic-engine calls, is an older parallel implementation of the same per-piece API, and its only
live constructor calls are the 18 fixtures in `testMaxMultiRegion`. Move those onto `plq_1p` and
75 calls leave the surface without porting anything.

**Tried, and reverted.** The API gap is one method (`biconjugateP`, ~20 lines, ported cleanly).
The BEHAVIOUR gap is not.

> **CORRECTION, same evening.** The line below originally read "four of six short tests that pass
> under `plq_1piece` fail under `plq_1p`". That was measured on only two of the four. Baselined
> properly afterwards: `testConjugate`, `testFractional` and `testConvex` pass under `plq_1piece`
> and fail under `plq_1p` -- THREE regressions -- while **`testPCE2` fails under BOTH** and is a
> pre-existing defect the new assertions exposed, not a swap regression. See the entry above it.

With the fixtures swapped:

    testPCE0        pass
    testBiconjugate pass
    testConvex      FAIL  the convex envelope EXCEEDS f by 12.92 on the domain
    testConjugate   FAIL  f*(1,0) = 0.25, sampled sup over the domain is 2
    testPCE2        FAIL  f*(0,0) = 0.0429, sup over the domain is 0
    testFractional  FAIL  symbolic:coeffs:NotAPolynomial in plq_1p.conjugateFunction

Two of these are not "different but defensible" -- an envelope ABOVE `f` and a conjugate BELOW the
sup are both violations of the definition. Note `testConvex` uses `PRect3`, a four-vertex domain
put through `convexEnvelope` WITHOUT triangulating first; `plq_1p` appears to require the
per-triangle input its own `triangulate` produces, where `plq_1piece` tolerated the polygon.
`testFractional` fails differently again: it hand-builds a RATIONAL envelope face, and
`plq_1p.conjugateFunction` calls `coeffs` on it, which refuses a non-polynomial.

**Why this was invisible until today.** Those four tests had NO ASSERTIONS until this session --
they ran the pipeline, printed, and returned. Under the old regime the swap would have looked
free, because every one of them "passed". The numeric checks added a few hours earlier are what
turned an invisible behaviour change into four specific, located failures. That is the clearest
return on that work so far.

**What it means for the port.** T6 is not a deletion, it is a MIGRATION with real defects to fix
first -- either `plq_1p` gains the polygon and rational-face handling `plq_1piece` has, or those
fixtures move to the triangulated form and `testFractional` gets a polynomial envelope. Until
then `plq_1piece`'s 75 symbolic calls stay, and a sym-free CCA2 must either port them or land
this migration properly. Reverted to green; nothing is left half-swapped.

## ~~2026-08-19 (night) — T1/T2: the number type for a sym-free CCA2 is Q(sqrt(d)), and it is built and tested~~ (the FIELD is refuted 2026-08-20, see the T1 entry above; the three design rules stand)

T8 (below) established that rationals alone cannot carry this pipeline. T1 is therefore decided
and T2 is implemented: `exactQ`, values `a + b*sqrt(d)` with `a`, `b` rational (int64
numerator/denominator, lowest terms) and `d` a squarefree positive integer, or 0 for a purely
rational value. 19 unit tests, green (`exactQTest`, in the FAST bucket, 0.3 s).

**Three design decisions, each of which is a refusal:**

1. **ONE quadratic extension, and mixing two RAISES.** `sqrt(2) + sqrt(3)` is an error, not a
   promotion to a tower. Silently building `Q(sqrt2, sqrt3)` is precisely how an exact type turns
   back into a symbolic engine by accident; making it loud means the caller learns WHICH operation
   needs the tower instead of every operation paying for one. A rational combines with anything.
2. **int64 multiply RAISES on overflow.** MATLAB's int64 SATURATES silently, and a wrong answer
   that looks exact is the one outcome worse than a slow one. Cross-cancellation before every
   multiply (`mulRat` divides out the two cross-gcds first) is what keeps the guard rare -- a
   40-term telescoping product of moderate fractions stays exact, and is a test.
3. **`fromDouble` REFUSES what it cannot represent.** `exactQ(pi)` errors rather than rounding.
   Converting at the boundary is the caller's job; a type that quietly approximates at its
   constructor has no exactness to offer.

**Sign is the operation everything downstream rests on**, and it is computed with no floating
point at all: when `a` and `b` share a sign the answer is that sign, otherwise `a^2*bd^2` is
compared to `b^2*d*ad^2` in integers. Tested on successive convergents of `sqrt(2)`, which
alternate about it with shrinking margins -- `1393/985` is below by 3.7e-7 and `665857/470832` is
above by 1.6e-12, decided by integers one apart (1940449 vs 1940450).

**Two tests caught their own author**, which is the argument for writing them first: the
convergent's sign was asserted backwards (1393/985 is BELOW sqrt(2), not above), and
`exactQ(2,1,1,1,2)` was read as `2*sqrt(2)` when the constructor's order `(an,ad,bn,bd,d)` makes it
`2 + sqrt(2)`. Both are now stated in the tests so the next reader does not repeat them.

**The two refuted number types are pinned as tests**, not just as prose: the A.5 cevian foot
`5/2 - sqrt(5)/2` is representable and its sign is exact, and `4 - 2*sqrt(2)` -- the value that
arrived as two doubles one ULP apart and made a shared facet invisible to `merge` -- compares
equal to itself and negates exactly.

## 2026-08-19 (night) — T8 ANSWERED, NEGATIVELY: the A.5 surd is intrinsic, so a sym-free CCA2 needs quadratic surds

`ALGORITHM.md` left one question open on the way to removing the Symbolic Toolbox: A.5's split
foot is irrational (`5/2 - sqrt(5)/2`), and a rational route would make the whole pipeline
representable over the rationals. It recorded that a NUMERIC single bounded triangle "was measured
exact at 7 of 7 probe points and never introduces a surd", and asked whether that generalises.

**First, the surd's actual source is narrower than the file suggests.** A.4 (TWO convex edges)
envelopes to a single rank-1 PSD quadratic -- no split, no surd. Only A.5 (THREE convex edges)
splits along the cevian. So the question is about the 3-convex-edge case alone.

**And the tightness objection does NOT apply to the conjugate.** `ALGORITHM.md` argues the A.5
foot must not be moved because A.4's closed form is only tight on the right side of it -- but that
is a statement about STEP 1, the envelope. The conjugate satisfies `f* = max_k (f|T_k)*` for ANY
cover, with no tightness condition, so a rational split IS sound provided each sub-triangle is
conjugated DIRECTLY (`conjPieceCPLQ`) rather than through A.4's envelope formula. That is the
opening this probe tested.

**It is useful only if the split drops every sub-triangle to <= 1 convex edge**, which is what
`conjPieceCPLQ` accepts. For `f = u1*u2` an edge is CONVEX exactly when its direction `d` has
`d1*d2 > 0`. Two natural rational splits were measured (`.claude/t8RationalSplitProbe.m`):

    centroid split            all three sub-triangles REFUSED (conjPieceCPLQ:notImplemented)
    one axis-parallel cut     REFUSED on every fixture, including the A.4 and 1-edge controls

The centroid fails for a reason that generalises: a cevian's slope lies BETWEEN the slopes of the
two edges it separates, so on a triangle whose edges all have positive slope every cevian has
positive slope too and the count never drops. One axis-parallel cut gives each sub-triangle only
ONE affine edge, leaving two convex ones.

**What would work, and why it is not worth it.** Only a decomposition into axis-aligned RIGHT
triangles gets there: two affine legs plus at most one convex hypotenuse. Two things kill it:

1. **It is rational only when the bilinear FRAME is.** A general indefinite `Q` needs
   `M = bilinearFrame(Q)`, built from `sqrt(lam1/2)` and `sqrt(-lam2/2)` (convEnvCPLQ:181) -- so
   the surd reappears in the frame itself, before any split. The route would serve the `b*x*y +
   linear` family only.
2. **It multiplies the cell count feeding a QUADRATIC-cost Step 3.** Each triangle becomes up to
   six, and Step 3 is the measured bottleneck (`maximumP`, ~n^2 in cells). Trading an exact surd
   for six times the cells is the wrong direction.

**Conclusion, and it decides T1: exact arithmetic over Q(sqrt(d)) is REQUIRED.** A sym-free CCA2
cannot be built on rationals alone; the number type has to carry one quadratic extension. Doubles
remain refused for the reason already recorded -- one ULP made a shared facet invisible to `merge`
and Step 3's cell count grew without bound -- and rational snapping remains refused for the reason
recorded under attempt 3 (vertex denominators do not bound downstream ones; 1e5 became 1e25).

## 2026-08-19 (evening) — the crash tests became tests, and the one red turned out to be an unguarded POLE in slopeAtVertex

Four things, and the fourth is the reason the other three were worth doing.

### The two big suites had 32 tests and ZERO assertions

`testcPLQ` (8) and `testMaxMultiRegion` (24) each ran a prefix of
`triangulate -> convexEnvelope -> conjugate -> maximum -> biconjugateF`, PRINTED the answer, and
returned. They passed whenever nothing threw. That was 6333 s of a 7207 s bucket -- 90% of it --
to establish that the functions return. For contrast, `conjCPLQTest` has 74 assertions and costs
238 s.

Every one now asserts against a DEFINITION rather than a golden value (`plqCheck`), so nothing
here needs re-pinning when a representation changes:

    convex envelope   co f <= f on the domain, and co f = f at the vertices
    conjugate         f*(s) = sup_{x in D} <s,x> - f(x), against a numeric sup over the domain
    biconjugate       f** <= f on the domain, and f** convex along sampled segments

The sup sampler is a LOWER bound (boundary at a fine step, every vertex, plus an interior grid),
so `f* < sup_sampled` is a definite defect while `f* > sup_sampled` within tolerance is expected;
both directions are checked with the tolerance each deserves.

**Two bugs in the CHECKS themselves, found by running them.** (a) Step 1's envelope is
quadratic-over-LINEAR and its denominator vanishes at some domain VERTICES -- exactly where an
underestimator check most wants to sample -- so `evalFunctionNDomain` raised
`symbolic:kernel:DivisionByZero` and the test died on the fixture's own geometry;
`plqCheck.safeEval` now treats a pole as "no value here". (b) `plq.conjugate` leaves each piece's
result in `.conjugates` (one cell per envelope FACE), while `.maxConjugate` -- what
`maximumConjugate` writes -- is the object that equals the sup; comparing the wrong one reported
every point "uncovered".

### The heavy tests are split by STAGE, with the intermediate cached

`plqStage.get(fixture, stage, fcn)` caches each stage keyed by (fixture, stage) and invalidates on
the mtime of any repository `.m` file. So a cold run costs what it always did, an edit recomputes,
and a re-run after a failure re-does only the broken stage onward. `testRectBiconj` was ONE test
of 3198 s that, when it began failing, produced nothing but an exception after 53 minutes; it is
now four tests that each say which stage broke.

The cache is derived data under `.claude/stagecache/`, is not tracked, and deleting it must change
nothing but runtime. **No test may assert against a cached value the same run did not verify** --
otherwise a stale cache becomes a passing test, which is worse than no cache.

### A `--verylong` bucket

`--slow` is now the five suites carrying the real assertions and finishes in minutes;
`--verylong` is `testcPLQ` + `testMaxMultiRegion`, the pipeline endurance run. Keeping them
together meant the bucket you run after touching `region.m` cost two hours, so it did not get run
-- which is precisely how the stale `conjCPLQTest` expectation went unnoticed for a day.

### THE RED: a pole, in a routine nobody has touched since the Phase 1 integration

`testcPLQ/testRectBiconj` passed in the morning's serial run and failed in the evening's. Stack:

    sym/subs                     Division by zero
    region/slopeAtVertex   1374  abs(subs(drx2.f, vars, pt)) < 1.0d-6
    region/simplifyUnboundedRegion
    functionNDomain/maximumP
    functionNDomain.maxOfList
    plq/biconjugateF

`git log -S slopeAtVertex` puts its last change at `92c9c96`, the Phase 1 integration. **So none of
the day's four rewrites broke it; they changed which regions reach it.**

The guard asks "is `dg/dy ~ 0` at pt", so that a vertical tangent gets slope `inf` instead of a
division by zero -- but it asked by EVALUATING `dg/dy` at pt, and on this path the constraint can
be RATIONAL (every second-pass conjugate is), so `dg/dy` has a denominator of its OWN which can
vanish at pt. **A pole is the opposite of what the guard is looking for**: `dg/dy` blowing up means
the tangent is HORIZONTAL. So an unevaluable guard must fall through to the ratio -- and the ratio
`-dg/dx / dg/dy` is exactly where the two denominators CANCEL. Fixed that way, with
`simplifyFraction` on the ratio and `intmax` (this routine's existing "vertical/undecided" marker)
kept for the case where nothing cancels.

## 2026-08-19 (later) — the symbolic-removal list is DONE for region.m: isFeasible was weak as well as slow, and both rootsIn fallbacks are dead

Items 5d/5e/5f/5g. Live `solve()` in `region.m` is now **3**, two of which are `region.rootsIn`'s
own fallback and one of which (`isFeasible`) is kept deliberately, for conic pairs only.

### 5f. The sweep outside region.m found ONE genuine instance, and one already-correct one

`functionNDomain.m` has six `solve()` sites, `plq_1p.m` one, and the pattern that produced the
`removeTangent` defect appears in exactly one of them:

* **`conjugateOfPiecePoly`'s isQuad chord rewrite** (line ~1296) already asked the right question
  -- "is root one feasible, else take the other" -- but asked it in a way that breaks twice:
  `my2(2)` is indexed UNCONDITIONALLY when root one fails, so a conic with a single root at the
  chord midpoint (a tangency) raises `Index exceeds array bounds`; and the second root is never
  itself tested for feasibility OR realness, so a complex root can become the probe and the
  orientation test then reads it. Now `probeOnConstraint`, which is that question stated once.
* **`pointBetweenOnCurve`** (line ~2210) was ALREADY property-based -- it takes
  `min(abs(sol))`, the root nearest the chord -- and needs nothing. Worth recording as the
  pattern to copy.
* The remaining sites solve SYSTEMS (`solve([eq1,eq2], vars)` for a gradient equation,
  `solve(ineq1,[x,y])` for a subdifferential row) rather than picking one root of one polynomial
  by position. Different mechanism; out of scope, and left alone.

### 5e. Both `rootsIn` fallbacks are DEAD on every path four suites exercise

Instrumented the five branches and ran `regionTest`, `functionNDomainTest`, `maxQuaParTest`,
`convEnvCPLQTest` (71 tests, all green):

    not a polynomial in v   0
    constant in v           2
    affine                  5
    QUADRATIC              23
    degree > 2              0

So every call takes the closed form; **the two `solve()` fallbacks were never entered**. They stay
anyway -- they cost nothing when unreached and they are the honest answer for a shape this file
does not currently produce -- but nothing is relying on them. Scope: these four suites. The slow
bucket was not instrumented, so a slow-only path could in principle differ.

### 5d. `isFeasible` was pairwise, and pairwise emptiness is a WEAK test

It asked `solve(g_i <= 0, g_j <= 0)` for every PAIR -- O(n^2) calls into the engine -- and pairwise
emptiness cannot see n-way infeasibility: `x >= 1, y >= 1, x + y <= 1` is infeasible with every
PAIR of its three constraints perfectly satisfiable. The affine constraints are now decided
TOGETHER by one LP (`region.maxLinear`, whose `st = -1` IS the certificate), and `solve` is kept
only for pairs a CONIC takes part in, which the linear relaxation drops. Cheaper AND strictly
stronger on the affine part.

Two things preserved deliberately: the `g_i == -g_j` refusal stays for EVERY pair (the LP would
report `g == 0` feasible, and changing that verdict is not what this change is for), and an EMPTY
region object now returns false through an explicit guard rather than erroring on
`size(obj.ineqs,2)` -- `region()` itself returns `region.empty()` for an input it can already see
is infeasible, so that path is reachable.

Nine cases, all correct (`.claude/isFeasibleCases.m`), including the three-way case above, a
tangency, a conic disjoint from a line, and an unbounded wedge.

### 5g. The `solveForY` dead end is struck

That entry's closing recommendation -- "make the callers order-INDEPENDENT" -- has now been
carried out at four sites and is `region.probeOnConstraint`. It is marked struck with the recipe,
so it reads as settled rather than as an open invitation to retry the naive substitution. One
refinement experience added: the property is not always feasibility. A probe that reads a curve's
SIGN near a vertex needs the root on the same BRANCH, which is the nearest one.

## 2026-08-19 — three more symbolic-removal sites: isconvex rewritten, removeTangent's probe was on the WRONG BRANCH, getNormalConeEdgeQ0 was dead

Items 5a/5b/5c of the symbolic-removal list, all measured before being changed. Live `solve()` in
`region.m`: **10 → 3**, and two of the three are `region.rootsIn`'s own fallback, so ONE genuine
call remains (`isFeasible`, which solves a pair of INEQUALITIES and is a different mechanism).

### 5a. `isconvex` — same rewrite as the normal cones, and it refuses exactly where A3 says it should

Four `solve()` calls, same shape as the ones fixed on 2026-08-18: substitute one abscissa, take
root ONE, and if that point is infeasible flip to the other ABSCISSA — never the other root. Now
`region.probeOnConstraint`, first FEASIBLE root.

**Differential over the 22 captured merge operands, against an independent midpoint test** (walks
each region's other edge at the shared vertex at three step sizes, choosing no root by position):
34 probes, **30 agree, 0 differ**, 4 undecided-by-the-oracle.

**On those 4 the old code said TRUE and the new code says FALSE, and the new answer is the better
one** — the oracle finds no point of the region on that edge at any step size, so the old `true`
came from a probe point that had failed its own feasibility test and was used anyway. It costs
nothing: they are `mg_0003/0004/0009/0013`, which are FOUR OF THE FIVE REFUSALS A3 PROVED CORRECT.
The pair is refused either way; it is now refused earlier, without the LP.

Refusing is the safe direction here by construction: `isconvex` is a LOCAL necessary condition,
its caller reads false as "no shared facet" (`mergeTally 'noSharedFacet'`), and `unionIsExact`
makes the real decision afterwards. A false costs compactness, never correctness.

### 5b. `removeTangent` — the probe was landing on the FAR BRANCH of the conic, 22 times in 34

This one is not a reproducibility fix, it is a DEFECT. The probe exists to read the conic's sign
just inside the region near a vertex, by forming the midpoint of the vertex and a nearby point of
the conic. At `sx = px + 0.1` a conic generally has TWO points, one near `py` and one far away,
and the code took the first REAL root. Measured over every conic-vertex pair in the captured
operands and the three curved fixtures:

    vertex/conic probes 34   single-root 2   FIRST-REAL-ROOT IS ON THE OTHER BRANCH: 22

and not marginally — the gaps are 1.26, 2.31, 2.45, 3.20, 5.99 units for a probe meant to sit
0.1 from the vertex. A midpoint built from a point 6 units away says nothing about the
neighbourhood it is supposed to describe.

The fix is the branch-continuous choice, which is also order-independent: the root NEAREST `py`.

### 5c. `getNormalConeEdgeQ0` — dead, so it was removed rather than rewritten

Its two `solve()` calls encode a TANGENCY condition (the line of a given slope meeting a conic in
a double root), which is the discriminant and has a closed form. But `grep` over the executed tree
finds exactly one occurrence of the name: its own definition. It came in verbatim with the Phase 1
cPLQ integration (`92c9c96`) and was never called; the archival original remains in `cPLQ/`, which
is never executed. Rewriting code with no caller adds a path no test can reach, so it was deleted
with a note in its place. Its live siblings (`getNormalConeEdge`, `...Q`, `...Q3`, `...QE`) are
untouched.

## 2026-08-19 — A3 ANSWERED: merge never over-claims, and every unresolved refusal is a CURVED certificate

`SCIP_READINESS.md` A3 asked whether `unionIsExact`'s ~43 refusals per fold are CORRECT or a
defect, and said not to optimise the gate before knowing which. Measured, and the answer splits
cleanly by REASON.

**The oracle is decisive, not a sampler, and it comes from `unionIsExact`'s own algebra.** merge
returns `M = A' ∩ B'` (each region with the shared facet deleted) and `M ⊆ A ∪ B` ALWAYS. So the
only way merge can be wrong is to LOSE a point: a point of `B` violating some constraint of `A'`,
or the mirror. Finding one PROVES the refusal correct; that is a witness, not evidence. Every
witness the numeric search reports is re-verified exactly against the region objects
(`.claude/a3score.m`). Note this also means the question "is `A ∪ B` convex" is a detour --
losing a point is the whole criterion.

**22 calls captured over three folds of the A.4/A.5 quadrilateral** (temporary probe in `merge`,
reverted):

    accepted, no lost point                      12     <- sound, every one
    refused, LOST POINT FOUND (proven correct)    5
    refused, none found in 2e6 samples            5

**DEFECTS: ZERO.** No accepted merge loses a point. That is the result that matters most, and it
says the cell counts are not hiding an over-claim.

**Every LP-decided refusal is CORRECT.** All five proven-correct refusals are
`exactAnotInB` / `exactBnotInA` / one `curved_positiveSomewhereOnP`, each with an exact witness:

    mg_0003, mg_0004  exactBnotInA                witness (-1.542,  5.240)  in B, violates A'
    mg_0009           curved_positiveSomewhereOnP witness (-2.067,  6.091)  in B, violates A'
    mg_0013           exactAnotInB                witness (-1.051,  4.442)  in A, violates B'
    mg_0018           exactBnotInA                witness (-0.722,  2.988)  in B, violates A'

**Every unresolved refusal is a CURVED certificate, and they are the whole conservatism budget.**
The five with no witness are all `curved_positiveSomewhereOnP` (four) or `curved_positiveAtAVertex`
(one) -- i.e. `region.certifiesNonPositive` could not certify a conic constraint. Re-probed with
2e6 candidates over a box padded to 20x each region's extent: still no lost point. Each of those
pairs carries exactly one or two curved constraints per operand.

**What this decides for Phase C.** The gate worth sharpening is `certifiesNonPositive`, NOT the
LP and NOT the refusal policy -- the LP's verdicts are provably right, and the refusal policy is
what keeps merge sound. The upper bound on what a perfect curved certificate could buy is those
five of 22 calls; whether that is worth doing is a Phase C cost question, not a correctness one.

**Honest limit:** "no lost point in 2e6 samples" is not a proof of convexity, only strong evidence
that those five merges were available. The five witnesses in the other direction ARE proofs.

## 2026-08-18 (evening, sixth) — the first symbolic-removal site LANDED: probes pick a root by feasibility, not by position

The third entry said the remaining `solve()` sites are a CALLER REWRITE, not a substitution,
because each computes a probe point on a curve and then picks one root BY POSITION, and `solve()`'s
ordering is an unreproducible MuPAD convention. `getNormalConeVertexQ` is now rewritten that way,
which is item 1 of that programme and the one the specification (fifth entry) was blocking.

**What replaced what.** Eight `py = solve(ey, y)` calls, each followed by "take root one; if the
point is infeasible try the other ABSCISSA" -- never the other ROOT -- became four calls to a new
`region.probeOnConstraint(cIdx, xs)`: *the first (x,y) with y a real root of constraint cIdx at x
that the region actually contains*, over the candidate abscissae in order. Roots come from
`region.rootsIn`, the quadratic formula in closed form, with `solve()` kept as the fallback for
anything not polynomial of degree <= 2 or with a leading coefficient it cannot show nonzero. Live
`solve()` calls in `region.m`: **16 -> 10**, and two of the ten are that fallback.

**MEASURED against the specification, which is the whole reason this could be attempted.** With
the edge list -- the contract of the fifth entry -- all three curved fixtures still score
**0 of 72 wrong directions at every vertex**. Without it, on the slot fallback applied to bounded
fixtures where it does not belong, the cones move CLOSER to the definition:

    piece 9    32 -> 29 wrong directions
    half lens  43 -> 29
    parabola    5 ->  5

**Three pinned values moved, all orientations**, and that is the mechanism working: piece 9's
(3,2) row and both of the half-lens's first column flip sign, because a different -- feasible --
probe now decides which side the cone is on. `regionTest.normalConesOnCurvedEdgesAreUnchanged` is
re-pinned with the reason written into it; the CONTRACT test is unchanged and still green.

**One defect went with the rewrite.** In the `cNext` block the second attempt read
`obj.ineqs(cj)` where the first read `obj.ineqs(cNext)` -- so a vertex whose probe failed on the
left was re-probed against the OTHER edge's constraint. Same copy-paste family as the two index
inversions already corrected in this routine. Asking one constraint on both sides is what
`probeOnConstraint` does by construction.

**Not claimed: a speed-up.** The point is reproducibility -- the answer no longer depends on
`solve()`'s root order -- and the closed form is a by-product. fast 217 / 0, regionTest 18 / 0.

## 2026-08-18 (evening, fifth) — getNormalConeVertexQ's specification, established. It IS the normal cone, and given the edge list it is exact.

The previous entry said the routine could not be replaced because there was no statement of what a
replacement must satisfy, and that the way to get one was to read the CONSUMER. Done, and the
answer is short. `functionNDomain.getSubdiffVertexT1` reads a row of `NC` in exactly three ways:

  * it takes `coef = getLinearCoeffs(row)` and branches on the SIGNS of the two coefficients;
  * it takes the slope as `m = d(row)/d(s1)` -- which is the row's own slope only when the
    coefficient of `s2` is `+-1`, so **normalisation is part of the contract**, not cosmetics;
  * it THROWS THE CONSTANT AWAY and re-anchors the line at `grad f(v_j)`.

The region it then builds is `region(subdV(j,:), ...)`, and region's convention is `expr <= 0`
feasible (`ptFeasible`). So the only observable object is the LINEAR PART under `<= 0`, and the
identity the vertex branch of the conjugate rests on is
`subdiff f(v_j) = grad f(v_j) + N_D(v_j)`. Hence:

> **Row j's linear parts, read as `<= 0`, must cut out the normal cone of the region at vertex j;
> the coefficient of `s2` must be `+-1` or `0` (and then `s1`'s must be `+-1`); the constant is
> free.**

**GIVEN THE EDGE LIST, THE COMMITTED IMPLEMENTATION SATISFIES THAT EXACTLY.** Same oracle as
yesterday, built from the definition (`u` in the cone at `v` iff `v` maximises `<u,.>` over the
region near `v`; local form, 72 directions, sampled to 5% of the vertex), now applied with the
`<= 0` reading and with `eIdx` supplied:

    piece 9   eIdx [2 1 3]   disagree 0   (per-vertex 0 0 0)     -- vertex 3 has THREE active constraints
    half lens eIdx [1 3]     disagree 0   (per-vertex 0 0)       -- vertex 1 is a CUSP
    parabola  eIdx [2 1]     disagree 0   (per-vertex 0 0)

The cusp and the concave-conic vertex are the two cases that killed the gradient rewrite. This
routine gets both right, which is the reason it exists: it builds the cone from the constraint's
own tangent, so a concave conic and a cusp are ordinary.

**What the earlier "4-30 of 72 disagree" was measuring.** The eIdx-less SLOT fallback, on fixtures
it does not apply to. Without `eIdx` the routine pairs vertex j with constraints `j` and `j+1`,
the layout of an UNBOUNDED region (slot 1 reserved for the ray, nv+1 slots for nv vertices). All
three fixtures are BOUNDED, so it reads the wrong constraint at some vertices -- the parabola's
vertex `(1,1)` is bounded by the parabola and `y = 1`, and the fallback used `y = 1` and `x = -1`,
giving a cone 5 directions too small. With the right pair the same code is exact. Totals without
`eIdx`: 32, 43 and 5 wrong directions.

    parabola v1 = (-1,1)  3 active constraints  fallback slots [1 2]  disagree 0
    parabola v2 = ( 1,1)  active [1 2]          fallback slots [2 3]  disagree 5 (all "too small")

**Also measured: with a lens, the edge list must carry the TRAVERSAL ORDER, not just the right
pair.** Both admissible lists for the half-lens use constraints {1,3} and both vertices see both,
yet `eIdx = [3 1]` is wrong at vertex 1 on 36 of 72 directions while `[1 3]` is exact. The routine
deduces each row's orientation from a probe on the ARRIVING edge, so which of the two is arriving
changes the answer. Mathematically it should not; that is the probe fragility of the previous
entry, now with a reproducible case.

**Kept:** `regionTest.vertexConesMatchTheDefinition` -- the specification as an executable test,
72 directions per vertex over the three fixtures, ~30 s, and GREEN against committed code (unlike
the property test reverted yesterday, which was red because it tested the fallback). regionTest is
now 17 / 0 in 44 s.

**AND THE FALLBACK IS SOUND ON ITS OWN LAYOUT -- settled the same evening, so it is not open.**
Two suites (`functionNDomainTest`, `unboundedFaceTest`, both green) captured ZERO eIdx-less calls,
so the pipeline was not going to answer it; the layout was rebuilt instead. `getEdgeNosInf` puts
the constraint carrying the ray at vertex 1 in slot 1 and edge j (vertex j to j+1) in slot j+1, so
at vertex j the arriving element is slot j and the leaving one slot j+1 -- exactly the pair the
fallback takes. Built that way (`removeInfV`, `poly2orderUnbounded`, the scatter) on
`{y >= x^2, -2 <= x <= 2}`: **0 of 72 directions wrong at both vertices**. On a BOUNDED region edge
j is slot j and the same pair is off by one, which is the whole of the earlier disagreement. Kept
as `regionTest.theSlotFallbackIsRightOnTheUnboundedLayout`; regionTest is 18 / 0 in 45 s.

**A REAL LIMIT, found by sharpening the oracle, and it is benign.** At 1 degree of boundary
sampling the oracle called a direction 0.96 degrees OUTSIDE a cone "inside" -- the nearest sampled
boundary direction landed exactly on the degenerate value. At 0.2 degrees (afforded by vectorising
the feasibility test) that artifact is gone, and one genuine difference surfaces: **where the
region is on the CONCAVE side of a conic the exact normal cone is not CLOSED.** At piece 9's vertex
2, locally `{x <= y^2/4}`, the tangent perpendicular `(1,0)` has region points strictly ahead of it
(`x = y^2/8 > 0`) while every direction just inside the cone does not. The routine returns the
CLOSED cone. The two differ by one ray, which the cell decomposition does not distinguish --
adjacent conjugate cells share their boundaries anyway -- so the test excludes directions ON the
cone's boundary from its verdict, and says why.

**What is left is a narrow one:** a BOUNDED region for which `edgeIndexList` refuses (returns
`ok = false`) still reaches the fallback, and there the pair is off by one. Nothing has been seen
to produce one; the oracle to check it is a reusable static, `regionTest.coneVsDefinition`.

## ~~2026-08-18 (evening, fourth) — getNormalConeVertexQ does NOT compute a normal cone. Its specification is unknown.~~

> **OVERTURNED the same evening** by the entry above: the specification IS the normal cone, and
> given `eIdx` this routine computes it exactly. What the measurement below caught is the
> eIdx-less slot fallback applied to bounded fixtures. The gradient rewrite's two failure modes
> stand -- they are why the tangent-based construction is the right one.

The gradient rewrite was built and measured, and the measurement is the useful part.

**The rewrite.** Replace slope-and-probe with `grad g`: the region is `{g_i <= 0}`, so at vertex v
the outward normal of an ACTIVE constraint is `+grad g_i(v)` -- direction and sign in one
evaluation. The active set is found by `g_i(v) = 0`, so no edge-slot convention is needed. 21.7 KB
of code became 5.3 KB, and `solve()` calls in `region.m` fell 19 -> 11.

**It is not correct, and neither is what it replaced.** Tested against the definition of a normal
cone, with an ORACLE built from that definition rather than from any implementation:

  * `u` in the normal cone at `v`  <=>  `v` maximises `<u,.>` over the region.
  * GLOBAL form for convex sets; LOCAL form (over the region intersected with a small ball) for
    non-convex ones.

Two mathematical failure modes, both real:

  * **cone too BIG** -- piece 9's vertex 2. The region is locally `{x <= y^2/4}`, the CONCAVE side
    of the conic (Hessian eigenvalues 0 and -10), so the gradient half-plane does not contain it:
    points with `x = y^2/8 > 0` are inside, and direction `(1,0)` is therefore NOT in the true
    cone though it is a gradient.
  * **cone too SMALL** -- the half-lens vertex 1. The parabola's gradient `(-4,0)` coincides with
    `-x`'s, but the region is a CUSP there, whose normal cone is strictly larger than the cone the
    gradients generate.

So the gradient cone equals the normal cone only under convexity plus a constraint qualification,
and this codebase's pieces satisfy neither.

**AND THE ORIGINAL FAILS THE SAME TEST.** Measured: the committed implementation disagrees with
the global definition on 4-30 of 72 directions per vertex, and also fails the local one. So
`getNormalConeVertexQ` is not computing "the normal cone of this region" in either sense.

**What that means, and it is the actionable part.** Either it computes something else BY DESIGN --
its own header describes an edge-SLOT convention ("the cone at vertex j is built from ineqs(j) and
ineqs(j+1)") tied to how `conjugateOfPiecePoly` consumes it, which is a different object from the
region's normal cone -- or it carries a latent defect that the pipeline compensates for elsewhere.
**Until that is settled the routine cannot be replaced**, because there is no statement of what a
replacement must satisfy. Establishing it means reading the CONSUMER (`getSubdiffVertexT1`, which
re-anchors these rows at `grad f(v)`), not the routine.

**Kept:** the characterization test, which pins the current output in ~13 s. **Reverted:** the
gradient implementation and the property test, the latter because it fails against the committed
code and would therefore be red from the start.

## 2026-08-18 (evening, third) — EVERY solve()-for-a-probe-point site is ORDER-DEPENDENT. Measured.

A differential oracle (`scratchpad/diff_roots.m`, seconds to run) compared `solve()` against the
closed-form quadratic on seven conics x six probe abscissae, asking the question the CALLERS
actually ask -- "the first REAL root":

    cases 42, disagree 16, and every disagreement is a case where BOTH roots are real.

So the picks differ whenever the probe line crosses the conic twice, which for these shapes is
the common case, not the exception. `removeTangent` loops for the first REAL root rather than
taking `solve()`'s first, and even that is not enough: with two real roots, position still decides.

**This generalises to the whole substitution programme.** Every remaining `solve()` in `region.m`
-- the normal cones, `isconvex`, `removeTangent`, `pointBetweenOnCurve` -- computes a PROBE POINT
on a curve and then picks one root by POSITION. `solve()`'s ordering is an internal MuPAD
convention (measured: `(x+y)^2-4x` at `x=1` ascending, `x^2+y^2-4` at `x=1` descending), so no
closed form reproduces it, and substituting changes which probe is used -- hence which orientation
is deduced, hence the answer.

**The work is therefore not a substitution, it is a caller rewrite.** Each site has to pick its
root by a PROPERTY -- first feasible, first satisfying the orientation test -- instead of by
position. That is strictly better than today (it also fixes sites that give up while a good probe
sits on the root they never examined), but it changes behaviour and so needs a property oracle per
site.

**And one such oracle exposed a second problem.** `getNormalConeVertexQ` called WITHOUT `eIdx`
falls back to "constraint j bounds edge j" by slot -- the convention `edgeIndexList` exists to
replace. Tested standalone against the definition (s is in the normal cone at v iff v maximises
<s,.> over the region), it disagrees on 4-30 of 72 directions per vertex. So a fast test for that
routine MUST supply the edge-index list the pipeline supplies; testing it bare pins a
configuration the pipeline does not use.

**Order of work for the next attempt**, which is now clear:
1. property oracle for the normal cone WITH `eIdx`, mirroring `conjugateOfPiecePoly`'s call;
2. rewrite the probe selection to be property-based, verified against it;
3. only then substitute closed forms, which at that point cannot change the answer.

## 2026-08-18 (evening, second attempt) — the closed-form normal cones are SLOWER. Measured, not guessed.

The first attempt was reverted because `solve()`'s root order is not reproducible. This attempt
removed the dependence instead of trying to match it: `region.probeOnCurve` tries BOTH sides of
the vertex and BOTH roots, keeping the first feasible probe, so ordering stops mattering.

**It works and it is slower.** Five calls of `getNormalConeVertexQ` on piece 9:

        solve()                 2.89 s
        closed form + probe     3.77 s      +30%

The reason is arithmetic, not subtlety: `conicCoefsSym` costs SIX substitutions per call, and the
probe then runs up to four `ptFeasible` tests, where the original did one or two `solve` calls on
a small expression. The closed form only wins when the thing replaced is expensive relative to the
setup -- true for `getVertices`, where one extraction served O(n^2) pairs, and false here, where
the extraction serves a single vertex.

**Also measured on the way:** converting only ONE of the four probe blocks left the cones
unchanged; converting the other three MOVED them. So the probe choice is genuinely load-bearing,
and any future attempt has to justify the new choice rather than assume equivalence.

**What was kept:** `regionTest/normalConesOnCurvedEdgesAreUnchanged`, a characterization test that
pins the curved-edge cones for three regions in about 13 seconds. It caught both failures in this
attempt -- a static/instance placement bug, and the probe-order behaviour change -- where the
suites that cover these paths (`conjCPLQTest`, `testMaxMultiRegion`) are in the ~100-minute slow
bucket. That test is the reusable result of the exercise; the optimisation is not.

**Rule this establishes for the remaining `solve()` sites:** extraction cost is amortised over the
number of uses. Convert where one extraction serves many operations; leave it where it serves one.

## ~~2026-08-18 (evening) — REVERTED: closed-form NORMAL CONES. solve()'s root ORDER is not reproducible.~~

> **SETTLED 2026-08-18/19 — this entry's own recommendation was carried out, and the recipe is
> now `region.probeOnConstraint`.** The diagnosis below stands and is why the naive substitution
> was right to revert; what is stale is reading it as an OPEN item. Four sites have since been
> rewritten the way the last bullet asks for -- pick the root by a PROPERTY, not by position:
>
>   * `getNormalConeVertexQ` (2026-08-18): eight `solve()` calls -> four `probeOnConstraint`
>     calls, first FEASIBLE root. Cones still exact against the definition given `eIdx`.
>   * `isconvex` (2026-08-19): four calls, same recipe; agrees with an independent oracle on
>     every decidable probe.
>   * `removeTangent` (2026-08-19): the property there is NOT feasibility but BRANCH CONTINUITY
>     -- the root nearest the vertex's own ordinate. First-real-root was on the wrong branch in
>     **22 of 34** probes, so this one was a defect, not just an ordering risk.
>   * `getNormalConeEdgeQ0` (2026-08-19): deleted instead of rewritten -- it had no caller.
>
> Live `solve()` in `region.m`: **16 -> 3**, two of the three being `region.rootsIn`'s own
> fallback. The "do not retry without the slow bucket" warning was honoured: every one of these
> went through it. Details in the 2026-08-18 (evening, sixth) and 2026-08-19 entries.

## 2026-08-18 (evening) — REVERTED: closed-form NORMAL CONES. solve()'s root ORDER is not reproducible.

- **Tried:** `region.solveForY`, a closed-form replacement for the ten
  `ey = subs(ineq, x, px); py = solve(ey, y)` pairs in `getNormalConeVertexQ` and
  `getNormalConeEdgeQ0`. Fixing `x = px` turns a conic into a quadratic in `y`, so the roots are
  the quadratic formula. Eight sites converted cleanly and the algebra is right.
- **Why it was reverted, and it is NOT the mathematics:** the callers take **`py(1)`** as a probe
  point, test it with `ptFeasible`, and if it fails flip `px` to the other side -- they never try
  the OTHER root. So which root comes first changes the probe point, and can change the cone.
- **`solve()`'s ordering cannot be reproduced.** Measured against it on six conics:

        (x+y)^2 - 4x  at x = 1   solve gives  [-3 ; 1]                   ASCENDING
        x^2 + y^2 - 4 at x = 1   solve gives  [3^(1/2) ; -3^(1/2)]       DESCENDING

  Sorting ascending fixes the first and breaks the second. It is an internal MuPAD convention,
  not a rule to reimplement.
- **What would make it safe**, and is the right fix when this is next picked up: make the callers
  order-INDEPENDENT -- try both roots and keep whichever is feasible, instead of taking `py(1)`
  and flipping `px`. That is strictly more robust than today's behaviour, but it CHANGES
  behaviour, so it needs the slow bucket, where `conjCPLQTest` and `testMaxMultiRegion` exercise
  the curved normal-cone paths. **DONE — that is `region.probeOnConstraint`; see the note above
  this entry. The one refinement experience added: the property is not always feasibility. A
  probe that reads a curve's SIGN near a vertex needs the root on the same BRANCH, which is the
  nearest one, not the first feasible one.**
- **Do not retry without the slow bucket.** Fast and normal do not cover those paths, so they
  would have passed either way -- which is exactly how a silent geometry regression gets in.

## 2026-08-18 (measurement) — `noSharedFacet` is MOSTLY HONEST; the open question moved to `unionIsExact`

`.claude/step3adjacency.m` with `CCA2_ADJ_FOLDS=3`, on the A.4/A.5 quadrilateral -- 137
same-function pairs over 36 cells carrying 7 distinct functions:

    merge sees | shared plane | how they meet | pairs
    no         | no           | do not touch  |  42   correctly refused
    no         | no           | a point       |   2   correctly refused
    yes        | yes          | do not touch  |  14   correctly refused
    yes        | yes          | a point       |  16   correctly refused
    NO         | no           | A SEGMENT     |   7   meet along a CONIC merge cannot see
    NO         | yes          | A SEGMENT     |   4   genuine facet-detection misses
    yes        | no           | a segment     |   6   seen via a quadratic facet
    yes        | yes          | a segment     |  46   reach unionIsExact

**54% of the pairs are correctly not merged** -- they do not touch, or touch at a single point,
where the union of two cells is not convex and merging would be wrong. So the raw
`noSharedFacet` count (346 at fold 3, 819 at fold 5) is NOT a defect count, which is what the
2026-08-17 entry suspected and this settles.

**Only 11 pairs are genuine detection misses**, and they split in two: 4 share an affine
hyperplane the symbolic test does not match, and 7 meet along a CONIC that neither the linear
facet search nor the quadratic branch identifies.

**The open question has moved.** 52 pairs reach `unionIsExact` and about 9 merges happen at that
fold, so roughly 43 are refused by the sound gate. Whether those refusals are RIGHT is not
established here -- two cells can share a facet, touch along a segment, and still have a
non-convex union. **Measure that before optimising it**: take a handful of the 46 and check
directly whether `A u B` is convex. Chasing `unionIsExact` without that risks the same mistake
this session made twice, of reading a correct refusal as a defect.

## 2026-08-18 (decision) — `RatPar.V` stays `{mustBeNumeric}`. The mesh route is a known-inexact FALLBACK.

**The question.** Mesh vertices are declared `V (:,2){mustBeNumeric}` on `RatPar`, for the whole
RatPar / RatPol / QuaPol / QuaPar lattice, so `convEnvCPLQ` returns `sqrt(2)` as `1.4142` and two
cells can carry two roundings of one number. That is the last double leak, and it costs 5 of 31
same-function pairs a real shared facet on the tri case.

**DECIDED: do not change it.** Three reasons, in order of weight.

1. **It costs cells and time, not correctness.** Every consequence measured is a merge refused and
   a cell not collapsed. No wrong value has been traced to it.
2. **The path that matters no longer goes through the mesh.** `plq_1p.triangulate` +
   `splitTightTriangleSym` is the DEFAULT since 2026-08-18 and keeps its geometry exact without a
   RatPol mesh; the general quadrilateral -- the case that motivated all of this -- is exact
   (worst denominator 56). `convEnvCPLQ` + `ratPolToPlq` survives as `conjCPLQ`'s FALLBACK for a
   rational envelope face, reached only after `conjPieceCPLQ` and the direct face routes decline.
3. **The change is lattice-wide and the classes say so on purpose.** `RatPar`'s header records
   that every property is declared in exactly one place because a property defined in two
   superclasses is fatal in MATLAB and would make `QuaPol < RatPol & QuaPar` impossible. Relaxing
   `mustBeNumeric` on `V` (and consistently on `f`, `den`, `P`, `Ec`) touches every consumer, and
   those consumers do numeric work -- `double(V)`, `norm`, indexing, plotting.

**What would reverse this decision:** a WRONG VALUE traced to mesh rounding, or the rational-face
fallback becoming common again. If it is only cost that bites, the cheaper move is to widen the
A.4/A.5 path so fewer inputs fall back at all -- that is where the code is already going.

## 2026-08-18 (last) — the parallelogram's "remaining 1%" was MY REFERENCE, not the code

Piece 9 of `f*` for `x*y` over `conv{(0,0),(2,0),(2.5,1),(0.5,1)}` is EXACT at all ten probe
points, one cell each -- `BAD 0 of 10`. The singular-quadratic overlap at `(1/2,1/4)` was the only
real defect, and `functionNDomain.singularEdgeCut` closed it.

**The three points I reported as a lingering ~1% over-claim were a broken reference.** The
brute-force sup was a grid over the piece, and for this piece the sup is attained AT A VERTEX --
`(1/4, 7/8)`, on the conic. A grid never lands there, so it reported `2.840439` where the true sup
is `2.875`, and the correct vertex cone read as over-claiming by 1%. Adding the exact vertices to
the candidate set makes all ten agree.

**The reusable part: a grid is not a reference for a sup that a VERTEX attains.** Every conjugate
in this codebase is a sup of an affine form, so its maximiser is at an extreme point far more
often than not -- put the region's own vertices in the reference before reading any disagreement
as a defect. This is the second time in one session that a measurement artefact was reported as a
code defect (the first is the correction entry about unmeasured Step 3 numbers).

## 2026-08-18 (later) — that revert was WRONG: the code was fine, the methods were in the wrong block

**Overturns the entry below, and the mistake is worth more than the fix.** The vertex-plus-arc-
tangency maximum was reverted after `testfunctionNDomain/testMerge` and
`cplqAdapterTest/twoTriangleSquareMaxMatchesNumericSup` went red, and it was written up as a
soundness failure with a guess about redundant conic facets. It was neither. The actual error:

    The class region has no Constant property or Static method named 'quadraticParts'.

`quadraticParts` and `lineMeetsConic` were appended next to `maxAffineOverRegion`, an INSTANCE
method, and then called as `region.quadraticParts(...)`. Moving them into the `methods (Static)`
block is the whole fix. `testMerge` has no assertions, so its "failure" was that exception --
which a single look at the report would have said, and which cost a revert instead.

**The lesson, and it is the reusable part: READ THE ERROR BEFORE THEORISING ABOUT THE MATH.** A
red test in this codebase is as likely to be a MATLAB scoping rule as a defect in the geometry,
and the report says which. The `unionIsExact` gate is soundness-critical, which is exactly why the
temptation was to assume the deep explanation.

**Measured with it back in**, same three folds of the A.4/A.5 quadrilateral:

    fold 3 cells                38 -> 36
    quadFacet_exactAnotInB      63 -> 41
    merges at fold 3             7 -> 9
    TOTAL                      828 -> 804 s

Green: fast 206/0, normal 11/0, testfunctionNDomain + regionTest 17/0, cplqAdapterTest 4/0.

## 2026-08-18 — the A.4/A.5 split is now the DEFAULT

Both objections recorded on 2026-08-16 are gone, measured:

  * `testcPLQ` with the split on: **8 passed / 0 failed in 2273 s**, against 4728 s and one ERROR.
    `testRectBiconj` is one of the eight -- that exception was a casualty of the double leaks and
    nothing in the test or the split was changed to fix it.
  * assembling `f*` for the general quadrilateral: 86 cells / 73 min -> 60 cells / 43 min.

2273 s against 1542 s off is a BUCKET question by the standing rule (2026-08-17), not a blocker:
`testcPLQ` is already in the slow bucket and finishes well inside its timeout. Verified with the
flip alone: fast 206/0, normal 11/0 in 230 s.

`plq_1p.appendTriangle` now gates on `CCA2_NO_A45_SPLIT` (opt OUT); `CCA2_A45_SPLIT` is still
honoured so the two quadrilateral tests that set it keep working.

## ~~2026-08-18 — REVERTED: an exact max over a region with a CURVED facet, for unionIsExact~~ (OVERTURNED, see below)

- **Tried:** replacing `region.impliedBy`'s LP-over-the-linear-relaxation with an exact maximum
  over the region itself (`maxAffineOverRegion` + `holdsOn`). The relaxation drops a conic facet,
  so it is sound but conservative exactly when that conic is what would have cut the violating
  part away -- `quadFacet_exactAnotInB`, 98 of fold 5's refusals and the largest NAMED gate left.
- **The argument, which still looks right:** a linear form on a compact set attains its max on the
  boundary; on a straight edge that is an endpoint (a vertex), on a conic arc it is an endpoint or
  a point where `grad h` is PARALLEL to the form. That parallel condition is affine, so it meets
  the conic in at most two points, in closed form. Vertices plus those points therefore cover
  every candidate.
- **Why it was reverted:** it FAILS. `testfunctionNDomain/testMerge` and
  `cplqAdapterTest/twoTriangleSquareMaxMatchesNumericSup` both go red with it in, and both go
  green the moment it comes out (11/0 with only the A.4/A.5 default flip live). Since the whole
  point of `unionIsExact` is soundness, a gate that admits merges the old one refused and then
  changes a VALUE is refuted, not debuggable by guesswork.
- **Where the hole most likely is, for whoever picks it up:** the argument needs the region's
  `vx`/`vy` to be every corner of its boundary, and needs the region to be the compact convex set
  it is assumed to be. A REDUNDANT conic facet is enough to break the first -- it makes `lin`
  false, so the refinement engages, while the vertex list still describes only the polyhedral
  part. Test that specific case first, on a region built by hand, before touching `merge` again.
- **Do not re-run this as-is.** The next attempt should be measured on
  `testfunctionNDomain/testMerge` FIRST, which is a unit-level merge test and takes seconds.

## 2026-08-17 (last) — item 1's root cause is the MESH VERTEX TYPE, and that is a design change

The 5 same-function pairs whose shared facet `merge` cannot see are two doubles of the same exact
number, one ULP apart. Traced to the source, and it is not `conjConvexOverPiece` (which now carries
whatever it is given exactly) and not `convEnvCPLQ` (which contains no `double(` at all):

    convEnvCPLQ on x*y over conv{(0,0),(3,3),(1,2)} returns envelope vertices
        (0,0)  (1.4142,1.4142)  (1,2)  (1.5,1.5)  (3,3)  (1.5858,1.5858)

`1.4142` is `sqrt(2)` and `1.5858` is `2 - sqrt(2)`, stored as DOUBLES. The reason is one line in
`RatPar.m`:

    V (:,2){mustBeNumeric} % nv x 2 matrix storing unique vertices

The mesh vertices are CONSTRAINED to be numeric, for the whole RatPar / RatPol / QuaPol / QuaPar
lattice, whose header records that every property lives there once and deliberately. So an exact
cevian point cannot survive being stored, no matter how exactly it was computed.

**This is therefore a design change, not a fix**, and it is why it was not attempted here: relaxing
`mustBeNumeric` on `V` (and consistently on `f`, `den`, `P`, `Ec`) touches every consumer of the
lattice, and the classes' own header explains why the properties are declared in exactly one place.
It also has a cheaper alternative worth pricing first: the A.4/A.5 path already keeps its geometry
exact WITHOUT going through a RatPol mesh (`splitTightTriangleSym` -> `plq_1p.triangulate`), and it
is the path the general quadrilateral uses. If the `ratPolToPlq` route is only needed for rational
envelopes, its inexactness costs cells and time but no correctness -- price that before changing the
lattice.

**What is measured, so the size of the prize is known:** on the tri case 5 of 31 same-function pairs
at fold 1 lose a real facet to this. On the A.4/A.5 quadrilateral, which does NOT go through the
mesh, the conjugate is exact (worst denominator 56).

## 2026-08-17 (latest) — testRectBiconj PASSES with the split on: the correctness blocker is gone

`testcPLQ/testRectBiconj` ERRORED with `CCA2_A45_SPLIT` set, and that -- not the runtime -- was
the stated reason the split could not be the default. Re-run against the exactness work
(`domain.mE`/`cE`, `region.limitOfFAtVertices`, `plq_1p.quadPartsOf`, `conjConvexOverPiece`):

    RESULT passed=1 failed=0 incomplete=0

So the exception was a casualty of the double leaks -- 145-digit coefficients and comparisons
`isAlways` could not decide -- and not an independent defect. Nothing was changed in that test or
in the split to achieve it.

**What remains for the default is therefore ONLY cost**, which by the standing rule
(2026-08-17, "DECIDED") is a bucket question rather than a blocker, unless it is on the way to not
finishing. Measure `testcPLQ` with the split on against its 1542 s off / 4728 s on, and note that
the machine is shared so a single timing decides nothing.

## 2026-08-17 (last) — the blow-up, finally COUNTED: two causes, and neither is the arithmetic

`.claude/step3adjacency.m` classifies every same-function pair of fold-1 cells three ways at once:
does `merge`'s own test see a shared facet; do the two carry the same hyperplane with opposite
orientation (numeric rows, so independent of how the constraint is written); and do they actually
MEET in a segment rather than a point. Run with Step 2 exact, on the control case:

    merge sees | same hyperplane | how they meet | pairs
    no         | no              | do not touch  |  5
    no         | yes             | a point only  |  7
    NO         | yes             | A SEGMENT     |  5     <-- facet detection FAILS
    yes        | yes             | do not touch  |  1
    yes        | yes             | a point only  | 14
    yes        | yes             | A SEGMENT     |  6     <-- reaches unionIsExact; 1 merged

**Only the 11 `segment` pairs are ones a merge could ever be right about.** 21 of the 38 pairs meet
at a POINT and must not merge -- the union of two cells touching at a corner is not convex, so
those refusals are correct and were never the problem. That alone retires "578 refusals" as a
number to reason from: most of them are right.

**The 11 split evenly into two DIFFERENT defects:**

1. **5 pairs share a real facet that `ineqs(i) == -ineqs(j)` does not find** -- with exact
   arithmetic, so this is not rounding. `symbolicFunction.eq` is `if (obj1.f == obj2.f)`, a
   STRUCTURAL comparison (its own comment says "change to isAlways"), so the same constraint
   written at a different positive SCALE does not match. `region.normalize1` is supposed to
   prevent that by dividing by `abs(coeffs(f,vars))(end)` -- check whether it picks the same
   term for both operands, which for two differently-written forms it need not.
2. **6 pairs are found and then refused by `unionIsExact`** -- and the fold-1 tally's
   `lin_exactCurvedTest = 6` matches exactly, so all six are `certifiesNonPositive` declining.
   It declines by design outside its hypothesis: a rational `h`, a non-convex quadratic, or a
   relaxation with no vertex. Find out WHICH of the three before extending it.

Both are small and both are now counted rather than guessed. That is the state to start from.

## 2026-08-17 (final for the session) — Step 2 is EXACT now, and Step 3 STILL does not merge

**Read this before spending another hour on exactness.** Three double leaks are fixed and the
whole conjugate is exact -- worst denominator across the six quadrilateral pieces went
`1.2e18 / 9.7e33 / 9.0e16 / 2.6e144 / 6.1e18 / 1.4e145` to `4 / 4 / 7 / 56 / 14 / 56`. On the
control case, Step 3 is essentially UNCHANGED:

              before exactness            after exactness
    FOLD 1    20 cells, 7 fns, 1 merge    21 cells, 7 fns, 1 merge
    FOLD 2    37 cells, 11 fns, 2 merges  37 cells, 10 fns, 2 merges
    FOLD 3    57 cells, 10 fns, 4 merges  53 cells,  9 fns, 4 merges
    fold 3 refusals: noSharedFacet 475 of 511.

**So the ULP-apart doubles were REAL and were NOT the cause.** The entry below identified two
constraints one ULP apart on `4 - 2*sqrt(2)` and concluded the facet test was blinded by
arithmetic. The arithmetic is now exact and the test still does not find the facets. That
conclusion is retracted; what remains true is that the leaks were genuine defects (one of them
produced a y-intercept of `-9.06e-72` where the exact answer is `0`) and are worth having fixed on
their own.

**What is NOT yet explained, and is where to start next.** 15 of 31 same-function pairs at fold 1
carry the same hyperplane with OPPOSITE orientation, measured off `linearForm`'s numeric rows,
while `merge`'s `ineqs(i) == -ineqs(j)` sees nothing. That measurement predates the exactness work
and should be REPEATED first -- if it still holds with exact numbers, the remaining candidates are:

  * `symbolicFunction.eq` is `if (obj1.f == obj2.f)`, a STRUCTURAL test (its own comment says
    "change to isAlways"), so two forms of the same constraint differing by a positive SCALE do
    not match. `region.normalize1` divides by `abs(coeffs(f,vars))(end)`, which is the highest
    term in MATLAB's ordering -- check that it picks the same term for both operands.
  * sharing a hyperplane is necessary for adjacency but NOT sufficient: the two cells may lie on
    the same line yet not touch. The probe above does not distinguish those, and should.

Do that measurement before writing any code. The three hypotheses this session tested were each
plausible, each partly right about a real defect, and each wrong about the cell count.

## 2026-08-17 (latest, and it corrects the entry below) — the facet test cannot match two doubles of the same number

**The entry below is WRONG about the control case, and the way it is wrong is worth keeping.** It
concluded that exactness was not the lever because the "all-rational" case blew up too. That case
is not rational: `conjConvexOverPiece` converts `Q, L, c` and the piece's vertices to DOUBLE by
design (its own lines 59 and 73), whatever it is handed. There was no control.

**What the refusals actually are.** With merge's two heuristics removed, every refusal on that case
became `noSharedFacet` -- 578 of fold 3's 612. A direct check of the fold-1 cells says that is a
FAILING TEST, not geometry: of 31 same-function pairs, **15 carry the same hyperplane with opposite
orientation and `ineqs(i) == -ineqs(j)` does not see it** (12 are seen by both, 4 share nothing).

**And here is why it cannot see it.** Cells 8 and 9 meet along one facet, and carry

    s_2 - 659536895553805/562949953421312       = 5276295164430440/4503599627370496
    s_2 - 5276295164430439/4503599627370496

-- two doubles of the same exact number, `4 - 2*sqrt(2)`, **one ULP apart**. No comparison can
identify them: not the structural `==` this code uses, and not `isAlways` either, because they are
genuinely different rationals. The facet is real and the arithmetic has destroyed the evidence.

**So the chain is one chain, and exactness is the lever after all:** a double enters Step 2 ->
the same quantity acquires two different values in two different cells -> merge cannot match the
shared facet -> nothing merges -> the cell count grows without bound. `domain.mE`/`cE` was one
source of doubles and is fixed; `conjConvexOverPiece` is the other, and it is the one left.

**What the heuristic removal was worth, honestly.** Nothing yet, on the measurements: cell counts
and successful merges are unchanged (20/37/57 cells, 1/2/4 merges). It is still the right code --
two heuristics replaced by `region.certifiesNonPositive`, a sound closed-form certificate -- and it
is what will do the work once the facets can be found again, since `unionIsExact` then becomes the
gate that actually decides. But it fixed nothing on its own and should not be described as if it
had.

## 2026-08-17 (latest) — The CONTROL case: the doubles are real but they are NOT the blow-up

**This corrects the entry below.** `domain.mE`/`cE` really were double arrays and that really was a
defect -- an exact slope arrived as `0.6` and an exact zero y-intercept as `-9.06e-72`, which is a
wrong value, not merely an expensive one, and fixing it took the quadrilateral's worst conjugate
coefficient from `1e144` to `1e33`. But it changed Step 3 **not at all**: cell counts and the merge
tally came back byte-identical.

**The control that decides it.** `.claude/step3cost.m` with `CCA2_STEP3_CASE=tri` runs `x*y` over
`conv{(0,0),(3,3),(1,2)}` through `convEnvCPLQ` + `ratPolToPlq` -- four pieces, ALL RATIONAL, no
A.4/A.5 split, no surds and no doubles anywhere:

    FOLD 1: paired=17 -> cells=20 distinctF= 7   merge: okLinear=1  noSharedFacet=14 quadCutsOther=14 quadMismatch=34
    FOLD 2: paired=36 -> cells=37 distinctF=11   merge: okLinear=2  noSharedFacet=54 quadCutsOther=26 quadMismatch=58
    FOLD 3: paired=60 -> cells=57 distinctF=10   merge: okLinear=4  noSharedFacet=266 quadCutsOther=50 quadMismatch=272 lin_exactCurvedTest=19

57 cells for 10 distinct functions, and **4 successful merges out of 612 attempts**. Same blow-up,
clean numbers. So exactness is not the lever: **`region.merge`'s own gates are**, and they are what
to fix.

**Ranked by what actually fires, measured:**

1. **`quadMismatch`** -- 272 of fold 3's 612. When BOTH regions carry a quadratic constraint and
   they do not share one as a facet, merge demands that EVERY quadratic of A equal EVERY quadratic
   of B, as a cross product, and refuses otherwise. Two adjacent cells each carrying a different
   parabolic arc elsewhere have a perfectly convex union; this refuses them outright.
2. **`noSharedFacet`** -- 266. With clean numbers a large part of this is honest: 10 functions over
   57 cells means groups of ~6, most of whose pairs are genuinely not adjacent. Check before
   attacking it.
3. **`quadCutsOther`** -- 50. The other heuristic: refuse if one region's quadratic meets the other
   anywhere but at a vertex.
4. **`exactCurvedTest`** -- 19. The SOUND gate: `unionIsExact` cannot certify `A subset B'` when a
   constraint it must test is non-affine, so it refuses.

**The fix these point to, and it is one fix.** 1 and 3 are heuristics standing in for the
certificate 4 cannot supply. Give `unionIsExact` that certificate -- "is `max h <= 0` over the other
region" for a non-affine `h`, which for a CONVEX quadratic over a polyhedron is decided exactly by
its vertices plus its recession directions, the region's own curved facets being droppable because
the linear relaxation is a superset -- and then the two heuristics can go, because `unionIsExact` is
the exact criterion (`M = A' n B'` equals `A u B` iff `A subset B'` and `B subset A'`) and refusing
only ever costs compactness, never correctness.

## 2026-08-17 (later) — Step 3's blow-up MEASURED: merge stops working entirely, and doubles are why

Two measurements, on `x*y` over `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}` with `CCA2_A45_SPLIT` on.
Both are reproducible with `.claude/step3cost.m`.

**Merge succeeds ZERO times from fold 2 on.** 190 attempts, no success, while the 29 surviving
cells carry only **7 distinct functions** -- and `distinctF` is 7 at fold 1 too, so the answer
never needs more than a handful of cells. The cell sequence 5, 14, 29, 45, 70, 86 is not cells
being created, it is merges being refused. Merging is also half the cost (182 s of fold 2's
223 s). Refusals by reason: `noSharedFacet` 70, `quadMismatch` 50, `quadCutsOther` 34,
`exactCurvedTest` 30, `exactAnotInB`/`exactBnotInA` 6.

**And the reason merge cannot decide anything is that Step 2 reintroduces DOUBLES.** Worst
denominator by stage: split sub-triangle domains **12**, Step 1 envelope faces **20**, Step 2
`conjugates` **1.2e18 / 9.7e33 / 2.6e144 / 1.4e145**, `maximumConjugate` unchanged. The split is
exact and Step 1 is exact; `plq_1p.conjugate` is where it goes wrong. Piece 4's conjugate
constraints carry `7307585874000779/9007199254740992` (denominator **2^53**) alongside the exact
`30^(1/2)/12 - 15^(1/2)/6 + 3/4` -- the same kind of quantity in two forms, one exact and one its
own double -- and one coefficient is 97 digits long.

**Why this reframes the work.** `merge` finds a shared facet by `ineqs(i) == -ineqs(j)` and
compares quadratics by `~=`. Neither test can succeed when one side is exact and the other is a
double of it, which is exactly what `noSharedFacet=70` and `quadMismatch=50` look like. So the
merge gates are not established as the defect yet: the numbers they are being asked to compare
are. **Fix the Step 2 leak first, then re-measure the tally**, and only judge `unionIsExact` and
the two quadratic pre-checks against clean numbers.

This is attempt 3's pathology surviving downstream of the fix that removed it from the split --
the same lesson the quadrilateral entry records, in a place nobody had looked.

## 2026-08-17 — DECIDED: the split stays opt-in until Step 3's cost is fixed, and the principle behind it

**The user's call, and the reasoning is worth keeping because it settles a whole class of future
questions.**

**Chosen: option (a)** — leave `CCA2_A45_SPLIT` opt-in and attack Step 3's cell blow-up first,
because that is also what the parallelogram's residual 4% needs. Options (b) (diagnose
`testRectBiconj` first) and (c) (flip now) were declined.

**The principle the user attached to it, which is general and outranks the flag:**

> Ultimately all computations have to be CORRECT even if they take a long time. In that case,
> split the unit tests between fast, medium and slow. Correctness is more important than speed —
> **unless** it is so slow that the user does not wait, because getting a timeout is not helpful
> either.

So the reason the split is still opt-in is **not** that it costs 1542 s → 4728 s. A correct path
that is slow gets its test moved down a bucket (`.claude/suite.sh` already has fast / normal /
slow) and keeps running. The only live objection is the **`testcPLQ/testRectBiconj` EXCEPTION**,
which is a correctness question, plus the risk that 4728 s is on the way to "does not finish",
which is the one failure mode the principle does not tolerate either. Both of those are what
Step 3's cost work addresses.

**What this rules out for good:** proposing that a proven-correct path be reverted or gated purely
because it is slower. That is not an available answer here.

## 2026-08-16 (last) — Making the A.4/A.5 split the DEFAULT: deliberately not done

- **Tried:** shipping the split on by default, since it fixes a documented crash and a documented
  wrong answer (`plq_1p.appendTriangle`, `splitTightTriangleSym`).
- **Why it was not kept:** with it on, `testcPLQ` runs 4728 s against 1542 s off (its historical
  time is 1427 s), **and `testcPLQ/testRectBiconj` ERRORS**. That test has no assertions — it runs
  `triangulate`, `maximum`, `biconjugateF` and nothing else — so the failure is an exception, not a
  wrong value. Undiagnosed. Turning the split on therefore trades a documented, LOUD failure on one
  domain shape for a new one on another, which is not a trade to make silently.
- **Before retrying:** fix Step 3's cell blow-up (the cost) and diagnose `testRectBiconj`'s
  exception (the correctness question), in that order. Both are in `TODO.md` with their numbers.
  The switch is `CCA2_A45_SPLIT`; the two quadrilateral tests set it themselves, so the fix stays
  exercised either way.
- **Evidence:** the runs behind those two timings are the `testcPLQ`-only runs of 2026-08-16,
  uncontended in both directions; commit `ba3457d`.

## 2026-08-16 (later) — The general quadrilateral, FIXED on the fourth attempt. What the three failures were each worth.

The method the user prescribed is what worked, and it is worth stating on its own: **do it
symbolically first, and only then reach for explicit formulas.** Attempt 3 failed for exactly the
reason attempt 4 succeeded — it computed the same geometry in double precision and let `sym` turn
the result into `2^53` denominators.

**The fix.** `splitTightTriangleSym` splits a triangle into sub-triangles on each of which cPLQ's
own closed form for THAT sub-triangle's convex-edge count IS the convex envelope, and
`plq_1p.triangulate` emits them as PIECES. A.4 gives one two-convex-edge half (its own form
unchanged and now tight) plus one one-convex-edge half (A.3's rational form, which the `nCE == 1`
branch derives analytically); A.5 gives two two-convex-edge halves that recurse into A.4. Exact
symbolic arithmetic throughout, so the cevian foot — irrational in general — stays a compact surd
(`5/2 − sqrt(5)/2`, `3/2 − 3·sqrt(5)/10`).

**Two things learned while writing it that are worth keeping.**

- **`twoEdgeQuadPlain`'s ± branch search is vestigial: it always returns `s = +1`.** Substituting
  `y = mh·x + qh` into the `s`-form gives `x²` coefficient `mh·denom/denom = mh`, `x` coefficient
  `qh·denom/denom = qh` and constant `0` — for BOTH signs. So the form touches `x·y` along edge `h`
  either way and the touching test cannot separate them; `s = −1` is merely undefined when
  `mh = mw`. cPLQ hard-codes `s = +1` and is right to.
- **The "no split needed" test has a one-line geometric reading.** With `s = +1` the curvature of
  `q1` along the weak edge `AB` is `(sqrt(mh·mw)·dx + dy)² / (sqrt(mh) + sqrt(mw))²` — a perfect
  square — so it vanishes exactly when `AB` is PARALLEL to the cevian direction `−sqrt(mh·mw)`.
  That is the same degeneracy `convEnvCPLQ`'s header describes as "both candidate cevians
  degenerate", arrived at from the other side.

**What the split costs, and where the cost actually is.** The no-split path — every input cPLQ
itself ever had — is 20 ms, which is why the fast bucket did not move. The split paths are 0.3 s
(A.4) and 1.2 s (A.5). What is expensive is what comes AFTER: six pieces instead of two, and
Step 3's cross-piece maximum then takes **73 minutes**, with the cell count running 5, 14, 29, 45,
70, 86. The answer is exact at both levels — 10 of 10 on the per-piece max, 8 of 8 on the assembled
one — so this is a scaling problem in `maxOfList`, not a defect in the split. It has its own
`TODO.md` item.

**And the real reason it is expensive is the ALGEBRAIC DEGREE, not the piece count — which is worth
knowing before optimising the wrong thing.** A.4's cevian foot is irrational, so a split
sub-triangle has SURD coordinates and every symbolic operation downstream works in a quadratic
extension instead of the rationals. Measured on `testcPLQ`, whose domains are general polygons
carrying `x·y`: **1542 s with the split off (matching its historical 1427 s), over 3100 s with it on
and still unfinished when stopped**, uncontended both times — while only two of its six domains
gain a piece at all (2 → 3 and 1 → 2). Two extra pieces cannot explain a 2×+ slowdown; surds can.
`CCA2_NO_A45_SPLIT` therefore exists as an opt-out, and correctness keeps the default.

## 2026-08-16 — The parallelogram's `emptyResult`: TWO defects, both fixed, and two more measured and not taken

Traced from `QuaParCPLQ:conj:emptyResult` down to a single piece, then to a unit reproducer that
runs in about a minute. Piece 9 of `f*` for `x·y` over `conv{(0,0),(2,0),(2.5,1),(0.5,1)}` went from
**6 of 10 probe points wrong or uncovered to 2**, against a brute-force sup over the piece.

**How it was found, which is the reusable part.** `f**` of a bounded domain is finite exactly ON
that domain, and `f**` is a MAX, so EVERY per-piece conjugate must be finite there. Evaluating all
12 groups at six points inside the parallelogram showed three with holes — 9, 11 and 12 — while the
other nine covered everything. That turned "the max comes out empty" into "these three pieces are
wrong", in one cheap measurement. The accumulator and group 11 covered DISJOINT halves of the
parallelogram, so the empty intersection was the honest answer to a wrong question.

**Defect 1 — `region.simplifyUnboundedRegion` declares any region with no finite vertex empty.** A
half-plane has none; so do a slab and the whole plane. And a half-plane is exactly what a TANGENT
vertex produces: the cone at a vertex is built from the two edges meeting there, and when those are
tangent — an arc and its chord touching, which is how a curvilinear piece ends — both half-planes
are the SAME one. Here the cone at `(1/4,7/8)` is `{2x/3 + y ≥ 4/3}` carrying `x/4 + 7y/8 − 1/2`,
and it was being deleted. Now refuted by a WITNESS: a feasible point certifies non-emptiness, while
failing to find one proves nothing, so the old verdict stands whenever no witness turns up.

**Defect 2 — the edge list, again, in its other form.** The piece has 3 vertices and 3 genuine
edges plus a conic touching one vertex: 4 constraints for 3 vertices, so `size(ineqs,2) == nv` calls
a BOUNDED region unbounded, `endNv` comes out `nv−1`, and it is built one edge cell short. Same
root cause as bug 1, without the lens collision — so `conjugateOfPiecePoly` now derives the edge
list for any bounded piece the count mislabels, with the count-matching branches keeping precedence.

**Measured and NOT taken, both at the `isQuad` chord rewrite.** `biconjugateTest`'s own comment
asked for new evidence before re-trying this; here it is, and it says leave it alone.
- Chording **the vertices the conic actually touches** (`region.vertexOfEdge`) instead of
  `vx(1),vx(2)`: makes this piece WORSE, 2 wrong of 10 → 3. The arc's cell grows to swallow a
  region the vertex cone was answering correctly.
- **Skipping the rewrite entirely**: changes this piece not at all, 2 of 10 either way.

**What the last 2 of 10 are.** The chord's edge cell and the arc's claim the SAME region, and the
chord's is checked first and is wrong there (0.0176 and −0.0138 against 0.0333 and 0.0363; the
arc's cell has the right value at both). `f` is a SINGULAR convex quadratic on this piece —
constant along one whole edge, `∇f = 0` at two of the three vertices — so `functionNDomain.getInterior`,
which eliminates `s` between `x = ∂₁f` and `y = ∂₂f`, returns the gradient map's image LINE rather
than a curve separating the two cells. That is the next defect, and it is not in the chord rewrite.

## 2026-08-16 — The A.4/A.5 split as a DOMAIN split: written, measured, REVERTED. The blocker is arithmetic, not structure.

Third attempt at the general quadrilateral, and the first one that failed for a reason none of the
previous analysis predicted. **The structure was right; the ARITHMETIC is what kills it.**

**What was built.** `plq_1p.triangulate` split each indefinite triangle into the A.4/A.5
sub-triangles before emitting pieces, so every sub-piece lands on a path Step 2 already has — a
2-convex-edge sub-triangle where cPLQ's own closed form IS tight, plus a 1-convex-edge one whose
A.3 rational envelope the `nCE == 1` branch derives analytically. The sub-triangles came from
`convEnvCPLQ` itself (its FACES are the split, in x coordinates, recursion done), so nothing had to
be moved out of that file. It worked at Step 1 exactly as designed.

**Why it was reverted: it turns the crash into a HANG.** The quadrilateral's first conjugate ran
for over 45 minutes without finishing, against ~3 minutes before, and had stopped producing output
— stuck inside a symbolic call, with 3.8 MB of `isAlways:TruthUnknown` warnings behind it carrying
coefficients like

    −609298613085773108668343859 / 14507109835375549828038656 .

A hang is worse than a crash here: this repository's own suite driver exists because "a wedged
suite no longer costs the whole run."

**The cause, and it is inherent to the approach.** [COAP] A.4's cevian has slope `−sqrt(mh·mw)`, so
its foot is IRRATIONAL in general. `convEnvCPLQ` is double precision throughout, and `sym` of a
double is EXACT — a denominator near `2^53`. Snapping the new vertices to the simplest rational
within `1e-10` (which the reverted code did, and which is the same numeric→symbolic hygiene
`ratPolToPlq` already applies to coefficients) bounds the VERTEX denominators but not the
downstream ones: the conjugate is a rational function of those coordinates, so a few squarings and
additions carry `1e5` denominators to `1e25`, and MuPAD's `isAlways` then cannot decide anything.

**So the split has to be carried SYMBOLICALLY, not snapped.** The cevian foot is an exact algebraic
number — intersect the line through the far vertex of slope `−sqrt(mh·mw)` with the other convex
edge — and `sqrt` is something the symbolic layer handles natively and keeps small. That means
implementing A.4 (and A.5's smooth-fit line) in symbolic form beside the existing double version,
rather than calling `convEnvCPLQ` for the geometry. It is a bounded piece of work — the formulas
are two line intersections and a curvature test — but it is new code, not wiring.

**What NOT to do next.** Do not re-try any of: installing A.4/A.5 faces as ENVELOPE faces (Step 2
has no rational-envelope branch — 2026-08-15 entry), calling `convEnvCPLQ` for the domains and
rounding (this entry), or leaving `nCE == 2` alone in the hope the error is elsewhere (it is the
envelope itself — 2026-08-15 entry, with the −0.2835 measurement).

## 2026-08-15 — cPLQ's `nCE == 2` envelope is NOT the envelope, and the obvious fix cannot work

**The defect, derived and then measured.** `plq_1p.convexEnvelope1`'s `nCE == 2` branch applies
[COAP]'s single-quadratic form to the WHOLE triangle. That form touches `x·y` along both convex
edges and is a valid convex MINORANT — but Appendix A.4 shows it is tight only over a sub-region,
which is exactly what `convEnvCPLQ`'s `splitTwoConvexEdges` tests for and splits on. This branch
never tests.

Measured on `conv{(2.5,1.5),(2,0),(0,0)}` carrying `x·y`, the piece `triangulate` produces from the
test quadrilateral. cPLQ returns

    0.954915·y − 0.572949·x + 0.427051·x·y + 0.286475·x² + 0.159153·y²

whose minimum over the triangle is **−0.2835**, at `(1,0)`. On that triangle `x ≥ 0` and `y ≥ 0`, so
`x·y ≥ 0` and the affine minorant `0` is admissible: **the true envelope is ≥ 0 everywhere**, and
this is strictly below it. A too-small envelope gives a too-large conjugate, and that is the whole
error — `f*(0,0) = 0.28647` for a truth of `0`, `f*(0.5,1) = 1.00464` for `1`. `convEnvCPLQ` on the
same triangle returns **2 faces** — it does apply the split — with minimum `0`.

The `0.28647` is not a coincidence worth chasing: it is the form's own `x²` coefficient,
`m_h·m_w/(√m_h+√m_w)²` with the two convex slopes `3` and `0.6`.

**AND ROUTING STEP 1 THROUGH `convEnvCPLQ` DOES NOT FIX IT.** Tried, measured, reverted:
`plq_1p.conjugateFunction`'s `nCE == 2` branch reads its envelope's coefficients with
`coeffs(envelope.f.f, vars)` and matches monomials. `convEnvCPLQ`'s A.4/A.5 faces are **RATIONAL**,
so it raises `symbolic:coeffs:NotAPolynomial` outright. **cPLQ's Step 2 has no rational-envelope
branch at all**, and the dispatch keys on the PIECE's `nCE` rather than on the envelope face in
hand — which the routine's own header already complains about for a different reason.

**So the split belongs in the DOMAIN, not in the envelope.** The route that already works for
rational faces is `conjCPLQ`'s `conjEnvelopeViaCPLQ`: it hands each rational face to cPLQ as its own
PIECE via `ratPolToPlq`, and lets `plq.maximum` take the max. The same idea applies here — have
`plq_1p.triangulate` split a 2- or 3-convex-edge triangle into the A.4/A.5 sub-triangles and emit
each as a piece, recursing while `splitTwoConvexEdges` still reports `needsSplit`. Every sub-piece
is then a triangle on which cPLQ's own closed form IS tight, so Step 2 is untouched and every
envelope stays polynomial.

**What that costs, and why it was not done unattended.** `splitTwoConvexEdges`, `splitThreeConvex`
and their helpers (`classifyConvexEdges`, `solveTriangleBF`, `envelopeFromClassified`,
`bilinearFrame`, …) are all file-local to `convEnvCPLQ.m`, so exposing them means moving a
connected web of functions out of a well-tested file — and `triangulate` feeds every Case C result,
so the blast radius is the whole symbolic pipeline. That is a design change with a full
re-verification behind it, not a fix.

## 2026-08-15 — The general quadrilateral's `nCE == 3` wiring: written, measured, REVERTED

**Do not re-land this until cPLQ's Step 2 is fixed.** The `nCE == 3` branch is not the reason the
general convex quadrilateral fails; it is only the reason it fails LOUDLY.

**What was written, and it is correct as far as it goes.** In `plq_1p.convexEnvelope1`, for
`nCE == 3`, build the triangle as a one-face `QuaPol` carrying `x·y` — safe, because reaching that
line with an indefinite `q` means `plq_1p.isCanonicalXY` held, so `q` is EXACTLY `x·y` with no
linear or constant part — call `convEnvCPLQ` (which has [COAP] A.5's `splitThreeConvex`, and
recurses so each half also gets A.4's tightness check), convert with `ratPolToPlq`, and install the
faces as envelope pieces. `plq_1p.conjugate` already loops over envelope pieces and accumulates
`conjfia`, so several per input piece is the normal shape, and Step 3's max over the results is
correct because a sup over a union is the max of the sups. **Measured: 4 envelope faces for the
offending triangle, two quadratic and two rational, all ≤ `x·y`, and `conj` stops raising
`MATLAB:badsubscript`.**

**Why it was reverted.** With the branch in, `f*` of the quadrilateral is WRONG at 4 of 8 probe
points — over-claiming `0.28647` at `(0,0)` where the truth is `0`, `1.00464` at `(0.5,1)` where it
is `1`, and uncovered at `(-1,0.5)` and `(-2,-2)`. **A silent wrong answer is worse than a loud
refusal**, so the crash stays until the thing underneath it is fixed.

**Where the fault actually is, separated by measurement rather than argued.** `triangulate` splits
the test quadrilateral into piece 1 `[2.5 1.5; 2 0; 0 0]` (`nE = 2`) and piece 2
`[2.5 1.5; 0 0; 0.5 1]` (`nE = 3`). Evaluating each piece's own Step 2 conjugate inside Case C:

- **piece 2 gets its 4 envelope faces and Step 2 returns ZERO conjugate cells for it.** The new
  envelope is computed and then discarded, so the wiring buys nothing today.
- **piece 1 — untouched by any of this — carries the whole error**: 6 cells, and every wrong value
  above is its.

**The measurement that localises it, and the trap it avoids.** That same `nE = 2` triangle
conjugated ON ITS OWN via `QuaPol.conj` is exact at 7 of 7 probes — which looks like an alibi and is
not: a single bounded triangle takes the NUMERIC route (`conjBoundedPolygon`), not cPLQ. **The
numeric Step 2 is right on this input and the vendored symbolic one is not.** Checking a suspect
piece "on its own" through the public API can silently change which implementation runs; evaluate
`p.pieces(k).maxConjugate` inside the pipeline instead. `assertStep3MatchesPieces` correctly does
not fire here, because Step 3 does agree with Step 2 — the gate is doing its job, and the fault is
one stage earlier than it looks at.

## 2026-08-15 (last) — BUG 1 FIXED: the edge list, and why the four earlier attempts could not have worked

**What the refactor turned out to be, and it is smaller than the note below predicted.** Four
attempts had all tried to make the LENS fit one of the two slot conventions — free a slot, park a
claim, drop a constraint. None can work, and the reason is worth stating plainly: *the conventions
cannot express a lens at all.* Both say "edge `j` is at `ineqs(j)`" or "at `ineqs(j+1)`", i.e. they
identify an edge by a VERTEX INDEX — and a lens's two edges join the SAME pair of vertices, so no
assignment of slots to constraints can name them apart. The question is not which constraint gets
which slot; it is that slots are the wrong addressing scheme for this shape.

So `functionNDomain.edgeIndexList` derives `eIdx(j)` — the constraint bounding the edge from vertex
`j` to vertex `j+1` — from the geometry, and the three edge-indexed readers take it as an argument.
Two things made this much less invasive than the earlier note predicted:

- **Both loops are in `conjugateOfPiecePoly` itself**, so the list is built in the first and
  consumed in the second with no class field, no signature change on the pieces, nothing to carry.
- **`getNormalConeEdgeQ` and `getNormalConeEdgeQ3` are the SAME routine.** Both build the
  perpendicular to the EDGE'S OWN constraint at each endpoint, oriented by the other endpoint; they
  differ only in which slot they believe that constraint occupies. Given the list they collapse to
  one routine (`getNormalConeEdgeQE`, over `region.coneNormalAt`). Likewise `getSubdiffVertexT2` and
  `getSubdiffVertexT2Q` are identical on these inputs — `T2Q`'s extra third column is only ever
  read when `NCE(:,3)` is non-zero, and neither edge routine sets it. **Three of the "four routines
  that move together" were two routines wearing four names.**

**How the list is decided, and the trap in the obvious version.** A constraint bounds the edge
between two consecutive vertices when both vertices lie on it AND its own curve between them stays
in the region. The second half is not a refinement: on the pipeline's own lens BOTH conics pass
through both vertices, and the redundant one meets the region nowhere else. The first version of
this preferred "the slot today's convention would use" before filtering on that, which handed edge
2 the redundant conic — a constraint bounding nothing. **Filter first, prefer second.** The
preference is what keeps every piece that works today on exactly the indices it has.

**Scope, and why it is safe.** The new path is entered only when two constraints that each bound a
genuine edge STILL share an edge number after `spreadCollidingEdges` has moved everything it can
(`edgesStillCollide`) — the lens signature, reached by nothing else. No constraint is dropped: the
two unsound drops recorded below stay ruled out, and the vertex cones' own feasibility probes
(`ptFeasible`, inside `getNormalConeVertexQ`) need the full constraint set.

**Measured.** Both half-lenses conjugate to 3 cells — two vertex cones plus the arc; the chord's
cell is a ray and drops out — exact against a brute-force sup over the lens at 12 points, 0 wrong.
The three the old code got wrong: `(0,0)` and `(-1,0.5)`, `+inf` where the truth is `0`, and
`(2,-1)`, `0` where the truth is `1/2`.

**Two defects were found by re-reading the finished code, not by a test**, and both are recorded
because the second is the same shape as three other bugs in this file. (i) The two halves of
`getNormalConeVertexQ` were given `eIdx(j)` and `eIdx(j+1)` where the routine's own probe points say
they mean `eIdx(j-1)` and `eIdx(j)` — invisible at `nv = 2`, where predecessor and successor
coincide, and wrong for any larger region. **When a routine's index convention is unclear, read its
PROBE POINTS: they are the geometry it is actually using.** (ii) `edgeIndexList` built `nv-1`
entries for an unbounded region while its consumer walks all `nv` vertices — an index used outside
the guard that produced it, for the fourth time in this codebase.

## 2026-08-15 — BUGS 3 and 4 FIXED: the curved envelope over an unbounded face

Recorded because the DERIVATIONS are the reusable part, and because the shape of the argument is
different from the affine cases already in `convEnvUnbounded`'s header.

- **Wedge, one flat ray `d1` and one convex ray `d2`:** `co q` is `q` with its CROSS TERM deleted,
  `q(v) + α·g1 + β·g2 + β²·A22/2`. The derivation is the file's own method — parametrise the affine
  minorants by their gap parameters and minimise the gap at an arbitrary target — but with one
  crucial difference: **the minimiser DEPENDS on the target point** (through `β0`), which is
  exactly why this envelope is not affine and why the affine argument does not extend. `α → ∞`
  forces `A12 ≥ 0`; a negative `A12` makes `d1 + t·d2` a recession direction of negative curvature,
  so the envelope is `−inf`, and that is now reported rather than answered.
- **Half-strip convex along the ray:** it separates **only when the base edge is Q-ORTHOGONAL to
  the ray**. Then `co q = q(v1) + s·(q(v2) − q(v1)) + t·⟨∇q(v1), d⟩ + t²·(d'Q d)/2`. `w'Q d ≠ 0` is
  refused loudly, with the value in the message — the two directions interact and the minimising
  minorant moves with the target in both coordinates, so there is no separable answer.
- **Checked against the fixtures, not fitted to them:** `co(x·y + I_K) = y²` on `{0 ≤ y ≤ x}` and
  `co(−x²+y²) = −x + y²` on `{0 ≤ x ≤ 1, y ≥ 0}` are what the general formulas produce, and both
  match what `unboundedFaceTest` derives by hand. `unboundedFaceTest` 18 / 0, from 16 / 2.
- **What is still not covered:** a wedge with BOTH rays convex, and a half-strip whose base edge is
  not Q-orthogonal to its ray. Both are refused, not approximated.

## 2026-08-15 — BUG 2 FIXED: a tangent built where the gradient vanishes

**A vanishing gradient at a cone's apex is a recurring failure mode in `region.m`, and this is the
second routine to fall into it.** `SUPPORT_MATRIX.md` §8.2(e) records the first:
`simplifyUnboundedRegion` decided emptiness from probe directions built out of constraint SLOPES at
a vertex where the split conic's gradient vanishes, and `region.witnessAwayFrom` was written for
it. Same input, same trap, different routine.

- **What it was:** `removeTangent` takes a quadratic constraint active at a vertex, builds the
  TANGENT LINE to it there, and deletes any constraint equal to that tangent as redundant. At the
  APEX OF A CONE the quadratic's gradient VANISHES — there is no tangent line, every direction is
  tangent — and whatever it computes is meaningless. It then deletes a constraint that matches.
- **Measured, on the 4-cone fan:** the assembled cell carrying `s1²/4 + s2²/2` lost its `−s1 ≤ 0`,
  keeping only `{s2 ≤ 0, s2²/2 − s1² ≤ 0, s1² − 2s2² ≤ 0}` — two constraints **blind to the sign of
  `s1`** — so the region became symmetric under `s1 → −s1` and claimed the mirror wedge.
  `f*(-3,-2.4)` came back `5.130` for a truth of `4.500`. It is now `4.5`; `conjCPLQTest` 25 / 0.
- **Fix:** refuse to conclude anything from a vanishing gradient.
- **Cleared on the way, so nobody re-checks them:** `region.redundantSubset` certifies nothing
  redundant on that constraint set (correct), and `simplifyUnboundedRegion` leaves the constraint
  alone. Step 2 is right too — each primal piece's own conjugate has the correct quadrant
  constraints, and the per-piece max at that point is the truth.
- **How it was found, and this is the transferable part:** by bisecting the pipeline rather than
  reading it. Dump the cells carrying the offending quadratic immediately before and after
  `mergeL`; the constraint is present before and gone after. Then feed that exact region to each
  simplification `mergeL` applies, in turn, and see which one drops it. Three routines, one run.

## 2026-08-15 (later) — BUG 1: three defects fixed, and one attempted fix that is UNSOUND

**Fixed, and each was necessary before the next was visible.** The lever that made this tractable
was a unit-level reproducer: the half-lens `{(s1+s2)² ≤ 4·s2, s2 ≤ s1}` carrying `s1`, built by
hand as a `region` and conjugated against a brute-force sup over its own boundary. No pipeline,
seconds per run, and it went from **2 identical wrong cells to 3 cells exact at all 10 probe
points**. Build that first next time; the pipeline runs took 10–40 minutes each.

1. `getEdgeNosInf` numbers an edge by one of its endpoint VERTICES, so a LENS — two edges joining
   the same pair — gets one number for both, and the last-write-wins scatter destroys one.
2. `getNormalConeVertexQ` (the routine that builds a vertex cone from the CONSTRAINT's own tangent
   rather than the chord) indexed its second constraint as `j+1` unwrapped, so it raised
   `badsubscript` on any BOUNDED region — which is why its only caller gated it behind a
   constraint COUNT and sent every bounded region to the polyhedral routine instead. Wrapped
   cyclically; identical to `j+1` for the unbounded layout, so nothing that worked changes.
3. `biconj` handed its second conjugation the curved MESH `conj` has returned since 2026-08-13,
   and `quaPolToPlq` refuses a curved domain. It now asks `conjCPLQ` for the symbolic form.

**UNSOUND, and reverted: freeing an edge slot by dropping the constraint holding it.**

- **Tried:** the lens's two edges need slots 1 and 2, which are held by constraints with a single
  vertex on them. Dropping constraints with `nOn ≤ 1` on a bounded region frees them.
- **Why it failed:** a constraint active at exactly one vertex of a convex region can still be
  ESSENTIAL. Removing it enlarges the piece, and an enlarged piece of `f*` has a SMALLER conjugate
  domain — so `f**` loses coverage somewhere else. Measured: with the drop, `f**` is exact at
  `(0.25,0.25)` and `(0.1,0.1)` and `+inf` at `(0.9,0.6)` and `(0.6,0.6)`; without it, exactly the
  other way round. Both are 5 of 7 and only one is sound.
- **A second, independent unsoundness found on the way, and worth its own note:** the boundedness
  test written for that drop read the vertex list AFTER `removeInfV`, which deletes the `±intmax`
  box-clip vertices that are the only mark of an unbounded region — so every region looked
  bounded. Fixed on its own merits (read it before `removeInfV`, which is the codebase's own
  convention), and it did NOT rescue the drop: the harmed piece is genuinely bounded.
- **A FOURTH attempt, also reverted: PARKING the slot claim instead of dropping the constraint.**
  Since only `endNv` constraints are ever read as edges — and on a lens that is ONE — the arc and
  its chord fight for one slot, and the scatter's last-write-wins hands it to the chord. Of the
  two, only the ARC has a two-dimensional dual cell (every point of an arc has its own normal, so
  it sweeps a cone; a straight edge has one normal and contributes a ray). So: give the arc the
  slot and PARK the chord's claim above the edge slots, dropping nothing — the region as a
  constraint set is unchanged, and only the label "this is edge number k" moves.
  The reasoning is sound and the result was worse: `conjugateOfPiecePoly` returned NO pieces at
  all, so it breaks something beyond the lens. Not diagnosed further.
- **Before retrying:** do not look for a better rule for dropping, and do not re-try parking
  without first finding out which OTHER piece the parking breaks. Give `conjugateOfPiecePoly` an
  explicit EDGE LIST instead of a count with two conventions (`endNv = nv` or `nv-1`; edge `j` at
  `ineqs(j)` or `ineqs(j+1)`). It cannot be done in that routine alone — `j` indexes
  `getNormalConeEdgeQ`/`Q3`'s output and `getSubdiffVertexT2`/`T2Q`'s `subdE` simultaneously, so
  all four move together. That is why the original "derive the edge list from the geometry" note
  underestimated the job.

## 2026-08-15 — Two of the five "remaining bugs" were described WRONG. Measure before fixing.

Both descriptions had been written from a symptom and carried forward as fact. Each cost an
attempt before measurement refuted it. Corrected shapes below.

> **Both are now resolved** — bug 5 fixed the same day (see its own entry), bug 1 taken from 0 to
> 5 of 7 probe points. This entry is kept for the corrected DIAGNOSES, which is what it is for.

### BUG 1 — "conjugateOfPiecePoly returns the conjugate of the chord"

- **What the record said:** the routine decides how many edges a piece has from
  `size(d.ineqs,2) == d.nv`, a COUNT standing in for "is this region unbounded", so the half-lens
  takes the unbounded convention and its arc is never read as an edge. Fix: derive the edge list
  from geometry.
- **What is actually true, measured** (instrumented dump of the pre-scatter constraint list):
  the half-lens arrives with `nv = 2` and FIVE constraints, and `edgeNo = [3 1 2 2 2]` —
  **three** constraints claim edge slot 2, all with both vertices on them:
      con 3: (s1+s2)^2 - 4*s1     con 4: (s1+s2)^2 - 4*s2     con 5: s2 - s1
  `getEdgeNosInf` numbers an edge by one of its ENDPOINT VERTICES, and a lens has two edges
  joining the SAME pair — so they are indistinguishable to it. The scatter
  `d.ineqs(edgeNo) = d.ineqs` is last-write-wins, so one edge is destroyed and another stored
  twice. Feed the arc first and both slots end up holding the CHORD; feed the chord first and both
  hold the arc. That is the whole of the "conjugate of the chord" symptom, and the count test is a
  consequence, not the cause.
- **Also measured, and this is why fixing the numbering is not enough:** hand-build the lens with
  ONLY its two genuine edges (`{(s1+s2)^2 <= 4*s2, s2 <= s1}` carrying `s1`, `nv = 2`,
  `nineq = 2`) and `conjugateOfPiecePoly` still returns **1 cell**, not the 4 the piece needs
  (2 vertex cones + 2 edge cells). So the downstream cell generation does not handle a 2-vertex
  CURVED region either.
- **Tried, twice, both REVERTED:** (i) spreading colliding genuine edges over distinct numbers;
  (ii) the same, plus dropping constraints that bound no edge — a constraint whose curve between
  the shared vertices leaves the region (which is what con 3 is, and which LP redundancy cannot
  see because it is not linear), and constraints active at a single vertex of a bounded region.
  Both make the lens's slots correct — measured, `nineq = 2` with the arc and the chord in
  distinct slots — and both take the second conjugation from WRONG VALUES to
  `QuaParCPLQ:conj:emptyResult`, no pieces at all.
- **Before retrying:** the edge numbering is necessary and not sufficient. Do the downstream half
  FIRST — make `conjugateOfPiecePoly` produce 4 cells for a bounded 2-vertex region with one
  curved edge, checked on the hand-built lens above, which needs no pipeline. Only then re-apply
  the numbering fix. Note also that `size(d.ineqs,2) == d.nv` is used a SECOND time, to choose
  between the polyhedral and the quadratic-aware normal-cone routines
  (`getNormalConeVertex`/`getNormalConeEdgeQ3` vs `getNormalConeVertexQ`/`getNormalConeEdgeQ`),
  so changing the count changes which of those runs.
- **Separately: the pinned test no longer fails the way it says it does.** Since `conj` began
  returning a MESHED QuaPar for a bounded multi-face domain (2026-08-13),
  `biconjugateOverATwoFaceSubdivisionIsTheEnvelope` ERRORS at `quaPolToPlq:curvedEdge` — the
  second conjugation is handed a curved QuaPar and `quaPolToPlq` requires a polyhedral domain. The
  symbolic route reaches the lens defect only when forced. Both need fixing; they are different.
- **Refuted on the way:** `convEnvCPLQ` is NOT a route to `f**` here. `f** = conv f` for a compact
  domain, but `convEnvCPLQ` is Step 1, a PER-TRIANGLE envelope with no cross-piece hull — measured
  on this input, it returns the per-piece envelopes (0.25 at (0.5,0.5) where the truth is 0).

### BUG 5 — "a piece that spans TWO sub-arcs of the SAME conic"

- **What the record said** (written 2026-08-14, from the piece's vertex list): piece 4 `src [1 6]`
  has g1's arc cut twice and spans both sub-arcs, so its one curve slot holds one and the other
  becomes a chord.
- **Measured, and it is false:** evaluating piece 4's own conic at its vertices gives
      V1 (-2.960609, 0.9606088) -> 0        V5 (-2.744821, 1.372827) -> 0
      V4 (-2, 2) -> +0.149        V2 -> +0.161        V3 -> +0.327
  with a tolerance of 8.8e-07. Only V1 and V5 are on it. The two curved edges are on **different**
  conics: piece 4's own slot holds the SPLITTING curve `{g1f1 = g2f6}`, and the flattened edge
  `(-2,2) -> (-2.744821,1.372827)` is a sub-arc of **g1's arc**, which its neighbour piece 5
  `src [2 6]` carries curved on `[-0.015625 -0.03125 -0.015625 -0.25 0.25 -1]`.
- **So it is the ordinary two-DIFFERENT-arcs case**, which `splitTwoArcPiece` exists for — and the
  cell does carry the arc going in: `src [1 6]`, `nv = 4`, `arcPos0 = 2`, arc
  `(-2,2) -> (-3.125,1.125)`. But `splitTwoArcPiece` is called exactly ONCE on this input and not
  for this cell, so the loss happens before it: either the crossed-arc restoration finds no
  position for the sub-arc (`findArcPosition` returning 0) or the piece is emitted before reaching
  it.
- **AND THE CAUSE IS NOW LOCATED, by that instrumentation.** For `src [1 6]` on that shift:
      HITS: 2 hits on edges [2 3] at (-2.744821,1.372827) and (-2.960609,0.9606088), arcPos0 = 2
      RESTORE half: nv=5 curveAfter=5 **p=4** straightCut=0
  So the crossed-arc restoration works perfectly: it FINDS the inherited sub-arc at edge 4, the
  splitting curve is genuinely curved (straightCut=0), and the half is correctly handed to
  `splitTwoArcPiece(half, 4, arcEc0)`. **`splitTwoArcPiece` then returns it UNSPLIT**, which its
  own header says it may do — "if neither works the piece is returned unsplit ... the assembly's
  own arrangement check is what would catch a dropped arc". It did, three stages later, as the
  orphan.
- **Why it finds no cut: THE TWO ARCS ARE ADJACENT.** They share vertex V5 (arc at edge 4, split
  curve at edge 5, `nv = 5`). Its two candidate chords are `arcPos+1 -> c+1` = V5->V1 and
  `arcPos -> c` = V4->V5 — but those ARE edges 5 and 4 themselves, so one chain comes out with 2
  vertices and the `numel(chain) < 3` guard skips both. The `nv == 3` fallback that handles the
  shared-vertex case (cut from the shared vertex to the midpoint of the opposite straight edge)
  does not apply at `nv = 5`.
- **FIXED 2026-08-15, exactly as prescribed:** the shared-vertex case is generalised to `nv >= 4`
  with the ordinary diagonal from S to a NON-ADJACENT vertex. On this piece S = V5 and `V5 -> V2`
  gives chains {2,3,4,5} (carrying edge 4, the inherited arc) and {5,1,2} (carrying edge 5, the
  split curve): one arc each. Guarded by `insideStraightHull` like the existing candidates, and
  each half passed through `splitAtReflexVertex`. **Seeded sweep 17 exact / 0 wrong / 1 errored →
  18 / 0 / 0**; `maxQuaParTest` 29 / 0; fast bucket 203 / 0. Pinned by
  `arcVsArcSplitsTwoADJACENTArcsOnAPieceWithADiagonal`.
- **Worth keeping:** the symptom pointed at the wrong stage twice. A straight edge facing an
  identical CURVED one at 8e-16 reads as a clip dropping a conic; the actual cause was a
  subdivision that declined to cut and announced it only by leaving an arc flattened. When a piece
  reaches assembly with an edge its neighbour calls curved, ask which producer returned it UNSPLIT
  before looking at the clip.
- Do NOT write a same-conic sub-arc splitter; that shape has never been observed, and the code for
  it was written and removed on 2026-08-15.

## 2026-08-14 — RESOLVED: the two-arc ray split, and why two attempts at it failed

- **Outcome:** `arcVsArcRefusesAnUnboundedTwoArcSplit` is GREEN and `maxQuaParTest` is 28 / 0 (from
  25 / 1). The entries below about the ray split stand as history; the construction was right in
  outline all along, and **neither reverted attempt failed for the reason it was thought to**.
- **What actually blocked it, both times.** Two defects, neither in the split:
    1. `pieceRecessionRays` took the parabola's axis from `arcNullDirs`, which solves `d·Q·d' = 0`
       *exactly* and returns **nothing** when `b²−4ac` comes out negative — which is what a
       floating-point parabola's `Q` does about half the time, being only semidefinite up to
       rounding (measured `−2.78e-17` on the pinned fixture's first arc). The derived chord was
       then silently never emitted, so the half's constraint region stayed a slab open at BOTH
       ends and `reccConeViolation` refused it. This is why **check (5) of the six did not
       separate the cases**: it was answering a question about a constraint set that was not the
       piece's. `parabolaArcFrame` has always taken the smallest-magnitude eigenvector; do that.
    2. `curveAfter ≠ 0` was read as "this edge is curved" in five places, when `boundedPiece` also
       sets it for a STRAIGHT splitting curve (`curveEc` all zeros). Two of those places call
       `parabolaArcFrame` on the zero conic and raise `degenerateAxis` — **so `MAXQP_ASSERT`
       crashed on three of the four arc-vs-arc fixtures, and the invariants that would have named
       all of this never ran on the inputs that needed them.** An invariant that errors is off.
- **And the seeded shift `[1.4979 3.6486]`, blamed on the split twice:** the split's halves pass
  all three exact invariants. The wrong 0.3531 came from a *different* piece, `src [2 4]`, whose
  only finite edge `polyConstraints` had dropped (defect 2) — the refusal had merely been masking
  it. It is now refused explicitly; see the next entry.
- **Lesson worth keeping:** when a gate refuses a construction you have independent reason to
  believe is right, suspect the gate. Both reverts were correct decisions on the evidence
  available, and the evidence was wrong because the tooling was broken in a way that was silent.

## 2026-08-14 — "`{f1=f2}` must be a PAIR of parallel lines here" (diagnosis, refuted by measuring)

- **Tried:** explaining the last wrong answer on seeded shift `[1.4979 3.6486]` — piece `src [2 4]`
  carrying g1 face 2's ZERO quadratic while g2 face 4 beat it by `+Inf` along its own ray — as a
  subdivision gap: `{f1=f2}` is a degenerate conic, so it can be a pair of parallel lines, and a
  half strictly on one side of the line `splitCell` cut along could be crossed by the other, which
  would leave through the recession cone and contribute no finite crossing to find. That story fits
  every symptom, and a guard was written for it (`assertWinnerHoldsAtInfinity`).
- **Why it failed:** the conic says otherwise. For that cell
  `diffRow = [0 0 0 0 0 0 0 −1.4979 −3.6486 5.4652]` — **its entire quadratic part is zero**, so
  `{f1=f2}` is a SINGLE straight line and nothing straddles it. `delta = 0`, `det3 = 0`,
  `eig(Q) = [0 0]`.
- **What it actually was:** the piece's only two vertices ARE the two crossing points, so both lie
  ON that line; `assignSide` had nothing to read the winner from and fell back to a centroid which
  is on the line too. The winner of a whole unbounded piece came out of floating-point noise.
  `assignSide` now reads such a piece in its RECESSION CONE, sharing the probe
  `assignSideFromCone` has used since it was written for the same problem in
  `splitUnboundedAtOneCrossing`. That shift now assembles CORRECTLY: the seeded sweep is 17 exact /
  0 wrong / 1 errored of 18, from 16 / 0 / 2.
- **Before retrying:** the guard is kept as a cheap exact backstop and currently fires on nothing.
  Do not read its existence as evidence that the pair-of-lines case occurs — no input has ever been
  observed to produce it. **Classify the conic before theorising about its shape**; one line of
  `eig(Q)` would have saved the detour.

## 2026-08-14 — Why check (5) may have failed to separate the two-arc ray split's cases

- **Tried (2026-08-13):** "each half's recession cone must equal the cone its own rays span" was
  one of the six checks applied to the unbounded two-arc ray split, and it did not separate the
  pinned fixture from the seeded shift that assembles to a wrong value. It was therefore recorded
  below as exhausted, and the search was pointed at WHERE the cut starts instead.
- **Why that conclusion may be premature:** the check is implemented by `pieceRecessionRays`, and
  until 2026-08-14 that routine derived the arc's chord — the constraint that makes a concave-side
  arc-piece compact at all — by asking which side the piece's other VERTICES fall on. That is the
  same rule `DECISIONS.md` already records as unsound one level up, in `QuaPar.chordCuts`, where it
  killed two green tests: a lens has both vertices ON the chord, so they say nothing. It also had
  no gate on *when* a chord may be emitted, so a piece that genuinely straddles the chord line
  could be handed one and be reported compact when it is not. Either way the check was answering a
  question about a constraint set that was not the piece's.
- **Corrected:** both questions are now settled by the conic itself. Along the chord `X0 + t·ch`
  the conic restricts to `q(t) = A·t·(t−1)` with `A = ch·Q·ch'`, because both endpoints are on it.
  `A ≤ 0` means the chord's interior is inside the kept side, the piece straddles the line and no
  chord may be emitted (the lens, and every convex-side face). `A > 0` means the piece touches the
  line only at the arc's endpoints and lies on the side the arc's own interior points are on,
  reached along the parabola's axis from the chord midpoint. The vertex test survives as a veto.
- **Status: CONFIRMED, and it was the whole of it.** Correcting the axis derivation is exactly what
  turned both halves of the pinned fixture's cut from "cannot be proved" into "proved". See the
  RESOLVED entry at the top.
- **Evidence:** `arcVsArcRefusesAnUnboundedTwoArcSplit` green; `maxQuaParTest` 28 / 0; fast bucket
  202 / 0.

## 2026-08-14 — A newly minted OUTGOING ray was given sign +1 (a live bug, not a dead end)

- **Recorded here because the reasoning is easy to re-derive wrongly.** `polyConstraints` reads a
  ray's outward normal as `sign · rot90ccw(direction)`, and both `dirIn` and `dirOut` store the
  direction pointing from the apex OUT to infinity. A piece is walked CCW with its interior on the
  LEFT, so the incoming ray is traversed along `−dirIn` and takes `+1`, and the outgoing ray is
  traversed along `+dirOut` and takes `−1`.
- **The bug:** both branches of `splitCell` that mint an escaping ray wrote `+1` for the OUTGOING
  one. That flips the kept half-plane to the far side of the ray's line, so the piece's constraint
  region is the reflection of its true region across that line — over-extended on one side, short
  on the other, which is the shape every far-field wrong answer here has had.
- **What NOT to do:** do not "fix" the INHERITED signs the same way. A sign is a property of the
  `P{k}` entry the ray came from, not of its role, and a face whose whole boundary is two rays
  sharing one apex can legitimately carry the same sign on both — that generalisation is itself a
  recorded fix, with its own test.
- **Evidence:** derived from `polyConstraints`' own HISTORY note; spotted independently while the
  ray split was being reverted and recorded in that commit. `newRaySign` now states it once.
  **Verified 2026-08-14**: no regression on any bucket, and an independent review confirmed the
  derivation (the previous `+1` gave the two halves of a split the *same* outward normal across
  their shared ray — overlap on one side, hole on the other).

## 2026-08-13 — Making the arc's chord a REAL EDGE by splitting the neighbour (option B)

- **Tried:** Nothing was built. It was proposed — and recommended — as the way to close a face
  whose arc is concave towards it, resolving the open choice left by the entry below.
- **Why it failed:** Geometry, checked before implementing. The chord runs through the **interior
  of the neighbour** on the other side of the arc, so making it a real edge splits that neighbour
  and leaves the offending face's own edge list unchanged. It cannot supply the missing constraint
  to the face that needs it. No arrangement can: the chord is not on that face's boundary.
- **Before retrying:** Do not. The chord must be **derived per face**, which is what
  `QuaPar.chordCuts` now does — and it resolved the whole far-field defect. Two soundness rules
  came out of building it, both learned by breaking tests: the side must be read from a point just
  INSIDE the face (the arc's midpoint stepped along the inward normal), never from the face's other
  vertices — a lens has both vertices ON the chord, so they say nothing and the wrong side is
  chosen; and a chord is emitted only when every other vertex and ray direction agrees, i.e. when
  it is redundant for the face, so it can never shrink a face below its true region.
- **Evidence:** `QuaPar.chordCuts` + `SUPPORT_MATRIX.md` §4.1. The vertex-based side rule broke
  `maxQuaParCurvedMatchesGroundTruthOnRandomQuadrilaterals` and
  `maxQuaParSplitsACellThatAlreadyCarriesAnArc`; the interior-point rule turns both green.

## 2026-08-13 — Splitting an UNBOUNDED two-arc piece with a ray (implemented, reverted)

- **Tried:** An unbounded piece carrying two arcs cannot be separated by a chord (no chord closes
  an unbounded piece), so: cut with a RAY from the vertex between the two arcs, along a direction
  the piece recedes in. Each half then keeps one arc, one original ray and the new one. It works —
  `arcVsArcRefusesAnUnboundedTwoArcSplit` assembles with it.
- **Why it failed:** It makes seeded shift `[1.4979 3.6486]` assemble to a value wrong by 0.3531.
  A silent wrong answer is worse than the refusal it replaces, which is what that test's own
  comment says, so it was reverted.
- **Before retrying, fix:** Six checks were tried and NONE separates the good case from the bad:
  (1) the ray recedes every straight constraint; (2) it stays inside BOTH arcs for its whole
  infinite length — necessary, a ray did leave through an arc; (3) each half admits a point just
  inside its own arc — necessary, it caught a mis-oriented pair; (4) the new ray's SIGN (an
  OUTGOING ray needs `-1` in `polyConstraints`' `sign*rot90ccw(dir)`, not the `+1` the
  escape-to-infinity branch uses — that branch looks wrong the same way and deserves its own
  check); (5) each half's recession cone equals the cone its rays span; (6) scoping to the
  half-strip shape. Start past all six: the open question is WHERE the cut starts, not which
  direction it takes — test whether the vertex between the arcs is the right starting point at
  all, or whether the cut must start on one of the arcs.
- **Evidence:** `maxQuaParTest/arcVsArcMatchesGroundTruthOverRandomShifts` (seed 20260803, N=18)
  versus `arcVsArcRefusesAnUnboundedTwoArcSplit`; commit "Revert the unbounded two-arc ray split".

## 2026-08-13 — "Equal areas and a shared polyline prove two halves tile" (refuted reasoning)

- **Tried:** Concluding, from a measurement, that a two-arc split was NOT responsible for a
  coverage hole: the two halves shared their cut polyline bit-exactly (17 digits), were both CCW,
  and their polygon areas summed to the parent's exactly. That sent the next session's search into
  assembly instead.
- **Why it failed:** Both facts were true and the conclusion was wrong. Area says nothing about
  which side of a **bent** boundary a point falls on. The halves did not tile: the cut polyline
  left one of them REFLEX at the bend, and every point-location test here reads a face as an
  intersection of half-planes, which is exact only for a CONVEX face.
- **Before retrying:** To decide whether pieces tile, use the per-edge cross product at the
  disputed point, not area. `splitAtReflexVertex` now splits such a half along a diagonal.
- **Evidence:** `arcVsArcDoesNotCrashOnSeededQuadSplits`, fixture 1, the point
  `(0.998629534754574, -0.0523359562429444)` — outside both halves, inside the parent.

## 2026-08-13 — Subdividing a bounded arc-piece whose constraint region is non-compact

- **Tried:** Twice, on the two quadrilateral-fixture pieces `src [1 2]` and `[1 6]`.
  Cutting the piece with a line parallel to the arc's chord, just past the arc's
  deepest point, to close the non-compact side.
- **Why it failed:** Both halves came back still non-compact, and the straight half
  wrongly kept the arc tag from `clipPolyHalfPlane`. The second attempt settled the
  reason: for a piece whose arc is **concave towards it**, the constraint set is a
  wedge intersected with the **outside** of a parabola. A cut parallel to the chord
  leaves the arc-side sliver still receding along the two side edges' own direction,
  which neither the cut nor the parabola blocks.
- **Before retrying, fix:** Do not retry any subdivision. This is a **representation**
  decision, not a subdivision bug. The constraint that would close the region is the
  arc's own **chord** — redundant for the region (which lies entirely on one side of
  it) but not one of its edges. So either a face may list a redundant bounding conic,
  or the chord becomes a real edge by splitting the neighbour across it. Pick one of
  those two before writing any more clipping code.
- **Evidence:** `maxQuaParTest`, quadrilateral fixture, pieces `src [1 2]` and `[1 6]`.
  (A separate `fixArcTag` fix has since stopped the losing half carrying the arc's
  constraint, but that does not make the cut work.)

## 2026-08-04 — Guarding against non-compact arc-pieces (four guards, all reverted)

- **Tried:** Four different guards, after establishing that arc-vs-arc results are
  correct near the arcs and wrong in the far field:
    1. Piece-level compactness guard on `triHalf` output.
    2. Piece-level compactness guard on every bounded curved piece.
    3. Post-assembly check: "a bounded face admits a far point", using QuaPar's own
       EC + P signs.
    4. Post-assembly check: "two faces disagree at a far point".
- **Why it failed:** (1) and (2) used the wrong sign convention and rejected valid
  faces — (2) errored on all three fixtures. (3) and (4) **detect correctly**, but
  they also fire on `[-1 0.75]` and `[2 -0.5]`, because those results genuinely *are*
  non-compact far out. A backstop that refuses non-compact results is therefore
  correct and would turn every arc-vs-arc result red.
- **Before retrying, fix:** Do not add a guard. The fix must be **upstream**: give an
  unbounded quadratic face's arc sub-pieces their **ray** boundaries, so they stay
  unbounded and compact-as-faces, instead of closing them with straight edges. That
  is a rework of the arc-vs-arc clip/split, not a guard.
- **Evidence:** Ring sweeps of `g.eval` against the pointwise max — `[-1 0.75]` 0/60
  wrong at radius 8 but 2/60 at radius 30 (worst 6.1); `[2 -0.5]` 2/60 at radius 8,
  11/60 at radius 30 (worst 33.8). `[0.5 0.5]`: assembled FACE 15 carries g1-face-1's
  quadratic over a tiny triangle near (-2,2) yet admits (-3.98,0.61) two units away,
  overlapping the correct zero face (FACE 9).

## 2026-08-03 — `decideWinner`: sampling `diff` at 1-D stationary points

- **Tried:** Making `decideWinner` sound on unbounded cells by additionally sampling
  `diff = f1 - f2` at the 1-D stationary point of each edge and each ray.
- **Why it failed:** Sound and strict, but it does not fix the defect and was
  confirmed neutral on all three fixtures. `{diff = 0}` is a parabola (rank-1 `Q`,
  which `splitCell` already asserts) whose branches run to infinity along `Q`'s null
  direction, strictly *between* the cell's two bounding rays. The sign change is off
  the 1-skeleton entirely, so no amount of sampling along edges and rays sees it.
- **Before retrying, fix:** Do not attempt a probe-based patch. A real fix needs
  **both** a sound sign-over-the-cone test in `decideWinner` **and** a `splitCell`
  that can cut a cell along a parabola entering and leaving at infinity — the current
  one asserts exactly two finite crossings, and here there are zero.
- **Evidence:** Cell `src[4 3]` stores `winner=g2` but `g1=0 > g2=-8.648` at the
  interior point `s=(-5.4843,1.5866)`; symmetrically `src[6 5]` stores `winner=g1`
  while g2 wins at `s=(-1.1298,4.1007)`. 4/68 samples, worst 8.648. Kept out of the
  commit.

## 2026-08-03 — Returning both components from a disconnecting curved cut

- **Tried:** Replacing `clipPolyByConic`'s blanket refusal of a disconnecting curved
  cut with "return both components".
- **Why it failed:** Two components are **unrepresentable**. A separated survivor
  would need the cutting conic running to infinity, i.e. an unbounded curved edge,
  which `QuaPar` cannot hold. Separately, the premise was wrong: the cut measured
  **connected** on all three fixtures, so the single-cell construction was right.
- **Before retrying, fix:** `Do not retry` unless `QuaPar` gains unbounded curved
  edges. The blanket refusal was instead replaced with a real connectivity test,
  which now refuses only the genuinely disconnecting case — and that case does not
  arise on these fixtures.
- **Evidence:** `partitionReport` OK on all three fixtures after the connectivity
  test landed.

## 2026-08-13 — Refuted diagnoses (the analysis was wrong, not just the fix)

Each of these was written up in detail and then disproved. Re-deriving any of them
from the same symptoms is the failure mode this entry exists to prevent.

- **"Piece 5 (`src[2 2]`) is over-extended and its ray should terminate."** Refuted
  by the sign data: `evalConic(Ecut)@V = [0, -0.046, -0.015]`, so the vertex `(-2,2)`
  sits on the g2-face-1 (discard) side and the far ray is legitimately g2-face-2. The
  clip is correct. The geometric measurements in that write-up are still good; the
  conclusion is not.
- **"The missing g1f6 ∩ g2f2 mirror is the defect."** The sharper diagnosis that
  replaced the above, itself now superseded — the mirror (`src[6 2]`) is built
  correctly.
- **"`decideWinner` wrongly declares an unbounded cell decided" (for `[2 -0.5]`).**
  The real cause was `assignSide` reading the winner at `piece.V(2,:)`, which for
  `splitCell`'s unbounded "rest" piece (`curveAfter=1`) is a **crossing point** where
  `diff ≈ 0` and its sign is noise, so both halves could get the same winner. Fixed
  in `53fc9fd` by reading the vertex farthest from `{diff=0}`. The parabola-to-infinity
  write-up was the wrong lead for this fixture.
- **"Piece 16's ray is piece 5's partner."** It is not — piece 16's ray lies on
  `x+y=0.5`, a parallel but different line from piece 5's `x+y=0`. They were never
  meant to match.

## 2026-08-03 — Retired hypotheses about the orphan ray — do not re-try

- **Tried / Why each failed:**
    - *The orphan ray is one physical ray covered to different extents* — nothing
      lies on it.
    - *The cut must be restricted to the arc's own span* — it must not; the straight
      half-planes are applied first.
    - *`clipPolyByConic` emits clockwise cells* — a guard on every bounded output
      never fires.
    - *`polyL` is non-convex* — all four curved faces measured convex.
    - *Pair `(6,1)` is spurious* — it is a genuine thin cell. The "zero intersection"
      reading was a resolution artifact: a 0.0625-spaced grid cannot see a cell of
      area 0.008.
- **Before retrying, fix:** `Do not retry.` Each was measured, not argued.

## 2026-08-24 — the exact integer layer has a PRECISION BUDGET, and it is the resultant that spends it

Measured while building `conicMeet`, the vertex layer of the sym-free conjugate. Recorded because
it is the one place the "faces and conics exactly, vertices approximately" design has a real cost,
and because the number is worth knowing before someone re-derives it.

**The budget.** `ratQ` is exact to 2^53. The Sylvester resultant of two conics is degree 2 in each
conic's coefficient vector, so a pair of conics with entries of size `M` produces a quartic with
entries of size about `M^2` -- and `M ~ 1e8` already exceeds the budget. Measured exactly: two
circles expressing a `1e-4` near-tangency need `M = 1e8` to be written exactly, and the resultant
reaches `2.3995600239996e17`, which raises `ratQ:overflow`.

**That is the designed outcome, and the alternative is the refuted one.** The only other thing the
layer could do is round the resultant, whose roots are then a plausible wrong answer -- which is
precisely the `double`-plus-tolerance design `DECISIONS.md` 2026-08-17 already refuted at a cost of
Step 3's cell count growing without bound. Raising is the correct behaviour and
`conicMeetTest/aTighterNearTangencyOVERFLOWSLOUDLYRatherThanReturningAWrongAnswer` pins it.

**What it does and does not bound.** It bounds how NEAR a degeneracy can be expressed exactly, not
the size of the problem: in the pipeline the conic coefficients are differences of face functions,
so their size is set by the input data, not by an arbitrarily fine perturbation. If it fires on
real input the answer is a wider integer type (int64 with checked arithmetic, or a small bignum),
never a tolerance.

**Two implementation facts that came out of the same measurement, both load-bearing.**

1. **The gcd chain overflows before the resultant does, and content removal fixes it.** The
   squarefree step divides the quartic by `gcd(p, p')` so that repeated roots become simple. A
   plain pseudo-remainder sequence squares the leading coefficients at every reduction: on the
   random sweep a quartic with six-digit entries reached `2.0006503375508298e23` in three steps.
   Taking the primitive part INSIDE the reduction loop -- not only at the end of each
   pseudo-remainder -- keeps the whole chain inside 2^53, and it is free because a gcd is
   unchanged by a constant factor.
2. **A repeated root is the normal case, not an exotic one, and it silently loses intersections.**
   `x^2 + 4y^2 = 4` against `x^2 = y^2` has resultant `(5x^2-4)^2`; `xy = 1` against
   `x^2 + y^2 = 5` has `(x^2-3)^2` after the shear. `roots` returns a double root to about
   `sqrt(eps)`, so the candidate misses the second conic by ~1e-7 and was REJECTED -- two of four
   intersections lost, with no error raised. Two changes were needed and both are kept: the exact
   squarefree reduction above, and polishing each candidate with 2D Newton on the pair BEFORE the
   acceptance test rather than after it. Correctness does not depend on the squarefree step (it is
   wrapped so that an overflow there falls back to the unreduced polynomial); accuracy does.

## 2026-08-24 — a test's `syms x y real` poisons the WHOLE MATLAB process, and the fast bucket shares one

Found while adding `ratQTest`. Recorded because the symptom points nowhere near the cause and the
next person to write a symbolic cross-check in a test will hit it.

**What happened.** `ratQTest` cross-checks `ratQ.diffConic` against the Symbolic Toolbox and opened
with `syms x y real`. Four tests in `regionTest` then failed with

    Unable to convert expression containing symbolic variables into double array.

They pass when `regionTest` is run alone, and fail when `ratQTest` runs first in the same MATLAB.
The fast bucket runs its whole bucket in ONE process deliberately (`.claude/suite.sh`: 26 MATLAB
startups is minutes of pure overhead), so adding one suite made four unrelated tests red.

**Why.** `syms x y real` does two things, and the second is not scoped to the function: it creates
the variables in the caller's workspace, AND it attaches a `real` ASSUMPTION to the symbols named
`x` and `y` in the shared MuPAD session. `region.m` builds its own `sym('x')` / `sym('y')`, which
are the same symbols, and the stray assumption changes what its expressions evaluate to.

**The rule.** In a test that needs symbols, use UNIQUELY NAMED ones and no assumptions --
`x = sym('ratQTest_x')` -- rather than `syms`. A name no production file uses cannot collide, and
without an assumption there is nothing to leak. `assume(x,'clear')` in a teardown would also work
and is strictly worse: it has to run, and a test that errors early does not reach it.

**Not a reason to give the suite its own process.** One process per suite is what the slow bucket
does and it costs ~8 s of startup each; the fix here is one line in the offending test.

## 2026-08-24 — the arc-bulge refusal in maxQuaPar is LIFTABLE, and the remaining blocker is ASSEMBLY, not representation

Worked most of the way; the item is scope-reduced rather than closed, and this entry exists so the
next attempt starts from the design instead of from the symptom.

**The claim that was wrong.** Two routines refused the same configuration -- a straight clip line
cutting a cell's parabolic edge TWICE, both arc endpoints on the same side -- on the grounds that
the result "is either disconnected or bounded by two separate arcs, neither of which is
representable as one QuaPar face". `insertPassthroughVertices` refused the sibling case (a
neighbour's vertex in the open interior of an arc) with the same reasoning, adding that lifting it
"means generalizing a piece to carry several arcs (curveAfter becoming a set)".

**It does not.** One arc per face is an invariant maintained by SUBDIVIDING, exactly as `splitCell`
already does for its own reason, and there is a canonical cut that always works:

> Cut the cell along the line through the split point PARALLEL TO THE PARABOLA'S AXIS.
> Such a line meets the conic EXACTLY ONCE -- that is precisely what `parabolaArcFrame`'s header
> already states makes `u` a global monotone parameter -- so it divides the arc into two sub-arcs,
> one per half, and creates no second arc anywhere. Both halves carry the same face function, so
> the FUNCTION is unchanged and only the subdivision is finer.

For the bulge case the cut point is forced, not chosen: along the arc, `nrm*x'-c` is the quadratic
`A2 u^2 + A1 u + A0` (`parabolaArcFrame.lineCoeffs`), so it has ONE stationary point
`u* = -A1/(2 A2)` -- the quantity `arcBulgesAcross` already computes to DETECT the case. Each
sub-arc lies on one side of it, so `nrm*x'-c` is monotone along each and the clip line crosses each
at most once. The recursion is therefore one level deep by construction: clipping either half by
the split line itself restricts to `u - u*`, which is affine, so `arcBulgesAcross` returns false.

**Implemented and green.** `clipPolyHalfPlane` returns a LIST, `bulgeSplit`/`splitAtArcU` do the
cut, `clipByFace` carries the list through, and `insertPassthroughVertices` splits instead of
raising. `maxQuaParTest` 29/0 and `addQuaParTest` 4/0 with it in.

**What still fails, and it is a DIFFERENT problem.** Both fixtures that used to hit the refusal now
get two stages further and die in `assemblePieces`:

    a boundary edge of piece 1 src [1 1] has no matching neighbour:
    (1.000000,1.000000)->(0.000000,0.000000), curved=1.
    Closest candidate: piece 4 src [2 4] (0.000000,0.000000)->(0.500000,0.500000) curved=0

So one side of a shared boundary carries the edge as CURVED and the other as STRAIGHT, and
`matchHalfEdges` pairs curved with curved. That is not a representation limit -- it is the
subdivision being inconsistent across the shared arc, which is the thing `insertGlobalPassthrough`
exists to fix and evidently does not yet fix for this case. **Do not re-attack the refusal; attack
the matching.** The likely shape of the answer is a global pass that collects, per CONIC, every
`u` at which any piece needs a split, and applies all of them to every piece carrying that conic --
so the subdivision is consistent by construction rather than per-pair.

**Behaviour is unchanged for callers.** `maxQuaPar:internal` still starts with `maxQuaPar:`, which
is what `conjCPLQ`'s fallback catch tests, so under `route='auto'` these inputs still fall through
to the symbolic Case C exactly as before. The change is kept because the primitives are correct and
tested and the blocker is now named; it does not by itself remove a fallback.

## 2026-08-24 — the numeric Step 3 can DROP a cell on an unbounded fold, and the fix is a cross-check, not a narrower gate

The most important finding of the session, because the failure is silent.

**What happened.** Removing `conjCPLQ`'s `isDomBounded` gate let a 4-cone fan take the numeric
route. At `s = (-2,-3)` it returned **2.0** where the definition sup is **4.5**. The assembled
result had **4** cells for a fold of four 4-cell conjugates: the cell carrying one face's strip
had been dropped. Every probe point of the OTHER orientation of the same fixture was exact -- which
is the signature of a dropped REGION rather than an arithmetic error: right almost everywhere,
wrong on a set, and silent.

It was caught by `conjCPLQTest/step3UnboundedAssemblyMatchesTheTruth`, whose closed-form truth is
independent of the pipeline. Two other tests in the same file failed only because they read the
result with `evalFunctionNDomain(g.fnd, ...)`, which a meshed result does not have -- a ROUTE
expectation, not a value one, and the file's own `evalConjResult` exists for exactly that.

**The fix is a verification, not a restriction.** `f* = max_k (q_k + I_{P_k})*` is an IDENTITY, and
the per-face conjugates that were folded are still in hand, so their pointwise max is f* exactly --
no reference implementation, no quadrature, no tolerance on the mathematics. `conjPolygonalDomain`
now checks the assembled result against that max at a spread of directions and magnitudes, at the
result's own vertices, and at the dual points of the input's vertices; on a disagreement it raises
a `maxQuaPar:`-prefixed error, which is what `conjCPLQ`'s fallback catch tests, so the domain routes
to the symbolic Case C instead of returning a number nobody checked.

The check is ONE-SIDED and is documented as such: sampling can miss a defect, it cannot invent one.

**Why not simply put the gate back.** The gate refused an entire shape family -- every unbounded
domain -- to avoid a defect that affects some of them. With the cross-check, the ones that assemble
correctly are answered in 0.01-0.1 s and the ones that do not fall back exactly as before. That is
strictly better than the gate, and it also now protects the BOUNDED path, which had no such check.

**What it does not do.** It does not fix the drop. `maxQuaPar`'s unbounded fold losing a cell is a
real defect with a reproducer (`conjCPLQTest/step3UnboundedAssemblyMatchesTheTruth`'s fixture, the
`F = [3 2;2 1;1 4;4 3]` orientation), and `PARTITION OVERLAP` diagnostics fire on it. That is the
next thing to fix, and it is separate from the arc-split work above.

## 2026-08-24 (last) — G2b FIXED: the split direction at a line pair's SINGULAR POINT

The dropped-cell defect recorded above is closed, and the cause is one line.

**The minimal reproducer is two functions of one variable each.**

    h1 = max(s1,0)^2/2      as two cells split by the s2-axis
    h2 = max(s2,0)^2/2      as two cells split by the s1-axis

Their max must split the FIRST QUADRANT along `s2 = s1` -- h1 wins below it, h2 above -- so the
result has 5 cells. `maxQuaPar` returned **4**, with the whole first quadrant carrying h1, and
`max(h1,h2)(2,3) = 2.0` against a truth of `4.5`. Confirmed pre-existing by running the same
input against `maxQuaPar.m` as of `8ea857b`.

**The chain, traced rather than guessed.** `{h1 = h2}` is `s1^2 = s2^2`: a degenerate conic, and a
PAIR of lines crossing at the origin -- which is the quadrant cone's apex and its only vertex.
`decideWinner` correctly refused to decide (the two rays give opposite asymptotic signs).
`splitCell` then found ONE boundary crossing, at that apex, and handed it to
`splitUnboundedAtOneCrossing`, which takes the branch direction as the perpendicular to the
GRADIENT of the difference -- and the gradient VANISHES at a line pair's crossing point. It
returned "no branch here"; the caller's tangency branch then resolved the winner at the cell's
centroid, which for a cone IS the apex, sitting exactly ON the curve. The winner came out of the
sign of a computed zero.

**The fix.** Where the gradient vanishes, `diffRow(X + t d) = (t^2/2) d' Q d`, so the branches
through `X` are exactly the NULL DIRECTIONS of the quadratic part `Q`. Those come off `Q`'s
eigen-decomposition in closed form -- `d = sqrt(-l2) u1 +- sqrt(l1) u2` -- with no root-finding and
no case analysis on which coefficient happens to vanish (`nullDirectionsOf`). Both candidates are
then put through the existing recession test, which picks the one that actually cuts the cone.

**Why it mattered so much.** A cone whose apex is the singular point is the dual of any face that
is a WEDGE, so every 4-cone conjugate hit it. It is the reason `conjCPLQ` had to gate unbounded
input away from the numeric route at all, and with it fixed the 4-cone fan assembles exactly
against the definition sup (8 cells, error 0 at every probe point) instead of returning 2.0 where
the truth is 4.5.

Pinned by `maxQuaParTest/twoHalfPlaneQuadraticsSplitTheirSharedQUADRANT` (the two-line reproducer,
asserting the cell COUNT as well as the values, because 4-vs-5 is the defect) and
`aFourConeFanFoldsExactlyAgainstItsOwnPieces`.

**The cross-check stays.** It cost nothing, it is what caught this, and it is one-sided.

## 2026-08-24 (later) — G1 is a LENS, not a matching bug: two operands' boundaries share endpoints but not the curve

Sharpening the earlier G1 entry with what the piece dump actually shows, so the next attempt does
not start from the error message.

**The symptom.** `f = x^2+y^2` on one triangle of the unit square and `f = xy` on the other. The
numeric route dies in `assemblePieces`:

    a boundary edge of piece 1 src [1 1] has no matching neighbour:
    (1,1)->(0,0), curved=1.  Closest candidate: piece 4 src [2 4] (0,0)->(0.5,0.5) curved=0

**What the pieces show.** Piece 1 is `V = [(2,0); (1,1); (0,0)]` with `curveAfter = 2`, so its arc
IS the edge `(1,1)->(0,0)` -- g2's parabola `(x+y)^2/4 = x`. The pieces offered as neighbours,
4 and 5, are `[(0,0),(0.5,0.5)]` and `[(0.5,0.5),(1,1)]`, both STRAIGHT, both `src [2 4]`. And
`(0.5,0.5)` is **not on the parabola**: `(0.5+0.5)^2/4 = 0.25`, not `0.5`.

**A CELL IS MISSING, and this is measured rather than argued.** A hand-rolled containment test was
tried first and is not trustworthy -- written against the chord it misclassifies exactly this
region, and rewritten against the arc it excluded a point that IS in piece 1. The test that settles
it uses each operand's OWN point location, captured in the same session as the pieces so the face
numbering cannot drift:

    g1 nf=7  g2 nf=6
    srcs produced (15): [1 1;1 2;1 6;2 4;3 2;3 5;3 6;4 5;5 1;5 3;5 4;5 6;5 7;6 4;6 7]

      (0.20,0.40) g1 face 4, g2 face 2 -> src present: 0
      (0.30,0.50) g1 face 4, g2 face 2 -> src present: 0
      (0.15,0.30) g1 face 4, g2 face 2 -> src present: 0
      (0.45,0.55) g1 face 4, g2 face 2 -> src present: 0
      (0.60,0.40) g1 face 1, g2 face 2 -> src present: 1     <- control, below the diagonal

Every point of the lens lies in g1's face 4 and g2's face 2, and the fold produced **no piece with
src [4 2]** -- while the control point one step away, in [1 2], is there. So
`clipByFace(g1.face4, g2.face2)` returned nothing for a pair whose intersection is not empty. The
orphaned arc in `assemblePieces` is the SYMPTOM; the missing cell is the defect.

**Where it is NOT.** `polyConstraints` already skips a curved edge's chord, and says why in its own
comment -- so the lens is not being clipped away by the chord half-plane. The remaining candidates,
in order: the operand SWAP at the top of `clipByFace` (which of the two carries the arc decides
whether it becomes the cell's own arc or the cutting conic); `clipPolyByConic`; and the three
reduction passes `dropDegeneratePieces` / `dedupPieces` / `dropSubsumedPieces`, any of which could
discard a thin sliver. Instrument `clipByFace` for that ONE pair rather than reading it.

**Why the arc-split work does not close this.** That work makes ONE cell's arc divisible; this needs
the OVERLAY of the two operands' boundaries to be complete, which is a different property. The
"consistent per pair rather than globally" phrasing in the earlier entry is right in outcome and
wrong in cause -- the missing thing is a cell, not an alignment.

**Related, and worth reading first.** `SUPPORT_MATRIX.md` 4.1 already lists "half-edges and boundary walks identifying
an edge by its endpoints alone, which four arcs between the same two points make ambiguous" among
the defects fixed on 2026-08-13. This is the same family and it is not fully closed: a chord and an
arc with the same endpoints are two different edges. If the coverage measurement confirms a hole,
the fix is to build that face when `clipByFace` produces a cell whose straight edge and the other
operand's arc share both endpoints, rather than letting one stand for both.

## 2026-08-25 — a RANDOMIZED definition check finds two PRE-EXISTING defects the fixtures missed

`checkConjAgainstDefinition` runs `conj` on random convex polygons carrying a random convex,
indefinite, concave or affine quadratic, and compares against `sup_{x in P} <s,x> - q(x)` computed
by a scan plus a pattern search. 24 cases, seed 20260824.

**Twenty-two of twenty-four are EXACT** -- `worst = +0.000e+00`, or `2e-16` on two convex cases.
That covers every sign class on 3-, 4- and 5-gons, and it is the broadest evidence this project has
that `conj` is right rather than merely green on the shapes someone thought to write down.

**Two are not, and BOTH ARE PRE-EXISTING.** Re-run against a pristine snapshot of `b9243d3` -- the
tree as it stood before any of 2026-08-24's work -- they reproduce with the same magnitude to every
digit printed:

    case 21  3-gon, indefinite xy   conj is 2.742e-02 BELOW the sup      (a MINORANT)
    case 17  5-gon, indefinite xy   MATLAB:badsubscript                  (a crash)

    W21 = [0.6057047151 0.9300751811; -0.3353947472 0.5251524293; -1.082499617 0.08448609744]
    f21 = [0 1 0 -0.7177913413 -0.6075645347 -0.6781835233]

    W17 = [0.9180323951 0.1778978365; 1.189169914 0.308421002; 0.407632005 1.314428266;
           -1.091935477 -0.4741970304; -0.8285840095 -0.5854877405]
    f17 = [0 1 0 1.03766872 0.8199895741 -0.2736502742]

Both are `xy` on a polygon with a general affine part, i.e. the A.3/A.4/A.5 machinery, which the
2026-08-24 work did not touch. Case 21 is a single TRIANGLE, so it has no cross-piece max and no
fold: the minorant comes out of the per-piece closed form itself. That is the same SHAPE of defect
as `DECISIONS.md` 2026-08-19's "`plq_1p`'s A.4 branch computes a MINORANT, not the envelope",
recorded there for the symbolic path.

**The other half of the differential is the good news.** Every other number is identical between the
two trees, and the five CONVEX cases changed exactly as intended -- `QuaPar` became `QuaPol` and
they got 2-3x faster (1.51 -> 0.40 s, 1.32 -> 0.27 s, 0.69 -> 0.32 s, 0.67 -> 0.20 s,
0.57 -> 0.23 s). So the night's changes altered the representation and the route, and altered no
value anywhere in the family.

**A cheap check that would have caught case 21, not yet built.** For any vertex `v` of `dom f`,
`f*(s) >= <s,v> - f(v)` -- a one-line lower bound, valid for every route including the single-piece
one, which has no "max of pieces" identity to test against. A conjugate that dips below the max of
those affine minorants is definitely wrong. Before wiring it as a refusal, MEASURE how many existing
fixtures violate it: a check that reroutes half the suite to the symbolic path is not an improvement.

## 2026-08-25 (later) — the VERTEX bound catches nothing; the EDGE bound catches exactly the defect

G6 was proposed as `f*(s) >= <s,v> - f(v)` at every vertex of `dom f` -- a one-line, always-valid
lower bound that would work on every route including the single-piece one, which has no
"max of pieces" identity to check against. **Measured before building it, and the proposal as
stated is useless:**

    vertex bound, 24 random cases:   0 violate    (worst -8.9e-16, machine noise -- including
                                                   case 21, the known 2.7e-2 minorant)

The reason is immediate once seen: case 21's sup is not attained at a vertex, so the vertex bound
is slack exactly where the answer is wrong.

**The EDGE bound is the right strengthening, and it is just as cheap.** Along the segment
`v_i + t(v_j - v_i)`, the objective `<s,x> - q(x)` is a QUADRATIC in `t`, so its maximum over
`t in [0,1]` is one closed-form expression -- two endpoint values plus, when the quadratic is
concave, the interior stationary point. It is a valid lower bound on `f*` for the same reason a
vertex is: the sup over the domain is at least the sup over any subset. Measured on the same 24:

    vertex + edge bound:             1 violates   -- case 21, at exactly -2.742e-02,
                                                     and every other case at ~1e-15

So it is SHARP on this family: it fires on the one known wrong answer and on nothing else. That is
the property a check has to have before it is wired as a refusal -- one that reroutes correct
answers to the symbolic path is not an improvement.

**Cost.** One quadratic per edge per probe point, no iteration, no engine. Cheaper than the fold
cross-check already in `conjPolygonalDomain`, and unlike that one it covers the SINGLE-PIECE route.

**What it still cannot see.** A sup attained strictly inside the domain. For a CONCAVE or
INDEFINITE piece the maximiser of `<s,x> - q(x)` is on the boundary, so the bound is tight there;
for a CONVEX piece it can be interior -- but that is the route `conjConvexPolygon` handles in
closed form and which the random sweep found exact in every case. Note the limit rather than
claiming the check is complete.

## 2026-08-25 (last) — `--verylong -j N` RACES on the stage cache, and it produced a spurious red

The daily gate was run twice on the current tree and once on a pristine snapshot of `b9243d3`, and
the three results settle both questions the night raised about it.

    baseline b9243d3, -j 2, uncontended     26 pass /  7 fail / 1 timeout
    current tree,     -j 2, uncontended     26 pass /  7 fail / 1 timeout    <- identical
    current tree,     -j 2, contended       25 pass /  8 fail / 1 timeout    <- one extra

**So the night's work changed the verylong bucket not at all.** The seven failures and the timeout
are pre-existing, `testMaxMultiRegion/testPCE2` among them (the handoff already listed it as the
known red).

**The eighth is a RACE, not a regression.** `testcPLQ/rectConvexEnvelopeUnderestimates` failed only
in the contended run, comes back `fail=1 incomplete=1` (i.e. it ERRORED rather than mis-asserted),
and **passes when run alone**. The mechanism is in the tooling, not the mathematics:

* `plqStage` caches each pipeline stage under `.claude/stagecache/`, keyed by (fixture, stage), and
  `get` LOADS that file when it is fresh;
* `suite.sh --verylong -j N` runs N MATLAB processes against that ONE directory;
* the job list is per TEST, and consecutive tests of one fixture are consecutive STAGES -- so with
  `-j 2`, `rectTriangulates` (stage 1) and `rectConvexEnvelopeUnderestimates` (stage 2) run at the
  same time. A MISSING cache is safe, because `get` then recomputes; a HALF-WRITTEN one is not,
  because `load` on a partial `.mat` throws.

The timing fits: stage 1 passed at 28 s and stage 2 errored at 29 s.

**The fix is an atomic write**, and it is three lines: `save` to a unique temporary name in the same
directory and then `movefile` onto the real one, so a reader sees either the old file or the new one
and never a partial. A `load` that throws should also fall back to recomputing rather than failing
the test. NOT DONE HERE: verifying it needs another contended `--verylong` run, and an unverified
change to a caching layer is exactly what should not be committed unattended.

**Until it is fixed, read a `--verylong -j N` red against a `-j 1` re-run of that one test before
believing it.** A spurious red in the daily gate costs a morning.

## 2026-08-25 (final) — the EDGE lower bound is now a DEFAULT refusal, and it is fully gated

Built, measured on every bucket, and turned on. `conjEdgeLowerBound(q, S)` returns
`max over the boundary of dom f of [<s,x> - f(x)]` -- one closed-form quadratic maximisation per
edge, no iteration, no engine -- and `conjCPLQ` raises `PLQ:conjCPLQ:belowEdgeBound` when the
conjugate it is about to return falls below it.

**Why it is allowed to be the default.** Three measurements, in the order they were taken:

    24 random polygons + quadratics    fires on exactly 1, the case wrong by 2.7e-2
    every FAST and SLOW suite, on      363 pass / 0 fail -- flags no correct answer
    fast / slow / verylong, on         298/0 · 88/0 · 26 pass 7 fail 1 timeout
                                       -- the verylong figure IDENTICAL to a pristine b9243d3

So it converts a class of silent wrong answer into a named refusal and touches nothing that works.

**Why it RAISES rather than falling back**, which was the one design question worth measuring: on
the known-bad triangle the SYMBOLIC route returns the same wrong value to six digits
(2.997553 against a bound of 3.013340). The defect is in the shared Step 1 / Step 2 closed form,
not in a route, so there is nothing to fall back to -- and a fallback that silently produced the
same minorant would be worse than an error.

**What it covers that nothing else did.** `conjPolygonalDomain`'s fold cross-check verifies an
assembled result against `f* = max_k (q_k + I_P_k)*`, which needs at least TWO pieces. A single
triangle has no such identity, and that is exactly where G4 lives.

**What it cannot see, stated rather than glossed.** A sup attained strictly INSIDE the domain. For
a concave or indefinite piece the maximiser is on the boundary and the bound is tight -- pinned by
`conjEdgeLowerBoundTest/theBoundIsTIGHTWhereTheSupIsOnTheBOUNDARY`; for a convex piece it can be
interior and the bound is genuinely slack, pinned as an inequality plus a witness by
`theBoundIsVALIDButSLACKWhereTheSupIsINTERIOR`.

**The consequence for G4.** `conj` of `xy` over some triangles now RAISES where it used to return a
minorant. That is the intended trade and the defect itself is still open.

## 2026-08-25 — REFUTED: branching to isolate one unattended run from another in the SAME working tree

- **Tried:** the overnight run created `overnight/2026-08-24` so its work would be separable from a
  parallel `proof/` session committing to `main`.
- **Why it failed:** two sessions in one repository share one working tree and therefore ONE HEAD.
  Branching did not isolate this run -- it silently MOVED the other session onto this run's branch,
  and both runs' commits landed there. It also made a `git add -A` here sweep two of that session's
  in-progress files into a commit, and it left `main` looking stale while all the work sat
  elsewhere.
- **Before retrying:** `Do not retry.` Branching cannot isolate sessions that share a checkout.
  **Stay on `main`** -- which is the standing preference for this repo anyway, and which makes the
  problem disappear rather than solving it: with both sessions on `main` there is no branch for the
  other one to be dragged onto, and interleaved linear history in a solo repo is not a defect. The
  branch was folded back in by fast-forward with nothing lost.
- **CORRECTED the same day.** This entry first recommended `git worktree add`, and that is WRONG
  for this setup -- it was generalised from the incident without being checked. Measured:
  `git worktree add <path> main` fails with `fatal: 'main' is already used by worktree at ...`, so
  worktrees FORCE two different branches, which is the very thing that caused the tangle and the
  thing the commit-on-`main` preference exists to avoid. Worktrees earn their keep when you
  genuinely want different branches checked out at once (`arch/co1d` runs six method-paths that
  way); that is a different problem from two sessions editing one project.
- **The other two symptoms had other causes**, and neither is fixed by any branching scheme: the
  `git add -A` sweep is fixed by staging explicit paths, and the `suite.sh`-edited-mid-run incident
  was one session editing a file its own `bash` was executing.
- **Evidence:** `proof/MORNING.md` reaches the same conclusion independently from the other side of
  the incident; this run's `MORNING.md` opens with it. A third symptom of the same cause: editing
  `.claude/suite.sh` while a `bash` process was executing it killed a slow-bucket run mid-flight
  with a bogus syntax error.

## 2026-08-25 (split) — the Lean proof is its own repository; the shared-working-tree problem is gone

- **Closes** the entry above (`REFUTED: branching to isolate one unattended run from another in the
  SAME working tree`). That entry ends at "stay on `main` and accept interleaved history", which is
  a way of living with the collision rather than removing it. Removing it is cheaper: the two
  sessions were only ever entangled because the Lean formalisation sat at `proof/` INSIDE this
  repository, so a `proof/` session and a MATLAB session shared one HEAD, one index and one
  `git add -A`.
- **Done:** `proof/` moved verbatim to `AI/CCA2proof` and given its own repo (initial commit
  `e629247`, the same 32 tracked files, `.lake/` still ignored). This repository keeps the
  pre-split history of `proof/`; the new one starts fresh. No history was rewritten here.
- **What this changes in practice:** two sessions, two working trees, two HEADs -- so neither the
  branch tangle, the cross-session `git add -A` sweep, nor a stale-looking `main` can recur between
  MATLAB work and proof work. The commit-on-`main` preference still stands inside EACH repo.
- **What it does not change:** references to `proof/...` in this file and in `MORNING.md` are the
  historical record of when the two lived together and are left as written. Only the live pointer
  in `.claude/SESSION_HANDOFF.md` was repointed.

## 2026-08-25 (G5) — the frame-change branch collapsed the conjugate's BLOCK INDEX, and the crash was standing in front of an intractable symbolic max

- **Found:** `MATLAB:badsubscript` in `plq_1p.maximumConjugate`, on `checkConjAgainstDefinition`
  (seed 20260824) case 29 and case 17. `TODO.md` recorded it as "SOME indefinite 5-gons"; case 29
  is a **4-gon**, so the shape class in that note is wrong and a fix aimed at 5-gons would have
  missed half of it.
- **The defect.** `conjugate`'s FRAME-CHANGE branch redoes an indefinite non-`x*y` piece in the
  z-frame and reads the answer back. It copied the z-frame object's ENVELOPE -- two faces on this
  input -- while replacing its `conjfia`, the per-face block boundaries into `conjugates`, with
  the single block `[1 nConj+1]`. Measured, both triangles of case 29:

        envelope faces = 2, numel(conjfia) = 2, nConjugates = 11

  `maximumConjugate` sizes its loop from `size(envelope,2)` and indexes `conjfia(i+1)`, so face 2
  asks for `conjfia(3)` of a 2-element array. The crash is the visible half; the other half is
  that had it not crashed, the max ACROSS the z-frame's faces would never have been taken and the
  piece would have returned a minorant.
- **Fixed, in two places.** The branch carries `objT.conjfia` (`substituteFrame` is a per-cell
  substitution, so order and count are preserved and objT's boundaries are the right ones), and
  `maximumConjugate` now takes its block count from `numel(conjfia)-1` rather than from the
  envelope. The two arrays are NOT interchangeable and only the ordinary path keeps them in step:
  the SEPARABLE branch legitimately returns the whole conjugate as one block -- a mesh that
  partitions the dual plane, already maximal -- over an envelope with several faces.
- **REFUTED, same session: taking the max in the z-frame first.** The argument is sound -- the max
  commutes with the read-back because `s -> M's - a` is an invertible affine map, and the z-frame
  coefficients are rational where the substituted ones are surds over forty-digit integers. It
  makes no difference: neither order finishes in **25 minutes** on case 29, because `maximumP` on
  this fixture's two envelope faces is itself the intractable step. Do not retry; take the simpler
  route and leave `conjugate`'s cost where it was (38 s, measured).
- **What that means.** G5's crash was in front of a symbolic max that does not terminate in
  practical time, so case 29 goes from a fast wrong CRASH to a correct computation nobody can wait
  for. The answer for that family is the NUMERIC route, which declines it today with
  `maxQuaPar:notImplemented` (clipPolyByConic separating an unbounded cell) -- that gap, not this
  branch, is what would make it fast.
- **Pinned by** `conjCPLQTest/frameChangedPieceKeepsItsEnvelopeBLOCKS`, which asserts the
  INVARIANT (one block per envelope face plus a terminator) rather than the value, because a value
  assertion here is a test nobody can run.

## 2026-08-25 (G4) — the minorant is in the FOLD, not the closed form, and it is 4x, not 2.7e-02

- **`TODO.md` said** "a single triangle, so there is no fold and no cross-piece max -- the
  per-piece closed form itself is short". **Measured: false on both counts.**
- **There IS a fold.** Step 1 splits sweep case 21's triangle into **four** faces -- the nCE==3
  cevian split with each half re-split -- two of them slivers of area 8.7e-05 against 2.7e-02.
  Every face's own Step 2 conjugate is EXACT at the bad point: 1.032507658472 to twelve digits,
  four times over, against the closed-form bilinear sup. Folding faces 2 and 3 keeps it; the
  **fourth fold** returns 1.005089907622. So the assembly is the defect and the closed form is not.
- **Pairwise is fine, accumulated is not.** `max(face1,face4)`, `max(face1,face2)`,
  `max(face1,face3)` and `max(face4,face1)` are each exact. Only folding a sliver into the
  ALREADY-ACCUMULATED mesh loses the value, which is why every earlier two-operand `maxQuaPar`
  test passed.
- **It is far worse than the recorded number.** The 2.742e-02 deficit was the worst on a probe grid
  of radius <= 6. The fold-vs-pieces cross-check finds `f*(-10,0) = 47.10181578` against a true
  10.86895777 -- an **over-estimate by a factor of four**. The G6 edge lower bound cannot see it:
  that bound is one-sided and catches minorants only.
- **With `MAXQP_ASSERT = 2`** the fold raises two invariant violations on this input, both real and
  neither fixed here: a piece whose declared rays are the NEGATIVE of the direction its constraint
  region recedes along, and a piece that carries one operand's quadratic where the other is larger
  by `Inf` along a ray. The second is the winner-domination failure that produces the wrong value.
- **What landed:** `conjMaxOfSubTriangles` now cross-checks the fold against the pieces it was
  built from, exactly as `conjPolygonalDomain` has since the 4-cone fan incident. The identity
  `f* = max_k (q_k + I_{T_k})*` makes the pieces their own oracle, so this needs no reference
  implementation.
- **It REFUSES rather than falling back**, and that is the same trade G6 made. A `maxQuaPar:`
  identifier would route the triangle to the symbolic Case C; measured, that did not finish in 25
  minutes, against **2.5 s** for the refusal, and nothing says Case C's answer on it is better. The
  new identifier is `PLQ:conjCPLQ:foldDroppedACell`.
- **Still open:** the `maxQuaPar` defect itself. What is closed is the misdiagnosis -- and the
  silence.

## 2026-08-25 (G10) — the accumulated fold: one defect FIXED, one REFUTED fix, and the real cause named

The `maxQuaPar` half of G4. Three things were measured; two of them are now in the code and one is
a refutation worth more than the attempt.

**FIXED — `pieceRecessionRays` decided the recession cone from floating-point noise.** Its header
claims "all sign tests are exact (sym over the pieces' rational coefficients)", and that is the
flaw: the coefficients arrive as DOUBLES and are lifted with `sym()`, so the comparisons are exact
on numbers whose last bits are rounding. Two of them compare against a mathematical ZERO:

    A(d) = d*Q*d'   is 0 along a parabola's null direction -- which for a half-strip piece is the
                    direction being tested
    cross(E,d)      is 0 when d lies along a ray edge -- which for that piece it does

Measured on G4's fixture, on a piece with `nv=2, curveAfter=1, unbounded` whose `dirIn` and
`dirOut` are the SAME direction stored as two different doubles (they differ in the last bit):

    A(dirIn)  = +4.37e-20     passes `A < 0` outright, and is not `== 0`, so the linear-term
                              TIE-BREAK -- the branch written for exactly this case -- never ran
    A(dirOut) = -2.49e-17     rejected by `A < 0`
    cross(-dirIn, +dirIn) = -5.55e-17   rejected the one direction every real constraint admits

The piece therefore reported receding along the exact NEGATIVE of its own declared ray, and the
whole decision was made by noise of order 1e-17. Both tests now compare against a tolerance scaled
by their operands. That widens the TIE case rather than weakening it: the tie is then settled by
the linear term, measured `|grad.d| = 8.1e-02`, four orders of magnitude clear of the tolerance.
This is not cosmetic -- `reccConeViolation` gates `halfIsProvedWellFormed`, so a false violation
makes a legitimate two-arc ray cut get refused.

**REFUTED — routing `splitCell`'s single-tangency branch through `assignSide`.** That branch reads
an unbounded cell's winner from the CENTROID of its finite vertices, which says nothing about which
operand wins far out; `assignSide` reads the recession cone instead and applies
`assertWinnerHoldsAtInfinity`. It fixes G4 -- `f*(-10,0)` stops answering 47.10181578 against a
true 10.86895777 -- and it BREAKS a case that was right: `conjSymFreeTest`'s A.5 three-convex-edge
triangle {(5/2,3/2),(0,0),(1/2,1)} carrying `x*y` matches the closed-form bilinear sup to 3.6e-15
over the whole probe grid under the centroid read, and REFUSES through `assignSide`. Its cell
straddles {f1=f2} as well, but the straddle never reaches the assembled answer. Trading a correct
answer for a refusal is not an improvement, so this is reverted, with the measurement kept at the
site. **Before retrying: do not.** The backstop is right that both cells straddle and wrong that
both are unusable, so a version of this that merely softens the backstop is the same idea again.

**The real cause, and where the fix belongs.** `splitCell` reaches that branch only when it found
ONE boundary crossing of the splitting curve. For an UNBOUNDED cell that is not a tangency at all:
the curve can enter through the boundary once and leave through the RECESSION CONE, so the cell has
one finite crossing and genuinely splits into two unbounded parts, each with its own winner.
`splitUnboundedAtOneCrossing` exists for exactly that and declined on both fixtures. Make it
succeed and neither cell needs a winner read off a centroid, and neither needs the backstop.

## 2026-08-25 (G10, second pass) — `splitUnboundedAtOneCrossing` conflated PARALLEL with ALONG, and declined a real split

- **Where the guess comes from.** `splitCell`'s tangency branch reads an unbounded cell's winner
  off the centroid of its finite vertices, which is the G10 wrong answer. It is reached only when
  `splitUnboundedAtOneCrossing` returns `[]`, so the question is not "how do we guess better" but
  "why did the split decline".
- **Measured, on `conjSymFreeTest`'s A.5 triangle {(5/2,3/2),(0,0),(1/2,1)} with `x*y`:**

        cell nv=2, dirIn = dirOut = [0.5145 -0.8575]      (a HALF-STRIP)
        |grad| at the crossing = 1.727e-01                 (not the singular-point case)
        candidate [ 0.5145 -0.8575]: max(n.w) = -5.551e-17  -> recedes the cell
        candidate [-0.5145  0.8575]: max(n.w) = +1.524e-01  -> does not

  The good candidate was then thrown away by `norm(w - dirIn) > 1e-7 && norm(w - dirOut) > 1e-7`,
  whose comment is "must not run along one of the cell's own rays, which would cut off nothing".
- **That test conflates PARALLEL with ALONG.** For a half-strip the escaping branch is parallel to
  both rays by construction, so EVERY candidate is rejected and the cell is reported as a tangency.
  A cut actually cuts off nothing only when the whole ray `{X + t w}` already lies in the boundary,
  i.e. when the crossing X sits ON that ray edge and w is that ray's direction. From any other
  boundary point the same direction slices the cell in two: a half-strip cut parallel to its rays
  is two half-strips. The test now asks that instead, using `hit.edge`.
- **Result.** The A.5 triangle splits and stays exact (3.552714e-15 against the closed-form
  bilinear sup over the whole probe grid, unchanged). fast 303/0/0.
- **What it does NOT fix, stated rather than glossed.** G4's own cell still assembles wrongly. With
  the definition checks off the numeric route now returns `Inf` at some probe points where it used
  to return a finite over-estimate -- the assembled mesh has a HOLE, so the split path has its own
  gap for this cell shape. Both are wrong and both are refused in production (`foldDroppedACell`,
  2.3 s), so nothing regressed that a caller can see; but the next step for G10 is that hole, not
  another winner heuristic.
- **`splitDeclined` now records the reason** under `MAXQP_ASSERT`. Every decline here ends in the
  centroid guess, and three separate investigations of that guess have begun by re-instrumenting
  this function.

## 2026-08-25 (B4/G3) — an unbounded non-convex face is TWO cases, and neither is "not implemented"

- **What the entry said:** "a non-convex face over an UNBOUNDED polygon declines by name today (the
  fan-triangulation route needs a BOUNDED face)". True of the NUMERIC route, and read ever since as
  though the case were unimplemented. Measured, it is two different things:
- **Envelope FINITE -> already answered, exactly.** `x*y` on the first quadrant is bounded below on
  its own recession cone, so conv f is finite. Case C returns the right function in ~14 s:

        f*(s) = sup_{x,y >= 0} s1 x + s2 y - x y  =  0 for s <= 0, +inf otherwise

  i.e. the INDICATOR of the third quadrant, checked at ten dual points against that closed form.
  So this half is a symbolic FALLBACK -- a cost, on the list with G1 and G2 -- and what it wants is
  a numeric route, not a construction.
- **Envelope -inf -> correctly refused.** `x^2 - y^2` on the same quadrant runs to -inf along the
  y-axis, which lies IN the recession cone, so conv f = -inf and f* = +inf everywhere.
  `convEnvUnbounded:minusInfinity` says exactly that, and it is the RIGHT answer: what blocks it is
  having nowhere to put an EMPTY domain, the same representation blocker as Case A's
  `conjugateHasEmptyDomain`. Not a missing construction.
- **The general rule**, worth having written down: for an unbounded face, conv f is finite iff f is
  bounded below on the face's recession cone; otherwise f* is +inf everywhere and dom f* is empty.
  That single test separates the two halves above and should be what any future numeric route
  dispatches on.
- Both halves are pinned by tests in `conjCPLQTest`; `SUPPORT_MATRIX.md` now carries them as two
  rows instead of one GAP.

## 2026-08-25 (G13) — `testPSqroot` is a LEGACY red; `conj` gets the same input right

- **Measured**, triangle {(-1,1),(-3,-3),(-4,-3)} carrying `x*y`, at `s = (-1,-1)`:

        truth (closed form; the sup is on the edge (-4,-3)->(-1,1) at t = 0.75)   1.75
        QuaPol.conj, numeric Case B                                              1.75    correct
        plq_1piece + plq.maximum, which is what testPSqroot runs                -5       WRONG

  `-5` is exactly the objective's value at the VERTEX (-4,-3), so the legacy pipeline takes a
  vertex where the sup is strictly inside an edge -- the same shape of defect as G4, in a different
  implementation.
- **It does not reach `conj`.** Forcing `conj` down the symbolic route on the same triangle raises
  `PLQ:conjCPLQ:cplqFailed` after 102 s rather than returning -5, because Case C never sees a
  single triangle (see below). So this red belongs to `plq_1piece`, i.e. T6's migration, and does
  not block the conjugate.

## 2026-08-25 (G12) — REFUTED: gating Case B on `~forceSymbolic`

- **Tried:** `route='symbolic'` is documented as "skips straight to Case C", and Case B ignored it
  -- the exact counterpart of the `'numeric'` hole closed earlier the same day. Gating Case B on
  `~forceSymbolic` looks like the obvious fix.
- **Why it failed:** Case C does not cover a single triangle. Its own header scopes it to `nf>1`
  and/or a non-triangular face, and sending the triangle above there raises
  `PLQ:conjCPLQ:cplqFailed` after 102 s. `'symbolic'` has nowhere to go for that shape, so the
  no-op is the honest behaviour and the note now sits at the branch.
- **What it did expose, and it is real:** `biconjCPLQ` asks for the symbolic form exactly when the
  numeric first conjugate `isMeshed && hasCurvedEdge`, because the second conjugation cannot take a
  curved mesh. For a single triangle that request is answered with the SAME curved mesh, so the
  escape does nothing. Filed as `TODO.md` G12. The fix is a symbolic form for a single triangle --
  Step 2's own cPLQ output -- not a reroute of the whole domain.

## 2026-08-25 (G11a) — `rectConjugateMatchesTheSup` was a FALSE red: the oracle was the coarse thing

- **One of the seven `verylong` reds, and it was never a defect.** The diagnostic reads
  `f*(0,0) = 2 EXCEEDS the sup 1.99996661` -- the conjugate exceeding the reference by 3.34e-05,
  asserted at 1e-05.
- **Measured against the closed form** on the piece in question, `testcPLQ`'s PRect piece 3, the
  triangle {(0,-4),(1,3),(1.8708,-0.2583)} carrying `x*y`:

        s = ( 0, 0)   closed form  2.0000000000    plq_1p  2               grid  1.9999666056
        s = (-1, 1)   closed form -0.4285714286    plq_1p -0.428571429     grid -0.4286032524

  The conjugate is exactly right to ten digits. The GRID is short by 3.3e-05.
- **The test was asserting the reference to an accuracy it does not have.**
  `plqCheck.supOverDomain` was a pure 200-per-edge / 900-interior sample, and
  `conjugateMatchesSup` checks BOTH directions at 1e-05 relative. `plqCheck`'s own header, and
  CCA2's `CLAUDE.md`, both say the sampled sup is a LOWER bound and that only `f* < sup` is a
  definite defect -- but the upper assertion at 1e-05 quietly demands the opposite.
- **Fixed at the oracle, not at the tolerance.** For a quadratic `f` the objective `<s,x> - f(x)`
  is quadratic, so its sup over a convex polygon is attained at a vertex, at a per-edge stationary
  point, or at the interior stationary point -- a finite list, computed exactly. `supOverDomain`
  now takes the max of the grid and that list. Both are lower bounds, so the result is still a
  lower bound and can only move toward the truth; for a quadratic on a bounded polygon it IS the
  truth. Measured after: 2.220e-16 on both points above, and the test passes in 28 s.
- **Expect the `f* < sup - tol` direction to fire more readily now.** That is the definite-defect
  direction, so anything it newly catches is a real minorant the coarse grid was hiding. The
  `--verylong` gate is the place that shows up.

## 2026-08-25 (G14) — every `verylong` red measured so far belongs to `plq_1piece`, and `conj` gets the same inputs right

Each red was taken to its own fixture and the SAME input put through `QuaPol.conj`:

    red                                class used     conj on the same input
    ---------------------------------  ------------   ---------------------------------------
    pce2* (was testPCE2)               plq_1piece     EXACT at all nine dual points
    testPSqroot                        plq_1piece     EXACT (1.75, where legacy returns -5)
    testOpenconvex                     plq_1piece     plq_1p REFUSES, and correctly
    rectConjugateMatchesTheSup         plq_1p         FALSE red -- the oracle; now green

- **`testPSqroot`:** legacy returns the objective's value at the VERTEX (-4,-3) where the sup is
  strictly inside the edge to (-1,1) at t = 0.75. G4's shape of defect, different implementation.
- **`testOpenconvex`:** legacy's envelope literally contains **2147483647** = `intmax('int32')`,
  the ray DIRECTION MARKER being read as an ordinary coordinate, and exceeds f by 2.147e+09 at
  (-1,0). `plq_1p` raises `convEnvUnbounded:minusInfinity` instead, which is the right answer --
  `x*y` on that half-strip runs to -inf along x = -1, so conv f = -inf. `quaPolToPlq`'s header
  records this defect as fixed *in plq_1p*; `plq_1piece` never got the fix.
- **`testPCE2`:** the one red everyone quotes. `QuaPol.conj` on {(0,0),(1,0),(2,1)} with `x*y` is
  exact at all nine dual points; the legacy route gives f*(0,0) = 0.0429 against a sup of 0.

**Consequence for the session's question.** "Fix any remaining bug for conj" is much closer to done
than the red count suggests: no measured red implicates `conj`. What the reds pin is `plq_1piece`,
the class T6 intends to retire, and `DECISIONS.md` 2026-08-19 already records that swapping the
fixtures to `plq_1p` is a MIGRATION with its own defects rather than a free change. Deciding
between migrating them and marking them as legacy pins is a call for Yves, not a bug hunt.

Still unmeasured: `testMax3` and `testMaxBiconjugate` (both `plq_1piece`, so expect the same), and
testcPLQ's `rectMaximumIsTheConjugateOfTheWholeDomain` and
`rectBiconjugateIsAConvexUnderestimator`, which DO use `plq_1p` and are the two worth measuring
next.

## 2026-08-25 (G10/G1) — the hole is `assemblePieces`, and G10 and G1 are the SAME defect

Two findings, both from probe tracking rather than argument.

**They are one defect.** G10's hole -- the point `(-0.866025, 0.5)` where `conj`'s numeric route
returns `Inf` on the G4 fixture -- was located in each operand's own faces, the way `TODO.md` G1
says to locate its own: the point lies in g1's face **21** and g2's face **6**, and the fold
produced **no piece with src [21 6]**. That is G1's signature verbatim, on a different fixture. One
fix closes both.

**But it is NOT `clipByFace`.** G1's entry names "the operand SWAP at the top of clipByFace,
`clipPolyByConic`, and the three reduction passes" as the candidates. Measured, with the probe
tracked through every one of them:

    enters clipByFace inside BOTH faces
    after the clip loop        -> covered by 1 piece
    dropDegeneratePieces       -> covered by 1 piece
    dedupPieces                -> covered by 1 piece
    dropSubsumedPieces         -> covered by 1 piece
    insertGlobalPassthrough    -> covered by 1 piece
    ...and the assembled QuaPar returns Inf there

So the cell is produced, survives every straight clip, the curved cut and all four reduction
passes, and is still covering the point when `assemblePieces` is called. **The loss is inside
`assemblePieces`.** That is exactly what `SUPPORT_MATRIX.md` section 4.6 already calls "the live
one": `matchHalfEdges` pairs curved with curved, and the subdivision is consistent per face PAIR
rather than globally.

**The instruments stay**, because three separate investigations this session began by adding them
by hand:
  * `MAXQP_PROBE` -- set a point and `clipByFace` names the exact straight clip, curved cut or
    reduction pass at which a pair that CONTAINS it stops containing it, with the half-plane and
    the probe's slack;
  * `dropDegeneratePieces` reports what it removes and why, at `MAXQP_ASSERT >= 0` -- deliberately
    not `>= 1`, because a violated invariant aborts before the drops and 0 is the setting you want
    when chasing a hole;
  * `splitDeclined` names why `splitUnboundedAtOneCrossing` gave up.

**Next step is `assemblePieces`, and it is the documented redesign** (globally-consistent
subdivision), not a small fix. Do not attempt it as a patch to `matchHalfEdges`.

## 2026-08-25 (overnight, G1/G4/G10) — the assembly hole is `matchHalfEdges` rejecting curved-vs-straight, and the SAGITTA is the discriminator

Four measurements, an attempt that got two stages further, and the blocker behind it. The attempt
is kept verbatim at `.claude/assembly_attempt_2026-08-25.diff` -- it is not committed because it
leaves the tree red, but it is 90% of the work and should be re-applied, not re-derived.

**1. The piece is produced; assembly loses it.** Tracked with `MAXQP_PROBE`, the cell covering the
uncovered point survives every straight clip, the curved cut and all four reduction passes, and is
still covering it when `assemblePieces` is called.

**2. WHICH piece, and why it is fragile.** Piece 28, src `[21 6]`, is an UNBOUNDED piece whose two
own vertices are

        [-0.20801507 -0.90722285]   and   [-0.20801318 -0.90721878]

**4.5e-06 apart**, against `matchHalfEdges`' `tolPos = 1e-3`. It is a cone whose apex arithmetic
split into two points. Its zero-length segment matched a neighbour it has no business matching,
`buildGlobalVertices` unified its vertices with distant ones, and face 28 came out spanning a
different region entirely.

**3. With those sub-tolerance edges collapsed FIRST, the real obstruction appears, and it is
`SUPPORT_MATRIX` 4.6's "live one" exactly.** `matchHalfEdges` then reports TWO orphans, each with
**zero candidates** -- rejected at generation, not lost to greedy consumption:

        piece 1 src [1 3]  straight  (-0.209275,-0.909936)->(-0.210638,-0.912444)
        piece 2 src [2 4]  CURVED    (-0.210638,-0.912444)->(-0.209301,-0.909885)   dist 5.76e-05

They are the same boundary, recorded as a chord by one piece and as an arc by the other, and
`matchHalfEdges` rejects the pair by rule: "a curved half-edge only ever matches another curved
one". That rule is CORRECT as far as its own reason goes -- an arc and its own chord share both
endpoints while bounding different sets -- so the fix is not to drop it.

**4. THE DISCRIMINATOR IS THE SAGITTA, and it separates the two cases by five orders of
magnitude.** A genuine chord departs from the conic in the MIDDLE; an edge that two pieces recorded
differently does not. Testing the straight edge's midpoint against the conic, measured on this
fixture:

        true pair   piece 1 src [1 3] vs piece 2 src [2 4]    residual 1.172e-06
        false pair  piece 2 src [2 4] vs piece 3 src [3 1]    residual 1.127e-01

Use it RELATIVE, not against a fixed tolerance: accept when the midpoint is no further off the
conic than the ENDPOINTS are (`rM <= max(1e-12, 10*max(rA,rB))`). A fixed 1e-06 threshold rejects
the true pair by a hair, because its endpoints already sit ~1e-06 off the conic having come from a
different arithmetic path.

**What the attempt achieves and where it stops.** With (a) sub-tolerance edges collapsed -- arcs
included, since an arc 5.76e-05 long has area ~3e-09 and is noise wearing a conic -- and (b) the
sagitta test, case 21 gets PAST assembly, and `case 8`'s A.5 triangle stays exact to 3.6e-15
throughout. It then dies one stage later in `clipArcByHalfPlane:internal`, "no crossing of the clip
line found within the arc's own u-span": exactly one endpoint outside, so a crossing must exist,
but on an arc that short the quadratic's discriminant comes out negative and no real root is found.

**Next step, and it is small and named:** make `clipArcByHalfPlane` handle the degenerate case --
one endpoint outside, no real root in span means the arc is tangent to the clip line within
tolerance, and the answer is to cut at the outside endpoint. Then re-apply the diff.

**Why it is not committed.** It turns `conjEdgeLowerBoundTest`'s two case-21 tests red, correctly:
they pin that the fixture is refused by a refusal we RECOGNISE, and `clipArcByHalfPlane:internal`
is a raw internal error. fast was 301/2 with it, 303/0 without.

## 2026-08-25 (overnight, continued) — `clipArcByHalfPlane` FIXED; and the assembly wall is the one its own history predicts

**FIXED, and it is independent of everything else.** `clipArcByHalfPlane` raised
`clipArcByHalfPlane:internal`, "no crossing of the clip line found within the arc's own u-span".
That error is not reachable geometrically: the branch runs only when exactly ONE endpoint is
outside, so `v(u) = nrm*X(u)' - c` has opposite signs at the two ends and is continuous along the
arc, and a root exists between them by the intermediate value theorem. What failed was the
closed-form SOLVE -- on a short arc (endpoints 5.76e-05 apart) the discriminant is a difference of
nearly equal quantities and comes out negative, so `solveQuadLocal` returns nothing.

Bisection cannot fail for the same reason the root must exist. It is now the closed form's
BACKSTOP: the quadratic still answers every case it can, and bisection runs only when it returns
nothing in span. `clipArcByHalfPlaneTest` 7/0, fast 303/0.

**The assembly wall, measured.** With the attempt re-applied on top of that fix, case 21 gets one
stage further again and then stops on a THIRD disagreement of the same kind:

    piece 11 src [8 3]  (-0.207912,-0.907429)->(-0.208593,-0.908683)   length 1.43e-03
    piece 18 src [10 4] (-0.209275,-0.909936)->(-0.208015,-0.907223)   length 2.99e-03
    closest candidate rejected on POSITION at dist 1.43e-03, against tolPos = 1e-3

Piece 11's edge is very nearly a SUB-SEGMENT of piece 18's -- a T-junction -- but its interior
endpoint lies **9.1e-05 off** piece 18's edge, so it is not a clean one, and three half-edges
orphan together (pieces 11, 12 and 18) because one long edge would have to pair with two short
ones.

**That is the wall `assemblePieces`' own HISTORY describes**, and the sequence of scales makes it
concrete: the same boundary is disagreed about at 5.76e-05, then 9.1e-05, then 1.43e-03, while
genuine features of this arrangement live at the same scales. No single tolerance separates them --
which is exactly what that header says, and why it proposes vertex PROVENANCE (tagging each edge
with which g1/g2 face-pair boundary produced it) as the resolution, noting it "turned out to be
unnecessary" for the cases known then. This fixture is the case that makes it necessary.

**So item 1's residue is now: implement edge provenance in `assemblePieces`.** The two fixes in
`.claude/assembly_attempt_2026-08-25.diff` are prerequisites and are correct as far as they go;
they are not committed only because without provenance the fixture still ends in a raw internal
error. Do not spend another night on tolerance tuning: three independent tolerances have now been
measured to be un-separable on this input.

## 2026-08-25 (overnight, item 2) — `rectMaximumIsTheConjugateOfTheWholeDomain` is a HOLE, and it does not reach `conj`

- **The red is a hole, not a wrong value:** `PRect f* uncovered at (-0.5,2)`. Same symptom class as
  item 1, in different code -- this is cPLQ's symbolic cross-piece `maximumConjugate`, not
  `maxQuaPar`.
- **`conj` gets the same domain RIGHT.** PRect is `x*y` over two adjacent polygons whose union is
  the hexagon `[-5,-4; 0,-4; 2,0; 2,1; 1,3; -5,5]`. Put that through `QuaPol.conj`: it returns a
  QuaParCPLQ in 151 s and `f*(-0.5,2) = 37.5`, which is the exact sup (attained at the vertex
  (-5,5): `-0.5*(-5) + 2*5 - (-5)(5) = 37.5`). Exact at every other probe point too.
- **Why the difference, and it is Step 0.** Both of PRect's pieces carry the SAME quadratic, so
  `mergeSameQuadFaces` deletes the shared edge and hands Case C ONE face. The cross-piece maximum
  that holes is then never performed. So this red cannot be reached through `conj` with this
  fixture -- it is reachable only by calling `plq.maximum` on a two-piece `plq` directly, which is
  what the test does.
- **The general question is still open and is NOT this fixture.** A QuaPol whose two faces carry
  DIFFERENT quadratics survives Step 0 and would exercise the same cross-piece max. Attempted with
  `x*y` and `x*y + x` over PRect's two polygons: `conj` did not finish in 40 minutes, so this is
  unmeasured rather than answered. A cheaper two-piece fixture is what that needs.
- **Split done.** The `max` stage was computed as `tri.maximum`, and `plq.maximum` re-runs
  convexEnvelope and conjugate on every piece before the cross-piece step -- work the `conj` stage
  has already cached and `rectConjugateMatchesTheSup` has already verified. It now builds from the
  cached `conj` stage via `crossPieceMax`, so attempts at this red stop paying for the per-piece
  conjugates. `rectBiconjugateIsAConvexUnderestimator` was building the same `plqStage` key by a
  DIFFERENT expression, so a cold `max` cost whichever path got there first; both now use the same
  one.

## 2026-08-25 (overnight, B3) — the POINT case is closed; the block was `QuaPar.eval`, not the mathematics

B3 was three refusals for a full-domain quadratic that is not strictly convex, all classified
earlier today with their closed forms: an EMPTY dom f* (a negative eigenvalue), a POINT (Q = 0,
f affine), and a LINE (PSD of rank 1). All three were recorded as blocked on the return type.

**One of the three was not.** dom f* being a single point is a NEEDLE -- `nv=1, ne=0, nf=0` -- and
`QuaPar`'s constructor has always anticipated exactly that shape:

    if ismember([obj.nv,obj.ne],[1,0;2,1],'rows'), obj.nf=0; end   %needle / segment / ray

What was missing was one branch of `QuaPar.eval`. Its first branch handles a full domain, its
second handles an EDGE CHAIN (`max(F,[],'all')==0`, which needs at least one edge), and a needle
matches neither: the face loop then runs over an empty `P` and every point comes back `+inf`,
including the needle's own vertex. Measured directly -- a hand-built needle at (3,-2) carrying -5
evaluated to `Inf` there.

So the domain existed in the representation and nowhere else. `eval` now has the branch, and
`conjCPLQ` Case A returns the answer instead of raising: for f = 3x - 2y + 5 it gives -5 at (3,-2)
and +inf elsewhere, which is the definition.

**The other two stay refused, and for different reasons.** A LINE is not a segment between two
vertices, which is all the `dim < 2` mesh can express; an EMPTY domain is not `nf=0` (that means
`dim < 2`, not empty) and has no encoding at all. Those two remain on the return-type work.

**Tests changed**, with the umbrella CLAUDE.md 8 sentence: `fullDomainNonStrictlyConvexIsClassified
NotLumped`'s table listed the affine row as raising `conjugateIsAPoint`, and
`theFullDomainRefusalsCarryTheClosedForm` used the affine case to check that a refusal carries its
closed form. Both pinned a gap that is now an answer, so the affine row leaves the table, the
message test is repointed at the rank-1 PSD case (still a refusal), and a new test asserts the
returned needle's value and its +inf elsewhere.

## 2026-08-25 (overnight, item 3) — the biconjugate timeout is DOWNSTREAM of item 2's hole

Where the >3600 s goes, from the two `--verylong -j 4` gates run earlier today plus tonight's split:

    tri + conj  (per-piece triangulate and conjugate)      ~28 s   -- rectConjugateMatchesTheSup
    cross-piece max                                      ~1534 s   -- the rest of the 1562 s 'max'
    biconjugateF                                         >2000 s   -- and does not finish

So the cost is not distributed: it is almost entirely the cross-piece maximum and then
`biconjugateF`, and only the last of those fails to terminate.

**And it is being fed a broken input.** Item 2 established that the `max` stage on this fixture
leaves a HOLE -- `PRect f* uncovered at (-0.5,2)`. `biconjugateF` computes f** = (f*)*, so it is
being asked to conjugate a mesh that does not cover the plane. That is the first thing to rule in
or out, and it is cheap to test now that the stages are split: repair or excise the uncovered
region in the cached `max` result and re-run `biconjugateF` alone. If it terminates, the timeout is
a symptom of item 2 and not an independent defect; if it does not, the two are separate.

A stage-by-stage timing run was started tonight to confirm the split above end to end. It had not
finished after several hours and was stopped when the run ended, which is itself consistent with
the numbers: the chain is dominated by two stages that together already exceed the budget.

## 2026-08-26 (G15) — REFUTED: the biconjugate timeout is NOT item 2's hole. It is the known scaling defect

Yesterday's G15 entry proposed that `rectBiconjugateIsAConvexUnderestimator`'s >3600 s is a symptom
of item 2's hole, on the grounds that `biconjugateF` computes f** = (f*)* and is therefore handed a
mesh that does not cover the plane. Measured today, and it is wrong.

    max stage (from cache)                       1016 s, 32 cells
    ...of which 1 of 10 probe points UNCOVERED   the hole is real and still there
    conjugateOfPiecePoly                          191 s, 32 cells -> 111 cells in 32 blocks

**The hole does not drive the cost.** The step that consumes the broken mesh -- and the only step
whose behaviour a hole could plausibly change -- takes 191 s of a budget of 3600, and it TERMINATES.
It does not choke, error, or spin on the uncovered region; it simply conjugates the cells that are
there.

**What is left is arithmetic, and it is the defect already on the list.** After that step there are
**111 cells in 32 blocks**, and `biconjugateF` must take the pointwise max across them. That is
`functionNDomain.maxOfList`, and `TODO.md`'s "STEP 3's CROSS-PIECE MAXIMUM DOES NOT SCALE" entry
measured exactly this shape on a different fixture: folds of 93, 294, 647, 1273 and 2087 s with the
cell count running 5, 14, 29, 45, 70, 86. At 111 cells, 3600 s is the expected cost rather than a
hang.

**So G15 folds into that entry rather than standing on its own**, and the fix that closes it is the
one already written there: merge same-function neighbours after each fold instead of only at the
end, starting with the POLYHEDRAL cells where `unionIsExact` already decides exactly. The prediction
this makes is testable and cheap -- if the surplus cells are the same "adjacent cells carrying the
same function, never merged" as on the quadrilateral, the 111 should collapse to of the order of a
dozen, and the timeout should disappear without anything in `biconjugateF` changing.

**Do not chase the hole for this.** It is a real defect (item 2 / G1 / G10) and worth fixing on its
own account, but it is not what makes this test time out.

## 2026-08-26 (B3) — the LINE case is closed, and it uncovered two real defects in `QuaPar`

B3's LINE case -- Q positive semidefinite of rank 1, so f is affine along null(Q) and f* is finite
only on a line -- was recorded as blocked on the return type: "QuaPar's dim<2 domain is a SEGMENT or
ray between two vertices, not a line". Measured, that was wrong twice over.

**A line IS representable**: two opposite RAYS from one point, `E(:,3) = 0`, which the constructor
has always accepted (it names `[nv,ne] = [2,1]` "segment / ray" beside the needle). What stopped it
were two defects, both in code that had nothing to do with B3:

**1. `QuaPar.eval`'s edge-chain branch treated every edge as a SEGMENT.** It called
`belongToEdge(V1,V2,x)` without the `isSegment` flag, so a RAY edge was read as the segment between
its apex and its DIRECTION marker, and the domain was silently truncated at a point with no
geometric meaning. Measured: the two-ray line through the origin gave `+inf` at (2,2).

**2. `belongToEdge`'s ray test only worked for rays pointing into the positive quadrant.** Its range
check was

        b(b) = all(V1 <= Xb + tol, 2);

a COORDINATE-WISE comparison. The contract is one dot product -- X is on `[V1,V2)` iff it is
collinear and `(X - V1)·(V2 - V1) >= 0` -- and that is direction-agnostic. With the flag fixed but
this not, the ray from (0,0) toward (-1,-1) still reported `+inf` at (-2,-2), a point on it. Pinned
by `QuaParTest/belongToEdgeHandlesARayInANYDirection`, which checks seven directions including the
axis-aligned and negative ones, plus the opposite ray in each case.

**Both fixed, so `conj` returns the answer.** Verified against the closed form
`f*(s) = 1/2 (s-L)' pinv(Q) (s-L) - kappa` on two fixtures, `x^2/2` and `(x+y)^2/2 + 2x`, at four
points along the line and two off it -- exact.

**Only the EMPTY case remains a refusal**, and it is the one with no encoding at all: dom f* empty
is not `nf = 0` (that means dim < 2). B3 is now one of three rather than three of three.

**Worth noting for the return-type work**: two of B3's three "blocked on representation" items
turned out to be blocked on `eval` instead. The representation was already there both times.

## 2026-08-26 (item 1 / scaling) — the proposed fix is ALREADY IN PLACE; the refusal is `noSharedFacet`

`TODO.md`'s "STEP 3's CROSS-PIECE MAXIMUM DOES NOT SCALE" entry (2026-08-16) ends with **"Where to
start: merge same-function neighbours after each fold, not only at the end."** Measured today, that
is already what happens, so the entry's advice is stale and following it would produce nothing.

`functionNDomain.maxOfList` calls `maximumP(true)` after EVERY fold, and `maximumP` calls `mergeL`
twice. The existing `CCA2_TRACE_MAXP` instrumentation shows it working, on PRect3 (`x*y` over a
quadrilateral):

    [maxP] in=11 afterSplit=11 merge1=8 merge2=7 (18 s)
    CROSS-PIECE max: 7 cells in 44 s

11 cells in, 7 out -- merging per fold, and succeeding.

**What the cost actually is.** `region.mergeTally`, which counts refusals by reason:

    cross-piece max -- 11 attempts:
       noSharedFacet    7
       okLinear         3
       emptyOperand     1

**`noSharedFacet` is 7 of 11**, and it is not a decision that the two cells cannot be merged -- it
is `sharesFacet` returning false because it could not PROBE the edge. Its own header says so: "this
routine is a LOCAL necessary condition ... a false here costs compactness, never correctness". So
the surplus cells are pairs that `unionIsExact` was never asked about, because a probe on
`probeOnConstraint` failed first.

**That is where to start instead**, and it is a much smaller target than the entry suggests: make
`sharesFacet` probeable on the edges it currently gives up on (the `probeOnConstraint` fallback when
`vx+0.1`/`vx-0.1` are both infeasible), or let a failed probe fall through to `unionIsExact` rather
than short-circuiting to "no". The second is strictly safer -- it can only ask a question that is
currently skipped.

**Scale note.** PRect3 is small enough to iterate on: 44 s for the cross-piece max and 68 s for
`biconjugateF`, against 1016 s and >3600 s for PRect. Use it, not PRect, when working this.

## 2026-08-26 (item 2, assemblePieces) — the wall is a T-JUNCTION whose vertex is 9.1e-05 off the edge

Not attempted tonight beyond locating it precisely; recorded so the next attempt starts here.

`insertGlobalPassthrough` already exists to eliminate exactly the failure the assembly attempt hits.
Its header describes the mechanism: two adjacent pieces share a boundary; one side gets split at a
crossing point P and the other does not, so a long ray cannot pair with a segment, and re-inserting
P as a vertex of the unsplit piece restores the pairing.

**And it states the tolerance regime it was built for: "P ... lies EXACTLY on the decided ray
(measured perp ~2e-15 on the arc-vs-arc fixtures)".**

On `TODO.md` G4's fixture the T-junction vertex is **9.1e-05** off the edge it should split -- ten
orders of magnitude worse than the case that routine was written against, and the same scale as this
arrangement's genuine features (5.76e-05, 1.43e-03). So the insertion is not failing because a
tolerance is a little too tight. **The two pieces genuinely disagree about where the shared boundary
is**, having computed it by different arithmetic paths, and no threshold can separate that
disagreement from a real feature -- which is the conclusion `assemblePieces`' own HISTORY already
reached and why it proposes vertex PROVENANCE.

**So the next attempt should tag half-edges with the (g1 face, g2 face, constraint) that produced
them and pair on that**, using geometry only to break ties among candidates that already agree on
provenance. Do not start by widening `insertGlobalPassthrough`'s tolerance: that is the third
tolerance on this fixture to be measured un-separable, after `matchHalfEdges`' tolPos and
`pieceRecessionRays`' sign tests.

## 2026-08-26 (item 1, second pass) — REFUTED: the surplus cells are NOT unprobed facets

The entry above proposed the fix: `noSharedFacet` is 7 of 11 merge refusals, `sharesFacet`'s
`isconvex` is only a LOCAL probe that returns false when it cannot probe the edge, so let a failed
probe fall through to `unionIsExact` -- "strictly safe, since it only asks a question currently
skipped".

**Built it and measured it. The reasoning is right and the conclusion is wrong.**

The fallback keeps a pair whose facet has ALREADY matched (opposite constraints, matching edge
vertices) but whose local probe failed, and hands it to `unionIsExact`; conservatively, only when no
probe-passing candidate exists and exactly one unprobed candidate does. On PRect3:

    before   11 attempts:  noSharedFacet 7, okLinear 3, emptyOperand 1
             7 cells, 44 s
    after    13 attempts:  noSharedFacet 5, okLinear 3, unprobedFacetTried 2,
                           lin_exactAnotInB 2, emptyOperand 1
             7 cells, 44 s

**Two pairs were newly asked, and `unionIsExact` REFUSED BOTH** (`exactAnotInB` -- A is not
contained in B'). The cell count and the runtime are unchanged.

**So the 7 refusals decompose as: 2 unprobeable facets that are genuinely unmergeable, and 5 pairs
with no shared facet at all** -- cells that are simply not adjacent, correctly not merged. There is
no surplus hiding behind the probe. `merge` is already achieving what it can on this fixture, and
7 cells is what it can.

**Reverted.** The change is sound and costs two extra LPs to learn nothing; keeping it would add a
code path with no measured benefit.

**What this means for the scaling defect.** Two of its three plausible causes are now eliminated:
the merge is not mis-SCHEDULED (it already runs per fold) and it is not being BLOCKED by the local
probe. So on a fixture where the cell count really does blow up -- the quadrilateral at 5, 14, 29,
45, 70, 86 -- the growth must be in the SPLIT, i.e. `maximumP` creating cells, not `mergeL` failing
to remove them. `[maxP] in=N afterSplit=M` is the number to read next, and on PRect3 it is
`in=11 afterSplit=11`, which creates nothing. **Measure afterSplit on the quadrilateral before
touching merge again.**

## 2026-08-26 (item 1, third pass) — the blow-up MEASURED: it is `mtimes` pairing and CURVED-facet refusals

The refutation above says: if the merge is neither mis-scheduled nor blocked by the local probe,
read `afterSplit` on the fixture where the count really does blow up. Done -- `x*y` over
conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}, 6 pieces after triangulation:

    [maxP] in=14 afterSplit=15 merge1=14 merge2=12 ( 45 s)
    [maxP] in=29 afterSplit=32 merge1=28 merge2=28 (126 s)
    [maxP] in=42 afterSplit=45 merge1=44 merge2=43 (294 s)
    [maxP] in=60 afterSplit=68 merge1=48 merge2=44 (551 s)
    [maxP] in=63 afterSplit=63 merge1=58 merge2=56 (571 s)
    CROSS-PIECE max: 56 cells in 2944 s

**First, the entry's own numbers are STALE and should be replaced.** It records 5, 14, 29, 45, 70,
86 cells and 73 minutes, measured 2026-08-16. Today the same fixture gives 14, 29, 42, 60, 63 -> 56
cells in 49 minutes. Intervening work has already taken a third off both. Still far too slow, but
re-measure before quoting it.

**The split is NOT the growth.** `afterSplit` exceeds `in` by 1, 3, 3, 8, 0 -- fifteen cells created
across the whole run. And the merge is working: fold 4 takes 68 down to 44.

**The growth is in the PAIRING.** `in` for each fold is what `objR * g` produced -- `mtimes`
intersects every cell of the accumulated result with every cell of the next operand -- and it runs
12 -> 29, 28 -> 42, 43 -> 60. That is the quadratic blow-up of an arrangement overlay, and it is
inherent to folding pairwise; what is supposed to contain it is the merge that follows.

**Where the merge loses, by count** (3236 attempts):

    noSharedFacet                   2743      pairs that are simply NOT ADJACENT -- with 63 cells
                                              there are ~2000 pairs, so this is expected, not a defect
    quadFacet_exactAnotInB           374      <-- THE REAL ONE
    okQuadFacet                       11
    okLinear                          26
    everything else                  <30

**`quadFacet_exactAnotInB` is the target: 374 refusals against 37 successes.** These are ADJACENT
pairs sharing a CURVED facet, put to `unionIsExact`, which answers "A is not contained in B'". That
is the conic case, and it is where the surplus cells live -- not in the probe (refuted above), not
in the schedule (refuted above), not in the split (refuted here).

**Next step, and now it is specific:** take one `quadFacet_exactAnotInB` pair from this fixture and
determine whether the refusal is CORRECT (the union really is not convex, so the cells must stay
separate and 56 is near-optimal) or CONSERVATIVE (the LP relaxation over a curved facet cannot
certify a union that is in fact exact). `region.m`'s own header says the curved certificate is
deliberately relaxed -- "makes the certificate harder to obtain, never wrong" -- so conservative is
the likelier of the two, and that relaxation is the thing to sharpen.

## 2026-08-27 (B3) — the EMPTY case closes too, and B3 is DONE. The pattern is worth naming.

dom f* being EMPTY was the last of B3's three, and the one that looked genuinely blocked: "nf = 0
means dim < 2, not empty, so there is no encoding at all".

**That is true of the DOMAIN and beside the point. The FUNCTION is representable.** f* = +inf
everywhere, and a full-domain mesh (`nv=0, nf=1`) whose constant is `+inf` evaluates to exactly
that. Checked before adopting it: `eval` returns Inf at every point including far out, `isMeshed`
and `kind` behave, `isConvex` says true (the indicator of the empty set IS convex), and addition
propagates +inf rather than producing NaN.

**ONE CONSEQUENCE, stated rather than left to be discovered.** Conjugating the result gives -inf
everywhere -- mathematically right, since (+inf)* = -inf, and NOT a PLQ function. So a caller
computing f** on a quadratic with a negative curvature direction now gets a -inf mesh instead of an
error. `biconjCPLQ` short-circuits convex inputs and never reaches it, and nothing else in the
repository conjugates twice without checking, so nothing is broken today. Pinned by
`conjCPLQTest/theEMPTYDomainConjugateRoundTripsToMinusInfinity` so it is deliberate, and that is
the site to revisit if a return type ever gains a way to say "not a PLQ function".

**THE PATTERN, which is the transferable part.** All three of B3 were recorded as blocked on the
RETURN TYPE. None of them was:

    POINT  (2026-08-25)  the needle existed; `QuaPar.eval` had no branch for it
    LINE   (2026-08-26)  two rays existed; `eval` ignored the ray marker, and `belongToEdge`'s
                         ray test was coordinate-wise so it only worked in the positive quadrant
    EMPTY  (2026-08-27)  a +inf full-domain mesh existed; nobody had tried it

Three for three, the representation was already there and the block was one routine that did not
know about it. **Before recording anything else as "blocked on the return type", build the object
by hand and evaluate it.** That is a five-minute check, and here it would have saved the entry
being written three times.

## 2026-08-27 (item 1) — the curved tightening is skipped on UNBOUNDED cells, which is most of them

The open question was whether `quadFacet_exactAnotInB` (374 refusals against 37 successes on the
quadrilateral) is a CORRECT refusal or a CONSERVATIVE one. Measured, and it is neither exactly:
**the tightening that would decide it does not run at all on most of these cells.**

`maxAffineOverRegion` has an exact closed-form bound for a region with conic facets -- vertices plus
the points where the arc's tangent is perpendicular to the objective -- added 2026-08-17 for exactly
this conservatism. It bails to the loose LP through several early returns. Instrumenting each,
on PRect3:

    tighten_polyhedral          6      the region has no conic facet: LP already exact, fine
    tighten_unboundedVertex     2      <-- the tightening GIVES UP

`tighten_unboundedVertex` is the guard "every vertex must be a finite point for the compactness
argument to apply" -- it returns as soon as any vertex is at `intmax`, the ray marker. **A conjugate
is full of cones**, so this fires on a large fraction of the cells that have a conic facet at all,
and those pairs fall back to the LP over the linear relaxation, which is precisely the bound whose
conservatism produced `exactAnotInB`.

**The guard is sound and the argument behind it is too narrow.** A linear form on an UNBOUNDED
region need not attain its max, so the vertex-plus-tangency candidate set is not complete -- correct.
But the max is still bounded, and decidable, whenever the form is bounded on the region's RECESSION
CONE, which is exactly the situation where `unionIsExact`'s question has a finite answer. The
candidate set then needs the recession directions added to it, not the whole region rejected.

**So the answer to the question is: CONSERVATIVE, and located.** Extend the tightening to unbounded
regions by testing the objective on the recession cone first -- unbounded above there means "no
certificate" as now, bounded means the finite candidates plus the cone's extreme rays decide it.
`pieceRecessionRays` already computes exactly those directions for `maxQuaPar`'s pieces.

**Not attempted here**, because the last two attempts on this defect were both refuted by
measurement after the code was written, and this one deserves a bounded experiment first: take ONE
`exactAnotInB` pair from PRect3, work its union out by hand, and confirm the union really is convex
before building anything.

## 2026-08-27 (item 2) — provenance is NOT the fix: the arrangement is finer than the tolerance

`assemblePieces`' own HISTORY proposes vertex PROVENANCE as the resolution for its residual
matching failures, and yesterday's entry accepted that and said the next attempt should tag
half-edges with the face pair and constraint that produced them. **Measure first**, because
threading a per-edge tag through thirteen producers (`clipPolyHalfPlane`, `clipPolyByConic`,
`splitCell`, `boundedPiece`, `splitTwoArcPiece`, `insertPassthroughVertices`, ...) is a large,
error-prone change. Measured on `TODO.md` G4's fixture, 30 pieces:

    edge lengths      min 4.484e-06   median 2.987e-03   max 9.982e-01
    edges shorter than matchHalfEdges' tolPos = 1e-3      9 of 42
    pieces whose whole DIAMETER is under tolPos           4 of 30
    pieces under 10x tolPos                               9 of 30

**The median edge is three times the matching tolerance, and a fifth of the edges are shorter than
it.** The orphaned piece 1 (src [1 3]) has its two vertices 3.0e-03 apart -- the entire piece is
barely above the resolution at which the assembly claims to identify points.

**That is not a pairing-strategy problem and provenance will not fix it.** Provenance answers "which
of several candidate neighbours is the right one" -- a disambiguation. Here the difficulty is prior:
at this scale the pieces' own COORDINATES, computed along different arithmetic paths, disagree by
more than the features they are meant to distinguish. Tagging an edge with the constraint that made
it does not make two disagreeing copies of that edge agree, and the fallbacks would still have to
compare positions.

**Two directions that do address it, neither cheap:**

  1. **Do not produce the sub-tolerance pieces.** They are slivers from a fan of cells meeting near
     one point; a snap-round of the arrangement, or collapsing sub-tolerance cells before assembly
     (the parked diff's first half), removes the input rather than teaching the matcher to cope.
     The parked attempt shows this gets two stages further before hitting the next disagreement.
  2. **Make the coordinates agree.** The scale is set by float noise across arithmetic paths, so
     exact rational vertices for the polyhedral part would remove the disagreement rather than
     bound it. That is `TODO.md`'s `ratQ`/`conicMeet` direction, and it is a large piece of work.

**What is now ruled out**, and recorded so it is not proposed a fourth time: tolerance tuning
(three separate tolerances measured un-separable on this fixture), and provenance tagging (this
entry). The measurement to start any future attempt from is the edge-length distribution above --
if a change does not move it, it will not fix the assembly.

## 2026-08-27 (item 4) — the cPLQ hole is the KNOWN `mergeL`/`removeTangent` gap, and it is NOT the assembly's problem

G11's remaining `plq_1p` red -- `rectMaximumIsTheConjugateOfTheWholeDomain`, the hole at
(-0.5, 2) -- was worth checking against item 2's finding: both are "the assembled result does not
cover a point", so the natural guess is one cause. Measured, and it is not.

    cells: 32;  eval at the hole = NaN
    nearest cell misses by 1.668523e-02

**1.67e-02 is a genuine GAP, not a near miss.** Item 2's assembly failures live at 1e-6 to 1e-3, at
or below the tolerances meant to resolve them. Here the nearest cell is four orders of magnitude
further away than that: no tolerance, no snapping and no provenance scheme would close it, because
there is nothing there to close. A region is simply absent.

**And it is already a named, documented gap.** The stack through the failing run is

    region.removeTangent <- functionNDomain.mergeL <- maximumP <- maxOfList <- plq.maximumConjugate

which is `DESIGN.md`'s "`mergeL`/`removeTangent` exact-tie-point gap" -- recorded there five times,
including "the remaining known bug ... which `QuaParCPLQ.conj` now inherits when composing".
`removeTangent` deletes a constraint it believes is tangent-redundant; at an exact tie point it
deletes one the region needs, and the region then fails to cover territory that is its own.

**So item 4 is not a new defect and not the same defect as item 2.** Two different "hole" symptoms
with two different causes, in two different implementations:

    maxQuaPar assembly   1e-6..1e-3   coordinates disagree below the matching tolerance   (item 2)
    cPLQ maximumConjugate   1.7e-2    a constraint is wrongly deleted as tangent-redundant (here)

The second has a name, a location, and a documented history. It should be worked as itself --
`removeTangent`'s tie-point handling -- rather than as part of the assembly work.

## 2026-08-27 (G17, corrected) — `removeTangent` DISCARDS whole regions, and that is the hole's mechanism

**First, a correction to the entry above.** It identified the hole as "the known
`mergeL`/`removeTangent` exact-tie-point gap" on the strength of a stack trace, and called the
mechanism "removeTangent deletes a constraint the region needs". **The deletion half of that is
impossible**, and the reasoning is one line: deleting a constraint makes a region BIGGER. A hole is
a point covered by NO region, and enlarging regions cannot create one. The identification was right
about the routine and wrong about how it does it.

**Measured, with every branch of `removeTangent` instrumented.** On the PRect max stage, 5 events:

    deletedLinear     at ( 3.1591, -1.0227)   s_1 - (9*s_2)/5 - 5
    DISCARDED_REGION  at ( 3.1591, -1.0227)   148*s_1 - 196*s_2 + 14*s_1*s_2 + s_1^2 + 49*s_2^2 - 684
    deletedLinear     at ( 0.7368, -2.3684)   s_1 - (9*s_2)/5 - 5
    deletedLinear     at ( 3.5447, -0.0307)   s_1 - (104*s_2)/43 + ...
    DISCARDED_REGION  at ( 0.8766,  1.3033)   4*s_1*s_2 - 40*s_2 + s_1^2 + 4*s_2^2 + 40

**The three deletions are all far from the hole at (-0.5, 2) and, per the argument above, could not
have caused it either way.** What can, and is not in that entry's account at all, is
`removeTangent`'s `lin & ~tin` branch:

        if lin & ~tin
            obj = region.empty;     % the WHOLE region, discarded
            return
        end

**Two regions were discarded, and the second one's quadratic is satisfied at the hole**:
at `(-0.5, 2)`, `4(-0.5)(2) - 40(2) + 0.25 + 16 + 40 = -27.75 <= 0`. So the uncovered point lies
inside the constraint of a region that `removeTangent` threw away. That is the mechanism, and it is
a much sharper target than "the tie-point gap".

**Stated as a candidate, not a proof**: the discarded region has other constraints that were not
captured, so containment is established for the one named constraint rather than for the region.
Confirming it is the next step and is cheap -- capture the full region at the discard and test the
point against all of its constraints.

**Why the branch exists at all**: `lin & ~tin` means the midpoint probe says the quadratic holds
there but the tangent does not. For a genuinely tangent pair that is contradictory, so the code
reads it as "this region is empty". At an exact tie point the probe is unreliable, and a region
that is merely thin gets destroyed rather than kept.

## 2026-08-27 (G17, third pass) — REFUTED: `removeTangent` is not the cause at all

The entry above found `removeTangent`'s `lin & ~tin` branch discarding whole regions, noted that one
discarded region's single named constraint was satisfied at the hole, and called it "the mechanism"
while flagging that containment was established for one constraint rather than the region.
**Checked the whole region. It does not contain the hole. Neither does the other.**

    discard 1 at (3.1591,-1.0227):  4 constraints, 2 violated, worst +9.100000e+00
    discard 2 at (0.8766, 1.3033):  5 constraints, 1 violated, worst +3.983315e+00

Both misses are of order 1 -- nowhere near the point. The single-constraint check that looked
promising was simply the one constraint of five that happens to hold there; the region as a whole is
far away.

**So `removeTangent` is exonerated on this fixture**, and the whole G17 identification collapses:

    claim 1  "removeTangent deletes a constraint the region needs"   -- impossible (deletion enlarges)
    claim 2  "removeTangent discards a region containing the hole"   -- measured false, both discards

The only evidence ever behind G17 was a stack trace, and a stack trace of a symbolic WARNING at
that -- `removeTangent` is simply where cPLQ's `isAlways` calls are noisiest, so it appears in
warning traces from any run. That is not evidence of causation, and I treated it as such twice.

**What is actually known about this hole, and it is little:** 32 cells, the point (-0.5,2) uncovered,
the nearest cell missing by 1.67e-02. That gap is four orders above the `maxQuaPar` assembly's
scale, so it is still a genuinely absent region rather than a tolerance artefact -- that part stands.

**The next step is a bisection, not another hypothesis.** `maxOfList` folds cell groups one at a
time and the fold count here is small. Evaluate the accumulated result at (-0.5,2) after EVERY fold:
the first fold at which it becomes NaN identifies the operand pair and the step, the way
`MAXQP_PROBE` did for the assembly. Do that before naming a routine again.

## 2026-08-27 (G17, LOCATED) — the hole appears in fold 4's `maximumP`, after the pairing

The bisection prescribed by the refutation above, run rather than reasoned about. Fold the
per-piece conjugates one at a time and count how many regions contain (-0.5,2) at each step --
before the pairing, after `mtimes`, and after `maximumP`:

    start:  3 cells, covered by 1
    fold 2: covered before=1   after mtimes=1   after maximumP=1   ( 8 cells)
    fold 3: covered before=1   after mtimes=1   after maximumP=1   (19 cells)
    fold 4: covered before=1   after mtimes=1   after maximumP=0   (25 cells)   <-- LOST
    fold 5: covered before=0   after mtimes=0   after maximumP=0   (32 cells)

**The point is covered entering fold 4 and by the paired object, and uncovered after
`maximumP`.** So the loss is in `maximumP` -- the split-and-merge step -- at one identified fold,
not in the pairing and not in any earlier fold.

**A note on the probe, since it cost two attempts.** After `mtimes` each region carries TWO
functions, so `evalFunctionNDomain` cannot read a paired object -- it errors in `subs`. Coverage has
to be counted from the CONSTRAINTS at that stage, which is what this version does. Anyone
instrumenting this pipeline will hit the same thing.

**Every operand covers the point** (measured: the five groups give 37.5, 37.5, 2.5, 2.5, 2.5 at it,
and 37.5 is the truth), so this is purely a fold defect -- no per-piece conjugate is at fault.

**Where that leaves the accusation of `removeTangent`.** It is called from `mergeL`, which
`maximumP` calls twice, so it remains *reachable* at the right step -- but the previous entry
measured both of its region discards as far from the point, so it is not the mechanism by which
fold 4 loses it. `maximumP` also SPLITS cells before merging (`in=N afterSplit=M`), and a split that
drops territory would produce exactly this. **Next: instrument `maximumP`'s own split/merge stages
for coverage at fold 4** -- the same before/after count, one level deeper. That isolates it to the
split loop or to `mergeL`, and only then is naming a routine justified.

## 2026-08-27 (G17, one level deeper) — the loss is in `mergeL`'s FIRST pass, not the split and not the second pass

The next step prescribed above, run: instrument `maximumP`'s own split loop and its two `mergeL`
calls separately, counting coverage of (-0.5,2) at each. `g17e.m` (scratchpad) replicates
`maximumP`'s body by hand for fold 4 only, reusing `g17d`'s fold-2/3 setup and its
constraint-based `coveredBy` (a paired object cannot be read by `evalFunctionNDomain`, per the
prior entry's note):

    entering fold 4:                19 cells, covered=1
    fold 4 after mtimes:            35 cells, covered=1
    fold 4 after SPLIT loop:        39 cells, covered=1
    fold 4 after mergeL #1:         27 cells, covered=0   <-- LOST
    fold 4 after mergeL #2:         25 cells, covered=0

**The split creates 4 cells and loses nothing** -- covered stays 1 through mtimes and the split
loop, ruling out `maximumP`'s own split-and-simplify path (the `d1 = d1.simplifyUnboundedRegion`
calls at both branches) as the mechanism. **The first `mergeL` call is where the point's covering
cell is merged away without carrying its coverage forward**, dropping 39 cells to 27. The second
`mergeL` call (39->27->25) does not touch the hole again.

So the accusation narrows from "somewhere in `maximumP`" to "somewhere inside `mergeL`'s first
pass" -- which still calls `removeTangent`, but the previous refutation (2026-08-27, third pass)
already ruled that routine out on this fixture by containment check. **Next: instrument `mergeL`
itself the same way** -- coverage before/after each pairwise merge attempt within the single
`mergeL(objSplit)` call, to name which of its ~39-choose-2 attempts is the one that drops the
covering cell (or merges two cells into one that no longer contains the point, which would be a
different and more specific defect than a discard).

## 2026-08-27 (item 3) — hand-checked with a minimal witness: CONFIRMED conservative, not correct

Item 1 (above) argued analytically that `quadFacet_exactAnotInB` is conservative because
`maxAffineOverRegion` gives up on an unbounded region, and prescribed a bounded experiment before
touching the code: take one such pair and confirm by hand that the union really is convex.

**The real fixture is too expensive to use for that check.** `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}`
needs the symbolic conjugate route (irrational sub-triangle vertices, sqrt(15)/sqrt(30)) before a
single fold runs, and the fold's own `isAlways` calls on those nested radicals are slow enough that
a cold run of `.maximum` did not reach a single `quadFacet_exactAnotInB` case in over two hours.
Killed; not worth waiting on for a question that does not need this fixture at all.

**Built the smallest possible witness instead** (`region`'s constructor takes raw inequalities
directly, no plq pipeline needed):

    A = {s1 >= s2^2, s2 >= 0}            (unbounded, one conic facet, one linear)
    B = {s1 <= s2^2, s1 >= 0, s2 >= 0}   (bounded above by the same arc)

shared on `s2^2 - s1 = 0`. `A.unionIsExact(B, 1, 1)` returns exactly `(false, 'exactAnotInB')`.

**By hand, and confirmed by a 200,000-point sample with zero mismatches: `A union B` is the
quarter-plane `{s1>=0, s2>=0}`, which is convex.** So the refusal is wrong to refuse -- CONFIRMED
conservative, not a correct rejection of a genuinely non-convex union.

**The mechanism is sharper than "gives up".** `A.maxAffineOverRegion([-1 0])` (i.e. the max of
`-s1`, testing whether `s1 >= 0` holds on A) returns `st = 1`, UNBOUNDED ABOVE -- not `st = 2`,
undecided. That is because "unbounded regions keep the LP answer" (`maxAffineOverRegion`'s own
header) and the LP relaxation DROPS THE CONIC FACET ENTIRELY: with only `s2 >= 0` left, nothing
bounds `s1` at all, so the relaxed LP reports unbounded even though the true region, respecting its
one curved facet, has `s1 >= s2^2 >= 0` everywhere and the true max of `-s1` is exactly 0 at the
origin. The recession cone of the TRUE region (only the trivial direction, since any direction with
`s1` decreasing violates the parabola for large t) is what the relaxation cannot see.

**This is exactly the fix item 1 named**, now with a concrete, cheap, reproducible counterexample
rather than an argument: extend `maxAffineOverRegion`'s unbounded-region path to test the objective
on the region's actual recession cone (`pieceRecessionRays` already computes it for `maxQuaPar`
pieces) before falling back to the linear relaxation's LP. `region.m`'s temporary dump
instrumentation used to reach this point through the real fixture was reverted -- it is not needed
and should not be committed. The witness above (`region([s2^2-s1;-s2],...)` etc.) is cheap enough to
turn directly into a unit test for `maxAffineOverRegion` and `unionIsExact` once the fix lands.

## 2026-08-27 (item 1, REFUTED) — the recession-cone fix is UNSOUND for a curved region; caught before landing

Built the fix item 1 and item 3 both pointed at: extend `maxAffineOverRegion`'s unbounded-region
path by testing cRow against the region's recession cone (linear facets' extreme rays, plus the
conic's null direction, mirroring `pieceRecessionRays.m`'s method) before falling back to the
loose LP. It made the item-3 witness exact (`unionIsExact` on `A={s1>=s2^2,s2>=0}`,
`B={s1<=s2^2,s1>=0,s2>=0}` went from `(false,'exactAnotInB')` to `(true,'ok')`, and
`maxAffineOverRegion(A,[-1 0])` from `Inf/unbounded` to the correct `0/decided`).

**Then a second, harder probe on the SAME region A refuted it before it was ever run against the
test suite.** `maxAffineOverRegion(A, [0 1])` -- max of `s2` over A -- came back `val=0, st=0`
(decided, finite). **The true answer is `+Inf`**: `(100,10)` satisfies `s1=100 >= s2^2=100`, and
`s2` can be made arbitrarily large the same way for any `t`. A confidently WRONG answer, not a
conservative one -- exactly the failure mode this file's every LP certificate exists to prevent.

**Why: the theorem the fix relied on is for POLYHEDRA, and A is not one.** "cRow is bounded above
on a closed convex set C iff cRow.d <= 0 for every d in C's recession cone" holds for a
POLYHEDRON (`region.maxLinear`'s own domain), where the recession cone's straight-line extreme
rays really do capture every direction the set extends in. It does NOT hold for a set with a
CURVED boundary: A's recession cone (straight-line asymptotic directions only) is exactly the
single ray `(1,0)` -- moving along any FIXED direction with a nonzero s2-component eventually
violates `s1>=s2^2` for large t, confirmed by hand for both `x0=(0,0)` and `x0=(M,0)` for large M.
`cRow=(0,1)` satisfies `cRow.(1,0) = 0 <= 0`, passing the recession-cone test, and yet `s2` is
still unbounded on A -- along the CURVE `(t^2,t)`, whose direction `(t,1)/norm` does not converge
to any FIXED ray as `t->infinity` in the way a recession-cone argument needs (it does converge in
ANGLE to `(1,0)`, which is exactly why the straight-line test is fooled: the set's extent in the
`s2` coordinate grows without a matching straight-line witness).

**REVERTED. Do not retry the recession-cone-only construction; it needs a strictly stronger
ingredient than `pieceRecessionRays`'s straight-line rays provide.** The `region.m` edit
(`recessionRaysGeneral`, `recedesAllGeneral`, `quadNullDirsNumeric`) was fully reverted before
being run against any test (`git checkout -- region.m`); nothing shipped.

**What a correct fix needs instead.** The straight-line recession-cone check is still valid and
necessary for the LINEAR facets (unchanged from the existing LP). What it cannot see is growth
ALONG the conic arc itself: a parabola parametrizes as one coordinate quadratic in the other (here
`s1 = s2^2`), and `cRow` evaluated along that parametrization is a QUADRATIC in the free parameter,
whose behaviour as the parameter -> +-infinity (bounded above iff the parameter's own coefficient
in cRow's composition is <= 0, exactly the classic "does this parabola open the right way for this
linear functional" check) is what actually decides boundedness on the arc's unbounded branch --
not the arc's asymptotic straight-line direction. This is a DIFFERENT and more specific
computation than a recession cone: it needs the conic's own explicit parametrization along each
unbounded branch, checked against cRow directly, in addition to (not instead of) the existing
straight-line recession test for the linear facets.

**Before attempting that:** hand-verify it on this SAME witness first (A's arc parametrizes as
`s2=t, s1=t^2`; `cRow=(0,1)` composed along it is `t`, unbounded as `t->infinity` -- correctly
predicts the counterexample above; `cRow=(-1,0)` composed is `-t^2`, bounded above by 0 at `t=0` --
correctly predicts the finite answer that WAS right) before writing it into `region.m`, exactly as
item 1 already prescribed for the first attempt and this entry now re-prescribes for the reason it
failed.

## 2026-08-27 (G17, ROOT CAUSE FOUND) — `certifiesNonPositive` accepts a WRONG `quadprog` exit code as proof the region is empty

Traced one level deeper than "mergeL's first pass" (previous entry): replicated `mergeL`'s two
accumulation loops by hand on the cached 39-cell fold-4 state (`fold4_objSplit.mat`, scratchpad).
The covering cell (12) merges cleanly with cell 17 (still covers `s=(-0.5,2)`), then merges with
cell 21 -- `region.merge` reports SUCCESS and the result no longer covers `s`. **This is not a
discard: it is a wrong "the union is exact" certificate.**

`r1217` (cells 12+17 merged) and `r21` share a LINEAR facet (`s1-sqrt(14)s2+4<=0` against its
negative). Deleting it and intersecting what's left needs `unionIsExact`'s two subset checks; the
one that fails silently is "does `r21`'s one CURVED constraint,
`h = 40 s2 - 4 s1 s2 - s1^2 - 4 s2^2 - 40 <= 0`, hold everywhere on `r1217`?", answered by
`certifiesNonPositive`. **It returned `tf=true, why='ok'`. It is wrong**: at `s`,
`h(s) = 27.75 > 0`, and `s` is confirmed feasible for `r1217`'s own linear relaxation directly
(`all(A*s<=b)` is true) -- `r1217` is not empty and `h` is not `<=0` on all of it.

**The mechanism, reproduced directly.** `h`'s Hessian is `Q=[[-2,-4],[-4,-8]]`, eigenvalues
`{-10, 0}` -- CONCAVE, rank 1 (not the convex branch). `certifiesNonPositive`'s concave branch
maximises `h` over `r1217`'s linear relaxation via `quadprog(-Q,-L,A,b,...)` and treats exit code
`ef==-2` as "the region is empty, so `h<=0` holds vacuously" (`tf=true`). Run by hand with
`Display,'iter'`: **`ef` comes back `-2`, and the iterates do not converge** -- `Fval` diverges
past `-6.4e7`, the returned point is `(-3.2e6, 1.6e6)`, and the solver's own message says
"unable to find a point that satisfies the constraints within the value of the constraint
tolerance." **`s` is feasible; the region is not empty.** `quadprog`'s `-2` here is a
NUMERICAL-FAILURE code on a QP that is unbounded above with a rank-DEFICIENT (semidefinite, not
definite) Hessian -- not the "primal infeasible" code the surrounding logic assumes it always
means. The two situations share an exit code and `certifiesNonPositive` cannot tell them apart.

**Why this shape is exactly where it bites.** `r1217` is an UNBOUNDED cell (this whole
investigation has been inside `maxQuaPar`'s unbounded-cell folds), and `h` is concave with a NULL
DIRECTION (rank 1) -- precisely the combination (unbounded feasible region, semidefinite not
definite Hessian) that interior-point QP solvers are known to struggle to certify correctly:
`quadprog` cannot always distinguish "truly infeasible" from "feasible but the objective runs away
along a flat direction, and the algorithm never converges" when the Hessian has a null space
aligned with the region's own recession cone.

**This is G17's actual root cause**, three levels deeper than where it was first bisected
(`maximumP` -> `mergeL`'s first pass -> the merge of cells 12+17+21 -> `certifiesNonPositive`'s
`ef==-2` branch). It is NOT `removeTangent` (exonerated twice already) and NOT a discard anywhere
in `maximumP`'s split loop (also exonerated, previous entry) -- it is a false "yes" from the ONE
place in this file that calls `quadprog` rather than `linprog`.

**Not fixed here** -- this took the whole session's remaining budget to isolate and deserves a
clean-headed fix rather than a rushed one, especially after item 1's recession-cone attempt was
just refuted by an analogous "trusted a plausible-looking answer without a second probe" mistake.
**What a fix needs:** never trust `ef==-2` alone as "empty" -- either verify emptiness
independently (e.g. `region.maxLinear` on the SAME `A,b` with an arbitrary objective, which uses
`linprog` and is not subject to this Hessian-driven failure mode) before accepting the vacuous-true
branch, or detect the "unbounded objective, degenerate Hessian" case directly (Q's null space
intersected with the region's recession cone, checking `L'd` there) and refuse rather than
guess when `ef` is not unambiguously convergent. Either way: reproduce this EXACT case
(`fold4_objSplit.mat`'s cells 12, 17, 21) as the regression test before touching the code, and
re-run `certifiesNonPositive` with `Display,'iter'` after any fix to confirm `quadprog` itself
now agrees the region is not empty rather than merely patching the interpretation.

## 2026-08-27 (overnight, G1/G10) — the parked assembly diff, MEASURED on the full case-21 `conj` call

Re-applied `.claude/assembly_attempt_2026-08-25.diff` (the `collapseTinyEdges` pass plus
`matchHalfEdges`'s sagitta test) and ran `checkConjAgainstDefinition`'s case 21 end to end
(`q.conj('cplq')`, the actual reproducer G4/G10 name), since `clipArcByHalfPlane`'s bisection
backstop -- the one piece TODO.md said was still needed before re-applying the diff -- already
landed on 2026-08-25 and is in the tree.

**Baseline (diff NOT applied), for comparison:** 2.4 s, `PLQ:conjCPLQ:foldDroppedACID` at
`(-10,0)`, `47.10181578` against `10.86895777` -- exactly TODO.md G4's documented numbers,
confirming the reproducer is right.

**With the diff applied:** 292.5 s (120x longer), and a DIFFERENT failure --
`PLQ:conjCPLQ:cplqFailed`, an internal MuPAD error
(`mupadengine.evalin2sym: Invalid return value 'undefined'`) inside Case C's symbolic fallback,
not the fold cross-check. **This is not a clean win.** The diff changes `maxQuaPar`'s assembly
enough that the fast `foldDroppedACell` refusal (2.4 s) no longer fires on this fixture -- the
numeric route runs much longer before eventually declining differently, or not declining at all
and falling through to Case C, which then hits an UNRELATED, pre-existing bug in the symbolic
pipeline (not documented anywhere in `DECISIONS.md`/`TODO.md`/`SUPPORT_MATRIX.md` before now).

**Reverted** (`git checkout -- maxQuaPar.m`) -- the diff stays parked, unchanged from
2026-08-25, not committed. Nothing here is a regression (the file was untouched on `main`), but it
is also not the improvement the diff's own header claims for THIS specific measurement: the prior
session's "gets case 21 past assembly" referred to a narrower, assembly-only test, not the whole
`conj` call, and the two do not agree on this fixture.

**What this adds, for whoever picks G1/G10 up next:** the diff is not simply safe to land as-is;
turning it on changes the FAILURE MODE on at least one fixture from a fast, well-understood,
documented refusal into a slow, unrelated, previously-unseen one. Any future attempt to land it
should first check case 21 (and the other reds it was measured against) do not regress from "fast
named refusal" to "slow unrelated crash" -- that is now a concrete acceptance test, not just "does
assembly complete."

**A separate, smaller finding worth its own line:** Case C's `mupadengine.evalin2sym` failure on
this fixture (with the diff applied) is new information -- Case C was known to not finish in 25
minutes on this input (TODO.md G4), but not previously known to fail outright with an internal
MuPAD error. Not investigated further tonight; it is downstream of a change that is not landing.

## 2026-08-27 (item 3, FIXED) — the arc-parametrization fix, with a second gap found and closed before landing

Item 1's recession-cone attempt was refuted the same day (see that entry) by a direct
counterexample: a recession cone cannot see growth ALONG a curved arc. The correct ingredient it
prescribed -- parametrize the conic's own unbounded branch and check `cRow` composed along it --
was built, and validated the same way the first attempt should have been: hand-derive a
counterexample BEFORE trusting a large random sweep.

**A second gap, found on paper before writing any code.** The arc-only check answers "does cRow
grow along the CURVE", but a region's boundary is straight edges PLUS one arc, and the true
maximiser can sit on a STRAIGHT edge instead. Constructed by hand:

    Region = {s1 <= s2^2, s1 >= s2, s2 >= 0}     (both the arc s1=s2^2 and the straight edge
                                                    s1=s2 reach infinity as s2 -> infinity)
    cRow = (-1, 2)

Along the straight edge (s1=s2): `cRow.z = -s2+2s2 = s2`, UNBOUNDED. Along the arc (s1=s2^2):
`cRow.z = -s2^2+2s2`, bounded above (critical point at s2=1, value 1). An arc-only check would
find Ac<0 (bounded) and wrongly report the region's sup as ~1, missing the straight edge's genuine
unboundedness entirely.

**Fix: combine two independent sufficient conditions for unboundedness**, checked before ever
tightening:
1. a straight recession ray -- admitted by every linear facet AND by the conic's own recession
   condition (`d'Qd<0`, or the tie `d'Qd==0` with `grad.d<=0`, valid for this codebase's rank<=1
   parabola facets) -- with `cRow.d>0`. Sound as a sufficient trigger even though (per item 1's
   refutation) it is NOT sufficient for concluding the converse.
2. arc growth via `parabolaArcFrame`'s parametrization, exactly as item 1 prescribed.

If NEITHER fires, the region is bounded in that direction, and the true max is the best of the
arc's own finite candidates (endpoints/critical point) and the region's genuine finite vertices --
reusing the SAME candidate logic the bounded-case (st==0) branch already had.

**A third, more subtle gap, caught by the validation harness itself.** The FIRST version of this
fix's own "is the arc's admissible u-range empty" check returned `st=-1` (empty region) when it
came up empty -- but an empty ARC does not mean an empty REGION: a convex region can lie entirely
on one side of the conic, never touching the curve, and still be a perfectly good nonempty
polygon. Caught by a large random sweep (many `closed=-Inf, bruteforce=<finite>` mismatches).
**Fixed by ABSTAINING** (keep the caller's original val/st) rather than concluding emptiness --
sound because, by construction, this code only ever runs after the caller's own `region.maxLinear`
has already ruled out `st==-1`.

**Validated** against a true 2D brute-force oracle (the arc, every straight edge parametrized and
clipped, plus a far-field grid -- the FIRST oracle version only sampled the arc, which is exactly
why it missed the straight-edge gap above; rebuilding it properly was itself part of closing that
gap): 5000+ cases at the prototype level (raw numeric, both `+` and `-` conic orientations), plus
~150 built as real `region` objects through the actual `maxAffineOverRegion`. Zero genuine
disagreements; the few flagged were either the oracle's own sampling threshold (values a few times
below the reporting cutoff, still clearly growing) or a PRE-EXISTING, unrelated fragility in
`region.linearForm` on degenerate random test input (the error fires before this fix's code is
ever reached, on the unmodified call path).

**Landed**: fast 309/0, normal 12/0, no regressions. `TODO.md`'s scaling-defect entry has the one
remaining open question: whether this closes the ORIGINAL 374-refusal fixture's cell-count
blowup, which is too expensive to re-measure quickly and has not been re-run yet -- this fix is
validated on the MECHANISM (`maxAffineOverRegion`'s correctness), not yet on that specific
fixture's cell count.

## 2026-08-27 (overnight, item 3 vs the ORIGINAL fixture) — measured against the real scaling target: no material win

Item 3's fix was validated on a minimal witness and a hand-derived counterexample, and against a
brute-force oracle on thousands of synthetic cases -- but never against the actual fixture the
scaling defect was named on. Measured it overnight with `.claude/step3cost.m` (`x*y` over
`conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}`, the A.4/A.5 quadrilateral):

    FOLD 1: paired= 14 -> cells= 12   (baseline: paired=14 -> cells=12, unchanged)
    FOLD 2: paired= 24 -> cells= 23   (baseline: paired=29 -> cells=28)
    FOLD 3: paired= 46 -> cells= 36   (baseline: paired=42 -> cells=43)
    FOLD 4: paired= 56 -> cells= 50   (baseline: paired=60 -> cells=44)
    FOLD 5: paired= 64 -> cells= 58   (baseline: paired=63 -> cells=56)
    TOTAL 2794 s (baseline: 2944 s / "49 min")

**Final cell count: 58, against the baseline's 56 -- no improvement, and if anything two cells
worse.** distinctF stayed at 8 both times (the true minimum the answer needs), so the surplus
(~7x) is essentially unchanged.

**The refusal-reason tally moved, though not enough to matter.** Summed across all 5 folds:
`quadFacet_exactAnotInB` fired 139 times this run (of 2010 total merge attempts, 6.9%) against
the previously-recorded 374 (of 3236, 11.6%) -- a real drop in both count and proportion, and
`okQuadFacet` successes are up too (11 here). But the intermediate PAIRED counts also differ
(e.g. fold 4: paired=56 here vs paired=60 before), so this is not a clean apples-to-apples
before/after on identical inputs -- Step 1/2's own numeric route, which the fix also touches
(`maxAffineOverRegion` is called from more places than just this merge), plausibly produces
slightly different intermediate cells upstream, so the two runs are not diffing the exact same
arrangement fold by fold.

**Conclusion: the fix is a genuine, validated correctness fix (confirmed sound by an independent
brute-force oracle and a hand-derived counterexample), and it is NOT what closes this scaling
defect.** Most likely reason, from the fix's own documented scope: `tightenUnboundedFacet`
requires EXACTLY ONE curved facet (`numel(qidx) ~= 1` aborts) -- a real, deliberate restriction,
not a bug, since the argument's proof needs the single-parabola structure. If a meaningful share
of this fixture's surplus cells carry TWO OR MORE conic facets (plausible: six pieces folded
pairwise, so an accumulated cell can inherit conics from multiple ancestors), those pairs are
still refused exactly as before, no matter how correct the one-facet case now is.

**Kept the fix** -- it is a real bug fix, independently justified, and the fast/normal suites are
green with it in. **The scaling defect itself stays OPEN**, and the next step for it is
measuring how many `quadFacet_exactAnotInB` refusals on THIS fixture actually carry 2+ curved
facets (instrument `unionIsExact` directly, not `mergeTally`'s aggregate reasons) before
attempting to extend the argument to that case -- which is a materially harder proof (two conics'
combined recession behaviour, not one).

## 2026-08-27 (overnight, item 3 follow-up) — CONFIRMED: 60% of the remaining refusals are the two-conic case, out of the fix's scope

The previous entry's hypothesis, measured directly. Instrumented `unionIsExact`'s bare
`exactAnotInB` branch (a TEMP probe, reverted after this measurement -- `git checkout -- region.m`
immediately after, nothing shipped) to count `objA`'s curved-facet count on every refusal, on the
same fixture (`x*y` over `conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}`):

    probeMultiQ_nq1: 55    (39%) -- exactly one curved facet: tightenUnboundedFacet's target
    probeMultiQ_nq2: 84    (60%) -- two curved facets: OUTSIDE the fix's scope, refused as before
    probeMultiQ_nq3:  2    ( 1%) -- three: also outside

**Confirmed.** The majority of what is left (60%) is the harder case the fix deliberately declined
to handle (`tightenUnboundedFacet` aborts outright when `numel(qidx) ~= 1`), which is exactly why
tightening the one-conic case barely moved this fixture's final cell count (58 vs 56 baseline,
previous entry). The one-conic slice (39%) is real and IS being helped (`okQuadFacet` successes
rose, `quadFacet_exactAnotInB`'s overall RATE dropped) -- it is just not the majority share on this
particular fixture, where cells accumulate conics from multiple ancestors after several folds.

**What a two-conic extension would need, sketched but not attempted.** The single-conic argument's
two mechanisms (a straight recession ray, or growth along ONE arc's own parametrization) do not
compose directly for two conics: the region's TRUE recession cone is the intersection of both
conics' individual recessive-direction sets (not their union, and not decidable from either alone),
and there is no single global parameter analogous to `parabolaArcFrame`'s `u` once two curves are
both active constraints -- the boundary can now alternate between the two arcs, not just one arc
plus straight edges. This is a materially different and harder proof, not a small generalisation of
the one just built, and should be scoped as its own item rather than folded into the existing one
if it is picked up.

## 2026-08-27 (overnight, G17 end-to-end) — abandoned after 12 hours; the isolated reproducer stands as verification

`g17_full` (the full `rectMaximumIsTheConjugateOfTheWholeDomain` fixture, `.maximum` end to end,
with the `certifiesNonPositive` fix in place) was left running overnight to confirm the hole at
(-0.5,2) actually closes, beyond the isolated three-cell reproducer already confirmed. **Killed
after 12 hours of continuous CPU-bound execution** (process start 22:14:17, CPU time 43344s at
kill time -- essentially 100% of one core the whole time, not starved, genuinely still computing).
That is a wall-clock cost this fixture never asked for before: the UNFIXED pipeline reached its
(wrong) answer in 1562s (~26 min). The fix trades a fast wrong merge for real verification work at
every affected merge, and this fixture apparently has enough of them that the full pipeline no
longer finishes in any reasonable overnight budget.

**Not a red flag about correctness** -- every other check passed: the exact isolated reproducer
(cells 12, 17, 21 from the cached fold-4 state) now correctly refuses the bad merge and still
covers the point; fast (309/0) and normal (12/0) are green; the item-3 sweep found zero unsound
answers across thousands of cases. This is a COST finding, not a correctness one, and it belongs
with item 3's own follow-up (this session's other finding): if 60% of the remaining merge
refusals on the OTHER scaling fixture need the two-conic case, the same is plausible here --
`certifiesNonPositive`'s newly-honest verification may simply be running the two-conic path's
existing (already correct, already slow) machinery far more often than it used to, now that it is
no longer short-circuited by a wrong quick answer.

**Left as `UNVERIFIED` at the full-pipeline level**, correctly reported as such rather than
claimed. Re-running this specific fixture end to end is not worth another overnight slot on its
own; if it matters, instrument it the same way as item 3's follow-up (count nq per refusal on
THIS fixture) before deciding whether to spend more time on it, or accept the isolated-reproducer
level of verification as sufficient given the mechanism is otherwise fully confirmed.

## 2026-08-28 — the two-conic straight-ray extension MEASURED zero effect on the scaling fixture, exactly as its own commit message predicted

This session's `tightenUnboundedFacet` extension (mechanism 1, straight recession ray, generalized
from one curved facet to two) landed with an honest caveat: "most remaining two-conic refusals are
bounded regions needing a genuine finite tightening, which this does not attempt." Measured
directly with `.claude/step3cost.m`'s full 5-fold run on the A.4/A.5 quadrilateral (the same
fixture 2026-08-27's "item 3 vs the ORIGINAL fixture" entry used):

    FOLD 5: paired= 64 -> cells= 58   (IDENTICAL to the item-3-only, single-conic-fix numbers)
    quadFacet_exactAnotInB summed over all 5 folds: 139   (IDENTICAL to the 139 recorded 2026-08-27)
    TOTAL 2546 s   (against 2794 s on 2026-08-27, and 2944 s / "49 min" on the original baseline --
                    a real speed improvement from OTHER work since, not from this extension)

**Zero change to cell count or refusal tally, exact match to the pre-extension numbers.** This is
the confirmation the caveat predicted, not a surprise: `tightenUnboundedFacet` only runs when the
LINEAR-only relaxation says unbounded, and mechanism 1 only fires when a genuine straight ray
recedes every facet -- a condition that is vacuously false on every region that is actually
BOUNDED, which is what most of this fixture's accumulated two-conic cells are. The fix is sound
(regionTest's dedicated case still pins it) and correctly scoped; it simply does not touch the
refusal population THIS fixture has. **Confirms, rather than refutes, the commit message's own
disclaimer** -- recorded so nobody re-measures this expecting a different answer without a new
mechanism (the genuinely open two-conic BOUNDED tightening, still unattempted).

## 2026-08-28 — the parked assembly diff's downstream MuPAD crash is the ALREADY-KNOWN Step 3 gap, not a fresh bug

2026-08-27's entry flagged the diff's `mupadengine.evalin2sym` crash on case 21 as "new information
... not investigated further." Investigated now, without re-landing the diff: applied it
uncommitted (`git apply`), re-ran case 21's `q.conj('cplq')` directly (not the whole suite), read
the full error with its own message text this time instead of just the identifier.

**It is not a fresh bug.** The error's own text: *"Step 2 fell back to cPLQ for a rational envelope
face, and cPLQ's symbolic pipeline failed at `mupadengine.evalin2sym`... Step 2 itself is known to
complete on these envelopes; what is not yet reliable is Step 3, cPLQ's cross-piece maximum
(`plq.maximumConjugate` -> `functionNDomain.maximumP` -> `region.maximum`). See `SUPPORT_MATRIX.md`
section 1.2."* That is `conjEnvelopeViaCPLQ`'s OWN pre-existing, documented error wrapper -- the
diff does not introduce a new failure mode, it just makes THIS fixture reach an old, already-named
one (the legacy Step 3 `region.maximum`/`mergeL`/`isAlways` pipeline, the exact machinery
`.claude/step3cost.m` profiles) instead of exiting earlier via `foldDroppedACell`.

**Consequence: the parked diff still has no clean path to landing**, for a reason distinct from
2026-08-27's "not a clean win" -- it isn't trading a known failure for an unknown one, it is
trading a FAST known failure for a SLOW instance of a DIFFERENT already-known failure (SUPPORT_MATRIX
1.2's Step 3 unreliability). Landing it would need Step 3's legacy reliability gap closed first,
which is a separate, larger, already-tracked item, not a consequence of the assembly diff itself.
Reverted (`git checkout -- maxQuaPar.m`); nothing committed. Acceptance test
(`sweepCase21FailsFastAndNamedNotSlowAndUnrelated`) re-confirmed green at baseline afterward.

## 2026-08-28 (later) — mechanism 3 (different-axis convex conics) ALSO measures zero on the scaling fixture

Landed a genuine, separate closed-form fix today (`region.m`'s `tightenUnboundedFacet`): two
CONVEX (PSD) curved facets on different axes are provably bounded, closing a real `Inf` bug
(`region([y^2-x,x^2-y],[x y])`, true max 0.25, previously answered `Inf`). Measured its effect on
the actual scaling fixture the same way as mechanism 1 was measured (`.claude/step3cost.m`, full
5-fold run): **paired=64->cells=58 and `quadFacet_exactAnotInB` sums to 139 -- IDENTICAL to both
the pre-mechanism-3 AND the pre-mechanism-1 numbers.**

**So this fixture's 84 `nq2` (two-curved-facet) refusals are not the different-axis-convex case
either.** They must be same-axis pairs (mechanism 1 already resolves those when a straight ray
exists, but a same-axis pair need not HAVE one) or involve at least one concave facet (TODO.md's
Step-4 entry: "a piece whose arc is CONCAVE towards it" is a documented, real configuration).
Neither mechanism 1 nor mechanism 3 was ever expected to cover those -- this is confirmation, not
a surprise, and it sharpens what a future attempt needs: same-axis two-conic tightening (a
DIFFERENT closed form -- the shared axis makes a 1-D reduction possible along it, not attempted)
and/or the mixed-concavity case (where recession is generically permissive, per mechanism-1's own
analysis, and the real question shifts to the linear facets combined with ONE conic, closer to
the already-solved one-conic case than to a genuinely new two-conic proof).

Total wall time 2226 s against 2546 s last time -- likely shared-machine variance (AI/CLAUDE.md
sec 3), not attributed to this fix.

## 2026-08-28 (last) — the mechanism-1 sign-gap fix ALSO measures zero: three real fixes, three zero-effect measurements, and that is itself the finding

Generalizing mechanism 3 exposed a real, previously-latent bug in mechanism 1 itself (only one
sign of each conic's null direction was ever tested, since `quadNullDirsNumeric` returns one
representative and the candidate list never added its negation) -- fixed alongside the same-axis
generalization. Re-measured the scaling fixture a third time: **`paired=64->cells=58`,
`quadFacet_exactAnotInB` sums to 139 -- IDENTICAL to all three prior measurements** (pre-mechanism-1,
post-mechanism-1, post-mechanism-3).

**Three consecutive real, independently-verified bug fixes; three consecutive exact-zero
measurements on this one fixture.** Each fix is real and correctly scoped (regionTest's dedicated
cases prove each one on its own witness, cross-checked against brute-force sampling) -- this is
not a case of chasing phantom bugs. It means this fixture's 84 `nq2` refusals are consistently
outside ALL FOUR mechanisms tried (straight-ray, different-axis-convex, same-axis-opposite-sense-
convex, and the sign fix that should have made the first three exhaustive for their own scope).
**The remaining population is very likely CONCAVE-facet-involving** (mixed convex/concave, or
concave/concave) -- the one class none of today's mechanisms touch, per this session's own
analysis (concave facets recede almost everywhere, shifting the real question to the linear
facets combined with the convex partner, which needs its own argument, not attempted).

**Stopping the fixture re-measurement here.** A fourth run without a new mechanism to test would
add no information (AI/CLAUDE.md sec 5) -- this is the point to change what is being tried, not
to keep measuring the same thing.

## 2026-08-28 (final) — the ORIGINAL DIAGNOSIS was wrong: `exactAnotInB` refusals are correctly-decided TRUE violations, not proof-completeness gaps

Rung 1 (observe, no editing): instrumented `holdsOn` directly (temp probes, `git checkout --
region.m` after -- nothing shipped) at its TWO distinct failure points, since `unionIsExact`'s
`why='exactAnotInB'` covers both and 2026-08-27's original nq1/nq2/nq3 count never distinguished
them:

    branch 1 (st ~= 0, i.e. maxAffineOverRegion could not DECIDE a value)  -- HOLDSON_PROBE
    branch 2 (val decided but val > bc(i), a genuine VIOLATION)             -- HOLDSON_PROBE2

Ran fold 2 alone (14 `exactAnotInB` refusals, matching every prior measurement exactly).
**HOLDSON_PROBE fired ZERO times. HOLDSON_PROBE2 fired all 14, with real excess margins**
(3.12, 1.85, 0.93, 0.666, 0.196 -- not near-zero numerical noise) **and mostly `nq=1`, not
`nq=2`** (12 of 14 single-curved-facet, only 2 two-curved-facet).

**This overturns the premise behind three of today's fixes.** 2026-08-27's "60% of refusals are
nq2, outside the fix's scope" measurement counted curved-facet numbers on EVERY `exactAnotInB`
refusal without asking WHY `holdsOn` returned false -- lumping genuine geometric violations
(branch 2: A really is not contained in B', so refusing IS correct) together with proof gaps
(branch 1: `maxAffineOverRegion` could not decide, which mechanisms 1/3 target). On this fixture,
EVERY refusal in the sample is branch 2. There is nothing here for a better boundedness proof to
fix, one-conic or two-conic: `maxAffineOverRegion` already computes the right, decided answer,
and that answer correctly says the merge is unsound. This is not a new defect and needs no fix --
it is `unionIsExact` doing exactly its job, on cells that genuinely do not tile convexly.

**Why this was not caught earlier.** The 2026-08-27 measurement was real and its raw counts
(nq1=55, nq2=84, nq3=2) are not in question -- what was never checked, by that session or by
this one until now, is whether ANY of those refusals were actually undecided rather than
correctly decided. Three real, independently-verified bug fixes came out of chasing this anyway
(the different-axis, same-axis, and sign-gap fixes all stand on their own witnesses and remain
correct and committed) -- but none of them could ever have moved THIS fixture's count, because
the count was never measuring what it was assumed to measure.

**What this means for the scaling defect.** It stays open, and the real target is different from
what TODO.md's G4/item-3 entries have said since 2026-08-25: the ~7x cell-count surplus (58 cells
against distinctF's true minimum of 8) is not primarily a merge-CERTIFICATE gap -- it is that the
cells being generated upstream (Step 1/2's split, or `maxQuaPar`'s own fold construction) are
genuinely too fine an arrangement for `unionIsExact` to recombine, because many of them really do
not pairwise-union convexly. The next attempt should look UPSTREAM of `unionIsExact` -- at why
the arrangement has cells this fragmented in the first place -- not at the merge certificate
again.

## 2026-08-28 (very last) — the upstream cause, first witness: a genuine SLIVER cell failing against THREE different neighbours

One more rung-1 step (temp probe on `holdsOn`, geometry only, reverted after -- nothing shipped):
dumped the actual cell behind three of fold 2's `exactAnotInB` refusals. All three are the SAME
cell A --

    A.nv=4, vertices (0.572949,1.80902) (0.818182,1.72727) (0.0919554,0.954022) (0.25,0.875)
    A.ineqs: 4 linear + 1 quadratic (a genuine parabola, (s1+0.6 s2)^2 - 4.8 s1 <= 0)
    area 0.0379, diameter 1.061  -- a SLIVER: about 15x thinner than a "round" shape of the
    same diameter (a disc of that diameter has area ~0.56)

-- refused against three DIFFERENT candidate neighbours (three separate merge attempts), each
with a real, non-numerical-noise excess (0.196, 0.93, 3.125 on the violated bound).

**Working hypothesis, not yet proven:** this matches the SAME topological signature
`maxQuaPar.m`'s `assemblePieces` documents for its own (different, numeric) pipeline -- "a
genuinely AMBIGUOUS 3-way (or more) cluster... there is no valid pairwise resolution however the
candidates are chosen" -- except here in the LEGACY `region`/`functionNDomain` pipeline, and about
CELL SHAPE rather than vertex identity: a thin sliver sitting where 3+ neighbours meet may need a
genuinely N-ARY union to recombine correctly, which `unionIsExact`'s inherently PAIRWISE test
(one A, one B, checked in isolation) cannot express even in principle -- no single B' can contain
a sliver whose true covering is the union of several neighbours together.

**Not confirmed, and deliberately not chased further this session** -- confirming it needs: (1)
identifying ALL of A's true geometric neighbours (not just the 3 sampled here) and checking
whether their COMBINED union (not any one pairwise attempt) is what A's own boundary actually
needs; (2) if so, deciding whether an N-ary merge test is buildable, or whether the real fix is
upstream (whatever fold step produced this sliver in the first place, avoiding it at the source
per `AI/CLAUDE.md` sec 5's rung 5 -- reduce scope / prevent, rather than reconcile after the
fact). Either is real, scoped, multi-step work for a fresh session, not a continuation of
today's `region.m` fixes.

## 2026-08-28 (truly last) — the sliver hypothesis is REFUTED; the real mechanism is same-CURVE, not same-FACET, false candidates

Tested the multi-way-junction hypothesis above directly rather than leaving it as a guess:
dumped the FULL geometry (not just the one violated constraint) for every distinct candidate `B`
`region.merge` tried against the sliver `A` on its curved edge (`ia=3`) -- four distinct `B`
shapes across six attempts. Built a brute-force oracle (2,000,000 samples inside `A`'s own
bounding box, Python, independent of MATLAB) and checked: does `A` lie inside `B_a`, `B_b`, `B_c`,
`B_d`, or their UNION?

**Zero coverage. Of 799,743 samples landing inside `A`, exactly 1 also landed in any of the four
candidates, combined.** The N-ary hypothesis is REFUTED, decisively -- these are not slices of a
group that together reconstruct `A`'s true neighbourhood.

**Why, found by comparing bounding boxes: every one of the four candidates shares exactly ONE
VERTEX with `A` -- `(0.572949, 1.80902)` -- and none of the other three.** They are different
ARCS of the SAME underlying parabola (`region.merge`'s candidate test,
`obj.ineqs(mq1(i)) == -obj2.ineqs(mq2(j))`, checks EQUATION equality only), meeting `A` at that
one point and extending away from it in unrelated directions (`B_c`/`B_d` run out to
`x approx -10.6, y approx 11.5`). They were never going to merge with `A` -- not because of a
missing N-ary proof, but because they do not share a boundary SEGMENT with it at all, only a
point on the same infinite curve.

**This is the real, mundane, and actionable finding.** `region.merge`'s "same conic equation"
match is a NECESSARY but not SUFFICIENT condition for two pieces to be genuine neighbours: many
pieces can lie on one shared parabola without ever touching along an arc. Every one of today's
`quadFacet_exactAnotInB` refusals sampled (6 of 6, this witness) is `unionIsExact` correctly
rejecting a candidate that was never a real neighbour, not a proof gap of any kind, one-conic or
N-ary. **A cheap pre-filter -- do `A` and `B` actually share two consecutive vertices (an edge),
not merely a vertex or the same curve equation -- would eliminate these attempts before paying
for an LP-backed `unionIsExact` call at all.** Whether it would also shrink the final CELL COUNT
(58 vs. `distinctF`=8) is a separate, still-open question: correctly-refused pairs that were
never real neighbours say nothing about whether `A`'s TRUE neighbour (whichever piece actually
shares its full curved edge) was tried and also failed, or was never generated as a candidate at
all. That is where a future session should look first -- not at `unionIsExact`'s certificate,
which this session has now shown twice over is not the defect.

## 2026-08-28 (correction) — the edge-adjacency pre-filter this session proposed is UNSAFE, and there is a direct precedent

Before implementing the "cheap pre-filter" the previous entry recommended (skip `unionIsExact`
when the matched curved facet's two operands do not share two matching endpoints, only a bare
equation), checked `region.merge`'s own HISTORY comment first -- and it names this exact class of
mistake already made and reverted. `quadCutsOther` (removed 2026-08-17) "refused when one
region's conic met the other anywhere but at a vertex" -- a DIFFERENT-looking but SAME-CLASS
heuristic (a plausible local-geometry pre-check standing in for the real certificate) -- and it
was wrong: "two adjacent cells each carrying a different parabolic arc ELSEWHERE have a perfectly
convex union and were refused by the first outright", costing 322 of 612 attempts on the control
fixture. The file's own conclusion: `unionIsExact` is "the EXACT criterion... so dropping a
conservative pre-check cannot make a wrong merge happen" -- which cuts both ways: ADDING one back
CAN.

**Why my proposed filter is not provably safe.** `unionIsExact`'s soundness argument (`A subset
B'` and `B subset A'`) is purely ALGEBRAIC -- it does not require the matched facet to be a
genuinely shared boundary SEGMENT, only that the inequality removed from each side is the same
equation with opposite orientation. Two operands whose curved facet only touches at one point
could, in principle, still satisfy the global containment conditions and merge validly, the same
way two cells sharing a DIFFERENT parabolic arc elsewhere can still union convexly. I have six
empirical samples where the arcs did not share a segment and `unionIsExact` correctly refused all
six -- that is evidence, not a proof, and `quadCutsOther`'s own history is a direct warning
against trusting six examples of a plausible-sounding geometric shortcut.

**Correction: do not implement the pre-filter as stated.** It would need either a proof that
sharing only one endpoint makes the algebraic containment conditions provably impossible (not
attempted), or it stays exactly what `unionIsExact` already is -- a call that is always made and
always correct, just possibly slower on candidates like this one. The actionable finding stands
(these six merges are correctly refused, and the sliver's TRUE neighbour, if any, is still
unidentified) -- the SPEEDUP idea does not, until someone proves it sound.

## 2026-08-28 (item 1, continued) — fold 2's 14 refusals trace to just TWO troublesome cells, and they look like each other's true neighbour

Captured EVERY `exactAnotInB` refusal in fold 2 (compact vertex-only dump, temp probe, reverted).
All 14 trace to exactly two cells:

    A  = (0.573,1.809)(0.818,1.727)(0.092,0.954)(0.25,0.875)   [the sliver, area 0.038]
    A2 = (0,0)(0.25,0.875)(0.092,0.954)

each repeatedly tried against the SAME small set of same-curve-different-arc false candidates
already identified (`B1`..`B4` from the earlier witness). **`A` and `A2` share TWO vertices --
`(0.25,0.875)` and `(0.092,0.954)` -- the exact pair that would make them genuine edge-adjacent
neighbours**, not just same-curve coincidences. Not confirmed whether they carry the same
function value or were ever tried against EACH OTHER (not captured by this dump, which only
logs the FAILING pair, not every attempted one) -- that is the concrete next check: if `A` and
`A2` have matching `.f` and were never paired against each other by `mergeL`'s same-function
grouping, that would be a real, narrow, checkable gap (not a proof gap, a GROUPING one -- two
cells with the same function value that `mergeL`'s `isAlways(simplifyFraction(...)==0)` test
failed to recognise as equal, for the usual symbolic-form reasons documented elsewhere in this
file). If they were tried and refused, the next question is why THAT specific pair fails
`unionIsExact` despite the apparent adjacency.

**Net effect: this fixture's whole 14-refusal count reduces to a two-cell, not a
twenty-three-cell, question.** Whoever picks this up should start by checking whether `A.f`
equals `A2.f` (one `isAlways` call) before doing anything else.

## 2026-08-28 (item 4) — the conic-conic closed form's PIECES already exist, unwired, one layer up the design

Checked before writing off item 4 as needing a from-scratch quartic implementation: `conicMeet.m`
and `ratQ.m` already exist, are tested (`conicMeetTest.m`, 12/0), and do EXACTLY what
`region.getVertices`' remaining `solve()` calls need -- the real intersection points of two
conics from an exact integer Sylvester resultant, Newton-polished and verified, with the
resultant itself kept as an exact certificate. **They are called from nowhere in the live
pipeline.** `grep -rln "conicMeet(" *.m` finds only `conicMeet.m` and its own test.

**Why they are not a drop-in for `region.m`.** Both require INTEGER-coefficient conics
(`ratQ.conic` asserts exact integers via `gcd`), and `region.m`'s facets are general `sym`
expressions that routinely carry surds (the A.4/A.5 split's whole reason for existing). More
fundamentally, `conicMeet`'s own header states its design premise: "vertices are stored as
intersections of conics with approximate coordinates as needed" -- i.e. it is built for
`QuaCon`, a NAMED-vertex H-representation (`CONJ_FIELD_PROOF.md` sec 6, `doc/QuaConExample.md`)
that does not yet exist as code (`QuaCon.m` is referenced by comment in three files and appears
only as design docs under `doc/`, never committed as a class). `region.m` stores vertices as
literal symbolic COORDINATES, used directly and exactly throughout (the same exactness this
session's other findings depend on) -- dropping a Newton-polished FLOAT into that model risks
exactly the "one ULP breaks merge" class of defect `ratQ.m`'s own header describes fixing
elsewhere.

**Corrected scope, from "needs a quartic implementation" (this session's earlier, too-pessimistic
read) to "needs QuaCon, or a proven-safe adapter into region's exact model, neither of which
exists yet."** Building `QuaCon` is the architecturally right answer per this codebase's own
design docs, and is a genuine multi-session rewrite (it touches how EVERY mesh class stores
vertices), not attempted here. Not implementing an adapter either -- the same soundness-first
discipline as the edge-adjacency finding above: dropping approximate coordinates into an exact
model needs a proof it cannot cost precision anywhere that matters, not a hopeful try.

## 2026-08-28 (item 1, resolved) — A and A2 are NOT the same function; the "adjacent sliver" theory is refuted, and this reframes the whole scaling defect

Checked the concrete question the previous entry raised (temp probe in `mergeL`'s own
function-equality test, reverted after): does `A` (the sliver) share a function value with
`A2` (its apparent edge-adjacent neighbour)? Caught all 4 times they were compared across
fold 2 -- **`sameF=0` every time.** They are genuinely different functions that happen to share
an edge. Not a merge failure; ordinary adjacency between two different pieces of the answer.

**This means `A`'s only SAME-FUNCTION partners in this fold are exactly the same-curve-
different-arc candidates already found (B1..B4)**, none of which share a real boundary segment
with it (2026-08-28 earlier entry, the brute-force refutation). So within fold 2, `A`'s function
genuinely has NO valid merge partner at all -- not because one was missed, but because none
exists yet: either the true completing piece has not been produced by this fold (arrives later,
if ever), or `A` is a genuinely separate connected component of that function's argmax region.

**Reframes the scaling defect's own premise, which has stood since 2026-08-25.**
`.claude/step3cost.m`'s header calls `distinctF` "the count the answer actually needs: any
excess over it is cells that ought to have merged" -- but a function whose true domain has (or
temporarily has, mid-fold) more than one connected convex piece breaks that equivalence. `A`'s
15x-thinner-than-round sliver, isolated from every other same-function cell sampled, is a
concrete instance where "cells > distinctF" does NOT mean "a merge was missed." **Whether this
generalises to most of the 58-vs-8 gap, or `A` is a rare outlier, is now the open question** --
and it is a question about the COUNT's own definition, not about `unionIsExact`,
`region.merge`'s candidate generation, or any boundedness proof. Recommended next step for
whoever picks this up: pick 2-3 more of fold 5's 50-cell surplus at random and check
whether each has a genuine, findable same-function neighbour it failed to merge with, or is
(like `A`) alone in its group -- that ratio is what tells you whether this is one outlier or the
main story.

## 2026-08-28 (item 2, attempted and reverted) — the parked assemblePieces diff has a SECOND regression beyond case 21

User authorized landing the parked diff and accepting case 21's known trade-off (fast refusal ->
the already-tracked Step 3 legacy crash, SUPPORT_MATRIX.md 1.2). Applied it, updated
`sweepCase21FailsFastAndNamedNotSlowAndUnrelated` to
`sweepCase21HitsTheKnownStep3LegacyGapNotANewFailure` (pinning `PLQ:conjCPLQ:cplqFailed`,
confirmed at 204 s), then ran the fast bucket as the standard regression check before committing.

**Found a SECOND, different, previously-uncharacterized regression: `conjEdgeLowerBoundTest`
2 failures + 1 incomplete**, on a fixture unrelated to case 21. `checkOrphanHalfEdges` now raises
`maxQuaPar:internal` ("piece 11 src [8 3] has no matching neighbour... closest candidate ...
dist 0.00143") where it previously didn't -- one test expects a refusal from a specific allowed
set and `maxQuaPar:internal` is not in it; the other test errors outright instead of completing.
0.00143 sits just ABOVE the diff's own `collapseTinyEdges` threshold (`tolPos = 1e-3`), so this
looks like the exact "genuinely small feature vs. arithmetic noise, ambiguous scale" territory
`TODO.md`'s G10 entry already flagged (median edge length ~2.987e-03, several edges near/under
`tolPos`) -- but on a DIFFERENT fixture than the one that analysis was done on.

**Reverted** (`git checkout -- maxQuaPar.m conjCPLQTest.m`) rather than deciding this is also an
acceptable trade-off -- the user's authorization covered the ONE known, already-documented case
21 trade-off specifically, not an open-ended "accept whatever else breaks." Confirmed clean at
baseline after revert (`conjEdgeLowerBoundTest` 5/0). **Not committed; nothing shipped.**

**For whoever picks this up:** the diff needs either a tolerance fix for edges just above
`tolPos` (risky without re-deriving why 1e-3 was chosen, and tolerance-tuning has its own bad
history in this file -- see `region.merge`'s HISTORY on `quadCutsOther`) or a real look at why
piece 11's boundary lands 0.00143 from its closest candidate instead of exactly on it. Re-run the
FULL fast bucket (not just the two fixtures already known) before concluding it's isolated to
this one.

## 2026-08-29 (item 1, generalized) — the whole fixture's 141 refusals concentrate into ~12 cell signatures, and the sliver persists across every fold to the end

Fold 2's finding (14 refusals, 2 cells) generalized to the FULL 5-fold run: a compact signature
probe (`nv`, first vertex, temp, reverted) over all 141 `exactAnotInB` refusals in the whole
fixture (14+43+26+58 across folds 2-5, matching the per-fold tallies exactly) found only
**12 DISTINCT cell signatures** account for every one of them:

    37x A[nv=4 v1=(0.25,1.4091)]      28x A[nv=5 v1=(0.57295,1.809)]
    25x A[nv=4 v1=(0.57295,1.809)]    21x A[nv=4 v1=(0.25,0.875)]
    14x A[nv=4 v1=(0.33768,1.5591)]    8x A[nv=3 v1=(0,0)]
    (6 more, 1-2 occurrences each)

**The top 5 signatures account for 125 of 141 (89%).** And the original sliver
(`v1=(0.57295,1.809)`) is not a fold-2 curiosity: it reappears with `nv=4` (25x) AND `nv=5` (28x)
-- meaning a cell rooted at that same vertex persists, unmerged, from fold 2 all the way through
fold 5, repeatedly tried against the same family of same-curve-different-arc false candidates
and never once finding a real partner.

**This is the closing finding for item 1, and it is good news for tractability.** The 58-vs-8
surplus is not 50 independent small problems -- it is dominated by a HANDFUL of specific,
persistently-stuck cells (best guess ~5-6 underlying cells, since some signatures are the same
cell surviving into a later fold with one more vertex from an intervening split). Whoever
attempts a real fix should trace these ~12 signatures' actual histories (which fold created each,
what its true intended neighbour would have to be) rather than treating this as a broad,
diffuse problem needing a general redesign. `unionIsExact`, candidate generation, and the
boundedness proofs have all been independently checked and cleared this session -- what remains
is specific to these few cells.

## 2026-08-29 (item 1, root cause) — the persistent sliver sits on a genuine HIGH-DEGREE HUB VERTEX, present before any folding at all

Traced the persistent sliver's vertex `(0.572949, 1.80902)` back to its origin: a standalone
script running just Steps 1+2 (no folding) on the quadrilateral fixture found this EXACT point
already shared by **4 separate cells in piece 1's own Step-2 output and 4 separate cells in
piece 2's**, independently -- 8 distinct cells meeting at one point before fold 1 even runs.

**This explains the whole pattern, not just this one witness.** `.claude/step3cost.m`'s own
header says `f*` of `x*y` over a convex quadrilateral has "a cone per vertex" -- so a dual point
like this one is very plausibly one of those cone apexes, where every piece whose primal
triangle touches the corresponding primal vertex contributes a cell that meets there. With 6
pieces (the A.4/A.5 split of a 4-vertex quadrilateral), a single hub vertex can accumulate cells
from MOST of them, and each pairwise fold (`region.merge`'s "same conic equation" match) finds
many candidates that touch at that one point without sharing a real edge -- exactly the
same-curve-different-arc false positives measured all session. **The refusal volume is a
structural consequence of PAIRWISE folding near a high-degree vertex, not a missing merge**: as
each new piece folds in, it re-tests the accumulating fan against every existing cell that
happens to touch the same point, and most of those tests were never going to succeed.

**This reframes the "fix" question one more time.** It is not "find cell X's missing partner"
(there may be none to find near a genuine fan vertex -- the correct decomposition there may
legitimately need several small pieces, one per angular wedge). It is "does the pairwise-fold
architecture generate MORE candidate pairs near a hub vertex than the final answer needs", which
is a question about the FOLD STRATEGY (perhaps: group cells by which hub vertex they touch and
resolve each hub's fan once, rather than re-testing it against every incoming piece) -- a
genuine architectural change, not a bug fix, and out of scope for anything short of a dedicated
redesign effort. Not attempted further; this is the natural stopping point for item 1's
diagnosis.

**Confirmatory check, 2026-08-29 (final):** worried this might instead signal a genuine
COVERAGE defect (a same-function group split by OTHER functions' territory sitting wrongly in
between) rather than the benign "correctly fragmented, never reconciled" reading above. Checked
the bounding box of the sliver's confirmed same-function siblings (B1-B4): it spans a huge area
(x in [-10.6, 0.82], y in [0.9, 11.5]) against the sliver's own tiny 0.038-area footprint --
consistent with one large, genuinely convex argmax region (plausible for a dominant face's
conjugate over much of the dual), fragmented by the incremental fold process. **Decisive point:
fold 5 is the LAST fold -- all 6 pieces are already in, and they still have not merged.** So
this is not "hasn't finished yet"; it is the final, persistent state, confirming rather than
contradicting the hub-vertex/different-arc diagnosis above, not a separate coverage bug.

## 2026-08-29 (landed) — `unionIsExact` memoized: a small, safe, MEASURED win, not the fold-strategy fix

Rather than the full fold-strategy redesign (still not attempted -- a genuine algorithm change),
found and landed a smaller, provably-safe optimization: `region.unionIsExact` is a PURE function
of its four arguments (nothing global or mutable feeds it), so caching by a canonical string key
cannot change any answer, only skip recomputing one already known. Verified the redundancy is
real before implementing: a partial (fold 1-3) run with a temporary key-logging probe found
**30 unique keys out of 46 total calls -- 34% exact duplicates**, consistent with the hub-vertex
finding above (a same-function cell that fails against a same-curve-different-arc candidate in
one fold is tried against that SAME candidate again in every later fold, since neither side
changes and `mergeL` re-derives its groups from scratch each fold).

**Landed as `unionIsExact` (memoized wrapper, persistent `containers.Map` cache) calling the
unchanged original logic, renamed `unionIsExactCompute`.** One bug caught before it shipped:
first attempt called the new helper methods with STATIC syntax (`region.uieKey(...)`) though
they are instance methods in the same block as `unionIsExact` -- `regionTest` caught it
immediately (2 failures, `MATLAB:subscripting:classHasNoPropertyOrMethod`), fixed by calling
them as `objA.uieKey(...)` instead.

**Correctness: IDENTICAL to every prior run.** Full 5-fold measurement on the reference
fixture: `paired=64->cells=58`, every merge-tally reason count byte-identical to all previous
runs (this session's and earlier). The cache changes nothing about what gets computed, only
skips recomputing it.

**Speed: real but modest, not a breakthrough.** TOTAL 2186 s against the 2226-2546 s range of
prior runs, with each fold's `maximumP` a bit faster than its own prior measurements (fold 5:
428 s vs 428-503 s previously). Consistent with `unionIsExact` being only PART of a fold's cost
(the 2026-08-18 profile's own breakdown: `getVertices` 125.6 s, region ctor 84.3 s, `merge` 66.1 s
of one ~200 s fold) -- eliminating a third of its redundant calls saves a real slice, not the
whole fold. **This is not the fold-strategy fix** (TODO.md G4, DECISIONS.md 2026-08-29 "item 1,
root cause") -- that remains open and is what would actually change the surplus cell count, not
just the wall-clock cost of computing it. Verified clean: fast 312/0, normal 12/0, regionTest
27/0.

## 2026-08-29 (landed, second) — `getVertices` memoized too: same argument, bigger win

Re-profiled the reference fixture (folds 1-4, MATLAB's own profiler) with `unionIsExact` already
memoized to find what dominates now. `region.getVertices` topped the list (407 s, 406 calls) --
`unionIsExact` had dropped out of the top 20 entirely, confirming the first memoization worked.

**Same soundness argument applies, and is even cleaner here.** `region`'s ONLY properties are
`ineqs, nv, vx, vy, vars` -- `nv`/`vx`/`vy` are EXACTLY what `getVertices` computes, and nothing
else in the object can affect the answer. The function's own long-standing HISTORY comment
already says it "is called more than once on the same (already-populated) object"
(`removeTangent` re-invokes it after deleting ineqs). Verified the redundancy before
implementing: a 3-fold key-logging probe (temp, reverted) found **312 unique keys of 439 calls
(29% duplicates)**, one region's ineqs recurring 10 times.

Landed as a memoized `getVertices` wrapping the unchanged original body (renamed
`getVerticesCompute`), same pattern as `unionIsExact` -- this time with no bug on the first try
(learned from the earlier mistake, used instance-call syntax throughout).

**Correctness: IDENTICAL to every prior run, again.** Full 5-fold measurement:
`paired=64->cells=58`, every merge-tally count byte-identical to every run this session and
before. **Speed: a real, bigger win this time** -- TOTAL **2008 s**, against 2186 s with only
`unionIsExact` memoized and the 2226-2546 s range before either landed. Every fold's `maximumP`
dropped substantially (fold 2: 109 s vs 138-160 s; fold 4: 233 s vs 316-386 s), consistent with
`getVertices` being the LARGER of the two redundant-computation sources this session found.
Verified clean: fast 312/0, normal 12/0, regionTest 27/0, slow bucket 96/0 (543 s, against 612 s
with only `unionIsExact` memoized and 798 s originally).

**Still not the fold-strategy fix.** Both memoizations reduce WASTED recomputation; neither
changes what gets computed or the final cell count (58, unchanged). The actual surplus (58 vs
`distinctF`=8) is still open, per the hub-vertex diagnosis above.

## 2026-08-29 (landed, third) — `simplifyUnboundedRegion` memoized too: same pattern, another real win

Re-profiled again (folds 1-3) with both prior memoizations in place. `region.getVertices`'
CACHE HITS were confirmed working (127 of 406 calls skipped, cost concentrated in
`getVerticesCompute`'s 279 genuine misses) -- `simplifyUnboundedRegion` (167 s, 256 calls) was
the next candidate. **Same argument, verified before implementing**: it reads only
`obj.ineqs/nv/vx/vy/vars`, has no global state or side effect anywhere in its ~400-line body
(checked directly), and `region` is a VALUE class (`classdef region`, no `< handle`), so caching
a copy of its result is exactly as safe as a local variable. A 3-fold key-logging probe found
**161 unique keys of 256 calls -- 37% duplicates**, the highest rate of the three found this
session.

Landed the same way: memoized wrapper, unchanged body renamed `simplifyUnboundedRegionCompute`.

**Correctness: IDENTICAL to every prior run, third time.** Full 5-fold measurement:
`paired=64->cells=58`, every merge-tally count byte-identical to every run this session.
**Speed: another real, substantial win** -- TOTAL **1830 s**, against 2008 s (two memoizations),
2186 s (one), and the 2226-2546 s range before any. Cumulative from the original 2944 s
baseline: roughly 38% faster, with zero change to any answer at any step.

**Cumulative summary of this session's three memoizations** (all pure-function caches on
`region.m`, all verified byte-identical output at every stage): reference fixture TOTAL time
2944 s -> 1830 s; slow bucket 798 s -> 493 s. Verified clean at every step: fast 312/0, normal
12/0, regionTest 27/0, slow bucket 96/0.

`mtimes` (the fold's initial pairing step,
NOT touched by any of these three) is now comparable to or larger than `maximumP` in the later,
larger folds (fold 5: 417 s vs 260 s) -- a natural next place to look, though it is a
genuinely different kind of cost (new cross-products each fold, not obviously redundant in the
same way) and has not been checked.

## 2026-08-29 (checked, ruled out) — `region.plus` has ZERO redundancy; `mtimes`'s cost is genuinely irreducible by caching

Checked `mtimes`'s own dominant cost, `region.plus` (region intersection, called once per (i,j)
pair in `mtimes`'s O(m x n) loop), for the same redundancy pattern the other three had. Read it
first: it already has a cheap early exit (a direct symbolic-equality check finds a trivially
contradictory pair before ever building the combined ineqs list or calling the expensive
constructor), so any remaining memoization win would have to come from repeated (obj1,obj2)
pairs across DIFFERENT fold's `mtimes` calls.

**A 3-fold key-logging probe (temp, reverted) found 246 unique keys of 246 total calls -- ZERO
duplicates.** Confirms the reasoning that predicted this: `mtimes`'s inner loop tests every
(i,j) pair exactly once by construction, and across folds both operands change (the accumulated
result grows, and each fold introduces a genuinely NEW piece), so the same pair essentially
never recurs. **`region.plus` is NOT a memoization candidate -- this is real, unique,
irreducible-by-caching work**, unlike the other three functions. A different kind of
optimization (a cheap geometric pre-filter to skip pairs that provably cannot intersect, before
paying for the symbolic equality checks and constructor) is the only remaining lever here, and
it carries real soundness risk (must never reject a genuinely non-empty pair) that a pure cache
does not -- not attempted, and this session's memoization thread stops here, on a clean negative
result rather than a guess.

## 2026-08-30 — extra `mergeL` passes on the FINAL fold's output find nothing: confirms the hub fragmentation is genuinely structural, not order-dependent

Tested the cheapest possible fix for the 58-vs-8 scaling defect before attempting any redesign:
what if the pairwise merge just needs to be tried MORE times, or in a different pass, rather than
a real architectural change? `maximumP` already calls `mergeL` twice per fold; ran the reference
fixture to completion (58 cells, 8 distinct functions, matches every prior measurement exactly)
and called `mergeL` on the FINAL result two MORE times, outside any fold.

**Result: 58 -> 58 -> 58. Zero additional merges, both times.** This rules out an order-of-pairing
artifact (mergeL's grouping is greedy over increasing index; a different traversal order finding
new merges was the remaining untested explanation) and confirms [[item 1's root-cause]] finding
directly: the same-function cells around a hub vertex do not merge PAIRWISE in any order because
most pairs never share more than the single hub point as a boundary -- not a missing candidate,
not a merge-order problem, a genuine structural mismatch between what pairwise `region.merge` can
prove and what the hub's fan actually needs (an N-ARY simultaneous union of all wedges meeting at
the point, which `unionIsExact` was never designed to attempt two cells at a time).

**Why this closes off the "cheap fix" avenue for good.** Any fix must therefore either (a) merge
an entire hub's fan as one N-ary operation with its own soundness proof (unionIsExact's pairwise
argument does not generalize automatically -- a genuinely new, harder correctness question), or
(b) avoid the fragmentation upstream in Steps 1/2 so the hub's cells are never split into separate
wedges in the first place. Both are real redesigns, not increments; neither attempted this
session, consistent with the standing conclusion. `.claude/step3cost_extramerge.m`-equivalent
probe script not committed (trivial to rebuild: run the reference fixture to completion, then
call `mergeL` on the result once or twice more).

## 2026-08-30 (later) — G16's "wait for G1" precondition still not met: the 3-piece elliptical-edge witness dies on an UNRELATED, earlier gap

Checked whether G1 landing (2026-08-28, the `assemblePieces` fix) changed reachability for
`doc/QuaConExample.md`'s minimal 3-piece elliptical-edge witness (triangle `(0,0),(60,10),(-5,10)`
cut by two cevians, pieces 1/3 non-adjacent, both positive definite -- the exact, math-verified
case where `f*` needs a genuinely elliptical edge). Ran `q.conj('cplq')` on it directly.

**Identical failure to the one this doc already recorded 2026-08-20/21, unaffected by G1:**
`PLQ:conjCPLQ:cplqFailed` at the same dual point `s=(-121,-121)`, assembled max `NaN` vs the
per-piece pointwise max `0`. This is `SUPPORT_MATRIX.md` 1.2's tracked pre-existing gap (the same
family as sweep case 21's Step 3 fallback, [[G1/G4/G10 landed]]) -- a coverage/assembly defect at
the far corner of the initial probe grid, unrelated to G1's `matchHalfEdges`/sagitta fix. The run
never gets far enough to test whether `maxQuaPar` would hit `QuaPar:notParabola` on a genuine
elliptical edge; the input dies on this earlier gap first, exactly as before.

**Conclusion: G16 ("build `QuaCon` when a three-piece input reaches Step 3 with an elliptical
edge") is still not triggerable, but for a different reason than its own text states.** It is not
waiting on G1 specifically -- it is waiting on `SUPPORT_MATRIX.md` 1.2's general Step-3 assembly
gap, which is a separate, already-tracked, unresolved item. Do not re-check this precondition
again after any future G1-adjacent fix; check it after 1.2 closes instead. `.claude`-equivalent
probe script not committed (trivial to rebuild from `doc/QuaConExample.md` section 2's exact
coefficients).

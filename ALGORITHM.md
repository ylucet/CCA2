# The algorithm — what runs, in what order, and why

_Written 2026-08-18. This is the ORDER of operations, not a list of routines. It exists because
the ordering is the design: the same mathematics costs seconds or hours depending on when each
question is asked._

There are **two operators**, and they are not the same computation:

* `conj` — the Fenchel conjugate `f*`.
* `biconj` — `f** = cl co f`, i.e. **the closed convex envelope**. CCA2 has a DIRECT algorithm
  for this (Step 1 of [COAP]/[JOGO]), so computing it as `conj(conj(f))` is a detour.

---

## The governing principle

**Classify first, split last, and prefer the direct operator over a composition.**

Every expensive path in this toolbox comes from violating one of those three. The measured
example: `f = (x²+y²)/2` over the unit square is convex, so `f** = f` and there is nothing to
compute — yet it cost **436 s**, because the code split it into triangles, conjugated each,
took a symbolic cross-piece maximum, and conjugated the result back.

---

## `biconj(f)` — the convex envelope

    0. NORMALISE      merge adjacent faces carrying the same quadratic   [not yet built]
    1. f CONVEX       -> return f                        (nothing to compute)      [built]
    2. f SEPARABLE on a box  -> 1-D envelope per axis, sum                         [built]
    3. f BILINEAR on a box   -> McCormick / Al-Khayyal-Falk, in closed form        [built]
    4. SINGLE bounded TRIANGLE -> Step 1 directly        (no conjugation at all)   [built]
    5. otherwise      -> conj(conj(f))                   (the exception, not the rule)

Case 3 is the closed form `max( b(xl*y + yl*x - xl*yl), b(xu*y + yu*x - xu*yu) )` for `b > 0`,
and the other affine pair for `b < 0` (write `b*x*y = -|b|*x*y` and take `-|b|` times the CONCAVE
envelope, which turns the min into a max). The two pieces meet on a DIAGONAL of the box -- the
anti-diagonal when `b > 0`, the main diagonal when `b < 0` -- which is what makes the result two
triangular faces. It does not help SCIP, which applies McCormick itself; it helps a CCA2 user,
for whom `co(x*y)` over a rectangle is the most elementary nonconvex envelope there is and should
not cost a minute.

Case 4 is deliberately NOT widened to "any single bounded piece". `convEnvCPLQ` classifies by
Hessian: the PSD branch returns `q` over any convex `P`, but the negative-semidefinite and
INDEFINITE branches are stated for a TRIANGLE and raise `convEnvCPLQ:notImplemented` otherwise.
Widening was tried on 2026-08-18 and turned all four bilinear box rows into that error.

Cases 1 and 4 are the same short-circuit at two scales, and the convex half of it **already
existed inside Step 1** — `convEnvCPLQ` returns `q` unchanged as soon as the Hessian is positive
semidefinite. The problem was never that the check was missing; it was that `biconj` only reached
Step 1 for a TRIANGLE, so a box (`nv == 4`) fell through to the double conjugation. Case 1 now
catches it before the dispatch, without calling Step 1 at all.

Step 0 matters more than it looks, and is the one piece not yet built. `biconjugateTest` hands the
unit square in as two triangles sharing a diagonal; merged, it is one face, cases 1–3 become
reachable, and the piece-coupling that forces case 5 disappears.

Case 5 is genuinely needed only when the envelope COUPLES several pieces — the convex hull of a
union is not determined piecewise.

---

## `conj(f)` — the conjugate

    0. NORMALISE      merge adjacent faces carrying the same quadratic
    1. f SEPARABLE on a box -> 1-D conjugate per axis, sum
    2. per piece, classify by the SIGN OF THE HESSIAN:
         convex / affine -> envelope = f (co f = f);
                            cells from normal cones / KKT, on the polygon AS IT STANDS
         concave         -> envelope is affine, from the vertex values;
                            the conjugate is the support function
         indefinite      -> change variables to x*y, triangulate,
                            [COAP] A.2 / A.3 / A.4 by convex-edge count,
                            A.5-split when three edges are convex
    3. Step 3: max ACROSS pieces -- only if more than one piece survived

---

## Why the triangulation exists, and why it is not universal

The whole apparatus of triangulating, counting convex edges and selecting among A.2/A.3/A.4/A.5
is **the indefinite case's algorithm**. It is per-triangle because those closed forms are. No
other sign class needs any of it:

* **convex** — `co f = f`. Step 1 is the identity. Step 2 is the KKT active set (a cell per
  vertex, per edge, one interior), which works on a polygon with ANY number of affine facets.
  A convex piece never needed a triangle; splitting it only forces Step 3 to glue back together
  what was never broken.
* **concave** — `co f` is the affine interpolant of the vertex values. On a many-sided polygon
  the envelope really is piecewise affine over the LOWER HULL, so a split is genuinely needed
  there; but any fan stays CORRECT, because conjugation turns a union into a max
  (`f* = maxₖ (f|Tₖ)*`). It is suboptimal, not wrong.
* **separable on a box** — drops to 1-D entirely, below all of the above.

Separability is tested FIRST and specifically BEFORE the change of variables, because an
indefinite DIAGONAL quadratic such as `x² − y²` is separable as written, and rotating it into
the `x*y` frame destroys both the separability and the box. That ordering is the whole reason
`x² − y²` over a box went from `MATLAB:badsubscript` to exact.

Separability is a property of the FUNCTION AND THE DOMAIN TOGETHER. `x*y` separates in rotated
coordinates (`u = x+y, v = x−y`), but the box is not a box there, so that does not help and is
deliberately not attempted.

---

## Two invariants worth exploiting

1. **`f*` is always convex.** So in `biconj = conj(conj(f))` the second conjugation always has a
   convex input and must never enter the indefinite machinery.
2. **`conj(f) = conj(co f)`.** Conjugating the envelope gives the same answer, so when Step 1 is
   cheap it is always worth doing first.

---

## The `isConvex` flag — trusted, with one free guard

Convexity of a SINGLE quadratic piece is free to check (`eig` of a 2×2 Hessian). Convexity of a
MULTI-PIECE PLQ is not: it needs per-piece convexity plus consistency of the gradient jump across
every shared edge. The flag therefore earns its keep exactly where verification is expensive, and
it unlocks the largest short-circuit there is (`biconj` of a convex `f` is `f`).

The hazard is equally specific: a wrong flag makes `biconj` return a NON-convex `f` as its own
envelope, silently. Two mitigations, both cheap:

* **The free necessary condition is still checked.** Every piece's Hessian must be positive
  semidefinite. That is necessary, not sufficient, and costs microseconds — so a flag that
  contradicts it is refused LOUDLY rather than trusted.
* **Full verification under an opt-in assert**, in the established `MAXQP_ASSERT` /
  `QUAPAR_VALIDATE` style.

The flag lives on the function object, not in an argument, so it travels with the data.

---

## Rational split points — where the principle applies, and where it must not

Split points divide into **chosen** and **determined**, and the distinction is the whole story.

* **CHOSEN — keep rational.** The fan triangulation of a polygon is free: any triangulation
  works, because `f* = maxₖ (f|Tₖ)*`. Its vertices are the polygon's own, so it is already
  rational whenever the input is. The genuinely free points that are NOT yet rational by
  construction are elsewhere — notably the midpoint `M` picked between two arcs in the two-arc
  split, chosen for convenience rather than determined by any tightness condition.
* **DETERMINED — must not be moved.** A.4/A.5's split is the cevian along which cPLQ's closed
  form stops being tight. Move it and one sub-triangle straddles the tightness boundary, so the
  closed form there is a strict MINORANT rather than the envelope. That is precisely the defect
  the `nCE == 2` branch carried before it was fixed: a minorant reaching `−0.2835` where the
  truth is `≥ 0`, and a correspondingly too-large conjugate.

**Prior art, and it is recorded.** Attempt 3 of the general-quadrilateral fix did snap this split
to rationals — taking the geometry from double-precision `convEnvCPLQ` and snapping vertices to
the simplest rational within `1e-10`. `DECISIONS.md` records why it failed: bounding the VERTEX
denominators does not bound the downstream ones, because the conjugate is a rational function of
those coordinates, so a few squarings carried `1e5` to `1e25` and the run hung. Attempt 4
succeeded by doing the opposite — carrying the irrational foot as a compact surd
(`5/2 − √5/2`).

**The cost is bounded, not runaway.** A.4's foot generates a QUADRATIC EXTENSION, and the
measured price is `testcPLQ` 1542 s → 2188 s. That buys a correct answer where the default
previously crashed.

**The better lever for the same goal** is not to relocate the split but to not need A.4 at all: a
single bounded triangle conjugated through the numeric route (`conjBoundedPolygon`) was measured
exact at 7 of 7 probe points and never introduces a surd. Whether that route can serve the
2-convex-edge case generally is unmeasured and is the open question here.

---

## Measured effect of getting the order right

| case | before | after |
|---|---|---|
| `(x²+y²)/2` on `[0,1]²`, first conjugation | 2 triangles, 14 cells, Step 3, **456 s**, err 1.6e-4 | 9 cells, **4.2 s**, exact 10/10 |
| `x² − y²` on `[0,1]²` | `MATLAB:badsubscript` | exact, 6 cells, 1.9 s |
| `f*` of `x*y` over a general quadrilateral | 86 cells, 73 min | 60 cells, 43 min |

And `biconj` over a box, the whole of `checkBoxEnvelopeForSCIP`, start of session to end:

| box case | before | after |
|---|---|---|
| `x*y` | 62 s, no mesh | **0 s**, `QuaPol` |
| `3xy+7x-2y+5` | 49 s, no mesh | **0 s**, `QuaPol` |
| `x^2-y^2` | **ERROR** | **0 s**, `QuaPol` |
| `(x^2+y^2)/2` | 456 s, no mesh | **0 s**, `QuaPol` |
| `x*y` on a sub-box | 43 s, no mesh | **0 s**, `QuaPol` |
| `x*y` on a wide box | 42 s, no mesh | **0 s**, `QuaPol` |

Every row also went from a mesh-less `QuaParCPLQ` to a meshed `QuaPol`, which is what a consumer
can actually read (`RETURN_TYPE.md` records that change).

---

## The 1-D UNBOUNDED case -- deliberately not built, and why

`conjSeparableOverBox` requires FINITE bounds, because every piece it builds is anchored on an
endpoint value `phi(l)` or `phi(h)`. Dropping a bound is easy algebra and awkward representation:

* `a > 0` -- unboundedness makes it EASIER: the maximiser `t* = (sigma-d)/(2a)` always exists, so
  each dropped bound removes a clamp. On all of `R` the conjugate is `(sigma-d)^2/(4a)`, ONE
  piece, no breakpoints.
* `a = 0` -- the conjugate's DOMAIN degenerates. On `[l,inf)` it is a half-line; on `R` it is the
  single point `sigma = d`.
* `a < 0` -- on any unbounded interval the sup is `+inf`, i.e. the envelope is `-inf`. Degenerate,
  and it must be REPORTED rather than answered.

The representational hazard is the reason to wait: a point-domain in one axis crossed with a full
line in the other gives a 2-D cell with EMPTY INTERIOR, and `region` reasons about vertices. That
is exactly the shape this codebase keeps getting wrong -- `simplifyUnboundedRegion` declared a
half-plane empty for want of a finite vertex, fixed 2026-08-16.

It is also no longer NEEDED for the box cases: `biconj` now short-circuits convex, separable and
bilinear inputs, so none of them conjugates twice. It would still pay for a genuinely unbounded
separable input (`x^2 + y^2` over a quadrant), but that needs TWO changes, not one -- `f*` is
piecewise, so the detector would have to test separability PER CELL, not on a single quadratic.

The 9 cells for `(x²+y²)/2` are the product structure the answer actually has: `f` separates and
the box is a product, so `f*(s) = g(s₁) + g(s₂)` with `g` the 1-D conjugate in three pieces.
Triangulation destroys that structure and Step 3 pays to rebuild it.

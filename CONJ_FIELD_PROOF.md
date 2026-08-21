# The number field of `conj` and `biconj`: answered, NEGATIVELY

_2026-08-20. Written in a session with read-only access to the rest of the repo; this file is new
and referenced by nothing. No `.m` file, no other `.md` file, and no commit was touched._

## Answer in one line

**No.** For a **continuous** PLQ with rational data — the real input class, since a PLQ is a
function on a closed polyhedral subdivision and is therefore continuous — the output of `conj` is
**not** always representable over `Q(sqrt(d))`, nor over any iterated quadratic extension
`Q(sqrt d1)(sqrt d2)...`. An explicit **five-piece continuous** counterexample has a vertex of
`f*` of algebraic degree 4 with Galois group **S4** (order 24, not a 2-group). The same
counterexample refutes the claim for `biconj` / `convEnv`.

The failure is narrow and completely characterised, which is the useful part:

| stored datum of `f*` | field |
|---|---|
| face functions `f` (the quadratic coefficients) | **always rational** — Theorem 1 |
| edge conics `Ec` (the parabolas / lines) | **always rational** — Theorem 1 |
| vertices `V` | degree `<= 4`; **a 2-group (so a tower works) whenever two of the three meeting pieces are ADJACENT** — Theorem 3; **S4, so no tower, when all three are pairwise non-adjacent** — section 5 |

So the whole problem lives in the vertex list, and there are **two** distinct shortfalls there, one
much commoner than the other:

* a vertex is generically of **degree 4** — measured at 10 of 12 continuous configurations (8.0) —
  and every `exactQ` element `a + b sqrt d` has degree at most 2. So `exactQ` already cannot hold
  the *ordinary* vertex. A two-level tower (a nested radical) can.
* when the three pieces meeting at the vertex are **pairwise non-adjacent**, the group is S4 and
  **no** tower of any depth can hold it (section 4). That needs at least four pieces in the input,
  and the three cannot be a chain. (That bound is for the S4 *field* result only. The separate
  *shape* finding — that an edge of `f*` can be an elliptical arc, 7.4b — needs just two
  non-adjacent pieces, so **three** pieces suffice for it: `doc/QuaConExample.md` §2.)

`exactQ` is therefore the right idea in the wrong coordinate system *and* one level short of even
the common case. Section 8 lists the data structure options and what each one buys.

---

## 1. Scope, and one correction

`f` is PLQ with rational data in full generality: a polyhedral subdivision with rational vertices,
`f = q_k` on the closed piece `T_k` with `q_k(x) = 1/2 x'Q_k x + beta_k'x + gamma_k` rational, and
`f = +inf` off the union. Convexity is **not** assumed, `Q_k` may be indefinite, singular or zero,
and pieces may be unbounded.

**Continuity is not an assumption, it is forced.** The pieces are closed, they overlap on their
shared edges, and `f` is a function — so `q_k = q_l` on `T_k ∩ T_l` and `f` is continuous. An
earlier version of this note used a discontinuous three-piece example; that is a legal `QuaPol`
*data structure* but not a PLQ *function*, so it does not answer the question. It is kept in
section 6 because it is the minimal illustration of the mechanism, and because the data structure
does accept it.

Continuity is the whole difficulty: it is exactly what makes the three-piece construction
impossible (Theorem 3) and forces five.

The claim under test is the **coefficient field only**, not the output *shape*. The shape claim is
separately false and already known to be: `biconjCPLQ` Case B returns a `RatPol`, a genuine
rational function, not a quadratic.

"Representable" is taken in the weakest useful sense: **any** finite tower of quadratic extensions
counts as success. That is strictly weaker than `exactQ`'s single `d`, so refuting it refutes
`exactQ` a fortiori.

The criterion used throughout: **an algebraic number lies in some iterated quadratic extension of
Q iff the Galois closure of its minimal polynomial has degree a power of 2.** For an irreducible
quartic this is decided by the resolvent cubic: reducible ⇒ the group is `D4`, `C4` or `V4`, all
2-groups, so a tower works; irreducible ⇒ the group is `S4` or `A4`, orders 24 and 12, and no
tower can reach the root.

---

## 2. The coefficients never leave Q, and the vertices have degree at most 4

### Lemma 1 — every candidate value function of one piece is a rational quadratic

Fix a piece and write `g(s) = sup{ <s,x> - q(x) : x in T }`. Where the sup is finite it is
attained, and the maximiser is one of three things.

* **a vertex `v` of `T`** (rational): `g_v(s) = <s,v> - q(v)` — rational **affine**.
* **the relative interior of an edge** `e = [v,w]`, `d = w - v` (both rational). Along the edge
  `phi(t) = g_v(s) + t*L(s) - t^2*alpha/2` with `alpha = d'Qd` and `L(s) = <s,d> - <grad q(v),d>`,
  both rational. A stationary point needs `alpha != 0`, and then `t* = L(s)/alpha` and

      g_e(s) = g_v(s) + L(s)^2 / (2*alpha)

  — a rational quadratic whose quadratic part is **rank one**.
* **an interior critical point** (needs `Q` invertible): `x* = Q^{-1}(s - beta)` and

      g_int(s) = 1/2 (s - beta)' Q^{-1} (s - beta) - gamma

  — a rational quadratic of full rank.

All three families are rational because `Q, beta, gamma, v, d` are. Degenerate cases delete a
branch; they never make one irrational. **QED**

### Theorem 1 — `f*` is assembled from those, and nothing else

`f = min_k (f + I_{T_k})` pointwise, so `f* = max_k (q_k + I_{T_k})*`. A max **selects** among the
functions of Lemma 1; it never manufactures a new one. Hence

* every **face function** of `f*` is one of the rational quadratics of Lemma 1;
* every **edge curve** of `f*`'s subdivision lies in `{A = B}` for two of them, so it is a conic
  with **rational** coefficients;
* every **vertex** of `f*`'s subdivision is a point where at least three of them agree. **QED**

### Theorem 2 — degree bound

> Every vertex of the subdivision of `f*` has coordinates algebraic of degree **at most 4** over Q.

At a vertex three rational quadratics `A, B, C` agree, so the vertex lies on the two rational
conics `{A = B}` and `{A = C}`. Two conics with no common component meet in at most 4 points;
eliminating `y` gives a rational quartic. If they share a component the vertex lies on a rational
line and the degree drops to at most 2. **QED**

Convexity of `f*` is what makes the "three quadratics agree" reading correct: a convex function
restricted to an open region is convex, so every face function of `f*` is a **convex** quadratic
and near a vertex `f*` is locally the max of the three. The vertex is therefore exactly a base
point of the pencil, not a lower-dimensional accident.

The A.5 cevian foot `5/2 - sqrt(5)/2` in `DECISIONS.md` is a **degree-2** instance — the easy end
of a range that reaches 4.

### Theorem 4 — the SAME bound for the envelope

> Every vertex of the subdivision of `f** = conv f` has coordinates of degree **at most 4** over Q.

Assume `f` proper with `conv f` closed (bounded `dom f` is the clean case; see the caveat below).
Write `S*` for the subdivision of `dom f*`: by Theorem 1 its 2-cells carry rational quadratics
`g`, its 1-cells are arcs of rational conics, and by Theorem 2 its 0-cells have degree at most 4.

For a closed proper convex `f*`, `x in d f*(s)` iff `s in d f**(x)`, so the relative interiors of
the sets

      R(C)  :=  union of  d f*(s)  over  s in relint C ,      C a cell of S*

are disjoint and cover `dom f**`: they **are** the subdivision of `f**`. Enumerating them by the
type of `C` gives every cell, and with it every corner.

| cell `C` of `S*` | `d f*(s)` on it | cell `R(C)` of `f**` | its corners |
|---|---|---|---|
| 2-cell, `f* = g`, `Hess g = H` nonsingular | the point `Hs + a` | 2-cell, `f** = g*` quadratic; `R(C) = H·C + a` | `H·(corners of C) + a` |
| 2-cell, `H` of rank 1 (an edge branch `g_e`) | a point on a fixed line | **1-cell**: a segment | images of the extremisers of an affine `L` over `C` |
| 2-cell, `H = 0` (a vertex branch `g_v`) | the constant `v` | **0-cell**: the single point `v` | `v` itself — **rational** |
| 1-cell, arc where `g_1, g_2` meet | the segment `[grad g_1(s), grad g_2(s)]` | 2-cell, the ruled (`RatPol`) piece | `grad g_i(endpoints of C)` |
| 0-cell, vertex `p` where `g_1,g_2,g_3` meet | the triangle `conv{grad g_k(p)}` | 2-cell on which `f**` is **affine** | `grad g_k(p)` |

So every vertex of `f**`'s subdivision is one of exactly four kinds:

1. a vertex `v` of the **input** subdivision — **rational**, degree 1;
2. `grad g_k(p) = H_k p + a_k` with `H_k, a_k` rational and `p` a vertex of `S*` — a rational
   affine image of a degree-`<=4` point, so **degree `<= 4`**;
3. `grad g_e(t)` where `t` extremises a rational affine functional over a cell bounded by rational
   conics: `t` is either a vertex of `S*` (kind 2) or a **tangency** point, the solution of a
   rational linear equation together with a rational conic — degree `<= 2`;
4. an intersection of one of those boundary arcs with the boundary of `conv(dom f)`, whose edges
   are rational lines — a rational conic met with a rational line, degree `<= 2`.

Each boundary arc above is `grad g_i` (rational affine) applied to an arc of a rational conic,
hence again an arc of a **rational** conic, which is what makes kinds 3 and 4 rational-vs-rational
intersections. All four kinds are degree at most 4. **QED**

Two remarks.

* **Sharpness, and it is kind 2.** In the counterexample of section 4 the three points
  `grad g_k(p) = x*_k(p)` are the corners of an affine cell of `f**` and all six coordinates are
  degree 4 with Galois group S4 (section 5). So the bound is attained on the envelope side too,
  and by the same mechanism.
* **Extra incidences only help.** If four faces of `f*` meet at one point there are three
  independent rational conics through it rather than two; the system is over-determined and the
  degree can only drop.

**Caveat.** The proof is written for `conv f` closed with the sup attained — bounded `dom f` is
the clean case. Unbounded pieces add recession cells; their data (recession cones of rational
polyhedra, and the directions in which `f*` is `+inf`) are rational, so the same argument applies,
but that case analysis was not written out here and is not claimed.

### Corollary to Theorem 4 — why `f**` is RATIONAL on a POLYHEDRAL subdivision

Theorem 4's table also explains, in one line each, the two things that make `biconjCPLQ` Case B
return a `RatPol`.

**Why the subdivision is polyhedral.** Look at what the gradient maps can produce. An **affine**
face of `f*` (a vertex branch `g_v`) has constant gradient `v` — a single point, and a vertex of
`T`. A **rank-one** face (an edge branch `g_e = g_v + L^2/(2 alpha)`) has gradient
`v + (L(s)/alpha) d`, which sweeps a **straight segment along the edge `e` of `T`** — no matter
how curved that face is in `s`-space. So every cell of `f**` is the convex hull of points and
segments lying on the vertices and edges of `T`: a polygon. The parabolic arcs of `f*` are
flattened by the gradient map, which is why curvature in `s`-space does not become curvature in
`x`-space.

**Why the function is rational and not quadratic.** On a ruled cell `f**` is affine along each
ruling, but the ruling's slope moves. Eliminating the ruling parameter is what creates a
denominator. Worked completely for `q(x) = x1 x2` on `T = (0,0), (1,1), (2,0)`, where
`alpha = d'Qd = 2 d1 d2` is positive on the single edge `(0,0)-(1,1)`:

    faces of f*:   g_(0,0) = 0     g_(1,1) = s1+s2-1     g_(2,0) = 2 s1
                   g_e = (s1+s2)^2/4,  valid where t* = (s1+s2)/2 in [0,1]

    the arc {g_e = g_(2,0)} parametrised by t:   s(t) = ( t^2/2 , 2t - t^2/2 )   (polynomial in t)
    the ruling at t joins the vertex (2,0) to the point (t,t) on the convex edge
    so the ruled cell is the TRIANGLE (2,0), (0,0), (1,1)  -- straight edges

    on it,  f**(x) = <s(t), x> - f*(s(t)) = lambda t^2 , and with
        mu = lambda = (2 + x2 - x1)/2      nu = lambda t = x2        (both AFFINE in x)

        f**(x) = 2 x2^2 / (2 + x2 - x1)

a **quadratic numerator over a linear denominator** — a `RatPol` piece, not a `QuaPol` one.
Checked against `f**(x) = sup_s <x,s> - f*(s)` computed numerically at 12 points of the cell:
agreement `<= 1.7e-16`, and on the convex edge itself `f**(t,t) = t^2 = q(t,t)` as it must be.

**The numerator is QUADRATIC, never cubic.** The general elimination above gives

      f**(x) = [ A mu^2 + (B + c0) mu nu + (C + alpha) nu^2 + q(w) mu ] / mu

with `mu, nu` affine in `x`. The only thing that could produce a cubic is `nu * <s(t), d>` with
`<s(t), d>` quadratic in `t` — but `<s(t), d> = L(s(t)) + <grad q(v), d> = alpha t + const` is
**linear** in `t` by the definition of the arc parameter. So the numerator stops at degree 2. That
matches `RatPol.m` exactly ("the numerator is a quadratic and the denominator is a nonzero linear
function", citing [COAP] Step 4 and [JOGO] Prop. 1), and the `1x10` cubic basis in which `f` is
stored is a storage convention shared with `QuaPol`, not a claim about the degree.

**Scope — this is the SINGLE-piece case, and it does not carry over.** The elimination is linear in
the ruling parameter only because one of the two faces bounding the arc is AFFINE (a vertex
branch), which is what a single piece gives. When both faces are full-rank — which needs several
pieces — the same elimination is a **quartic** in the ruling parameter and `f**` is algebraic of
degree 4 on that cell, not rational at all. See 5.1.

The other ruled-cell type behaves the same way. When the arc separates two **edge** branches that
share a base vertex — the A.4 two-convex-edge case — the affine part cancels,
`(2s1+s2)^2/8 - (s1+2s2)^2/8 = (3/8)(s1+s2)(s1-s2)` is a pair of lines, and the elimination gives
`f** = (2/9)(x1+x2)^2` on `T = (0,0),(2,1),(1,2)`: a pure rank-one PSD **quadratic** with no
denominator at all, which is the A.4 result already recorded in `DECISIONS.md`. Both cases checked
against a numerically computed `f**`, agreement `<= 1.1e-16` (`degreeF2.m`).

---

## 3. What continuity buys: adjacent pieces can never do it

### Theorem 3 — an adjacent pair gives a DEGENERATE conic

> Let pieces `i` and `j` share an edge, so continuity forces `q_j - q_i = l*m` with `l` the affine
> form vanishing on the shared edge and `m` affine. Then the conic `{g_i = g_j}` is **degenerate**
> — its 3x3 matrix is singular, and it is a real pair of lines.

Writing `l(x) = u'x` (translate the shared edge line to the origin) and `m(x) = v'x + e`,

      Q_j = Q_i + u v' + v u' ,        beta_j = beta_i + e u

— both perturbations are *in the direction `u`*, which is what forces the degeneracy. Verified as
a **symbolic identity** in `adjProof.m`: with `Q_i`, `u`, `v`, `e`, `beta_i` all free symbols, the
determinant of the 3x3 conic matrix of `g_i - g_j` simplifies to exactly `0`. The determinant of
its quadratic part is `-(u x v)^2 / (4 det Q_i det Q_j) <= 0`, so the pair of lines is real, not a
conjugate complex pair. Independently sampled on 946 random rational adjacent pairs: relative
`|det|` at most `1.4e-15`, against `9.6e-2` for a control with the continuity constraint removed.

### Corollary — such a vertex is always in a quadratic tower

The pencil spanned by the two conics through a vertex of `f*` contains all three differences
`g_i - g_j`, `g_i - g_k`, `g_j - g_k`. A **degenerate rational member of the pencil is a rational
root of the resolvent cubic.** So if *any two* of the three pieces meeting at that vertex are
adjacent, the resolvent cubic is reducible, the Galois group is a 2-group, and the vertex **does**
lie in an iterated quadratic extension of Q.

This is exactly what the search saw. Of twelve feasible three-piece continuous configurations
(fan of three triangles, pieces 1-2 and 2-3 adjacent), ten had an irreducible quartic and **all
ten** had a reducible resolvent cubic. None is a counterexample; they are Theorem 3 in action.

**So an S4 counterexample needs three pairwise NON-adjacent pieces**, which needs at least four
pieces in the subdivision. Two non-adjacent pieces — three in total — are already enough for the
weaker *shape* statement of 7.4b, since that needs one non-degenerate edge rather than a triple
point (`doc/QuaConExample.md` §2).

---

## 4. The continuous counterexample: five pieces, S4

### 4.1 The input

Domain: the triangle `(0,0), (60,10), (-60,10)` cut by four cevians from the apex, giving five
triangles fanned from `(0,0)`. Active pieces are **1, 3 and 5** — pairwise non-adjacent, since 2
and 4 separate them.

    V = [0 0; 60 10; 15 10; 5 10; -5 10; -15 10; -60 10]
    E = [1 2 1; 1 3 1; 1 4 1; 1 5 1; 1 6 1; 1 7 1;
         2 3 1; 3 4 1; 4 5 1; 5 6 1; 6 7 1]
    F = [1 0; 2 1; 3 2; 4 3; 5 4; 0 5;
         1 0; 2 0; 3 0; 4 0; 5 0]

    f-rows. CCA2 stores the RAW HESSIAN, not the monomial coefficients: QuaPol.matrixForm reads
    a row as  f(x) = 1/2 x'Qx + L'x + c10  with Q = [c5 c6; c6 c7] and L = [c8;c9], so a row is
    [Q11 Q12 Q22 beta1 beta2 gamma].

      piece 1   T = (0,0),(60,10),(15,10)      [  3   -1    5    0    1   0 ]
      piece 2   T = (0,0),(15,10),(5,10)       [  3   -5   17    6   -8   0 ]
      piece 3   T = (0,0),(5,10),(-5,10)       [ 15   -2   11    2   -6   0 ]
      piece 4   T = (0,0),(-5,10),(-15,10)     [  7    2   17    8   -3   0 ]
      piece 5   T = (0,0),(-15,10),(-60,10)    [ 11   11   35   10    0   0 ]

Equivalently `q_k(x) = 1/2 x'Q_k x + beta_k'x`, constant term 0 throughout (continuity at the apex
forces all five to share it):

      Q1 = [ 3  -1;  -1   5]  det  14    beta1 = ( 0,  1)
      Q2 = [ 3  -5;  -5  17]  det  26    beta2 = ( 6, -8)
      Q3 = [15  -2;  -2  11]  det 161    beta3 = ( 2, -6)
      Q4 = [ 7   2;   2  17]  det 115    beta4 = ( 8, -3)
      Q5 = [11  11;  11  35]  det 264    beta5 = (10,  0)

All five are positive definite, so every piece is a convex quadratic — but `f` is **not** convex:
the gradient jumps across all four internal edges. The construction is a C^0 spline,
`q_{k+1} = q_k + l_k m_k` with `l_k` vanishing on the k-th cevian, so continuity holds by
construction and was checked: the largest `|q_k - q_{k+1}|` over sample points on the four shared
edges is **exactly 0**.

**It is a PLQ in the sense of Rockafellar-Wets Definition 10.20**, clause by clause
(`reference/ROCKAFELLAR-98a.pdf`): *"dom f can be represented as the union of finitely many
polyhedral sets, relative to each of which f(x) is given by an expression of the form
1/2<x,Ax> + <a,x> + alpha"*.

* `dom f` is the triangle `(0,0), (60,10), (-60,10)`, the union of the **five** triangles above —
  finitely many, all polyhedral;
* on each, `f = 1/2 x'Q_k x + beta_k'x + 0` with `Q_k` symmetric, `beta_k` a vector, `alpha = 0`;
* `f` is single-valued: the expressions agree on every shared face, checked at 201 points on each
  of the four cevians, max `|q_k - q_{k+1}| = 4.5e-13` (floating point only; by construction
  `q_{k+1} - q_k = l_k m_k` with `l_k` vanishing on the shared cevian, so the difference is
  identically zero there).

The definition asks for no more than this — in particular it does not require the polyhedral sets
to form a subdivision, which this one does anyway.

**CCA2's own test agrees.** With the rows entered in the raw-Hessian convention above,
`QuaPol.isEdgeContinuous` returns `1` on all four interior edges (and `p.isConvex()` is false, as
intended — `f` is deliberately nonconvex).

### 4.2 The vertex

Pieces 1, 3 and 5 attain the sup simultaneously at

      p = (3.907778936001628211405642364 ,
           1.414595978550451447003477036)

and `p_x` is a root of the **irreducible** quartic

      278339044 u^4 - 4292927872 u^3 + 49300591700 u^2 - 92964619904 u - 198299511119

whose **resolvent cubic is irreducible** over Q and whose discriminant is not a perfect square:

      Galois group = S4, order 24 -- not a power of 2

so **`p` lies in no iterated quadratic extension of Q**, hence in no `Q(sqrt d)`. `p_y` is
independently degree 4 with the same certificate:

      278339044 u^4 + 2972370528 u^3 - 59118592492 u^2 - 808937576592 u + 1253092464913

None of the three pencil generators is degenerate — the 3x3 determinants of `{g1=g3}`, `{g1=g5}`,
`{g3=g5}` are `-73/2254`, `-29/264`, `-3/7084` — which is precisely the escape from Theorem 3 that
pairwise non-adjacency buys.

### 4.3 Why `p` really is a vertex of `f*`

All checked in `cont5verify.m`:

1. **The active set at `p` is exactly `{1,3,5}`.** With `V = 2.861058907266`:

       piece 1  q*(p) = 2.861058907266  = V   argmax ( 1.4252,  0.3680)  inside T_1
       piece 2  q*(p) = 2.756640477602  < V
       piece 3  q*(p) = 2.861058907266  = V   argmax ( 0.2225,  0.7145)  inside T_3
       piece 4  q*(p) = 2.145083735693  < V   (argmax not even in T_4)
       piece 5  q*(p) = 2.861058907266  = V   argmax (-0.8666,  0.3128)  inside T_5

   Each `Q_k` is positive definite, so `sup over T_k <= q_k*(p)` whatever `T_k` is; pieces 2 and 4
   are therefore strictly below `V` no matter how the far edge is placed.
2. **An independent oracle agrees exactly.** `f*(s)` recomputed by exact constrained maximisation
   over each of the five triangles — unconstrained optimum when feasible, else per-edge 1-D optima
   and the vertices, using none of the CCA2 pipeline — matches `max_k g_k(s)` to **0.000e+00** at
   12 points on a circle of radius `1e-4` around `p`, and pieces 2 and 4 never win.
3. **All three regions are full-dimensional.** The three argmax points span a triangle of area
   `0.4303 != 0`, so `p` is a genuine vertex. Their angular widths at `p` are
   **162.55 / 36.32 / 161.13 degrees** — no thin sliver.

### 4.2b  `p` written in radicals

`S4` is **solvable**, so `p` IS expressible in radicals — every quartic is, by Ferrari. What `S4`
forbids is a tower of SQUARE roots; the radical expression necessarily contains a **cube root**,
and that is precisely the step no nesting of square roots can reproduce.

**The minimal polynomials**, both irreducible over Q:

      p_1 :  278339044 u^4 - 4292927872 u^3 + 49300591700 u^2 - 92964619904 u - 198299511119
      p_2 :  278339044 u^4 + 2972370528 u^3 - 59118592492 u^2 - 808937576592 u + 1253092464913

and `p_2` is a rational function of `p_1`, so the whole point lies in ONE quartic field
`Q(p) = Q(p_1)`:

      p_2 = (13482 p_1^2 - 182718 p_1 + 620099) / (110 (38 p_1 + 571))            [28 digits]

**The radical form.** With `u = t + 268307992/69584761` the quartic depresses to
`t^4 + P t^2 + Q t + R`,

      P = 425707401224338541 / 4842038963427121
      Q = 193167905847776214801877200 / 336932124022763955703081
      R = -2814575554967315391400224471082919 / 93781365293385553668053913394564

Ferrari's resolvent `8w^3 + 8P w^2 + (2P^2 - 8R) w - Q^2 = 0` clears to the integer cubic

      113523256198491191891880187027339833985032892561 w^3
    +   9980855325579385301469850064861835302110211272781 w^2
    + 222783867488801978187818427546464643933788563320900 w
    -  4664229981201917045367715995791328757606080484980000  =  0

which is irreducible over Q and has **negative discriminant** — one real root, two complex. So this
is *not* casus irreducibilis and REAL radicals suffice. Cardano gives, with

      sigma = -425707401224338541 / 14526116890281363
      A     =  219506412795617619231936130909 / 9097167348614626803983187
      B     =  195362540618591648880971999102170823763342365178607257500
              / 340569768595473575675640561082019501955098677683

      w_0 = sigma + cbrt(A + sqrt(B)) + cbrt(A - sqrt(B)) = 12.6859208594788251902336146784

and then, with `rho = sqrt(2 w_0) = 5.03704692443475757426537478177`,

      p_1 = 268307992/69584761 + ( -rho + sqrt( -2 w_0 - 2P + 2Q/rho ) ) / 2
          = 3.90777893600162821140564236354

Substituted back into the quartic the residual is `1.3e-29`. The other sign choices give the
second real root `-1.23313760731441271458` and the complex pair `6.3743675887 +- 10.3543330671 i`.
`p_2` follows from the rational relation above.

The **cube root is unavoidable**: an expression in square roots alone would place `p_1` in a
2-tower, and the Galois closure has order 24.

---

## 5. The biconjugate falls with it

At `p` three pieces attain the sup, so

      d f*(p) = conv{ x*_1(p), x*_3(p), x*_5(p) }

and `f** = conv f` is **affine exactly on that triangle**, with gradient `p`. Those three points
are therefore vertices of `f**`'s own subdivision.

*That the triangle is exactly the affine cell is a proof, not an observation.* Let
`A(x) = <p,x> - f*(p)`. By the definition of `f*`, `A <= f` everywhere, so `A <= f**`. At each
`x*_k(p)` the sup defining `f*(p)` is attained, so `A(x*_k) = f(x*_k) >= f**(x*_k) >= A(x*_k)`,
giving equality at all three corners. `f**` is convex and equals the affine `A` at the three
corners, so `f** <= A` on the triangle by convexity — hence `f** = A` there. And the region on
which `f**` has gradient `p` is exactly `d f*(p) = conv{x*_1,x*_3,x*_5}`, so the cell is maximal
and the three points are genuinely its corners. Checked numerically as well
(`envCheck.m`): `max |f** - A|` over interior samples is `2.5e-13`, while just outside each edge
`f** - A` is strictly positive and grows with the offset (`+2.5e-05` to `+3.8e-03` at offsets
`0.02` to `0.10`), with `f**` computed as `sup_s <x,s> - f*(s)` from the exact `f*` oracle.

All six coordinates were computed exactly by elimination and every one is degree 4 with an
irreducible resolvent cubic — S4 again:

      x*_1 = ( 1.4252493327541852 ,  0.36796906226092733 )
      x*_3 = ( 0.22245192703800505,  0.71449998478422378 )
      x*_5 = (-0.86662232198484083,  0.3127840434395343  )

      x*_1(1)   54554452624 u^4 - 243304777184 u^3 + 1301838417096 u^2 - 2074588008584 u + 791631204841
      x*_1(2)   54554452624 u^4 + 111499531360 u^3 -  316769385352 u^2 - 1590088667272 u + 621438933889
      x*_3(1)   7214826359524 u^4 - 4853168507936 u^3 + 5459472239060 u^2 - 1239400609552 u + 41302210321
      x*_3(2)   7214826359524 u^4 - 9619577562560 u^3 - 8947336752844 u^2 - 2536524010720 u + 8008553034241
      x*_5(1)   19399118010624 u^4 + 54575693111808 u^3 + 51689921408928 u^2 + 71148714376608 u + 47417344633633
      x*_5(2)   160323289344 u^4 - 92838320640 u^3 - 97223324256 u^2 - 3202210848 u + 11819738137

So the answer for `biconj` is the same **No**, on the same continuous input, and it is reached
without conjugating twice — which matters because `biconjCPLQ` deliberately does not do that.

### 5.1  The explicit formula for `f**`, including on `R`

`f = min_k (q_k + iota_{T_k})`, and each `epi(q_k + iota_{T_k})` is **convex**, so `f**` is the
convex hull of five convex sets and

      f**(x) = min { sum_k lam_k q_k(z_k) :  sum_k lam_k z_k = x,  z_k in T_k,  lam in the simplex }

which is explicit and valid on all of `dom f** = dom f`. The cells are indexed by the support of the
optimal `lam`.

**|supp lam| = 1.** `f** = q_k`; these are the `Omega_k`, with `q_k` the original rational
quadratics.

**|supp lam| = 2**, pieces `i` and `j`. At the optimum the two pieces share ONE tangent plane:
`grad q_i(z_i) = grad q_j(z_j) = s` with equal intercepts, which is exactly `g_i(s) = g_j(s)` —
`s` on that arc of `f*`. Then

      f**(x) = <s,x> - f*(s),     x = lam z_i(s) + (1-lam) z_j(s),   z_k(s) = Q_k^{-1}(s - beta_k)

so `f**` is AFFINE along each ruling.

**|supp lam| = 3.** `f**` is affine: `f**(x) = <s_0,x> - f*(s_0)` with `s_0` the corresponding
VERTEX of `f*`'s subdivision. Over `p` this is the S4 triangle of section 5.

**The degree on a two-piece cell — it is NOT rational.** Eliminating `s` from the two-piece system
leaves a polynomial in `lam` of **degree 4**, whose coefficients are polynomials of degree at most
2 in `x` (computed for pieces 1 and 5; the leading coefficient is the constant `-144055296`, and
the `lam^3` coefficient is constant too). So on a two-piece cell `f**` is an **algebraic function
of degree at most 4** in `x` — not a rational function, and in particular *not* quadratic over
linear. That is a multi-piece phenomenon: for a SINGLE piece one of the two faces is affine, the
elimination is linear in `lam`, and the result is the quadratic-over-linear that `RatPol` stores
(the Corollary to Theorem 4).

**Checked.** On three genuine edges of `f*` — pieces 5|3, 1|2 and 3|2 — taking `s_0` onto the arc
by Newton, forming `x` on the ruling and evaluating `<s,x> - f*(s)` reproduces a direct
`sup_s <x,s> - f*(s)` to between `4.3e-07` and `2.2e-06` (the residual is the reference sup's own
grid), and every such `x` is strictly below every `q_k` that contains it — so those points really
are in `R`, not in any `Omega_k`.

---

## 6. The three-piece discontinuous example (data structure only)

Kept because it is the minimal illustration of the mechanism and because `QuaPol` accepts it, but
it is **not a PLQ function** — it is discontinuous across its three internal edges, so it answers
the question only for the data structure, not for `conj` of a PLQ.

    V = [4 0; 24 10; -16 10; 4 -20]
    E = [1 2 1; 1 3 1; 1 4 1; 2 3 1; 3 4 1; 4 2 1]
    F = [2 1; 3 2; 1 3; 2 0; 3 0; 1 0]

      piece 1   T = (4,0),(4,-20),(24,10)      [ 1/3  -1/3   4/3   -2     4     8  ]
      piece 2   T = (4,0),(24,10),(-16,10)     [ 1/4   0     1     -1/2  -2     5/2]
      piece 3   T = (4,0),(-16,10),(4,-20)     [ 2     1     1      1     1    -5/2]

(raw-Hessian convention, as in 4.1)

Triple point `p = (0.608050881512364091, 0.358525944978487752)`, `p_x` a root of the irreducible
`3t^4 - 24t^3 + 10t^2 + 160t - 96`, resolvent cubic `t^3 - (10/3)t^2 - (896/3)t - 11008/9`
irreducible, discriminant `98559664128` not a square
(`313942^2 = 98559579364 < 98559664128 < 98560207249 = 313943^2`) — Galois group S4. Verified the
same way: argmaxes strictly interior, oracle agreement `0.000e+00`, active-piece counts `[5 3 4]`
on a circle of radius `1e-3`.

---

## 7. What this means for the number type

### 7.1 No tower can work, but the gap is small and nameable

The refutation is not "we need a bigger `d`" and not "we need `Q(sqrt2, sqrt3)`". A degree-4 S4
number is not reachable by **any** number of nested square roots, so growing `exactQ` into a tower
type buys nothing. But Theorems 2 and 3 bound the damage precisely:

* degree is never worse than 4;
* the *only* failing configuration is a vertex of `f*` where three **pairwise non-adjacent** pieces
  of the input attain the sup simultaneously;
* every vertex involving an adjacent pair — which includes every vertex of `f*` arising from a
  one-piece or two-piece input, and every vertex of a three-piece *chain* — is a 2-group vertex
  and **is** representable in a quadratic tower;
* a single convex piece is entirely rational (7.3).

That makes `exactQ` + a documented `notImplemented` on the non-adjacent-triple case a defensible
`SUPPORT_MATRIX.md` GAP rather than a silent wrong answer — provided the case is **detected**, and
it is detectable: factor the resolvent cubic of the pencil and refuse when it has no rational root.

### 7.2 The surd is only in `V`, and `V` is the one thing you could stop storing

Theorem 1 says exactly where the irrationality is and is not: face functions rational, edge conics
rational, `E`/`F`/`P` combinatorial, and only `V` of degree up to 4. A `QuaPar` in **H-form** —
faces as sign conditions on rational conics, vertices named by the *pair of conics they solve*
rather than by coordinates — is entirely rational for any rational-data input, with no extension
field at all.

The honest caveat is not small: you still have to decide **predicates** about those points — is
this intersection inside that face, which crossing comes first along an arc, do two faces share a
facet. Those are sign tests on degree-4 algebraic numbers: decidable exactly (resultants, Sturm
sequences, real-algebraic-number arithmetic) and never requiring a coordinate to be written down,
but not `exactQ` comparisons and not free. The design choice this result forces is between

    (a) exact real-algebraic-number arithmetic, degree <= 4    -- correct, heavier than exactQ
    (b) H-form storage + exact predicates only                 -- correct, biggest rewrite
    (c) exactQ + a detected, documented refusal (7.1)          -- a real, nameable GAP

### 7.3 A single convex piece is completely safe

If `q` is convex on a rational polygon `T`, the objective is concave, so the sup is at the
unconstrained optimum when feasible and on the boundary otherwise, and the regions are

* the interior branch on `{s : Q^{-1}(s-beta) in T}` — an affine preimage of a rational polygon,
  hence a **rational polygon**;
* edge and vertex branches, separated by the rational lines `x*(s) in edge` and `t* in [0,1]`.

The subdivision is **polyhedral and rational**, with no parabola and no surd anywhere. Parabolas —
and the whole question — begin with indefinite pieces and with the cross-piece max.

### 7.4 Whether ONE indefinite piece already suffices is OPEN

Both counterexamples live in the cross-piece max (Step 3). The sharper question — can a **single
indefinite** piece do it, which would put the defect inside `conjPieceCPLQ` (Step 2) — is not
settled here. The algebra permits it: for indefinite `q` the faces of `f*` are the vertex branches
`g_v` and the edge branches `g_e = g_v + L_e^2/(2 alpha_e)`, and a triple point of two edge
branches and a vertex branch is again the intersection of two rational parabolas. What is not
automatic is feasibility. Note that when `w` is an endpoint of `e`,

      g_e - g_v = L^2/(2 alpha)        and        g_e - g_w = (L - alpha)^2/(2 alpha)

are perfect squares, so those boundaries are the **lines** `t* = 0` and `t* = 1`, not parabolas —
the single-piece analogue of Theorem 3. A parabolic boundary needs the argmax to jump from the
interior of edge `e` to a **non-adjacent** vertex, and a degree-4 triple point needs two such
jumps to meet. A bounded search found no feasible instance; that is inconclusive, not a proof.

### 7.4b The section-4 example's `f*` has an ELLIPTICAL edge, so it is not a `QuaPar`

Not part of the field question, but it surfaced while checking the closure argument, and it is a
second and independent finding about the same example. **Traced against the exact oracle, not
inferred from the conic algebra.**

> **ANSWERED — the theorem is at fault, not the example.** This subsection closed by asking
> whether the construction violates an unidentified hypothesis of the `QuaPar` structure theorem,
> or whether the theorem's cross-piece step needs a condition. It is the second.
> `doc/QuaConExample.md` §1 locates the break: [JOGO] Theorem 6's proof assumes the two functions
> compared in the maximum always include a linear one *"because the pieces compared are adjacent
> in the primal"*, and Step 3b assumes the same — but Step 3b iterates over **all** pieces of `f`,
> so once there are three or more, some compared pair is non-adjacent and both functions can be
> elliptic quadratics. The adjacent case is sound and is exactly Theorem 3 above, so the theorem's
> conclusion holds for every comparison its proof actually covers. [COAP] inherits the gap in its
> route; its conclusion is not refuted. Section 4 therefore stands as written.

For a boundary between two faces of `f*`, the conic `{g_i = g_j}` has quadratic part
`(H_i - H_j)/2`, so its discriminant is `B^2 - 4AC = -det(H_i - H_j)`. For the three faces meeting
at `p`:

| pair | adjacent? | 3x3 det | quadratic-part eigenvalues | type |
|---|---|---|---|---|
| `{g1=g3}` | no | `-73/2254` | `+0.102437`, `+0.307500` | **ellipse**, centre `(-1.366197, 6.845070)`, semi-axes `8.960823`, `5.171951` |
| `{g1=g5}` | no | `-29/264` | `+0.082554`, `+0.314633` | **ellipse**, centre `(-11.083333, 10.916667)`, semi-axes `20.244388`, `10.369803` |
| `{g3=g5}` | no | `-3/7084` | `-0.085593`, `+0.072841` | **hyperbola** |
| `{g1=g2}`, `{g2=g3}`, `{g3=g4}`, `{g4=g5}` | yes | `0` (Theorem 3) | — | degenerate: real **pairs of lines**, legal `QuaPar` edges |

The adjacent pairs behave exactly as Theorem 3 predicts and are `QuaPar`-legal. The non-adjacent
ones are non-degenerate and **none is a parabola or a line**.

**Traced.** Parametrising the ellipse `{g1 = g3}` so that `theta = theta_0` is `p`, and computing
at each point which pieces attain `f*` with the exact per-piece oracle:

    dtheta   s                          pieces attaining f*
    -0.300   (4.249226,  3.659893)      2
    -0.200   (4.190551,  2.869235)      2
    -0.120   (4.103589,  2.264943)      [1 3]   <- genuine 1|3 edge
    -0.060   (4.015368,  1.830747)      [1 3]   <- genuine 1|3 edge
    -0.030   (3.963972,  1.620320)      [1 3]   <- genuine 1|3 edge
    -0.010   (3.927040,  1.482637)      [1 3]   <- genuine 1|3 edge
     0.000   (3.907779,  1.414596)      [1 3 5] <- the triple point p
    +0.010   (3.887991,  1.347098)      5

every traced point lying on the conic to `<= 2e-15`. So there is an arc of **positive length**
(`dtheta` from about `-0.12` to `0`, on an ellipse of semi-axes `~9` and `~5`) along which exactly
pieces 1 and 3 attain the supremum — a genuine 1-cell of `f*`'s subdivision, and it is an
**elliptical arc**.

`QuaPar.m` requires every edge conic to satisfy `b^2 - 4ac = 0` (a parabola, or a line as the
all-zero degenerate case). An elliptical arc is neither, so **this `f*` cannot be stored as a
`QuaPar`**, and it is also the case `clipArcByConic.m`'s header rules out ("no ellipse or
hyperbola can arise as a dual boundary").

**How much weight to put on this.** The trace uses only the definition of the conjugate — exact
constrained maximisation over each of the five triangles — so it is independent of the `QuaPar`
theory and of the CCA2 pipeline. The input is a bona fide continuous PLQ: a triangle cut by four
cevians, five pieces, all `Q_k` positive definite, continuity verified exactly (section 4.1). If
the `QuaPar` structure theorem is nevertheless correct as stated, then the error is in the input
or in my reading of what a PLQ is, not in the trace — so **the input above is the thing to check
first**, not the arithmetic.

**The field results do not depend on this.** Theorems 1, 2 and 4 use only that `f*` is a max of
rational quadratics; they never assume the subdivision is parabolic. Sections 4 and 5 stand either
way.

### 7.5 The minimum piece count — THREE for the shape result, open for S4

Five pieces is what the search produced. This subsection guessed four and reasoned that three
*pairwise* non-adjacent pieces are needed; **that reasoning was too strong**, because it conflated
the two results.

* **Shape result** (an edge of `f*` that is not a parabola or a line, 7.4b): needs one
  non-degenerate edge, i.e. just **two** non-adjacent pieces — so **three pieces suffice**.
  `doc/QuaConExample.md` §2 exhibits it: delete pieces 4 and 5 of 4.1, leaving the triangle
  `(0,0), (60,10), (-5,10)` cut by the two cevians to `(15,10)` and `(5,10)`. All three Hessians
  positive definite, pieces 1 and 3 separated by piece 2, and `{g_1 = g_3}` is the ellipse

      93 s1^2 + 38 s1 s2 - 6 s1 + 39 s2^2 - 482 s2 - 1003 = 0,   disc -13064,  3x3 det -8650208

  with pieces 1 and 3 tying at `s* = (-17/62 + sqrt(612030)/186, 3/2)` and piece 2 strictly below.
  Every number there is rational or the square root of one — the degree-4 quartic plays no part.
* **S4 field result** (section 4): does need a triple point, hence three pairwise non-adjacent
  pieces, hence at least four in the subdivision. Whether four is achievable — a central triangle
  with three glued to its edges — is **still not attempted**.

### 7.6 Iterating

`conj` maps rational data to degree-`<=4` data; a second `conj` applied to that gives degree at
most 16 over Q. One more reason `biconjCPLQ`'s refusal to compute `f**` as `conj(conj(f))` is the
right call, beyond the unbounded-domain reason it documents.

---

## 8. Data structure options

### 8.0 The measurement that drives the choice

`exactQ` was built for the A.5 cevian foot `5/2 - sqrt5/2`, which is a **vertex**. Theorem 1 says
face coefficients and edge conics are rational, so the quadratic extension is not needed for them
at all — **the vertex layer is `exactQ`'s entire reason to exist**. And the vertex layer needs
more than it offers:

> Of the twelve feasible **continuous** three-piece configurations (all in the "safe" adjacent case
> of Theorem 3, all with a 2-group Galois group), the vertex quartic factors over Q as
>
>     degree 4 (irreducible) : 10 of 12
>     2 + 2                  :  2 of 12
>     with a rational root   :  0 of 12
>
> Every element of `exactQ` has degree at most 2 over Q. **So `exactQ` cannot hold the generic
> vertex — not in the exotic S4 case, but in the ordinary one.** It is the right type one level
> short. (The 2-group structure does mean each of those ten lies in a two-level tower, i.e. a
> nested radical `a + b*sqrt(c + d*sqrt(e))`, not a single `a + b*sqrt(d)`.)

The operations a data structure has to support are worth naming, because they, not storage, decide
the cost:

    (i)   build a vertex        intersect two rational conics
    (ii)  sign of a rational polynomial at a vertex   (inside/outside a face; which piece wins)
    (iii) equality and ordering of vertices           (merge shared facets; order crossings on an arc)
    (iv)  evaluate f* at a point
    (v)   arithmetic on face and conic coefficients   -- RATIONAL only, by Theorem 1

Only (i)-(iii) ever leave Q, and none of them needs a coordinate to be *written down*.

### 8.1 The options

**1. Symbolic (`sym`) — the status quo.**
*Pro:* no gaps at all; handles degree 4, S4, and any future operator without special-casing;
already implemented.
*Con:* the "symbolic always works" premise is false in practice — a comparison can return
**unknown**, not merely slow. This session hit it: `symbolic:sym:isAlways:TruthUnknown` ("Unable to
prove `x^2 == x*y`", "`x^2 == y^2`", "`y^2 == 1`") at `plq_1p.conjugateFunction:571/574`, followed
immediately by the crash at `:579`. No canonical form, so equal quantities can compare unequal and
`merge` misses shared facets — the recorded mechanism by which the cell count grows. Expression
swell; `solve`/`simplify` cannot be interrupted from M-code (why each slow suite gets its own
process). Toolbox licence plus VPN. And it blocks the whole sym-free port (T1-T5).

**2. `double` + tolerance.**
*Pro:* fastest, no dependency.
*Con:* already refuted and recorded — one ULP made a shared facet invisible to `merge` and Step 3's
cell count grew without bound. Degree-4 vertices make near-coincidences the norm, not the
exception. Listed only for completeness.

**3. `exactQ`, a single `Q(sqrt d)` — where T2 is now.**
*Pro:* exact, int64, sign with no floating point, and the "mixing raises" rule is a genuinely good
discipline. Sufficient for everything Theorem 1 makes rational, and for degree-2 vertices such as
the A.5 foot.
*Con:* **insufficient for the generic vertex** (8.0). Additionally, two vertices carrying different
`d` cannot be compared at all under the mixing rule, and a single `f*` will routinely contain
several.
*Verdict:* right idea, one level short. Keep it for what is rational; do not make it the vertex
type.

**4. Iterated quadratic tower — nested radicals, `exactQ` over `exactQ`.**
*Pro:* by Theorem 3 this covers **every** vertex except the pairwise-non-adjacent triple, which is
to say every vertex arising from a one-, two-, or chained-three-piece input. Still exact, still
algebraic, sign tests still cheap.
*Con:* has a **proved, reachable gap** (section 4). Comparing two towers needs a common tower, so
the "mixing raises" discipline gets much harder to enforce, and the number of distinct towers in
one output is not small. The work is close to option 5, which has no gap.
*Verdict:* a defensible interim, a poor destination.

**5. General real algebraic numbers of degree <= 4** (each vertex stored as a rational quartic plus
an isolating interval; comparisons by resultants / Sturm sequences).
*Pro:* **complete for `conj` and for the envelope** — Theorems 2 and 4 both cap the degree at 4, so
this is not an open-ended algebraic-number package but a fixed, small one. Uniform: no per-vertex
field bookkeeping, no "mixing raises", no compositum. Exactly the classical solution for conic
arrangements.
*Con:* heavier per operation than `exactQ`; needs an isolating-interval + refinement layer; the
coefficients want bignum rationals rather than int64 — the counterexample's quartics already carry
14-digit coefficients and they will grow. The degree-4 cap is proved for bounded domains; the
unbounded/recession case is argued but not written out (Theorem 4's caveat).
*Verdict:* the smallest complete answer if vertices continue to be stored.

### The difference between 6 and 7

They are not competing designs — **6 is a storage decision and 7 is 6 plus an evaluation
strategy.** 7 cannot be done without 6, and 6 without 7 is correct but slow. Stated separately
because they are separate pieces of work and 6 delivers the correctness on its own.

**6. H-representation — store the rational conics, never the vertex coordinates.**

What changes: a face is stored as a list of sign conditions on rational conics; an edge is a
rational conic plus its two bounding vertices *named by which conics they solve*; a vertex is not a
pair of numbers at all but the label "the intersection of conic `i` and conic `j` that lies in
this branch". Nothing in the mesh is ever an irrational number.

*Pro:* by Theorem 1 **every stored number is rational** — no extension field appears in the data
structure at all. Canonicalisation becomes easy (normalise a rational conic vector to a primitive
integer vector), which is exactly the operation `merge` needs and the one symbolic form cannot
provide. And the degree does not grow along a chain of operators: if two objects have rational
face functions and rational edge conics, so do their sum, their max, and the overlay of their
subdivisions, so their vertices are again intersections of *rational* conics — degree stays `<= 4`
over Q however many operators are composed, instead of compounding.
*Con:* it relocates the algebra rather than removing it. Every predicate in (ii)-(iii) — is this
intersection inside that face, which of two crossings comes first along an arc, do these two faces
share a facet — is still a sign test at a degree-4 algebraic point, and under 6 alone **every one
of them** pays the full exact cost (a resultant or a Sturm sequence). Largest rewrite: V-form is
baked into `RatPar`'s `V/E/F/P`, `eval`, `createP`, `orderEdges`, plotting, and every test that
names a vertex. Plotting and user-facing output still need coordinates, so a numeric "realise this
vertex" path is required anyway — but it is then an *output* convenience, not the stored truth.
*Verdict:* the right long-term shape; correct but with no performance story of its own.

**7. The same thing, with *filtered* predicates.**

What 7 adds is only this: before running the exact kernel, evaluate the predicate in interval (or
double-with-error-bound) arithmetic. If the resulting interval excludes zero, the sign is already
*certain* and the answer is returned at double speed. Only when the interval straddles zero — a
genuine near-degeneracy — does the exact degree-4 kernel run.

*Pro:* this is the standard, proven arrangement in exact geometric computing, and the filter
succeeds on the overwhelming majority of calls, so the exact path becomes a cold path. Stored data
stay rational and canonical, exactly as in 6; the answers are identical to 6's, never approximate.
*Con:* two paths that must agree, and the filter has to be **certified** — a real error bound on
the interval evaluation, not a hand-tuned tolerance — plus either a separation bound or a
certified refinement loop so the fallback terminates. Get that wrong and you have option 2 wearing
a disguise, which is the refuted design.
*Verdict:* the recommended target. 6 buys correctness and canonical form; 7 buys back the speed.

**8. `exactQ` (or option 4) plus a *detected* refusal.**
*Pro:* smallest step from where T2 is; keeps int64 speed; and the failure is **detectable** rather
than silent — factor the vertex quartic (or the pencil's degeneracy cubic) and refuse when it does
not split into rational or quadratic factors. That turns the whole issue into one nameable
`SUPPORT_MATRIX.md` GAP.
*Con:* with plain `exactQ` the refusal rate is high (8.0: ten of twelve), so this only makes sense
paired with option 4, and then it still refuses the section-4 configuration. It is a gap, not a
solution.
*Verdict:* legitimate interim, and much better than being silently wrong — but only on top of 4.

### What option 9 actually is

**9. A per-problem number field `Q(theta)` with one primitive element.**

The idea. Instead of letting each irrational number carry its own minimal polynomial (option 5),
choose **one** algebraic number `theta` such that *every* number appearing in this particular
output lies in `Q(theta)`. Fix `theta`'s minimal polynomial `M(t)`, of degree `n`. Then every
value in the computation is written as

      c_0 + c_1 theta + ... + c_{n-1} theta^{n-1} ,   c_i rational

— an `n`-vector of rationals. Addition is componentwise; multiplication is polynomial
multiplication followed by reduction modulo `M`; inversion is the extended Euclidean algorithm
modulo `M`; and the sign of an element is decided by evaluating at the isolating interval of the
chosen real root of `M`. This is what Pari/GP, Sage and Magma call *number field arithmetic*, and
it is the textbook way to compute exactly with algebraic numbers.

Why it looks attractive here. It is uniform and closed: one type, one representation, `+ - * /`
all total, and none of `exactQ`'s "mixing raises" discipline, because by construction there is
nothing to mix — every number already lives in the same field.

Why it fails for this problem. The premise is the problem: the vertices of one `f*` do **not** all
lie in a common small field. Each vertex generates its own degree-4 field, and those fields are
generically unrelated. To hold them all you must take the **compositum**, whose degree multiplies:
two unrelated quartic fields give up to 16, three up to 64, `k` of them up to `4^k`. Since
arithmetic cost grows at least quadratically in `n`, and since the cell count in Step 3 is already
the measured bottleneck and grows with the input, the degree of `theta` grows with exactly the
quantity you most need to scale. Computing a primitive element for the compositum is itself
expensive, and it must be redone whenever a new vertex appears — that is, continuously, during the
sweep.

Option 5 avoids all of this by never combining fields at all: each vertex keeps its own quartic,
and a comparison between two numbers from different fields is done by a resultant, which is a
bounded local computation on two quartics rather than a global change of representation. So option
9 is strictly worse than 5 here, and it is worth naming only so it is not rediscovered as a good
idea later.

*Verdict:* not recommended — the right tool when all your numbers genuinely share one field, which
is exactly what this problem does not give you.

### 8.1b Composition: what grows, and what does not

The natural worry about any exact type is that using CCA2 *as a toolbox* — chaining operators —
makes the representation grow without bound. That worry is justified, but it splits into two axes
with completely different answers, and the one that actually bites is **not** the degree-4
vertices.

**Bit size: grows, but slowly, and it is the manageable axis.** Along a chain that stays inside the
class (one `conj` from a `QuaPol`, then `add` / `max` / scalar multiply, each step introducing new
small rational data), measured over six rounds with an independent rational perturbation of each
face (`growthB.m`):

    round            0     1     2     3     4     5     6
    face coeffs      4     4     6     6     5     5     5     digits
    conic coeffs     5     6     6     6     6     6     6     digits
    vertex quartic  12    18    19    18    19    20    20     digits
    vertex degree    4     4     4     4     4     4     4

Growth is **additive**, not multiplicative — denominators accumulate but reduce by `gcd`, and the
vertex quartic is a fixed-degree polynomial in the conic coefficients, so its size stays a small
multiple of theirs. Bignum rationals absorb this. Symbolic is strictly worse on the same axis.

**Degree: bounded inside the class, and unbounded the moment you leave it.** The degree-4 cap
rests on Theorem 1 — face functions and edge conics rational. That survives `add`, `max`, scalar
multiplication and subdivision overlay, because each of those produces face functions that are
rational combinations of rational ones and boundary curves that are differences of them. It does
**not** survive conjugating an object that has a genuine **parabolic** edge.

Here is why, and it is elementary. Let `h` be convex, equal to `g1` on one side of the parabola
`{y = x^2}` and to `g2` on the other. Continuity forces `g2 - g1` to vanish on the whole parabola,
so `g2 = g1 + c (y - x^2)` with `c` rational, and hence `grad g2 - grad g1 = c(-2 s1, 1)`. On the
arc the subgradient is the segment between the two gradients, so an argmax of `<x,s> - h(s)`
sitting on the arc, at `s = (t, t^2)`, satisfies `x = grad g1(t,t^2) + lambda c (-2t, 1)`.
Eliminating `lambda` gives, with `g1 = 1/2(a1 s1^2 + 2 a2 s1 s2 + a3 s2^2) + a4 s1 + a5 s2 + a6`,

      2 a3 t^3 + 3 a2 t^2 + (a1 + 2 a5 - 2 x2) t + (a4 - x1)  =  0

— a rational **cubic**, irreducible for generic data. Instance found and checked (`closure.m`):
`a = (-1,-5,1,-1,-1)`, `c = -2`, `x = (-3,1)` gives `2t^3 - 15t^2 - 5t + 2`, irreducible over Q,
with the relevant root `t* = -0.54085397978781` carrying `lambda = 0.4984` strictly inside `(0,1)`
— a genuine interior subgradient weight, residuals `8.9e-16` and `0`. If `a3 = 0` the cubic
degenerates to a quadratic, so the escape is generic, not special.

**CRUCIALLY, this does NOT happen for `f*` itself.** The coefficient that drives the cubic is
`a3`, the `s2^2` coefficient of the adjacent face *in the frame where the edge is* `y = x^2`. For a
genuine conjugate of a `QuaPol` it is **always zero**, so the cubic degenerates to a quadratic (or
lower) and the degree never escapes. The reason is structural: a parabolic edge of `f*` separates
an **edge branch** `g_e = g_v + L^2/(2 alpha)` (rank-one quadratic part) from a **vertex branch**
`g_v` (affine). The affine one has `a3 = 0` outright; and continuity forces
`Hess g_e - Hess g_v = c * Hess(y - x^2) = c * diag(-2,0)`, so `Hess g_e = diag(2c, 0)` — rank one
with its kernel along the parabola's axis, giving `a3 = 0` again. Either way the cubic's leading
term vanishes. The `a3 = 0` control row in `closure.m` is exactly this case, and it lands on a
quadratic.

So the escape needs an object whose faces are *general* quadratics across a parabolic edge, which
a `QuaPol` conjugate never produces. The instance above is a legitimate convex piecewise quadratic
with a parabolic seam, but it is **not** in the image of `conj` on `QuaPol`s.

**Where this leaves the derived operators — stated as a question, not a claim.** `moreau` calls
`conj` once, on `mu*f + 1/2||.||^2`, still polyhedral when `f` is — safe. `infConv.m:26-28` is
`conj(toQuaPar(conj f).add(toQuaPar(conj g)))`, and adding `g*`'s face `B` to both faces of `f*`
across one of `f*`'s parabolic edges leaves the edge unchanged but changes `a3` from `0` to `B`'s
own `s2^2` coefficient — so the structural reason for `a3 = 0` is destroyed. Whether that actually
produces a degree-3 point depends on whether such a configuration is reachable, i.e. on whether
`f box g` is always PLQ for convex PLQ `f, g` — since `f* + g* = (f box g)*`, it is exactly that.
**Not verified either way here.** An earlier draft of this section asserted that `infConv`,
`proxAverage` and `lasryLions` leave the class; that was inferred from the call signature alone
and is withdrawn.

**What this means for the choice.** The class *is* closed on the path that matters: `QuaPol ->
conj -> QuaPar`, and `conj(conv f) = (conv f)* = f*` returns the same `QuaPar`, so the `RatPol`
envelope conjugates straight back into the class. Within `QuaPol -> conj`, plus `add`, `max` and
scalar multiplication, faces and conics stay rational and vertices stay degree `<= 4`, so the
degree-`<=4` kernel of option 5 is complete for it permanently, no matter how long the chain.
The one open boundary is the `f* + g*` route above.

### 8.2 Summary and recommendation

| option | exact | complete for `conj` | stored data | main cost |
|---|---|---|---|---|
| 1 symbolic | yes* | yes | expressions | comparisons can return *unknown*; no canonical form; slow; blocks the port |
| 2 double + tol | no | no | doubles | refuted: lost facets, unbounded cell growth |
| 3 `exactQ` | yes | **no** (deg 4 vertices) | `a + b sqrt d` | one level short |
| 4 quadratic tower | yes | **no** (S4 case) | nested radicals | proved gap; tower bookkeeping |
| 5 real algebraic, deg <= 4 | yes | **yes** (Thms 2 and 4) | quartic + interval | needs an algebraic-number layer, bignums |
| 6 H-form | yes | yes | **all rational** | biggest rewrite; predicates still need 5 |
| 7 H-form + filtered predicates | yes | yes | **all rational** | two paths; termination bound |
| 8 `exactQ`/tower + detection | yes | no, but **honest** | as 3/4 | refuses a reachable case |
| 9 single number field | yes | yes | `Q(theta)` | compositum blow-up |

*\* "yes" for symbolic is qualified: it is exact in principle, but an undecided `isAlways` is a
wrong answer or a crash in practice, which is what was observed.*

**Recommendation.** Fix the closed sub-API of 8.1b first — one `conj` from a `QuaPol`, then
`add`/`max`/scale — because that is what makes any bounded number type a *permanent* answer rather
than a temporary one. Within it, target **7**, reach it via **5**. Concretely:

1. Keep `exactQ` — but re-scope it. Theorem 1 means the face and conic layers need only exact
   **rationals**, so `Rat`/int64 covers them; `exactQ`'s quadratic extension is not load-bearing
   there. Do not extend it into a tower as an end state (option 4's gap is proved).
2. Add a **degree-<=4 real algebraic number** kernel for the vertex layer only: a rational quartic
   with an isolating interval, sign-of-a-rational-polynomial-at-it, and comparison. That single
   component makes `conj` **and the envelope** exact and complete, and it is bounded work because
   Theorems 2 and 4 both cap the degree at 4.
3. Then move the mesh to **H-form** so no vertex is ever stored, and add the interval filter. At
   that point every number in a `QuaPar` is rational and the algebraic kernel is a cold path.
4. Until step 2 lands, pair whatever is in place with option **8**'s detection, so the reachable
   failure is a named refusal instead of a wrong answer.

---

## 9. Two ledgers: intrinsic vs incidental surds

**Intrinsic — the number type must handle these or refuse them by name.**

| where | what | degree / group |
|---|---|---|
| vertex of `f*`, three pairwise NON-adjacent pieces | base point of a pencil of rational conics, no degenerate member | 4, **S4** — no tower |
| vertex of `f*`, some two pieces adjacent | pencil has a rational degenerate member (Theorem 3) | `<= 4`, 2-group — tower works |
| vertex of `f**`'s affine face over such a point | the argmax points `x*_k(p)` | same as the `f*` vertex above it |
| splitting an adjacent pair's degenerate conic into its two lines | one `sqrt` of a rational discriminant | 2 |
| A.5 cevian foot `5/2 - sqrt5/2` | already in `DECISIONS.md` | 2 |

**Incidental — an artefact of the route, not of the answer.**

`convEnvCPLQ:181` builds `bilinearFrame(Q)` from `sqrt(lam1/2)` and `sqrt(-lam2/2)`, so a surd
enters the **coefficients** of the intermediate representation. Theorem 1 proves the face functions
of the *final* `f*` are rational, so any irrational coefficient surviving into `f*` is a property
of the path taken, not of the object computed. That does not make the frame wrong — it is a
legitimate change of variables — but an irrational **coefficient** and an irrational **vertex** are
different findings and must not be pooled. Only the second is forced.

(This entry is about `f*`. The intermediate `conv f` is a different object — `RatPol` pieces,
genuine rational functions — and nothing is claimed here about whether *its* coefficients must be
rational.)

---

## 10. What was actually run

MATLAB R2024b, Symbolic Math Toolbox, on the shared machine. Scripts are in the session
scratchpad, deliberately not in the repo; none imports CCA2 except the two `runCCA2*` ones.

| script | what it establishes | result |
|---|---|---|
| `galoisSearch.m` | small integer convex quadratics with irreducible quartic **and** irreducible resolvent cubic | 2 certificates in 400 trials |
| `buildCE.m`, `verifyCE.m` | the three-piece (discontinuous) example, oracle, transversality | agreement `0.000e+00`, counts `[5 3 4]` |
| `degrees.m` | exact minimal polynomials of all 8 stored coordinates, three-piece case | all degree 4, all resolvent-irreducible |
| `contCE.m`, `contExact.m` | twelve feasible CONTINUOUS three-piece configurations | 10 irreducible quartics, **all 10** resolvent-reducible |
| `adjDegen.m` | is adjacency the reason? 946 random rational adjacent pairs vs a control | relative `\|det\|` `1.4e-15` vs `9.6e-2` |
| `adjProof.m` | Theorem 3 as a **symbolic identity** in free symbols | `det == 0` identically |
| `cont5.m` | 5-piece continuous fan, active set `{1,3,5}`, all feasibility constraints | 1 candidate in 400 000 trials |
| `cont5exact.m`, `cont5verify.m`, `cont5poly.m` | the counterexample: S4 certificate, continuity, active set, oracle, region widths, exact minimal polynomials | all as reported in section 4 |
| `onePiece.m` | is one indefinite piece enough (7.4)? | no certificate found; **inconclusive** |
| `envCheck.m` | the envelope cell structure on the counterexample | `f** = A` inside to `2.5e-13`, strictly above outside |
| `growth.m`, `growthB.m` | bit-size growth along a legal operator chain (8.1b) | additive: quartic 12 -> 20 digits over 6 rounds, degree stays 4 |
| `closure.m` | is the class closed under `conj` of a parabolic-edge object (8.1b)? | **no** — irreducible rational cubic, degree 3 |
| `conicType.m` | conic type of each `f*` boundary (7.4b) | adjacent: degenerate line pairs; non-adjacent: ellipse / ellipse / hyperbola |
| `edgeType.m` | trace the `U_3\|U_6` edge against the exact oracle (7.4b) | positive-length arc of a bounded **ellipse**; `f*` is not a `QuaPar` |
| `edgeReal.m` | are the edges at `p` real? circle walk at radii 1e-3 / 1e-2 / 1e-1 | exactly 3 transitions each, same angles, argmax jumps 1.25 / 1.16 / 2.29 |
| `ratpol.m`, `degreeF2.m` | the `RatPol` structure for ONE piece, worked and verified | `f** = 2x2^2/(2 + x2 - x1)`; matches the numeric `f**` to `1.7e-16` |
| `radicals2.py` | `p` in REAL radicals, Ferrari + Cardano (4.2b) | resolvent cubic irreducible with disc < 0; residual `1.3e-29` |
| `fpp_formula.py` | `p_2` as a rational function of `p_1`; degree of `f**` on a two-piece cell (5.1) | quartic in `lam` — algebraic of degree 4, NOT rational |
| `verify_ruled2.py` | the two-piece formula on three genuine edges of `f*` (5.1) | matches a direct sup to `4.3e-07 .. 2.2e-06`; every point in `R` |
| `quacon.py` | the figure `doc/QuaCon.svg` | exact rational face data; edges classified exactly |
| `runCCA2.m`, `runCCA2cont.m` | feed both counterexamples to CCA2's own `conj` | **FAILED** — see below |

Four self-corrections, recorded because each would otherwise have produced a wrong answer
silently:

* the CCA2 `f`-rows first published here used the **monomial** convention
  `[x^2 xy y^2 x y const]` and so halved the two diagonal entries. `QuaPol.matrixForm` stores the
  **raw Hessian** (`f = 1/2 x'Qx + L'x + c10`), so the mesh handed to CCA2 encoded a different
  function, and `isEdgeContinuous` correctly reported the interior edges as discontinuous. With
  the rows corrected (4.1) `isEdgeContinuous` returns 1 on all four interior edges. **No
  mathematical result is affected**: every theorem, oracle check and trace in this note is
  computed from `Q_k`, `beta_k` directly and never reads a `QuaPol`. The two `conj` failures below
  were re-run on the corrected inputs and both reproduce unchanged.

* the original framing assumed a PLQ could be discontinuous. It cannot; the three-piece example is
  a legal data structure but not a PLQ, and the whole continuous case (sections 3-5) exists
  because of that correction.
* an LLL-based integer-relation finder written for `degrees.m` returned "no relation of degree
  <= 6" even for a number whose quartic was already known. Discarded rather than debugged, and
  replaced by exact elimination (substitute, resultant, factor over Q, pick the factor vanishing
  at the numeric value). Every degree quoted here comes from the exact route.
* the first `onePiece.m` run computed the sup by sampling 4001 points per edge. A sampled sup
  **under-reports** by `O(h^2)`, the same size as the acceptance tolerance, so genuine triple
  points could be rejected. Replaced by the exact closed-form max over vertices and feasible edge
  optima.

**CCA2's `conj` cannot compute either input today**, in two different ways.

* The **three-piece** one ran 69 s and died with `symbolic:kernel:InvalidReturnType` inside
  `symbolicFunction.getDen` <- `symbolicFunction.subsF` <- `region.funcVertices:2731`, after
  printing an `intmax` sentinel coordinate (`2147483647`) for a vertex.
* The **five-piece continuous** one is accepted as a mesh (`nv=7 ne=11 nf=5`, `PARTITION OK`
  twice) and then dies after 8 s at `plq_1p.conjugateFunction:579` with
  `MATLAB:UndefinedFunction / Unrecognized function or variable 'd'`, via
  `plq_1p.conjugate:384` <- `plq.maximum:167` <- `conjCPLQ:164`. It is immediately preceded by
  four `symbolic:sym:isAlways:TruthUnknown` warnings at `conjugateFunction:571` and `:574`
  (`Unable to prove 'x^2 == x*y'`, `'x^2 == y^2'`, `'y^2 == 1'`), so the `isAlways` tests fall
  through to a branch that never assigns `d`. That is a plain code defect on a reachable path,
  not a numeric or representational limit — but it is a **separate finding** from this
  document's subject and was not investigated further.

Neither failure is evidence for or against anything here:
the result is about the mathematical object `f*`, proved in section 2 and exhibited in section 4
against an oracle that uses none of the pipeline. It does mean neither counterexample can be
observed end-to-end inside CCA2 today, so this is a **design-level** finding, not a reproduction
of a wrong answer. `region.m` and `exactQ.m` were under active edit in another session and were
not touched or investigated.

Nothing in the repository was modified, no test bucket was run (nothing changed that could turn a
test red), and nothing was committed.

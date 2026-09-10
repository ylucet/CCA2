# Why the `I3 | I5` arc really is a hyperbola

_Written 2026-09-10. A self-contained proof for the reader of `doc/QuaConExample.md` who is not
convinced that a **non-degenerate** hyperbolic arc can occur in the conjugate of a bivariate PLQ
function, because every example she built came out as a **degenerate** one — a pair of straight
lines._

**The short version.** Her examples are not wrong and the document is not wrong. She has been
comparing pieces that share an edge, and for such a pair the conic is *forced* to degenerate into
two lines — that is Proposition 5 below, and it is the reason a two-piece example can never
produce a curved arc. The document's arc comes from two pieces that share **no** edge, and for
those nothing is forced. Whether a given pair gives a genuine hyperbola, a genuine ellipse, or two
lines is decided by **two determinants of the input data**, computable by hand before any conjugate
is calculated (Theorem 6). For the pair in question they are `-265` and `144`, and that is the
whole proof.

Everything here is exact rational arithmetic. Section 9 lists what was recomputed, and with which
independent tool.

---

## 1. Setting and notation

Let `f : R^2 -> R u {+oo}` be piecewise linear-quadratic in the sense of Rockafellar–Wets
Definition 10.20: `dom f` is a finite union of polyhedra `T_1, ..., T_N` and on each of them `f`
agrees with one quadratic expression

        q_k(x) = 1/2 <x, Q_k x> + <beta_k, x> + gamma_k ,

with `Q_k` symmetric. Throughout we assume every `Q_k` is **positive definite**, which is the case
in `QuaConExample.md` §4.1; `f` itself is *not* convex, because the gradient jumps across the
internal edges.

Since `dom f` is the union of the `T_k`,

        f*(s) = sup_{x in dom f} [ <s,x> - f(x) ]  =  max_{k}  h_k(s) ,
        h_k(s) := sup_{x in T_k} [ <s,x> - q_k(x) ] ,                              (1.1)

so `f*` is a finite max of `N` functions, each of which is the conjugate of one quadratic
*restricted to one polyhedron*. Write `g_k` for the conjugate of the same quadratic taken over all
of `R^2`:

        g_k(s) := sup_{x in R^2} [ <s,x> - q_k(x) ] .

**Lemma 1 (the interior cell).** `g_k(s) = 1/2 <s - beta_k, Q_k^{-1}(s - beta_k)> - gamma_k`, and
the (unique) maximiser is `x_k(s) = Q_k^{-1}(s - beta_k)`.

*Proof.* The objective `x |-> <s,x> - q_k(x)` is strictly concave (because `Q_k > 0`), so it has a
unique maximiser, characterised by first-order optimality `s - Q_k x - beta_k = 0`, i.e.
`x = x_k(s)`. Substituting, and writing `x = x_k(s)` for brevity,

        <s,x> - 1/2<x,Q_k x> - <beta_k,x> - gamma_k
          = <s - beta_k, x> - 1/2<x, Q_k x> - gamma_k
          = <Q_k x, x> - 1/2 <x, Q_k x> - gamma_k          (since s - beta_k = Q_k x)
          = 1/2 <x, Q_k x> - gamma_k
          = 1/2 <s - beta_k, Q_k^{-1}(s - beta_k)> - gamma_k .                     ∎

Define the **interior cell** of piece `k`,

        C_k := { s : x_k(s) in int T_k } ,                                          (1.2)

an open convex polyhedron (it is the preimage of `int T_k` under the affine bijection
`s |-> Q_k^{-1}(s - beta_k)`). On `C_k` the unconstrained maximiser is feasible for `T_k`, hence

        s in C_k   =>   h_k(s) = g_k(s) .                                           (1.3)

These are the faces `I1, ..., I5` of `QuaConExample.md` §3.2. If, on some open set, the max in
(1.1) is attained by piece `i` on one part and by piece `j` on another, and both parts meet the
respective interior cells, then the boundary between them is contained in the level set

        Z_{ij} := { s : g_i(s) = g_j(s) } ,                                          (1.4)

which is a conic, since `g_i - g_j` is a quadratic polynomial in `s`. **The question "is the arc
curved?" is the question "is `Z_{ij}` a non-degenerate conic?"** — a degenerate conic is a union of
lines, so an arc of one is a straight segment.

---

## 2. The two invariants of a conic, and what they decide

Write a quadratic polynomial in `s = (s1,s2)` as

        P(s) = A s1^2 + 2B s1 s2 + C s2^2 + 2p1 s1 + 2p2 s2 + R                     (2.1)

and attach to it the symmetric `3x3` matrix and its leading `2x2` block

              [ A  B  p1 ]                       [ A  B ]
        M  =  [ B  C  p2 ] ,             A2  =   [ B  C ] ,       p = (p1, p2) ,
              [ p1 p2 R  ]

so that `P(s) = <(s,1), M (s,1)>`. What follows are three elementary facts, proved here rather
than cited because the precise form of the second is what the rest of the note leans on.

**Proposition 2 (centre and canonical form).** If `det A2 != 0` then `P` has a unique critical
point ("the centre")

        c = -A2^{-1} p ,        and       P(c + u) = <u, A2 u> + k ,
        with       k = R - <p, A2^{-1} p> = det M / det A2 .                        (2.2)

*Proof.* Expanding, `P(s) = <s, A2 s> + 2<p, s> + R`, so with `s = c + u`,

        P(c+u) = <u, A2 u> + 2<A2 c + p, u> + [ <c, A2 c> + 2<p,c> + R ] .

The linear term vanishes exactly when `A2 c + p = 0`, i.e. `c = -A2^{-1}p`, which is the unique
solution because `A2` is invertible. The constant is then

        <A2^{-1}p, p> - 2<p, A2^{-1}p> + R = R - <p, A2^{-1}p> = k .

For the second expression for `k`, use the Schur-complement determinant formula on the block matrix
`M = [[A2, p],[p^T, R]]`:  `det M = det A2 * (R - p^T A2^{-1} p) = det A2 * k`.                ∎

The **type** of the conic is read off `det A2`. Writing (2.1) in the more familiar
`a s1^2 + b s1 s2 + c s2^2 + ...` form, `a = A`, `b = 2B`, `c = C`, so

        b^2 - 4ac = 4B^2 - 4AC = -4 det A2 ,                                        (2.3)

and `det A2 < 0` (equivalently `b^2 - 4ac > 0`, equivalently `A2` **indefinite**) is the hyperbolic
case. Diagonalising `A2 = O diag(lam1, lam2) O^T` by a rotation `O`, Proposition 2 turns `P = 0`
into

        lam1 u^2 + lam2 v^2 + k = 0 ,      lam1 * lam2 = det A2 .                   (2.4)

**Corollary 3 (in the hyperbolic case, `det M` is exactly the "two lines or not" switch).** Suppose
`det A2 < 0`, so `lam1 < 0 < lam2`. Then:

* if `k != 0` (equivalently `det M != 0`) the solution set of (2.4) is a **genuine hyperbola**: two
  disjoint unbounded branches, centre `c`, semi-axes `sqrt(|k/lam1|)` and `sqrt(|k/lam2|)`, and two
  distinct real asymptotes through `c` with directions given by the two real roots of the leading
  form `A + 2B m + C m^2 = 0`;
* if `k = 0` (equivalently `det M = 0`) then (2.4) reads `lam1 u^2 = -lam2 v^2`, i.e.
  `v = +- sqrt(-lam1/lam2) u`: exactly **two distinct real lines through `c`** — a *degenerate*
  hyperbola.

*Proof.* Immediate from (2.4) after dividing by `-k` in the first case (`lam1 < 0 < lam2` gives the
signs `u^2/a^2 - v^2/b^2 = 1` after possibly swapping the two axes) and by factoring the difference
of two squares in the second. `k = 0 <=> det M = 0` because `det M = k det A2` and
`det A2 != 0`.                                                                                  ∎

So for the pair we care about, *everything* hinges on the sign of `det A2` and the vanishing of
`det M`. It is worth recording the one direction of the factorisation criterion that we use later
in general form.

**Lemma 4 (a product of affine forms has `det M = 0`).** If `P = l * m` with
`l(s) = a1 s1 + a2 s2 + a3` and `m(s) = b1 s1 + b2 s2 + b3` affine, then, writing
`alpha = (a1,a2,a3)` and `bet = (b1,b2,b3)` in `R^3`,

        M = 1/2 ( alpha bet^T + bet alpha^T ) ,       hence     det M = 0 ,

and moreover

        det A2 = a1 b1 a2 b2 - ((a1 b2 + a2 b1)/2)^2 = -((a1 b2 - a2 b1)/2)^2 <= 0 ,  (2.5)

with strict inequality iff `l` and `m` are non-parallel. In words: **a real pair of lines is always
of hyperbolic (or, if parallel, parabolic) type — never elliptic.**

*Proof.* Homogenise: `l(s) = <alpha, (s,1)>` and `m(s) = <bet, (s,1)>`, so
`P(s) = <alpha,(s,1)> <bet,(s,1)> = <(s,1), alpha bet^T (s,1)>`; symmetrising the matrix (which
does not change the quadratic form) gives `M = (alpha bet^T + bet alpha^T)/2`. That matrix is a sum
of two rank-one matrices, hence has rank at most `2 < 3`, hence `det M = 0`. For (2.5), expand `l*m`
to read off `A = a1 b1`, `C = a2 b2`, `2B = a1 b2 + a2 b1`, and use the elementary identity
`((P+Q)/2)^2 - PQ = ((P-Q)/2)^2` with `P = a1 b2`, `Q = a2 b1`.                                 ∎

> **Aside, and it matters for the reader's confusion.** By Lemma 4 the degenerate conics she has
> been getting — two crossing lines — *are* hyperbolas. They are the degenerate members of the
> family. So the question was never "hyperbola or something else"; it was always "degenerate or
> not", i.e. `det M = 0` or `det M != 0`.

---

## 3. The conic `Z_ij` and its primal companion

Fix two pieces `i != j`. Introduce

        D    := Q_i^{-1} - Q_j^{-1} ,
        w    := Q_i^{-1} beta_i - Q_j^{-1} beta_j ,
        kap  := 1/2 ( <beta_i, Q_i^{-1} beta_i> - <beta_j, Q_j^{-1} beta_j> ) - (gamma_i - gamma_j) .

Expanding the two expressions of Lemma 1,

        g_i(s) - g_j(s) = 1/2 <s, D s> - <w, s> + kap ,

so, doubling to clear the `1/2` and comparing with (2.1),

        2 (g_i - g_j)(s) = <s, D s> - 2<w, s> + 2 kap    has matrix

               [  D     -w   ]
        M_d =  [ -w^T   2 kap ]        (the DUAL matrix; A2 = D, p = -w, R = 2 kap).   (3.1)

Do the same on the primal side. With `delta := beta_i - beta_j`,

        2 (q_i - q_j)(x) = <x, (Q_i - Q_j) x> + 2 <delta, x> + 2(gamma_i - gamma_j)

               [ Q_i - Q_j        delta        ]
        M_p =  [   delta^T   2(gamma_i-gamma_j) ]      (the PRIMAL matrix).             (3.2)

Both matrices are built from the same twelve numbers. The next two lemmas say they carry the same
information.

---

## 4. Lemma A — the type transfers

**Lemma A.** `det D * det Q_i * det Q_j = det(Q_i - Q_j)`, and consequently, since
`det Q_i det Q_j > 0`,

        sign( det D ) = sign( det(Q_i - Q_j) ) .

*Proof.* Factor the difference of inverses:

        D = Q_i^{-1} - Q_j^{-1} = Q_i^{-1} ( Q_j - Q_i ) Q_j^{-1} ,

which is checked by multiplying out: `Q_i^{-1}(Q_j - Q_i)Q_j^{-1} = Q_i^{-1} - Q_j^{-1}`. Taking
determinants, `det D = det(Q_j - Q_i) / (det Q_i det Q_j)`. Finally in dimension `2`,
`det(-N) = (-1)^2 det N = det N`, so `det(Q_j - Q_i) = det(Q_i - Q_j)`.                        ∎

Combining Lemma A with (2.3) applied to (3.1) gives a formula worth displaying, because it lets one
compute the discriminant of the arc from the two input Hessians and nothing else. For the
*unscaled* difference `g_i - g_j` (matrix `M_d / 2`, so `A2 = D/2`):

        b^2 - 4ac  =  -4 det(D/2)  =  - det D  =  - det(Q_i - Q_j) / (det Q_i det Q_j) .  (4.1)

**Interpretation.** `Q_i - Q_j` indefinite means neither piece is uniformly more curved than the
other — neither `Q_i >= Q_j` nor `Q_j >= Q_i` in the Loewner order. *That* is the geometric content
of "the arc is hyperbolic rather than elliptic", and it is a statement purely about the input.

---

## 5. Lemma B — the degeneracy transfers

**Lemma B.** `det M_d * det Q_i * det Q_j = - det M_p`. In particular
`det M_d = 0 <=> det M_p = 0`.

*Proof.* Both sides are polynomials in the twelve entries `(Q_i, Q_j, beta_i, beta_j, gamma_i,
gamma_j)` — on the left, after clearing the inverses through the factor `det Q_i det Q_j` — so it is
enough to prove the identity on the dense open set where in addition `S := Q_j - Q_i` is invertible;
it then extends to all data by continuity (a polynomial vanishing on a dense set vanishes
identically).

Abbreviate `A := Q_i`, `B := Q_j`, `a := beta_i`, `b := beta_j`, `c := gamma_i - gamma_j`,
`delta := a - b`, `S := B - A`. Apply Proposition 2's Schur formula to each matrix:

        det M_d = det D  * ( 2 kap - <w, D^{-1} w> ) ,                                  (5.1)
        det M_p = det(A-B) * ( 2c - <delta, (A-B)^{-1} delta> ) .                       (5.2)

By Lemma A, `det D = det(A-B) / (det A det B)`. So the claim `det M_d det A det B = -det M_p`
reduces, after cancelling `det(A-B) != 0`, to the scalar identity

        2 kap - <w, D^{-1} w>  =  - ( 2c - <delta, (A-B)^{-1} delta> ) ,

i.e., since `2 kap = <a, A^{-1}a> - <b, B^{-1}b> - 2c` and `(A-B)^{-1} = -S^{-1}`, to

        <a, A^{-1} a> - <b, B^{-1} b> - <w, D^{-1} w>  =  - <delta, S^{-1} delta> .     (5.3)

Now compute. From `D = A^{-1} S B^{-1}` (Lemma A's factorisation) we get
`D^{-1} = B S^{-1} A`. Put

        u := A^{-1} a ,      v := B^{-1} b ,      z := u - v ,

so that `a = A u`, `b = B v`, `w = u - v = z`, and `delta = A u - B v`. In these variables (5.3)
reads

        <u, A u> - <v, B v> - <z, B S^{-1} A z> + <A u - B v, S^{-1}(A u - B v)> = 0 .   (5.4)

Handle the third and fourth terms separately, using `A = B - S` in the third and
`A u - B v = B z - S u` in the fourth (check: `B z - S u = B u - B v - S u = (B - S)u - Bv = Au -
Bv`).

* Third term. `B S^{-1} A = B S^{-1} (B - S) = B S^{-1} B - B`, hence
  `<z, B S^{-1} A z> = <z, B S^{-1} B z> - <z, B z>`.
* Fourth term. With `S^{-1}` symmetric,

        <Bz - Su, S^{-1}(Bz - Su)> = <z, B S^{-1} B z> - 2 <Bz, S^{-1} S u> + <Su, S^{-1} S u>
                                   = <z, B S^{-1} B z> - 2 <Bz, u> + <Su, u> .

Substituting both into (5.4), the `<z, B S^{-1} B z>` terms cancel and we are left with

        <u,Au> - <v,Bv> + <z,Bz> - 2<Bz,u> + <u,Su> .

Use `A + S = B` on the first and last terms: `<u,Au> + <u,Su> = <u,Bu>`. Expand
`<z,Bz> = <u,Bu> - 2<u,Bv> + <v,Bv>` and `-2<Bz,u> = -2<u,Bu> + 2<u,Bv>`. The total is

        <u,Bu> - <v,Bv> + <u,Bu> - 2<u,Bv> + <v,Bv> - 2<u,Bu> + 2<u,Bv>  =  0 ,

which is (5.4).                                                                                 ∎

Lemmas A and B together say something clean enough to remember:

> **The Legendre transform preserves the projective type of the difference conic.** For two
> strictly convex quadratics, the primal conic `{q_i = q_j}` and the dual conic `{q_i* = q_j*}` have
> discriminants of the same sign and `3x3` determinants that vanish together. Ellipse goes to
> ellipse, hyperbola to hyperbola, line pair to line pair.

Two convenient closed forms, obtained by combining Lemma B with the scaling rules (multiplying a
conic by `t` scales `det A2` by `t^2` and `det M` by `t^3`); they are what makes the whole question
a hand computation. For the unscaled `g_i - g_j`,

        b^2 - 4ac    = - det(Q_i - Q_j) / (det Q_i det Q_j) ,                     [= (4.1)]
        det M(g_i-g_j) = - det M_p / ( 8 det Q_i det Q_j ) .                             (5.5)

---

## 6. Lemma C — why an adjacent pair can never give a curved arc

**Proposition 5 (adjacency forces degeneracy).** Suppose `f` is (finite and) continuous on
`T_i u T_j` and `T_i n T_j` contains a segment of positive length, lying on the line
`L = { x : l(x) = 0 }` for some non-zero affine `l`. Then

        q_i - q_j = l * m       for some affine  m ,

hence `det M_p = 0` (Lemma 4), hence `det M_d = 0` (Lemma B): the conic `Z_ij` is **degenerate**.
Moreover by (2.5) it is either two distinct real lines (the generic case) or a double line — never
a genuine hyperbola and never a genuine ellipse.

*Proof.* Continuity of `f` on `T_i u T_j` means `f` is single-valued there, so `q_i(x) = q_j(x)` at
every `x in T_i n T_j`, in particular on a segment `[x0, x1] c L` with `x0 != x1`. Parametrise
`L` by `t |-> x0 + t (x1 - x0)`. The composition `t |-> (q_i - q_j)(x0 + t(x1-x0))` is a univariate
polynomial of degree at most `2` vanishing for every `t in [0,1]`, hence infinitely often, hence
identically. So `q_i - q_j` vanishes on **all** of `L`, not just on the shared segment.

Now choose affine coordinates `(y1, y2)` in which `l(y) = y1`, so `L = {y1 = 0}`. In these
coordinates `P := q_i - q_j` is again a polynomial of degree `<= 2`:

        P = c20 y1^2 + c11 y1 y2 + c02 y2^2 + c10 y1 + c01 y2 + c00 .

`P` vanishes identically on `{y1 = 0}`, i.e. `c02 y2^2 + c01 y2 + c00 = 0` for all `y2`, so
`c02 = c01 = c00 = 0`. Therefore

        P = y1 ( c20 y1 + c11 y2 + c10 ) = l * m ,      m affine,

as claimed. The rest is Lemma 4 and Lemma B.                                                    ∎

Two consequences, which between them explain the reader's experience completely.

**Corollary 5a (two pieces are never enough).** Let `f` have two pieces whose intersection contains
a segment — which is the normal situation for a two-piece subdivision of a connected planar domain,
since the two polyhedra have to meet along their common boundary. Then the single pair available is
adjacent, so by Proposition 5 the only conic in sight is a pair of lines. **No such two-piece
example can exhibit a curved arc, whatever the two Hessians are.** This is exactly the family the
reader built, which is why the outcome was always the same.

**Corollary 5b (what is needed instead).** A non-degenerate arc requires a pair `(i,j)` with
`T_i n T_j` containing no segment. Note that sharing a single *point* is harmless: it imposes the
one scalar equation `q_i = q_j` at that point, not the two-parameter family of equations that
forces a factorisation. The first configuration that admits such a pair is **three pieces with the
middle one separating the other two** — which is precisely `QuaConExample.md` §2. (In §4.1's
five-piece fan, pieces `3` and `5` share only the apex `(0,0)`, where both quadratics vanish.)

---

## 7. The theorem, and the whole example in one table

Collecting Lemmas A, B and Corollary 3:

**Theorem 6.** Let `T_i, T_j` be two pieces with positive definite Hessians `Q_i, Q_j`, and let
`M_p` be the primal `3x3` matrix (3.2). Then the conic `Z_ij = {g_i = g_j}` separating the two
interior faces is

* of **hyperbolic** type iff `det(Q_i - Q_j) < 0`, of elliptic type iff `> 0`, parabolic iff `= 0`;
* **non-degenerate** iff `det M_p != 0`;

and therefore a **genuine hyperbola** iff `det(Q_i - Q_j) < 0` **and** `det M_p != 0`. Both tests
are integer determinants of the input data; no conjugate need be computed.

Applied to `QuaConExample.md` §4.1, whose data are

        Q1 = [ 3 -1; -1  5]  det  14   beta1 = ( 0,  1)   gamma = 0
        Q2 = [ 3 -5; -5 17]  det  26   beta2 = ( 6, -8)   gamma = 0
        Q3 = [15 -2; -2 11]  det 161   beta3 = ( 2, -6)   gamma = 0
        Q4 = [ 7  2;  2 17]  det 115   beta4 = ( 8, -3)   gamma = 0
        Q5 = [11 11; 11 35]  det 264   beta5 = (10,  0)   gamma = 0

the theorem gives the following for all **ten** pairs of the five pieces, where the last two columns are the closed forms (4.1) and
(5.5) and may be compared line by line with the tables in `QuaConExample.md` §3.3 and §4.2:

| pair | adjacent? | `det(Q_i-Q_j)` | `det M_p` | verdict | `b^2-4ac` | `det M` of `g_i-g_j` |
|---|---|---|---|---|---|---|
| `1 \| 2` | yes | `-16` | **`0`** | two lines | `4/91` | `0` |
| `2 \| 3` | yes | `-81` | **`0`** | two lines | `81/4186` | `0` |
| `3 \| 4` | yes | `-64` | **`0`** | two lines | `64/18515` | `0` |
| `4 \| 5` | yes | `-9` | **`0`** | two lines | `3/10120` | `0` |
| `1 \| 3` | no | `+71` | `584` | **ellipse** | `-71/2254` | `-73/2254` |
| `1 \| 4` | no | `+39` | `1024` | **ellipse** | `-39/1610` | `-64/805` |
| `1 \| 5` | no | `+96` | `3248` | **ellipse** | `-2/77` | `-29/264` |
| `2 \| 4` | no | `-49` | `-40` | **hyperbola** | `49/2990` | `1/598` |
| `2 \| 5` | no | `-112` | `-224` | **hyperbola** | `7/429` | `7/1716` |
| `3 \| 5` | no | `-265` | `144` | **HYPERBOLA** | `265/42504` | `-3/7084` |

Read the table twice.

1. **Every one of the four adjacent pairs has `det M_p = 0` exactly, and hyperbolic type.** That is
   Proposition 5 in action: adjacent comparisons give *degenerate hyperbolas*, i.e. pairs of
   straight lines, and they can give nothing else. This is what the reader has been observing, and
   it is a theorem about her construction, not evidence about ours.
2. **Every one of the six non-adjacent pairs has `det M_p != 0`** — non-degeneracy is the generic
   situation the moment adjacency stops forcing the issue — and their types are then decided by
   `det(Q_i - Q_j)`: three ellipses and three hyperbolas. For the three pairs that meet at the
   triple point `p`, the rows reproduce exactly the six numbers printed in `QuaConExample.md`
   (`-71/2254`, `-2/77`, `265/42504` for `b^2-4ac`, and `-73/2254`, `-29/264`, `-3/7084` for the
   `3x3` determinants), recomputed here from the primal data through (4.1) and (5.5) without
   touching a conjugate. That agreement is a strong consistency check on the document in both
   directions.
3. **A caveat on how much the table claims.** It classifies the conic `Z_ij` for every pair. It
   does *not* claim that every pair actually contributes an edge to `f*`: whether the two interior
   cells `C_i, C_j` meet, and whether some third piece takes the maximum first, is separate
   information, and `QuaConExample.md` §3.1 records that six of the 23 formulas never win anywhere.
   Section 9 below is what establishes that the `3|5` conic is genuinely realised as an edge of
   `f*`. The table's role is the converse direction, which is the one under dispute: for an
   adjacent pair no amount of geometry could produce a curved edge, because there is no curve to
   produce.

---

## 8. The witness, spelled out

Pieces `3` and `5`. They are separated by pieces `4` and `2` and meet only at the apex, so
Proposition 5 does not apply and nothing constrains `q_3 - q_5`.

**Type.** `Q3 - Q5 = [ 4 -13 ; -13 -24 ]`. Its indefiniteness needs no eigenvalues, just two
vectors:

        <e1, (Q3-Q5) e1> =   4 > 0 ,        <e2, (Q3-Q5) e2> = -24 < 0 ,

so `Q3 - Q5` is indefinite; concretely `det(Q3-Q5) = 4*(-24) - (-13)^2 = -96 - 169 = -265 < 0`.
Piece 3 is the stiffer one along `x1`, piece 5 along `x2` — neither dominates the other.
Hyperbolic type. By (4.1) the arc's discriminant is `b^2 - 4ac = 265/(161*264) = 265/42504 > 0`.

**Non-degeneracy.** `beta3 - beta5 = (-8,-6)` and `gamma3 - gamma5 = 0`, so

              [   4  -13  -8 ]
        M_p = [ -13  -24  -6 ] ,
              [  -8   -6   0 ]

        det M_p = 4(0-36) - (-13)(0-48) + (-8)(78-192) = -144 - 624 + 912 = 144 != 0 .

By Lemma B, `det M_d = -144/42504 != 0`. **By Theorem 6 the arc `{g_3 = g_5}` is a non-degenerate
hyperbola.** ∎

**The curve itself.** Written out from Lemma 1 and cleared of denominators by
`2 * 42504 = 2 det Q3 det Q5`, `Z_35` is

        -2731 s1^2 + 4598 s1 s2 + 2189 s2^2 + 107420 s1 + 9988 s2 - 421996 = 0 ,

which is the equation in `QuaConExample.md` §3.3. Its invariants are

        A2      = [ -2731  2299 ; 2299  2189 ] ,     det A2 = -11 263 560 = -265 * 42504
        det M   = -260 148 962 304 = -144 * 42504^2
        centre  c = ( 2496/265 , -3226/265 ) = ( 9.41886..., -12.17358... )
        k       = det M / det A2 = 6 120 576 / 265 = 23 096.5132...

— the two Lemma A / Lemma B predictions confirmed on the nose (the powers of `42504` are the
scaling rules for `det A2` and `det M`). Since `k != 0`, Corollary 3 applies: eigenvalues of `A2`
are `-271 - sqrt(11 337 001)` and `-271 + sqrt(11 337 001)`, i.e. `-3638.0463...` and `+3096.0463...`, and the canonical form is

        u^2 / a^2  -  v^2 / b^2  =  1 ,      a = 2.51964... ,  b = 2.73130... ,  e = 1.47481...

with asymptote slopes the roots of `2189 m^2 + 4598 m - 2731 = 0`, namely `0.4829249...` and
`-2.5834274...` (the discriminant `45 054 240` is not a perfect square, so the slopes are irrational
conjugates). Had `det M` been `0`, `k` would have been `0` and these two slopes would have *been*
the curve.

**Why it looks straight, and why no plot can settle this.** The realized edge of `f*` — verified
below — runs from the triple point outwards, and every point of it sits far out on one branch:
the sampled endpoints are at distances `15.287` and `55.970` from the centre, i.e. `6.1` and `22.2`
semi-axes `a`. A hyperbola is asymptotically straight, so over that stretch the arc is nearly a
line segment: with the two endpoints joined by a chord of length `40.68`, the maximum perpendicular
deviation of the arc from that chord is `0.0479`, about **1 part in 850**. Its curvature at the
midpoint is `1.18e-3`. An eyeball check, or a plot at any reasonable scale, will call this a
straight line. **The determinant will not.**

---

## 9. What was recomputed, and with what

Independent recomputation on 2026-09-10 (MATLAB R2024b on the UBC VPN; Python 3 / SymPy 1.14;
SNAT at `../SNAT`; Lean 4.33 + Mathlib).

* **CCA2, its own code path.** `QuaPol` accepts §4.1's five rows and reports
  `isEdgeContinuous = 1` on all four cevians and `isConvex = 0` — the input is a continuous
  nonconvex PLQ, as claimed. `conjConvexOverPiece` (the KKT active-set construction, not the
  `conj('cplq')` pipeline, which still dies at its Step-3 gate — a separate finding) rebuilds all
  five interior faces, and they agree character for character with §3.2's `I1 ... I5`. From those
  two formulas, the conic above, `b^2-4ac = 45 054 240` and `det M = -260 148 962 304` for the
  cleared form; both Lemma A and Lemma B verified symbolically (`-265` vs `-265`, `-144` vs `-144`).
  At ten points of the conic, `f* = max_k g_k` with `g_3 - g_5 = 0` to machine zero and a margin of
  `0.18` to `0.97` over pieces `1, 2, 4` — so the arc really is an edge of `f*` and both faces
  attain the max there.
* **An oracle sharing no code with any of the above** (40-digit `mpmath`; per piece, the sup over
  the closed triangle by enumeration — interior stationary point when feasible, each closed edge as
  a clamped 1-D concave quadratic, plus the vertices). The `I3|I5` boundary of `f*` was located by
  **bisecting on the argmax**, never using the conic; the located points satisfy the conic equation
  to `1e-35`. Three of them span a triangle of area `0.0085` (`2 * area = 0.017 != 0`): the edge
  contains **three non-collinear points**, which a pair of lines could not do. Boundary points were
  found for `s2 = 2, 2.5, ..., 6, 10, 20, 40` and none at `s2 = 80`.
* **SNAT** (its own convexity-agnostic 1-D conjugate engine, no CCA2 code): slicing `f` at fixed
  height `y` and conjugating the 1-D *nonconvex* PLQ slice, then maximising over `y`. At
  `s = (17/5, 5/2)` the outer maximisers sit at `y = 1303/1610` (piece 3) and `y = 91/240`
  (piece 5) and SNAT returns `115291/32200 = g_3(s)` and `3557/960 = g_5(s)` **exactly, as
  rationals**; at `s = (18/5, 5/2)` it returns `116631/32200 = g_3(s)` and `37087/10560 = g_5(s)`.
  The winner swaps between the two points, which is the arc being crossed. Both face formulas
  therefore hold in a second, independent engine.
* **Lean 4 + Mathlib**, compiled clean with no `sorry` (axioms only `propext`, `Classical.choice`,
  `Quot.sound`): Lemma 4 (`det3_of_line_pair`, `det2_of_line_pair` with the identity
  `-((a1 b2 - a2 b1)/2)^2`), Lemmas A and B in division-free polynomial form (`transfer_det2`,
  `transfer_det3`), the derivation of the printed integer conic from the primal entries
  (`dual_coeffs_are_the_printed_conic`), and the two non-factorisation statements
  (`primal_not_a_line_pair`, `dual_not_a_line_pair`).

---

## 10. If you want to build your own

Theorem 6 is a recipe. Pick two positive definite `2x2` integer matrices with

        det(Q_i - Q_j) < 0        (equivalently: neither Q_i >= Q_j nor Q_j >= Q_i)

— easiest by making one stiffer in `x` and the other stiffer in `y`, as `Q3, Q5` are — then pick
`beta_i, beta_j, gamma_i, gamma_j` and check `det M_p != 0`, which is generic. Finally arrange the
geometry so the two pieces **share no segment**: a fan of three or more triangles from a common
apex, with the pair you care about separated by at least one other piece, is the smallest such
arrangement. Then verify, as §9 does, that the two interior cells really do meet — non-degeneracy
of the conic is about the curve, and one still has to check that `f*` actually switches between
those two faces somewhere rather than handing the maximum to a third piece first.

A warning worth repeating: if you compare the two Hessians and get `det(Q_i - Q_j) > 0` you will
find a perfectly genuine **ellipse** arc, as pairs `1|3`, `1|4` and `1|5` do. That is also a refutation of
a parabolic subdivision, and it is in fact what `QuaConExample.md` §2's three-piece minimal example
produces. The hyperbola needs the extra condition that the two Hessians be Loewner-incomparable.

---

_Companion documents: `doc/QuaConExample.md` (the computed conjugate, the 17 faces and the 41
edges), `CONJ_FIELD_PROOF.md` §4 (the input and the `S4` vertex) and its Theorem 3 (the adjacent
case, of which Proposition 5 here is a self-contained proof)._

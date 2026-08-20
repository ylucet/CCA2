# The explicit example: one triangle, two square roots, one coordinate

_2026-08-20. Produced by `.claude/t1RadicandProbe.m` and the run below; kept because the T1
decision (`DECISIONS.md`) turns on it and re-deriving it costs a MATLAB session._

`exactQ` carries **one** quadratic extension `Q(sqrt(d))` and raises when two are mixed. That rule
is right; the field is too small. Here is the smallest input that shows it, and it is not a corner
case — it is one of the two triangles of the A.4/A.5 quadrilateral every Step 3 measurement in this
project is taken on.

## The input

    T = conv{ (5/2, 3/2), (0, 0), (1/2, 1) },   f = x*y

An edge is CONVEX for `x*y` exactly when its slope is positive, and all three are:

    edge (5/2,3/2)-(0,0)     slope 3/5
    edge (0,0)-(1/2,1)       slope 2
    edge (1/2,1)-(5/2,3/2)   slope 1/4

Three convex edges is [COAP] A.5: split along the middle cevian into two two-convex-edge halves,
each of which recurses into A.4. A.4's cevian has slope `-sqrt(mh*mw)` for the two convex edges it
separates, so each half contributes the square root of a DIFFERENT product of slopes — and nothing
makes two such products differ by a square.

## What comes out

`splitTightTriangleSym(T)` returns four sub-triangles:

    T1 = { (0, 0),      (1/2, 1),   (sqrt(30)/6, sqrt(30)/10) }
    T2 = { (1/2, 1),    (sqrt(30)/6, sqrt(30)/10),
                        (sqrt(30)/12 - sqrt(15)/6 + 5/4,  sqrt(30)/20 - sqrt(15)/10 + 3/4) }
    T3 = { (5/2, 3/2),  (1/2, 1),   (5/2 - sqrt(15)/3, 3/2 - sqrt(15)/5) }
    T4 = { (1/2, 1),    (5/2 - sqrt(15)/3, 3/2 - sqrt(15)/5),
                        (sqrt(30)/12 - sqrt(15)/6 + 5/4,  sqrt(30)/20 - sqrt(15)/10 + 3/4) }

**Read T2's third vertex.** Its x-coordinate is

    sqrt(30)/12 - sqrt(15)/6 + 5/4

— a SINGLE NUMBER that needs both extensions. This is not "two cells that happen to live in
different fields and could be kept apart"; there is no representation of that one vertex over any
`Q(sqrt(d))`. `Q(sqrt(15), sqrt(30)) = Q(sqrt(2), sqrt(3), sqrt(5))`, of degree 8 over Q.

And the neighbouring triangle of the same quadrilateral brings a third one:

    conv{(0,0),(2,0),(5/2,3/2)}  ->  S1 = { (5/2, 3/2), (2, 0), (5/2 - sqrt(5)/2, 3/2 - 3*sqrt(5)/10) }
                                     S2 = { (2, 0), (5/2 - sqrt(5)/2, 3/2 - 3*sqrt(5)/10), (0, 0) }

so Step 3's cross-piece max subtracts a `sqrt(5)` cell from a `sqrt(15)` one.

## What `exactQ` does with them

    a = exactQ(0,1,1,1,15)        %  sqrt(15)
    b = exactQ(0,1,1,1,30)        %  sqrt(30)

    a + b   RAISES  exactQ:fieldMismatch
    a * b   RAISES  exactQ:fieldMismatch
            "cannot combine sqrt(15) with sqrt(30) -- exactQ carries ONE quadratic
             extension by design. Whatever needs both is the operation to look at."

The refusal did its job: it named the operation. `a*b` is the interesting one — the product is
`15*sqrt(2)`, perfectly representable, in a THIRD field again.

## The reading

Square roots of squarefree integers are closed under multiplication up to a rational factor
(`sqrt(15)*sqrt(30) = 15*sqrt(2)`), so the family this pipeline generates is not a tower of
arbitrary algebraic numbers — it is the MULTIQUADRATIC field `Q(sqrt(p1), ..., sqrt(pk))` over the
primes actually seen, an element being a rational combination of `sqrt(m)` for squarefree
`m | p1...pk`. Two properties survive the generalisation, and they are the two the pipeline needs:

* those `sqrt(m)` are linearly INDEPENDENT over Q, so an element is zero exactly when every
  coefficient is zero — zero-testing stays exact and cheap, with no floating point;
* sign follows by the same elimination `exactQ.sign` already uses for `k = 1`, applied recursively:
  write `x = a + b*sqrt(p_k)` with `a, b` in the smaller field; if `a` and `b` share a sign that is
  the answer, otherwise compare `a^2` against `b^2 * p_k` in the smaller field.

`exactQ` is a correct implementation of the `k = 1` case. It is the base case of the type to
build, not a dead end.

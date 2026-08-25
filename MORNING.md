# Morning report — 2026-08-24 overnight run

Branch: `overnight/2026-08-24` (it fast-forwards onto `main`). It ran on a branch
despite the standing "commit on main" preference because a parallel `proof/`
session is committing to `main` tonight.

Task: steps 1–8 of the sym-free `conj` plan. `biconj` untouched throughout.

## The finding that re-shaped the plan

**The numeric conjugate path was already 100% sym-free before this run started.**
Counting non-comment call sites, `conjCPLQ`, `conjPieceCPLQ`, `convEnvCPLQ`,
`maxQuaPar`, `clipArcByConic` and `mergeSameQuadFaces` contain **zero**
`sym` / `subs` / `simplify` / `solve` / `isAlways` / `coeffs` between them. Every
symbolic call on the conj route sits behind ONE dispatch: the fallback to Case C.

So the job was never "rewrite Step 3 without the symbolic engine". It was
**shrink the set of inputs that fall back** — and that is what the night did.
`checkConjSymFree.m` (new) measures that set and names the reason for each.

    baseline, 16 fixtures:  SYMBOLIC 2 of 16, both maxQuaPar, 86–112 s each
                            (the numeric route answers in 0.01–1 s)
    end of run, 17:         SYMBOLIC 3 of 17 — the same two, plus max(0,x,y),
                            which the baseline family did not contain because
                            its unbounded fixture was malformed. The unbounded
                            CONVEX family went from "no numeric route at all"
                            to 0.16 s.

## What changed

1. **`Conic` level of the lattice** — `Conic.m`, `RatCon.m`, and `RatPar.m`
   reduced to the `Par` marker. The subdivision axis is now `Pol < Par < Conic`;
   every existing type is unchanged. `Par.m`'s claim that non-parabolic edges
   never arise is corrected in place: it holds for every comparison
   [COAP]/[JOGO] Theorem 6's proof covers, and fails in general.
   The file is `Conic.m`, not `Con.m`: `CON` is a Windows device name, and a
   plain read of `Con.m` hangs on the console.
2. **`ratQ.m`** — exact rational arithmetic on coefficient vectors, two canonical
   forms because a face function has a scale and a conic does not. Integer-valued
   doubles with a 2^53 guard, not int64 (which saturates silently). 17 tests.
3. **`conicMeet.m`** — two conics intersected through their exact integer
   Sylvester resultant; the certificate is exact, only its root is floating
   point. 12 tests.
4. **`conjConvexPolygon.m`** — the main piece. The conjugate of a CONVEX
   quadratic over ANY convex polygon, bounded or not, in closed form, with no
   triangulation and no curved edge; returns a `QuaPol`. 10 tests, against closed
   forms and a refined definition-sup.
5. **`conjCPLQ` rewired** — convex faces go through (4) whole, and the
   `isDomBounded` gate that sent every unbounded domain to the symbolic path is
   gone. Plus `route='numeric'` (refuses to fall back, so a test can pin the
   ROUTE) and a `CCA2_CONJ_FALLBACK` global that records why.
6. **A cross-check that makes the numeric route SAFE** — see "the one that
   nearly got away" below. This is the most important change of the night.
7. **`maxQuaPar` arc splits** — a clip line cutting a cell's parabola twice, and
   a neighbour's vertex inside an arc, both used to raise `notImplemented`. They
   now SUBDIVIDE, by cutting along the line through the split point parallel to
   the parabola's axis, which meets the conic exactly once. Two documented GAPs
   became OK; `maxQuaParTest` 29/0.
8. **`conjSymFreeTest.m`** — pins which ROUTE each shape takes, including the
   remaining fallbacks, so a regression onto the symbolic path is visible.
9. **Retired** — `exactQ` is off the plan (nothing but its own test uses it);
   `TODO.md` marks the degree-≤4 algebraic kernel and the "detected refusal"
   item CANCELLED by the H-form premise rather than deferred.

## The one that nearly got away

Opening the unbounded route exposed a **silent wrong answer**: on a 4-cone fan
the numeric fold returned `2.0` at `s = (-2,-3)` where the definition sup is
`4.5`. The assembly had dropped a cell — 4 cells for a fold that needs many more
— and every probe point of the same fixture in its other orientation was exact,
which is what a dropped *region* looks like: right almost everywhere, wrong on a
set, and silent. It was caught by an existing test whose truth is a closed form.

The fix is a **verification, not a narrower gate**. `f* = max_k (q_k + I_P_k)*`
is an identity and the per-face conjugates are still in hand, so their pointwise
max *is* f\* exactly. The fold is now checked against it and DECLINES on a
disagreement, which routes that input to the symbolic fallback. Inputs that
assemble correctly keep the 0.01–0.1 s numeric answer; the rest are as slow as
they were before, and none is wrong. The check is one-sided by construction: it
can miss a defect, it cannot invent one. It now protects the bounded path too,
which had none.

## What is broken

Nothing red. Four things are **narrowed or named, not closed** — each pinned by a
test or a reproducer, and listed in `TODO.md` as G1–G3:

- **G1** — an indefinite quadratic over a multi-piece POLYGON still falls back.
  The arc split in (7) works and is tested, but the pieces then fail to pair up
  in `assemblePieces`: the subdivision is consistent per face PAIR, not globally.
  Iterating the global pass to a fixed point was tried and did not fix it.
  Attack the MATCHING, not the refusal — `DECISIONS.md` has the design.
- **G2** — an AFFINE face over an UNBOUNDED polygon (`max(0,x,y)`). Its
  conjugate is a support function, `+inf` off a cone, and `maxQuaPar` refuses an
  operand that is not finite everywhere. `TODO.md` prices a third route that
  would cover every piecewise-LINEAR input in one construction.
- **G2b** — `maxQuaPar` drops a cell on some unbounded folds. This is the defect
  the cross-check catches; it is contained, not fixed, and every input it hits
  pays the symbolic path. Reproducer named in `TODO.md`.
- **G3** — a non-convex face over an unbounded polygon; declines by name.

One bounded limit found and pinned rather than fixed: the exact integer layer
cannot express a near-tangency closer than about 1e-3 (the resultant wants
~1e17, the budget is 2^53). It RAISES rather than rounding, by design.

## Needs a decision

Nothing blocking. Two notes:

- **Three tests in `conjCPLQTest` were changed**, and each change is justified in
  a comment at the site. Two asserted an EXACT return class (`QuaPar`) where the
  result is now a `QuaPol` — a `QuaPar` with every edge conic pinned to zero,
  which satisfies the property those tests state they are about ("closed-form
  numerics, not the symbolic fallback") strictly better. Three read the result
  via `g.fnd`, which only the symbolic form has, and now use the file's own
  route-agnostic `evalConjResult`. **No value assertion was weakened.**
- `conjCPLQTest` got slower (38 s → 154 s) because the fixture that exposed G2b
  now takes the symbolic fallback. That is the honest cost of the cross-check.

## Where I stopped

Everything is committed on the branch; `git log --oneline main..HEAD` is the
list. The next item is **G2b** (a wrong answer, contained but not fixed), then
**G1**. `DECISIONS.md` 2026-08-24 has four entries written so that each starts
from the design rather than from the symptom.

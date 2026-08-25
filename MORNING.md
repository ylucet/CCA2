# Morning report — 2026-08-24 overnight run

Branch: **`main`** — and the branch this run started on was a mistake worth
recording. It created `overnight/2026-08-24` to isolate itself from the parallel
`proof/` session, but the two sessions share one working tree and therefore one
HEAD: branching did not isolate this run, it **moved the other session onto this
run's branch**. Both runs' commits are on it. It has been fast-forwarded back
into `main` (a pointer move, no history rewritten, nothing lost) and both
sessions are on `main` again, which is also the standing preference. The branch
pointer is kept.

The real fix, which `proof/MORNING.md` reaches independently: two unattended runs
in one repository need **separate git worktrees**, not separate folders.

The same shared tree produced two other incidents, both resolved: a `git add -A`
here swept two of the proof session's in-progress files into a commit (un-tracked
one commit later, files on disk untouched, and they re-added them properly), and
editing `.claude/suite.sh` while a `bash` process was executing it killed one
slow-bucket run mid-flight with a bogus syntax error.

Task: steps 1–8 of the sym-free `conj` plan. `biconj` untouched throughout.

---

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
    end of run, 17:         SYMBOLIC 3 of 17 — the same two (G1), plus
                            max(0,x,y) (G2), which the baseline family did not
                            contain because its unbounded fixture was malformed.

That "3 of 17" understates the change, so read it with the two lines that matter:

* the unbounded CONVEX family went from **no numeric route at all** to 0.16 s;
* and it was returning **wrong values** when the gate came off, which is the
  next section.

## The one that nearly got away — and it is now fixed

Opening the unbounded route exposed a **silent wrong answer**: a 4-cone fan
returned `2.0` at `s = (-2,-3)` where the definition sup is `4.5`. The assembly
had dropped a cell, and every probe of the same fixture in its other orientation
was exact — the signature of a dropped *region*: right almost everywhere, wrong
on a set, silent. An existing test whose truth is a closed form caught it.

Two things were done, and both are worth keeping:

1. **A cross-check, so the numeric route cannot return an unchecked number.**
   `f* = max_k (q_k + I_P_k)*` is an identity and the per-face conjugates are
   still in hand, so their pointwise max *is* f\* exactly. The fold is compared
   against it and DECLINES on a disagreement, routing that input to the symbolic
   fallback. One-sided by construction: it can miss a defect, it cannot invent
   one. It protects the bounded path too, which had none.
2. **The defect itself, traced to one line.** The minimal reproducer is two
   functions of one variable:

       h1 = max(s1,0)^2/2   (two cells, split by the s2-axis)
       h2 = max(s2,0)^2/2   (two cells, split by the s1-axis)

   Their max must split the first quadrant along `s2 = s1` — 5 cells. It gave 4.
   `{h1=h2}` is `s1²=s2²`, a **pair of lines crossing at the origin**, which is
   the cone's apex and its only vertex. `splitUnboundedAtOneCrossing` took the
   branch direction as the perpendicular to the gradient — and the gradient
   *vanishes* at a line pair's crossing point. It gave up; the caller's tangency
   branch then read the winner at the cell's centroid, which for a cone **is**
   that apex, lying on the curve. The winner came out of the sign of a computed
   zero. Where the gradient vanishes the branches are the **null directions of
   the quadratic part**, available in closed form from its eigen-decomposition —
   that is the fix. A cone whose apex is that singular point is the dual of any
   *wedge* face, so every 4-cone conjugate hit it. The fan now assembles to 8
   cells with error 0 against the definition sup.

## What changed

1. **`Conic` level of the lattice** — `Conic.m`, `RatCon.m`, `RatPar.m` reduced
   to the `Par` marker. The subdivision axis is `Pol < Par < Conic`; every
   existing type is unchanged. `Par.m`'s claim that non-parabolic edges never
   arise is corrected in place. The file is `Conic.m` because `CON` is a Windows
   device name and a plain read of `Con.m` hangs on the console.
2. **`ratQ.m`** — exact rational arithmetic on coefficient vectors, two canonical
   forms because a face has a scale and a conic does not. 17 tests.
3. **`conicMeet.m`** — two conics intersected through their exact integer
   Sylvester resultant; only the root is floating point. 12 tests.
4. **`conjConvexPolygon.m`** — a CONVEX quadratic over ANY convex polygon,
   bounded or not, closed form, no triangulation, returns a `QuaPol`. 10 tests.
5. **`conjCPLQ` rewired** — convex faces go through (4) whole; the `isDomBounded`
   gate is gone; `route='numeric'` refuses to fall back so a test can pin the
   ROUTE; `CCA2_CONJ_FALLBACK` records why.
6. **The fold cross-check** and **the singular-point fix**, above.
7. **`maxQuaPar` arc splits** — a clip line cutting a cell's parabola twice, and
   a neighbour's vertex inside an arc, both used to raise `notImplemented`. They
   now subdivide along the line through the split point parallel to the
   parabola's axis, which meets the conic exactly once. Two more documented GAPs
   became OK.
8. **`conjSymFreeTest.m`** pins which route each shape takes, fallbacks included.
9. **Retired** — `exactQ` is off the plan; the degree-≤4 algebraic kernel and the
   "detected refusal" item are marked CANCELLED by the H-form premise.

## What is broken

Nothing red that this run caused. **fast 298/0, slow 88/0** (the slow bucket
identical to its pre-change baseline), and the `verylong` daily gate was run
three times to settle it:

    baseline b9243d3, -j 2, uncontended     26 pass / 7 fail / 1 timeout
    current tree,     -j 2, uncontended     26 pass / 7 fail / 1 timeout   <- identical
    current tree,     -j 2, contended       25 pass / 8 fail / 1 timeout   <- one extra

So **the night's work changed `verylong` not at all**; its seven failures and its
timeout are pre-existing (`testMaxMultiRegion/testPCE2` among them, which the
handoff already listed). The eighth appeared only in the contended run, comes
back as an ERROR rather than a mis-assertion, and passes standalone: it is a race
on `plqStage`'s shared cache, filed as G7 with the three-line fix. Until that
lands, **read a `--verylong -j N` red against a `-j 1` re-run of that one test
before believing it.**

Two gaps remain, both named in `TODO.md` and pinned by tests that will go green
when they are fixed:

- **G1 — a CELL is missing, and `clipByFace` is where.** An indefinite quadratic
  over a multi-piece polygon. Measured with each operand's own point location,
  captured alongside the pieces so the numbering cannot drift: every point of the
  lens between g1's straight diagonal and g2's arc lies in **g1 face 4 and g2
  face 2**, and the fold produces **no piece with src [4 2]** — while the control
  point one step away, in `[1 2]`, is there. So `clipByFace` returned nothing for
  a pair whose intersection is not empty, and the orphaned arc in
  `assemblePieces` is the symptom rather than the defect. `polyConstraints`
  already skips a curved edge's chord, so the candidates are the operand SWAP at
  the top of `clipByFace`, `clipPolyByConic`, and the three reduction passes.
  **Instrument that one pair.**
- **G2 — an AFFINE face over an UNBOUNDED polygon** (`max(0,x,y)`). Its
  conjugate is a support function, `+inf` off a cone. `TODO.md` prices a route
  that would cover every piecewise-LINEAR input in one construction and never
  enter `maxQuaPar`; it is a new operator, so it was not started.

One bounded limit found and pinned rather than fixed: the exact integer layer
cannot express a near-tangency closer than about 1e-3 (the resultant wants
~1e17, the budget is 2^53). It RAISES rather than rounding, by design.

## Needs a decision

Nothing blocking. Two notes:

- **Three tests in `conjCPLQTest` were changed**, each justified at the site. Two
  asserted an EXACT return class (`QuaPar`) where the result is now a `QuaPol` —
  a `QuaPar` with every edge conic pinned to zero, which satisfies the property
  those tests state they are about strictly better. Three read the result via
  `g.fnd`, which only the symbolic form has, and now use the file's own
  route-agnostic `evalConjResult`. **No value assertion was weakened.**
- `conjCPLQTest` got slower (38 s → 179 s) because the two G1 fixtures take the
  symbolic fallback. That is the honest cost of the cross-check.

## Where I stopped

All of it is committed on `main`. `DECISIONS.md` has six entries dated
2026-08-24, each written so the next attempt starts from the design rather than
the symptom; `SUPPORT_MATRIX.md` §1.2 and §4 are refreshed and §§4.4–4.5 are new;
`TODO.md` opens with the measured gap list.

Next: **G1** (the lens), then **G2**.

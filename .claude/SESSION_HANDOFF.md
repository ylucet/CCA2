# Session Handoff

_Last updated: 2026-07-18T04:30:00Z_

## What happened this session

Attempted the highest-value TODO from the prior handoff: deriving `conjPieceCPLQ`'s closed-form
rational-piece conjugate. Made substantial, numerically-verified mathematical progress but did
**not** land it — the geometric region-assembly has a confirmed remaining bug. Since shipping it
would violate this codebase's own standard of thorough validation before landing (and returns
`Inf` for some dual points, not just a wrong number), **the code was reverted**; `conjPieceCPLQ.m`
is back to its last-committed state (`083de95`), all 147 tests still pass. This handoff exists to
capture the derivation so the next session doesn't have to redo it.

## The math that IS verified (safe to reuse without re-deriving)

Setup: a RatPol triangle face `p(x,y) = N(x,y)/D(x,y)` (N quadratic numerator, `D=[g h k]` a
nonconstant affine denominator) produced by `convEnvCPLQ`'s Appendix A.3 one-convex-edge envelope
(`envelopeFromClassified` case 1). **Structural fact** (proved algebraically from that formula and
confirmed numerically through `convEnvCPLQ`'s own non-orthogonal `bilinearFrame` rotation): `N`
and `D` vanish **together at exactly one triangle vertex** — call it `W`, the vertex opposite the
classified convex edge; the other two vertices `A,B` (CCW: `W -> A -> B -> W`) are the classified
edge's own endpoints, where `p` is smooth and finite.

1. **Decomposition** `p(x,y) = psi0(x,y) + psi1(x,y)^2/psi2(x,y)` (psi_i affine, psi2 = ±D): divide
   N by D eliminating whichever of x,y has the larger-magnitude coefficient in D; the remainder is
   *exactly* a perfect square of an affine function of the other variable (discriminant ~0 — this
   is what the "N,D vanish together at one vertex" structure forces). Closed-form numeric formulas
   for psi0,psi1,psi2's coefficients from N,D are written out (verified: matches `RatPol.evalRational`
   exactly, including through the real `convEnvCPLQ` pipeline with an indefinite, non-bilinear
   quadratic input).
2. **Gradient confined to a parabola**: `grad p(x,y) = s(t(x,y))` where `t=psi1/psi2`,
   `s(t)=(a0+2a1t-a2t^2, b0+2b1t-b2t^2)` (a_i,b_i = psi_i's x,y-coefficients) — p has rank-1 Hessian
   everywhere, like `conjPSDRank1QuadTriangle` but with a genuinely curved (parabolic, not linear)
   gradient image, since `t` is a *projective* pencil-of-lines-through-W parameter, not a rotated
   linear coordinate.
3. **`p` is literally affine with CONSTANT gradient along the WHOLE edges W-A and W-B** (not just
   at their endpoints) — `s(tA)` along all of W-A, `s(tB)` along all of W-B (tA=t(A), tB=t(B)) —
   because those edges are exactly the pencil's own two extreme members.
4. **Value along any cevian**: for s = s(t) exactly, `F*(s(t)) = c2*t^2 - 2*c1*t - c0` (c_i = psi_i's
   constant terms) — verified against dense sup-sampling over the *whole* 2D triangle, not just the
   line.
5. **`Erow` (the "E" face's 2D value formula, for s off the parabola but with x* on edge A-B)**:
   since edge A-B's own direction is **always perpendicular to psi2's gradient (a2,b2)** (an
   affine-invariant incidence fact: psi2's zero-LINE is parallel to edge A-B by construction —
   verified through the non-orthogonal rotation too), the edge-tangent elimination
   `t* = [ex*s1+ey*s2 - (ex*a0+ey*b0)] / [2*(ex*a1+ey*b1)]` (e=(ex,ey)=B-A) is AFFINE in s (not
   quadratic), and `F*(s) = [c2*t*^2-2*c1*t*-c0] + lambda*betaAB` where
   `lambda = <s-s(t*), nAB>` (nAB = edge A-B's unit outward normal) and `betaAB = <nAB,A>`. Verified
   exactly against `g.eval` on the real `QuaPar` object built from this Erow formula (matches to
   1e-6 or better at dozens of test points, including through the real pipeline).
6. **Vertex values**: `fA = N(A)/D(A)`, `fB = N(B)/D(B)` (ordinary, finite), `fW = psi0(W)` (the
   finite removable-singularity limit — `psi1(W)=psi2(W)=0` exactly).
7. **`ConeW = {s : <s,A> - fA <= <s,W> - fW} ∩ {s : <s,B> - fB <= <s,W> - fW}`** — i.e. TWO CLEAN
   HALF-PLANES (`L_A<=L_W` and `L_B<=L_W`), bounded by the FULL lines through A and B respectively
   (not rays) — confirmed against 199/200 random dual points (the one "mismatch" was a test-
   tolerance artifact at a near-tie point). This is the correct, *complete* description of where
   `w` wins — no curved boundary needed for it at all.

## The bug that is NOT fixed (why nothing was landed)

The four regions (`E`, `ConeA`, `ConeB`, `ConeW`) don't assemble into a single valid `QuaPar`
polygon the way `conjPSDRank1QuadTriangle`'s five faces do. Specifically:

- A parabola conic ("crs", from eliminating `t` between `s1(t)` and the `(b2,-a2)`-combination
  `Traw=2*K*t`) was built expecting it to be the E/ConeW divider (mirroring `rank1EdgeQuad`'s own
  boundary-conic trick). It is NOT that — empirically (9+ test points, exact to machine precision)
  **`crs<=0` is instead exactly the condition "t\* (edge-tangent-based) lies in [tA,tB]"**, i.e. it
  duplicates what `edge2`/`edge4` (the straight t\*=tA / t\*=tB lines through sA=s(tA), sB=s(tB))
  already test. It is NOT the E-vs-ConeW divider.
- The genuine E-vs-ConeW divider, derived directly as `Erow(s) - L_W(s)` (both already-verified
  formulas), IS a valid parabola (discriminant exactly 0, vanishes at sA and sB) and DOES correctly
  separate E from ConeW **inside** the `t*` range (verified: positive at a confirmed-E point,
  negative at a confirmed-ConeW point within range). But its sign becomes MEANINGLESS outside the
  range (Erow isn't even a valid candidate there), so it cannot be reused as ConeW's *outer*
  boundary either — at points like `s=(4,4)` (using the session's worked m=0.6,q=-0.3 example
  triangle), this divider's sign flips even though `{LA<=LW,LB<=LW}` (fact 7 above) correctly and
  unambiguously says ConeW.
- So **ConeW's true region is a union of two differently-shaped pieces** — "beaten E,ConeA,ConeB
  within the strip" (bounded by the Erow-minus-LW parabola) and "outside the strip entirely"
  (bounded just by the two straight half-planes from fact 7) — and QuaPar's face membership test is
  a pure AND of half-plane/conic inequalities (`P{i}` chains), with no OR support. Concretely
  reproduced: `conjPieceCPLQ(RatPol(...))` for this piece type returns `Inf` at points like
  `s=(4,4)`/`s=(6,2)` in the worked example (confirmed via `g.eval`), i.e. no face's constraint
  intersection covers them, even though the true `f*(s)` there is finite and equals `L_W(s)`.
  A short-lived attempt to add a 5th duplicate "ConeW_outer" face (same `L_W` formula, bounded only
  by the two half-planes from fact 7) was not completed: giving it a proper closed polygonal
  boundary without re-touching (and risking silently overlapping/overwriting) the already-correct
  E/ConeA/ConeB/ConeW-inner faces needs care not yet worked out.

## Next steps

- **Resume here, don't re-derive**: facts 1-7 above are solid and reusable as-is. The remaining
  work is purely the polygon-assembly puzzle (see previous section) — likely needs either (a) a
  genuinely different topological decomposition than the 4-face `conjPSDRank1QuadTriangle`-style
  guess this session tried, or (b) a 5th/6th duplicate-formula face for `ConeW`'s outer region,
  carefully bounded so it cannot overlap `E`'s territory (which the two half-planes provably don't,
  per fact 7's 199/200 test — the missing piece is just the closed-polygon bookkeeping, e.g.
  whether `edge3`/`edge5` as FULL lines from `sA`/`sB` naturally close up via their own crossing
  point, and whether that crossing point needs to be a genuine 5th dual vertex or can be avoided).
  A dense angular/grid sweep around candidate apex points (as this session did) is the fastest way
  to sanity-check any new hypothesis before writing code.
  See `conjPieceCPLQ.m`'s own file-header TODO (unchanged this session) for the original recipe
  pointer (`cPLQ/plq_1piece.m`'s `conjugateFunction` "type 1" branch) — this session's derivation
  reached considerably further than that recipe (closed numeric formulas throughout, no live `sym`
  needed) but the recipe's own vertex/edge subdifferential construction (`getSubdiffVertexT1` etc.)
  might still hold a hint for the topology puzzle above if consulted again.
- Once a correct construction exists, validate with the SAME rigor as `conjPSDRank1QuadTriangle`
  (numeric sup-sampling across many random triangles/quadratics, an oracle-based `pickRep`, and a
  dedicated `conjPieceCPLQTest` case) before landing — do not skip this given how many false leads
  this session had that "worked" on the first few test points checked.
- Unchanged from before: exact `[LOCATELLI]` citation in DESIGN.md; 2/741 (0.3%) residual
  `maxQuaPar:internal` crashes; `QuaPar.orderEdges`/`createP`'s near-degenerate-triangle error;
  `partialConj` for `'cplq'`/`'pqp'`; `add` for `RatPol`/`RatPar`; conjugate engines `'pqp'`/`'graph'`;
  standalone `RatPol.conj` gap; a proper random-sample stress test of 3-convex-edge split frequency.

## Relevant files

- `conjPieceCPLQ.m` — **unchanged this session** (reverted back to `083de95`'s state after the
  attempt above did not pan out). The rational-piece branch still raises
  `conjPieceCPLQ:notImplemented` exactly as before.
- No other files touched. `cPLQ/` — untracked reference clone, read again this session (Appendix
  B/C.2 in the COAP paper PDF `doc/s10589-026-00781-5.pdf`, pages 28-34, were also read this
  session — confirmed Appendix B is actually about conjugating the *original* bilinear/quadratic
  piece directly (`conjBilinearXYoneCE` etc., already implemented), not the rational envelope
  output — so it does not contain the missing recipe either).

/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.Shapes

/-!
# A computable classifier over ℚ, and the census that validates it

`Quad.kind` is noncomputable: deciding `disc < 0` over `ℝ` needs `Classical`. But
`../CONJ_FIELD_PROOF.md` Theorem 1 says that for **rational** input data the face
functions and the edge conics of `f*` stay rational — only the *vertices* of the
subdivision need algebraic numbers of degree up to 4. So the conic *type* is
genuinely decidable on the real input class, and this file provides the decision
procedure.

## What is here

* `RatQuad` — the same six coefficients over `ℚ`, with `disc`, `det3` and a
  **computable** `kind` returning the seven-way `ConicKind`.
* `RatQuad.kind_toQuad` — the bridge: the rational classifier agrees with the
  real `Quad.kind` on the image of `ℚ`. Without this the computation would be a
  separate object from the theorem.
* The three branches over `ℚ`, and `ratCand`, so a whole `QuaPol`'s candidate
  list can be computed.
* **The census.** The five-piece example of `../CONJ_FIELD_PROOF.md` §4.1 is
  encoded and checked against numbers computed outside Lean, in a different
  language, by a different route.

## The census, and why it is a real check

Three independent agreements, all `#eval`-checked below:

1. `ratCand` of the five-piece example has **7 branches per piece and 23 after
   deduplication** — exactly `../doc/QuaConExample.md` §3.1's count, which was
   obtained by an entirely separate pipeline. Deduplication is doing real work
   here: 35 branches collapse to 23 because adjacent pieces share edge and vertex
   branches, and that only happens if the branch formulas are right.
2. All **ten** curved edges tabulated in §3.3 classify as that table says —
   five ellipses and five hyperbolas.
3. The **four adjacent pairs** of §4.1 have `det3 = 0` exactly, so their interior
   branches tie along crossing lines. That is `Shapes.lean`'s Theorem 3 happening
   in a concrete instance.
-/

namespace QuaConProof

/-! ### The classifier over ℚ -/

/-- A quadratic in two variables with rational coefficients. -/
@[ext]
structure RatQuad where
  /-- coefficient of `s₁²` -/
  a : ℚ
  /-- coefficient of `s₁s₂` -/
  b : ℚ
  /-- coefficient of `s₂²` -/
  c : ℚ
  /-- coefficient of `s₁` -/
  d : ℚ
  /-- coefficient of `s₂` -/
  e : ℚ
  /-- constant coefficient -/
  f : ℚ
  deriving DecidableEq, Repr, Inhabited

namespace RatQuad

instance : Zero RatQuad := ⟨⟨0, 0, 0, 0, 0, 0⟩⟩

@[simp] lemma zero_a : (0 : RatQuad).a = 0 := rfl
@[simp] lemma zero_b : (0 : RatQuad).b = 0 := rfl
@[simp] lemma zero_c : (0 : RatQuad).c = 0 := rfl
@[simp] lemma zero_d : (0 : RatQuad).d = 0 := rfl
@[simp] lemma zero_e : (0 : RatQuad).e = 0 := rfl
@[simp] lemma zero_f : (0 : RatQuad).f = 0 := rfl

/-- `b² - 4ac`, over `ℚ`. -/
def disc (q : RatQuad) : ℚ := q.b ^ 2 - 4 * q.a * q.c

/-- Four times the `3×3` conic determinant, over `ℚ`. -/
def det3 (q : RatQuad) : ℚ :=
  4 * q.a * q.c * q.f + q.b * q.d * q.e - q.a * q.e ^ 2 - q.c * q.d ^ 2 - q.f * q.b ^ 2

/-- The determinant of the Hessian, `4ac - b²`. -/
def hessDet (q : RatQuad) : ℚ := 4 * q.a * q.c - q.b ^ 2

def sub (p q : RatQuad) : RatQuad :=
  ⟨p.a - q.a, p.b - q.b, p.c - q.c, p.d - q.d, p.e - q.e, p.f - q.f⟩

instance : Sub RatQuad := ⟨sub⟩

@[simp] lemma sub_a (p q : RatQuad) : (p - q).a = p.a - q.a := rfl
@[simp] lemma sub_b (p q : RatQuad) : (p - q).b = p.b - q.b := rfl
@[simp] lemma sub_c (p q : RatQuad) : (p - q).c = p.c - q.c := rfl
@[simp] lemma sub_d (p q : RatQuad) : (p - q).d = p.d - q.d := rfl
@[simp] lemma sub_e (p q : RatQuad) : (p - q).e = p.e - q.e := rfl
@[simp] lemma sub_f (p q : RatQuad) : (p - q).f = p.f - q.f := rfl

/-- **The decision procedure.** Computable, unlike `Quad.kind`. -/
def kind (q : RatQuad) : ConicKind :=
  if q = 0 then .wholePlane
  else if q.det3 = 0 then
    if q.disc < 0 then .pointOrEmpty
    else if q.disc = 0 then .linesParallelOrRepeated
    else .crossingLines
  else
    if q.disc < 0 then .ellipse
    else if q.disc = 0 then .parabola
    else .hyperbola

/-- The real quadratic with the same coefficients. -/
def toQuad (q : RatQuad) : Quad := ⟨q.a, q.b, q.c, q.d, q.e, q.f⟩

@[simp] lemma toQuad_disc (q : RatQuad) : (toQuad q).disc = (q.disc : ℝ) := by
  simp only [toQuad, Quad.disc, disc]; push_cast; ring

@[simp] lemma toQuad_det3 (q : RatQuad) : (toQuad q).det3 = (q.det3 : ℝ) := by
  simp only [toQuad, Quad.det3, det3]; push_cast; ring

@[simp] lemma toQuad_eq_zero_iff (q : RatQuad) : toQuad q = 0 ↔ q = 0 := by
  rw [Quad.eq_zero_iff, RatQuad.ext_iff]
  simp only [toQuad, zero_a, zero_b, zero_c, zero_d, zero_e, zero_f, Rat.cast_eq_zero]

/-- **The bridge.** The computable rational classifier and the noncomputable real
one give the same answer. This is what makes the `#eval`s below statements about
the objects the theorems are about. -/
theorem kind_toQuad (q : RatQuad) : (toQuad q).kind = q.kind := by
  unfold Quad.kind kind
  simp only [toQuad_disc, toQuad_det3, toQuad_eq_zero_iff, Rat.cast_eq_zero,
    Rat.cast_lt_zero]

end RatQuad

/-! ### The three branches over ℚ -/

/-- A point of the rational plane. -/
abbrev RatPt := ℚ × ℚ

namespace RatQuad

def eval (q : RatQuad) (x : RatPt) : ℚ :=
  q.a * x.1 ^ 2 + q.b * x.1 * x.2 + q.c * x.2 ^ 2 + q.d * x.1 + q.e * x.2 + q.f

def gradAt (q : RatQuad) (x : RatPt) : RatPt :=
  (2 * q.a * x.1 + q.b * x.2 + q.d, q.b * x.1 + 2 * q.c * x.2 + q.e)

def alongCurv (q : RatQuad) (D : RatPt) : ℚ :=
  2 * q.a * D.1 ^ 2 + 2 * q.b * D.1 * D.2 + 2 * q.c * D.2 ^ 2

end RatQuad

/-- The rational vertex branch. -/
def ratVertexBranch (q : RatQuad) (v : RatPt) : RatQuad :=
  ⟨0, 0, 0, v.1, v.2, -(q.eval v)⟩

/-- The curvature of `q` along the line from `v` to `w`. -/
def ratEdgeCurv (q : RatQuad) (v w : RatPt) : ℚ := q.alongCurv (w.1 - v.1, w.2 - v.2)

/-- The rational edge branch. -/
def ratEdgeBranch (q : RatQuad) (v w : RatPt) : RatQuad :=
  let d1 := w.1 - v.1
  let d2 := w.2 - v.2
  let al := ratEdgeCurv q v w
  let g := q.gradAt v
  let K := g.1 * d1 + g.2 * d2
  ⟨d1 ^ 2 / (2 * al), d1 * d2 / al, d2 ^ 2 / (2 * al),
   v.1 - K * d1 / al, v.2 - K * d2 / al, -(q.eval v) + K ^ 2 / (2 * al)⟩

/-- The rational interior branch. -/
def ratInteriorBranch (q : RatQuad) : RatQuad :=
  let D := q.hessDet
  ⟨q.c / D, -q.b / D, q.a / D, (-2 * q.c * q.d + q.b * q.e) / D,
   (q.b * q.d - 2 * q.a * q.e) / D,
   (q.c * q.d ^ 2 - q.b * q.d * q.e + q.a * q.e ^ 2) / D - q.f⟩

/-- One piece: a vertex list and a quadratic. -/
structure RatPiece where
  /-- the vertices -/
  verts : List RatPt
  /-- the quadratic carried on the piece -/
  q : RatQuad
  deriving Repr

/-- All candidate branches of one piece. -/
def ratBranches (p : RatPiece) : List RatQuad :=
  (p.verts.map (ratVertexBranch p.q))
  ++ ((p.verts.flatMap fun v => p.verts.map fun w => (v, w)).filter
        (fun vw => 0 < ratEdgeCurv p.q vw.1 vw.2)).map (fun vw => ratEdgeBranch p.q vw.1 vw.2)
  ++ (if p.q.hessDet = 0 then [] else [ratInteriorBranch p.q])

/-- The deduplicated candidate list of a whole `QuaPol`. -/
def ratCand (ps : List RatPiece) : List RatQuad :=
  (ps.flatMap ratBranches).dedup

/-! ### The five-piece example of `../CONJ_FIELD_PROOF.md` §4.1

A triangle cut by four cevians. Rows are entered in the raw-Hessian convention
`[Q₁₁ Q₁₂ Q₂₂ β₁ β₂ γ]` of §4.1 and converted here, so `a = Q₁₁/2`, `b = Q₁₂`,
`c = Q₂₂/2`. All five Hessians are positive definite. -/

namespace Census

/-- Mesh vertices, numbered as in §4.1. -/
def mv : List RatPt := [(0, 0), (60, 10), (15, 10), (5, 10), (-5, 10), (-15, 10), (-60, 10)]

def q1 : RatQuad := ⟨3/2, -1, 5/2, 0, 1, 0⟩
def q2 : RatQuad := ⟨3/2, -5, 17/2, 6, -8, 0⟩
def q3 : RatQuad := ⟨15/2, -2, 11/2, 2, -6, 0⟩
def q4 : RatQuad := ⟨7/2, 2, 17/2, 8, -3, 0⟩
def q5 : RatQuad := ⟨11/2, 11, 35/2, 10, 0, 0⟩

def fivePiece : List RatPiece :=
  [⟨[(0,0), (60,10), (15,10)], q1⟩,
   ⟨[(0,0), (15,10), (5,10)],  q2⟩,
   ⟨[(0,0), (5,10),  (-5,10)], q3⟩,
   ⟨[(0,0), (-5,10), (-15,10)], q4⟩,
   ⟨[(0,0), (-15,10), (-60,10)], q5⟩]

/-! #### Check 1 — the branch count  (compiler-checked)

`../doc/QuaConExample.md` §3.1: "Each piece contributes seven cells to its own
conjugate — one where the unconstrained maximiser lies inside the triangle, one
per edge, one per vertex: 35 cells over the five pieces ... 35 cells collapse to
**23 distinct formulas**."

Deduplication is the load-bearing part: 35 collapses to 23 only if adjacent
pieces really do produce *identical* branch coefficients on their shared cevians
and shared vertices, which is a nontrivial consequence of the three formulas.

These two use `native_decide`, so they rest on the Lean **compiler** rather than
the kernel, and carry the extra axiom `Lean.ofReduceBool`. That is confined to
this section; nothing in `Shapes.lean` or `QuaCon.lean` depends on it, and
`#print axioms conj_isQuaCon` remains clean. Everything below Check 1 is
kernel-checked by `norm_num`. -/

/-- Seven branches per piece, deduplicated within the piece. -/
theorem census_sevenPerPiece :
    (fivePiece.map fun p => (ratBranches p).dedup.length) = [7, 7, 7, 7, 7] := by
  native_decide

/-- **23 distinct candidate quadratics**, matching §3.1. -/
theorem census_twentyThree : (ratCand fivePiece).length = 23 := by
  native_decide

/-! #### Check 2 — `interiorBranch q₁` reproduces the published face `I₁`

`../doc/QuaConExample.md` §3.2 lists
`I₁ = 5s₁²/28 + s₁s₂/14 - s₁/14 + 3s₂²/28 - 3s₂/14 + 3/28`. Recomputing it from
the primal `(Q₁, β₁, γ₁)` is an end-to-end check of the interior branch formula
against a number produced by a different pipeline. -/
theorem census_interiorBranch_q1 :
    ratInteriorBranch q1 = ⟨5/28, 1/14, 3/28, -1/14, -3/14, 3/28⟩ := by
  norm_num [ratInteriorBranch, RatQuad.hessDet, q1, RatQuad.ext_iff]

/-! #### Check 3 — the ten curved edges of §3.3

That table gives an exact integer equation and a type for each. The coefficients
below are transcribed from it; the **types** are what this file computes. -/

def I1I3   : RatQuad := ⟨93, 38, 39, -6, -482, -1003⟩
def I1I5   : RatQuad := ⟨415, 418, 319, 4636, -2332, -24104⟩
def I3I5   : RatQuad := ⟨-2731, 4598, 2189, 107420, 9988, -421996⟩
def I4V3   : RatQuad := ⟨17, -4, 7, -3734, -2226, 104172⟩
def E45I4  : RatQuad := ⟨-28, 12, -21, 1680, 6678, -334389⟩
def E45I5  : RatQuad := ⟨-87, 110, -55, 6668, 25300, -1282588⟩
def E12I5  : RatQuad := ⟨5969, 5390, -847, 67532, -22748, -353236⟩
def E13I3  : RatQuad := ⟨152, 256, 17, -176, -1044, -2588⟩
def E34I3  : RatQuad := ⟨128, -12, -45, 14228, 9144, -433732⟩
def E56I5  : RatQuad := ⟨19, 154, -77, -9884, 35420, -2848244⟩

/-- **All ten classify exactly as §3.3's table says**: five ellipses, five
hyperbolas, and **no parabola**. Kernel-checked. -/
theorem census_tenCurvedEdges :
    RatQuad.kind I1I3 = .ellipse ∧ RatQuad.kind I1I5 = .ellipse ∧
    RatQuad.kind I4V3 = .ellipse ∧ RatQuad.kind E45I4 = .ellipse ∧
    RatQuad.kind E45I5 = .ellipse ∧ RatQuad.kind I3I5 = .hyperbola ∧
    RatQuad.kind E12I5 = .hyperbola ∧ RatQuad.kind E13I3 = .hyperbola ∧
    RatQuad.kind E34I3 = .hyperbola ∧ RatQuad.kind E56I5 = .hyperbola := by
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;>
    norm_num [RatQuad.kind, RatQuad.disc, RatQuad.det3, RatQuad.ext_iff,
      I1I3, I1I5, I4V3, E45I4, E45I5, I3I5, E12I5, E13I3, E34I3, E56I5]

/-! #### Check 4 — Theorem 3 on the four adjacent pairs

Pieces `k` and `k+1` share a cevian and `f` is continuous across it, so
`Shapes.lean`'s `det3_interiorBranch_sub_of_factorisation` predicts a degenerate
tie conic. Here it is, computed: `det3 = 0` exactly, in all four cases. -/
theorem census_adjacentPairsDegenerate :
    (ratInteriorBranch q2 - ratInteriorBranch q1).det3 = 0 ∧
    (ratInteriorBranch q3 - ratInteriorBranch q2).det3 = 0 ∧
    (ratInteriorBranch q4 - ratInteriorBranch q3).det3 = 0 ∧
    (ratInteriorBranch q5 - ratInteriorBranch q4).det3 = 0 := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
    norm_num [RatQuad.det3, ratInteriorBranch, RatQuad.hessDet, q1, q2, q3, q4, q5]

/-- and they are **crossing line pairs**, not points: the discriminant is
positive in each case. -/
theorem census_adjacentPairsCrossingLines :
    RatQuad.kind (ratInteriorBranch q2 - ratInteriorBranch q1) = .crossingLines ∧
    RatQuad.kind (ratInteriorBranch q3 - ratInteriorBranch q2) = .crossingLines ∧
    RatQuad.kind (ratInteriorBranch q4 - ratInteriorBranch q3) = .crossingLines ∧
    RatQuad.kind (ratInteriorBranch q5 - ratInteriorBranch q4) = .crossingLines := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;>
    norm_num [RatQuad.kind, RatQuad.disc, RatQuad.det3, ratInteriorBranch,
      RatQuad.hessDet, RatQuad.ext_iff, q1, q2, q3, q4, q5]

/-! #### Check 5 — the parabola witness of panel (3b)

`q = x₁x₂` on `(0,0), (1,1), (2,0)` — an indefinite piece. Its edge branch along
`(0,0)..(1,1)` against the vertex branch at `(2,0)` is a genuine parabola. By
`Shapes.lean` a parabola is always a flat-versus-flat comparison, and
`det3_vertexBranch_sub_edgeBranch_self` says it is non-degenerate exactly when
the vertex is off the edge's line — `(2,0)` is off the line `x₁ = x₂`. -/
def qxy : RatQuad := ⟨0, 1, 0, 0, 0, 0⟩

theorem census_parabolaWitness :
    RatQuad.kind (ratVertexBranch qxy (2, 0) - ratEdgeBranch qxy (0, 0) (1, 1))
      = .parabola := by
  norm_num [RatQuad.kind, RatQuad.disc, RatQuad.det3, ratVertexBranch, ratEdgeBranch,
    ratEdgeCurv, RatQuad.eval, RatQuad.gradAt, RatQuad.alongCurv, RatQuad.ext_iff, qxy]

end Census

end QuaConProof

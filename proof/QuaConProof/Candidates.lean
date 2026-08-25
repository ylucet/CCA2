/-
Copyright (c) 2026 Yves Lucet. All rights reserved.
-/
import QuaConProof.QuaPol

/-!
# The candidate family, and the cells it defines

`PROJECT_PLAN.md` §0.4 and §0.5.

For one piece `T` carrying the quadratic `q`, write `ψ(x) = ⟨s,x⟩ - q(x)`, the
objective whose supremum over `T` the conjugate takes. Where that supremum is
attained, the maximiser is a stationary point of `ψ` restricted to the affine
hull of some face, and there are only three possibilities:

| branch | attained at | shape of the resulting function of `s` |
|---|---|---|
| `vertexBranch q v` | a vertex `v` | affine |
| `edgeBranch q v w` | the stationary point along the line `v..w` | rank-one quadratic |
| `interiorBranch q` | the stationary point in the plane | full-rank quadratic |

Each is a `Quad`, i.e. a quadratic **in `s`**. This file defines them, proves
each really is `ψ` evaluated at the corresponding stationary point
(`vertexBranch_eval`, `edgeBranch_eval`, `interiorBranch_eval` — kernel-checked
identities, so a sign error cannot survive), assembles them into `cand`, and
defines the cells by activity pattern.

## Why pairs of vertices, not edges of the polygon

`branches` ranges over **all ordered pairs** of vertices, not over the edges of
the polytope. That is what lets the whole development avoid polytope face theory,
which mathlib does not have: Carathéodory puts the maximiser in the convex hull
of an affinely independent subset of at most three vertices, and the pair case is
then a pair of *vertices*, with no need to know it is an edge. See
`DECISIONS.md`, 2026-08-21, "Polytope face theory avoided entirely".

## Why `0 < α` on the edge branch

Along the line through `v` and `w`, `ψ` is `g_v(s) + tL(s) - t²α/2` with
`α = δᵀHδ`. For `α > 0` this is a downward parabola with the maximum at
`t* = L/α`, which is the branch. For `α < 0` the restriction of `ψ` to the line
is unbounded above — checked numerically over a thousand random instances — and
for `α = 0` it is affine, so its value on the segment is attained at an endpoint
and is already a vertex branch. Excluding `α ≤ 0` therefore loses nothing.

## Main definitions

* `psi` — the objective `⟨s,x⟩ - q(x)`.
* `vertexBranch`, `edgeBranch`, `interiorBranch` and their attaining points.
* `cand` — the candidate family of a `QuaPol`, a `Finset Quad`.
* `active`, `cell` — the activity pattern at a point, and the cell it cuts out.

## Main results

* `*_eval` — each branch is `ψ` at the stated point.
* `psi_le_edgeBranch` — the edge branch really is the maximum along the line.
* `cell_disjoint`, `iUnion_cell_eq_univ`, `conj_eq_of_mem_cell`,
  `finite_nonempty_cells`, `cell_subset_eqLocus` — everything about the cells
  that does *not* need the selection lemma. These are the conjuncts of the main
  theorem that come out of the construction.
-/

namespace QuaConProof

open scoped Classical

/-- The objective of the conjugate for a single quadratic: `ψ(x) = ⟨s,x⟩ - q(x)`.
`conj` is the supremum of this over the pieces. -/
def psi (q : Quad) (s x : Plane) : ℝ := dot s x - q.eval x

namespace Quad

/-- The gradient of `q` at `x`. With `q = a s₁² + b s₁s₂ + c s₂² + d s₁ + e s₂ + f`
the Hessian is `[[2a, b], [b, 2c]]`, so `∇q(x) = (2a x₁ + b x₂ + d, b x₁ + 2c x₂ + e)`. -/
def gradAt (q : Quad) (x : Plane) : Plane :=
  (2 * q.a * x.1 + q.b * x.2 + q.d, q.b * x.1 + 2 * q.c * x.2 + q.e)

/-- The second derivative of `q` along the direction `δ`, i.e. `δᵀHδ`. -/
def alongCurv (q : Quad) (δ : Plane) : ℝ :=
  2 * q.a * δ.1 ^ 2 + 2 * q.b * δ.1 * δ.2 + 2 * q.c * δ.2 ^ 2

/-- The determinant of the Hessian, `4ac - b²`. It is `-disc`, so the Hessian is
singular exactly when the conic `{q = 0}` is parabolic. -/
def hessDet (q : Quad) : ℝ := 4 * q.a * q.c - q.b ^ 2

lemma hessDet_eq_neg_disc (q : Quad) : q.hessDet = -q.disc := by
  simp only [hessDet, disc]; ring

end Quad

/-! ### The vertex branch -/

/-- The value of `ψ` at a fixed point `v`, as a quadratic in `s`: it is affine,
`s ↦ ⟨s,v⟩ - q(v)`. -/
def vertexBranch (q : Quad) (v : Plane) : Quad :=
  ⟨0, 0, 0, v.1, v.2, -(q.eval v)⟩

@[simp] lemma vertexBranch_eval (q : Quad) (v s : Plane) :
    (vertexBranch q v).eval s = psi q s v := by
  simp only [vertexBranch, Quad.eval, psi, dot]; ring

/-- The vertex branch has no quadratic part: it is affine in `s`. -/
@[simp] lemma vertexBranch_a (q : Quad) (v : Plane) : (vertexBranch q v).a = 0 := rfl
@[simp] lemma vertexBranch_b (q : Quad) (v : Plane) : (vertexBranch q v).b = 0 := rfl
@[simp] lemma vertexBranch_c (q : Quad) (v : Plane) : (vertexBranch q v).c = 0 := rfl

/-! ### The edge branch -/

/-- The curvature of `q` along the line from `v` to `w`. The edge branch exists
exactly when this is positive. -/
def edgeCurv (q : Quad) (v w : Plane) : ℝ := q.alongCurv (w - v)

/-- Along a recession direction the curvature is that of the direction itself. -/
@[simp] lemma edgeCurv_add (q : Quad) (v r : Plane) :
    edgeCurv q v (v + r) = q.alongCurv r := by
  simp only [edgeCurv, add_sub_cancel_left]

/-- The slope of `ψ` at `v` along the line towards `w`. -/
def edgeSlope (q : Quad) (v w s : Plane) : ℝ := dot s (w - v) - dot (q.gradAt v) (w - v)

/-- The stationary point of `ψ` on the **line** through `v` and `w`. -/
noncomputable def edgePoint (q : Quad) (v w s : Plane) : Plane :=
  v + (edgeSlope q v w s / edgeCurv q v w) • (w - v)

/-- The value of `ψ` at the stationary point along the line `v..w`, as a
quadratic in `s`. Its quadratic part has rank one. -/
noncomputable def edgeBranch (q : Quad) (v w : Plane) : Quad :=
  ⟨(w.1 - v.1) ^ 2 / (2 * edgeCurv q v w),
   (w.1 - v.1) * (w.2 - v.2) / edgeCurv q v w,
   (w.2 - v.2) ^ 2 / (2 * edgeCurv q v w),
   v.1 - dot (q.gradAt v) (w - v) * (w.1 - v.1) / edgeCurv q v w,
   v.2 - dot (q.gradAt v) (w - v) * (w.2 - v.2) / edgeCurv q v w,
   -(q.eval v) + dot (q.gradAt v) (w - v) ^ 2 / (2 * edgeCurv q v w)⟩

/-- **Along the line through `v` and `w`, `ψ` is a downward parabola in the
parameter `t`**, with constant term the vertex branch, slope `edgeSlope` and
curvature `edgeCurv`. Everything about the edge branch follows from this one
identity plus a scalar fact about parabolas. -/
lemma psi_along_line (q : Quad) (v w s : Plane) (t : ℝ) :
    psi q s (v + t • (w - v))
      = psi q s v + t * edgeSlope q v w s - t ^ 2 * edgeCurv q v w / 2 := by
  simp only [psi, edgeSlope, edgeCurv, Quad.alongCurv, Quad.gradAt, Quad.eval, dot,
    Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, Prod.fst_sub, Prod.snd_sub,
    smul_eq_mul]
  ring

/-- The scalar fact: a downward parabola `g + tL - t²α/2` with `α > 0` peaks at
`g + L²/(2α)`. The gap is `(L - tα)²/(2α) ≥ 0`. -/
lemma parabola_le_peak {g L α t : ℝ} (hα : 0 < α) :
    g + t * L - t ^ 2 * α / 2 ≤ g + L ^ 2 / (2 * α) := by
  have hpos : (0 : ℝ) < 2 * α := by linarith
  have key : g + L ^ 2 / (2 * α) - (g + t * L - t ^ 2 * α / 2)
      = (L - t * α) ^ 2 / (2 * α) := by
    field_simp; ring
  have hnn : 0 ≤ (L - t * α) ^ 2 / (2 * α) :=
    div_nonneg (sq_nonneg _) hpos.le
  linarith

/-- **The edge branch, in closed form**: the vertex branch plus `L²/(2α)`.

A kernel-checked identity, so the six coefficients of `edgeBranch` cannot be off
by a sign or a factor of two. -/
theorem edgeBranch_eval_formula (q : Quad) (v w s : Plane) (hα : edgeCurv q v w ≠ 0) :
    (edgeBranch q v w).eval s
      = psi q s v + edgeSlope q v w s ^ 2 / (2 * edgeCurv q v w) := by
  have hd1 : (w - v).1 = w.1 - v.1 := rfl
  have hd2 : (w - v).2 = w.2 - v.2 := rfl
  simp only [edgeBranch, Quad.eval, psi, edgeSlope, dot, hd1, hd2]
  field_simp
  ring

/-- **The edge branch is `ψ` at the stationary point along the line.** -/
theorem edgeBranch_eval (q : Quad) (v w s : Plane) (hα : edgeCurv q v w ≠ 0) :
    (edgeBranch q v w).eval s = psi q s (edgePoint q v w s) := by
  rw [edgeBranch_eval_formula q v w s hα, edgePoint, psi_along_line]
  field_simp
  ring

/-- **The edge branch is the maximum of `ψ` along the whole line**, when the
curvature is positive.

This is why `branches` keeps only the pairs with `0 < edgeCurv`: for `α < 0` the
restriction of `ψ` to the line is unbounded above, and for `α = 0` it is affine,
so its maximum over the segment is at an endpoint and is already a vertex
branch. -/
theorem psi_le_edgeBranch (q : Quad) (v w s : Plane) (hα : 0 < edgeCurv q v w) (t : ℝ) :
    psi q s (v + t • (w - v)) ≤ (edgeBranch q v w).eval s := by
  rw [psi_along_line, edgeBranch_eval_formula q v w s (ne_of_gt hα)]
  exact parabola_le_peak hα

/-! ### The interior branch -/

/-- The stationary point of `ψ` in the plane: the solution of `∇q(x) = s`. -/
noncomputable def interiorPoint (q : Quad) (s : Plane) : Plane :=
  ((2 * q.c * (s.1 - q.d) - q.b * (s.2 - q.e)) / q.hessDet,
   (2 * q.a * (s.2 - q.e) - q.b * (s.1 - q.d)) / q.hessDet)

/-- The value of `ψ` at the unconstrained stationary point, as a quadratic in
`s`. Its quadratic part has full rank. -/
noncomputable def interiorBranch (q : Quad) : Quad :=
  ⟨q.c / q.hessDet, -q.b / q.hessDet, q.a / q.hessDet,
   (-2 * q.c * q.d + q.b * q.e) / q.hessDet, (q.b * q.d - 2 * q.a * q.e) / q.hessDet,
   (q.c * q.d ^ 2 - q.b * q.d * q.e + q.a * q.e ^ 2) / q.hessDet - q.f⟩

/-- The interior point really is the stationary point: `∇q` there is `s`. -/
theorem gradAt_interiorPoint (q : Quad) (s : Plane) (hD : q.hessDet ≠ 0) :
    q.gradAt (interiorPoint q s) = s := by
  rw [Prod.ext_iff]
  constructor <;>
  · simp only [Quad.gradAt, interiorPoint]
    rw [← mul_div_assoc, ← mul_div_assoc, ← add_div, div_add' _ _ _ hD,
      div_eq_iff hD, Quad.hessDet]
    ring

/-- **The interior branch is `ψ` at the unconstrained stationary point.** -/
theorem interiorBranch_eval (q : Quad) (s : Plane) (hD : q.hessDet ≠ 0) :
    (interiorBranch q).eval s = psi q s (interiorPoint q s) := by
  have h2 : q.hessDet ^ 2 ≠ 0 := pow_ne_zero 2 hD
  simp only [interiorBranch, interiorPoint, Quad.eval, psi, dot]
  field_simp
  ring_nf
  rw [Quad.hessDet] at *
  ring

/-! ### The candidate family -/

namespace QuaPiece

/-- The directions along which an edge branch is taken: from a vertex towards
another vertex, or from a vertex along a recession direction.

The second family is Phase 7. A maximiser whose face is the line `v + ℝr` has
value `edgeBranch q v (v + r)`, exactly as a maximiser on the line through two
vertices has value `edgeBranch q v w` — `edgeBranch` depends on `v` and on the
direction `w - v` only, so one definition serves both. -/
noncomputable def edgeCandidates (p : QuaPiece) : Finset (Plane × Plane) :=
  (p.verts ×ˢ p.verts) ∪ (p.verts ×ˢ p.rays).image (fun vr : Plane × Plane => (vr.1, vr.1 + vr.2))

/-- The ordered pairs along which `q` is strictly concave-down, so that the edge
branch exists. -/
noncomputable def edgePairs (p : QuaPiece) : Finset (Plane × Plane) :=
  (p.edgeCandidates).filter (fun vw : Plane × Plane => 0 < edgeCurv p.q vw.1 vw.2)

/-- The interior branch, when the Hessian is nonsingular; otherwise nothing. -/
noncomputable def interiorPart (p : QuaPiece) : Finset Quad :=
  haveI := Classical.dec (p.q.hessDet = 0)
  if p.q.hessDet = 0 then ∅ else {interiorBranch p.q}

/-- All candidate branches contributed by one piece: one per vertex, one per
ordered pair of vertices with positive curvature along it, and the interior one
when the Hessian is nonsingular. -/
noncomputable def branches (p : QuaPiece) : Finset Quad :=
  p.verts.image (vertexBranch p.q)
  ∪ (p.edgePairs).image (fun vw : Plane × Plane => edgeBranch p.q vw.1 vw.2)
  ∪ p.interiorPart

lemma vertexBranch_mem_branches {p : QuaPiece} {v : Plane} (hv : v ∈ p.verts) :
    vertexBranch p.q v ∈ p.branches := by
  simp only [branches, Finset.mem_union, Finset.mem_image]
  exact Or.inl (Or.inl ⟨v, hv, rfl⟩)

lemma mem_edgeCandidates_vert {p : QuaPiece} {v w : Plane} (hv : v ∈ p.verts)
    (hw : w ∈ p.verts) : (v, w) ∈ p.edgeCandidates :=
  Finset.mem_union.2 (Or.inl (Finset.mem_product.2 ⟨hv, hw⟩))

lemma mem_edgeCandidates_ray {p : QuaPiece} {v r : Plane} (hv : v ∈ p.verts)
    (hr : r ∈ p.rays) : (v, v + r) ∈ p.edgeCandidates :=
  Finset.mem_union.2 (Or.inr (Finset.mem_image.2
    ⟨(v, r), Finset.mem_product.2 ⟨hv, hr⟩, rfl⟩))

/-- The edge branch of any admissible pair with positive curvature is a
candidate. Both families of `edgeCandidates` go through this one lemma. -/
lemma edgeBranch_mem_branches_of_cand {p : QuaPiece} {v w : Plane}
    (h : (v, w) ∈ p.edgeCandidates) (hcurv : 0 < edgeCurv p.q v w) :
    edgeBranch p.q v w ∈ p.branches :=
  Finset.mem_union.2 (Or.inl (Finset.mem_union.2 (Or.inr
    (Finset.mem_image.2 ⟨(v, w), Finset.mem_filter.2 ⟨h, hcurv⟩, rfl⟩))))

lemma edgeBranch_mem_branches {p : QuaPiece} {v w : Plane} (hv : v ∈ p.verts)
    (hw : w ∈ p.verts) (hcurv : 0 < edgeCurv p.q v w) :
    edgeBranch p.q v w ∈ p.branches :=
  edgeBranch_mem_branches_of_cand (mem_edgeCandidates_vert hv hw) hcurv

/-- **Phase 7.** The edge branch taken from a vertex along a recession
direction is a candidate too. -/
lemma edgeBranch_ray_mem_branches {p : QuaPiece} {v r : Plane} (hv : v ∈ p.verts)
    (hr : r ∈ p.rays) (hcurv : 0 < edgeCurv p.q v (v + r)) :
    edgeBranch p.q v (v + r) ∈ p.branches :=
  edgeBranch_mem_branches_of_cand (mem_edgeCandidates_ray hv hr) hcurv

lemma interiorBranch_mem_branches {p : QuaPiece} (h : p.q.hessDet ≠ 0) :
    interiorBranch p.q ∈ p.branches := by
  refine Finset.mem_union.2 (Or.inr ?_)
  simp only [interiorPart, if_neg h, Finset.mem_singleton]

end QuaPiece

/-- The **candidate family** of a `QuaPol`: every branch of every piece,
deduplicated by `Finset`. Deduplication is load-bearing — two *distinct* members
have a nonzero difference, which is what makes their equality locus a conic
rather than the whole plane. -/
noncomputable def cand (f : QuaPol) : Finset Quad :=
  f.pieces.biUnion QuaPiece.branches

lemma mem_cand {f : QuaPol} {p : QuaPiece} (hp : p ∈ f.pieces) {q : Quad}
    (hq : q ∈ p.branches) : q ∈ cand f :=
  Finset.mem_biUnion.2 ⟨p, hp, hq⟩

/-! ### Cells -/

/-- The **activity pattern** at `s`: which candidates attain the conjugate there. -/
noncomputable def active (f : QuaPol) (s : Plane) : Finset Quad :=
  (cand f).filter (fun q : Quad => ((q.eval s : ℝ) : EReal) = f.conj s)

lemma mem_active_iff {f : QuaPol} {s : Plane} {q : Quad} :
    q ∈ active f s ↔ q ∈ cand f ∧ ((q.eval s : ℝ) : EReal) = f.conj s :=
  Finset.mem_filter

lemma active_subset_cand (f : QuaPol) (s : Plane) : active f s ⊆ cand f :=
  Finset.filter_subset _ _

/-- The **cell** of a given activity pattern. Cells are pairwise disjoint and
cover the plane by construction; `cell f ∅` is the region where the conjugate is
`⊤`, and needs no separate treatment. -/
def cell (f : QuaPol) (S : Finset Quad) : Set Plane := {s | active f s = S}

lemma mem_cell_iff {f : QuaPol} {S : Finset Quad} {s : Plane} :
    s ∈ cell f S ↔ active f s = S := Iff.rfl

/-- Every point lies in the cell of its own activity pattern. -/
lemma mem_cell_active (f : QuaPol) (s : Plane) : s ∈ cell f (active f s) := rfl

/-- **Cells are pairwise disjoint.** Immediate: a point determines its activity
pattern. -/
theorem cell_disjoint (f : QuaPol) {S T : Finset Quad} (h : S ≠ T) :
    Disjoint (cell f S) (cell f T) := by
  rw [Set.disjoint_left]
  intro s hs ht
  exact h (hs ▸ ht ▸ rfl)

/-- **Cells cover the plane.** -/
theorem iUnion_cell_eq_univ (f : QuaPol) : (⋃ S, cell f S) = Set.univ :=
  Set.eq_univ_of_forall fun s => Set.mem_iUnion.2 ⟨active f s, mem_cell_active f s⟩

/-- **On a cell, the conjugate agrees with each of its active quadratics.** -/
theorem conj_eq_of_mem_cell {f : QuaPol} {S : Finset Quad} {q : Quad} (hq : q ∈ S)
    {s : Plane} (hs : s ∈ cell f S) : f.conj s = ((q.eval s : ℝ) : EReal) := by
  rw [mem_cell_iff] at hs
  subst hs
  exact (mem_active_iff.1 hq).2.symm

/-- **Only finitely many cells are nonempty**, since every activity pattern is a
subset of the finite candidate family. -/
theorem finite_nonempty_cells (f : QuaPol) : {S : Finset Quad | (cell f S).Nonempty}.Finite := by
  apply Set.Finite.subset (cand f).powerset.finite_toSet
  rintro S ⟨s, hs⟩
  rw [mem_cell_iff] at hs
  subst hs
  simpa using active_subset_cand f s

/-- **A cell with two active quadratics lies inside their equality locus**, and
that locus is a conic. This is the sixth conjunct of the main theorem, and it
does not depend on the selection lemma. -/
theorem cell_subset_eqLocus {f : QuaPol} {S : Finset Quad} {q₁ q₂ : Quad}
    (h₁ : q₁ ∈ S) (h₂ : q₂ ∈ S) :
    cell f S ⊆ {s : Plane | q₁.eval s = q₂.eval s} := by
  intro s hs
  have e₁ := conj_eq_of_mem_cell h₁ hs
  have e₂ := conj_eq_of_mem_cell h₂ hs
  exact_mod_cast e₁.symm.trans e₂

/-- The same, packaged with the conic conclusion. -/
theorem isConic_eqLocus_of_mem {S : Finset Quad} {q₁ q₂ : Quad}
    (_h₁ : q₁ ∈ S) (_h₂ : q₂ ∈ S) (hne : q₁ ≠ q₂) :
    IsConic {s : Plane | q₁.eval s = q₂.eval s} := isConic_eqLocus hne

/-! ### Sanity checks -/

namespace Sanity

/-- For `q(x) = x₁² + x₂²` the interior branch is `s ↦ (s₁² + s₂²)/4`, the
classical conjugate of `|x|²`. Hand-check: `sup_x ⟨s,x⟩ - |x|²` is at `x = s/2`
with value `|s|²/2 - |s|²/4 = |s|²/4`. -/
example : interiorBranch ⟨1, 0, 1, 0, 0, 0⟩ = ⟨1/4, 0, 1/4, 0, 0, 0⟩ := by
  simp only [interiorBranch, Quad.hessDet]
  norm_num

/-- The vertex branch at `v = (1,2)` of `q = x₁² + x₂²` is `s ↦ s₁ + 2s₂ - 5`. -/
example : vertexBranch ⟨1, 0, 1, 0, 0, 0⟩ (1, 2) = ⟨0, 0, 0, 1, 2, -5⟩ := by
  simp only [vertexBranch, Quad.eval]
  norm_num

/-- Along the `s₁`-axis from `(0,0)` to `(1,0)`, `q = x₁² + x₂²` has curvature
`2`, and the edge branch is `s ↦ s₁²/4` — rank one, as it must be. -/
example : edgeCurv (⟨1, 0, 1, 0, 0, 0⟩ : Quad) (0, 0) (1, 0) = 2 := by
  simp only [edgeCurv, Quad.alongCurv]
  norm_num

example : edgeBranch (⟨1, 0, 1, 0, 0, 0⟩ : Quad) (0, 0) (1, 0) = ⟨1/4, 0, 0, 0, 0, 0⟩ := by
  simp only [edgeBranch, edgeCurv, Quad.alongCurv, Quad.gradAt, Quad.eval, dot]
  norm_num

/-- The Hessian determinant of `q = x₁x₂` is `-1`: indefinite, and nonsingular,
so the interior branch exists. This is the piece that produces a parabolic edge
(`Conic.lean`, `Parab`). -/
example : (⟨0, 1, 0, 0, 0, 0⟩ : Quad).hessDet = -1 := by
  simp only [Quad.hessDet]; norm_num

end Sanity

end QuaConProof

# Session Handoff

_2026-08-24_

## Blocked

Nothing.

## State

- Branch `main` @ `f807af0` — "Phase 7 done: Frank-Wolfe, and conj_isQuaCon
  loses its boundedness hypothesis"
- Pushed: **no** — four commits ahead of `origin/main`, needs `git push`
- Verification: `lake build` green 2026-08-24, no warnings · **0 sorry** ·
  `#print axioms conj_isQuaCon` = `[propext, Classical.choice, Quot.sound]`
- Known reds: none. 13 files in `QuaConProof/`.
- Two `native_decide` counts in `Rational.lean` carry an extra axiom — audited
  in `SORRY_LEDGER.md`, nothing upstream depends on them.

## What changed this session

**Phase 7 is done.** `conj_isQuaCon` no longer assumes bounded pieces.
`QuaPiece` carries `rays : Finset Plane`, `T = convexHull verts + coneHull rays`;
`IsDirRep` in `Bary.lean` runs one induction over vertex and ray supports at
once; `FrankWolfe.lean` proves attainment for a quadratic bounded above on
`conv V + cone R`. Only `conj_ne_top` and `dom_conj_eq_univ` keep `f.Bounded`,
and they are false without it.

## Next

All optional. Ask Yves which, if any.

1. **Write-up** — `PROOF_NOTES.md` mapping Lean lemmas to `../CONJ_FIELD_PROOF.md`;
   then decide paper vs CCA2 note. Do this first if it is heading anywhere.
2. **Conic normal forms** — the only thing that would let the write-up say
   "ellipse" as a proved classification rather than a discriminant sign.
3. **Realisation, remaining half** — that a *single* cell is infinite.

## Files

- `TODO.md` — the three open items, costed
- `DECISIONS.md` — 21 entries; read before re-attempting anything. The two new
  ones say why the H-representation/Farkas route to Frank-Wolfe is still shut,
  and why the curvature dichotomy must be over `conv R`, not the generators
- `SORRY_LEDGER.md` — zero, plus the axiom audit
- `QuaConProof/FrankWolfe.lean` — Phase 7's real content
- `QuaConProof/QuaCon.lean` — `conj_isQuaCon`, plus `Sanity.rayPol`, the
  unbounded witness whose `⊤` cell is inhabited

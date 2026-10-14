# Session Handoff

_2026-10-13_

## Blocked
- Erratum text changes: Yves to review erratum-lean.tex / erratum-diff.pdf
  (erratum repo D:\AI\67191168115b2b8bc24afaab).

## State
- CCA2 `main` @ `292002c` (overnight/2026-10-10 merged)
- Erratum repo `main` @ `93ca39b` (overnight/2026-10-13 merged)
- Pushed: yes, both (GitHub; Overleaf)
- Lean: `lake build` green, 267 theorems axiom-clean (2026-10-13)
- CCA2 MATLAB tests: not run this session (needs VPN)
- Known reds: none

## Next
1. If wanted: planar Minkowski-Weyl (H->V unbounded pieces, V->H with cone)
   to drop the cited Rockafellar 19.6; see lean/README "Not formalized".
2. CCA2: run fast/normal on VPN (merged commits not re-run).

## Files
- `erratum-lean.tex` — revised erratum after Lean
- `erratum-diff.pdf` — latexdiff vs erratum.tex
- `lean/ErratumProofs/Theorem2Subdivision.lean` — main subdivision theorem
- `lean/README.md` — statement map, not-formalized list

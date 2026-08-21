# Session Handoff — conj number-field question

_2026-08-21. Read-only session alongside the port session; its own
`SESSION_HANDOFF.md` was deliberately NOT overwritten._

## Blocked

Nothing. The QuaPar-shape question is answered: the theorem is at fault, not the
example — `doc/QuaConExample.md` §1, folded into `CONJ_FIELD_PROOF.md` §7.4b.

## State

- Branch `main` @ `a548da8` — "The conjugate of a rational-data PLQ needs
  degree-4 algebraic vertices, not Q(sqrt d)"
- Pushed: yes — 2 commits, `4f39478..d42a2f9`, 2026-08-21
- Tests: none run — no `.m` file was touched this session
- Known reds: unchanged; nothing here affects them

## Next

1. Fold the decided parts into `DECISIONS.md` — deliberately not edited here.
   Candidates: degree-4 vertices, Theorem 3, [JOGO] Thm 6 gap, exactQ one short.
2. Vertex-layer type: degree-4 real algebraic kernel (§8 option 5), then H-form.
3. Open: the S4 minimum piece count (§7.5, four not attempted) and whether ONE
   indefinite piece suffices (§7.4).

## Files

- `CONJ_FIELD_PROOF.md` — the whole result; §8 is the data-structure options
- `doc/QuaCon.svg` — the five-row pipeline figure
- §4.1 — the five-piece input, raw-Hessian rows `[Q11 Q12 Q22 b1 b2 gamma]`
- §7.4, §7.5, §5.1 — the other open items, each scoped

# Session Handoff

_Last updated: 2026-08-13_

## What happened this session

**The arc-vs-arc far-field defect is fixed.** It came down to one sentence: *a curved edge is a
bounded ARC and its conic is not*, so `QuaPar.eval`'s "every bounding conic, sign-oriented, ≤ 0"
rule is exact on a parabola's convex side and admits points arbitrarily far away on its concave
side. `QuaPar.chordCuts` derives the missing constraint (the arc's chord) per face. Six other
defects were found and fixed alongside it, each of which had been producing a silently wrong
answer. The work was merged to `main` and pushed.

## Where things stand

- Branch: `main` @ `634c60d` — "Task 3: the slow bucket's new failures are ROUTE changes, not
  regressions" (plus a handoff commit on top)
- Pushed: yes
- **Fast bucket: 200 pass / 1 fail. Slow bucket: 111 pass / 4 fail (was 106/9).** All five reds are
  documented open items, not unknowns:
  - `maxQuaParTest/arcVsArcRefusesAnUnboundedTwoArcSplit` — see next steps
  - `biconjugateTest/biconjugateOverATwoFaceSubdivisionIsTheEnvelope` — §7 defect
  - `conjCPLQTest/step3UnboundedAssemblyMatchesTheTruth` — §8.2(e) over-claim
  - `unboundedFaceTest` ×2 — curved envelope over an unbounded face
- Measured: `sweepMaxQuaParCurvedSplit(20260802,200)` **30 → 59 assembled** of 142 sampled, 0 of
  1031 vertices / 571 midpoints / 3540 interior points wrong. A 397-quadrilateral far-field sweep:
  **7 of 64 wrong → 0 of 64**.
- **`conj` on a multi-face bounded domain now returns a MESHED `QuaPar`**, not the symbolic
  `QuaParCPLQ` — the numeric route completes where it used to fall back. Five tests pinned the old
  routing and were updated after verifying the new answers are exact (error 0 against the
  closed-form sup, where the old route was checked to 2e-3). This matters for SCIP: the mesh is the
  `V/E/F` the bridge reads.
- New tools, all off by default: `MAXQP_ASSERT` (three exact piece invariants), `QUAPAR_VALIDATE`
  (eval raises when two faces admit a point with different values), `verifyMaxIsExactSymbolically`
  (proves `g = max(g1,g2)` over whole regions by closed-form minimisation, no sampling).

## Next steps

1. **`arcVsArcRefusesAnUnboundedTwoArcSplit`** — the last `maxQuaPar` red. A ray split was
   implemented and REVERTED (it turns one seeded shift into a silent wrong answer). **Read
   `DECISIONS.md` first**: six checks are recorded that do NOT separate the good case from the bad,
   so start past them. Open question: where the cut starts, not which direction it takes.
2. **Covering proof** for `verifyMaxIsExactSymbolically` — every per-face check now passes on all
   four arc-vs-arc fixtures, so "no holes" is the only part of Phase 4 still resting on
   `partitionReport`'s sampling. The hole fixed this session is the argument for doing it.
3. **Docs**: `SUPPORT_MATRIX.md` §7/§8 still list the far-field defect as open and don't mention
   `chordCuts` or the new tools (§4.1 is written and current). `FARFIELD_FIX_PLAN.md` describes a
   plan now largely executed.
4. **Then SCIP/QPLIB**, in the order that bites: wire `biconj` into `SCIP/src/cca2ConvexEnvelope.m`
   → expose value+subgradient off `QuaParCPLQ` → fix diagonal terms over a box (`x²−y²`,
   `(x²+y²)/2` on `[0,1]²` still error in the second conjugation) → performance (~40–60 s/term).

## Relevant files

- `DECISIONS.md` — dead ends and refuted reasoning. Read before attacking item 1.
- `TODO.md` — live action items (its header is current; older sections are history).
- `SUPPORT_MATRIX.md` §4.1 — what the defect was and the re-measured sweeps.
- `maxQuaPar.m` — `clipPolyByConic` (corner cuts, `clipAtVertexCrossings`), `splitCell`
  (`splitUnboundedAtOneCrossing`), `splitTwoArcPiece`, `splitAtReflexVertex`, `matchHalfEdges`,
  `assertPiecesWellFormed`.
- `QuaPar.m` — `chordCuts`, `eval`'s validate mode, `orderEdges` (face-local walk).
- `verifyMaxIsExactSymbolically.m`, `pieceRecessionRays.m` — the acceptance and invariant tools.

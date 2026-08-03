# TODO

_Seeded 2026-08-02 at the start of the overnight run, from the task given when it
was launched: "return both components from clipPolyByConic's disconnecting conic
cut, then continue fixing splitCell and assemblePieces until the three arc-vs-arc
tests pass". The repository had no TODO.md; the acceptance criterion is precise
(three named tests green), so the run works from this list._

## Next up

- [ ] `clipPolyByConic`: return BOTH components when a curved cut disconnects an
      unbounded cell. The branch is reached on all three fixtures and currently
      refuses (`maxQuaPar:notImplemented`); before the refusal it silently built
      one cell from two components, which put a hole and an overlap in the
      arrangement.
- [ ] Re-run the three curved fixtures and fix whatever surfaces next inside
      `clipPolyByConic` (the bounded branch has the same one-component
      assumption).
- [ ] `splitCell`: fix whatever the disconnection fix surfaces there.
- [ ] `assemblePieces`: orphan ray / half-edge matching, once the arrangement is
      a genuine partition.
- [ ] All three arc-vs-arc tests green:
      `twoCurvedWhereTheSplitCurveCrossesAnArc`,
      `twoCurvedThatMustAssembleAcrossRays`,
      `twoCurvedAssembleAcrossRaysSecondFixture`.
- [ ] No regression: `maxQuaParTest` full, then the fast and normal buckets.
- [ ] Decide what to do with the investigation diagnostics added while locating
      the defect (`coverageReport`, the enriched orphan error, the orientation
      guard) — keep the ones that earn their place, drop the rest.

## Done recently

## Blocked

# Morning report — 2026-08-02 overnight run

Branch: `overnight/2026-08-02`

Task as given: return both components from `clipPolyByConic`'s disconnecting
conic cut, then continue fixing `splitCell` and `assemblePieces` until the three
arc-vs-arc tests pass.

Acceptance: `maxQuaParTest/twoCurvedWhereTheSplitCurveCrossesAnArc`,
`twoCurvedThatMustAssembleAcrossRays`, `twoCurvedAssembleAcrossRaysSecondFixture`
all green, with no regression elsewhere.

## What changed

- Seeded `TODO.md`. The repository had none, which is normally a stop condition
  for this mode; the task and its acceptance criterion were given explicitly at
  launch, so the run works from a list derived from them rather than stopping.

## What is broken

- At start of run: the three arc-vs-arc tests fail. `maxQuaParTest` 16 pass /
  3 fail. Everything else in the repository is green apart from two pre-existing
  deliberate reds (`unboundedFaceTest`'s two curved-envelope pins) and one
  by-design failure (`biconjugateTest/biconjugateOverATwoFaceSubdivisionIsTheEnvelope`).

## Needs a decision

## Where I stopped

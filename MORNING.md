# Morning report — 2026-11-04 overnight run

Branch: overnight/2026-11-04

Tasks: (A) the triangle's ~20 faces that own no sample point; (B) Example 2 pair 2's remaining
extra faces. TODO.md "Next up — 2026-11-04 overnight run".

## What changed
- A DONE (ed21ccc): the triangle's 20 extra faces were 18 EXACT DUPLICATES (same value, same rows;
  now kept once at assembly) and 2 quadratic faces touching a bridge at the corner (0,0), certified
  with two more subdivision levels beside a certified corner square. 56 -> 36 faces, every one the
  maximum somewhere, values identical on 400000+ points, no holes. The dedupe keeps each face's
  repeated rows (the boundary tracer reads them by position; uniquing them broke the kink's
  boundary check). The c = -2 kink AlgCrv also loses its copies: 22 -> 4 faces.
## What is broken
Nothing known. Start: main 58331a4, fast 826/0, normal 35/0, slow 127/0 (+2 known), Ex2Subdivision 3/0.
## Needs a decision
## Where I stopped
A: suites running on the final version. B: pair 2's 34 open calls are tangent touches at the
cells' corners; a prototype (edge factors divided out of big-integer face rows + corner blow-ups
for them) certifies the first calls, ~2 min each -- measuring all 34 and the kink's cost next.

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
- B (partial): the exact corner chart rewritten on raw BigInteger, ~50x faster, identical
  (d4506cd). Pair 2's remaining faces diagnosed: tangent touches at cell corners. A certificate
  for them was built and measured, NOT adopted (see Needs a decision).
- Example 2 pairs now 20/62/34 (pair 3 lost a face to A's deeper corner pass).

## What is broken
Nothing known. Start: main 58331a4, fast 826/0, normal 35/0, slow 127/0 (+2 known), Ex2Subdivision 3/0.
## Needs a decision
- **B: adopt the corner certificate for Example 2 pair 2?** It drops 62 -> 60 faces (the 30
  non-winning faces are harmless: values identical), but pair 2 builds in ~820 s instead of ~220 s
  and the kink suite runs ~2x longer. Options: (a) leave it (recommended -- the faces cost only
  overlay size, and the overlay is already declined by default); (b) adopt it as an opt-in for
  overlay builds only; (c) research a separator through the touching point (exact value and
  gradient at dyadic corners). Patch: .claude/proto/2026-11-05-top-row-corners.patch.
## Where I stopped
Run finished: A done, B diagnosed and parked on the decision above. Tree green on 5 commits since
main: fast 826/0, normal 35/0, slow 128/0 (+2 known), Ex2Subdivision 3/0.

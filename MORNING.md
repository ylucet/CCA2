# Morning report — 2026-09-19 overnight run

Branch: overnight/2026-09-19

## What changed
- **The bridge face's region, measured rather than constructed** (item 1). Its membership is
  CONSTANT on the cells of the arrangement -- 50000 points, 40 cells, 0 mixed, each bridge entirely
  inside exactly one cell -- so the region half of the problem that defeated two constructions is
  settled. DECISIONS.md 2026-09-19.
- **`cellInterior`'s two-pair case is decided, not skipped** (item 3). It was reachable: four
  concurrent lines drawn from two conjugate pairs leave 8 of 16 sign patterns undecidable in either
  quadratic field alone, and every one of them is an EMPTY cell that was being kept. The certificate
  now runs in the compositum Q(t1,t2) (`ratQ.qs2*`), with an exact zero test and interval bisection
  on exact rational midpoints. 8 of 16 live, none unsure. cellInteriorTest 7/0.

## What is broken

## Needs a decision

**Does the SUBDIVISION axis grow?** A bridge face carries the envelope's value, and where it
competes with a quadratic face of another pair the mesh must draw `FP(x, y, q(x,y)) = 0`, a rational
curve of degree well above two. Measured: 2 of 40 cells on the c = -2 fixture, and the competitors
are faces of different pairs, so it is a genuine max between pairs and not a partition inside one.
A conic subdivision cannot hold that edge. The options:

1. **Grow the subdivision axis** to carry a higher-degree curve -- the `AlgAlg` level set aside on
   2026-09-13 because only the function axis seemed to need it. Biggest change; gives a real mesh,
   a boundary, and the paper's figures. **This is what I would pick**, because every other route
   leaves `biconj` returning something that cannot draw its own domain.
2. **Return the face list** for a curved tie -- what `envelopeEval` already answers from -- and keep
   the mesh for straight/tangential ties only. Cheap, correct, and no boundary.
3. **Mesh only where no bridge competes** and refuse by name otherwise. Narrow and honest, but on
   this fixture it still refuses, so it buys little.

Nothing was built on this tonight; the item is scope-reduced in TODO.md with the region half done.

## Where I stopped (running)
Item 1 is scope-reduced and its residue is back on the list; item 2 (compaction) is next.

Item 1 (the bridge face's region) is under measurement rather than construction, because two
constructions over sign conditions are already refuted. What is measured so far, on the c = -2
fixture:

- the bridge predicate is CONSTANT on the cells of the arrangement -- 50000 points, 40 cells with
  two or more points, 0 mixed, and each bridge face is entirely inside exactly ONE of them and
  outside the other 39. So its region is a union of arrangement cells, which is what neither failed
  construction could get from sign conditions alone.
- the outer-bound rows add NO new curves: the 23 curves are the ones the quadratic faces and corner
  wedges already contribute.

Still to measure before building: whether the WINNER changes inside a cell once bridges are in the
max. If it does, the mesh needs the curve {bridge value = quadratic value}, which is not a conic.
Run started. Items taken from TODO.md "Next up": the bridge face's region, mesh
compaction, and `cellInterior`'s skipped two-pair certificate case.

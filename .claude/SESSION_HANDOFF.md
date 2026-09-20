# Session Handoff

_2026-09-19 (after the overnight run)_

## Blocked
- **Does the subdivision axis grow?** The only item left, and a type decision — see
  MORNING.md. A bridge face carries the envelope's value, so where it competes
  with a quadratic face the mesh must draw a curve of degree well above two.
- `cPLQ/` is gitignored AND shadows tracked `functionNDomain.m` — your call.
- `maxQuaPar:notDegenerate` reachability — needs the symbolic pipeline (hours).

## State
- Branch `overnight/2026-09-19`, 7 commits ahead of `main`, nothing pushed
- Tests: fast 688/0 over 48 suites · slow 98/0 · verylong 37/0 with 8 documented
  `assumeFail` quarantines (both gates re-run 09-19)
- Known reds: none

## Where the envelope stands
The mesh now carries BOTH curved-tie face kinds for what they decide: a corner for
its wedge (coverage), a bridge for its predicate (coverage), with membership asked
once per cell because the predicate is constant on cells. Compaction halved the
face count, 321 -> 160. What remains is the bridge/quadratic BOUNDARY, refused by
name as `envelopeMesh:bridgeCompetes`.

## Next
1. The decision above. With option 1 the work is the `AlgAlg` subdivision level;
   with option 2 `biconj` returns the face list for a curved tie and the mesh
   stays for straight and tangential ones.
2. `dropRedundantRows` is dead code now — delete after review.
3. Pin `envelopeMesh:bridgeCompetes` with a normal-bucket test (it is on the guard
   inventory as `todo`).

## Files
- `compactCells.m` — drop implied rows, merge to a fixpoint; thin faces keep every row
- `cellInterior.m` + `ratQ.qs2*` — emptiness, now including the compositum of two fields
- `envelopeMesh.m` — corners and bridges as coverage; `bridgeCompetes` is the boundary
- `.claude/bridgeCellConstant.m`, `.claude/bridgeWinnerConstant.m` — the night's measurements

# Morning report — 2026-08-02 overnight run

Branch: `overnight/2026-08-02`

Task as given: return both components from `clipPolyByConic`'s disconnecting
conic cut, then continue fixing `splitCell` and `assemblePieces` until the three
arc-vs-arc tests pass.

**Headline: the tests still fail, but the clip and split stages now produce a
VALID PARTITION where they previously overlapped, and the failure is confined to
`assemblePieces`. The task as stated — "return both components" — turned out to
be the wrong fix, and measuring settled that rather than argument.**

## What changed

- Seeded `TODO.md`. The repository had none, which is normally a stop condition
  for this mode; the task and its acceptance criterion were given explicitly at
  launch, so the run worked from a list derived from them rather than stopping.
- **Retracted the hole/overlap evidence** that drove the previous session's
  conclusions. `pieceContainsPt` omitted the arc edge, so for a 3-vertex curved
  piece it tested only two half-planes — an unbounded wedge. A sliver whose
  vertices all lie within 0.3 of (-2, 2.2) was reported as containing
  (-4.875, -0.325). The check now refuses to answer for a curved piece instead
  of producing confident nonsense.
- **Fixed a real bug** in `clipPolyByConic`'s no-crossing decision: it evaluated
  the cut conic at `interiorSample`, the centroid of the VERTICES, which is not
  necessarily inside a cell bounded by an inward-bulging arc — the same
  non-convexity `splitTwoArcPiece` already warns about. It now decides from the
  boundary signs, which are constant when there is no crossing.
- **Replaced the blanket refusal of a disconnecting curved cut with a real
  connectivity test** (flood fill over the survivor). Measured CONNECTED on all
  three fixtures — 708, 344 and 125 grid cells, one component each — so the
  single-cell construction was right all along and the refusal was
  over-conservative.
- **The partition is now OK on all three fixtures**, where before it reported an
  overlap. That is the substantive gain of the night.
- Pieces are tagged with their source `(k,l)` face pair; the orphan-ray error
  lists the other rays and flags any sharing the apex or containing it. Both
  earned their place — the diagnosis below is entirely from them.

## What is broken

- The three arc-vs-arc tests still fail: `maxQuaParTest` 16 pass / 3 fail,
  unchanged in count but for a much smaller reason.
- All three now fail in the SAME place: `assemblePieces`, a boundary ray with no
  matching neighbour.
- Diagnosis so far: the arrangement legitimately contains very thin cells (area
  ~0.008 beside cells of area ~10), and the orphan rays' apexes are exactly those
  cells' corners. The apexes of the unmatched rays differ by ~0.044, which is far
  above `matchHalfEdges`' 1e-3 apex tolerance, so they are genuinely DIFFERENT
  rays rather than one ray seen twice — this is not the near-degenerate cluster
  case that `checkOrphanHalfEdges` already handles for segments. Something about
  the thin-cell neighbourhood is producing a ray whose partner is not there.

## Needs a decision

Nothing blocks on a decision. One finding is worth your attention because it
contradicts the task as given:

- **"Return both components" is not the right fix, and for a genuine separation
  it is not even representable.** If a curved cut separates an unbounded cell,
  each component contains one ray AND is bounded in part by the cutting conic
  running to infinity — an unbounded curved edge, which QuaPar cannot hold
  (`assertCurvedEdgesAreArcs`). So the options were: (a) return two components
  (impossible for the separated case), (b) refuse always (over-conservative —
  measured, these fixtures are connected), (c) test connectivity and refuse only
  a real separation. I took (c), which is reversible and is what the code now
  does.

## Where I stopped

Mid-way through the last `TODO.md` item, `assemblePieces`' ray matching. The next
concrete step is in `TODO.md` under **Next up**, and **Retired hypotheses** there
lists five things already refuted by measurement so they are not retried:

1. the orphan is one physical ray covered to different extents,
2. the cut must be restricted to the arc's span,
3. `clipPolyByConic` emits clockwise cells,
4. `polyL` is non-convex,
5. pair `(6,1)` is spurious — it is a genuine thin cell, and the "zero
   intersection" reading was a grid-resolution artifact (0.0625 spacing cannot
   see a cell of area 0.008). That one is worth remembering: it nearly led to a
   fix that would have deleted legitimate cells.

Everything is committed on `overnight/2026-08-02`; `main` is untouched. No
regression anywhere: the only failures in the repository are the three
arc-vs-arc pins, the two pre-existing curved-envelope pins in
`unboundedFaceTest`, and the one by-design `biconjugateTest` failure.

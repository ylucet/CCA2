# Morning report — 2026-09-04/05 overnight run

Branch: `overnight/2026-09-04` (12 commits, not merged)

Everything green: **fast 477/0/0**, **normal 13/0**, **slow 98/0** (the night's daily gate, run once
at the start). No red tests, nothing quarantined.

## What changed

Ten of the twelve commits are defect fixes, and **every defect was found by a differential test
against something independent** — the legacy pipeline, an LP, or a closed form — not by reading code.

**`conjQ` — five real defects, all in the same corner of the input space.**

- **Three in HALF-PLANE pieces** (a piece bounded by a line: two opposite rays from one point, where
  every vertex and every recession direction lies ON that line, so the piece's own geometry says
  nothing about which side it is on).
  1. The edge's SIDE could not be derived, so no half-plane was recorded and the domain silently
     became the whole plane. Only `F` knows; it is now consulted there.
  2. The OUTWARD normal was re-derived from those same vertices, so one of two opposite rays got the
     wrong side and its KKT multiplier condition was flipped. That is what left cells overlapping.
  3. A half-plane's RECESSION CONE is two-dimensional and no two rays generate it, so a null
     direction into it was missed and a `+infinity` sup came out finite.

  Result: on `QuaPol.examples{12}`, **28 value and 91 domain disagreements with the legacy
  conjugate went to 0 and 0**.
- **A dual domain that is a single POINT was discarded**, so `examples{19}` answered `+infinity`
  everywhere where `f*(0,0) = 0`. Now recovered exactly and emitted with equality sides.
- **A four-piece input crashed** (`MATLAB:badsubscript`): the corner-naming loop used the pre-merge
  cell count. It could not fire until the merge began removing cells, which only started after an
  earlier restructure.

**`biconjQ` — the same audit, three more defects, plus two features.**

- **A non-convex face was answered WRONGLY** (`+Inf` at four points that are in the domain). The
  envelope's DOMAIN is `conv(P)`, not `P`. Unlike the conjugate this could not be fixed by
  splitting: the convex hull of a union is not the union of hulls.
- **FULL-PLANE inputs raised `noFace`** — `QuaPol.energy`, the most elementary convex input there is.
- **UNBOUNDED CONVEX inputs were refused**, 11 of the 13 fixtures the legacy envelope answered and
  this one did not. A third bug surfaced fixing it: the convexity guard compared hull sizes, which
  is meaningless for a cone, so every unbounded piece looked non-convex.
- **NEW: multi-piece envelopes.** The envelope couples its pieces, so it is not a fold — but that
  coupling is just ONE hull over all the lifted vertices. Verified against an LP at 2.2e-15.
- **The unbounded refusal now says WHICH KIND it is** — `minusInfinity` (a correct answer with
  nowhere to be stored) versus `unbounded` (a finite answer simply not computed). Measured: **all 8
  fixtures the legacy envelope answers and this one refuses are the FINITE kind**, so all 8 are
  reachable rather than blocked on a representation. That is the single most useful number for
  whoever takes the item next.
- **NEW: domains of dimension < 2.** A needle is its own envelope; a segment is a 1-D problem, so
  what decides it is the curvature ALONG the segment, not H in the plane.

**New verification machinery, all checked in.**

- `checkQuaConConsistent.m` — **no two cells may overlap while carrying different functions.** Such
  a mesh does not define a function at all: `eval` resolves it by first match, so the answer depends
  on cell ORDER, and the failure is silent. **2 of 29 conjugates violated it before tonight's fixes,
  0 of 29 after** — and both were caught by its EXACT half, not by sampling. Now asserted across the
  corpus for both operators.
- `.claude/legacy-diff.m` and `.claude/biconj-legacy-diff.m` — differential tests against the legacy
  pipeline over the library's own fixtures.
- `.claude/consistency-sweep.m`, `.claude/biconj-consistency-sweep.m`, `.claude/stress-probe.m`.

## What is broken

**Nothing known-wrong, and nothing quarantined.** The one quarantine opened during the night was
traced and fixed before morning.

One fixture still disagrees with the legacy conjugate — `examples2(9)`, 109 points — and **there the
legacy is the wrong one**: it returns finite values that a sampled sup already exceeds, and the
samples keep climbing as the reach grows. Recorded, not chased.

## Needs a decision

- **I branched, and your memory says not to.** `/overnight` creates a dated branch so the morning
  review is one `git log`; your recorded preference is to commit on `main` for solo research repos.
  I followed the skill since you invoked it explicitly. To take it as if it had been on `main`:
  `git switch main && git merge overnight/2026-09-04 --ff-only`.
- **`TODO.md` has no `## Next up` section**, which is what the skill works from — its live items are
  under dated `##` headings. I worked the top-of-file sections and triaged as I went. Say the word
  and I will restructure it to the convention.
- **Which of `biconjQ`'s three remaining refusals to take first.** My recommendation is the order in
  `TODO.md`: unbounded non-convex (9 fixtures, and the machinery is closest), then multi-piece
  convexity verification, then edge-convex — that last one being where `AlgAlg` finally earns its
  place, so it should not be started before something needs it.

## Where I stopped

Consolidated rather than starting a large feature late in the run: re-ran every probe and sweep after
the night's changes, updated the coverage tables in `TODO.md`, and left both operators' status
written down there rather than in this file.

**Two method notes worth keeping**, both about tests rather than code:

- **A shared assumption between code and oracle is invisible to differential testing.** Every
  multi-face fixture in the repo was malformed the same way (`F(j,:)` is `[left, right]`, so a
  diagonal's two sides need different columns). The `conj` test passed anyway, because the code and
  the oracle both read a face as the polygon its edges bound. Only going through a THIRD route —
  `eval` and `P`'s signs — exposed it.
- **Four of tonight's first five "defects" were my own fixtures or oracles**, not the code: a
  rotated face assignment, an oracle sampling the convex hull rather than the face, an oracle
  evaluating a multi-piece `f` with one piece's formula, and an oracle asserting equality at a
  collinear vertex that is genuinely dominated. Checking a hand-built fixture against `f.eval` before
  trusting it would have saved most of that.

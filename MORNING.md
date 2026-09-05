# Morning report — 2026-09-04 overnight run

Branch: overnight/2026-09-04

## What changed

- **`biconjQ` now does MULTI-PIECE envelopes** (`biconjQ.m`). The envelope couples its pieces, so
  it is not a fold like the conjugate -- but the coupling IS just one hull: each face's graph has
  its lower hull spanned by that face's lifted vertices, so the whole envelope is the lower hull of
  ALL of them at once. Verified against an LP computing the envelope straight from its definition:
  2.2e-15. `biconjQTest` 18/0, fast 471/0.
- **Differential test against the LEGACY conjugate at scale** (`.claude/legacy-diff.m`, new). Over
  the library's own fixture corpus: 16 fixtures both routes answer, 6 the EXACT route reaches and
  the legacy one cannot, 1 the legacy accepts and the exact deliberately refuses (inexact
  coefficients). One real defect found, quarantined -- see below.
- **Fixed: a piece bounded BY A LINE had no half-plane at all.** A half-plane piece has every vertex
  and every recession direction ON its boundary, so its own geometry cannot say which side it is;
  `F` is the only source and is now consulted there. Took `examples(12)` from 91 domain
  disagreements to 26.
- **Found and fixed a malformed-mesh convention in every multi-face fixture in the repo.** `F(j,:)`
  is `[left, right]` of edge `j`, so a diagonal's two sides get different columns: `F(5,:)` must be
  `[2 1]`, not `[1 2]`. Written the other way the mesh describes an unbounded region and
  `QuaPol.eval` reports the wrong piece. The `conj` multi-face test passed anyway because the code
  and the oracle shared the same (correct) reading -- invisible to differential testing, exposed
  only by going through `eval`/`P` as a third route. `conjQTest` still 38/38 after the fix.

## What is broken

- **`QuaPol.examples{12}`: `conjQ`'s fold leaves two cells OVERLAPPING with different values**, and
  `eval`'s first-match returns the smaller (0.20026 where the truth is 0.64930). Quarantined by name
  as `conjQTest/examples12FoldsTwoHalfPlanesAndOVERLAPS_quarantined` (an `assumeFail`, so the bucket
  reports 1 incomplete rather than a red). Ruled out already: it is NOT the merge (overlap survives
  with it disabled) and NOT the half-plane side (fixed separately tonight, moved the domain
  disagreements only). The two overlapping faces carry functions from DIFFERENT pieces, so it is
  `maxQ`'s fold. Full trace in `TODO.md`.

## Needs a decision

- **I branched, and your memory says not to.** `/overnight` creates a dated branch so the morning
  review is one `git log` and abandoning the night is one `git switch`; your recorded preference is
  to commit on `main` for solo research repos. I followed the skill because you invoked it
  explicitly, and it is reversible either way. If you would rather I had stayed on `main`:
  `git switch main && git merge overnight/2026-09-04 --ff-only`.
- **`TODO.md` has no `## Next up` section**, which is what the skill works from. Its live items are
  under dated `##` headings instead. I worked the top-of-file dated sections in order and triaged
  stale items as I went; if you want the `Next up` convention, say so and I will restructure it.

## Where I stopped

_(in progress)_

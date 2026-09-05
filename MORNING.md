# Morning report — 2026-09-04 overnight run

Branch: overnight/2026-09-04

## What changed

- **`biconjQ` now does MULTI-PIECE envelopes** (`biconjQ.m`). The envelope couples its pieces, so
  it is not a fold like the conjugate -- but the coupling IS just one hull: each face's graph has
  its lower hull spanned by that face's lifted vertices, so the whole envelope is the lower hull of
  ALL of them at once. Verified against an LP computing the envelope straight from its definition:
  2.2e-15. `biconjQTest` 18/0, fast 471/0.
- **Fixed THREE defects in HALF-PLANE pieces**, found by chasing the overlap quarantine: the edge
  side could not be derived so the domain became the whole plane; the outward normal was re-derived
  and flipped one of two opposite rays; and a half-plane's recession cone is two-dimensional, so a
  null direction into it was missed and a +infinity sup came out finite. `examples(12)` went from
  28 value + 91 domain disagreements against the legacy route to 0 and 0.
- **New invariant `checkQuaConConsistent.m`**: no two cells may overlap carrying different
  functions. 2 of 29 conjugates violated it before the fixes, 0 of 29 after -- and both were caught
  by its EXACT half, not by sampling. Now asserted across the corpus in `conjQTest`.
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

- **Nothing quarantined, no reds.** The overlap defect found earlier tonight was traced and FIXED
  (three separate bugs, all in half-plane pieces -- see *What changed*).
- Two fixtures still disagree with the LEGACY conjugate about the DOMAIN only (`examples(19)`, 1
  point; `examples2(9)`, 109 points), with zero value disagreements anywhere. Not traced yet, and
  the legacy route is itself unreliable on this family -- it raises its own Step-3-vs-Step-2
  cross-check on `examples(11)` -- so which side is right is an open question, not a known defect.

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

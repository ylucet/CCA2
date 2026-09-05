# Morning report — 2026-09-04 overnight run

Branch: overnight/2026-09-04

## What changed

- **`biconjQ` now does MULTI-PIECE envelopes** (`biconjQ.m`). The envelope couples its pieces, so
  it is not a fold like the conjugate -- but the coupling IS just one hull: each face's graph has
  its lower hull spanned by that face's lifted vertices, so the whole envelope is the lower hull of
  ALL of them at once. Verified against an LP computing the envelope straight from its definition:
  2.2e-15. `biconjQTest` 18/0, fast 471/0.
- **Found and fixed a malformed-mesh convention in every multi-face fixture in the repo.** `F(j,:)`
  is `[left, right]` of edge `j`, so a diagonal's two sides get different columns: `F(5,:)` must be
  `[2 1]`, not `[1 2]`. Written the other way the mesh describes an unbounded region and
  `QuaPol.eval` reports the wrong piece. The `conj` multi-face test passed anyway because the code
  and the oracle shared the same (correct) reading -- invisible to differential testing, exposed
  only by going through `eval`/`P` as a third route. `conjQTest` still 38/38 after the fix.

## What is broken

_(nothing yet)_

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

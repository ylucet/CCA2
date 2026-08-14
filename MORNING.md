# Morning report — 2026-08-14 overnight run

Branch: `overnight/2026-08-13` (created before midnight; the run continued past it)

Task as given: items 1–4 from `TODO.md` — the last `maxQuaPar` red, the covering proof, the
non-compact arc-pieces, and the doc refresh. Explicit instruction: do not pause, and keep attacking
a bug past the usual three-strikes rule.

## Headline

**`maxQuaParTest` is 28 pass / 0 fail, from a `main` baseline of 25 / 1.** All four items are done.
Every bucket re-measured: fast **203 / 0** (was 200 / 1), normal **6 / 0**, slow **111 / 4** —
the slow bucket identical to its recorded baseline, its four failures being the four documented
open items. **No regression anywhere.**

The last red, `arcVsArcRefusesAnUnboundedTwoArcSplit`, had survived two previous attempts. Neither
of them failed for the reason it was thought to:

> **The tooling that judged them was itself broken, in two ways, and silently.** `MAXQP_ASSERT` was
> *crashing* on three of the four arc-vs-arc fixtures, so the invariants that eventually named
> three defects had never run on the inputs that needed them. An invariant that errors is an
> invariant that is off, and nothing was noticing.

MATLAB was unusable for roughly the first half of the run — the UBC licence server is on an
internal domain and the VPN was down, so `matlab -batch` failed with License Manager Error -96.
That half produced documentation and code written blind. The VPN came back; **everything here has
now been run**, and two of the blind changes turned out to be wrong and were fixed on evidence.

## What changed

**Item 1 — `arcVsArcRefusesAnUnboundedTwoArcSplit`, GREEN.** Four defects:

1. `pieceRecessionRays` took the parabola's axis from an exact discriminant. `arcNullDirs` solves
   `d·Q·d' = 0` exactly and returns *nothing* when `b²−4ac` comes out negative — which is what a
   floating-point parabola's `Q` does about half the time, being only semidefinite up to rounding
   (measured `−2.78e-17`). The derived chord was then never emitted, the piece's constraint region
   stayed a slab open at both ends, and the recession-cone invariant refused it. **This is why
   check (5) of the six recorded heuristics never separated the good case from the bad.**
2. `curveAfter ≠ 0` was read as "this edge is curved" in five places, when `boundedPiece` also sets
   it for a *straight* splitting curve. `polyConstraints` emitted **no half-plane at all** for an
   ordinary straight edge — a piece was admitted two units outside itself and answered `0` where
   the truth was `0.35310191`; `pieceStraightEdges` skipped it, blinding every boundary
   minimisation; and two sites called `parabolaArcFrame` on an all-zero conic, which is the crash
   above.
3. With the invariants finally able to run, they named the last defect outright: a whole unbounded
   piece could have its **winner decided by floating-point noise**, because `splitCell`'s unbounded
   "rest" piece can have exactly the two crossing points as its vertices — both, by construction,
   *on* `{f1=f2}` — so neither they nor their centroid can decide it. Now read in the piece's
   recession cone. This is what took the seeded sweep from 16 exact / 0 wrong / 2 errored to
   **17 / 0 / 1**.
4. The ray split itself, restored and gated on all three exact invariants per half rather than on
   heuristics.

**Item 2 — the covering proof, DONE.** `verifyFacesCoverThePlane` decides from the constraint data,
not from probe points, that the faces leave no hole. It took three corrections after its first run
(a row-vs-column indexing slip; a tolerance read off the curve *parameter* rather than the *point*,
which on an arc reached ~0.1 and produced four false over-extensions; and three routes to passing
*vacuously* found by an independent review). `coverProofRejectsBrokenArrangements` now breaks a
certified result three ways and requires a finding each time — a prover that cannot be made to fail
is not a prover.

**Item 3 — was stale, and what was left of it is fixed.** Its representation question had been
answered on 2026-08-13 by deriving the chord per face. What remained was `pieceRecessionRays` still
reading the chord's side off the piece's *vertices* — the rule `DECISIONS.md` records as unsound
one level up — with no gate on when a chord may be emitted at all. Both now settled by the conic.

**Item 4 — docs, DONE.** `SUPPORT_MATRIX.md` §4's table was stale in both directions and every line
number in it predated the arc-vs-arc work by ~1400 lines; all re-derived. New §4.2 (the
verification tools) and §4.3 (the covering proof). §7 and §8 updated. `FARFIELD_FIX_PLAN.md` gets
an OUTCOME section scoring its five phases against what actually happened — its own diagnosis was
wrong, and that gap is the useful part of the record.

**Also fixed:** a newly minted OUTGOING ray was given sign `+1` where `polyConstraints`' convention
needs `−1`, giving the two halves of a split the same outward normal across their shared ray.

## What is broken

**One open case, and it is an error rather than a wrong answer** — a piece that spans TWO sub-arcs
of the same conic, which its single curve slot cannot hold, so the second is represented by its
chord and assembly finds an unpairable half-edge. It was masked until now by the two-arc refusal
upstream.

It is fully diagnosed, and the diagnosis is the useful part, because the symptom misleads: the
error names a straight edge facing an identical *curved* one at distance 8e-16, which reads like a
clip dropping a conic. It is not — on seeded shift `[-2.6434 -1.8066]`, g1's arc is cut **twice**
and one piece spans both sub-arcs. The neighbouring cell's halves are both correct, so
`splitCell`'s crossed-arc restoration is not at fault. `TODO.md` has the captured piece list.

I stopped there rather than writing the fix: it needs a subdivision that `splitTwoArcPiece` cannot
supply as-is (it assumes the two arcs lie on different conics), and that is a design decision, not
a patch.

## Needs a decision

Nothing blocks. Two things worth your attention:

1. **Turn `MAXQP_ASSERT=1` on in the test suite.** Level 1 is cheap, and this run is the argument
   for it: the invariants were off *and* crashing, and three defects lived behind that. I did not
   turn it on myself because it changes the suite's cost and that is your call.
2. **The branch is not merged and not pushed**, per this mode's rules. `main` is untouched.

## Where I stopped

All four items done, committed, and every bucket re-measured — nothing is outstanding.
`.claude/SESSION_HANDOFF.md` is rewritten for the next session; `TODO.md` and `DECISIONS.md` carry
the one open case with its full trace, and the lesson — *when a gate refuses a construction you
have independent reason to believe is right, suspect the gate.*

One retraction worth flagging, since it is the kind of thing that gets copied forward: partway
through I explained the winner defect as `{f1=f2}` being a **pair of parallel lines**, a real
subdivision gap, and wrote a guard for it. Classifying the conic refuted that in one line — its
whole quadratic part is zero, so it is a single straight line. The guard is kept as a cheap exact
backstop and its header says plainly that it now fires on nothing, so nobody reads its existence as
evidence that the case occurs. `DECISIONS.md` records the refutation.

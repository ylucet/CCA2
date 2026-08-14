# Morning report — 2026-08-14 overnight run

Branch: `overnight/2026-08-13` (created before midnight; the run continued past it)

Task as given: items 1–4 from `TODO.md` — the last `maxQuaPar` red, the covering proof, the
non-compact arc-pieces, and the doc refresh. Explicit instruction: do not pause, and keep attacking
a bug past the usual three-strikes rule.

## Headline

**`maxQuaParTest` is 28 pass / 0 fail, from a `main` baseline of 25 / 1.** All four items are done.
The fast bucket is **202 / 0** (was 200 / 1) and the normal bucket **6 / 0**. The slow bucket was
still running when this was written — **check it before merging**; nothing on this branch touches
the symbolic path, but it is the one thing not re-measured.

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
3. With the invariants finally able to run, they named the last defect outright — and it is a new
   **refusal**, not a fix. See "What is broken".
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

**One open case, and it is a refusal rather than a wrong answer.** An unbounded piece can straddle
`{f1 = f2}`: that curve is a degenerate conic, so it can be a **pair of parallel lines**, and a
half lying strictly on one side of the line `splitCell` cut along can still be crossed by the
other — which, if it leaves through the recession cone, contributes no finite boundary crossing for
`splitCell` to find. `assignSide` now detects this exactly (the asymptotic sign along each ray) and
errors instead of returning a silently wrong winner.

Reproducer: seeded shift `[1.4979 3.6486]` of the two-curved fixture, piece `src [2 4]`. The repair
is a `splitCell` that can cut a cell along a second line entering and leaving at infinity;
`DECISIONS.md` (2026-08-03) describes the same shape for a parabola and warns against patching it
with probes, so it was not attempted.

## Needs a decision

Nothing blocks. Two things worth your attention:

1. **Turn `MAXQP_ASSERT=1` on in the test suite.** Level 1 is cheap, and this run is the argument
   for it: the invariants were off *and* crashing, and three defects lived behind that. I did not
   turn it on myself because it changes the suite's cost and that is your call.
2. **The branch is not merged and not pushed**, per this mode's rules. `main` is untouched.

## Where I stopped

All four items done and committed. The slow bucket is the only measurement outstanding.
`.claude/SESSION_HANDOFF.md` is rewritten for the next session; `TODO.md` and `DECISIONS.md` carry
the open case and the lesson — *when a gate refuses a construction you have independent reason to
believe is right, suspect the gate.*

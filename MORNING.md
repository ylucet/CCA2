# Morning report — 2026-09-06 overnight run (second)

Branch: `overnight/2026-09-06-2`

Started from `main` @ `080f878`, tree clean. The first run of the night is already merged to `main`
and pushed, so this branch starts from a published base.

**Daily gates.** `--slow -j 4` and `--verylong -j 4` were started once, chained, in the background
at the top of this run. Both had last been run BEFORE the BigInteger Sturm change landed in
`7c66200`, so neither had seen `ratQ.polyPseudoRem`'s new fallback — a shared low-level routine, so
the gate is worth its cost.

**The three items** come from `TODO.md`'s Next up, in its order:

1. rank-1 EDGE cells — quadratic in one direction, a support function in the other
2. rank-2 INTERIOR cell — bounded and already Case-B-ready
3. confirm the two reclassified `N/R` guards through the wired symbolic path

## What changed

## What is broken

## Needs a decision

## Where I stopped

# Morning report — 2026-09-20 overnight run

Branch: overnight/2026-09-20 (1 commit; the run stopped almost immediately)

## What changed

Nothing was committed beyond the run's own setup. One edit is in the WORKING TREE, uncommitted and
unverified — see below.

## What is broken

Nothing in the repository. The RUN is what broke: both background jobs were killed by the system for
low memory, within minutes of each other —

- the daily SLOW gate, started at the top of the run as rule 9 asks;
- the verification of the `biconjQ` change, which is the only way to know whether that change works.

Claude Code stops background shells when the machine is critically short of memory, and says not to
restart them unattended. Nothing about either job is suspect; the machine simply has no room. Other
projects share it (umbrella rule 3), so this is likely something else's peak rather than ours.

## Needs a decision

**Only whether to re-run.** When the machine has room, the two commands are

    bash .claude/suite.sh --slow
    matlab -batch "addpath('.claude'); biconjRouteProbe"

and the second is what decides the uncommitted edit below. No design question is open: the route's
pieces were each verified this afternoon (`envelopeMeshAlg` against `envelopeEval`, `AlgCrv` against
both), and what is untested is only the wiring that joins them.

## Where I stopped

`biconjQ.m` carries an UNCOMMITTED edit: `assembledEnvelope` now tries `envelopeMesh` first and falls
back to `assembleAlgCrv(envelopeMeshAlg(env))` when — and only when — it raises
`envelopeMesh:bridgeCompetes`. So the RETURN TYPE would say which case the input was: a QuaCon when
every edge is a conic, an AlgCrv when a bridge face's value competes with a quadratic one.

It is left uncommitted deliberately (umbrella rule 10: never commit what has not been run). Note that
`biconjQTest` line 671 currently PINS the old behaviour -- the kinked fixture refusing with
`PLQ:biconjQ:convexPieceCoupling` -- so that test is expected to go red when this lands, and it is
the test that must then change: the contract it asserts is the one being deliberately replaced.

Two items after it, untouched and in TODO.md: the VERTEX table for an AlgCrv (two degree-4 curves
meeting needs a bigger kernel than `conicMeet`'s), and the QuaCon -> AlgCon property move that was
refuted this afternoon.

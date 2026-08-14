# Morning report — 2026-08-13 overnight run

Branch: `overnight/2026-08-13`

Task as given: items 1–4 from `TODO.md` / the session handoff —
1. `arcVsArcRefusesAnUnboundedTwoArcSplit` (the last `maxQuaPar` red),
2. the covering proof for `verifyMaxIsExactSymbolically`,
3. bounded arc-pieces whose constraint region is non-compact,
4. refresh `SUPPORT_MATRIX.md` §7/§8 and `FARFIELD_FIX_PLAN.md`.

Explicit instruction: do not pause for input, and keep attacking a bug past the
usual three-strikes rule.

## THE HEADLINE — MATLAB was unavailable for the whole run

The licence server `SLMS-SMATLABP1.ead.ubc.ca` does not resolve: it is on UBC's internal domain
(`ead.ubc.ca` fails too, `www.ubc.ca` resolves fine), so the **VPN is down**. `matlab -batch`
fails with License Manager Error -96 on every attempt. Connecting the VPN needs an interactive
login, which this mode cannot and should not do.

**So nothing in this run has been executed.** No test was run, no sweep was re-measured, no number
in any document was re-derived from a run. A poller is running in the background and will pick the
licence up the moment it comes back; if it did, the sections below say so.

The run therefore did the work that does not need MATLAB, and stopped short of committing any
change to the maths that would need a test to justify it.

## What changed

- **Item 4 (docs) — DONE.** `SUPPORT_MATRIX.md` §4's table was stale in both directions: three rows
  marked GAP are now handled, several real guards were missing from it, and every line-number
  citation predated the arc-vs-arc work by ~1400 lines, so each pointed at unrelated code. All
  re-derived from the source. §4.1's closing sentence described the ray split that was *reverted*
  rather than the code that is there. New §4.2 (the verification tools, which existed only in the
  session handoff) and §4.3 (the covering proof). §7 now says the far-field defect is closed; §8
  promotes the one remaining `maxQuaPar` case to its own entry and marks the far-field blocker
  resolved. `FARFIELD_FIX_PLAN.md` gets an OUTCOME section scoring its five phases against what
  actually happened — its own diagnosis was wrong, and that gap is the useful part of the record.
- **Item 2 (covering proof) — WRITTEN, NOT RUN.** `verifyFacesCoverThePlane.m`, plus
  `maxQuaParTest/arcVsArcResultsCoverThePlane`. Four checks on the constraint data, each about a
  whole curve rather than probe points; together they force the boundary of the union of the faces
  to be empty. The argument is in the file header. **This code has never been executed.**

## What is broken

- Everything that was red at the start of the run is still red; nothing was fixed and nothing was
  broken, because nothing ran.
- `verifyFacesCoverThePlane.m` and its test are unverified and should be assumed to fail on first
  run until proved otherwise. They are on this branch only.

## Needs a decision

1. **Nothing else can proceed until the VPN is up.** Items 1 and 3 are both "change the geometry,
   then measure" — the reverted ray split was reverted precisely *because* a measurement caught it
   silently returning a wrong value, so writing a replacement without being able to run the seeded
   sweep would repeat the mistake this repository has recorded twice.
2. If MATLAB is likely to be unavailable often, it is worth asking whether a **borrowed licence**
   (MathWorks network licences support offline borrowing for a fixed period) is available on this
   machine. That is a change to the environment, not to the project, so it was not attempted.

## Where I stopped

Items 1 and 3 are **not started as code**, deliberately. What was done for them instead is
analysis, recorded below and in `TODO.md`, so the next session with a working MATLAB starts from a
narrowed search rather than from the six checks already refuted.

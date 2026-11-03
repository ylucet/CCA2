# Morning report — 2026-11-02 overnight run

Branch: overnight/2026-11-02

Tasks: (1) the tangent-face certificate for erratum Example 2 (drop overlay faces tangent to the
winner); (2) guard exercises -- the `todo` ones, then the `symbolic-input` ones as a slow file.

## What changed
- Item 1, step 1 (189a48b): a new exact check that one face lies below another face, used where
  no quadratic fits between them (faces touching along a curve). Also fixed: the squared-factor
  reduction could not remove a 4th power, which is what a tangency produces. Example 2's pair
  face counts 42/98/59 -> 36/96/51. Fast suite 815/0.
- Item 1, step 2 (c0d5064): a face is only compared where it actually carries the value.
  Example 2 pairs -> 33/91/48 faces; the overlay's size bound 11567 -> 7551 (default limit 4000,
  so Example 2 is still declined by default, as decided).
- Item 1, measured end to end: Example 2 fully subdivided (limit lifted) has 6109 faces on this
  branch against 7949 on main, 23% fewer (one run each, shared machine: 1537 s vs 1794 s).
- Item 1, what is left (a92887d, DECISIONS): the remaining extra faces touch the winner exactly at
  cell corners. Two corner fixes were tried and reverted (no gain, or +1 face at 4x the time).
- Item 2 DONE: 20 more error checks now have a test that triggers them (8 in the fast/normal
  suites, 12 in a new slow file), and 6 that no valid input can reach are labelled with the
  reason. Two checks that could never fire were FIXED: a curved edge with no quadratic part
  crashed on a NaN, and an empty region crashed the convex-piece conjugate before its own message.
  4 remain untested (they need big fixtures that don't exist yet).

## What is broken
Nothing known. Tests on this branch: fast 821/0, normal 34/0, slow 122/0 (+2 known-incomplete),
Example 2 long check (limit lifted) 3/0.

## Needs a decision
None.

## Where I stopped
Both items done; the run is finished. Branch not merged or pushed. What could come next on item 1:
the extra faces touch the winner exactly at cell corners, along the cell's edge (DECISIONS
2026-11-02 (third)); a certificate along the edge itself is the idea not yet tried.

# Morning report — 2026-11-02 overnight run

Branch: overnight/2026-11-02

Tasks: (1) the tangent-face certificate for erratum Example 2 (drop overlay faces tangent to the
winner); (2) guard exercises -- the `todo` ones, then the `symbolic-input` ones as a slow file.

## What changed
- Item 1, step 1 (189a48b): a new exact check that one face lies below another face, used where
  no quadratic fits between them (faces touching along a curve). Also fixed: the squared-factor
  reduction could not remove a 4th power, which is what a tangency produces. Example 2's pair
  face counts 42/98/59 -> 36/96/51. Fast suite 815/0.
## What is broken
## Needs a decision
## Where I stopped
Item 1: measuring what Example 2's pair cells still keep after the new check.

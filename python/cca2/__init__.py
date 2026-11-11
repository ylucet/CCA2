"""cca2 -- exact conjugate / biconjugate of bivariate piecewise linear-quadratic functions.

Python port of the EXACT engine of the MATLAB toolbox CCA2 (conjQ, biconjQ). See PORT_PLAN.md and
PORT_NOTES.md in the parent folder for scope and for every deliberate divergence from MATLAB.

Conventions shared by every module (PORT_NOTES.md has the reasons):
  * exact data are Python ``int`` (unbounded) -- never numpy integers, never floats;
  * index VALUES stored in data (E, F, FC, Vname, ...) are 1-BASED with 0 meaning "none",
    exactly as in MATLAB, so fixtures and goldens are verbatim; only list subscripts shift by one;
  * function names inside a module keep the MATLAB spelling (``rat_q.feasible2`` is
    ``ratQ.feasible2``) so the two code bases can be grepped against each other.
"""
from .errors import CCA2Error  # noqa: F401

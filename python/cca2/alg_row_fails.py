"""alg_row_fails -- a CERTIFICATE that a set of algebraic rows cannot all hold anywhere in the
part of a box that lies in a region (algRowFails.m).

STATUS OF THE PORT (PORT_NOTES.md section 5): NOT YET PORTED. algRowFails.m is 1100 lines of
exact sign proofs (resultant-based row polynomials, Yun squarefree parts, blow-ups at corners,
quadratic surds). This stand-in answers False -- "no proof found" -- which is always a SOUND
answer: the callers use a True only to DROP a face or a piece as empty, so with False they keep
it. The resulting AlgCrv evaluates to the same values (an empty face owns no point); it can
carry faces MATLAB proves empty and drops.
"""

PORTED = False


def algRowFails(a, box, inBox, depth=6, factors=None, corners=None):
    return False

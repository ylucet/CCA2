"""Diagnostic: dump the live-cell ordering for a failing conj record, in Python.

Runs conjQ on a specific record index and prints, in order:
  - the final fN (canonical) + each cell's con rows
  - and instruments assembleQuaConCells to also dump the live-cell order
    immediately after the dead-cell drop, and after compactCells.
"""
import json
import sys
from fractions import Fraction

from cca2 import io_matlab
from cca2.conj_q import conjQ
from cca2.assemble_qua_con_cells import assembleQuaConCells
from cca2 import assemble_qua_con_cells as _a

IDX = int(sys.argv[1]) if len(sys.argv) > 1 else 36

recs = io_matlab.load('/home/ylucet/AI/CCA2/python/golden/conj.jsonl')
rec = recs[IDX]
P = io_matlab.quaPolFrom(rec['in'])

# Instrument compactCells to capture the live-cell ordering (fN/fD/FC in order)
_orig_compact = _a.compactCells
def spy_compact(fN, fD, FC, EcQ):
    print("### compactCells input: %d cells" % len(FC))
    for k in range(len(FC)):
        print("  live[%d] fN=%s fD=%s FC=%s" % (k, fN[k], fD[k], FC[k]))
    out = _orig_compact(fN, fD, FC, EcQ)
    fN2, fD2, FC2, info = out
    print("### compactCells output: %d cells" % len(FC2))
    for k in range(len(FC2)):
        print("  out[%d] fN=%s fD=%s FC=%s" % (k, fN2[k], fD2[k], FC2[k]))
    return out
_a.compactCells = spy_compact

g, dom = conjQ(P)
print("\n### FINAL: nf=%d" % g.nf)
for k in range(g.nf):
    print("  f[%d] fN=%s fD=%s" % (k, g.fN[k], g.fD[k]))
print("### expected (MATLAB):")
exp_fN = io_matlab._mat(rec['out']['fN'])
exp_fD = io_matlab._mat(rec['out']['fD'])[0]
for k in range(len(exp_fN)):
    print("  e[%d] fN=%s fD=%s" % (k, exp_fN[k], exp_fD[k]))

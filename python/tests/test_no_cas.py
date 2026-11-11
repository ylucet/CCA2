"""The compute path is CAS-free (the Python twin of conjSymFreeTest / checkConjSymFree.m):
no module of cca2 imports sympy, flint or any other computer-algebra system -- except
_poly.py, the one optional FLINT facade the port plan allows (exact integer kernels with a
pure-Python twin; never symbolic)."""
import ast
import os

import cca2

BANNED = {'sympy', 'flint', 'cypari2', 'sage', 'mpmath'}


def test_cca2_imports_no_cas():
    root = os.path.dirname(cca2.__file__)
    bad = []
    for name in sorted(os.listdir(root)):
        if not name.endswith('.py') or name == '_poly.py':
            continue
        tree = ast.parse(open(os.path.join(root, name)).read())
        for node in ast.walk(tree):
            mods = []
            if isinstance(node, ast.Import):
                mods = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                mods = [node.module]
            for m in mods:
                if m.split('.')[0] in BANNED:
                    bad.append('%s imports %s' % (name, m))
    assert bad == []

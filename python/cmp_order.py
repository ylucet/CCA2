"""Compare per-assemble live-cell fN ordering between Python (/tmp/py_conj36.txt) and
Rust (/tmp/rs_order36.txt). Both group by assemble call, in call order.
"""
import re
from collections import Counter

def parse_py(path):
    calls = []
    cur = None
    for line in open(path):
        m = re.match(r'### compactCells input: (\d+) cells', line)
        if m:
            cur = []
            calls.append(cur)
            continue
        m = re.match(r'\s*live\[\d+\] fN=\[(.*)\] fD=(\d+)', line)
        if m and cur is not None:
            cur.append(tuple(int(x) for x in m.group(1).split(',')))
        if line.startswith('### compactCells output') or line.startswith('### FINAL'):
            cur = None
    return calls

def parse_rs(path):
    calls = []
    cur = None
    for line in open(path):
        m = re.match(r'\[dbg_order\] assemble#\d+: (\d+) live cells', line)
        if m:
            cur = []
            calls.append(cur)
            continue
        m = re.match(r'\[dbg_order\]\s+#\d+ fN=(.*) fD=(\d+)', line)
        if m and cur is not None:
            cur.append(tuple(int(x) for x in m.group(1).split(',')))
    return calls

py = parse_py('/tmp/py_conj36.txt')
rs = parse_rs('/tmp/rs_order36.txt')
print("py calls:", len(py), " rs calls:", len(rs))
for ci in range(max(len(py), len(rs))):
    p = py[ci] if ci < len(py) else None
    r = rs[ci] if ci < len(rs) else None
    if p is None or r is None:
        print("call %d: MISSING one side (py=%s rs=%s)" % (ci, p is not None, r is not None))
        continue
    cp, cr = Counter(p), Counter(r)
    same = (cp == cr)
    order_same = (p == r)
    print("call %d: py_n=%d rs_n=%d multiset_same=%s order_same=%s" % (ci, len(p), len(r), same, order_same))
    if not order_same:
        # find first positional diff
        for i in range(min(len(p), len(r))):
            if p[i] != r[i]:
                print("   FIRST diff at %d: py=%s rs=%s" % (i, p[i], r[i]))
                break

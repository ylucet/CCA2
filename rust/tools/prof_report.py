#!/usr/bin/env python3
"""Report on a stack file written by the `prof` binary (src/bin/prof.rs).

    python3 tools/prof_report.py out.stacks [--tree DEPTH] [--min PCT] [--top N]

Addresses are resolved with `addr2line -f -C -i`; inlined frames are expanded. Prints inclusive and
self percentages per function, and (with --tree) the call tree of the crate's own functions, which
is what shows WHICH LOOP holds the time -- the question a parallelisation study asks.
"""
import collections
import re
import subprocess
import sys


def main():
    a = sys.argv[1:]
    path = a[0]
    depth = int(a[a.index("--tree") + 1]) if "--tree" in a else 0
    minpct = float(a[a.index("--min") + 1]) if "--min" in a else 1.0
    top = int(a[a.index("--top") + 1]) if "--top" in a else 40
    lines = open(path).read().splitlines()
    head = lines[0]
    exe = head.split()[2]
    stacks = [l.split() for l in lines[1:] if l.strip()]
    # return addresses: look one byte back so the CALL line is reported, except for the leaf
    want = set()
    for s in stacks:
        for i, x in enumerate(s):
            v = int(x, 16)
            if v < 1 << 40:
                want.add(v if i == 0 else v - 1)
    want = sorted(want)
    out = subprocess.run(["addr2line", "-f", "-C", "-i", "-a", "-e", exe] + ["%x" % v for v in want], capture_output=True, text=True).stdout.splitlines()
    res = {}
    cur = None
    i = 0
    while i < len(out):
        l = out[i]
        if l.startswith("0x"):
            cur = int(l, 16)
            res[cur] = []
            i += 1
            continue
        fn = re.sub(r"::h[0-9a-f]{16}$", "", l)
        loc = out[i + 1].split("/")[-1].split(" ")[0]
        res[cur].append((fn, loc))
        i += 2
    incl = collections.Counter()
    selft = collections.Counter()
    tree = collections.Counter()
    n = len(stacks)
    for s in stacks:
        names = []  # leaf first
        for i, x in enumerate(s):
            v = int(x, 16)
            if v >= 1 << 40:
                continue
            for fn, _loc in res.get(v if i == 0 else v - 1, []):
                names.append(fn)
        if not names:
            continue
        selft[names[0]] += 1
        for fn in set(names):
            incl[fn] += 1
        if depth:
            own = [f for f in reversed(names) if f.startswith("cca2::") and "{{closure}}" not in f]
            ded = []
            for f in own:
                if not ded or ded[-1] != f:
                    ded.append(f)
            for d in range(1, min(depth, len(ded)) + 1):
                tree[tuple(ded[:d])] += 1
    print(head)
    print("\ninclusive % (crate functions)")
    k = 0
    for fn, c in incl.most_common():
        if fn.startswith("cca2::") and 100.0 * c / n >= minpct:
            print("  %5.1f  %s" % (100.0 * c / n, fn))
            k += 1
            if k >= top:
                break
    print("\nself % (all functions)")
    for fn, c in selft.most_common(top):
        if 100.0 * c / n >= minpct:
            print("  %5.1f  %s" % (100.0 * c / n, fn))
    if depth:
        print("\ncall tree (crate functions, inclusive %)")
        for key in sorted(tree, key=lambda k: tuple((-tree[k[: i + 1]], k[i]) for i in range(len(k)))):
            p = 100.0 * tree[key] / n
            if p >= minpct:
                print("  %s%5.1f  %s" % ("  " * (len(key) - 1), p, key[-1].replace("cca2::", "")))


main()

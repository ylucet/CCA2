"""Delete a named MATLAB function (and its preceding comment block) from a .m file.

Cut = [first line of the function's own leading comment block] .. [its closing `end`].
The closing `end` is found by walking backward from the next `function` line (or EOF) to the
last line that is exactly `end` at the same indentation as the `function` keyword.
"""
import io, re, sys

def find_functions(lines):
    out = []
    for i, l in enumerate(lines):
        m = re.match(r'^(\s*)function\b', l)
        if m:
            out.append((i, len(m.group(1))))
    return out

def delete(path, names, dry=False):
    src = io.open(path, encoding='utf-8', errors='replace').read().split('\n')
    for name in names:
        fns = find_functions(src)
        start = indent = None
        for k, (i, ind) in enumerate(fns):
            if re.match(r'^\s*function\s+(?:\[[^\]]*\]\s*=\s*|[\w~]+\s*=\s*)?' + re.escape(name) + r'\s*[\(\s]', src[i]):
                start, indent, ki = i, ind, k
                break
        if start is None:
            print(f'  !! {path}: {name} NOT FOUND'); continue
        nxt = fns[ki + 1][0] if ki + 1 < len(fns) else len(src)
        endline = None
        for j in range(nxt - 1, start, -1):
            m2 = re.match(r'^(\s*)end\s*[;,]?\s*(%.*)?$', src[j])
            if m2 and len(m2.group(1)) >= indent - 2:
                endline = j; break
        if endline is None:
            print(f'  !! {path}: no closing end found for {name}'); continue
        # absorb a contiguous leading comment block belonging to this function
        top = start
        while top - 1 > 0 and re.match(r'^\s*%', src[top - 1]) and not re.match(r'^\s*end\b', src[top - 1]):
            top -= 1
        n = endline - top + 1
        print(f'  {path}: {name}  lines {top+1}-{endline+1}  ({n} source lines)')
        if not dry:
            del src[top:endline + 1]
    if not dry:
        io.open(path, 'w', encoding='utf-8', newline='').write('\n'.join(src))

if __name__ == '__main__':
    dry = '--dry' in sys.argv
    R = 'C:/Users/ylucet/AI/CCA2/'
    jobs = [
        ('region.m', ['simplifyOpenRegion', 'adjustNormalConeUnB', 'splitmax2', 'slopes2',
                      'getIntersectingFeasiblePts', 'getEdgeNosInf2']),
        ('functionNDomain.m', ['times', 'conjugateExprEdgesT1Poly']),
        ('QuaPar.m', ['evalMatrixForm', 'examples', 'examples2', 'examplesNonconvex', 'examplesDiscontinuous']),
        ('RatPol.m', ['evalMatrixForm', 'examples', 'examples2', 'examplesNonconvex', 'examplesDiscontinuous']),
        ('symbolicFunction.m', ['isConst', 'solveF', 'isParabolic', 'tangentOfSlope']),
    ]
    for f, names in jobs:
        print('==', f)
        delete(R + f, names, dry)

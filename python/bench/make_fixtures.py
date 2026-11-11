"""Write bench/fixtures.json: the exact-engine calls the MATLAB suites make (deduplicated), as
plain data both MATLAB and Python can rebuild. Refusals, inexact inputs and caller-flagged
fixtures are left out: the benchmark times the computation, not the guards."""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
from cca2 import io_matlab  # noqa: E402

out = []
seen = set()
files = [('fast', 'conj.jsonl'), ('fast', 'biconj.jsonl')] + \
        [('slow', os.path.basename(p)) for p in sorted(glob.glob(os.path.join(HERE, '..', 'golden', 'slow_*.jsonl')))]
for bucket, name in files:
    for r in io_matlab.load(os.path.join(HERE, '..', 'golden', name)):
        j = r['in']
        if r['err'] or not j['exact'] or j.get('fIsConvex'):
            continue
        key = (r['fn'], json.dumps([j['V'], j['E'], j['f'], j['F']]))
        if key in seen:
            continue
        seen.add(key)
        out.append({'id': len(out) + 1, 'bucket': bucket, 'src': name, 'fn': r['fn'],
                    'V': j['V'], 'E': j['E'], 'f': j['f'], 'F': j['F'], 'nf': j['nf'],
                    'kind': r['out']['kind']})
json.dump(out, open(os.path.join(HERE, 'fixtures.json'), 'w'))
from collections import Counter
print(len(out), Counter((x['bucket'], x['fn'], x['kind']) for x in out))

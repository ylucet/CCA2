#!/usr/bin/env bash
# export_trace.sh -- record the calls MATLAB makes to a list of top-level functions while a set
# of test suites runs: arguments and results (or the identifier raised), one JSON line per call.
#
#   bash python/tools/export_trace.sh OUTNAME "fn1 fn2 ..." suite1 [suite2 ...]
#       -> python/golden/trace_OUTNAME.jsonl
#
# Same shadow trick as export_golden.sh: for each function a wrapper of the same name, in the
# folder MATLAB runs from, calls a renamed copy of the real file. Nothing in the repo is modified.
set -eu
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$1"; FNS="$2"; shift 2
SH="$(mktemp -d)"
trap 'rm -rf "$SH"' EXIT
for fn in $FNS; do
    sed -E "0,/^function/{/^function/{s/\\b$fn\\b *\\(/${fn}_real(/}}" "$ROOT/$fn.m" > "$SH/${fn}_real.m"
    grep -q "^function.*${fn}_real(" "$SH/${fn}_real.m" || { echo "could not rename $fn"; exit 1; }
    cat > "$SH/$fn.m" <<M
function varargout = $fn(varargin)
    nout = max(nargout, 1);
    try
        [varargout{1:nout}] = ${fn}_real(varargin{:});
    catch e
        try, golden_any('$fn', varargin, {}, e.identifier); catch, end
        rethrow(e)
    end
    try, golden_any('$fn', varargin, varargout, ''); catch e2, disp(e2.message); end
end
M
done
mkdir -p "$ROOT/python/golden"
export CCA2_TRACE_FILE="$ROOT/python/golden/trace_$OUT.jsonl"
rm -f "$CCA2_TRACE_FILE"
list=$(printf "'%s'," "$@"); list="{${list%,}}"
cd "$SH"
matlab -batch "addpath('$ROOT/MATLAB'); addpath('$ROOT/python/tools'); p = strsplit(path, pathsep); m = p(contains(p, [filesep 'toolbox' filesep 'maple'])); if ~isempty(m), rmpath(m{:}); addpath(fullfile(matlabroot,'toolbox','symbolic','symbolic')); end; r = runtests($list); fprintf('COUNTS %d %d %d\n', sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]));" 2>&1 | tail -5
wc -l "$CCA2_TRACE_FILE"

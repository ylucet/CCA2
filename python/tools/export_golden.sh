#!/usr/bin/env bash
# export_golden.sh -- record every call the MATLAB test suites make to the exact-engine entry
# points (conjQ, biconjQ): input QuaPol, exact output, or the named refusal. One JSON line per call.
#
#   bash python/tools/export_golden.sh OUTNAME suite1 [suite2 ...]
#       -> python/golden/OUTNAME.jsonl
#
# HOW. A SHADOW directory holds logging wrappers named conjQ.m / biconjQ.m that call renamed copies
# of the real ones (conjQ_real.m / biconjQ_real.m, the same source with the first line renamed).
# MATLAB runs with the shadow as its CURRENT folder, which outranks the path, so every caller --
# QuaPol.conj, the tests, biconjQ itself -- reaches the wrapper. Nothing in the repo is modified.
set -eu
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="$1"; shift
SH="$(mktemp -d)"
trap 'rm -rf "$SH"' EXIT
sed '1s/conjQ(/conjQ_real(/' "$ROOT/MATLAB/conjQ.m" > "$SH/conjQ_real.m"
sed '1s/biconjQ(/biconjQ_real(/' "$ROOT/MATLAB/biconjQ.m" > "$SH/biconjQ_real.m"
cat > "$SH/conjQ.m" <<'M'
function [g, dom] = conjQ(obj)
    try
        [g, dom] = conjQ_real(obj);
    catch e
        try, golden_log('conjQ', obj, {}, e.identifier); catch, end
        rethrow(e)
    end
    try, golden_log('conjQ', obj, {g, dom}, ''); catch e2, disp(e2.message); end
end
M
cat > "$SH/biconjQ.m" <<'M'
function h = biconjQ(obj, varargin)
    try
        h = biconjQ_real(obj, varargin{:});
    catch e
        try, golden_log('biconjQ', obj, {}, e.identifier); catch, end
        rethrow(e)
    end
    try, golden_log('biconjQ', obj, {h}, ''); catch e2, disp(e2.message); end
end
M
mkdir -p "$ROOT/python/golden"
export CCA2_GOLDEN_FILE="$ROOT/python/golden/$OUT.jsonl"
rm -f "$CCA2_GOLDEN_FILE"
list=$(printf "'%s'," "$@"); list="{${list%,}}"
cd "$SH"
matlab -batch "addpath('$ROOT/MATLAB'); addpath('$ROOT/python/tools'); p = strsplit(path, pathsep); m = p(contains(p, [filesep 'toolbox' filesep 'maple'])); if ~isempty(m), rmpath(m{:}); addpath(fullfile(matlabroot,'toolbox','symbolic','symbolic')); end; r = runtests($list); fprintf('COUNTS %d %d %d\n', sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]));" 2>&1 | tail -5
wc -l "$CCA2_GOLDEN_FILE"

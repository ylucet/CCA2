% A3: are unionIsExact's REFUSALS correct?  Score every captured pair against the definition.
%
% THE DECISIVE TEST, from unionIsExact's own algebra. merge returns M = A' n B' (both regions
% with the shared facet deleted), and M is contained in A u B always. So a point that lies in
% A u B but NOT in M is a point merge would LOSE -- and finding one PROVES the refusal correct.
% Concretely: a point of B violating some constraint of A' (or of A violating one of B').
%
%   refused  + lost point found -> CORRECT       (merge would have been wrong)
%   refused  + none found       -> CONSERVATIVE? (a merge that may have been available)
%   accepted + lost point found -> DEFECT        (merge over-claims; must be zero)
%
% NUMERIC sampling, deliberately: this is a search for a WITNESS, and every witness it reports is
% re-verified exactly before being believed (verifyExact below). Symbolic ptFeasible over tens of
% thousands of candidates is minutes per pair; vectorised double evaluation is milliseconds.
cd(getenv('CCA2DIR'));
addpath(getenv('SPDIR'));
warning('off','symbolic:sym:isAlways:TruthUnknown');
D = dir(fullfile(getenv('CCA2_DUMP_MERGE'), 'mg_*.mat'));
fprintf('captured %d unionIsExact calls\n', numel(D));

tally = containers.Map();
rows = {};
for i = 1:numel(D)
    S = load(fullfile(D(i).folder, D(i).name));
    [z, side] = lostPoint(S.A, S.B, S.ia, S.ib);
    if isempty(z)
        verdict = 'noWitness';
    elseif verifyExact(S.A, S.B, S.ia, S.ib, z)
        verdict = 'LOST';
    else
        verdict = 'witnessUnverified';
    end
    if S.ok, cls = 'accepted'; else, cls = 'refused'; end
    key = sprintf('%-8s %-18s why=%s', cls, verdict, S.why);
    if isKey(tally, key), tally(key) = tally(key) + 1; else, tally(key) = 1; end
    rows{end+1} = struct('file', D(i).name, 'ok', S.ok, 'why', S.why, ...
                         'verdict', verdict, 'z', z, 'side', side); %#ok<AGROW>
end

k = keys(tally);
fprintf('\n--- classification -------------------------------------------------\n');
for t = 1:numel(k)
    fprintf('  %-56s : %d\n', k{t}, tally(k{t}));
end

fprintf('\n--- verdicts -------------------------------------------------------\n');
nDefect = 0; nCorrect = 0; nConserv = 0;
for i = 1:numel(rows)
    r = rows{i};
    if strcmp(r.verdict, 'LOST')
        if r.ok
            nDefect = nDefect + 1;
            fprintf('  DEFECT       %s: merge ACCEPTED but would lose (%g,%g) [%s]\n', ...
                    r.file, r.z(1), r.z(2), r.side);
        else
            nCorrect = nCorrect + 1;
            fprintf('  CORRECT      %s: refused (%s); witness (%g,%g) [%s]\n', ...
                    r.file, r.why, r.z(1), r.z(2), r.side);
        end
    elseif ~r.ok
        nConserv = nConserv + 1;
        fprintf('  CONSERVATIVE? %s: refused (%s), no lost point found\n', r.file, r.why);
    end
end
fprintf('\nDEFECTS %d   CORRECT refusals %d   CONSERVATIVE? %d   of %d calls\n', ...
        nDefect, nCorrect, nConserv, numel(D));

% ------------------------------------------------------------------------------------------
function [z, side] = lostPoint(A, B, ia, ib)
    [z, side] = deal([], '');
    z = probe(B, A, ia);
    if ~isempty(z), side = 'in B, violates a constraint of A'''; return, end
    z = probe(A, B, ib);
    if ~isempty(z), side = 'in A, violates a constraint of B'''; end
end

function z = probe(R, Other, drop)
% A point of R violating some constraint of Other with `drop` removed. Numeric; verified later.
    z = [];
    hR = handles(R);
    keep = true(1, size(Other.ineqs,2)); keep(drop) = false;
    hO = handles(Other, find(keep));
    P = boxSample(R, 60000);
    inR = true(size(P,1),1);
    for t = 1:numel(hR), inR = inR & (hR{t}(P(:,1), P(:,2)) <= 1e-12); end
    P = P(inR,:);
    if isempty(P), return, end
    for t = 1:numel(hO)
        bad = hO{t}(P(:,1), P(:,2)) > 1e-7;
        if any(bad)
            idx = find(bad);
            [~, w] = max(hO{t}(P(idx,1), P(idx,2)));    % the deepest violation, most robust
            z = P(idx(w),:); return
        end
    end
end

function tf = verifyExact(A, B, ia, ib, z)
% Re-ask the region objects themselves: is z genuinely in one region and genuinely outside the
% other's retained constraints? Symbolic, once, on a single point.
    tf = false;
    v = A.vars;
    for pair = {{A, B, ib}, {B, A, ia}}
        R = pair{1}{1}; O = pair{1}{2}; drop = pair{1}{3};
        if ~R.ptFeasible(v, z), continue, end
        for t = 1:size(O.ineqs,2)
            if t == drop, continue, end
            val = double(O.ineqs(t).subsF(v, z).f);
            if isfinite(val) && val > 1e-7, tf = true; return, end
        end
    end
end

function h = handles(R, idx)
    if nargin < 2, idx = 1:size(R.ineqs,2); end
    v = R.vars;
    h = cell(1, numel(idx));
    for k = 1:numel(idx)
        h{k} = matlabFunction(R.ineqs(idx(k)).f, 'Vars', {v(1), v(2)});
    end
end

function P = boxSample(R, n)
% Candidates over a box around R's finite vertices, padded so an unbounded region gets points
% out along its rays.
    vx = double(R.vx); vy = double(R.vy);
    ok = isfinite(vx) & isfinite(vy) & abs(vx) < 1e6 & abs(vy) < 1e6;
    if ~any(ok)
        lo = [-5 -5]; hi = [5 5];
    else
        lo = [min(vx(ok)) min(vy(ok))]; hi = [max(vx(ok)) max(vy(ok))];
        pad = max(1, 3*max(hi - lo)); lo = lo - pad; hi = hi + pad;
    end
    rng(20260819);
    P = lo + (hi - lo) .* rand(n, 2);
end

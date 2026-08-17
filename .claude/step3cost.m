% Step 3's cost, fold by fold, on the case that measures it: f = x*y over the general convex
% quadrilateral conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}, which the A.4/A.5 split turns into 6 pieces.
%
% WHAT THIS IS FOR. Steps 1 and 2 take ~25 s for all six pieces; functionNDomain.maxOfList then
% takes 73 minutes, with the cell count running 5, 14, 29, 45, 70, 86 -- roughly ten times what
% the answer needs (f* of x*y over a convex quadrilateral has a cone per vertex and a cell per
% edge). This script reports, per fold:
%
%   * the cell count after mtimes (the pairing), after maximumP's split loop, and after each of
%     the two mergeL passes -- so growth can be attributed to cells being MADE or merges REFUSED;
%   * the number of DISTINCT functions among the surviving cells, which is the count the answer
%     actually needs: any excess over it is cells that ought to have merged;
%   * region.mergeTally's refusal reasons for that fold.
%
% Run:  CCA2_STEP3_FOLDS=3 matlab -batch "run('.claude/step3cost.m')"
% CCA2_STEP3_FOLDS bounds the number of folds (default: all of them).

DIR = getenv('CCA2DIR');
if isempty(DIR), DIR = pwd; end
cd(DIR);
fprintf('=== repo: %s ===\n', DIR);

setenv('CCA2_A45_SPLIT', '1');
setenv('CCA2_TRACE_MAXP', '1');

% The pipeline emits megabytes of symbolic:sym:isAlways:TruthUnknown warnings on this input --
% enough to push the numbers this script exists to report out of any capped log. They are not
% the measurement; silence them here and nowhere else.
warning('off', 'symbolic:sym:isAlways:TruthUnknown');
warning('off', 'symbolic:isAlways:TruthUnknown');
% Every reported line also goes to its own file, so a truncated console never loses the result.
logPath = fullfile(DIR, '.claude', 'step3cost.log');
logFid = fopen(logPath, 'w');

maxFolds = str2double(getenv('CCA2_STEP3_FOLDS'));
if isnan(maxFolds), maxFolds = inf; end

% CCA2_STEP3_CASE picks the input. 'quad' (default) is the general convex quadrilateral the
% A.4/A.5 split turns into 6 pieces -- surd coordinates throughout. 'tri' is the all-RATIONAL
% reference triangle, which needs no split: it is the control that says whether merge's refusals
% are about exactness or about the gates themselves.
switch getenv('CCA2_STEP3_CASE')
    case 'tri'
        setenv('CCA2_A45_SPLIT', '');
        V = [0 0; 3 3; 1 2];
        E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
    otherwise
        V = [0 0; 2 0; 2.5 1.5; 0.5 1];
        E = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; F = [1 0; 1 0; 1 0; 1 0];
end
q = QuaPol(V, E, [0 1 0 0 0 0], F);

tAll = tic;
p = quaPolToPlq(q);
p = p.triangulate;
both(logFid, 'pieces: %d\n', p.nPieces);
for i = 1:p.nPieces
    tt = tic;
    p.pieces(i) = p.pieces(i).convexEnvelope;
    p.pieces(i) = p.pieces(i).conjugate;
    p.pieces(i) = p.pieces(i).maximumConjugate;
    both(logFid, 'piece %d: %d cells (%.0f s)\n', i, size(p.pieces(i).maxConjugate,2), toc(tt));
end
both(logFid, 'Steps 1+2 total %.0f s\n', toc(tAll));

acc = p.pieces(1).maxConjugate;
both(logFid, 'FOLD 0 (seed): cells=%d distinctF=%d\n', size(acc,2), nDistinct(acc));

for j = 2:min(p.nPieces, 1 + maxFolds)
    region.mergeTally('reset');
    t1 = tic;  acc = acc * p.pieces(j).maxConjugate;              tMul = toc(t1);
    nPaired = size(acc,2);
    t2 = tic;  acc = acc.maximumP(true);                          tMax = toc(t2);
    both(logFid, 'FOLD %d: paired=%3d -> cells=%3d distinctF=%3d   mtimes %.0f s, maximumP %.0f s\n', ...
        j-1, nPaired, size(acc,2), nDistinct(acc), tMul, tMax);
    both(logFid, '  merge: %s\n', tallyString(region.mergeTally('get')));
end
both(logFid, 'TOTAL %.0f s\n', toc(tAll));
fclose(logFid);
fprintf('summary also written to %s\n', logPath);

function n = nDistinct(acc)
% How many DIFFERENT functions the surviving cells carry. mergeL groups by exactly this test, so
% it is also the number of cells a perfect merge would leave.
    n = 0; reps = {};
    for i = 1:numel(acc)
        f = acc(i).f(1).f;
        isNew = true;
        for k = 1:numel(reps)
            if isAlways(simplifyFraction(f - reps{k}) == 0, 'Unknown', 'false')
                isNew = false; break
            end
        end
        if isNew, n = n + 1; reps{end+1} = f; end %#ok<AGROW>
    end
end

function s = tallyString(T)
    k = fieldnames(T);
    if isempty(k), s = 'no attempts'; return, end
    s = '';
    for i = 1:numel(k)
        s = [s sprintf('%s=%d ', k{i}, T.(k{i}))]; %#ok<AGROW>
    end
    s = strtrim(s);
end

function both(fid, varargin)
% Print to the console AND to the run's own log, so a capped console never loses the numbers.
    fprintf(varargin{:});
    if fid > 0, fprintf(fid, varargin{:}); end
end
